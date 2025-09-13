#!/usr/bin/env bash
# shellcheck disable=SC2086
set -ue
export DEBUILD_LINTIAN_OPTS="-i -I --show-overrides --pedantic"

SRC_DIR="$(dirname "$(readlink -f "${0}")")/../"
cd "${SRC_DIR}" || exit 1
rm -fr opt/build_deb
mkdir -p opt/build_deb

echo "DEB: Building for version ${PACKAGE_VERSION}"

# Generate orig tarball from git files
echo "DEB: Generate orig tarball from git files"
git ls-files | grep -v '^debian/' > opt/build_deb/git-files.txt
mkdir -p opt/build_deb/art-modern-${PACKAGE_VERSION}
while read -r line; do
    if [ -f "${line}" ]; then
        install -D "${line}" "opt/build_deb/art-modern-${PACKAGE_VERSION}/${line}"
    fi
done < opt/build_deb/git-files.txt
env -C opt/build_deb/ tar -czf art-modern_${PACKAGE_VERSION}.orig.tar.gz art-modern-${PACKAGE_VERSION}

# Copy debian files
echo "DEB: Copy debian files"
cp -r debian opt/build_deb/art-modern-${PACKAGE_VERSION}/debian

# Build the deb
echo "DEB: Build the deb"
env -C opt/build_deb/art-modern-${PACKAGE_VERSION} debuild -us -uc \
    > >( sed 's/^/DEBUILD-INFO: /' | tee opt/build_deb/debuild-info.log) \
    2> >( sed 's/^/DEBUILD-ERR: /' | tee opt/build_deb/debuild-err.log >&2 )
env -C opt/build_deb/art-modern-${PACKAGE_VERSION} debuild -- clean
# TODO: Supress -Wdate-time
