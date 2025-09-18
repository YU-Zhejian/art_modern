FROM almalinux/9-base
# Install necessary packages
RUN dnf update -y && \
    dnf install -y \
    gcc \
    gcc-c++ \
    make \
    cmake \
    python3 \
    boost-devel \
    zlib-devel \
    pkgconf-pkg-config \
    coreutils \
    sed \
    grep \
    fmt-devel \
    fmt \
    moodycamel-concurrentqueue-devel \
    abseil-cpp-devel \
    abseil-cpp \
    rpm-build \
    rpm-devel \
    rpmlint \
    diffutils \
    patch \
    rpmdevtools && \
    dnf clean all
# HTSLib have to be bundled. Ridiculous.
