#!/usr/bin/env bash
# Shell script that counts how many lines of code was written by us.
# shellcheck disable=SC2086

builtin set -ue
NAME="scc.sh"
VERSION=0.1

if builtin type -p "${SCC:-}" &>/dev/null; then
    true
elif builtin type -p scc &>/dev/null; then
    SCC=scc
elif builtin type -p cloc &>/dev/null; then
    SCC=cloc
else
    CLOC_INFO="scc or cloc required!"
    builtin exit 1
fi

# SHDIR="$(dirname "$(readlink -f "${0}")")"

LAST_COMMIT=$(git log --pretty=oneline --abbrev-commit --graph --branches -n 1)
AUTHOR_INFO=$(git shortlog --numbered --summary --email)
AD_MINUS=$(
    git log --numstat --pretty="%an$(echo -e "\t")%H" |
        awk '
    BEGIN{
        FS="\t"
    }
    {
        if (NF == 2){
            name = $1
        };
        if(NF == 3) {
            plus[name] += $1; minus[name] += $2
        }
    }
    END {
        for (name in plus) {
            print name":\t+"plus[name]"\t-"minus[name]
        }
    }' |
        sort -k2 -gr |
        sed 's;^;\t;'
)

SOURCES=$(
    git ls-files |
        grep -v '\.idea' |
        grep -v 'deps' |
        grep -v 'Illumina_profiles' |
        grep -v 'benchmark_other_simulators/src' |
        xargs
)
echo "Enumerating sources FIN"

CLOC_INFO=$("${SCC}" ${SOURCES})

cat <<EOF
${NAME} ver. ${VERSION}
Called by: ${0} ${*}
Repository version information:
	The last commit is: ${LAST_COMMIT}
Author Information:
${AUTHOR_INFO}
Author changes:
${AD_MINUS}
Code count:
${CLOC_INFO}
EOF

builtin exit 0
