#!/usr/bin/env bash
set -ue
cd "$(readlink -f "$(dirname "${0}")")"
if [ -f "sh.d/profile.sh.d/${1:-}.sh" ]; then
    bash --norc --noprofile "sh.d/profile.sh.d/${1}.sh"
else
    echo "Unknown profile: ${1:-}" >&2
    echo "Available profilers:"
    find "sh.d/profile.sh.d" -type f -name '*.sh' -printf '%f\n' | sed 's/\.sh$//'
fi
