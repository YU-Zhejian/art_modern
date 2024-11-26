#!/usr/bin/env bash
set +ue
. /opt/intel/oneapi/setvars.sh
set -ue
env -C opt/art_modern_build/ ninja
