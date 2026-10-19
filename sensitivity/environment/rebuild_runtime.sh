#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$ROOT"
python3 prepare_runtime.py
./run-r-local.sh -f rebuild_from_preserved_sources.R > runtime_rebuild.log 2>&1
./run-r-local.sh -e "stopifnot(getRversion()=='4.6.1', packageVersion('fwildclusterboot')=='0.14.3', packageVersion('dqrng')=='0.4.1'); print(sessionInfo())"
