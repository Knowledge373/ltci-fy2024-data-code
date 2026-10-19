#!/bin/sh
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
export R_LIBS_USER="$ROOT/r-library"
export R_MAKEVARS_USER="$ROOT/Makevars.low-memory"
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export TZ=UTC
export PATH="$ROOT/toolchain/usr/bin:$PATH"
export LIBRARY_PATH="$ROOT/toolchain/lib:${LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec "$ROOT/r-runtime/usr/lib/R/bin/R" --vanilla --slave "$@"
