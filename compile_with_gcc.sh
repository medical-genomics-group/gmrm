#!/bin/bash

set -e

module purge
module load gcc mvapich2 boost
module list

export ARDYH_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." >/dev/null 2>&1 && pwd )"

cd $ARDYH_ROOT

ARDYH_EXE=ardyh_g

EXEC=$ARDYH_EXE make -j8 $1 || exit 1

export ARDYH_EXE=$ARDYH_EXE

cd -
