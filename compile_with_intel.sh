#!/bin/bash

set -e

# before loading modules!!!
export ARDYH_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo ARDYH_ROOT = $ARDYH_ROOT
[[ "$ARDYH_ROOT" == *ardyh ]] || (echo "ARDYH_ROOT ($ARDYH_ROOT) expected to end with \"ardyh\"" && exit 1)

module purge
module load intel intel-mpi boost
module list

cd $ARDYH_ROOT
ARDYH_EXE=ardyh_i
EXEC=$ARDYH_EXE make -j8 $1 || exit 1
export ARDYH_EXE=$ARDYH_EXE
cd -
