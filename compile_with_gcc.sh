#!/bin/bash

set -e

# before loading modules!!!
export ARDYH_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo ARDYH_ROOT = $ARDYH_ROOT
[[ "$ARDYH_ROOT" == *ardyh ]] || (echo "ARDYH_ROOT ($ARDYH_ROOT) expected to end with \"ardyh\"" && exit 1)

module purge
module load gcc/8.3.0 mvapich2 boost
module list

cd $ARDYH_ROOT
ARDYH_EXE=ardyh_g
EXEC=$ARDYH_EXE make -j 8 $1 || exit 1
export ARDYH_EXE=$ARDYH_EXE
cd -
