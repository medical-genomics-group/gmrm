#!/bin/bash

set -e
#set -x

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo DIR = $DIR

source $DIR/../compile_with_gcc.sh $1

PGM_ROOT=$ARDYH_ROOT
PGM_EXE=$ARDYH_EXE

BENCH_DIR=/work/ext-unil-ctgg/etienne/data_bench
[ -d $BENCH_DIR ] || (echo "fatal: bench directory not found! $BENCH_DIR" && exit)
echo BENCH_DIR = $BENCH_DIR

OUT_DIR=/scratch/orliac/bench_${PGM_EXE}

NTASKS=2
NTHREADS_PER_TASK=18
TPS="--ntasks-per-socket 1"

export OMP_NUM_THREADS=$NTHREADS_PER_TASK
export MV2_ENABLE_AFFINITY=0

CMD_BASE="srun -p build -n $NTASKS $TPS --cpus-per-task $NTHREADS_PER_TASK -t 00:10:00 --mem=0 --cpu-bind=verbose \
$PGM_ROOT/bin/$PGM_EXE \
--bed-file $BENCH_DIR/test.bed \
--dim-file $BENCH_DIR/test.dim \
--group-index-file $BENCH_DIR/test.gri \
--group-mixture-file $BENCH_DIR/test.grm \
--shuffle-markers 0 \
--seed 123 \
--trunc-markers 10000 \
--verbosity 2 \
--iterations 12 \
--mimic-hydra"

PHEN1=$BENCH_DIR/test1.phen
PHEN2=$BENCH_DIR/test2.phen
PHEN3=$BENCH_DIR/test3.phen
PHEN4=$BENCH_DIR/test4.phen
PHEN5=$BENCH_DIR/test5.phen

PHENS1="--phen-files $PHEN1"
PHENS5="--phen-files $PHEN5"
PHENS1_2="--phen-files $PHEN1,$PHEN2"
PHENS1_3="--phen-files $PHEN1,$PHEN2,$PHEN3"
PHENS1_4="--phen-files $PHEN1,$PHEN2,$PHEN3,$PHEN4"
PHENS1_5="--phen-files $PHEN1,$PHEN2,$PHEN3,$PHEN4,$PHEN5"
PHENS15="--phen-files $PHEN1,$PHEN5"


CMD=${CMD_BASE}" "${PHENS1}; echo CMD = $CMD; $CMD; exit 0
CMD=${CMD_BASE}" "${PHENS5}; echo CMD = $CMD; $CMD; #exit 0
CMD=${CMD_BASE}" "${PHENS15}; echo CMD = $CMD; $CMD; exit 0
#CMD=${CMD_BASE}" "${PHENS3}; echo CMD = $CMD; $CMD;
#CMD=${CMD_BASE}" "${PHENS3}; echo CMD = $CMD; $CMD;
#CMD=${CMD_BASE}" "${PHENS4}; echo CMD = $CMD; $CMD;
