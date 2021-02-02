#!/bin/bash

set -e
#set -x

source ../compile_with_intel.sh $1

echo ARDYH_ROOT = $ARDYH_ROOT
echo ARDYH_EXE  = $ARDYH_EXE

BENCH_DIR=/work/ext-unil-ctgg/etienne/data_bench
[ -d $BENCH_DIR ] || (echo "fatal: bench directory not found! $BENCH_DIR" && exit)
echo BENCH_DIR = $BENCH_DIR

OUT_DIR=/scratch/orliac/bench_ardyh

NTASKS=6
NTHREADS_PER_TASK=6

export OMP_NUM_THREADS=$NTHREADS_PER_TASK


CMD_BASE="srun -p build -n $NTASKS --cpus-per-task $NTHREADS_PER_TASK -t 00:10:00 --mem=0 --cpu-bind=verbose \
../bin/$ARDYH_EXE \
--bed-file $BENCH_DIR/test.bed \
--dim-file $BENCH_DIR/test.dim \
--group-index-file $BENCH_DIR/test.gri \
--group-mixture-file $BENCH_DIR/test.grm \
--shuffle-markers 0 \
--seed 123 \
--trunc-markers 10000 \
--verbosity 2 \
--iterations 5"

PHENS1="--phen-files $BENCH_DIR/test1.phen"
PHENS2="--phen-files $BENCH_DIR/test1.phen,$BENCH_DIR/test2.phen"
PHENS3="--phen-files $BENCH_DIR/test1.phen,$BENCH_DIR/test2.phen,$BENCH_DIR/test3.phen"
PHENS4="--phen-files $BENCH_DIR/test1.phen,$BENCH_DIR/test2.phen,$BENCH_DIR/test3.phen,$BENCH_DIR/test4.phen"


CMD=${CMD_BASE}" "${PHENS1}
echo CMD = $CMD
$CMD

exit 0

CMD=${CMD_BASE}" "${PHENS2}
echo CMD = $CMD
$CMD

CMD=${CMD_BASE}" "${PHENS3}
echo CMD = $CMD
$CMD

CMD=${CMD_BASE}" "${PHENS4}
echo CMD = $CMD
$CMD




