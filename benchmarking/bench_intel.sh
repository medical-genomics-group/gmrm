#!/bin/bash

set -e
#set -x

source ../compile_with_intel.sh $1

echo ARDYH_ROOT = $ARDYH_ROOT
echo ARDYH_EXE  = $ARDYH_EXE

BENCH_DIR=/work/ext-unil-ctgg/etienne/data_bench/
[ -d $BENCH_DIR ] || (echo "fatal: bench directory not found! $BENCH_DIR" && exit)
echo BENCH_DIR = $BENCH_DIR

OUT_DIR=/scratch/orliac/bench_ardyh

NTASKS=7
NTHREADS_PER_TASK=4

export OMP_NUM_THREADS=$NTHREADS_PER_TASK


CMD_BASE="srun -p build -n $NTASKS --cpus-per-task $NTHREADS_PER_TASK -t 00:10:00 --mem=0 --cpu-bind=verbose \
../bin/$ARDYH_EXE \
--bedfile $BENCH_DIR/test.bed \
--dimfile $BENCH_DIR/test.dim \
--phenfiles $BENCH_DIR/test.phen \
--shuffle-markers 0 \
--seed 123 \
--trunc-markers 10 \
--S 0.0001,0.001,0.01 \
--verbosity 2 \
--iterations 20"

#BED
CMD=${CMD_BASE}
echo CMD = $CMD
$CMD




