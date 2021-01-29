#!/bin/bash

set -e
#set -x

source ../compile_with_gcc.sh $1

echo ARDYH_ROOT = $ARDYH_ROOT
echo ARDYH_EXE  = $ARDYH_EXE

BENCH_DIR=/work/ext-unil-ctgg/etienne/data_bench/
[ -d $BENCH_DIR ] || (echo "fatal: bench directory not found! $BENCH_DIR" && exit)
echo BENCH_DIR = $BENCH_DIR

OUT_DIR=/scratch/orliac/bench_ardyh

NTHREADS=4

export OMP_NUM_THREADS=$NTHREADS
export MV2_ENABLE_AFFINITY=0


CMD_BASE="srun -n 4 --cpus-per-task $NTHREADS -p build --mem 40G \
../bin/$ARDYH_EXE \
--bedfile $BENCH_DIR/test.bed \
--dimfile $BENCH_DIR/test.dim \
--phenfiles $BENCH_DIR/test.phen \
--shuffle-markers 0 \
--seed 123 \
--trunc-markers 10000 \
--S 0.0001,0.001,0.01 \
--verbosity 2 \
--iterations 3"

#BED
CMD=${CMD_BASE}
echo CMD = $CMD
$CMD




