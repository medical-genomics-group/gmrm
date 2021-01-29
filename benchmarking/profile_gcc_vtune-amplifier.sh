#!/bin/bash

#SBATCH --partition build
#SBATCH --mem 10GB
#SBATCH --cpus-per-task 2

set -e
#set -x

COMPILER=gcc

print_help() {
cat <<-HELP

Call the script like this:

$0 --compiler {gcc / intel}, default is gcc

HELP
exit 0
}

# Parse Command Line Arguments
while [ "$#" -gt 0 ]; do
    case "$1" in
        --compiler*|-c*)
            if [[ "$1" != *=* ]]; then shift; fi
            COMPILER="${1#*=}"
            ;;
        --help|-h)
            print_help;;
        *)
            printf "************************************************************\n"
            printf "* Error: Invalid argument, run --help for valid arguments. *\n"
            printf "************************************************************\n"
            exit 1
            ;;
    esac
    shift
done

if [ $COMPILER == "gcc" ]; then
    echo "GCC"
    source ../compile_with_gcc.sh $1
    #export GOMP_CPU_AFFINITY=verbose
elif [ $COMPILER == "intel" ]; then
    echo "INTEL"
    source ../compile_with_intel.sh $1
    #export KMP_AFFINITY=verbose
else 
    echo "fatal: unknown compiler $COMPILER. Check --help."
fi

echo HYDRA_ROOT = $HYDRA_ROOT
echo HYDRA_EXE  = $HYDRA_EXE

AVARS=/ssoft/spack/external/intel/2018.4/vtune_amplifier/amplxe-vars.sh
[ -f $AVARS ] || (echo "fatal: $AVARS not found." && exit 1)
source $AVARS
export MODULEPATH=/ssoft/spack/humagne/v1/share/spack/lmod/linux-rhel7-x86_S6g1_Mellanox/intel/18.0.5:$MODULEPATH


BENCH_DIR=/work/ext-unil-ctgg/etienne/data_bench/
[ -d $BENCH_DIR ] || (echo "fatal: bench directory not found! $BENCH_DIR" && exit)
echo BENCH_DIR = $BENCH_DIR

OUT_DIR=/scratch/orliac/bench_hydra2
[ -d $OUT_DIR ] && rm -rv $OUT_DIR

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#SEARCH_DIRS="--search-dir src:=$HYDRA_ROOT/src --search-dir sym:=$HYDRA_ROOT/build_gcc --search-dir bin:=$HYDRA_ROOT/bin"

ADVIXE1="amplxe-cl -c hotspots    -r $OUT_DIR/vtune  -data-limit=0 -- " #-no-auto-finalize $SEARCH_DIRS"
ADVIXE1="amplxe-cl -c concurrency -r $OUT_DIR/vtune  -data-limit=0 -- " #-no-auto-finalize $SEARCH_DIRS"
#ADVIXE2="amplxe-cl --collect=tripcounts --flop --project-dir=$OUT_DIR/vtune -no-auto-finalize $SEARCH_DIRS -data-limit=0 --"

CMD_BASE="srun -n 1 -t 00:03:00 --cpu-bind=verbose,sockets"

CMD_TAIL="$HYDRA_ROOT/bin/$HYDRA_EXE \
    --bfile $BENCH_DIR/test \
    --pheno $BENCH_DIR/test.phen \
    --mcmc-out-dir $OUT_DIR \
    --mcmc-out-name bench_hydra_epfl_$COMPILER \
    --mpibayes bayesMPI \
    --shuf-mark 0 \
    --seed 123 \
    --number-markers 20 \
    --number-individuals 500000 \
    --S 0.0001,0.001,0.01 \
    --verbosity 0 \
    --chain-length 1"

#BED
env | grep SLURM_
env | grep OMP_

CMD="${CMD_BASE} ${ADVIXE1} ${CMD_TAIL} \
    --bfile $BENCH_DIR/test"
echo CMD = $CMD
$CMD

#CMD=${CMD_BASE}" "${ADVIXE2}" "${CMD_TAIL}" \
#    --bfile $BENCH_DIR/test \
#echo CMD = $CMD
#$CMD

