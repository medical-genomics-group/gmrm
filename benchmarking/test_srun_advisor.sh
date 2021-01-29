#!/bin/bash

#SBATCH --partition debug

set -e
#set -x 

source ../compile_with_gcc.sh $1
#module load intel

ADVIXE_VARS=/ssoft/spack/external/intel/2018.4/advisor/advixe-vars.sh
[ -f $ADVIXE_VARS ] || (echo "fatal: $ADVIXE_VARS not found." && exit 1)
source $ADVIXE_VARS

which advixe-cl

export MODULEPATH=/ssoft/spack/humagne/v1/share/spack/lmod/linux-rhel7-x86_S6g1_Mellanox/intel/18.0.5:$MODULEPATH
#prepend_path("PATH","/ssoft/spack/external/intel/2018.4/compilers_and_libraries_2018.5.274/linux/bin")
#prepend_path("MANPATH","/ssoft/spack/external/intel/2018.4/compilers_and_libraries_2018.5.274/linux/man")
#prepend_path("CMAKE_PREFIX_PATH","/ssoft/spack/external/intel/2018.4/compilers_and_libraries_2018.5.274/linux/")
#setenv("INTEL_ROOT","/ssoft/spack/external/intel/2018.4/compilers_and_libraries_2018.5.274/linux")
#setenv("CC","icc")
#setenv("CXX","icpc")
#setenv("FC","ifort")
#setenv("F90","ifort")
#setenv("F77","ifort")
#setenv("AMPLXE_ROOT","/ssoft/spack/external/intel/2018.4/vtune_amplifier_2018.4.0.573462/bin64")
#setenv("ICCCFG","/ssoft/spack/external/intel/config/2018.4/compilers_and_libraries_2018.5.274/icl.cfg")
#setenv("ICPCCFG","/ssoft/spack/external/intel/config/2018.4/compilers_and_libraries_2018.5.274/icl.cfg")
#setenv("IFORTCFG","/ssoft/spack/external/intel/config/2018.4/compilers_and_libraries_2018.5.274/icl.cfg")
#prepend_path("PATH","/ssoft/spack/external/intel/2018.4/vtune_amplifier_2018.4.0.573462/bin64:/ssoft/spack/external/intel/2018.4/compilers_and_libraries_2018.5.274/linux/bin/intel64")
#prepend_path("INTEL_LICENSE_FILE","/ssoft/spack/external/intel/License")



#srun -p debug advixe-cl -- $HYDRA_ROOT/bin/$HYDRA_EXE

srun -n 1 -t 00:03:00 --mem 5GB --cpus-per-task 1  advixe-cl --collect=survey --project-dir=/scratch/orliac/bench_hydra/advisor -no-auto-finalize -- date

#srun -n 1 -t 00:03:00 --mem 5GB --cpus-per-task 1 --cpu-bind=verbose,sockets advixe-cl --collect=survey --project-dir=/scratch/orliac/bench_hydra/advisor -no-auto-finalize -- /home/orliac/DCSR/CTGG/hydra/bin/hydra_g --bfile /work/ext-unil-ctgg/etienne/data_bench//test --pheno /work/ext-unil-ctgg/etienne/data_bench//test.phen --mcmc-out-dir /scratch/orliac/bench_hydra --mcmc-out-name bench_hydra_epfl_gcc --mpibayes bayesMPI --shuf-mark 0 --seed 123 --number-markers 20 --number-individuals 500000 --S 0.0001,0.001,0.01 --verbosity 0 --chain-length 1 --bfile /work/ext-unil-ctgg/etienne/data_bench//test --mcmc-out-name profile_hydra_epfl_gcc_bed
