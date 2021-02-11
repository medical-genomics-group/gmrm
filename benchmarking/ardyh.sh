#!/bin/bash

module load intel intel-mpi boost

sh comp_intel.sh || exit 1

srun -n 4 --cpus-per-task 8 -p build --mem 40G ./bin/ardyh_i \
--bedfile /work/ext-unil-ctgg/marion/benchmark_simulation/data/ukb_chr2_N_QC.bed \
--dimfile /work/ext-unil-ctgg/marion/benchmark_simulation/data/ukb_chr2_N_QC.dim \
--phenfiles /work/ext-unil-ctgg/marion/benchmark_simulation/phen/sim_1/data.noHEAD.phen \
--group-index-file /work/ext-unil-ctgg/marion/benchmark_simulation/ukb_chr2_N_QC_groups7.group \
--group-mixture-file /work/ext-unil-ctgg/marion/benchmark_simulation/00001_0001_001_01_groups7.cva \
--shuffle-markers 0 \
--seed 123 \
--trunc-markers 57700 \
--S 0.0001,0.001,0.01 \
--verbosity 2 \
--iterations 20
