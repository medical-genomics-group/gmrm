#!/bin/bash

module load intel intel-mpi boost

sh comp_intel.sh

srun -n 2 -p build --mem 4G ./bin/ardyh_i \
--bedfile /work/ext-unil-ctgg/marion/benchmark_simulation/data/ukb_chr2_N_QC.bed \
--dimfile /work/ext-unil-ctgg/marion/benchmark_simulation/data/ukb_chr2_N_QC.dim \
--phenfiles /work/ext-unil-ctgg/marion/benchmark_simulation/phen/sim_1/data.noHEAD.phen,/work/ext-unil-ctgg/marion/benchmark_simulation/phen/sim_1/data.noHEAD.phen \
--shuffle-markers 0 \
--seed 123 \
--trunc-markers 200 \
--S 0.0001,0.001,0.01 \
--verbosity 2 \
--iterations 20


