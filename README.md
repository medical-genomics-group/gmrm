# ardyh



## Compilation 

Move into `ardyh` cloned repository:

```
cd /path/to/ardyh
```

Then, assuming a standard HPC environment:

```
module load intel intel-mpi boost eigen
make
```

or

```
module load gcc mvapich boost eigen
make
```



## Running

Assuming Slurm as job scheduler, create a file `sbatch_ardyh.sh` along these lines:

```
#!/bin/bash

#SBATCH --partition debug
#SBATCH --tasks 2
#SBATCH --cpus-per-task 2
#SBATCH --memory 4G

module load gcc mvapich boost eigen 

srun ./bin/ardyh_g \
--bedfile /work/ext-unil-ctgg/marion/benchmark_simulation/data/ukb_chr2_N_QC.bed \
--dimfile /work/ext-unil-ctgg/marion/benchmark_simulation/data/ukb_chr2_N_QC.dim \
--phenfiles /work/ext-unil-ctgg/marion/benchmark_simulation/phen/sim_1/data.noHEAD.phen,/work/ext-unil-ctgg/marion/benchmark_simulation/phen/sim_1/data.noHEAD.phen

```

Then submit it:

`sbatch sbatch_ardyh.sh`



## Considerations on implementing multi-traits processing

The markers' statistics are phenotype dependent (due to the NAs). When reading in the phenotypes, a bitmask is created, indicating NA with 0. A byte therefore contains the information for 8 individuals.