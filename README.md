# ardyh



## Compilation 

Move into `ardyh` cloned repository:

```
cd /path/to/ardyh
```

Then, assuming a standard HPC environment:

```
module load intel intel-mpi boost
make
```

or

```
module load gcc mvapich boost
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
--phenfiles /work/ext-unil-ctgg/marion/benchmark_simulation/phen/sim_1/data.noHEAD.phen,/work/ext-unil-ctgg/marion/benchmark_simulation/phen/sim_1/data.noHEAD.phen \
--shuffle-markers 1 \
--seed 12345 \
--iterations 10
--trunc-markers 5
--S 0.0001,0.001,0.01
--verbosity 2

```

Then submit it:

`sbatch sbatch_ardyh.sh`



## Processing options

| options           |      | values                                  | description                                                  |
| ----------------- | ---- | --------------------------------------- | ------------------------------------------------------------ |
| --bedfile         |      | /path/to/gen                            | Path to the PLINK bed file to be processed.                  |
| --phenfiles       |      | /path/to/phen1,<br />/path/to/phen2,... | Comma separated (**no space!**) list of phenotype files to be processed. <br />At least one phenotype file is expected. |
| --dimfile         |      | /path/to/dim                            | Path to the file containing the dimensions of the genotype: expected <br />to be a single line file containing 2 integers: N and M, N being the number<br />of individuals and M the number of markers. |
| --verbosity       |      | 0, 1, 2, 3                              | default is 0, printing the bare minimum information on the processing.<br />Level 3 is the debugging mode, printing extensively. |
| --shuffle-markers |      | 0, 1                                    | Shuffling (1) or not (0) the markers.                        |
| --seed            |      | unsigned integer                        | Seed for pseudo-random number generator (boost).             |
| --iterations      |      | unsigned integer >= 1                   | Number of iterations to run.                                 |
| --trunc-markers   |      | unsigned integer >= 1                   | Truncate the number of markers to process.                   |
| --S               |      | 0.0001,0.001,0.01                       | Component assignment; all > 0.0; in increasing order!        |
|                   |      |                                         |                                                              |



## Considerations on implementing multi-traits processing

1. The markers' statistics are phenotype dependent (due to the NAs). A mask is applied using a lookup tables with 64 x 4 entries. Note that mostly the same information will be reused (0b000011111) as phenotypes contains few NAs.
2. The same base seed is used when processing the different phenotype, that is, the markers for each phenotype are treated in the same way. Therefore, this multi-trait version can not be used to process several times the same phenotype.