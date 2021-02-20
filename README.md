# ardyh



## Compilation 

```
git clone https://github.com/medical-genomics-group/ardyh.git
cd ardyh
```

Then, assuming a standard HPC environment, e.g.

```
module load intel intel-mpi boost
make
```

or

```
module load gcc mvapich boost
make
```

Edit the Makefile if necessary to fit your environment.



## Running

Assuming Slurm as job scheduler, create a file `sbatch_ardyh.sh` along these lines:

```
#!/bin/bash

#SBATCH --ntasks 4
#SBATCH --cpus-per-task 6
#SBATCH --memory 10G
#SBATCH --time 01:00:00

module load gcc mvapich2 boost

srun $ardyh/bin/ardyh_g \
--bed-file /path/to/file.bed \
--dim-file /path/to/file.dim \
--phen-files /path/to/phen1.phen,/path/to/phen2.phen \
--shuffle-markers 1 \
--seed 12345 \
--iterations 10 \
--trunc-markers 5 \
--verbosity 2 \

```

Then submit it:

`sbatch sbatch_ardyh.sh`



## Processing options

| options              | values                                            | description                                                  |
| -------------------- | ------------------------------------------------- | ------------------------------------------------------------ |
| --bed-file           | /path/to/file.bed                                 | Path to the PLINK bed file to be processed.                  |
| --phen-files         | /path/to/file1.phen,<br />/path/to/file2.phen,... | Comma separated (**no space!**) list of phenotype files to be processed. <br />At least one phenotype file is expected. |
| --dim-file           | /path/to/file.dim                                 | Path to the file containing the dimensions of the genotype: expected <br />to be a single line file containing 2 integers: N and M, N being the number<br />of individuals and M the number of markers. |
| --group-index-file   | /path/to/file.gri                                 | Path to group index file. One line for each marker. On each line a label<br />and the group to which the markers belongs to. Each marker must<br />belong to a group. Groups start at 0 and are successive. Like 0,1,2,3 if<br />you have 4 groups. At least a single group with all markers. |
| --group-mixture-file | /path/to/file.grm                                 | Path to group mixture file. For each group, a line describing the component<br />mixture. Must start with 0.0. Next entry must be strictly greater. For all<br />groups, the same number of components is required. |
| --verbosity          | 0, 1, 2, 3                                        | Default is 0, printing the bare minimum information on the processing.<br />Level 3 is the debugging mode, printing extensively. |
| --shuffle-markers    | 0, 1                                              | Shuffling (1) or not (0) the markers. Default is shuffling.  |
| --seed               | unsigned integer                                  | Seed for pseudo-random number generator (Mersenne twister).  |
| --iterations         | unsigned integer >= 1                             | Number of iterations to run.                                 |
| --trunc-markers      | unsigned integer >= 1                             | Truncate the number of markers to process.                   |
| --mimic-hydra        |                                                   | If this option is passed a single PRNG is used as per in *hydra*. Can only process<br />a single phenotype. |
|                      |                                                   |                                                              |



## Considerations on implementing multi-traits processing

1. Because identical results are expected for a given phenotype whether it is processed alone or in multi-traits mode, we had to decouple the shuffling of the markers from the pseudo-random number generator (PRNG) used to sample distributions. This is imposed by the fact that we want to process the markers against all phenotypes in the same order on each iteration, in order to optimize the synchronization between tasks if a marker triggers an update on several phenotypes. Hence the option --mimic-hydra to reproduce the results from *hydra*.
2. The PRNG for the shuffling of the markers is seeded with --seed + rank.
3. The PRNG used for sampling distributions for each phenotype is seeded with --seed + (rank + 1)  x 1000.
4. The markers' statistics are phenotype dependent (due to the NAs). A mask is applied using a lookup tables with 64 x 4 entries. Note that mostly the same information will be reused (0b000011111) as phenotypes contains few NAs.
5. The same base seed is used when processing the different phenotype, that is, the markers for each phenotype are treated in the same way. Therefore, this multi-trait version can not be used to process several times the same phenotype.