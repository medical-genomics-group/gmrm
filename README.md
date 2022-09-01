The contents of the branch:

- All the previous files from gmrm ./src repository with changes made
- boost_1_79_0 library should be downloaded
- ./example folder - folder with original test data. File test_gri.gri was added, it specifies the division of markers into 2 groups.
- 4 folders with results_{num of groups}_{num of iters} (for the original test data from ./example folder) with different iterations number and different number of groups
- When the gmrm.exe runs and produces the results, the ./results folder is created with 3 files: result_betas (trace for 2000 betas for each phenotype, two neighbor cols correspond to one marker and 2 phenotypes), result_rescov (trace for values of residual covariance matrix), result_betascov (trace for effects covariance matrix for each group)
- ./src/simulation_stats.R - R script, taking data from ./results folder and computing medians for it.
