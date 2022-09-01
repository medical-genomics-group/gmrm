# Read simulation files
rescov <- read.csv(file = '../results/result_rescov.csv')
betascov <- read.csv(file = '../results/result_betascov.csv')
betas <- read.csv(file = '../results/result_betas.csv')

# burnin
burnin = 1000

#Find medians of all simulations for all columns
rescov_median = apply(rescov[c(burnin:nrow(rescov)),c(2:ncol(rescov))],2,median)
betascov_median = apply(betascov[c(burnin:nrow(betascov)),c(2:ncol(betascov))],2,median)
betas_medians = apply(betas[c(burnin:nrow(betas)),c(2:ncol(betas))],2,median)


# Put estimated beta medians of betas into dataframe
df_betas_medians = data.frame(matrix(NA, nrow = (ncol(betas)-1)/2, ncol = 2))
nrow = 1
for (i in seq(1, (ncol(betas)-1), 2)){
  df_betas_medians[nrow,] = betas_medians[i:(i+1)]
  nrow = nrow+1
}

# Number of betas non-zero medians
colSums(df_betas_medians != 0)
