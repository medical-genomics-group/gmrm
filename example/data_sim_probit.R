## simple R script to simulate example genotype data
## MRR 14.07.21
## This requires the software plink: https://www.cog-genomics.org/plink2

set.seed(171014)
require(MASS)

## set sample size, N
N = 10000

## total number of covariates, M
M = 20000

## simulate marker data
X <- matrix(rbinom(N*M,2,0.4),N,M)

## total variance explained by marker effects
h2 = 0.5

## variance of the marker effects
## only 500 causal markers
cm = 500
sigma_b = h2 / cm

## sample marker effects
index <- sample(1:M, cm)

## generate genetic values
beta <- rep(0,M)
b <- mvrnorm(cm, 0.0, sigma_b)
beta[index] = b#rep(0,cm)#b
g <- scale(X) %*% beta

## generate residuals
sigma_e = 1 - var(g)
e <- mvrnorm(N, 0, sigma_e)

## output phenotype
z = g + e
y <- z
y[y > 0] = 1
y[y <= 0] = 0

## output genetic data
X[X == 2] <- "AA"
X[X == 1] <- "AG"
X[X == 0] <- "GG"

## output to plink .ped/.map format
ped <- data.frame("FID" = 1:N,
                  "IID" = 1:N,
                  "PID" = rep(0,N),
                  "MID" = rep(0,N),
                  "Sex" = rep(1,N),
                  "phen" = rep(0,N))
ped <- cbind(ped,X)
write.table(ped,"test.ped", row.names=FALSE, col.names=FALSE, quote=FALSE)

map <- data.frame("chr" = rep(1,M),
                  "rs" = paste("rs",1:M, sep=''),
                  "dist" = rep(0,M),
                  "bp" = 1:M)
write.table(map,"test.map", row.names=FALSE, col.names=FALSE, quote=FALSE)

## convert from .ped/.map to plink binary format
system("plink --file test --make-bed --out test")

## remove .ped/.map files
system("rm *.ped")
system("rm *.map")
system("rm *.log")

## output phenotype files
phen <- data.frame("FID" = 1:N,
                   "IID" = 1:N,
                   "phen" = y[,1],
                   "latent" = z[,1])
write.table(phen[,c(1,2,3)],"test.phen", row.names=FALSE, col.names=FALSE, quote=FALSE)
#write.table(phen[,c(1,2,4)],"test.latent", row.names=FALSE, col.names=FALSE, quote=FALSE)