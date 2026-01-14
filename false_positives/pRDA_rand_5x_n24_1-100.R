# run pRDA on random variables 3 at a time 30x to achieve 100 lists of random candidate SNPs

# Run in R in directory /srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/04_rda/26_pRDA_rand_5x_n24_1-100
setwd("...")

# in R
# Load packages
library(tidyverse)
library(psych)
library(adespatial)
library(vegan)
library(rlist)
library(adegenet)
library(Hmisc)
library(tidyr)

# Load SNPs
snps_unlinked_5x <-read.PLINK(".../5x_TOA_only_filtered_SNPs_all_unlinked.raw")

# Impute missing SNPs
snps_unlinked_imp_5x <- apply(snps_unlinked_5x, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
# may take a few minutes

# Load environmental variables
envs <- read.csv(".../pRDA_rand_5x_n24_1-100.csv")

# RANDOM VARIABLES 1 - 3 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand1, rand2, rand3, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand1 + rand2 + rand3 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_1-3.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2965
# [1] 3126
# [1] 6486

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand1", "rand2", "rand3")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand1 rand2 rand3 
# 3684  2669  6224 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand1-3.csv")

# RANDOM VARIABLES 4 - 6 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand4, rand5, rand6, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand4 + rand5 + rand6 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_4-6.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1844
# [1] 3247
# [1] 1808

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand4", "rand5", "rand6")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand4 rand5 rand6 
# 2550  2293  2056 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand4-6.csv")

# RANDOM VARIABLES 7 - 9 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand7, rand8, rand9, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand7 + rand8 + rand9 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_7-9.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1523
# [1] 3789
# [1] 2411

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand7", "rand8", "rand9")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand7 rand8 rand9 
# 3825  1471  2427 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand7-9.csv")

# RANDOM VARIABLES 10 - 12 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand10, rand11, rand12, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand10 + rand11 + rand12 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_10-12.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1619
# [1] 2959
# [1] 2308

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand10", "rand11", "rand12")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand10 rand11 rand12 
# 2255   2967   1664 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand10-12.csv")

# RANDOM VARIABLES 13 - 15 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand13, rand14, rand15, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand13 + rand14 + rand15 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_13-15.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1161
# [1] 2453
# [1] 5391

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand13", "rand14", "rand15")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand13 rand14 rand15 
# 5168   1352   2485 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand13-15.csv")

# RANDOM VARIABLES 16 - 18 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand16, rand17, rand18, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand16 + rand17 + rand18 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_16-18.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2965
# [1] 3126
# [1] 6486

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand16", "rand17", "rand18")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand16 rand17 rand18 
# 3684   2669   6224

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand16-18.csv")

# RANDOM VARIABLES 19 - 21 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand19, rand20, rand21, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand19 + rand20 + rand21 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_19-21.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1219
# [1] 2237
# [1] 6520

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand19", "rand20", "rand21")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand19 rand20 rand21 
# 4202   1280   4494 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand19-21.csv")

# RANDOM VARIABLES 22 - 24 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand22, rand23, rand24, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand22 + rand23 + rand24 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_22-24.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1232
# [1] 1696
# [1] 3714

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand22", "rand23", "rand24")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand22 rand23 rand24 
# 1616   3744   1282 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand22-24.csv")

# RANDOM VARIABLES 25 - 27 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand25, rand26, rand27, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand25 + rand26 + rand27 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_25-27.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2481
# [1] 3195
# [1] 2830

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand25", "rand26", "rand27")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand25 rand26 rand27 
# 3124   2978   2404 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand25-27.csv")

# RANDOM VARIABLES 28 - 30 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand28, rand29, rand30, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand28 + rand29 + rand30 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_28-30.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1784
# [1] 1408
# [1] 3165

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand28", "rand29", "rand30")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand28 rand29 rand30 
# 1433   2452   2472 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand28-30.csv")

# RANDOM VARIABLES 31 - 33 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand31, rand32, rand33, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand31 + rand32 + rand33 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_31-33.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 3934
# [1] 4789
# [1] 2992

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand31", "rand32", "rand33")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand31 rand32 rand33 
# 3435   4659   3621 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand31-33.csv")

# RANDOM VARIABLES 34 - 36 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand34, rand35, rand36, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand34 + rand35 + rand36 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_34-36.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2867
# [1] 3292
# [1] 3941

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand34", "rand35", "rand36")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand34 rand35 rand36 
# 3198   3336   3566 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand34-36.csv")

# RANDOM VARIABLES 37 - 39 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand37, rand38, rand39, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand37 + rand38 + rand39 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_37-39.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1133
# [1] 2306
# [1] 4003

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand37", "rand38", "rand39")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand37 rand38 rand39 
# 3882   1191   2369 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand37-39.csv")

# RANDOM VARIABLES 40 - 42 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand40, rand41, rand42, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand40 + rand41 + rand42 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_40-42.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1000
# [1] 2070
# [1] 2650

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand40", "rand41", "rand42")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand40 rand41 rand42 
# 2540   1093   2087 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand40-42.csv")

# RANDOM VARIABLES 43 - 45 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand43, rand44, rand45, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand43 + rand44 + rand45 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_43-45.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2973
# [1] 3371
# [1] 3928

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand43", "rand44", "rand45")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand43 rand44 rand45 
# 4171   3340   2761 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand43-45.csv")

# RANDOM VARIABLES 46 - 48 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand46, rand47, rand48, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand46 + rand47 + rand48 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_46-48.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 928
# [1] 1867
# [1] 3184

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand46", "rand47", "rand48")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand46 rand47 rand48 
# 1331   3343   1305 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand46-48.csv")

# RANDOM VARIABLES 49 - 51 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand49, rand50, rand51, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand49 + rand50 + rand51 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_49-51.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 3931
# [1] 3597
# [1] 5675

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand49", "rand50", "rand51")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand49 rand50 rand51 
# 5246   4317   3640 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand49-51.csv")

# RANDOM VARIABLES 52 - 54 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand52, rand53, rand54, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand52 + rand53 + rand54 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_52-54.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1886
# [1] 3452
# [1] 3310

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand52", "rand53", "rand54")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand52 rand53 rand54 
# 3095   2601   2952 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand52-54.csv")

# RANDOM VARIABLES 55 - 57 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand55, rand56, rand57, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand55 + rand56 + rand57 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_55-57.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2353
# [1] 3643
# [1] 5818

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand55", "rand56", "rand57")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand55 rand56 rand57 
# 3616   5893   2305

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand55-57.csv")

# RANDOM VARIABLES 58 - 60 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand58, rand59, rand60, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand58 + rand59 + rand60 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_58-60.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1564
# [1] 3885
# [1] 6222

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand58", "rand59", "rand60")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand58 rand59 rand60 
# 2896   4616   4159 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand58-60.csv")

# RANDOM VARIABLES 61 - 63 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand61, rand62, rand63, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand61 + rand62 + rand63 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_61-63.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 4673
# [1] 3633
# [1] 3257

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand61", "rand62", "rand63")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand61 rand62 rand63 
# 4752   2601   4210 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand61-63.csv")

# RANDOM VARIABLES 64 - 66 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand64, rand65, rand66, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand64 + rand65 + rand66 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_64-66.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1399
# [1] 3089
# [1] 2679

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand64", "rand65", "rand66")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand64 rand65 rand66 
# 2498   2073   2596 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand64-66.csv")

# RANDOM VARIABLES 67 - 69 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand67, rand68, rand69, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand67 + rand68 + rand69 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_67-69.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1915
# [1] 3237
# [1] 4017

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand67", "rand68", "rand69")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand67 rand68 rand69 
# 3228   3001   2940 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand67-69.csv")

# RANDOM VARIABLES 70 - 72 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand70, rand71, rand72, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand70 + rand71 + rand72 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_70-72.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2929
# [1] 2210
# [1] 1873

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand70", "rand71", "rand72")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand70 rand71 rand72 
# 1932   2092   2988 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand70-72.csv")

# RANDOM VARIABLES 73 - 75 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand73, rand74, rand75, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand73 + rand74 + rand75 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_73-75.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1195
# [1] 3667
# [1] 3029

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand73", "rand74", "rand75")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand73 rand74 rand75 
# 3242   2248   2401 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand73-75.csv")

# RANDOM VARIABLES 76 - 78 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand76, rand77, rand78, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand76 + rand77 + rand78 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_76-78.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1643
# [1] 4340
# [1] 3873

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand76", "rand77", "rand78")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand76 rand77 rand78 
# 4268   1684   3904 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand76-78.csv")

# RANDOM VARIABLES 79 - 81 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand79, rand80, rand81, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand79 + rand80 + rand81 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_79-81.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1961
# [1] 2333
# [1] 2785

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand79", "rand80", "rand81")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand79 rand80 rand81 
# 2750   1858   2471 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand79-81.csv")

# RANDOM VARIABLES 82 - 84 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand82, rand83, rand84, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand82 + rand83 + rand84 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_82-84.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 1952
# [1] 4535
# [1] 4327

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand82", "rand83", "rand84")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand82 rand83 rand84 
# 4153   3646   3015 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand82-84.csv")

# RANDOM VARIABLES 85 - 87 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand85, rand86, rand87, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand85 + rand86 + rand87 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_85-87.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2460
# [1] 3695
# [1] 7164

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand85", "rand86", "rand87")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand85 rand86 rand87 
# 7161   2460   3698 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand85-87.csv")

# RANDOM VARIABLES 88 - 90 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand88, rand89, rand90, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand88 + rand89 + rand90 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_88-90.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 767
# [1] 2132
# [1] 4749

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand88", "rand89", "rand90")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand88 rand89 rand90 
# 3066   1805   2777 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand88-90.csv")

# RANDOM VARIABLES 91 - 93 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand91, rand92, rand93, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand91 + rand92 + rand93 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_91-93.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 4271
# [1] 3139
# [1] 2932

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand91", "rand92", "rand93")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand91 rand92 rand93 
# 4366   2771   3205 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand91-93.csv")

# RANDOM VARIABLES 94 - 96 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand94, rand95, rand96, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand94 + rand95 + rand96 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_94-96.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2632
# [1] 3324
# [1] 2775

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand94", "rand95", "rand96")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand94 rand95 rand96 
# 2918   2699   3114 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand94-96.csv")

# RANDOM VARIABLES 97 - 99 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand97, rand98, rand99, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand97 + rand98 + rand99 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_97-99.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 2026
# [1] 3318
# [1] 4526

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand97", "rand98", "rand99")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand97 rand98 rand99 
# 2280   3346   4244 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand97-99.csv")

# RANDOM VARIABLES 98 - 100 ###################################################################################################################################################################################################################################

# Subset 3 random variables + fishing data 

retained_preds <- subset(envs, select=c(rand98, rand99, rand100, fishing_ind_sum))

# pRDA_constrained 3 rands + fishing (ind sum) as condition

prda_1 <- rda(snps_unlinked_imp_5x ~ rand98 + rand99 + rand100 + Condition(fishing_ind_sum) , data=retained_preds, scale=T)

vif.cca(prda_1)

list.save(prda_1, "pRDA_rand_5x_n24_prda_98-100.rds")

# subset env variables to exclude fishing data
envs_retained <- subset(retained_preds, select=-c(fishing_ind_sum))

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per RDA axis using z = 3 standard deviations as a cut off
# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)

# [1] 718
# [1] 1986
# [1] 3316

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)

# organize  results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with each environmental predictors
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors

colnames(foo) <- c("rand98", "rand99", "rand100")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where n is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)

# rand100  rand98  rand99 
# 2366    1864    1790 

write.csv(cand,"pRDA_rand_5x_n24_prda_putative_adaptive_SNPs_rand98-100.csv")
