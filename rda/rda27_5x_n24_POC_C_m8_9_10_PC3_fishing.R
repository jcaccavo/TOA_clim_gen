# helpful links
# https://github.com/Tom-Jenkins/seascape_rda_tutorial/blob/master/4.Redundancy_analysis/4.redundancy_analysis.R
# https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples

# run on server (yangtze)
# install necessary packages using conda
# conda install -c conda-forge r-tidyverse
# conda install -c conda-forge r-psych
# conda install -c conda-forge r-adespatial
# conda install -c conda-forge r-vegan
# conda install -c conda-forge r-rlist
# conda install -c conda-forge r-adegenet

# Run in R in directory

# in R
# Load packages
library(tidyverse)
library(psych)
library(adespatial)
library(vegan)
library(rlist)
library(adegenet)

# Load SNPs
snps_unlinked_5x <-read.PLINK(".../5x_TOA_only_filtered_SNPs_all_unlinked.raw")
# Reading PLINK raw format into a genlight object... 
# 
# 
# Reading loci information... 
# 
# Reading and converting genotypes... 
# .
# Building final object... 
# 
# ...done.

dim(snps_unlinked_5x)
# [1]      24 2359009

sum(is.na(snps_unlinked_5x))
# [1] 0
# Warning message:
#   In is.na(snps_unlinked_5x) :
#   is.na() applied to non-(list or vector) of type 'S4'

# Impute missing SNPs
snps_unlinked_imp_5x <- apply(snps_unlinked_5x, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
# may take a few minutes

dim(snps_unlinked_imp_5x)
# [1]      24 2359009

sum(is.na(snps_unlinked_imp_5x))
# [1] 0

#FORWARD SELECTION#########################################################################################################

# Load predictor variables
envs <- read.csv(".../rda27_5x_n24_POC_C_m8_9_10_PC3_fishing.csv")

str(envs)
# 'data.frame':   24 obs. of  4 variables:
# $ PC1            : num  -1.924 0.558 1.648 1.648 1.648 ...
# $ PC2            : num  -3.69 -1.233 0.921 0.921 0.921 ...
# $ PC3            : num  2.477 -1.315 0.957 0.957 0.957 ...
# $ fishing_ind_sum: int  398 0 0 0 0 0 398 398 398 356 ...

#Forward selection environmental variables
sel_fs_envs <- forward.sel(Y = snps_unlinked_imp_5x, X = envs)
# Testing variable 1
# Testing variable 2
# Testing variable 3
# Procedure stopped (alpha criteria): pvalue for variable 3 is 0.301000 (> 0.050000)

sel_fs_envs
# variables order         R2      R2Cum    AdjR2Cum        F pvalue
# 1       PC1     1 0.04856934 0.04856934 0.005322487 1.123072  0.009
# 2       PC3     3 0.04824548 0.09681481 0.010797177 1.121758  0.019

pval.adj <- p.adjust(sel_fs_envs$pval, method = 'holm', n = 18)
sel_fs_envs$pval.adj <- pval.adj
sel_fs_envs
# variables order         R2      R2Cum    AdjR2Cum        F pvalue pval.adj
# 1       PC1     1 0.04856934 0.04856934 0.005322487 1.123072  0.009    0.162
# 2       PC3     3 0.04824548 0.09681481 0.010797177 1.121758  0.019    0.323

#COLLINNEARITY#########################################################################################################

# check the level of correlation between predictors (should be |r| < 0.7)
pdf("rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_pair_panel_collinearity_plot_envs.pdf")
pairs.panels(envs, scale=T)
dev.off()

# no cross-correlation between PCs or fishing data whatsover

#CONSTRAINED RDAs#########################################################################################################

#pRDA_constrained 1 all PCs (n = 3) + fishing (ind sum) as condition
prda_1 <- rda(snps_unlinked_imp_5x ~ PC1 + PC2 + PC3 + Condition(fishing_ind_sum) , data=envs, scale=T)
# takes <1 minute
prda_1
# Call: rda(formula = snps_unlinked_imp_5x ~ PC1 + PC2 + PC3 +
#             Condition(fishing_ind_sum), data = envs, scale = T)
# 
# -- Model Summary --
#   
#   Inertia Proportion Rank
# Total         2.359e+06  1.000e+00     
# Conditional   1.077e+05  4.565e-02    1
# Constrained   3.367e+05  1.427e-01    3
# Unconstrained 1.915e+06  8.116e-01   19
# 
# Inertia is correlations
# 
# -- Eigenvalues --
#   
#   Eigenvalues for constrained axes:
#   RDA1   RDA2   RDA3 
# 137194 101415  98098 
# 
# Eigenvalues for unconstrained axes:
#   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
# 141334 110926 105742 104500 103132 102498 100587  99825 
# (Showing 8 of 19 unconstrained eigenvalues)

# VIF - variance inflation factor
vif.cca(prda_1)
# fishing_ind_sum     PC1             PC2             PC3 
# 1.055950        1.054025        1.001227        1.000698 

# VIF all <3! will proceed 

summary(eigenvals(prda_1, model = "constrained"))
# Importance of components:
#                         RDA1      RDA2      RDA3
# Eigenvalue            1.372e+05 1.014e+05 9.810e+04
# Proportion Explained  4.075e-01 3.012e-01 2.913e-01
# Cumulative Proportion 4.075e-01 7.087e-01 1.000e+00

RsquareAdj(prda_1)
# $r.squared
# [1] 0.1427322
# 
# $adj.r.squared
# [1] 0.01524518

list.save(prda_1, "rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_prda_1.rds")

anova.cca(prda_1, permutations = 1000, by="terms")
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# Model: rda(formula = snps_unlinked_imp_5x ~ PC1 + PC2 + PC3 + Condition(fishing_ind_sum), data = envs, scale = T)
#           Df Variance      F  Pr(>F)
# PC1       1   116493 1.1560 0.02198 *
# PC2       1   104246 1.0345 0.30969
# PC3       1   115967 1.1508 0.03397 *
#   Residual 19  1914609
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#PLOT RDAs#########################################################################################################

#Set colors for plotting
# all points by CCAMLR areas 48, 58, 88
bg_pop1_all <- c("#00441b",
                 "#762a83",
                 "#762a83",
                 "#762a83",
                 "#762a83",
                 "#762a83",
                 "#00441b",
                 "#00441b",
                 "#00441b",
                 "#a6dba0",
                 "#a6dba0",
                 "#00441b",
                 "#a6dba0",
                 "#00441b",
                 "#00441b",
                 "#762a83",
                 "#762a83",
                 "#a6dba0",
                 "#00441b",
                 "#00441b",
                 "#00441b",
                 "#00441b",
                 "#00441b",
                 "#00441b") 
# population categories for legend
bg_pop1 <- c("#762a83","#a6dba0","#00441b")

# all points by cohort categories pre- and post-2000
bg_cohort2_all <- c("#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#980043",
                    "#d7b5d8",
                    "#d7b5d8",
                    "#d7b5d8",
                    "#d7b5d8",
                    "#d7b5d8",
                    "#d7b5d8",
                    "#d7b5d8",
                    "#d7b5d8") 
# population categories for legend
bg_cohort2 <- c("#d7b5d8","#980043")

# plot prda_1 by pop1 
pdf("rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_prda_1_pop1.pdf")
plot(prda_1, type="n", scaling=3, xlab="RDA1 40.8%", ylab="RDA2 30.1%")
points(prda_1, display="sites", pch=21, cex=1, scaling=3, bg=bg_pop1_all)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("48","58","88"), bty="o", col="gray32", pch=21, cex=1, pt.bg=bg_pop1)
dev.off()

# plot prda_1 by cohort2 
pdf("rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_prda_1_cohort2.pdf")
plot(prda_1, type="n", scaling=3, xlab="RDA1 40.8%", ylab="RDA2 30.1%")
points(prda_1, display="sites", pch=21, cex=1, scaling=3, bg=bg_cohort2_all)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("pre2000","post2000"), bty="o", col="gray32", pch=21, cex=1, pt.bg=bg_cohort2)
dev.off()

# screeplot of prda_1
pdf("rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_prda_1_screeplot.pdf")
screeplot(prda_1)
dev.off()

# CANDIDATE SNPS#########################################################################################################

# useful links
# https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html#conclusions

# Load predictor variables
envs <- read.csv(".../rda27_5x_n24_POC_C_m8_9_10_PC3_fishing.csv")

str(envs)
# 'data.frame':   24 obs. of  4 variables:
# $ PC1            : num  -1.924 0.558 1.648 1.648 1.648 ...
# $ PC2            : num  -3.69 -1.233 0.921 0.921 0.921 ...
# $ PC3            : num  2.477 -1.315 0.957 0.957 0.957 ...
# $ fishing_ind_sum: int  398 0 0 0 0 0 398 398 398 356 ...

# subset predictor variables to only include variables identified as significant from the ANOVA
# for POC-C, this was PC1 and PC3
envs_retained <- subset(envs, select=c(PC1,PC3))

str(envs_retained)
# 'data.frame':   24 obs. of  2 variables:
# $ PC1: num  -1.924 0.558 1.648 1.648 1.648 ...
# $ PC3: num  2.477 -1.315 0.957 0.957 0.957 ...

# if continuing directly from RDAs, skip loading the RDA and SNPs

# ***********************************************************************************************************************

# if coming back to an old analysis in a new R workspace, you must load the RDA, as well as the SNPs, and impute them

# Load selected RDA
prda_1 <- read_rds("rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_prda_1.rds")

# Load SNPs
snps_unlinked_5x <-read.PLINK(".../5x_TOA_only_filtered_SNPs_all_unlinked.raw")
# Reading PLINK raw format into a genlight object... 
# 
# 
# Reading loci information... 
# 
# Reading and converting genotypes... 
# .
# Building final object... 
# 
# ...done.

# Impute missing SNPs
snps_unlinked_imp_5x <- apply(snps_unlinked_5x, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
# may take a few minutes
length(snps_unlinked_imp_5x)
# [1] 56,616,216

# ***********************************************************************************************************************

# Extract SNP loadings for all axes (RDA 1 - 3, PC 1- 19)
load.rda <- scores(prda_1, choices=c(1:22), display="species")

# Plot histograms of SNP loadings from RDA axes
pdf("rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_prda_1_load_constrained_1.pdf")
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")
dev.off()

# Identify SNPs in the tails of the distribution
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Candidate SNPs per axis using z = 3 standard deviations as a cut off

# RDA 1
cand1 <- outliers(load.rda[,1],3)
length(cand1)
# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)
# PC 1
cand4 <- outliers(load.rda[,4],3)
length(cand4)
# PC 2
cand5 <- outliers(load.rda[,5],3)
length(cand5)
# PC 3
cand6 <- outliers(load.rda[,6],3)
length(cand6)
# PC 4
cand7 <- outliers(load.rda[,7],3)
length(cand7)
# PC 5
cand8 <- outliers(load.rda[,8],3)
length(cand8)
# PC 6
cand9 <- outliers(load.rda[,9],3)
length(cand9)
# PC 7
cand10 <- outliers(load.rda[,10],3)
length(cand10)
# PC 8
cand11 <- outliers(load.rda[,11],3)
length(cand11)
# PC 9
cand12 <- outliers(load.rda[,12],3)
length(cand12)
# PC 10
cand13 <- outliers(load.rda[,13],3)
length(cand13)
# PC 11
cand14 <- outliers(load.rda[,14],3)
length(cand14)
# PC 12
cand15 <- outliers(load.rda[,15],3)
length(cand15)
# PC 13
cand16 <- outliers(load.rda[,16],3)
length(cand16)
# PC 14
cand17 <- outliers(load.rda[,17],3)
length(cand17)
# PC 15
cand18 <- outliers(load.rda[,18],3)
length(cand18)
# PC 16
cand19 <- outliers(load.rda[,19],3)
length(cand19)
# PC 17
cand20 <- outliers(load.rda[,20],3)
length(cand20)
# PC 18
cand21 <- outliers(load.rda[,21],3)
length(cand21)
# PC 19
cand22 <- outliers(load.rda[,22],3)
length(cand22)

# all length() output values
# [1] 648
# [1] 5591
# [1] 4445
# [1] 2006
# [1] 5483
# [1] 8673
# [1] 10253
# [1] 2569
# [1] 3875
# [1] 1893
# [1] 5388
# [1] 1853
# [1] 7502
# [1] 8438
# [1] 6488
# [1] 2570
# [1] 2777
# [1] 3945
# [1] 6366
# [1] 13157
# [1] 4312
# [1] 5070

# Candidate total from both RDA (3) and PC (20) axes
ncand_all <- length(cand1) + length(cand2) + length(cand3) + length(cand4) + length(cand5) + length(cand6) + length(cand7) + length(cand8) + length(cand9) + length(cand10) + length(cand11) + length(cand12) + length(cand13) + length(cand14) + length(cand15) + length(cand16) + length(cand17) + length(cand18) + length(cand19) + length(cand20) + length(cand21) + length(cand22)
ncand_all
# [1] 113302

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)
ncand_rda3
# [1] 10684

# Candidate total from RDA axes 1 and 2 only
ncand_rda2 <- length(cand1) + length(cand2)
ncand_rda2
# [1] 6239

# Candidate total from RDA axis 1 only
ncand_rda1 <- length(cand1)
ncand_rda1
# [1] 648

# Candidate SNP names from RDA axes 1 and 2
prda1cand2 <- c(names(cand1), names(cand2))

# Candidate SNP names from RDA axes 1 and 2
prda1cand3 <- c(names(cand1), names(cand2), names(cand3))

# Candidate SNPs detected on multiple RDA axes
length(prda1cand3[duplicated(prda1cand3)])
# [1] 0
# No duplicates

# Candidate SNPs detected on RDA axes 1 and 2
length(prda1cand2[duplicated(prda1cand2)])
# [1] 0
# No duplicates

# In the event that there are duplicates, use the following to remove
# prda1cand3 <- prda1cand3[!duplicated(prda1cand3)] 
# prda1cand2 <- prda1cand2[!duplicated(prda1cand2)] 

# Set up the color scheme for plotting:
bgcol  <- ifelse(colnames(snps_unlinked_imp_5x) %in% prda1cand2, 'gray32', '#00000000')
snpcol <- ifelse(colnames(snps_unlinked_imp_5x) %in% prda1cand2, 'red', '#00000000')

## axes 1 & 2 - zooming in to just the SNPs here...
pdf("rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_prda_1_SNP_zoom_axes1_2.pdf")
plot(prda_1, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="rda24.prda1, axes 1 and 2", xlab="RDA1 40.8%", ylab="RDA2 30.1%")
points(prda_1, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(prda_1, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
dev.off()

pdf("rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_prda_1_SNP_zoom_axes2_3.pdf")
plot(prda_1, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="rda24.prda1, axes 2 and 3", xlab="RDA1 40.8%", ylab="RDA2 30.1%")
points(prda_1, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(prda_1, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
dev.off()

# Let’s see which environmental predictors are most strongly correlated with the first three RDA axes
intersetcor(prda_1)[,1:3]
#         RDA1        RDA2       RDA3
# PC1 -0.6294475  0.69423465 -0.2954998
# PC2 -0.3758424  0.01096054  0.9198150
# PC3 -0.6132545 -0.72210922 -0.2637790

# candidate SNPs on RDA1 represent multilocus sets of SNP genotypes associated most strongly with PC1
# candidate SNPs on RDA2 represent multilocus sets of SNP genotypes associated most strongly with PC2
# candidate SNPs on RDA3 represent multilocus sets of SNP genotypes associated most strongly with PC3

# organize results into one data frame with the axis, SNP name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with predictor (n = 2) PC1 & PC3
foo <- matrix(nrow=(ncand_rda3), ncol=2)  # 2 columns for 2 predictors
colnames(foo) <- c("PC1","PC3")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  
head(cand)
#    axis                     snp    loading        PC1        PC3
# 1    1   HiC_scaffold_10:554904_G 0.03525069 -0.7382876 -0.2712823
# 2    1  HiC_scaffold_10:6538681_T 0.03586705 -0.5586857 -0.3205270
# 3    1  HiC_scaffold_10:6539276_G 0.03564588 -0.2509156 -0.5940660
# 4    1  HiC_scaffold_10:9856833_C 0.03561058 -0.3904337 -0.4243519
# 5    1 HiC_scaffold_10:10078005_A 0.03586705 -0.5586857 -0.3205270
# 6    1 HiC_scaffold_10:10221119_T 0.03614859 -0.2986261 -0.5885249

length(cand$snp[duplicated(cand$snp)])
# takes <1 minute
# [1] 0

foo <- cbind(cand$axis, duplicated(cand$snp))
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
# 0 
# 648
table(foo[foo[,1]==2,2]) # no duplicates on axis 2
# 0 
# 5591
table(foo[foo[,1]==3,2]) # no duplicates on axis 3
# 0 
# 4445

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,6] <- names(which.max(abs(bar[4:5]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where k is the new column you will add, so k = n + 1
  cand[i,7] <- max(abs(bar[4:5]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[6] <- "predictor"
colnames(cand)[7] <- "correlation"

table(cand$predictor)
# PC1  PC3 
# 7237 3447 

# this is the total candidate SNPs explained all RDA axes

write.csv(cand,"rda27_5x_n24_POC_C_m8_9_10_PC3_fishing_prda_1_putative_adaptive_SNPs.csv")
