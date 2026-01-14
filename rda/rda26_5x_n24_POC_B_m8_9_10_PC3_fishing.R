# helpful links
# https://github.com/Tom-Jenkins/seascape_rda_tutorial/blob/master/4.Redundancy_analysis/4.redundancy_analysis.R
# https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/rda26_5x_n24_POC_B_m8_9_10_PC3_fishing")

# Load packages
library(tidyverse)
library(psych)
library(adespatial)
library(vegan)
library(rlist)
library(adegenet)

# Load SNPs
snps_unlinked_5x <-read.PLINK("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/genomic_data/5x_TOA_only_filtered_SNPs_all_unlinked.raw")
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
envs <- read.csv("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/rda26_5x_n24_POC_B_m8_9_10_PC3_fishing/rda26_5x_n24_POC_B_m8_9_10_PC3_fishing.csv")
# if running on server
envs <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/rda26_5x_n24_POC_B_m8_9_10_PC3_fishing.csv")

str(envs)
# 'data.frame':	24 obs. of  4 variables:
# $ PC1            : num  1.294 2.173 0.676 0.676 0.676 ...
# $ PC2            : num  -3.114 -0.269 0.912 0.912 0.912 ...
# $ PC3            : num  -1.843 0.379 -0.211 -0.211 -0.211 ...
# $ fishing_ind_sum: int  398 0 0 0 0 0 398 398 398 356 ...

#Forward selection environmental variables
sel_fs_envs <- forward.sel(Y = snps_unlinked_imp_5x, X = envs)
# Testing variable 1                                                                                                                                                                                                                                        
# Procedure stopped (alpha criteria): pvalue for variable 1 is 0.060000 (> 0.050000)                                                                                                                                                                        
# Error in forward.sel(Y = snps_unlinked_imp_5x, X = envs) :                                                                                                                                                                                                  
#   No variables selected. Please change your parameters.

sel_fs_envs
# Error: object 'sel_fs_envs' not found

#COLLINNEARITY#########################################################################################################

# check the level of correlation between predictors (should be |r| < 0.7)
pdf("rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_pair_panel_collinearity_plot_envs.pdf")
pairs.panels(envs, scale=T)
dev.off()

# no cross-correlation between PCs or fishing data whatsover

#CONSTRAINED RDAs#########################################################################################################

#pRDA_constrained 1 all PCs (n = 3) + fishing (ind sum) as condition
prda_1 <- rda(snps_unlinked_imp_5x ~ PC1 + PC2 + PC3 + Condition(fishing_ind_sum) , data=envs, scale=T)
# takes <1 minute
prda_1
# Call: rda(formula = snps_unlinked_imp_5x ~ PC1 + PC2 + PC3 + Condition(fishing_ind_sum), data = envs, scale = T)
# 
# Inertia Proportion Rank
# Total         2.359e+06  1.000e+00     
# Conditional   1.077e+05  4.565e-02    1
# Constrained   3.074e+05  1.303e-01    3
# Unconstrained 1.944e+06  8.240e-01   19
# Inertia is correlations 
# 
# Eigenvalues for constrained axes:
#   RDA1   RDA2   RDA3 
# 111456  98564  97429 
# 
# Eigenvalues for unconstrained axes:
#   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
# 166781 110796 105937 104552 104060 102535 100334  99773 
# (Showing 8 of 19 unconstrained eigenvalues)

# VIF - variance inflation factor
vif.cca(prda_1)
# fishing_ind_sum     PC1             PC2             PC3 
# 1.021108        1.007904        1.004509        1.008695 

# VIF all <3! will proceed 

summary(eigenvals(prda_1, model = "constrained"))
# Importance of components:
#                         RDA1      RDA2      RDA3
# Eigenvalue            1.115e+05 9.856e+04 9.743e+04
# Proportion Explained  3.625e-01 3.206e-01 3.169e-01
# Cumulative Proportion 3.625e-01 6.831e-01 1.000e+00

RsquareAdj(prda_1)
# $r.squared
# [1] 0.1303299
# 
# $adj.r.squared
# [1] 0.0002318511

list.save(prda_1, "rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_prda_1.rds")

anova.cca(prda_1, permutations = 1000, by="terms")
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# Model: rda(formula = snps_unlinked_imp_5x ~ PC1 + PC2 + PC3 + Condition(fishing_ind_sum), data = envs, scale = T)
#           Df Variance      F Pr(>F)
# PC1       1   102865 1.0054 0.5854
# PC2       1   105727 1.0334 0.3397
# PC3       1    98857 0.9663 0.7473
# Residual 19  1943866

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
pdf("rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_prda_1_pop1.pdf")
plot(prda_1, type="n", scaling=3, xlab="RDA1 36.3%", ylab="RDA2 32.1%")
points(prda_1, display="sites", pch=21, cex=1, scaling=3, bg=bg_pop1_all)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("48","58","88"), bty="o", col="gray32", pch=21, cex=1, pt.bg=bg_pop1)
dev.off()

# plot prda_1 by cohort2 
pdf("rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_prda_1_cohort2.pdf")
plot(prda_1, type="n", scaling=3, xlab="RDA1 36.3%", ylab="RDA2 32.1%")
points(prda_1, display="sites", pch=21, cex=1, scaling=3, bg=bg_cohort2_all)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("pre2000","post2000"), bty="o", col="gray32", pch=21, cex=1, pt.bg=bg_cohort2)
dev.off()

# much less clustering in the POC-B dataset with respect to the SG-B

# screeplot of prda_1
pdf("rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_prda_1_screeplot.pdf")
screeplot(prda_1)
dev.off()

# CANDIDATE SNPS#########################################################################################################

# useful links
# https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html#conclusions

# Load predictor variables
envs <- read.csv("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/rda26_5x_n24_POC_B_m8_9_10_PC3_fishing/rda26_5x_n24_POC_B_m8_9_10_PC3_fishing.csv")

str(envs)
# 'data.frame':	24 obs. of  4 variables:
# $ PC1            : num  1.294 2.173 0.676 0.676 0.676 ...
# $ PC2            : num  -3.114 -0.269 0.912 0.912 0.912 ...
# $ PC3            : num  -1.843 0.379 -0.211 -0.211 -0.211 ...
# $ fishing_ind_sum: int  398 0 0 0 0 0 398 398 398 356 ...

# subset predictor variables to only include variables identified as significant from the ANOVA
# given that none were identified as significant, moving forward with all 3 PCs
envs_retained <- subset(envs, select=c(PC1, PC2, PC3))

str(envs_retained)
# 'data.frame':   24 obs. of  3 variables:
# $ PC1: num  1.294 2.173 0.676 0.676 0.676 ...
# $ PC2: num  -3.114 -0.269 0.912 0.912 0.912 ...
# $ PC3: num  -1.843 0.379 -0.211 -0.211 -0.211 ...

# if continuing directly from RDAs, skip loading the RDA and SNPs

# ***********************************************************************************************************************

# if coming back to an old analysis in a new R workspace, you must load the RDA, as well as the SNPs, and impute them

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/rda26_5x_n24_POC_B_m8_9_10_PC3_fishing")

# Load selected RDA
prda_1 <- read_rds("rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_prda_1.rds")

# Load SNPs
snps_unlinked_5x <-read.PLINK("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/genomic_data/5x_TOA_only_filtered_SNPs_all_unlinked.raw")
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
pdf("rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_prda_1_load_constrained_1.pdf")
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
# [1] 6091
# [1] 3701
# [1] 5150
# [1] 15880
# [1] 6426
# [1] 3928
# [1] 9391
# [1] 5342
# [1] 2341
# [1] 7840
# [1] 7030
# [1] 3352
# [1] 9182
# [1] 9956
# [1] 3492
# [1] 2322
# [1] 3345
# [1] 2924
# [1] 9319
# [1] 7325
# [1] 12855
# [1] 3859

# Candidate total from both RDA (3) and PC (20) axes
ncand_all <- length(cand1) + length(cand2) + length(cand3) + length(cand4) + length(cand5) + length(cand6) + length(cand7) + length(cand8) + length(cand9) + length(cand10) + length(cand11) + length(cand12) + length(cand13) + length(cand14) + length(cand15) + length(cand16) + length(cand17) + length(cand18) + length(cand19) + length(cand20) + length(cand21) + length(cand22)
ncand_all
# [1] 141051

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)
ncand_rda3
# [1] 14942

# Candidate total from RDA axes 1 and 2 only
ncand_rda2 <- length(cand1) + length(cand2)
ncand_rda2
# [1] 9792

# Candidate total from RDA axis 1 only
ncand_rda1 <- length(cand1)
ncand_rda1
# [1] 6091

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
pdf("rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_prda_1_SNP_zoom_axes1_2.pdf")
plot(prda_1, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="rda24.prda1, axes 1 and 2", xlab="RDA1 36.3%", ylab="RDA2 32.1%")
points(prda_1, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(prda_1, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
dev.off()

pdf("rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_prda_1_SNP_zoom_axes2_3.pdf")
plot(prda_1, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="rda24.prda1, axes 2 and 3", xlab="RDA1 36.3%", ylab="RDA2 32.1%")
points(prda_1, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(prda_1, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
dev.off()

# Let’s see which environmental predictors are most strongly correlated with the first three RDA axes
intersetcor(prda_1)[,1:3]
#         RDA1        RDA2       RDA3
# PC1  0.6027838 -0.02623656  0.7809446
# PC2  0.7387110  0.22202354 -0.6061996
# PC3 -0.1433424  0.97477924  0.1607781

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

# add in the correlations of each candidate SNP with each predictor (n = 3)
# given that none were significant, working with all 3
foo <- matrix(nrow=(ncand_rda3), ncol=3)  # 3 columns for 3 predictors
colnames(foo) <- c("PC1","PC2","PC3")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  
head(cand)
#     axis                    snp     loading        PC1        PC2         PC3
# 1    1  HiC_scaffold_10:204373_T -0.03871754 -0.8880169 -0.2120474 -0.04737198
# 2    1  HiC_scaffold_10:467004_C -0.03385615 -0.4544402 -0.4225881  0.08496964
# 3    1  HiC_scaffold_10:554904_G -0.03345369 -0.5557538 -0.3427917  0.04225816
# 4    1  HiC_scaffold_10:893491_C -0.03948807 -0.3433457 -0.7052147 -0.26405203
# 5    1 HiC_scaffold_10:1209683_G -0.03652717 -0.5124733 -0.4343583  0.12076904
# 6    1 HiC_scaffold_10:1210310_T -0.03506982 -0.5138328 -0.3317365  0.03501947

length(cand$snp[duplicated(cand$snp)])
# takes <1 minute
# [1] 0

foo <- cbind(cand$axis, duplicated(cand$snp))
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
# 0 
# 6091
table(foo[foo[,1]==2,2]) # no duplicates on axis 2
# 0 
# 3701
table(foo[foo[,1]==3,2]) # no duplicates on axis 3
# 0 
# 5150

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,7] <- names(which.max(abs(bar[4:6]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where k is the new column you will add, so k = n + 1
  cand[i,8] <- max(abs(bar[4:6]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[7] <- "predictor"
colnames(cand)[8] <- "correlation"

table(cand$predictor)
# PC1  PC2  PC3 
# 7057 4325 3560

# this is the total candidate SNPs explained all RDA axes

write.csv(cand,"rda26_5x_n24_POC_B_m8_9_10_PC3_fishing_prda_1_putative_adaptive_SNPs.csv")

# ****UP TO HERE**** ####
# on server
# jc_screen_1

# need to compare this candidate SNP list (n = 14,942 candidates) with the LFMM list for POC-B, to whittle down to the final overlapping candidate list
# then can proceed with adptively enriched pRDA based on subset vcf using only overlapping candidates from pRDA and LFMM run on POC-B environmental data

# Adaptively enriched pRDA #########################################################################################################

# https://www.rdocumentation.org/packages/adegenet/versions/2.1.10/topics/extract.PLINKmap

# convert vcf file to bed on server using plink - run in prompt
# plink2 --vcf 5x_TOA_only_filtered_SNPs_unlinked_candidate_SNPs_SG_B_854.vcf --make-bed --allow-extra-chr --out 5x_TOA_only_filtered_SNPs_unlinked_candidate_SNPs_SG_B_854

# convert vcf file to raw on server using plink - run in prompt
# plink --vcf 5x_TOA_only_filtered_SNPs_unlinked_candidate_SNPs_SG_B_854.vcf --chr-set 95 --allow-extra-chr --recode A --out 5x_TOA_only_filtered_SNPs_unlinked_candidate_SNPs_SG_B_854

# Run in R in directory /srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/04_rda/#RELEVANT RUN FOLDER

# in R
# Load packages
library(tidyverse)
library(psych)
library(adespatial)
library(vegan)
library(rlist)
library(adegenet)

setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/04_rda/26_rda26_5x_n24_POC_B_m8_9_10_PC3_fishing")

# load candidate SNPs (n = 854)
candidate_snps_5x <-read.PLINK("5x_TOA_only_filtered_SNPs_unlinked_candidate_SNPs_SG_B_854.raw")
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

dim(candidate_snps_5x)
# [1]  24 854

sum(is.na(candidate_snps_5x))
# [1] 0
# Warning message:
#   In is.na(snps_unlinked_5x) :
#   is.na() applied to non-(list or vector) of type 'S4'

# Impute missing SNPs
candidate_snps_imp_5x <- apply(candidate_snps_5x, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
# takes <1 second

dim(candidate_snps_imp_5x)
# [1]  24 854

sum(is.na(candidate_snps_imp_5x))
# [1] 0

candidate_AF <- makefreq(candidate_snps_imp_5x)

# candidate FORWARD SELECTION#########################################################################################################

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/rda24_5x_n24_SG_B_m8_9_10_PC3_fishing")
# Load predictor variables
envs <- read.csv("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing.csv")

str(envs)
# 'data.frame':   24 obs. of  4 variables:
# $ PC1            : num  -1.32 2.44 2.26 2.26 2.26 ...
# $ PC2            : num  -2.275 -0.739 -0.327 -0.327 -0.327 ...
# $ PC3            : num  0.038 -0.425 -0.244 -0.244 -0.244 ...
# $ fishing_ind_sum: int  398 0 0 0 0 0 398 398 398 356 ...

#Forward selection environmental variables
cand_sel_fs_envs <- forward.sel(Y = candidate_snps_imp_5x, X = envs)
# Testing variable 1
# Testing variable 2
# Testing variable 3
# Testing variable 4
# Procedure stopped (alpha criteria): pvalue for variable 4 is 0.215000 (> 0.050000)
# takes <5 seconds

cand_sel_fs_envs
#           variables order         R2     R2Cum  AdjR2Cum         F pvalue
# 1             PC1     1 0.52165381 0.5216538 0.4999108 23.991795  0.001
# 2             PC3     3 0.06733777 0.5889916 0.5498479  3.440545  0.001
# 3 fishing_ind_sum     4 0.04125404 0.6302456 0.5747825  2.231429  0.008

pval.adj <- p.adjust(cand_sel_fs_envs$pval, method = 'holm', n = 18)
cand_sel_fs_envs$pval.adj <- pval.adj
cand_sel_fs_envs
#           variables order       R2     R2Cum  AdjR2Cum         F pvalue
# 1             PC1     1 0.52165381 0.5216538 0.4999108 23.991795  0.001
# 2             PC3     3 0.06733777 0.5889916 0.5498479  3.440545  0.001
# 3 fishing_ind_sum     4 0.04125404 0.6302456 0.5747825  2.231429  0.008
#     pval.adj
# 1    0.018
# 2    0.018
# 3    0.128

# candidate COLLINNEARITY#########################################################################################################

# check the level of correlation between predictors (should be |r| < 0.7)
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_candidate_pair_panel_collinearity_plot_envs.pdf")
pairs.panels(envs, scale=T)
dev.off()

# no cross-correlation between PCs or rands whatsover

# candidate CONSTRAINED RDAs#########################################################################################################

#pRDA_constrained 1 all PCs (n = 3) + fishing (ind sum) as condition
prda_can_1 <- rda(candidate_snps_imp_5x ~ PC1 + PC2 + PC3 + Condition(fishing_ind_sum) , data=envs, scale=T)
# takes <1 second

prda_can_1
# Call: rda(formula = candidate_snps_imp_5x ~ PC1 + PC2 + PC3 +
#             Condition(fishing_ind_sum), data = envs, scale = T)
# 
# -- Model Summary --
#   
#                   Inertia Proportion Rank
# Total         854.00000    1.00000     
# Conditional    32.99839    0.03864    1
# Constrained   531.04111    0.62183    3
# Unconstrained 289.96050    0.33953   19
# 
# Inertia is correlations
# 
# -- Eigenvalues --
#   
#   Eigenvalues for constrained axes:
#   RDA1  RDA2  RDA3 
# 490.0  25.9  15.2 
# 
# Eigenvalues for unconstrained axes:
#   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
# 56.05 34.35 26.04 23.50 21.31 15.24 12.44 11.50 
# (Showing 8 of 19 unconstrained eigenvalues)

# VIF - variance inflation factor
vif.cca(prda_can_1)
# fishing_ind_sum     PC1             PC2             PC3 
# 1.037914        1.016953        1.003203        1.017757 

# VIF all <3! will proceed 

summary(eigenvals(prda_can_1, model = "constrained"))
# Importance of components:
#                          RDA1     RDA2     RDA3
# Eigenvalue            489.9963 25.86450 15.18029
# Proportion Explained    0.9227  0.04871  0.02859
# Cumulative Proportion   0.9227  0.97141  1.00000

RsquareAdj(prda_can_1)
# $r.squared
# [1] 0.621828
# 
# $adj.r.squared
# [1] 0.5940457

list.save(prda_can_1, "rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_candidate_prda_can_1.rds")

anova.cca(prda_can_1, permutations = 1000, by="terms")
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# Model: rda(formula = candidate_snps_imp_5x ~ PC1 + PC2 + PC3 + Condition(fishing_ind_sum), data = envs, scale = T)
#           Df Variance       F   Pr(>F)    
# PC1       1   460.73 30.1899 0.000999 ***
# PC2       1    15.71  1.0296 0.292707    
# PC3       1    54.60  3.5776 0.000999 ***
#   Residual 19   289.96                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# candidate PLOT RDAs#########################################################################################################

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

# plot prda_can_1 by pop1 
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_candidate_prda_can_1_pop1.pdf")
plot(prda_can_1, type="n", scaling=3, xlab="RDA1 92.3%", ylab="RDA2 4.9%")
points(prda_can_1, display="sites", pch=21, cex=1, scaling=3, bg=bg_pop1_all)
text(prda_can_1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("topright", legend=c("48","58","88"), bty="o", col="gray32", pch=21, cex=1, pt.bg=bg_pop1)
dev.off()

# plot prda_can_1 by cohort2 
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_candidate_prda_can_1_cohort2.pdf")
plot(prda_can_1, type="n", scaling=3, xlab="RDA1 92.3%", ylab="RDA2 4.9%")
points(prda_can_1, display="sites", pch=21, cex=1, scaling=3, bg=bg_cohort2_all)
text(prda_can_1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("topright", legend=c("pre2000","post2000"), bty="o", col="gray32", pch=21, cex=1, pt.bg=bg_cohort2)
dev.off()

# screeplot of prda_can_1
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_candidate_prda_can_1_screeplot.pdf")
screeplot(prda_can_1)
dev.off()