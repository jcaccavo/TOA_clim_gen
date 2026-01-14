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

# Run in R in directory /srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/04_rda/#RELEVANT RUN FOLDER

# in R
# Load packages
library(tidyverse)
library(psych)
library(adespatial)
library(vegan)
library(rlist)
library(adegenet)

# Load SNPs
snps_unlinked_5x <-read.PLINK("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/03_gen_data/5x_TOA_only_filtered_SNPs_all_unlinked.raw")
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
envs <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/rda24_5x_n24_SG_B_m8_9_10_PC3_fishing.csv")

str(envs)
# 'data.frame':   24 obs. of  4 variables:
# $ PC1            : num  -1.32 2.44 2.26 2.26 2.26 ...
# $ PC2            : num  -2.275 -0.739 -0.327 -0.327 -0.327 ...
# $ PC3            : num  0.038 -0.425 -0.244 -0.244 -0.244 ...
# $ fishing_ind_sum: int  398 0 0 0 0 0 398 398 398 356 ...

#Forward selection environmental variables
sel_fs_envs <- forward.sel(Y = snps_unlinked_imp_5x, X = envs)
# Testing variable 1
# Testing variable 2
# Procedure stopped (alpha criteria): pvalue for variable 2 is 0.450000 (> 0.050000)

sel_fs_envs
# variables order         R2      R2Cum    AdjR2Cum        F pvalue
# 1       PC1     1 0.05148626 0.05148626 0.008372001 1.194182  0.001

pval.adj <- p.adjust(sel_fs_envs$pval, method = 'holm', n = 18)
sel_fs_envs$pval.adj <- pval.adj
sel_fs_envs
# variables order         R2      R2Cum    AdjR2Cum        F pvalue pval.adj
# 1       PC1     1 0.05148626 0.05148626 0.008372001 1.194182  0.001    0.018

#COLLINNEARITY#########################################################################################################

# check the level of correlation between predictors (should be |r| < 0.7)
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_pair_panel_collinearity_plot_envs.pdf")
pairs.panels(envs, scale=T)
dev.off()

# no cross-correlation between PCs or rands whatsover

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
# Constrained   3.242e+05  1.374e-01    3
# Unconstrained 1.927e+06  8.169e-01   19
# 
# Inertia is correlations
# 
# -- Eigenvalues --
#   
#   Eigenvalues for constrained axes:
#   RDA1   RDA2   RDA3 
# 129987  99249  94940 
# 
# Eigenvalues for unconstrained axes:
#   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8 
# 152455 109451 105187 104481 103011 101979 100439 100092 
# (Showing 8 of 19 unconstrained eigenvalues)

# VIF - variance inflation factor
vif.cca(prda_1)
# fishing_ind_sum    PC1             PC2             PC3 
# 1.037914        1.016953        1.003203        1.017757 

# VIF all <3! will proceed 

summary(eigenvals(prda_1, model = "constrained"))
# Importance of components:
#                         RDA1      RDA2      RDA3
# Eigenvalue            1.30e+05 9.925e+04 9.494e+04
# Proportion Explained  4.01e-01 3.062e-01 2.929e-01
# Cumulative Proportion 4.01e-01 7.071e-01 1.000e+00

RsquareAdj(prda_1)
# $r.squared
# [1] 0.1374206
# 
# $adj.r.squared
# [1] 0.008815339

list.save(prda_1, "rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_prda_1.rds")

anova.cca(prda_1, permutations = 1000, by="terms")
# Permutation test for rda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 1000
# 
# Model: rda(formula = snps_unlinked_imp_5x ~ PC1 + PC2 + PC3 + Condition(fishing_ind_sum), data = envs, scale = T)
# Df Variance      F   Pr(>F)
# PC1       1   128254 1.2645 0.001998 **
# PC2       1    99279 0.9788 0.596404
# PC3       1    96644 0.9528 0.898102
# Residual 19  1927139
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
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_prda_1_pop1.pdf")
plot(prda_1, type="n", scaling=3, xlab="RDA1 40.1%", ylab="RDA2 30.6%")
points(prda_1, display="sites", pch=21, cex=1, scaling=3, bg=bg_pop1_all)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("48","58","88"), bty="o", col="gray32", pch=21, cex=1, pt.bg=bg_pop1)
dev.off()

# plot prda_1 by cohort2 
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_prda_1_cohort2.pdf")
plot(prda_1, type="n", scaling=3, xlab="RDA1 40.1%", ylab="RDA2 30.6%")
points(prda_1, display="sites", pch=21, cex=1, scaling=3, bg=bg_cohort2_all)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("pre2000","post2000"), bty="o", col="gray32", pch=21, cex=1, pt.bg=bg_cohort2)
dev.off()

# screeplot of prda_1
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_prda_1_screeplot.pdf")
screeplot(prda_1)
dev.off()

# CANDIDATE SNPS#########################################################################################################

# useful links
# https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html#conclusions

# Load predictor variables
envs <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/rda24_5x_n24_SG_B_m8_9_10_PC3_fishing.csv")

str(envs)
# 'data.frame':   24 obs. of  4 variables:
# $ PC1            : num  -1.32 2.44 2.26 2.26 2.26 ...
# $ PC2            : num  -2.275 -0.739 -0.327 -0.327 -0.327 ...
# $ PC3            : num  0.038 -0.425 -0.244 -0.244 -0.244 ...
# $ fishing_ind_sum: int  398 0 0 0 0 0 398 398 398 356 ...

# subset predictor variables to only include variables identified as significant from the ANOVA, as well as fake variables
envs_retained <- subset(envs, select=c(PC1))

str(envs_retained)
# 'data.frame':   24 obs. of  1 variable:
# $ PC1: num  -1.32 2.44 2.26 2.26 2.26 ...

# if continuing directly from RDAs, skip loading the RDA and SNPs

# ***********************************************************************************************************************

# if coming back to an old analysis in a new R workspace, you must load the RDA, as well as the SNPs, and impute them

# Load selected RDA
prda_1 <- read_rds("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_prda_1.rds")

# Load SNPs
snps_unlinked_5x <-read.PLINK("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/03_gen_data/5x_TOA_only_filtered_SNPs_all_unlinked.raw")
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
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_prda_1_load_constrained_1.pdf")
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
# [1] 1010

# RDA 2
cand2 <- outliers(load.rda[,2],3)
length(cand2)
# [1] 2254

# RDA 3
cand3 <- outliers(load.rda[,3],3)
length(cand3)
# [1] 14664

# PC 1
cand4 <- outliers(load.rda[,4],3)
length(cand4)
# [1] 21

# PC 2
cand5 <- outliers(load.rda[,5],3)
length(cand5)
# [1] 1811

# PC 3
cand6 <- outliers(load.rda[,6],3)
length(cand6)
# [1] 2315

# PC 4
cand7 <- outliers(load.rda[,7],3)
length(cand7)
# [1] 927

# PC 5
cand8 <- outliers(load.rda[,8],3)
length(cand8)
# [1] 2587

# PC 6
cand9 <- outliers(load.rda[,9],3)
length(cand9)
# [1] 2966

# PC 7
cand10 <- outliers(load.rda[,10],3)
length(cand10)
# [1] 14229

# PC 8
cand11 <- outliers(load.rda[,11],3)
length(cand11)
# [1] 5210

# PC 9
cand12 <- outliers(load.rda[,12],3)
length(cand12)
# [1] 3285

# PC 10
cand13 <- outliers(load.rda[,13],3)
length(cand13)
# [1] 6348

# PC 11
cand14 <- outliers(load.rda[,14],3)
length(cand14)
# [1] 2091

# PC 12
cand15 <- outliers(load.rda[,15],3)
length(cand15)
# [1] 11733

# PC 13
cand16 <- outliers(load.rda[,16],3)
length(cand16)
# [1] 9315

# PC 14
cand17 <- outliers(load.rda[,17],3)
length(cand17)
# [1] 2482

# PC 15
cand18 <- outliers(load.rda[,18],3)
length(cand18)
# [1] 5348

# PC 16
cand19 <- outliers(load.rda[,19],3)
length(cand19)
# [1] 8523

# PC 17
cand20 <- outliers(load.rda[,20],3)
length(cand20)
# [1] 10412

# PC 18
cand21 <- outliers(load.rda[,21],3)
length(cand21)
# [1] 3225

# PC 19
cand22 <- outliers(load.rda[,22],3)
length(cand22)
# [1] 3671

# all output values (for future comparison)
# [1] 1010
# [1] 2254
# [1] 14664
# [1] 21
# [1] 1811
# [1] 2315
# [1] 927
# [1] 2587
# [1] 2966
# [1] 14229
# [1] 5210
# [1] 3285
# [1] 6348
# [1] 2091
# [1] 11733
# [1] 9315
# [1] 2482
# [1] 5348
# [1] 8523
# [1] 10412
# [1] 3225
# [1] 3671

# Candidate total from both RDA (3) and PC (20) axes
ncand_all <- length(cand1) + length(cand2) + length(cand3) + length(cand4) + length(cand5) + length(cand6) + length(cand7) + length(cand8) + length(cand9) + length(cand10) + length(cand11) + length(cand12) + length(cand13) + length(cand14) + length(cand15) + length(cand16) + length(cand17) + length(cand18) + length(cand19) + length(cand20) + length(cand21) + length(cand22)
ncand_all
# [1] 114427

# Candidate total from RDA axes only
ncand_rda3 <- length(cand1) + length(cand2) + length(cand3)
ncand_rda3
# [1] 17928

# Candidate total from RDA axes 1 and 2 only
ncand_rda2 <- length(cand1) + length(cand2)
ncand_rda2
# [1] 3264

# Candidate total from RDA axis 1 only
ncand_rda1 <- length(cand1)
ncand_rda1
# [1] 1010

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
pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_prda_1_SNP_zoom_axes1_2.pdf")
plot(prda_1, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), main="rda24.prda1, axes 1 and 2", xlab="RDA1 40.1%", ylab="RDA2 30.6%")
points(prda_1, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(prda_1, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
dev.off()

pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_prda_1_SNP_zoom_axes2_3.pdf")
plot(prda_1, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="rda24.prda1, axes 2 and 3", xlab="RDA1 40.1%", ylab="RDA2 30.6%")
points(prda_1, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(prda_1, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
text(prda_1, scaling=3, display="bp", col="#0868ac", cex=1)
dev.off()

# Let’s see which environmental predictors are most strongly correlated with the first three RDA axes
intersetcor(prda_1)[,1:3]
#         RDA1       RDA2        RDA3
# PC1  0.93832263 0.04197081 -0.21886628
# PC2 -0.03965053 0.99761377  0.04084409
# PC3 -0.22766396 0.02261746 -0.97124622

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

# add in the correlations of each candidate SNP with predictor (n = 1) PC1
foo <- matrix(nrow=(ncand_rda3), ncol=1)  # 1 column for 1 predictor
colnames(foo) <- c("PC1")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps_unlinked_imp_5x[,nam]
  foo[i,] <- apply(envs_retained,2,function(x) cor(x,snp.gen))
}
# may take a few minutes

cand <- cbind.data.frame(cand,foo)  
head(cand)
#     axis                   snp     loading        PC1
# 1    1  HiC_scaffold_10:564698_G  0.04553554  0.7923422
# 2    1 HiC_scaffold_10:1503610_G  0.04172289  0.6737594
# 3    1 HiC_scaffold_10:2543556_T  0.04413212  0.7915214
# 4    1 HiC_scaffold_10:2711464_A  0.04205947  0.7606572
# 5    1 HiC_scaffold_10:3231239_G -0.03366827 -0.6494701
# 6    1 HiC_scaffold_10:3720955_G -0.03402031 -0.5867432

length(cand$snp[duplicated(cand$snp)])
# takes <1 minute
# [1] 0

foo <- cbind(cand$axis, duplicated(cand$snp))
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
# 0 
# 1010
table(foo[foo[,1]==2,2]) # no duplicates on axis 2
# 0 
# 2254
table(foo[foo[,1]==3,2]) # no duplicates on axis 3
# 0 
# 14664

# to which predictor is each candidate SNP most strongly correlated 
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,5] <- names(which.max(abs(bar[4:4]))) # gives the variable | column 4 starts the variables, so 4:n where n depends on how many variables you have | and i,k where k is the new column you will add, so k = n + 1
  cand[i,6] <- max(abs(bar[4:4]))              # gives the correlation | in the brackets [] include the same columns indicated in the line above in brackets | and i,m where m is the second new column you will add, so m = k + 1
}
colnames(cand)[5] <- "predictor"
colnames(cand)[6] <- "correlation"

table(cand$predictor)
# PC1 
# 17928 

# this is the total candidate SNPs explained all RDA axes
# it remains an open question if it makes sense to include the non-significant variables (predictors) to isolate which SNPs go with PC1 versus the other 2..

write.csv(cand,"rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_prda_1_putative_adaptive_SNPs.csv")

# ****UP TO HERE**** ####
# running in jc_screen_1

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

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/rda24_5x_n24_SG_B_m8_9_10_PC3_fishing")

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

# BOURRET et al. analyses ####

# packages ####
library(tidyverse)
library(here)
library(ggpubr)
library(ggrepel)
library(ggcorrplot)
library(ggridges)
library(parallel)

library("ggVennDiagram")
library(eulerr)
library(RColorBrewer)

library(multcompView)

library(adegenet)
library(hierfstat)
# BiocManager::install("dartR")
# BiocManager::install("SNPRelate")
library(SNPRelate)
library(dartR)

library(pcadapt)
#BiocManager::install("qvalue")
library("qvalue")
library(robust)

`%nin%` = Negate(`%in%`)

# BiocManager::install("QuickPop")
# library(QuickPop) # COULD NOT INSTALL

library(raster)
library(corrplot)
library(vegan)

#Libraries that we will need
library(codep)
library(adespatial)
library(adegraphics)
library(ape)
library(car)

library(raster)
# BiocManager::install("rgdal")
# library(rgdal) # COULD NOT INSTALL locally - this package has been replaced by "sf" and "terra"
library(sf)
# BiocManager::install("rgeos")
# library(rgeos) # COULD NOT INSTALL locally
library(stringr)
library(rnaturalearth)
library(rnaturalearthdata)

library(terra)
library(tidyterra)

# BiocManager::install("rEEMSplots")
# library(rEEMSplots) # COULD NOT INSTALL

#detach("package:ggVennDiagram", character.only = T)
library(sp)
library(rworldmap)
library(rworldxtra)

# Genotype-environement association : identifying loci under selec --------

# skipping 1st 2 steps as we've already done this with the pRDA
# RDA_env.constrained   <- rda(freq.MAF.final.Gen_ZONE_FG ~   Tbtm.Winter + Tbtm.Summer + Sbtm.ann +  SSS.ann + SST.Larval + Condition(PC1 + PC2 + PC3),  Variables)
# plot(RDA_env.constrained, scaling = 2)

rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

rdadapt_env.constrained   <-rdadapt(prda_can_1, 2)

## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_env.constrained$p.values)

## Identifying the loci that are below the p-value threshold
outliers.constrained <- data.frame(Loci = colnames(candidate_snps_imp_5x)[which(rdadapt_env.constrained$p.values<thres_env)], 
                                   p.value = rdadapt_env.constrained$p.values[which(rdadapt_env.constrained$p.values<thres_env)])
outliers.constrained    
outliers.constrained  %>% nrow()
# [1] 155

## List of outlier names
outliers_rdadapt_env <- as.character(outliers.constrained$Loci)

RDA_env <- prda_can_1

pdf("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing_bourret_prda_can_1_screeplot.pdf")
screeplot(RDA_env, main="Eigenvalues of constrained axes")
dev.off()

## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores) %>% 
  mutate(type = ifelse(names %in% outliers_rdadapt_env, "Outlier", "Neutral"))

TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) %>% 
  mutate(ID = row.names(scores(RDA_env, choices=c(1,2), display="bp")) ) %>% left_join(envs)
TAB_var  # pull the biplot scores

# Error in `left_join()`:                                                                                                                                                                                                                                   
# ! `by` must be supplied when `x` and `y` have no common variables.                                                                                                                                                                                        
# ℹ Use `cross_join()` to perform a cross-join.                                                                                                                                                                                                             
# Run `rlang::last_trace()` to see where the error occurred.

# load Bourret's raster data to troubleshoot issue with TAB_var

# Env data from raster directly -------------------------------------------

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/genomic_offset/examples/bourret_scripts_data/00_Data/99_SIG")

## Loading the climatic rasters ####
ras_all <- stack(list.files(".", pattern = ".tif", full.names = T))
# Error in .local(x, ...) : 
#   lazy-load database '/srv/public/users/jcaccavo/miniconda3/lib/R/library/raster/R/raster.rdb' is corrupt
# In addition: Warning messages:
#   1: In get(var, frame, inherits = FALSE) :
#   internal error -3 in R_decompress1
# 2: In .local(x, ...) : restarting interrupted promise evaluation
# 3: In .local(x, ...) : internal error -3 in R_decompress1

# cannot import raster data - will try to skip the TAB_var step, hopefully it's not critical..

# this ^ was on the server; locally the raster loading worked
names(ras_all) <- names(ras_all) %>% str_remove("mon_")
names(ras_all) 
# [1] "Sbtm.ann_2075"    "Sbtm.ann"         "SSS.ann_2075"     "SSS.ann"          "SST.Larval_2075" 
# [6] "SST.Larval"       "Tbtm.Summer_2075" "Tbtm.Summer"      "Tbtm.Winter_2075" "Tbtm.Winter"    
ras_current <- ras_all[[str_subset(names(ras_all), "2075", negate = T)]]
names(ras_current)
# [1] "Sbtm.ann"    "SSS.ann"     "SST.Larval"  "Tbtm.Summer" "Tbtm.Winter"
ras_future <- ras_all[[str_subset(names(ras_all), "2075", negate = F)]]
names(ras_future) <- names(ras_future) %>% str_remove("_2075")
names(ras_future)
# [1] "Sbtm.ann"    "SSS.ann"     "SST.Larval"  "Tbtm.Summer" "Tbtm.Winter"

remove.NAs.stack<-function(rast.stack){
  nom<-names(rast.stack)
  test1<-calc(rast.stack, fun=sum)
  test1[!is.na(test1)]<-1
  test2<-rast.stack*test1
  test2<-stack(test2)
  names(test2)<-nom
  return(test2)
}

ras_current <- remove.NAs.stack(ras_current)
ras_current
# class      : RasterStack 
# dimensions : 333, 444, 147852, 5  (nrow, ncol, ncell, nlayers)
# resolution : 0.09, 0.09  (x, y)
# extent     : -79.99998, -40.03998, 40.02999, 69.99999  (xmin, xmax, ymin, ymax)
# crs        : +proj=longlat +datum=WGS84 +no_defs 
# names      :   Sbtm.ann,    SSS.ann, SST.Larval, Tbtm.Summer, Tbtm.Winter 
# min values :  0.9844299,  0.3367375, -0.8128369,  -1.7610348,  -2.1381688 
# max values :   35.03660,   35.67757,   19.54368,    23.49642,    10.16929 
ras_future  <- remove.NAs.stack(ras_future)
ras_future
# class      : RasterStack 
# dimensions : 333, 444, 147852, 5  (nrow, ncol, ncell, nlayers)
# resolution : 0.09, 0.09  (x, y)
# extent     : -79.99998, -40.03998, 40.02999, 69.99999  (xmin, xmax, ymin, ymax)
# crs        : +proj=longlat +datum=WGS84 +no_defs 
# names      :   Sbtm.ann,    SSS.ann, SST.Larval, Tbtm.Summer, Tbtm.Winter 
# min values :  0.8313547,  0.1447514, -3.2504046,  -1.6631174,  -2.0874052 
# max values :   35.23787,   35.69218,   20.32519,    24.73203,    13.40221 

## Standardization of the variables
## Recovering scaling coefficients
scale_env <- attr(scale(Env.data[,-1], center=TRUE, scale=TRUE), 'scaled:scale')
center_env <- attr(scale(Env.data[,-1], center=TRUE, scale=TRUE), 'scaled:center')

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/genomic_offset/examples/bourret_scripts_data/00_Data/00_FIleInfos")
pop.data <- read_csv("Projet_Infos.csv") 

region.df <- tibble(Gen_ZONE = pop.data$Gen_ZONE %>% unique()) %>% 
  mutate(RegionAssesment = ifelse(Gen_ZONE %in% c("SFA-13", "SFA-14", "SFA-15", "SFA-16"), "Scotian",
                                  ifelse(Gen_ZONE %in% c("SFA-8", "SFA-9", "SFA-10", "SFA-12"), "GSL",
                                         ifelse(Gen_ZONE %in% c("SFA-4","SFA-5", "SFA-6", "SFA-7"), "NFL",
                                                ifelse(Gen_ZONE %in% c("SFA-0", "WAZ", "EAZ"), "Northern",
                                                       ifelse(Gen_ZONE %in% "NAFO-3M", "NAFO-3M",NA))))))

region.df
# A tibble: 18 × 2
# Gen_ZONE RegionAssesment
# <chr>    <chr>          
#   1 SGSLZ    NA             
# 2 SFA-16   Scotian        
# 3 SFA-7    NFL            
# 4 SFA-15   Scotian        
# 5 SFA-13   Scotian        
# 6 SFA-14   Scotian        
# 7 NAFO-3M  NAFO-3M        
# 8 SFA-9    GSL            
# 9 SFA-12   GSL            
# 10 SFA-10   GSL            
# 11 SFA-8    GSL            
# 12 WAZ      Northern       
# 13 SFA-4    NFL            
# 14 EAZ      Northern       
# 15 SFA-0    Northern       
# 16 SFA-5    NFL            
# 17 SFA-6    NFL            
# 18 NA       NA 

pop.data <- pop.data %>% left_join(region.df) %>% 
  mutate(Gen_ZONE = ifelse(Gen_ZONE == "SFA-0", "EAZ", Gen_ZONE),
         Gen_ZONE = factor(Gen_ZONE, levels = c("WAZ", "EAZ", "SFA-4",
                                                "SFA-5", "SFA-6", "SFA-7", "SFA-8",
                                                "SFA-9", "SFA-10", "SFA-12",
                                                "SFA-13", "SFA-14", "SFA-15", "SFA-16",
                                                "NAFO-3M")),
         Gen_ZONE_FG = factor(Gen_ZONE_FG, levels = c("SFA-0-1",
                                                      paste(rep(levels(Gen_ZONE), each = 10), 1:10 , sep = "-"))),
         Gen_ZONE_FG = factor(Gen_ZONE_FG),
         Transcripto = ifelse(Gen_ZONE_FG %in% c("SFA-6-3", "SFA-8-4", "SFA-12-2", "SFA-12-3", "SFA-15-1"), "Yes", "No")
  )

pop.data
# A tibble: 1,920 × 20
# ID_GQ       Plaque Puit  Espece Espece_revision Cat_Sample Barcode ID_labo  Year Sexe  Multipart Egg     Lat
# <chr>         <dbl> <chr> <chr>  <chr>           <chr>      <chr>   <chr>   <dbl> <chr> <chr>     <chr> <dbl>
# 1 Dp_2064       15 B07   Dp     Dp              Sample     GCCAGT  Dp_2064  2018 NA    NA        NA     45.9
# 2 Dp_2065       18 F05   Dp     Dp              Sample     GCTCTA  Dp_2065  2018 M     N         N      45.9
# 3 Dp_2066       18 A07   Dp     Pm              Sample     CTATTA  Dp_2066  2018 NA    NA        NA     45.9
# 4 Dp_2067       14 C01   Dp     Dp              Sample     ACTA    Dp_2067  2018 NA    NA        NA     45.9
# 5 Dp_2074       18 A11   Dp     Dp              Sample     TGCAAG… Dp_2074  2019 M     N         N      44.8
# 6 Dp_2075       12 H09   Dp     Dp              Sample     CTACGGA Dp_2075  2019 F     NA        N      44.8
# 7 Dp_2077       12 C12   Dp     Dp              Sample     CGTGTG… Dp_2077  2019 M     N         N      44.8
# 8 Dp_2086       18 D08   Dp     Dp              Sample     GAATTCA Dp_2086  2019 M     N         N      44.8
# 9 Dp_2088       17 C01   Dp     Dp              Sample     ACTA    Dp_2088  2019 M     N         N      44.8
# 10 Dp_2089      12 A11   Dp     Dp              Sample     TGCAAG… Dp_2089  2019 F     NA        N      44.8
# ℹ 1,910 more rows
# ℹ 7 more variables: Long <dbl>, Gen_ZONE <fct>, Gen_ZONE_FG <fct>, Plaque_ID <chr>, POP <chr>,
#   RegionAssesment <chr>, Transcripto <chr>
# ℹ Use `print(n = ...)` to see more rows

head(pop.data$Gen_ZONE_FG)
# [1] <NA>     <NA>     <NA>     <NA>     SFA-16-2 SFA-16-2
# 60 Levels: SFA-0-1 WAZ-1 WAZ-2 WAZ-3 EAZ-1 EAZ-2 EAZ-3 EAZ-4 SFA-4-1 SFA-4-2 SFA-5-1 SFA-5-2 SFA-5-3 SFA-5-4 SFA-6-1 SFA-6-2 SFA-6-3 SFA-6-4 SFA-6-5 ... NAFO-3M-3

pop.too.small <-c("SFA-0-1", "WAZ-3", "SFA-7-7")

# Genetic Data ------------------------------------------------------------

load("00_Data/06_Filtering/06g_UniqueFinal/populations.14331snps.1513ind.n13HW.DP.r5.single.final.adegenet.Rdata")
# not possible to do as the original data are not part of the zenodo data assembly

# Adaptative index - Constrained ------------------------------------------

## Loading the climatic rasters

newproj  = "+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/genomic_offset/examples/bourret_scripts_data/00_Data/99_SIG")
## Species range shapefile (download information at beginning of tutorial)
range <- read_sf("range2.shp") 
range <- read_sf("add_coastline_high_res_line_v7_10.shp")
crs(range) <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

res_RDA_proj_current <- adaptive_index(RDA = prda_can_1, K = 2, env_pres = ras_current,  range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
# cannot run the adaptive index analyses because the Env.data object, which was used to create scale_env and center_env used the genomic data files, which were not included in the online data

# trying to convert my own .nc files to .tif files to have raster versions

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/environmental_data/CMEMS/ORA5_reanalysis/1985")
# sla = raster("sosaline_control_monthly_highres_2D_198512_CONS_v0.1.nc", level=

sal198512 <- terra::rast("sosaline_control_monthly_highres_2D_198512_CONS_v0.1.nc")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
sal198512
# class       : SpatRaster 
# dimensions  : 1021, 1442, 1  (nrow, ncol, nlyr)
# resolution  : 1, 1  (x, y)
# extent      : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
# source      : sosaline_control_monthly_highres_2D_198512_CONS_v0.1.nc:sosaline 
# varname     : sosaline (Sea Surface Salinity) 
# name        : sosaline 
# unit        :      PSU 
sit198512 <- terra::rast("iicethic_control_monthly_highres_2D_198512_CONS_v0.1.nc")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
sit198512
# class       : SpatRaster 
# dimensions  : 1021, 1442, 1  (nrow, ncol, nlyr)
# resolution  : 1, 1  (x, y)
# extent      : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
# source      : iicethic_control_monthly_highres_2D_198512_CONS_v0.1.nc:iicethic 
# varname     : iicethic (Ice thickness) 
# name        : iicethic 
# unit        :        m 
sic198512 <- terra::rast("ileadfra_control_monthly_highres_2D_198512_CONS_v0.1.nc")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
sic198512
# class       : SpatRaster 
# dimensions  : 1021, 1442, 1  (nrow, ncol, nlyr)
# resolution  : 1, 1  (x, y)
# extent      : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
# source      : ileadfra_control_monthly_highres_2D_198512_CONS_v0.1.nc:ileadfra 
# varname     : ileadfra (Ice concentration) 
# name        : ileadfra 
# unit        : Fraction 
mld1_198512 <- terra::rast("somxl010_control_monthly_highres_2D_198512_CONS_v0.1.nc")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
mld1_198512
# class       : SpatRaster 
# dimensions  : 1021, 1442, 1  (nrow, ncol, nlyr)
# resolution  : 1, 1  (x, y)
# extent      : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
# source      : somxl010_control_monthly_highres_2D_198512_CONS_v0.1.nc:somxl010 
# varname     : somxl010 (Mixed Layer Depth 0.01) 
# name        : somxl010 
# unit        :        m 
mld3_198512 <- terra::rast("somxl030_control_monthly_highres_2D_198512_CONS_v0.1.nc")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
mld3_198512
# class       : SpatRaster 
# dimensions  : 1021, 1442, 1  (nrow, ncol, nlyr)
# resolution  : 1, 1  (x, y)
# extent      : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
# source      : somxl030_control_monthly_highres_2D_198512_CONS_v0.1.nc:somxl030 
# varname     : somxl030 (Mixed Layer Depth 0.03) 
# name        : somxl030 
# unit        :        m 
sst_198512 <- terra::rast("sosstsst_control_monthly_highres_2D_198512_CONS_v0.1.nc")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
sst_198512
# class       : SpatRaster 
# dimensions  : 1021, 1442, 1  (nrow, ncol, nlyr)
# resolution  : 1, 1  (x, y)
# extent      : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
# source      : sosstsst_control_monthly_highres_2D_198512_CONS_v0.1.nc:sosstsst 
# varname     : sosstsst (Sea Surface temperature) 
# name        : sosstsst 
# unit        :        C 
cm_198512 <- terra::rast("vomecrty_control_monthly_highres_3D_198512_CONS_v0.1.nc")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
cm_198512
# class       : SpatRaster 
# dimensions  : 1021, 1442, 75  (nrow, ncol, nlyr)
# resolution  : 1, 1  (x, y)
# extent      : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
# source      : vomecrty_control_monthly_highres_3D_198512_CONS_v0.1.nc:vomecrty 
# varname     : vomecrty (Meridional Current) 
# names       : vomec~76001, vomec~58553, vomec~76817, vomec~62799, vomec~03613, vomec~30336, ... 
# unit        :         m/s,         m/s,         m/s,         m/s,         m/s,         m/s, ... 
cz_198512 <- terra::rast("vozocrtx_control_monthly_highres_3D_198512_CONS_v0.1.nc")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
cz_198512
# class       : SpatRaster 
# dimensions  : 1021, 1442, 75  (nrow, ncol, nlyr)
# resolution  : 1, 1  (x, y)
# extent      : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# coord. ref. : +proj=longlat +datum=WGS84 +no_defs 
# source      : vozocrtx_control_monthly_highres_3D_198512_CONS_v0.1.nc:vozocrtx 
# varname     : vozocrtx (Zonal Current) 
# names       : vozoc~76001, vozoc~58553, vozoc~76817, vozoc~62799, vozoc~03613, vozoc~30336, ... 
# unit        :         m/s,         m/s,         m/s,         m/s,         m/s,         m/s, ... 

ras_1985_all <- stack(sal198512,sit198512,sic198512,mld1_198512,mld3_198512,sst_198512,cm_198512,cz_198512)
# Error in .local(x, ...) : 
  # unused arguments (new("SpatRaster", ptr = new("Rcpp_SpatRaster", .xData = <environment>)), 
#                     new("SpatRaster", ptr = new("Rcpp_SpatRaster", .xData = <environment>)), 
#                     new("SpatRaster", ptr = new("Rcpp_SpatRaster", .xData = <environment>)), 
#                     new("SpatRaster", ptr = new("Rcpp_SpatRaster", .xData = <environment>)), 
#                     new("SpatRaster", ptr = new("Rcpp_SpatRaster", .xData = <environment>)), 
#                     new("SpatRaster", ptr = new("Rcpp_SpatRaster", .xData = <environment>)), 
#                     new("SpatRaster", ptr = new("Rcpp_SpatRaster", .xData = <environment>)))

# turn nc files into rasterstack ####

library(ncdf4)
setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/environmental_data/CMEMS/ORA5_reanalysis/1985")

fname <- file.choose()
sal198512_2 <- nc_open("sosaline_control_monthly_highres_2D_198512_CONS_v0.1.nc")
sal198512_2brick <- brick(fname, varname="sosaline")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(sal198512_2brick, filename="sal198512_2.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
sit198512_2 <- nc_open("iicethic_control_monthly_highres_2D_198512_CONS_v0.1.nc")
sit198512_2brick <- brick(fname, varname="iicethic")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(sit198512_2brick, filename="sit198512_2.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
sic198512_2 <- nc_open("ileadfra_control_monthly_highres_2D_198512_CONS_v0.1.nc")
sic198512_2brick <- brick(fname, varname="ileadfra")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(sic198512_2brick, filename="sic198512_2.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
mld1_198512_2 <- nc_open("somxl010_control_monthly_highres_2D_198512_CONS_v0.1.nc")
mld1_198512_2brick <- brick(fname, varname="somxl010")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(mld1_198512_2brick, filename="mld1_198512_2.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
mld3_198512_2 <- nc_open("somxl030_control_monthly_highres_2D_198512_CONS_v0.1.nc")
mld3_198512_2brick <- brick(fname, varname="somxl030")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(mld3_198512_2brick, filename="mld3_198512_2.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
sst_198512_2 <- nc_open("sosstsst_control_monthly_highres_2D_198512_CONS_v0.1.nc")
sst_198512_2brick <- brick(fname, varname="sosstsst")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(sst_198512_2brick, filename="sst_198512_2.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
cm_198512_2 <- nc_open("vomecrty_control_monthly_highres_3D_198512_CONS_v0.1.nc")
cm_198512_2brick <- brick(fname, varname="vomecrty")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(cm_198512_2brick, filename="cm_198512_2_3.tif", format="GTiff",overwrite=TRUE, bylayer=TRUE)

fname <- file.choose()
cz_198512_2 <- nc_open("vozocrtx_control_monthly_highres_3D_198512_CONS_v0.1.nc")
cz_198512_2brick <- brick(fname, varname="vozocrtx")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(cz_198512_2brick, filename="cz_198512_2_3.tif", format="GTiff",overwrite=TRUE, bylayer=TRUE)

# load 1985 raster stack ####

# if running locally
setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/environmental_data/CMEMS/ORA5_reanalysis/1985")
# if running on the server
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/29_genomic_offset/01_env_rasters/a_1985")

ras_all_1985 <- stack(list.files(".", pattern = ".tif", full.names = T))

names(ras_all_1985)
# [1] "cm_198512_2_3_1" "cz_198512_2_3_1" "mld1_198512_2"   "mld3_198512_2"   "sal198512_2"     "sic198512_2"     "sit198512_2"     "sst_198512_2" 

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/environmental_data/CMEMS/ORA5_reanalysis/2022")

fname <- file.choose()
sal2022 <- nc_open("sosaline_control_monthly_highres_2D_202212_OPER_v0.1.nc")
sal2022brick <- brick(fname, varname="sosaline")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(sal2022brick, filename="sal2022.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
sit2022 <- nc_open("iicethic_control_monthly_highres_2D_202212_OPER_v0.1.nc")
sit2022brick <- brick(fname, varname="iicethic")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(sit2022brick, filename="sit2022.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
sic2022 <- nc_open("ileadfra_control_monthly_highres_2D_202212_OPER_v0.1.nc")
sic2022brick <- brick(fname, varname="ileadfra")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(sic2022brick, filename="sic2022.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
mld1_2022 <- nc_open("somxl010_control_monthly_highres_2D_202212_OPER_v0.1.nc")
mld1_2022brick <- brick(fname, varname="somxl010")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(mld1_2022brick, filename="mld1_2022.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
mld3_2022 <- nc_open("somxl030_control_monthly_highres_2D_202212_OPER_v0.1.nc")
mld3_2022brick <- brick(fname, varname="somxl030")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(mld3_2022brick, filename="mld3_2022.tif", format="GTiff",overwrite=TRUE)

fname <- file.choose()
sst_2022 <- nc_open("sosstsst_control_monthly_highres_2D_202212_OPER_v0.1.nc")
sst_2022brick <- brick(fname, varname="sosstsst")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(sst_2022brick, filename="sst_2022.tif", format="GTiff",overwrite=TRUE)

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/environmental_data/CMEMS/ORA5_reanalysis/2022")

fname <- file.choose()
cm_2022 <- nc_open("vomecrty_control_monthly_highres_3D_202212_OPER_v0.1.nc")
cm_2022brick <- brick(fname, varname="vomecrty")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(cm_2022brick, filename="old_tgifs/cm_2022_3.tif", format="GTiff",overwrite=TRUE, bylayer=TRUE)

fname <- file.choose()
cz_2022 <- nc_open("vozocrtx_control_monthly_highres_3D_202212_OPER_v0.1.nc")
cz_2022brick <- brick(fname, varname="vozocrtx")
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named x BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
# [1] "vobjtovarid4: **** WARNING **** I was asked to get a varid for dimension named y BUT this dimension HAS NO DIMVAR! Code will probably fail at this point"
writeRaster(cz_2022brick, filename="old_tgifs/cz_2022_3.tif", format="GTiff",overwrite=TRUE, bylayer=TRUE)

# load 2022 raster stack ####

# if running locally
setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/environmental_data/CMEMS/ORA5_reanalysis/2022")
# if running on the server
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/29_genomic_offset/01_env_rasters/b_2022")

ras_all_2022 <- stack(list.files(".", pattern = ".tif", full.names = T))

names(ras_all_2022)
# [1] "cm_2022_3_1" "cz_2022_3_1" "mld1_2022"   "mld3_2022"   "sal2022"     "sic2022"     "sit2022"     "sst_2022"   

remove.NAs.stack<-function(rast.stack){
  nom<-names(rast.stack)
  test1<-calc(rast.stack, fun=sum)
  test1[!is.na(test1)]<-1
  test2<-rast.stack*test1
  test2<-stack(test2)
  names(test2)<-nom
  return(test2)
}

ras_all_1985 <- remove.NAs.stack(ras_all_1985)
ras_all_2022 <- remove.NAs.stack(ras_all_2022)

ras_all_1985
# class      : RasterStack 
# dimensions : 1021, 1442, 1472282, 8  (nrow, ncol, ncell, nlayers)
# resolution : 1, 1  (x, y)
# extent     : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# crs        : NA 
# names      : cm_198512_2_3_1, cz_198512_2_3_1, mld1_198512_2, mld3_198512_2,  sal198512_2,  sic198512_2,  sit198512_2, sst_198512_2 
# min values :     -2.13900256,     -1.86371887,    3.25130916,    3.25130916,   0.02503924,   0.00000000,   0.00000000,  -2.08541298 
# max values :       1.7417066,       1.8068730,  1507.2650146,  2591.7431641,   41.3624916,    0.9999983,    4.6301308,   31.5603733 

ras_all_2022
# class      : RasterStack 
# dimensions : 1021, 1442, 1472282, 8  (nrow, ncol, ncell, nlayers)
# resolution : 1, 1  (x, y)
# extent     : 0.5, 1442.5, 0.5, 1021.5  (xmin, xmax, ymin, ymax)
# crs        : NA 
# names      :  cm_2022_3_1,  cz_2022_3_1,    mld1_2022,    mld3_2022,      sal2022,      sic2022,      sit2022,     sst_2022 
# min values :  -1.68924916,  -1.69077957,   3.25130916,   3.25130916,   0.05512097,   0.00000000,   0.00000000,  -2.04820776 
# max values :    2.0137737,    1.9313731,  774.7990112, 1128.0826416,   43.5273705,    0.9999958,    4.6246314,   31.3035774 

Env.data <- data.frame(Gen_ZONE_FG = unique(pop(gi.final.Gen_ZONE_FG))) %>% 
  left_join(pop.data %>% dplyr::select(Gen_ZONE_FG, Lat, Long) %>% distinct(.keep_all = T))

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/rda24_5x_n24_SG_B_m8_9_10_PC3_fishing")
Env.data <- read.csv("rda24_5x_n24_SG_B_m8_9_10_PC3_fishing.csv")

Env.data
# 'data.frame':	24 obs. of  7 variables:
#   $ Lat            : num  -64.8 -60.7 -60.7 -60.7 -60.7 -60.7 -64.8 -64.8 -64.8 -62.2 ...
# $ Long           : num  167 -27.8 -27.8 -27.8 -27.8 -27.8 167 167 167 76.6 ...
# $ Gen_ZONE_FG    : int  88 48 48 48 48 48 88 88 88 58 ...
# $ PC1            : num  -1.32 2.44 2.26 2.26 2.26 ...
# $ PC2            : num  -2.275 -0.739 -0.327 -0.327 -0.327 ...
# $ PC3            : num  0.038 -0.425 -0.244 -0.244 -0.244 ...
# $ fishing_ind_sum: int  398 0 0 0 0 0 398 398 398 356 ...

Env.ras <- data.frame(raster::extract(ras_all_1985, Env.data[, c("Lat", "Long")]))

Env.r  <- cor(Env.ras, method = "pearson")

final.env.names <- data.frame(ID = c("cm_198512_2_3_1","cz_198512_2_3_1","mld1_198512_2","mld3_198512_2","sal198512_2","sic198512_2","sit198512_2","sst_198512_2"),
                              new.ID = c("cm_198512","cz_198512","mld1_198512","mld3_198512","sal198512","sic198512","sit198512","sst_198512"))

dimnames(Env.r)[[1]] <- data.frame(ID = dimnames(Env.r)[[1]]) %>% 
  left_join(final.env.names) %>% pull(new.ID)

dimnames(Env.r)[[2]] <- data.frame(ID = dimnames(Env.r)[[2]]) %>% 
  left_join(final.env.names) %>% pull(new.ID)

gg.r.env <- ggcorrplot::ggcorrplot(Env.r,   hc.order = TRUE, type = "full", lab = TRUE) +
  theme_bw(base_size = 8) + theme(axis.title = element_blank(),
                                  strip.background = element_rect(fill="white"),
                                  legend.position = "bottom") +
  facet_wrap(~"Spearman's correlations")

gg.r.env 

# RDAdapt #####

# pacakges ####

library(devtools)
# install.packages("geometry")
library(geometry)
# install_github("landscape-genomics/rdadapt")
library(rdadapt)
library(vegan) ## rda

## Other libraries essentially used for plotting purposes
library(terra)
library(reshape2)
library(ggplot2)
library(viridis)
# install.packages("wesanderson")
library(wesanderson)

rdadapt_env <- rdadapt(prda_can_1, 2)
str(rdadapt_env)
# 'data.frame':	854 obs. of  2 variables:
# $ p.values: num  3.58e-01 9.07e-01 2.37e-35 4.58e-01 7.45e-05 ...
# $ q.values: num  4.91e-01 5.05e-01 8.86e-35 4.91e-01 2.10e-04 ...

res_RDA_proj_current <- adaptive_index(prda_can_1, 2, ras_all_1985, env_mask = NULL, method = "loadings")
# Error in (function (classes, fdef, mtable)  : 
            # unable to find an inherited method for function ‘adaptive_index’ for signature ‘"rda", "numeric", "RasterStack", "NULL"’

genomic_offset(prda_can_1, 2, ras_all_1985, ras_all_2022, env_mask = NULL, method = "loadings")
# Error in (function (classes, fdef, mtable)  : 
            # unable to find an inherited method for function ‘genomic_offset’ for signature ‘"rda", "numeric", "RasterStack", "RasterStack", "NULL", "missing"’

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