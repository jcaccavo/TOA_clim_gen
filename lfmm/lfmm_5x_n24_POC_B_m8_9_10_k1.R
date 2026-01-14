# LFMM (Univariate latent-factor linear mixed model) using lfmm package function lfmm_ridge
# using K = 1 based on snmf analysis k with lowest cross-entropy value
# POINT OF CAPTURE - BIRTH YEAR (POC-B) environmental data

# helpful links
# https://rdrr.io/cran/lfmm/man/lfmm_ridge.html

# run on server
# Run in R in directory

# in R
# install packages (if needed)
### BiocManager::install(version = "3.18")
### BiocManager::install("qvalue")

# load packages ####
library(vegan)
library(qvalue)
library(adegenet)
library(lfmm)

# load SNPs ####
snps_unlinked_5x <-read.PLINK(".../5x_TOA_only_filtered_SNPs_all_unlinked.raw")

# impute missing SNPs
snps_unlinked_imp_5x <- apply(snps_unlinked_5x, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

# current_m #############################################################################################################################
setwd(".../09_lfmm1_5x_n24_POC_B_m8_9_10_k1_cm")

# import environmental data
cm <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_current_m.csv", header=FALSE)

str(cm)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  -0.01235 0.0226 -0.00549 -0.00549 -0.00549 ...

# run LFMM
toa.lfmm.cm.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=cm, K=1)
# takes <10 seconds

# get p-values
toa.pv.cm.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=cm, lfmm=toa.lfmm.cm.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_cm_k1_unadjusted_p.pdf")
hist(toa.pv.cm.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_cm_k1_GIF-adjusted_p.pdf")
hist(toa.pv.cm.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.cm.k1 <- qvalue(toa.pv.cm.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.cm.k1 < 0.1))
# 4742

# which SNPs with a FDR > 10%
toa.FDR.cm.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.cm.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.cm.k1,"lfmm1_5x_n24_POC_B_m8_9_10_cm_k1_putative_adaptive_SNPs_0.1.csv")

# current_z #############################################################################################################################
setwd(".../10_lfmm1_5x_n24_POC_B_m8_9_10_k1_cz")

# import environmental data
cz <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_current_z.csv", header=FALSE)

str(cz)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  0.1284 0.0412 -0.0672 -0.0672 -0.0672 ...

# run LFMM
toa.lfmm.cz.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=cz, K=1)
# takes <10 seconds

# get p-values
toa.pv.cz.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=cz, lfmm=toa.lfmm.cz.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_cz_k1_unadjusted_p.pdf")
hist(toa.pv.cz.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_cz_k1_GIF-adjusted_p.pdf")
hist(toa.pv.cz.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.cz.k1 <- qvalue(toa.pv.cz.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.cz.k1 < 0.1))
# 1482

# which SNPs with a FDR > 10%
toa.FDR.cz.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.cz.k1 < 0.1)]

#Save adaptive SNPs
write.csv(toa.FDR.cz.k1,"lfmm1_5x_n24_POC_B_m8_9_10_cz_k1_putative_adaptive_SNPs_0.1.csv")

# mxl_0.01 #############################################################################################################################
setwd(".../11_lfmm1_5x_n24_POC_B_m8_9_10_k1_mxl1")

# import environmental data
mxl1 <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_mxl0.01.csv", header=FALSE)

str(mxl1)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  67.7 49.3 142.7 142.7 142.7 ...

# run LFMM
toa.lfmm.mxl1.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=mxl1, K=1)
# takes <10 seconds

# get p-values
toa.pv.mxl1.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=mxl1, lfmm=toa.lfmm.mxl1.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_mxl1_k1_unadjusted_p.pdf")
hist(toa.pv.mxl1.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_mxl1_k1_GIF-adjusted_p.pdf")
hist(toa.pv.mxl1.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.mxl1.k1 <- qvalue(toa.pv.mxl1.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.mxl1.k1 < 0.1))
# 12639

# which SNPs with a FDR > 10%
toa.FDR.mxl1.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.mxl1.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.mxl1.k1,"lfmm1_5x_n24_POC_B_m8_9_10_mxl1_k1_putative_adaptive_SNPs_0.1.csv")

# mxl_0.03 #############################################################################################################################
setwd(".../12_lfmm1_5x_n24_POC_B_m8_9_10_k1_mxl3")

# import environmental data
mxl3 <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_mxl0.03.csv", header=FALSE)

str(mxl3)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  75.8 51 185.2 185.2 185.2 ...

# run LFMM
toa.lfmm.mxl3.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=mxl3, K=1)
# takes <10 seconds

# get p-values
toa.pv.mxl3.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=mxl3, lfmm=toa.lfmm.mxl3.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_mxl3_k1_unadjusted_p.pdf")
hist(toa.pv.mxl3.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_mxl3_k1_GIF-adjusted_p.pdf")
hist(toa.pv.mxl3.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.mxl3.k1 <- qvalue(toa.pv.mxl3.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.mxl3.k1 < 0.1))
# 19735

# which SNPs with a FDR > 10%
toa.FDR.mxl3.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.mxl3.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.mxl3.k1,"lfmm1_5x_n24_POC_B_m8_9_10_mxl3_k1_putative_adaptive_SNPs_0.1.csv")

# salinity #############################################################################################################################
setwd(".../13_lfmm1_5x_n24_POC_B_m8_9_10_k1_sal")

# import environmental data
sal <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_sal.csv", header=FALSE)

str(sal)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  33.9 34.1 34.2 34.2 34.2 ...

# run LFMM
toa.lfmm.sal.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=sal, K=1)
# takes <10 seconds

# get p-values
toa.pv.sal.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=sal, lfmm=toa.lfmm.sal.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_sal_unadjusted_p.pdf")
hist(toa.pv.sal.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_sal_k1_GIF-adjusted_p.pdf")
hist(toa.pv.sal.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.sal.k1 <- qvalue(toa.pv.sal.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.sal.k1 < 0.1))
# 10300

# which SNPs with a FDR > 10%
toa.FDR.sal.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.sal.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.sal.k1,"lfmm1_5x_n24_POC_B_m8_9_10_sal_k1_putative_adaptive_SNPs_0.1.csv")

# sea_ice_concentration #############################################################################################################################
setwd(".../14_lfmm1_5x_n24_POC_B_m8_9_10_k1_sic")

# import environmental data
sic <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_sea_ice_conc.csv", header=FALSE)

str(sic)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  0.619 0.889 0.934 0.934 0.934 ...

# run LFMM
toa.lfmm.sic.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=sic, K=1)
# takes <10 seconds

# get p-values
toa.pv.sic.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=sic, lfmm=toa.lfmm.sic.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_sic_k1_unadjusted_p.pdf")
hist(toa.pv.sic.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_sic_k1_GIF-adjusted_p.pdf")
hist(toa.pv.sic.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.sic.k1 <- qvalue(toa.pv.sic.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.sic.k1 < 0.1))
# 3474

# which SNPs with a FDR > 10%
toa.FDR.sic.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.sic.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.sic.k1,"lfmm1_5x_n24_POC_B_m8_9_10_sic_k1_putative_adaptive_SNPs_0.1.csv")

# sea_ice_thickness #############################################################################################################################
setwd(".../15_lfmm1_5x_n24_POC_B_m8_9_10_k1_sit")

# import environmental data
sit <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_sea_ice_thick.csv", header=FALSE)

str(sit)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  0.549 1.298 1.121 1.121 1.121 ...

# run LFMM
toa.lfmm.sit.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=sit, K=1)
# takes <10 seconds

# get p-values
toa.pv.sit.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=sit, lfmm=toa.lfmm.sit.k1, calibrate="gif")
# takes <1 minute

# Create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_sit_k1_unadjusted_p.pdf")
hist(toa.pv.sit.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_sit_k1_GIF-adjusted_p.pdf")
hist(toa.pv.sit.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.sit.k1 <- qvalue(toa.pv.sit.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.sit.k1 < 0.1))
# 12902

# which SNPs with a FDR > 10%
toa.FDR.sit.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.sit.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.sit.k1,"lfmm1_5x_n24_POC_B_m8_9_10_sit_k1_putative_adaptive_SNPs_0.1.csv")

# sea_surface_temperature #############################################################################################################################
setwd(".../16_lfmm1_5x_n24_POC_B_m8_9_10_k1_sst")

# import environmental data
sst <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_sst.csv", header=FALSE)

str(sst)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  -1.69 -1.84 -1.82 -1.82 -1.82 ...

# run LFMM
toa.lfmm.sst.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=sst, K=1)
# takes <10 seconds

# get p-values
toa.pv.sst.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=sst, lfmm=toa.lfmm.sst.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_sst_k1_unadjusted_p.pdf")
hist(toa.pv.sst.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_POC_B_m8_9_10_sst_k1_GIF-adjusted_p.pdf")
hist(toa.pv.sst.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.sst.k1 <- qvalue(toa.pv.sst.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.sst.k1 < 0.1))
# 3018

# which SNPs with a FDR > 10%
toa.FDR.sst.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.sst.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.sst.k1,"lfmm1_5x_n24_POC_B_m8_9_10_sst_k1_putative_adaptive_SNPs_0.1.csv")

# PC1 #############################################################################################################################
setwd(".../29_lfmm1_5x_n24_POC_B_m8_9_10_k1_PC1")

# import environmental data
PC1 <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_PC1.csv", header=FALSE)

str(PC1)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  1.294 2.173 0.676 0.676 0.676 ...

# run LFMM
toa.lfmm.PC1.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=PC1, K=1)
# takes <10 seconds

# get p-values
toa.pv.PC1.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=PC1, lfmm=toa.lfmm.PC1.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm_5x_n24_POC_B_m8_9_10_PC1_k1_unadjusted_p.pdf")
hist(toa.pv.PC1.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm_5x_n24_POC_B_m8_9_10_PC1_k1_GIF-adjusted_p.pdf")
hist(toa.pv.PC1.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.PC1.k1 <- qvalue(toa.pv.PC1.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.PC1.k1 < 0.1))
# 2919

# which SNPs with a FDR > 10%
toa.FDR.PC1.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.PC1.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.PC1.k1,"lfmm_5x_n24_POC_B_m8_9_10_PC1_k1_putative_adaptive_SNPs_0.1.csv")

# PC2 #############################################################################################################################
setwd(".../30_lfmm1_5x_n24_POC_B_m8_9_10_k1_PC2")

# import environmental data
PC2 <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_PC2.csv", header=FALSE)

str(PC2)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  -3.114 -0.269 0.912 0.912 0.912 ...

# run LFMM
toa.lfmm.PC2.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=PC2, K=1)
# takes <10 seconds

# get p-values
toa.pv.PC2.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=PC2, lfmm=toa.lfmm.PC2.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm_5x_n24_POC_B_m8_9_10_PC2_k1_unadjusted_p.pdf")
hist(toa.pv.PC2.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm_5x_n24_POC_B_m8_9_10_PC2_k1_GIF-adjusted_p.pdf")
hist(toa.pv.PC2.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.PC2.k1 <- qvalue(toa.pv.PC2.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.PC2.k1 < 0.1))
# 1747

# which SNPs with a FDR > 10%
toa.FDR.PC2.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.PC2.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.PC2.k1,"lfmm_5x_n24_POC_B_m8_9_10_PC2_k1_putative_adaptive_SNPs_0.1.csv")

# PC3 #############################################################################################################################
setwd(".../31_lfmm1_5x_n24_POC_B_m8_9_10_k1_PC3")

# import environmental data
PC3 <- read.csv(".../lfmm_5x_n24_POC_B_m8_9_10_PC3.csv", header=FALSE)

str(PC3)
# 'data.frame':   24 obs. of  1 variable:
#   $ V1: num  -1.843 0.379 -0.211 -0.211 -0.211 ...

# run LFMM
toa.lfmm.PC3.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=PC3, K=1)
# takes <10 seconds

# get p-values
toa.pv.PC3.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=PC3, lfmm=toa.lfmm.PC3.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm_5x_n24_POC_B_m8_9_10_PC3_k1_unadjusted_p.pdf")
hist(toa.pv.PC3.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm_5x_n24_POC_B_m8_9_10_PC3_k1_GIF-adjusted_p.pdf")
hist(toa.pv.PC3.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.PC3.k1 <- qvalue(toa.pv.PC3.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.PC3.k1 < 0.1))
# 4587

# which SNPs with a FDR > 10%
toa.FDR.PC3.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.PC3.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.PC3.k1,"lfmm_5x_n24_POC_B_m8_9_10_PC3_k1_putative_adaptive_SNPs_0.1.csv")
