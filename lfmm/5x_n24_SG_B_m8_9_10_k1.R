# LFMM (Univariate latent-factor linear mixed model) using lfmm package function lfmm_ridge
# using K = 1 based on snmf analysis k with lowest cross-entropy value
# SPAWNING GROUND - BIRTH YEAR (SG-B) environmental data

# helpful links
# https://rdrr.io/cran/lfmm/man/lfmm_ridge.html

# run on server (yangtze)
# Run in R in directory /srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/#RELEVANT RUN FOLDER

# in R
# install packages (if needed)
### BiocManager::install(version = "3.18")
### BiocManager::install("qvalue")

# load packages
library(vegan)
library(qvalue)
library(adegenet)
library(lfmm)

# load SNPs
snps_unlinked_5x <-read.PLINK("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/03_gen_data/5x_TOA_only_filtered_SNPs_all_unlinked.raw")
# impute missing SNPs
snps_unlinked_imp_5x <- apply(snps_unlinked_5x, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
# may take a few minutes

# current_m #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/01_lfmm1_5x_n24_SG_B_m8_9_10_k1_cm")

# import environmental data
cm <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_current_m.csv")

# run LFMM
toa.lfmm.cm.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=cm, K=1)
# takes <10 seconds

# get p-values
toa.pv.cm.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=cm, lfmm=toa.lfmm.cm.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_cm_k1_unadjusted_p.pdf")
hist(toa.pv.cm.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_cm_k1_GIF-adjusted_p.pdf")
hist(toa.pv.cm.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.cm.k1 <- qvalue(toa.pv.cm.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.cm.k1 < 0.1))
# 8109

# which SNPs with a FDR > 10%
toa.FDR.cm.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.cm.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.cm.k1,"lfmm1_5x_n24_SG_B_m8_9_10_cm_k1_putative_adaptive_SNPs_0.1.csv")

# current_z #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/02_lfmm1_5x_n24_SG_B_m8_9_10_k1_cz")

# import environmental data
cz <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_current_z.csv")

# run LFMM
toa.lfmm.cz.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=cz, K=1)
# takes <10 seconds

# get p-values
toa.pv.cz.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=cz, lfmm=toa.lfmm.cz.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_cz_k1_unadjusted_p.pdf")
hist(toa.pv.cz.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_cz_k1_GIF-adjusted_p.pdf")
hist(toa.pv.cz.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.cz.k1 <- qvalue(toa.pv.cz.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.cz.k1 < 0.1))
# 8340

# which SNPs with a FDR > 10%
toa.FDR.cz.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.cz.k1 < 0.1)]

#Save adaptive SNPs
write.csv(toa.FDR.cz.k1,"lfmm1_5x_n24_SG_B_m8_9_10_cz_k1_putative_adaptive_SNPs_0.1.csv")

# mxl_0.01 #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/03_lfmm1_5x_n24_SG_B_m8_9_10_k1_mxl1")

# import environmental data
mxl1 <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_mxl0.01.csv")

# run LFMM
toa.lfmm.mxl1.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=mxl1, K=1)
# takes <10 seconds

# get p-values
toa.pv.mxl1.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=mxl1, lfmm=toa.lfmm.mxl1.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_mxl1_k1_unadjusted_p.pdf")
hist(toa.pv.mxl1.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_mxl1_k1_GIF-adjusted_p.pdf")
hist(toa.pv.mxl1.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.mxl1.k1 <- qvalue(toa.pv.mxl1.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.mxl1.k1 < 0.1))
# 6854

# which SNPs with a FDR > 10%
toa.FDR.mxl1.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.mxl1.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.mxl1.k1,"lfmm1_5x_n24_SG_B_m8_9_10_mxl1_k1_putative_adaptive_SNPs_0.1.csv")

# mxl_0.03 #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/04_lfmm1_5x_n24_SG_B_m8_9_10_k1_mxl3")

# import environmental data
mxl3 <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_mxl0.03.csv")

# run LFMM
toa.lfmm.mxl3.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=mxl3, K=1)
# takes <10 seconds

# get p-values
toa.pv.mxl3.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=mxl3, lfmm=toa.lfmm.mxl3.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_mxl3_k1_unadjusted_p.pdf")
hist(toa.pv.mxl3.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_mxl3_k1_GIF-adjusted_p.pdf")
hist(toa.pv.mxl3.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.mxl3.k1 <- qvalue(toa.pv.mxl3.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.mxl3.k1 < 0.1))
# 6922

# which SNPs with a FDR > 10%
toa.FDR.mxl3.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.mxl3.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.mxl3.k1,"lfmm1_5x_n24_SG_B_m8_9_10_mxl3_k1_putative_adaptive_SNPs_0.1.csv")

# salinity #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/05_lfmm1_5x_n24_SG_B_m8_9_10_k1_sal")

# import environmental data
sal <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_sal.csv")

# run LFMM
toa.lfmm.sal.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=sal, K=1)
# takes <10 seconds

# get p-values
toa.pv.sal.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=sal, lfmm=toa.lfmm.sal.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_sal_unadjusted_p.pdf")
hist(toa.pv.sal.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_sal_k1_GIF-adjusted_p.pdf")
hist(toa.pv.sal.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.sal.k1 <- qvalue(toa.pv.sal.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.sal.k1 < 0.1))
# 10436

# which SNPs with a FDR > 10%
toa.FDR.sal.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.sal.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.sal.k1,"lfmm1_5x_n24_SG_B_m8_9_10_sal_k1_putative_adaptive_SNPs_0.1.csv")

# sea_ice_concentration #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/06_lfmm1_5x_n24_SG_B_m8_9_10_k1_sic")

# import environmental data
sic <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_sea_ice_conc.csv")

# run LFMM
toa.lfmm.sic.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=sic, K=1)
# takes <10 seconds

# get p-values
toa.pv.sic.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=sic, lfmm=toa.lfmm.sic.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_sic_k1_unadjusted_p.pdf")
hist(toa.pv.sic.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_sic_k1_GIF-adjusted_p.pdf")
hist(toa.pv.sic.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.sic.k1 <- qvalue(toa.pv.sic.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.sic.k1 < 0.1))
# 428

# which SNPs with a FDR > 10%
toa.FDR.sic.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.sic.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.sic.k1,"lfmm1_5x_n24_SG_B_m8_9_10_sic_k1_putative_adaptive_SNPs_0.1.csv")

# sea_ice_thickness #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/07_lfmm1_5x_n24_SG_B_m8_9_10_k1_sit")

# import environmental data
sit <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_sea_ice_thick.csv")

# run LFMM
toa.lfmm.sit.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=sit, K=1)
# takes <10 seconds

# get p-values
toa.pv.sit.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=sit, lfmm=toa.lfmm.sit.k1, calibrate="gif")
# takes <1 minute

# Create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_sit_k1_unadjusted_p.pdf")
hist(toa.pv.sit.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_sit_k1_GIF-adjusted_p.pdf")
hist(toa.pv.sit.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.sit.k1 <- qvalue(toa.pv.sit.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.sit.k1 < 0.1))
# 474

# which SNPs with a FDR > 10%
toa.FDR.sit.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.sit.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.sit.k1,"lfmm1_5x_n24_SG_B_m8_9_10_sit_k1_putative_adaptive_SNPs_0.1.csv")

# sea_surface_temperature #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/08_lfmm1_5x_n24_SG_B_m8_9_10_k1_sst")

# import environmental data
sst <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_sst.csv")

# run LFMM
toa.lfmm.sst.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=sst, K=1)
# takes <10 seconds

# get p-values
toa.pv.sst.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=sst, lfmm=toa.lfmm.sst.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_sst_k1_unadjusted_p.pdf")
hist(toa.pv.sst.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_sst_k1_GIF-adjusted_p.pdf")
hist(toa.pv.sst.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.sst.k1 <- qvalue(toa.pv.sst.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.sst.k1 < 0.1))
# 209

# which SNPs with a FDR > 10%
toa.FDR.sst.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.sst.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.sst.k1,"lfmm1_5x_n24_SG_B_m8_9_10_sst_k1_putative_adaptive_SNPs_0.1.csv")

# PC1 #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/25_lfmm1_5x_n24_SG_B_m8_9_10_k1_PC1")

# import environmental data
PC1 <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_PC1.csv", header=FALSE)

str(PC1)
# 'data.frame':	24 obs. of  1 variable:
# $ V1: num  -1.32 2.44 2.26 2.26 2.26 ...

# run LFMM
toa.lfmm.PC1.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=PC1, K=1)
# takes <10 seconds

# get p-values
toa.pv.PC1.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=PC1, lfmm=toa.lfmm.PC1.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_PC1_k1_unadjusted_p.pdf")
hist(toa.pv.PC1.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_PC1_k1_GIF-adjusted_p.pdf")
hist(toa.pv.PC1.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.PC1.k1 <- qvalue(toa.pv.PC1.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.PC1.k1 < 0.1))
# 12562

# which SNPs with a FDR > 10%
toa.FDR.PC1.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.PC1.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.PC1.k1,"lfmm1_5x_n24_SG_B_m8_9_10_PC1_k1_putative_adaptive_SNPs_0.1.csv")

# PC2 #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/26_lfmm1_5x_n24_SG_B_m8_9_10_k1_PC2")

# import environmental data
PC2 <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_PC2.csv", header=FALSE)

# run LFMM
toa.lfmm.PC2.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=PC2, K=1)
# takes <10 seconds

# get p-values
toa.pv.PC2.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=PC2, lfmm=toa.lfmm.PC2.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_PC2_k1_unadjusted_p.pdf")
hist(toa.pv.PC2.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_PC2_k1_GIF-adjusted_p.pdf")
hist(toa.pv.PC2.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.PC2.k1 <- qvalue(toa.pv.PC2.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.PC2.k1 < 0.1))
# 716

# which SNPs with a FDR > 10%
toa.FDR.PC2.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.PC2.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.PC2.k1,"lfmm1_5x_n24_SG_B_m8_9_10_PC2_k1_putative_adaptive_SNPs_0.1.csv")

# PC3 #############################################################################################################################
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/27_lfmm1_5x_n24_SG_B_m8_9_10_k1_PC3")

# import environmental data
PC3 <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/lfmm_5x_n24_SG_B_m8_9_10_PC3.csv", header=FALSE)

# run LFMM
toa.lfmm.PC3.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=PC3, K=1)
# takes <10 seconds

# get p-values
toa.pv.PC3.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=PC3, lfmm=toa.lfmm.PC3.k1, calibrate="gif")
# takes <1 minute

# create plots
# unadjusted p-values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_PC3_k1_unadjusted_p.pdf")
hist(toa.pv.PC3.k1$pvalue[,1], main="Unadjusted p-values")
dev.off()
# GIF-adjusted p values
pdf("lfmm1_5x_n24_SG_B_m8_9_10_PC3_k1_GIF-adjusted_p.pdf")
hist(toa.pv.PC3.k1$calibrated.pvalue[,1], main="GIF-adjusted p-values")
dev.off()

# calculate LFMM statistics
toa.qv.PC3.k1 <- qvalue(toa.pv.PC3.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.PC3.k1 < 0.1))
# 130088

# which SNPs with a FDR > 10%
toa.FDR.PC3.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.PC3.k1 < 0.1)]

# save adaptive SNPs
write.csv(toa.FDR.PC3.k1,"lfmm1_5x_n24_SG_B_m8_9_10_PC3_k1_putative_adaptive_SNPs_0.1.csv")

# RANDOM VALUE RUNS FOR FALSE POSITIVE ANALYSIS ####
setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/28_lfmm1_5x_n24_SG_B_m8_9_10_k1_rand")
# import random values 1 -3
envs <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/rda25_5x_n24_rand3_PCs.csv")
# import random values 4 -100
envs <- read.csv("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/02_env_data/pRDA_rand_5x_n24_4-100.csv")

# rand1 #############################################################################################################################

# subset rand1
retained_preds <- subset(envs, select=c(rand1))

# run LFMM
toa.lfmm.rand1.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand1.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand1.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand1.k1 <- qvalue(toa.pv.rand1.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand1.k1 < 0.1))
# 25378

# which SNPs with a FDR > 10%
toa.FDR.rand1.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand1.k1 < 0.1)]

# rand2 #############################################################################################################################

# subset rand2
retained_preds <- subset(envs, select=c(rand2))

# run LFMM
toa.lfmm.rand2.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand2.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand2.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand2.k1 <- qvalue(toa.pv.rand2.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand2.k1 < 0.1))
# 224

# which SNPs with a FDR > 10%
toa.FDR.rand2.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand2.k1 < 0.1)]

# rand3 #############################################################################################################################

# subset rand3
retained_preds <- subset(envs, select=c(rand3))

# run LFMM
toa.lfmm.rand3.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand3.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand3.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand3.k1 <- qvalue(toa.pv.rand3.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand3.k1 < 0.1))
# 5028

# which SNPs with a FDR > 10%
toa.FDR.rand3.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand3.k1 < 0.1)]

# combine outputs of candidate SNPs from random value LFMM runs 1 -3
max_length <- max(length(toa.FDR.rand1.k1), length(toa.FDR.rand2.k1), length(toa.FDR.rand3.k1))
length(toa.FDR.rand1.k1) <- max_length
length(toa.FDR.rand2.k1) <- max_length
length(toa.FDR.rand3.k1) <- max_length
rand1_3 <- cbind(toa.FDR.rand1.k1, toa.FDR.rand2.k1, toa.FDR.rand3.k1)

# save false candidate SNPs from random value LFMM runs 1-3 to .csv file
write.csv(rand1_3,"lfmm1_5x_n24_SG_B_m8_9_10_rand1-3_k1_putative_adaptive_SNPs_0.1.csv")

# rand4 #############################################################################################################################

# subset rand4
retained_preds <- subset(envs, select=c(rand4))

# run LFMM
toa.lfmm.rand4.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand4.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand4.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand4.k1 <- qvalue(toa.pv.rand4.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand4.k1 < 0.1))
# 441

# which SNPs with a FDR > 10%
toa.FDR.rand4.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand4.k1 < 0.1)]

# rand5 #############################################################################################################################

# subset rand5
retained_preds <- subset(envs, select=c(rand5))

# run LFMM
toa.lfmm.rand5.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand5.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand5.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand5.k1 <- qvalue(toa.pv.rand5.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand5.k1 < 0.1))
# 3582

# which SNPs with a FDR > 10%
toa.FDR.rand5.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand5.k1 < 0.1)]

# rand6 #############################################################################################################################

# subset rand6
retained_preds <- subset(envs, select=c(rand6))

# run LFMM
toa.lfmm.rand6.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand6.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand6.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand6.k1 <- qvalue(toa.pv.rand6.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand6.k1 < 0.1))
# 12040

# which SNPs with a FDR > 10%
toa.FDR.rand6.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand6.k1 < 0.1)]

# rand7 #############################################################################################################################

# subset rand7
retained_preds <- subset(envs, select=c(rand7))

# run LFMM
toa.lfmm.rand7.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand7.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand7.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand7.k1 <- qvalue(toa.pv.rand7.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand7.k1 < 0.1))
# 247

# which SNPs with a FDR > 10%
toa.FDR.rand7.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand7.k1 < 0.1)]

# rand8 #############################################################################################################################

# subset rand8
retained_preds <- subset(envs, select=c(rand8))

# run LFMM
toa.lfmm.rand8.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand8.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand8.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand8.k1 <- qvalue(toa.pv.rand8.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand8.k1 < 0.1))
# 102

# which SNPs with a FDR > 10%
toa.FDR.rand8.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand8.k1 < 0.1)]

# rand9 #############################################################################################################################

# subset rand9
retained_preds <- subset(envs, select=c(rand9))

# run LFMM
toa.lfmm.rand9.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand9.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand9.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand9.k1 <- qvalue(toa.pv.rand9.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand9.k1 < 0.1))
# 167

# which SNPs with a FDR > 10%
toa.FDR.rand9.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand9.k1 < 0.1)]

# rand10 #############################################################################################################################

# subset rand10
retained_preds <- subset(envs, select=c(rand10))

# run LFMM
toa.lfmm.rand10.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand10.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand10.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand10.k1 <- qvalue(toa.pv.rand10.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand10.k1 < 0.1))
# 135

# which SNPs with a FDR > 10%
toa.FDR.rand10.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand10.k1 < 0.1)]

# rand11 #############################################################################################################################

# subset rand11
retained_preds <- subset(envs, select=c(rand11))

# run LFMM
toa.lfmm.rand11.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand11.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand11.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand11.k1 <- qvalue(toa.pv.rand11.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand11.k1 < 0.1))
# 16664

# which SNPs with a FDR > 10%
toa.FDR.rand11.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand11.k1 < 0.1)]

# rand12 #############################################################################################################################

# subset rand12
retained_preds <- subset(envs, select=c(rand12))

# run LFMM
toa.lfmm.rand12.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand12.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand12.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand12.k1 <- qvalue(toa.pv.rand12.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand12.k1 < 0.1))
# 1102

# which SNPs with a FDR > 10%
toa.FDR.rand12.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand12.k1 < 0.1)]

# rand13 #############################################################################################################################

# subset rand13
retained_preds <- subset(envs, select=c(rand13))

# run LFMM
toa.lfmm.rand13.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand13.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand13.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand13.k1 <- qvalue(toa.pv.rand13.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand13.k1 < 0.1))
# 176

# which SNPs with a FDR > 10%
toa.FDR.rand13.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand13.k1 < 0.1)]

# rand14 #############################################################################################################################

# subset rand14
retained_preds <- subset(envs, select=c(rand14))

# run LFMM
toa.lfmm.rand14.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand14.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand14.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand14.k1 <- qvalue(toa.pv.rand14.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand14.k1 < 0.1))
# 21654

# which SNPs with a FDR > 10%
toa.FDR.rand14.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand14.k1 < 0.1)]

# rand15 #############################################################################################################################

# subset rand15
retained_preds <- subset(envs, select=c(rand15))

# run LFMM
toa.lfmm.rand15.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand15.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand15.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand15.k1 <- qvalue(toa.pv.rand15.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand15.k1 < 0.1))
# 6948

# which SNPs with a FDR > 10%
toa.FDR.rand15.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand15.k1 < 0.1)]

# rand16 #############################################################################################################################

# subset rand16
retained_preds <- subset(envs, select=c(rand16))

# run LFMM
toa.lfmm.rand16.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand16.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand16.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand16.k1 <- qvalue(toa.pv.rand16.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand16.k1 < 0.1))
# 25378

# which SNPs with a FDR > 10%
toa.FDR.rand16.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand16.k1 < 0.1)]

# rand17 #############################################################################################################################

# subset rand17
retained_preds <- subset(envs, select=c(rand17))

# run LFMM
toa.lfmm.rand17.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand17.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand17.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand17.k1 <- qvalue(toa.pv.rand17.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand17.k1 < 0.1))
# 224

# which SNPs with a FDR > 10%
toa.FDR.rand17.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand17.k1 < 0.1)]

# rand18 #############################################################################################################################

# subset rand18
retained_preds <- subset(envs, select=c(rand18))

# run LFMM
toa.lfmm.rand18.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand18.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand18.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand18.k1 <- qvalue(toa.pv.rand18.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand18.k1 < 0.1))
# 5028

# which SNPs with a FDR > 10%
toa.FDR.rand18.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand18.k1 < 0.1)]

# rand19 #############################################################################################################################

# subset rand19
retained_preds <- subset(envs, select=c(rand19))

# run LFMM
toa.lfmm.rand19.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand19.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand19.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand19.k1 <- qvalue(toa.pv.rand19.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand19.k1 < 0.1))
# 1157

# which SNPs with a FDR > 10%
toa.FDR.rand19.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand19.k1 < 0.1)]

# rand20 #############################################################################################################################

# subset rand20
retained_preds <- subset(envs, select=c(rand20))

# run LFMM
toa.lfmm.rand20.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand20.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand20.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand20.k1 <- qvalue(toa.pv.rand20.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand20.k1 < 0.1))
# 22003

# which SNPs with a FDR > 10%
toa.FDR.rand20.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand20.k1 < 0.1)]

# rand21 #############################################################################################################################

# subset rand21
retained_preds <- subset(envs, select=c(rand21))

# run LFMM
toa.lfmm.rand21.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand21.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand21.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand21.k1 <- qvalue(toa.pv.rand21.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand21.k1 < 0.1))
# 308

# which SNPs with a FDR > 10%
toa.FDR.rand21.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand21.k1 < 0.1)]

# rand22 #############################################################################################################################

# subset rand22
retained_preds <- subset(envs, select=c(rand22))

# run LFMM
toa.lfmm.rand22.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand22.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand22.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand22.k1 <- qvalue(toa.pv.rand22.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand22.k1 < 0.1))
# 12734

# which SNPs with a FDR > 10%
toa.FDR.rand22.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand22.k1 < 0.1)]

# rand23 #############################################################################################################################

# subset rand23
retained_preds <- subset(envs, select=c(rand23))

# run LFMM
toa.lfmm.rand23.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand23.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand23.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand23.k1 <- qvalue(toa.pv.rand23.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand23.k1 < 0.1))
# 219

# which SNPs with a FDR > 10%
toa.FDR.rand23.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand23.k1 < 0.1)]

# rand24 #############################################################################################################################

# subset rand24
retained_preds <- subset(envs, select=c(rand24))

# run LFMM
toa.lfmm.rand24.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand24.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand24.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand24.k1 <- qvalue(toa.pv.rand24.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand24.k1 < 0.1))
# 21044

# which SNPs with a FDR > 10%
toa.FDR.rand24.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand24.k1 < 0.1)]

# rand25 #############################################################################################################################

# subset rand25
retained_preds <- subset(envs, select=c(rand25))

# run LFMM
toa.lfmm.rand25.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand25.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand25.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand25.k1 <- qvalue(toa.pv.rand25.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand25.k1 < 0.1))
# 967

# which SNPs with a FDR > 10%
toa.FDR.rand25.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand25.k1 < 0.1)]

# rand26 #############################################################################################################################

# subset rand26
retained_preds <- subset(envs, select=c(rand26))

# run LFMM
toa.lfmm.rand26.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand26.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand26.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand26.k1 <- qvalue(toa.pv.rand26.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand26.k1 < 0.1))
# 92

# which SNPs with a FDR > 10%
toa.FDR.rand26.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand26.k1 < 0.1)]

# rand27 #############################################################################################################################

# subset rand27
retained_preds <- subset(envs, select=c(rand27))

# run LFMM
toa.lfmm.rand27.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand27.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand27.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand27.k1 <- qvalue(toa.pv.rand27.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand27.k1 < 0.1))
# 232

# which SNPs with a FDR > 10%
toa.FDR.rand27.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand27.k1 < 0.1)]

# rand28 #############################################################################################################################

# subset rand28
retained_preds <- subset(envs, select=c(rand28))

# run LFMM
toa.lfmm.rand28.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand28.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand28.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand28.k1 <- qvalue(toa.pv.rand28.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand28.k1 < 0.1))
# 591

# which SNPs with a FDR > 10%
toa.FDR.rand28.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand28.k1 < 0.1)]

# rand29 #############################################################################################################################

# subset rand29
retained_preds <- subset(envs, select=c(rand29))

# run LFMM
toa.lfmm.rand29.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand29.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand29.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand29.k1 <- qvalue(toa.pv.rand29.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand29.k1 < 0.1))
# 927

# which SNPs with a FDR > 10%
toa.FDR.rand29.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand29.k1 < 0.1)]

# rand30 #############################################################################################################################

# subset rand30
retained_preds <- subset(envs, select=c(rand30))

# run LFMM
toa.lfmm.rand30.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand30.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand30.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand30.k1 <- qvalue(toa.pv.rand30.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand30.k1 < 0.1))
# 228

# which SNPs with a FDR > 10%
toa.FDR.rand30.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand30.k1 < 0.1)]

# rand31 #############################################################################################################################

# subset rand31
retained_preds <- subset(envs, select=c(rand31))

# run LFMM
toa.lfmm.rand31.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand31.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand31.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand31.k1 <- qvalue(toa.pv.rand31.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand31.k1 < 0.1))
# 1148

# which SNPs with a FDR > 10%
toa.FDR.rand31.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand31.k1 < 0.1)]

# rand32 #############################################################################################################################

# subset rand32
retained_preds <- subset(envs, select=c(rand32))

# run LFMM
toa.lfmm.rand32.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand32.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand32.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand32.k1 <- qvalue(toa.pv.rand32.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand32.k1 < 0.1))
# 922

# which SNPs with a FDR > 10%
toa.FDR.rand32.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand32.k1 < 0.1)]

# rand33 #############################################################################################################################

# subset rand33
retained_preds <- subset(envs, select=c(rand33))

# run LFMM
toa.lfmm.rand33.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand33.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand33.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand33.k1 <- qvalue(toa.pv.rand33.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand33.k1 < 0.1))
# 178

# which SNPs with a FDR > 10%
toa.FDR.rand33.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand33.k1 < 0.1)]

# rand34 #############################################################################################################################

# subset rand34
retained_preds <- subset(envs, select=c(rand34))

# run LFMM
toa.lfmm.rand34.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand34.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand34.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand34.k1 <- qvalue(toa.pv.rand34.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand34.k1 < 0.1))
# 991

# which SNPs with a FDR > 10%
toa.FDR.rand34.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand34.k1 < 0.1)]

# rand35 #############################################################################################################################

# subset rand35
retained_preds <- subset(envs, select=c(rand35))

# run LFMM
toa.lfmm.rand35.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand35.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand35.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand35.k1 <- qvalue(toa.pv.rand35.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand35.k1 < 0.1))
# 7146

# which SNPs with a FDR > 10%
toa.FDR.rand35.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand35.k1 < 0.1)]

# rand36 #############################################################################################################################

# subset rand36
retained_preds <- subset(envs, select=c(rand36))

# run LFMM
toa.lfmm.rand36.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand36.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand36.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand36.k1 <- qvalue(toa.pv.rand36.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand36.k1 < 0.1))
# 360

# which SNPs with a FDR > 10%
toa.FDR.rand36.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand36.k1 < 0.1)]

# rand37 #############################################################################################################################

# subset rand37
retained_preds <- subset(envs, select=c(rand37))

# run LFMM
toa.lfmm.rand37.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand37.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand37.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand37.k1 <- qvalue(toa.pv.rand37.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand37.k1 < 0.1))
# 937

# which SNPs with a FDR > 10%
toa.FDR.rand37.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand37.k1 < 0.1)]

# rand38 #############################################################################################################################

# subset rand38
retained_preds <- subset(envs, select=c(rand38))

# run LFMM
toa.lfmm.rand38.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand38.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand38.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand38.k1 <- qvalue(toa.pv.rand38.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand38.k1 < 0.1))
# 27168

# which SNPs with a FDR > 10%
toa.FDR.rand38.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand38.k1 < 0.1)]

# rand39 #############################################################################################################################

# subset rand39
retained_preds <- subset(envs, select=c(rand39))

# run LFMM
toa.lfmm.rand39.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand39.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand39.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand39.k1 <- qvalue(toa.pv.rand39.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand39.k1 < 0.1))
# 424

# which SNPs with a FDR > 10%
toa.FDR.rand39.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand39.k1 < 0.1)]

# rand40 #############################################################################################################################

# subset rand40
retained_preds <- subset(envs, select=c(rand40))

# run LFMM
toa.lfmm.rand40.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand40.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand40.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand40.k1 <- qvalue(toa.pv.rand40.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand40.k1 < 0.1))
# 2394

# which SNPs with a FDR > 10%
toa.FDR.rand40.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand40.k1 < 0.1)]

# rand41 #############################################################################################################################

# subset rand41
retained_preds <- subset(envs, select=c(rand41))

# run LFMM
toa.lfmm.rand41.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand41.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand41.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand41.k1 <- qvalue(toa.pv.rand41.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand41.k1 < 0.1))
# 30

# which SNPs with a FDR > 10%
toa.FDR.rand41.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand41.k1 < 0.1)]

# rand42 #############################################################################################################################

# subset rand42
retained_preds <- subset(envs, select=c(rand42))

# run LFMM
toa.lfmm.rand42.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand42.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand42.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand42.k1 <- qvalue(toa.pv.rand42.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand42.k1 < 0.1))
# 8044

# which SNPs with a FDR > 10%
toa.FDR.rand42.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand42.k1 < 0.1)]

# rand43 #############################################################################################################################

# subset rand43
retained_preds <- subset(envs, select=c(rand43))

# run LFMM
toa.lfmm.rand43.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand43.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand43.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand43.k1 <- qvalue(toa.pv.rand43.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand43.k1 < 0.1))
# 1600

# which SNPs with a FDR > 10%
toa.FDR.rand43.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand43.k1 < 0.1)]

# rand44 #############################################################################################################################

# subset rand44
retained_preds <- subset(envs, select=c(rand44))

# run LFMM
toa.lfmm.rand44.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand44.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand44.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand44.k1 <- qvalue(toa.pv.rand44.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand44.k1 < 0.1))
# 3198

# which SNPs with a FDR > 10%
toa.FDR.rand44.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand44.k1 < 0.1)]

# rand45 #############################################################################################################################

# subset rand45
retained_preds <- subset(envs, select=c(rand45))

# run LFMM
toa.lfmm.rand45.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand45.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand45.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand45.k1 <- qvalue(toa.pv.rand45.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand45.k1 < 0.1))
# 14263

# which SNPs with a FDR > 10%
toa.FDR.rand45.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand45.k1 < 0.1)]

# rand46 #############################################################################################################################

# subset rand46
retained_preds <- subset(envs, select=c(rand46))

# run LFMM
toa.lfmm.rand46.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand46.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand46.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand46.k1 <- qvalue(toa.pv.rand46.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand46.k1 < 0.1))
# 19270

# which SNPs with a FDR > 10%
toa.FDR.rand46.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand46.k1 < 0.1)]

# rand47 #############################################################################################################################

# subset rand47
retained_preds <- subset(envs, select=c(rand47))

# run LFMM
toa.lfmm.rand47.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand47.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand47.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand47.k1 <- qvalue(toa.pv.rand47.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand47.k1 < 0.1))
# 5781

# which SNPs with a FDR > 10%
toa.FDR.rand47.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand47.k1 < 0.1)]

# rand48 #############################################################################################################################

# subset rand48
retained_preds <- subset(envs, select=c(rand48))

# run LFMM
toa.lfmm.rand48.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand48.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand48.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand48.k1 <- qvalue(toa.pv.rand48.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand48.k1 < 0.1))
# 63

# which SNPs with a FDR > 10%
toa.FDR.rand48.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand48.k1 < 0.1)]

# rand49 #############################################################################################################################

# subset rand49
retained_preds <- subset(envs, select=c(rand49))

# run LFMM
toa.lfmm.rand49.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand49.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand49.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand49.k1 <- qvalue(toa.pv.rand49.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand49.k1 < 0.1))
# 287

# which SNPs with a FDR > 10%
toa.FDR.rand49.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand49.k1 < 0.1)]

# rand50 #############################################################################################################################

# subset rand50
retained_preds <- subset(envs, select=c(rand50))

# run LFMM
toa.lfmm.rand50.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand50.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand50.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand50.k1 <- qvalue(toa.pv.rand50.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand50.k1 < 0.1))
# 5134

# which SNPs with a FDR > 10%
toa.FDR.rand50.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand50.k1 < 0.1)]

# rand51 #############################################################################################################################

# subset rand51
retained_preds <- subset(envs, select=c(rand51))

# run LFMM
toa.lfmm.rand51.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand51.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand51.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand51.k1 <- qvalue(toa.pv.rand51.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand51.k1 < 0.1))
# 983

# which SNPs with a FDR > 10%
toa.FDR.rand51.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand51.k1 < 0.1)]

# rand52 #############################################################################################################################

# subset rand52
retained_preds <- subset(envs, select=c(rand52))

# run LFMM
toa.lfmm.rand52.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand52.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand52.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand52.k1 <- qvalue(toa.pv.rand52.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand52.k1 < 0.1))
# 114

# which SNPs with a FDR > 10%
toa.FDR.rand52.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand52.k1 < 0.1)]

# rand53 #############################################################################################################################

# subset rand53
retained_preds <- subset(envs, select=c(rand53))

# run LFMM
toa.lfmm.rand53.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand53.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand53.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand53.k1 <- qvalue(toa.pv.rand53.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand53.k1 < 0.1))
# 205

# which SNPs with a FDR > 10%
toa.FDR.rand53.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand53.k1 < 0.1)]

# rand54 #############################################################################################################################

# subset rand54
retained_preds <- subset(envs, select=c(rand54))

# run LFMM
toa.lfmm.rand54.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand54.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand54.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand54.k1 <- qvalue(toa.pv.rand54.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand54.k1 < 0.1))
# 567

# which SNPs with a FDR > 10%
toa.FDR.rand54.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand54.k1 < 0.1)]

# rand55 #############################################################################################################################

# subset rand55
retained_preds <- subset(envs, select=c(rand55))

# run LFMM
toa.lfmm.rand55.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand55.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand55.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand55.k1 <- qvalue(toa.pv.rand55.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand55.k1 < 0.1))
# 2618

# which SNPs with a FDR > 10%
toa.FDR.rand55.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand55.k1 < 0.1)]

# rand56 #############################################################################################################################

# subset rand56
retained_preds <- subset(envs, select=c(rand56))

# run LFMM
toa.lfmm.rand56.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand56.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand56.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand56.k1 <- qvalue(toa.pv.rand56.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand56.k1 < 0.1))
# 2599

# which SNPs with a FDR > 10%
toa.FDR.rand56.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand56.k1 < 0.1)]

# rand57 #############################################################################################################################

# subset rand57
retained_preds <- subset(envs, select=c(rand57))

# run LFMM
toa.lfmm.rand57.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand57.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand57.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand57.k1 <- qvalue(toa.pv.rand57.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand57.k1 < 0.1))
# 226

# which SNPs with a FDR > 10%
toa.FDR.rand57.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand57.k1 < 0.1)]

# rand58 #############################################################################################################################

# subset rand58
retained_preds <- subset(envs, select=c(rand58))

# run LFMM
toa.lfmm.rand58.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand58.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand58.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand58.k1 <- qvalue(toa.pv.rand58.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand58.k1 < 0.1))
# 15166

# which SNPs with a FDR > 10%
toa.FDR.rand58.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand58.k1 < 0.1)]

# rand59 #############################################################################################################################

# subset rand59
retained_preds <- subset(envs, select=c(rand59))

# run LFMM
toa.lfmm.rand59.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand59.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand59.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand59.k1 <- qvalue(toa.pv.rand59.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand59.k1 < 0.1))
# 12391

# which SNPs with a FDR > 10%
toa.FDR.rand59.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand59.k1 < 0.1)]

# rand60 #############################################################################################################################

# subset rand60
retained_preds <- subset(envs, select=c(rand60))

# run LFMM
toa.lfmm.rand60.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand60.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand60.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand60.k1 <- qvalue(toa.pv.rand60.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand60.k1 < 0.1))
# 893

# which SNPs with a FDR > 10%
toa.FDR.rand60.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand60.k1 < 0.1)]

# rand61 #############################################################################################################################

# subset rand61
retained_preds <- subset(envs, select=c(rand61))

# run LFMM
toa.lfmm.rand61.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand61.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand61.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand61.k1 <- qvalue(toa.pv.rand61.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand61.k1 < 0.1))
# 2797

# which SNPs with a FDR > 10%
toa.FDR.rand61.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand61.k1 < 0.1)]

# rand62 #############################################################################################################################

# subset rand62
retained_preds <- subset(envs, select=c(rand62))

# run LFMM
toa.lfmm.rand62.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand62.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand62.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand62.k1 <- qvalue(toa.pv.rand62.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand62.k1 < 0.1))
# 154

# which SNPs with a FDR > 10%
toa.FDR.rand62.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand62.k1 < 0.1)]

# rand63 #############################################################################################################################

# subset rand63
retained_preds <- subset(envs, select=c(rand63))

# run LFMM
toa.lfmm.rand63.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand63.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand63.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand63.k1 <- qvalue(toa.pv.rand63.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand63.k1 < 0.1))
# 626

# which SNPs with a FDR > 10%
toa.FDR.rand63.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand63.k1 < 0.1)]

# rand64 #############################################################################################################################

# subset rand64
retained_preds <- subset(envs, select=c(rand64))

# run LFMM
toa.lfmm.rand64.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand64.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand64.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand64.k1 <- qvalue(toa.pv.rand64.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand64.k1 < 0.1))
# 7469

# which SNPs with a FDR > 10%
toa.FDR.rand64.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand64.k1 < 0.1)]

# rand65 #############################################################################################################################

# subset rand65
retained_preds <- subset(envs, select=c(rand65))

# run LFMM
toa.lfmm.rand65.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand65.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand65.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand65.k1 <- qvalue(toa.pv.rand65.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand65.k1 < 0.1))
# 1527

# which SNPs with a FDR > 10%
toa.FDR.rand65.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand65.k1 < 0.1)]

# rand66 #############################################################################################################################

# subset rand66
retained_preds <- subset(envs, select=c(rand66))

# run LFMM
toa.lfmm.rand66.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand66.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand66.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand66.k1 <- qvalue(toa.pv.rand66.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand66.k1 < 0.1))
# 232

# which SNPs with a FDR > 10%
toa.FDR.rand66.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand66.k1 < 0.1)]

# rand67 #############################################################################################################################

# subset rand67
retained_preds <- subset(envs, select=c(rand67))

# run LFMM
toa.lfmm.rand67.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand67.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand67.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand67.k1 <- qvalue(toa.pv.rand67.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand67.k1 < 0.1))
# 2522

# which SNPs with a FDR > 10%
toa.FDR.rand67.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand67.k1 < 0.1)]

# rand68 #############################################################################################################################

# subset rand68
retained_preds <- subset(envs, select=c(rand68))

# run LFMM
toa.lfmm.rand68.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand68.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand68.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand68.k1 <- qvalue(toa.pv.rand68.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand68.k1 < 0.1))
# 451

# which SNPs with a FDR > 10%
toa.FDR.rand68.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand68.k1 < 0.1)]

# rand69 #############################################################################################################################

# subset rand69
retained_preds <- subset(envs, select=c(rand69))

# run LFMM
toa.lfmm.rand69.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand69.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand69.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand69.k1 <- qvalue(toa.pv.rand69.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand69.k1 < 0.1))
# 62

# which SNPs with a FDR > 10%
toa.FDR.rand69.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand69.k1 < 0.1)]

# rand70 #############################################################################################################################

# subset rand70
retained_preds <- subset(envs, select=c(rand70))

# run LFMM
toa.lfmm.rand70.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand70.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand70.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand70.k1 <- qvalue(toa.pv.rand70.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand70.k1 < 0.1))
# 273

# which SNPs with a FDR > 10%
toa.FDR.rand70.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand70.k1 < 0.1)]

# rand71 #############################################################################################################################

# subset rand71
retained_preds <- subset(envs, select=c(rand71))

# run LFMM
toa.lfmm.rand71.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand71.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand71.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand71.k1 <- qvalue(toa.pv.rand71.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand71.k1 < 0.1))
# 159

# which SNPs with a FDR > 10%
toa.FDR.rand71.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand71.k1 < 0.1)]

# rand72 #############################################################################################################################

# subset rand72
retained_preds <- subset(envs, select=c(rand72))

# run LFMM
toa.lfmm.rand72.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand72.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand72.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand72.k1 <- qvalue(toa.pv.rand72.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand72.k1 < 0.1))
# 542

# which SNPs with a FDR > 10%
toa.FDR.rand72.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand72.k1 < 0.1)]

# rand73 #############################################################################################################################

# subset rand73
retained_preds <- subset(envs, select=c(rand73))

# run LFMM
toa.lfmm.rand73.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand73.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand73.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand73.k1 <- qvalue(toa.pv.rand73.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand73.k1 < 0.1))
# 335

# which SNPs with a FDR > 10%
toa.FDR.rand73.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand73.k1 < 0.1)]

# rand74 #############################################################################################################################

# subset rand74
retained_preds <- subset(envs, select=c(rand74))

# run LFMM
toa.lfmm.rand74.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand74.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand74.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand74.k1 <- qvalue(toa.pv.rand74.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand74.k1 < 0.1))
# 7345

# which SNPs with a FDR > 10%
toa.FDR.rand74.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand74.k1 < 0.1)]

# rand75 #############################################################################################################################

# subset rand75
retained_preds <- subset(envs, select=c(rand75))

# run LFMM
toa.lfmm.rand75.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand75.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand75.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand75.k1 <- qvalue(toa.pv.rand75.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand75.k1 < 0.1))
# 15010

# which SNPs with a FDR > 10%
toa.FDR.rand75.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand75.k1 < 0.1)]

# rand76 #############################################################################################################################

# subset rand76
retained_preds <- subset(envs, select=c(rand76))

# run LFMM
toa.lfmm.rand76.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand76.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand76.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand76.k1 <- qvalue(toa.pv.rand76.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand76.k1 < 0.1))
# 862

# which SNPs with a FDR > 10%
toa.FDR.rand76.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand76.k1 < 0.1)]

# rand77 #############################################################################################################################

# subset rand77
retained_preds <- subset(envs, select=c(rand77))

# run LFMM
toa.lfmm.rand77.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand77.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand77.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand77.k1 <- qvalue(toa.pv.rand77.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand77.k1 < 0.1))
# 243

# which SNPs with a FDR > 10%
toa.FDR.rand77.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand77.k1 < 0.1)]

# rand78 #############################################################################################################################

# subset rand78
retained_preds <- subset(envs, select=c(rand78))

# run LFMM
toa.lfmm.rand78.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand78.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand78.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand78.k1 <- qvalue(toa.pv.rand78.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand78.k1 < 0.1))
# 1609

# which SNPs with a FDR > 10%
toa.FDR.rand78.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand78.k1 < 0.1)]

# rand79 #############################################################################################################################

# subset rand79
retained_preds <- subset(envs, select=c(rand79))

# run LFMM
toa.lfmm.rand79.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand79.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand79.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand79.k1 <- qvalue(toa.pv.rand79.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand79.k1 < 0.1))
# 781

# which SNPs with a FDR > 10%
toa.FDR.rand79.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand79.k1 < 0.1)]

# rand80 #############################################################################################################################

# subset rand80
retained_preds <- subset(envs, select=c(rand80))

# run LFMM
toa.lfmm.rand80.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand80.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand80.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand80.k1 <- qvalue(toa.pv.rand80.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand80.k1 < 0.1))
# 607

# which SNPs with a FDR > 10%
toa.FDR.rand80.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand80.k1 < 0.1)]

# rand81 #############################################################################################################################

# subset rand81
retained_preds <- subset(envs, select=c(rand81))

# run LFMM
toa.lfmm.rand81.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand81.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand81.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand81.k1 <- qvalue(toa.pv.rand81.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand81.k1 < 0.1))
# 872

# which SNPs with a FDR > 10%
toa.FDR.rand81.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand81.k1 < 0.1)]

# rand82 #############################################################################################################################

# subset rand82
retained_preds <- subset(envs, select=c(rand82))

# run LFMM
toa.lfmm.rand82.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand82.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand82.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand82.k1 <- qvalue(toa.pv.rand82.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand82.k1 < 0.1))
# 850

# which SNPs with a FDR > 10%
toa.FDR.rand82.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand82.k1 < 0.1)]

# rand83 #############################################################################################################################

# subset rand83
retained_preds <- subset(envs, select=c(rand83))

# run LFMM
toa.lfmm.rand83.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand83.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand83.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand83.k1 <- qvalue(toa.pv.rand83.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand83.k1 < 0.1))
# 2767

# which SNPs with a FDR > 10%
toa.FDR.rand83.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand83.k1 < 0.1)]

# rand84 #############################################################################################################################

# subset rand84
retained_preds <- subset(envs, select=c(rand84))

# run LFMM
toa.lfmm.rand84.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand84.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand84.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand84.k1 <- qvalue(toa.pv.rand84.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand84.k1 < 0.1))
# 148

# which SNPs with a FDR > 10%
toa.FDR.rand84.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand84.k1 < 0.1)]

# rand85 #############################################################################################################################

# subset rand85
retained_preds <- subset(envs, select=c(rand85))

# run LFMM
toa.lfmm.rand85.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand85.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand85.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand85.k1 <- qvalue(toa.pv.rand85.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand85.k1 < 0.1))
# 8003

# which SNPs with a FDR > 10%
toa.FDR.rand85.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand85.k1 < 0.1)]

# rand86 #############################################################################################################################

# subset rand86
retained_preds <- subset(envs, select=c(rand86))

# run LFMM
toa.lfmm.rand86.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand86.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand86.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand86.k1 <- qvalue(toa.pv.rand86.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand86.k1 < 0.1))
# 9329

# which SNPs with a FDR > 10%
toa.FDR.rand86.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand86.k1 < 0.1)]

# rand87 #############################################################################################################################

# subset rand87
retained_preds <- subset(envs, select=c(rand87))

# run LFMM
toa.lfmm.rand87.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand87.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand87.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand87.k1 <- qvalue(toa.pv.rand87.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand87.k1 < 0.1))
# 3499

# which SNPs with a FDR > 10%
toa.FDR.rand87.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand87.k1 < 0.1)]

# rand88 #############################################################################################################################

# subset rand88
retained_preds <- subset(envs, select=c(rand88))

# run LFMM
toa.lfmm.rand88.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand88.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand88.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand88.k1 <- qvalue(toa.pv.rand88.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand88.k1 < 0.1))
# 1183

# which SNPs with a FDR > 10%
toa.FDR.rand88.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand88.k1 < 0.1)]

# rand89 #############################################################################################################################

# subset rand89
retained_preds <- subset(envs, select=c(rand89))

# run LFMM
toa.lfmm.rand89.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand89.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand89.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand89.k1 <- qvalue(toa.pv.rand89.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand89.k1 < 0.1))
# 14967

# which SNPs with a FDR > 10%
toa.FDR.rand89.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand89.k1 < 0.1)]

# rand90 #############################################################################################################################

# subset rand90
retained_preds <- subset(envs, select=c(rand90))

# run LFMM
toa.lfmm.rand90.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand90.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand90.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand90.k1 <- qvalue(toa.pv.rand90.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand90.k1 < 0.1))
# 9306

# which SNPs with a FDR > 10%
toa.FDR.rand90.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand90.k1 < 0.1)]

# rand91 #############################################################################################################################

# subset rand91
retained_preds <- subset(envs, select=c(rand91))

# run LFMM
toa.lfmm.rand91.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand91.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand91.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand91.k1 <- qvalue(toa.pv.rand91.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand91.k1 < 0.1))
# 2796

# which SNPs with a FDR > 10%
toa.FDR.rand91.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand91.k1 < 0.1)]

# rand92 #############################################################################################################################

# subset rand92
retained_preds <- subset(envs, select=c(rand92))

# run LFMM
toa.lfmm.rand92.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand92.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand92.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand92.k1 <- qvalue(toa.pv.rand92.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand92.k1 < 0.1))
# 506

# which SNPs with a FDR > 10%
toa.FDR.rand92.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand92.k1 < 0.1)]

# rand93 #############################################################################################################################

# subset rand93
retained_preds <- subset(envs, select=c(rand93))

# run LFMM
toa.lfmm.rand93.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand93.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand93.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand93.k1 <- qvalue(toa.pv.rand93.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand93.k1 < 0.1))
# 758

# which SNPs with a FDR > 10%
toa.FDR.rand93.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand93.k1 < 0.1)]

# rand94 #############################################################################################################################

# subset rand94
retained_preds <- subset(envs, select=c(rand94))

# run LFMM
toa.lfmm.rand94.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand94.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand94.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand94.k1 <- qvalue(toa.pv.rand94.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand94.k1 < 0.1))
# 1393

# which SNPs with a FDR > 10%
toa.FDR.rand94.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand94.k1 < 0.1)]

# rand95 #############################################################################################################################

# subset rand95
retained_preds <- subset(envs, select=c(rand95))

# run LFMM
toa.lfmm.rand95.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand95.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand95.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand95.k1 <- qvalue(toa.pv.rand95.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand95.k1 < 0.1))
# 103

# which SNPs with a FDR > 10%
toa.FDR.rand95.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand95.k1 < 0.1)]

# rand96 #############################################################################################################################

# subset rand96
retained_preds <- subset(envs, select=c(rand96))

# run LFMM
toa.lfmm.rand96.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand96.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand96.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand96.k1 <- qvalue(toa.pv.rand96.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand96.k1 < 0.1))
# 945

# which SNPs with a FDR > 10%
toa.FDR.rand96.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand96.k1 < 0.1)]

# rand97 #############################################################################################################################

# subset rand97
retained_preds <- subset(envs, select=c(rand97))

# run LFMM
toa.lfmm.rand97.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand97.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand97.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand97.k1 <- qvalue(toa.pv.rand97.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand97.k1 < 0.1))
# 18238

# which SNPs with a FDR > 10%
toa.FDR.rand97.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand97.k1 < 0.1)]

# rand98 #############################################################################################################################

# subset rand98
retained_preds <- subset(envs, select=c(rand98))

# run LFMM
toa.lfmm.rand98.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand98.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand98.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand98.k1 <- qvalue(toa.pv.rand98.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand98.k1 < 0.1))
# 1631

# which SNPs with a FDR > 10%
toa.FDR.rand98.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand98.k1 < 0.1)]

# rand99 #############################################################################################################################

# subset rand99
retained_preds <- subset(envs, select=c(rand99))

# run LFMM
toa.lfmm.rand99.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand99.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand99.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand99.k1 <- qvalue(toa.pv.rand99.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand99.k1 < 0.1))
# 9595

# which SNPs with a FDR > 10%
toa.FDR.rand99.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand99.k1 < 0.1)]

# rand100 #############################################################################################################################

# subset rand100
retained_preds <- subset(envs, select=c(rand100))

# run LFMM
toa.lfmm.rand100.k1 <- lfmm_ridge(Y=snps_unlinked_imp_5x, X=retained_preds, K=1)
# takes <10 seconds

# get p-values
toa.pv.rand100.k1 <- lfmm_test(Y=snps_unlinked_imp_5x, X=retained_preds, lfmm=toa.lfmm.rand100.k1, calibrate="gif")
# takes <1 minute

# calculate LFMM statistics
toa.qv.rand100.k1 <- qvalue(toa.pv.rand100.k1$calibrated.pvalue)$qvalues

# how many SNPs with a FDR > 10%
length(which(toa.qv.rand100.k1 < 0.1))
# 93

# which SNPs with a FDR > 10%
toa.FDR.rand100.k1 <- colnames(snps_unlinked_imp_5x)[which(toa.qv.rand100.k1 < 0.1)]

# combine outputs of candidate SNPs from random value LFMM runs 1 - 100
max_length <- max(length(toa.FDR.rand1.k1),
                  length(toa.FDR.rand2.k1),
                  length(toa.FDR.rand3.k1),
                  length(toa.FDR.rand4.k1),
                  length(toa.FDR.rand5.k1),
                  length(toa.FDR.rand6.k1),
                  length(toa.FDR.rand7.k1),
                  length(toa.FDR.rand8.k1),
                  length(toa.FDR.rand9.k1),
                  length(toa.FDR.rand10.k1),
                  length(toa.FDR.rand11.k1),
                  length(toa.FDR.rand12.k1),
                  length(toa.FDR.rand13.k1),
                  length(toa.FDR.rand14.k1),
                  length(toa.FDR.rand15.k1),
                  length(toa.FDR.rand16.k1),
                  length(toa.FDR.rand17.k1),
                  length(toa.FDR.rand18.k1),
                  length(toa.FDR.rand19.k1),
                  length(toa.FDR.rand20.k1),
                  length(toa.FDR.rand21.k1),
                  length(toa.FDR.rand22.k1),
                  length(toa.FDR.rand23.k1),
                  length(toa.FDR.rand24.k1),
                  length(toa.FDR.rand25.k1),
                  length(toa.FDR.rand26.k1),
                  length(toa.FDR.rand27.k1),
                  length(toa.FDR.rand28.k1),
                  length(toa.FDR.rand29.k1),
                  length(toa.FDR.rand30.k1),
                  length(toa.FDR.rand31.k1),
                  length(toa.FDR.rand32.k1),
                  length(toa.FDR.rand33.k1),
                  length(toa.FDR.rand34.k1),
                  length(toa.FDR.rand35.k1),
                  length(toa.FDR.rand36.k1),
                  length(toa.FDR.rand37.k1),
                  length(toa.FDR.rand38.k1),
                  length(toa.FDR.rand39.k1),
                  length(toa.FDR.rand40.k1),
                  length(toa.FDR.rand41.k1),
                  length(toa.FDR.rand42.k1),
                  length(toa.FDR.rand43.k1),
                  length(toa.FDR.rand44.k1),
                  length(toa.FDR.rand45.k1),
                  length(toa.FDR.rand46.k1),
                  length(toa.FDR.rand47.k1),
                  length(toa.FDR.rand48.k1),
                  length(toa.FDR.rand49.k1),
                  length(toa.FDR.rand50.k1),
                  length(toa.FDR.rand51.k1),
                  length(toa.FDR.rand52.k1),
                  length(toa.FDR.rand53.k1),
                  length(toa.FDR.rand54.k1),
                  length(toa.FDR.rand55.k1),
                  length(toa.FDR.rand56.k1),
                  length(toa.FDR.rand57.k1),
                  length(toa.FDR.rand58.k1),
                  length(toa.FDR.rand59.k1),
                  length(toa.FDR.rand60.k1),
                  length(toa.FDR.rand61.k1),
                  length(toa.FDR.rand62.k1),
                  length(toa.FDR.rand63.k1),
                  length(toa.FDR.rand64.k1),
                  length(toa.FDR.rand65.k1),
                  length(toa.FDR.rand66.k1),
                  length(toa.FDR.rand67.k1),
                  length(toa.FDR.rand68.k1),
                  length(toa.FDR.rand69.k1),
                  length(toa.FDR.rand70.k1),
                  length(toa.FDR.rand71.k1),
                  length(toa.FDR.rand72.k1),
                  length(toa.FDR.rand73.k1),
                  length(toa.FDR.rand74.k1),
                  length(toa.FDR.rand75.k1),
                  length(toa.FDR.rand76.k1),
                  length(toa.FDR.rand77.k1),
                  length(toa.FDR.rand78.k1),
                  length(toa.FDR.rand79.k1),
                  length(toa.FDR.rand80.k1),
                  length(toa.FDR.rand81.k1),
                  length(toa.FDR.rand82.k1),
                  length(toa.FDR.rand83.k1),
                  length(toa.FDR.rand84.k1),
                  length(toa.FDR.rand85.k1),
                  length(toa.FDR.rand86.k1),
                  length(toa.FDR.rand87.k1),
                  length(toa.FDR.rand88.k1),
                  length(toa.FDR.rand89.k1),
                  length(toa.FDR.rand90.k1),
                  length(toa.FDR.rand91.k1))
max_length
# [1] 27168

max_length <- max(length(toa.FDR.rand92.k1),
                  length(toa.FDR.rand93.k1),
                  length(toa.FDR.rand94.k1),
                  length(toa.FDR.rand95.k1),
                  length(toa.FDR.rand96.k1),
                  length(toa.FDR.rand97.k1),
                  length(toa.FDR.rand98.k1),
                  length(toa.FDR.rand99.k1),
                  length(toa.FDR.rand100.k1))
max_length
# [1] 18238

# had to do in 2 parts - it wasn't working to do it all at once for some reason
# obviously the first with the 92 values has the max we ultimately want to use, so re-run 1-92 and store that value for the next steps

length(toa.FDR.rand1.k1) <- max_length
length(toa.FDR.rand2.k1) <- max_length
length(toa.FDR.rand3.k1) <- max_length
length(toa.FDR.rand4.k1) <- max_length
length(toa.FDR.rand5.k1) <- max_length
length(toa.FDR.rand6.k1) <- max_length
length(toa.FDR.rand7.k1) <- max_length
length(toa.FDR.rand8.k1) <- max_length
length(toa.FDR.rand9.k1) <- max_length
length(toa.FDR.rand10.k1) <- max_length
length(toa.FDR.rand11.k1) <- max_length
length(toa.FDR.rand12.k1) <- max_length
length(toa.FDR.rand13.k1) <- max_length
length(toa.FDR.rand14.k1) <- max_length
length(toa.FDR.rand15.k1) <- max_length
length(toa.FDR.rand16.k1) <- max_length
length(toa.FDR.rand17.k1) <- max_length
length(toa.FDR.rand18.k1) <- max_length
length(toa.FDR.rand19.k1) <- max_length
length(toa.FDR.rand20.k1) <- max_length
length(toa.FDR.rand21.k1) <- max_length
length(toa.FDR.rand22.k1) <- max_length
length(toa.FDR.rand23.k1) <- max_length
length(toa.FDR.rand24.k1) <- max_length
length(toa.FDR.rand25.k1) <- max_length
length(toa.FDR.rand26.k1) <- max_length
length(toa.FDR.rand27.k1) <- max_length
length(toa.FDR.rand28.k1) <- max_length
length(toa.FDR.rand29.k1) <- max_length
length(toa.FDR.rand30.k1) <- max_length
length(toa.FDR.rand31.k1) <- max_length
length(toa.FDR.rand32.k1) <- max_length
length(toa.FDR.rand33.k1) <- max_length
length(toa.FDR.rand34.k1) <- max_length
length(toa.FDR.rand35.k1) <- max_length
length(toa.FDR.rand36.k1) <- max_length
length(toa.FDR.rand37.k1) <- max_length
length(toa.FDR.rand38.k1) <- max_length
length(toa.FDR.rand39.k1) <- max_length
length(toa.FDR.rand40.k1) <- max_length
length(toa.FDR.rand41.k1) <- max_length
length(toa.FDR.rand42.k1) <- max_length
length(toa.FDR.rand43.k1) <- max_length
length(toa.FDR.rand44.k1) <- max_length
length(toa.FDR.rand45.k1) <- max_length
length(toa.FDR.rand46.k1) <- max_length
length(toa.FDR.rand47.k1) <- max_length
length(toa.FDR.rand48.k1) <- max_length
length(toa.FDR.rand49.k1) <- max_length
length(toa.FDR.rand50.k1) <- max_length
length(toa.FDR.rand51.k1) <- max_length
length(toa.FDR.rand52.k1) <- max_length
length(toa.FDR.rand53.k1) <- max_length
length(toa.FDR.rand54.k1) <- max_length
length(toa.FDR.rand55.k1) <- max_length
length(toa.FDR.rand56.k1) <- max_length
length(toa.FDR.rand57.k1) <- max_length
length(toa.FDR.rand58.k1) <- max_length
length(toa.FDR.rand59.k1) <- max_length
length(toa.FDR.rand60.k1) <- max_length
length(toa.FDR.rand61.k1) <- max_length
length(toa.FDR.rand62.k1) <- max_length
length(toa.FDR.rand63.k1) <- max_length
length(toa.FDR.rand64.k1) <- max_length
length(toa.FDR.rand65.k1) <- max_length
length(toa.FDR.rand66.k1) <- max_length
length(toa.FDR.rand67.k1) <- max_length
length(toa.FDR.rand68.k1) <- max_length
length(toa.FDR.rand69.k1) <- max_length
length(toa.FDR.rand70.k1) <- max_length
length(toa.FDR.rand71.k1) <- max_length
length(toa.FDR.rand72.k1) <- max_length
length(toa.FDR.rand73.k1) <- max_length
length(toa.FDR.rand74.k1) <- max_length
length(toa.FDR.rand75.k1) <- max_length
length(toa.FDR.rand76.k1) <- max_length
length(toa.FDR.rand77.k1) <- max_length
length(toa.FDR.rand78.k1) <- max_length
length(toa.FDR.rand79.k1) <- max_length
length(toa.FDR.rand80.k1) <- max_length
length(toa.FDR.rand81.k1) <- max_length
length(toa.FDR.rand82.k1) <- max_length
length(toa.FDR.rand83.k1) <- max_length
length(toa.FDR.rand84.k1) <- max_length
length(toa.FDR.rand85.k1) <- max_length
length(toa.FDR.rand86.k1) <- max_length
length(toa.FDR.rand87.k1) <- max_length
length(toa.FDR.rand88.k1) <- max_length
length(toa.FDR.rand89.k1) <- max_length
length(toa.FDR.rand90.k1) <- max_length
length(toa.FDR.rand91.k1) <- max_length
length(toa.FDR.rand92.k1) <- max_length
length(toa.FDR.rand93.k1) <- max_length
length(toa.FDR.rand94.k1) <- max_length
length(toa.FDR.rand95.k1) <- max_length
length(toa.FDR.rand96.k1) <- max_length
length(toa.FDR.rand97.k1) <- max_length
length(toa.FDR.rand98.k1) <- max_length
length(toa.FDR.rand99.k1) <- max_length
length(toa.FDR.rand100.k1) <- max_length

rand1_100 <- cbind(toa.FDR.rand1.k1,
                   toa.FDR.rand2.k1,
                   toa.FDR.rand3.k1,
                   toa.FDR.rand4.k1,
                   toa.FDR.rand5.k1,
                   toa.FDR.rand6.k1,
                   toa.FDR.rand7.k1,
                   toa.FDR.rand8.k1,
                   toa.FDR.rand9.k1,
                   toa.FDR.rand10.k1,
                   toa.FDR.rand11.k1,
                   toa.FDR.rand12.k1,
                   toa.FDR.rand13.k1,
                   toa.FDR.rand14.k1,
                   toa.FDR.rand15.k1,
                   toa.FDR.rand16.k1,
                   toa.FDR.rand17.k1,
                   toa.FDR.rand18.k1,
                   toa.FDR.rand19.k1,
                   toa.FDR.rand20.k1,
                   toa.FDR.rand21.k1,
                   toa.FDR.rand22.k1,
                   toa.FDR.rand23.k1,
                   toa.FDR.rand24.k1,
                   toa.FDR.rand25.k1,
                   toa.FDR.rand26.k1,
                   toa.FDR.rand27.k1,
                   toa.FDR.rand28.k1,
                   toa.FDR.rand29.k1,
                   toa.FDR.rand30.k1,
                   toa.FDR.rand31.k1,
                   toa.FDR.rand32.k1,
                   toa.FDR.rand33.k1,
                   toa.FDR.rand34.k1,
                   toa.FDR.rand35.k1,
                   toa.FDR.rand36.k1,
                   toa.FDR.rand37.k1,
                   toa.FDR.rand38.k1,
                   toa.FDR.rand39.k1,
                   toa.FDR.rand40.k1,
                   toa.FDR.rand41.k1,
                   toa.FDR.rand42.k1,
                   toa.FDR.rand43.k1,
                   toa.FDR.rand44.k1,
                   toa.FDR.rand45.k1,
                   toa.FDR.rand46.k1,
                   toa.FDR.rand47.k1,
                   toa.FDR.rand48.k1,
                   toa.FDR.rand49.k1,
                   toa.FDR.rand50.k1,
                   toa.FDR.rand51.k1,
                   toa.FDR.rand52.k1,
                   toa.FDR.rand53.k1,
                   toa.FDR.rand54.k1,
                   toa.FDR.rand55.k1,
                   toa.FDR.rand56.k1,
                   toa.FDR.rand57.k1,
                   toa.FDR.rand58.k1,
                   toa.FDR.rand59.k1,
                   toa.FDR.rand60.k1,
                   toa.FDR.rand61.k1,
                   toa.FDR.rand62.k1,
                   toa.FDR.rand63.k1,
                   toa.FDR.rand64.k1,
                   toa.FDR.rand65.k1,
                   toa.FDR.rand66.k1,
                   toa.FDR.rand67.k1,
                   toa.FDR.rand68.k1,
                   toa.FDR.rand69.k1,
                   toa.FDR.rand70.k1,
                   toa.FDR.rand71.k1,
                   toa.FDR.rand72.k1,
                   toa.FDR.rand73.k1,
                   toa.FDR.rand74.k1,
                   toa.FDR.rand75.k1,
                   toa.FDR.rand76.k1,
                   toa.FDR.rand77.k1,
                   toa.FDR.rand78.k1,
                   toa.FDR.rand79.k1,
                   toa.FDR.rand80.k1,
                   toa.FDR.rand81.k1,
                   toa.FDR.rand82.k1,
                   toa.FDR.rand83.k1,
                   toa.FDR.rand84.k1,
                   toa.FDR.rand85.k1,
                   toa.FDR.rand86.k1,
                   toa.FDR.rand87.k1,
                   toa.FDR.rand88.k1,
                   toa.FDR.rand89.k1,
                   toa.FDR.rand90.k1,
                   toa.FDR.rand91.k1,
                   toa.FDR.rand92.k1,
                   toa.FDR.rand93.k1,
                   toa.FDR.rand94.k1,
                   toa.FDR.rand95.k1,
                   toa.FDR.rand96.k1,
                   toa.FDR.rand97.k1,
                   toa.FDR.rand98.k1,
                   toa.FDR.rand99.k1,
                   toa.FDR.rand100.k1)

# save false candidate SNPs from random value LFMM runs 1-100 to .csv file
write.csv(rand1_100,"lfmm1_5x_n24_SG_B_m8_9_10_rand1-100_k1_putative_adaptive_SNPs_0.1.csv")