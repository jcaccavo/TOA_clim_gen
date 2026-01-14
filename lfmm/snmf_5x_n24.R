# helpful links
# http://membres-timc.imag.fr/Olivier.Francois/LEA/tutorial.htm
# https://rdocumentation.org/packages/LEA/versions/1.4.0/topics/lfmm
# https://rdrr.io/bioc/LEA/man/lfmm2.html
# https://www.bioconductor.org/packages/release/bioc/html/LEA.html
# https://bioconductor.org/packages/release/bioc/html/qvalue.html

# run on server (yangtze)
# Run in R in directory /srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/05_lfmm/#RELEVANT RUN FOLDER

# in R
# for some reason, R package LEA not installing using usual conda method
### conda install -c conda-forge r-<pkg>
# or
### conda install -c r r-<pkg>
# instead, had to install within R using BiocManager using the following command
### BiocManager::install("LEA")
# but R was not happy with the version of Bioconductor that I had, giving the following error
### Error: Bioconductor version '3.16' requires R version '4.2'; use `version = '3.18'`
### with R version 4.3; see https://bioconductor.org/install
# so I updated the version of BiocManager to 3.18 in R with the following command
### BiocManager::install(version = "3.18")
# after having done this, I was then able to install LEA with the original suggested command
### BiocManager::install("LEA")

# in R
# needed to install package qvalue as well
# here too, the usual conda method in the prompt did not work
# had to go into R and use
### BiocManager::install("qvalue")

# in R
# Load packages
library(vegan)
library(LEA)
library(qvalue)
library(adegenet)
library(lfmm)

# currently in jc_screen_5

# use snmf to inform choice of k for lfmm
# using 20 repetitions per Varas-Myrik et al. 2022 https://doi.org/10.1016/j.foreco.2021.119856
toa.snmf <- snmf("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/24_gea/03_gen_data/5x_TOA_only_filtered_SNPs_all.ped",
                 K = 1:10, 
                 entropy = TRUE, 
                 repetitions = 20,
                 project = "new")
# run was aborted

# trying again - but because the folders created in this analysis are now being used for the lfmm that is ongoing, will run this on the ped file in the plink folder
toa.snmf <- snmf("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/21_plink/5X/5x_TOA_only_filtered_SNPs_all.ped",
                 K = 1:10, 
                 entropy = TRUE, 
                 repetitions = 20,
                 project = "new")
# running on jc_screen_2, began 31/08 11:45
# finished by 01/09 16:45 
# took <30 hours

# show the project
show(toa.snmf)

# ***** run ***** 
#   snmf class
# 
# file directory:                   K8/run20/ 
#   Q output file:                    5x_TOA_only_filtered_SNPs_all_r20.8.Q 
# G output file:                    5x_TOA_only_filtered_SNPs_all_r20.8.G 
# snmfClass file:                   5x_TOA_only_filtered_SNPs_all_r20.8.snmfClass 
# number of ancestral populations:  8 
# run number:                       20 
# regularization parameter:         10 
# number of CPUs:                   1 
# seed:                             2008861373 
# maximal number of iterations:     200 
# tolerance error:                  1e-05 
# Q input file:                      
#   cross-Entropy:                    1.235395 
# 
# ***** run *****
#   snmf class
# 
# file directory:                   K9/run20/ 
#   Q output file:                    5x_TOA_only_filtered_SNPs_all_r20.9.Q 
# G output file:                    5x_TOA_only_filtered_SNPs_all_r20.9.G 
# snmfClass file:                   5x_TOA_only_filtered_SNPs_all_r20.9.snmfClass 
# number of ancestral populations:  9 
# run number:                       20 
# regularization parameter:         10 
# number of CPUs:                   1 
# seed:                             2008861373 
# maximal number of iterations:     200 
# tolerance error:                  1e-05 
# Q input file:                      
#   cross-Entropy:                    1.335084 
# 
# ***** run *****
#   snmf class
# 
# file directory:                   K10/run20/ 
#   Q output file:                    5x_TOA_only_filtered_SNPs_all_r20.10.Q 
# G output file:                    5x_TOA_only_filtered_SNPs_all_r20.10.G 
# snmfClass file:                   5x_TOA_only_filtered_SNPs_all_r20.10.snmfClass 
# number of ancestral populations:  10 
# run number:                       20 
# regularization parameter:         10 
# number of CPUs:                   1 
# seed:                             2008861373 
# maximal number of iterations:     200 
# tolerance error:                  1e-05 
# Q input file:                      
#   cross-Entropy:                    1.449414 

# summary of the project
summary(toa.snmf)
# $repetitions
# K = 1 K = 2 K = 3 K = 4 K = 5 K = 6 K = 7 K = 8 K = 9
# with cross-entropy       20    20    20    20    20    20    20    20    20
# without cross-entropy     0     0     0     0     0     0     0     0     0
# total                    20    20    20    20    20    20    20    20    20
# K = 10
# with cross-entropy        20
# without cross-entropy      0
# total                     20
# 
# $crossEntropy
#       K = 1     K = 2     K = 3     K = 4     K = 5    K = 6    K = 7
# min  0.7185032 0.7447634 0.8062941 0.8877192 0.9752125 1.036098 1.131819
# mean 0.7196043 0.7491623 0.8134001 0.9036899 0.9968014 1.099785 1.185482
# max  0.7205038 0.7504007 0.8380674 0.9202292 1.0251004 1.142485 1.234553
#       K = 8    K = 9   K = 10
# min  1.235395 1.317830 1.407069
# mean 1.284267 1.381388 1.469744
# max  1.358882 1.450633 1.529008

# plot cross-entropy criterion ####
# of all runs of the toa.snmf
pdf("snmf_5x_n24_cross-entropy_criterion_1.pdf")
plot(toa.snmf, cex = 1.2, col = "lightblue", pch = 19)
dev.off()

# K = 4 ####

# get the cross-entropy of all runs for K = 4
ce4 <- cross.entropy(toa.snmf, K = 4)
ce4
# K = 4
# run 1  0.9069158
# run 2  0.8914107
# run 3  0.8877192
# run 4  0.9202292
# run 5  0.9078812
# run 6  0.8978856
# run 7  0.9128234
# run 8  0.9115645
# run 9  0.8923741
# run 10 0.9120675
# run 11 0.9047712
# run 12 0.9048179
# run 13 0.9066553
# run 14 0.8996530
# run 15 0.8972984
# run 16 0.8962061
# run 17 0.9188126
# run 18 0.8990435
# run 19 0.8983245
# run 20 0.9073451

# select the run with the lowest cross-entropy for K = 4
best4 <- which.min(ce4)
# [1] 3
# run 3  0.8877192

# K = 3 ####

# get the cross-entropy of all runs for K = 3
ce3 <- cross.entropy(toa.snmf, K = 3)
ce3
# run 1  0.8153636
# run 2  0.8135562
# run 3  0.8125079
# run 4  0.8090393
# run 5  0.8380674
# run 6  0.8123587
# run 7  0.8121740
# run 8  0.8148451
# run 9  0.8107250
# run 10 0.8143378
# run 11 0.8142780
# run 12 0.8078566
# run 13 0.8132141
# run 14 0.8145915
# run 15 0.8096345
# run 16 0.8062941
# run 17 0.8067720
# run 18 0.8141958
# run 19 0.8142348
# run 20 0.8139550

# select the run with the lowest cross-entropy for K = 3
best3 <- which.min(ce3)
# [1] 16
# run 16 0.8062941

# K = 2 ####

# get the cross-entropy of all runs for K = 2
ce2 <- cross.entropy(toa.snmf, K = 2)
ce2
# run 1  0.7504007
# run 2  0.7496951
# run 3  0.7493028
# run 4  0.7491740
# run 5  0.7503453
# run 6  0.7474912
# run 7  0.7498498
# run 8  0.7495604
# run 9  0.7493185
# run 10 0.7493563
# run 11 0.7496026
# run 12 0.7492387
# run 13 0.7466095
# run 14 0.7500824
# run 15 0.7501500
# run 16 0.7502350
# run 17 0.7484280
# run 18 0.7499824
# run 19 0.7447634
# run 20 0.7496608

# select the run with the lowest cross-entropy for K = 2
best2 <- which.min(ce2)
# [1] 19
# run 19 0.7447634

# K = 1 ####

# get the cross-entropy of all runs for K = 1
ce1 <- cross.entropy(toa.snmf, K = 1)
ce1
# run 1  0.7202784
# run 2  0.7197083
# run 3  0.7190252
# run 4  0.7194705
# run 5  0.7202093
# run 6  0.7185032
# run 7  0.7197199
# run 8  0.7195682
# run 9  0.7189918
# run 10 0.7197307
# run 11 0.7194435
# run 12 0.7191003
# run 13 0.7192647
# run 14 0.7201777
# run 15 0.7205038
# run 16 0.7198976
# run 17 0.7190375
# run 18 0.7198019
# run 19 0.7199969
# run 20 0.7196558

# select the run with the lowest cross-entropy for K = 1
best1 <- which.min(ce1)
# [1] 6
# run 6  0.7185032

# K = 5 ####

# get the cross-entropy of all runs for K = 5
ce5 <- cross.entropy(toa.snmf, K = 5)
ce5
# run 1  0.9828280
# run 2  1.0030936
# run 3  1.0031706
# run 4  0.9901861
# run 5  1.0251004
# run 6  0.9833961
# run 7  1.0093750
# run 8  0.9818957
# run 9  0.9982080
# run 10 0.9995602
# run 11 0.9919352
# run 12 0.9824775
# run 13 1.0071128
# run 14 0.9930475
# run 15 0.9981470
# run 16 0.9948296
# run 17 0.9752125
# run 18 1.0015416
# run 19 1.0207446
# run 20 0.9941656

# select the run with the lowest cross-entropy for K = 5
best5 <- which.min(ce5)
# [1] 17
# run 17 0.9752125

# K = 7 ####

# get the cross-entropy of all runs for K = 7
ce7 <- cross.entropy(toa.snmf, K = 7)
ce7
# run 1  1.213736
# run 2  1.131819
# run 3  1.185517
# run 4  1.190535
# run 5  1.234553
# run 6  1.164109
# run 7  1.159286
# run 8  1.164687
# run 9  1.204726
# run 10 1.175700
# run 11 1.155375
# run 12 1.206660
# run 13 1.187173
# run 14 1.191512
# run 15 1.180441
# run 16 1.195152
# run 17 1.211951
# run 18 1.200090
# run 19 1.172672
# run 20 1.183941

# select the run with the lowest cross-entropy for K = 7
best7 <- which.min(ce7)
# [1] 2
# run 2  1.131819

# k = 1 has the lowest cross-entropy, so should be used in subsequent lfmm analyses

# create labels for snmf barplot ####
# for different categorizations

pop_k4 <- c("48.5_2013",
            "58.4.2_2019",
            "48.5_2013",
            "48.5_2013",
            "48.5_2013",
            "88.1_2020",
            "48.1_2001",
            "88.1_2001",
            "88.1_2001",
            "88.1_2001",
            "88.1_2020",
            "48.2_2019",
            "88.1_2020",
            "88.1_2020",
            "58.4.2_2019",
            "58.4.2_2019",
            "88.1_2020",
            "88.1_2020",
            "88.1_2020",
            "48.2_2019",
            "58.4.2_2001",
            "88.1_2001",
            "88.1_2001",
            "88.1_2001")

pop_k7 <- c("88.1_2020",
            "58.4.2_2019",
            "88.1_2020",
            "48.1_2001",
            "48.2_2019",
            "48.5_2013",
            "48.5_2013",
            "48.5_2013",
            "48.2_2019",
            "88.1_2020",
            "88.1_2001",
            "48.5_2013",
            "58.4.2_2019",
            "88.1_2020",
            "88.1_2001",
            "88.1_2001",
            "88.1_2001",
            "88.1_2020",
            "88.1_2020",
            "58.4.2_2019",
            "88.1_2020",
            "58.4.2_2001",
            "88.1_2001",
            "88.1_2001")

area_k4 <- c("48",
             "58",
             "48",
             "48",
             "48",
             "88",
             "48",
             "88",
             "88",
             "88",
             "88",
             "48",
             "88",
             "88",
             "58",
             "58",
             "88",
             "88",
             "88",
             "48",
             "58",
             "88",
             "88",
             "88")

area_k3 <- c("48",
             "48",
             "48",
             "48",
             "48",
             "58",
             "58",
             "48",
             "88",
             "88",
             "88",
             "88",
             "58",
             "88",
             "88",
             "88",
             "48",
             "58",
             "88",
             "88",
             "88",
             "88",
             "88",
             "88")

cohort_k4 <- c("post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "pre2000",
               "pre2000",
               "pre2000",
               "pre2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "pre2000",
               "pre2000",
               "pre2000",
               "pre2000")

cohort_k2 <- c("post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "post2000",
               "pre2000",
               "pre2000",
               "pre2000",
               "pre2000",
               "pre2000",
               "pre2000",
               "pre2000",
               "pre2000")

# colorblind-safe 7 diverging color palette for barplot from colorbrewer
# https://colorbrewer2.org/#type=diverging&scheme=BrBG&n=7
my.colors <- c("#8c510a","#c7eae5","#d8b365","#5ab4ac","#f6e8c3","#01665e","#f5f5f5")

# axis colors
# col.axis parameter from axis can only take one value (e.g. black) and not a vector of color assignments

# create a barplots of the snmf analysis ####

# by pop 
# k = 4
pdf("snmf_5x_n24_k4_barplot_pop_1.pdf")
barchart(toa.snmf, K = 4, run = best4,
         border = NA, space = 0,
         col = my.colors,
         ylab = "Ancestry proportions",
         main = "Ancestry matrix by population") -> bp
axis(1, at = 0.5:length(bp$order),
     labels = pop_k4, las=2,
     cex.axis = .7)
dev.off()

# k = 7
pdf("snmf_5x_n24_k7_barplot_pop_1.pdf")
barchart(toa.snmf, K = 7, run = best7,
         border = NA, space = 0,
         col = my.colors,
         ylab = "Ancestry proportions",
         main = "Ancestry matrix by population") -> bp
axis(1, at = 0.5:length(bp$order),
     labels = pop_k7, las=2,
     cex.axis = .7)
dev.off()

bp$order
# 1 10 15 17 16  3  4  6  2  9 23  5 13 12 21 22 24  7  8 11 14 18 19 20

# by area 
# k = 4
pdf("snmf_5x_n24_k4_barplot_area_1.pdf")
barchart(toa.snmf, K = 4, run = best4,
         border = NA, space = 0,
         col = my.colors,
         ylab = "Ancestry proportions",
         main = "Ancestry matrix by CCAMLR area") -> bp
axis(1, at = 0.5:length(bp$order),
     labels = area_k4, las=2,
     cex.axis = .7)
dev.off()

# k = 3
pdf("snmf_5x_n24_k3_barplot_area_1.pdf")
barchart(toa.snmf, K = 3, run = best3,
         border = NA, space = 0,
         col = my.colors,
         ylab = "Ancestry proportions",
         main = "Ancestry matrix by CCAMLR area") -> bp
axis(1, at = 0.5:length(bp$order),
     labels = area_k3, las=2,
     cex.axis = .7)
dev.off()

bp$order
# 3  4  5  6  2 11 13 16  1  7  8  9 10 12 14 15 17 18 19 20 21 22 23 24

# by cohort 
# k = 4
pdf("snmf_5x_n24_k4_barplot_cohort_1.pdf")
barchart(toa.snmf, K = 4, run = best4,
         border = NA, space = 0,
         col = my.colors,
         ylab = "Ancestry proportions",
         main = "Ancestry matrix by cohort") -> bp
axis(1, at = 0.5:length(bp$order),
     labels = cohort_k4, las=2,
     cex.axis = .7)
dev.off()

# by cohort 
# k = 2
pdf("snmf_5x_n24_k2_barplot_cohort_1.pdf")
barchart(toa.snmf, K = 2, run = best2,
         border = NA, space = 0,
         col = my.colors,
         ylab = "Ancestry proportions",
         main = "Ancestry matrix by cohort") -> bp
axis(1, at = 0.5:length(bp$order),
     labels = cohort_k2, las=2,
     cex.axis = .7)
dev.off()

bp$order
# 3  4  5  6  1  2  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

# display the Q-matrix
# Q.matrix <- as.qmatrix(Q(toa.snmf, K = 4, run = best))
# Error in as.qmatrix(Q(toa.snmf, K = 4, run = best)) : 
#   could not find function "as.qmatrix"

# get the ancestral genotype frequency matrix, G, for the 2nd run for K = 4. 
G.matrix <- G(toa.snmf, K = 4, run = 2)
# takes <30 seconds

# snmf runs are automatically saved into an snmf project directory.
# The name of the snmf project file is the same name as 
# the name of the input file with a .snmfProject extension ("genotypes.snmfProject").
# The name of the snmf project directory is the same name as
# the name of the input file with a .snmf extension ("genotypes.snmf/")
# A unique project includes all runs for each input file.

setwd("/srv/public/users/jcaccavo/11_CCGA_full_seq/02_NovaSeq/02_WG/21_plink/5X")

# re-load a snmf project ####
# An snmf project can be load in a different session.
toa.snmf <- load.snmfProject("5x_TOA_only_filtered_SNPs_all.snmfProject")

# An snmf project can be removed.
# Caution: All files associated with the project will be removed.
remove.snmfProject("5x_TOA_only_filtered_SNPs_all.snmfProject")