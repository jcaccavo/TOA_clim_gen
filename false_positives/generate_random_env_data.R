# generate random environmental data
# pull n values for each random environmental variable from a standard normal distribution
# rnorm (n, mean = x, sd = y)
# n – number of observations we want rnorm to return
# mean – mean value of the normal distribution we are using
# sd – standard deviation of the normal distribution we are using
# for a standard normal distribution, mean = 0 and sd = 1
# do this 8 times for each environmental variable

current_m <- rnorm(24, 0, 1)
current_z <- rnorm(24, 0, 1)
sst <- rnorm(24, 0, 1)
sal <- rnorm(24, 0, 1)
sea_ice_thick <- rnorm(24, 0, 1)
sea_ice_conc <- rnorm(24, 0, 1)
mxl0.01 <- rnorm(24, 0, 1)
mxl0.03 <- rnorm(24, 0, 1)

rand_envs <- cbind(current_m, current_z, sst, sal, sea_ice_thick, sea_ice_conc, mxl0.01, mxl0.03)

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/rda23_5x_n24_rand")

write.csv(rand_envs,"rda23_5x_n24_rand_envs.csv")

# repeat process now 3 times for each PC

rand1 <- rnorm(24, 0, 1)
rand2 <- rnorm(24, 0, 1)
rand3 <- rnorm(24, 0, 1)

rand_envs <- cbind(rand1, rand2, rand3)

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/rda26_5x_n24_rand3")

write.csv(rand_envs,"rda26_5x_n24_rand3_PCs.csv")

# repeat process not 97 times for each PC
rand4 <- rnorm(24, 0, 1)
rand5 <- rnorm(24, 0, 1)
rand6 <- rnorm(24, 0, 1)
rand7 <- rnorm(24, 0, 1)
rand8 <- rnorm(24, 0, 1)
rand9 <- rnorm(24, 0, 1)
rand10 <- rnorm(24, 0, 1)
rand11 <- rnorm(24, 0, 1)
rand12 <- rnorm(24, 0, 1)
rand13 <- rnorm(24, 0, 1)
rand14 <- rnorm(24, 0, 1)
rand15 <- rnorm(24, 0, 1)
rand16 <- rnorm(24, 0, 1)
rand17 <- rnorm(24, 0, 1)
rand18 <- rnorm(24, 0, 1)
rand19 <- rnorm(24, 0, 1)
rand20 <- rnorm(24, 0, 1)
rand21 <- rnorm(24, 0, 1)
rand22 <- rnorm(24, 0, 1)
rand23 <- rnorm(24, 0, 1)
rand24 <- rnorm(24, 0, 1)
rand25 <- rnorm(24, 0, 1)
rand26 <- rnorm(24, 0, 1)
rand27 <- rnorm(24, 0, 1)
rand28 <- rnorm(24, 0, 1)
rand29 <- rnorm(24, 0, 1)
rand30 <- rnorm(24, 0, 1)
rand31 <- rnorm(24, 0, 1)
rand32 <- rnorm(24, 0, 1)
rand33 <- rnorm(24, 0, 1)
rand34 <- rnorm(24, 0, 1)
rand35 <- rnorm(24, 0, 1)
rand36 <- rnorm(24, 0, 1)
rand37 <- rnorm(24, 0, 1)
rand38 <- rnorm(24, 0, 1)
rand39 <- rnorm(24, 0, 1)
rand40 <- rnorm(24, 0, 1)
rand41 <- rnorm(24, 0, 1)
rand42 <- rnorm(24, 0, 1)
rand43 <- rnorm(24, 0, 1)
rand44 <- rnorm(24, 0, 1)
rand45 <- rnorm(24, 0, 1)
rand46 <- rnorm(24, 0, 1)
rand47 <- rnorm(24, 0, 1)
rand48 <- rnorm(24, 0, 1)
rand49 <- rnorm(24, 0, 1)
rand50 <- rnorm(24, 0, 1)
rand51 <- rnorm(24, 0, 1)
rand52 <- rnorm(24, 0, 1)
rand53 <- rnorm(24, 0, 1)
rand54 <- rnorm(24, 0, 1)
rand55 <- rnorm(24, 0, 1)
rand56 <- rnorm(24, 0, 1)
rand57 <- rnorm(24, 0, 1)
rand58 <- rnorm(24, 0, 1)
rand59 <- rnorm(24, 0, 1)
rand60 <- rnorm(24, 0, 1)
rand61 <- rnorm(24, 0, 1)
rand62 <- rnorm(24, 0, 1)
rand63 <- rnorm(24, 0, 1)
rand64 <- rnorm(24, 0, 1)
rand65 <- rnorm(24, 0, 1)
rand66 <- rnorm(24, 0, 1)
rand67 <- rnorm(24, 0, 1)
rand68 <- rnorm(24, 0, 1)
rand69 <- rnorm(24, 0, 1)
rand70 <- rnorm(24, 0, 1)
rand71 <- rnorm(24, 0, 1)
rand72 <- rnorm(24, 0, 1)
rand73 <- rnorm(24, 0, 1)
rand74 <- rnorm(24, 0, 1)
rand75 <- rnorm(24, 0, 1)
rand76 <- rnorm(24, 0, 1)
rand77 <- rnorm(24, 0, 1)
rand78 <- rnorm(24, 0, 1)
rand79 <- rnorm(24, 0, 1)
rand80 <- rnorm(24, 0, 1)
rand81 <- rnorm(24, 0, 1)
rand82 <- rnorm(24, 0, 1)
rand83 <- rnorm(24, 0, 1)
rand84 <- rnorm(24, 0, 1)
rand85 <- rnorm(24, 0, 1)
rand86 <- rnorm(24, 0, 1)
rand87 <- rnorm(24, 0, 1)
rand88 <- rnorm(24, 0, 1)
rand89 <- rnorm(24, 0, 1)
rand90 <- rnorm(24, 0, 1)
rand91 <- rnorm(24, 0, 1)
rand92 <- rnorm(24, 0, 1)
rand93 <- rnorm(24, 0, 1)
rand94 <- rnorm(24, 0, 1)
rand95 <- rnorm(24, 0, 1)
rand96 <- rnorm(24, 0, 1)
rand97 <- rnorm(24, 0, 1)
rand98 <- rnorm(24, 0, 1)
rand99 <- rnorm(24, 0, 1)
rand100 <- rnorm(24, 0, 1)

rand_envs <- cbind(rand4,
                   rand5,
                   rand6,
                   rand7,
                   rand8,
                   rand9,
                   rand10,
                   rand11,
                   rand12,
                   rand13,
                   rand14,
                   rand15,
                   rand16,
                   rand17,
                   rand18,
                   rand19,
                   rand20,
                   rand21,
                   rand22,
                   rand23,
                   rand24,
                   rand25,
                   rand26,
                   rand27,
                   rand28,
                   rand29,
                   rand30,
                   rand31,
                   rand32,
                   rand33,
                   rand34,
                   rand35,
                   rand36,
                   rand37,
                   rand38,
                   rand39,
                   rand40,
                   rand41,
                   rand42,
                   rand43,
                   rand44,
                   rand45,
                   rand46,
                   rand47,
                   rand48,
                   rand49,
                   rand50,
                   rand51,
                   rand52,
                   rand53,
                   rand54,
                   rand55,
                   rand56,
                   rand57,
                   rand58,
                   rand59,
                   rand60,
                   rand61,
                   rand62,
                   rand63,
                   rand64,
                   rand65,
                   rand66,
                   rand67,
                   rand68,
                   rand69,
                   rand70,
                   rand71,
                   rand72,
                   rand73,
                   rand74,
                   rand75,
                   rand76,
                   rand77,
                   rand78,
                   rand79,
                   rand80,
                   rand81,
                   rand82,
                   rand83,
                   rand84,
                   rand85,
                   rand86,
                   rand87,
                   rand88,
                   rand89,
                   rand90,
                   rand91,
                   rand92,
                   rand93,
                   rand94,
                   rand95,
                   rand96,
                   rand97,
                   rand98,
                   rand99,
                   rand100)

setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/rda/pRDA_rand_5x_n24_4-100")

write.csv(rand_envs,"pRDA_rand_5x_n24_4-100.csv")
