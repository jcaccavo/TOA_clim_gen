# conduct a PCA to reduce the dimensionality of the environmental variables for the rda

# helpful links
# https://medium.com/@victorallan/an-intuitive-guide-to-pca-in-r-a-step-by-step-tutorial-with-beautiful-visualization-examples-73bce5ee02e5
# http://www.sthda.com/english/wiki/get-pca-extract-the-results-for-individuals-variables-in-principal-component-analysis-r-software-and-data-mining # factoextra
# https://easystats.github.io/parameters/articles/parameters_reduction.html # parameters - feature reduction
# https://easystats.github.io/parameters/reference/principal_components.html # parameters - PCA
# https://easystats.github.io/parameters/reference/get_scores.html # parameters - get PCA scores
# https://drlee.io/secrets-of-pca-a-comprehensive-guide-to-principal-component-analysis-with-python-and-colab-6f7f3142e721 # PCA do's and don'ts

# if they are not already installed
# install.packages("factoextra") #Easy PCA plots and visualization
# install.packages("FactoMineR") #PCA and factor analysis functions  
# install.packages("ggplot2") #Enhanced data visualization
# install.packages("dplyr") #Data manipulation
# install.packages("gridExtra") #Arranging multiple plot grids
# install.packages("gt") # For creating beautiful tables
# install.packages("parameters")

# load packages
library(factoextra) #Easy PCA plots and visualization
library(FactoMineR) #PCA and factor analysis functions  
library(ggplot2) #Enhanced data visualization
library(dplyr) #Data manipulation
library(gridExtra) #Arranging multiple plot grids
library(gt) # For creating beautiful tables
library(parameters)

# prepare data ####

# set working directory
setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/environmental_data/pca/POC-C_8-10")

# import data
envs_pocc8910 <- read.csv("pca_env_data_POC-C_8-10.csv")

str(envs_pocc8910)
# 'data.frame':	24 obs. of  14 variables:
# $ sample       : chr  "1-425" "168" "312" "313" ...
# $ year         : int  2020 2019 2013 2013 2013 2013 2020 2020 2020 2019 ...
# $ lat          : num  -64 -61.4 -74.6 -74.6 -74.6 ...
# $ lon          : num  178 -34.7 -27.8 -27.8 -27.8 ...
# $ pop1         : int  88 48 48 48 48 48 88 88 88 58 ...
# $ cohort_2     : chr  "post2000" "post2000" "post2000" "post2000" ...
# $ current_m    : num  0.0507 0.0131 -0.0155 -0.0155 -0.0155 ...
# $ current_z    : num  0.1245 0.0204 -0.0477 -0.0477 -0.0477 ...
# $ sst          : num  -1.6 -1.8 -1.84 -1.84 -1.84 ...
# $ sal          : num  33.9 34.3 34.1 34.1 34.1 ...
# $ sea_ice_thick: num  0.824 1.045 1.427 1.427 1.427 ...
# $ sea_ice_conc : num  0.724 0.925 0.949 0.949 0.949 ...
# $ mxl0.01      : num  52.9 59.6 132.8 132.8 132.8 ...
# $ mxl0.03      : num  68.2 63.5 155 155 155 ...

# make year, pop1, and cohort_2 into factors
envs_pocc8910$year <- factor(envs_pocc8910$year, levels=unique(envs_pocc8910$year))
envs_pocc8910$pop1 <- factor(envs_pocc8910$pop1, levels=unique(envs_pocc8910$pop1))
envs_pocc8910$cohort_2 <- factor(envs_pocc8910$cohort_2, levels=unique(envs_pocc8910$cohort_2))

str(envs_pocc8910)
# 'data.frame':	24 obs. of  14 variables:
# $ sample       : chr  "1-425" "168" "312" "313" ...
# $ year         : Factor w/ 4 levels "2020","2019",..: 1 2 3 3 3 3 1 1 1 2 ...
# $ lat          : num  -64 -61.4 -74.6 -74.6 -74.6 ...
# $ lon          : num  178 -34.7 -27.8 -27.8 -27.8 ...
# $ pop1         : Factor w/ 3 levels "88","48","58": 1 2 2 2 2 2 1 1 1 3 ...
# $ cohort_2     : Factor w/ 2 levels "post2000","pre2000": 1 1 1 1 1 1 1 1 1 1 ...
# $ current_m    : num  0.0507 0.0131 -0.0155 -0.0155 -0.0155 ...
# $ current_z    : num  0.1245 0.0204 -0.0477 -0.0477 -0.0477 ...
# $ sst          : num  -1.6 -1.8 -1.84 -1.84 -1.84 ...
# $ sal          : num  33.9 34.3 34.1 34.1 34.1 ...
# $ sea_ice_thick: num  0.824 1.045 1.427 1.427 1.427 ...
# $ sea_ice_conc : num  0.724 0.925 0.949 0.949 0.949 ...
# $ mxl0.01      : num  52.9 59.6 132.8 132.8 132.8 ...
# $ mxl0.03      : num  68.2 63.5 155 155 155 ...

# create correlation matrix ####

# for correlation matrix, first column needs to be sample IDs (all unique)
# second column needs to be groupings, and the rest of the columns need to be quantitative variables
# subset data frame to have this structure before creating correlation matrix
# data frame should be called "data"
data <- subset(envs_pocc8910, select=-c(year, lat, lon, cohort_2))

str(data)
# 'data.frame':	24 obs. of  10 variables:
# $ sample       : chr  "1-425" "168" "312" "313" ...
# $ pop1         : Factor w/ 3 levels "88","48","58": 1 2 2 2 2 2 1 1 1 3 ...
# $ current_m    : num  0.0507 0.0131 -0.0155 -0.0155 -0.0155 ...
# $ current_z    : num  0.1245 0.0204 -0.0477 -0.0477 -0.0477 ...
# $ sst          : num  -1.6 -1.8 -1.84 -1.84 -1.84 ...
# $ sal          : num  33.9 34.3 34.1 34.1 34.1 ...
# $ sea_ice_thick: num  0.824 1.045 1.427 1.427 1.427 ...
# $ sea_ice_conc : num  0.724 0.925 0.949 0.949 0.949 ...
# $ mxl0.01      : num  52.9 59.6 132.8 132.8 132.8 ...
# $ mxl0.03      : num  68.2 63.5 155 155 155 ...

# code to create correlation matrix from website linked above

# Function for correlation matrix with stars. 
cor_mat_allan <-function(df,
                         type = "pearson",
                         digits = 3,
                         decimal.mark = ".",
                         use = "all",
                         show_significance = TRUE,
                         replace_diagonal = FALSE,
                         replacement = ""){
  
  # check arguments
  stopifnot({
    is.numeric(digits)
    digits >= 0
    use %in% c("all", "upper", "lower")
    is.logical(replace_diagonal)
    is.logical(show_significance)
    is.character(replacement)
  })
  # we need the Hmisc package for this
  require(Hmisc)
  
  # retain only numeric and boolean columns
  isNumericOrBoolean = vapply(df, function(x) is.numeric(x) | is.logical(x), logical(1))
  if (sum(!isNumericOrBoolean) > 0) {
    cat('Dropping non-numeric/-boolean column(s):', paste(names(isNumericOrBoolean)[!isNumericOrBoolean], collapse = ', '), '\n\n')
  }
  df = df[isNumericOrBoolean]
  
  # transform input data frame to matrix
  x <- as.matrix(df)
  
  # run correlation analysis using Hmisc package
  correlation_matrix <- Hmisc::rcorr(x, type = )
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value
  
  # transform correlations to specific character format
  Rformatted = formatC(R, format = 'f', digits = digits, decimal.mark = decimal.mark)
  
  # if there are any negative numbers, we want to put a space before the positives to align all
  if (sum(R < 0) > 0) {
    Rformatted = ifelse(R > 0, paste0(' ', Rformatted), Rformatted)
  }
  
  # add significance levels if desired
  if (show_significance) {
    # define notions for significance levels; spacing is important.
    stars <- ifelse(is.na(p), "   ", ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "*  ", "   "))))
    Rformatted = paste0(Rformatted, stars)
  }
  # build a new matrix that includes the formatted correlations and their significance stars
  Rnew <- matrix(Rformatted, ncol = ncol(x))
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep =" ")
  
  # replace undesired values
  if (use == 'upper') {
    Rnew[lower.tri(Rnew, diag = replace_diagonal)] <- replacement
  } else if (use == 'lower') {
    Rnew[upper.tri(Rnew, diag = replace_diagonal)] <- replacement
  } else if (replace_diagonal) {
    diag(Rnew) <- replacement
  }
  
  return(Rnew)
}


correlation_matrix= cor_mat_allan(data[,-c(1,2)])%>%as.data.frame()%>%gt( rownames_to_stub = TRUE) %>%tab_header(
  
  title = md("**Correlation Matrix**"),
  subtitle = "Package used : agricolae; PB-Perfect"
)%>%tab_source_note(
  source_note = "Source: PCA - Data analysis"
)%>%
  tab_options(
    heading.subtitle.font.size = 12,
    heading.align = "left",
    table.border.top.color = "red",
    column_labels.border.bottom.color = "red",
    column_labels.border.bottom.width= px(3)
  )%>%opt_stylize(style = 6, color = "cyan")%>%
  tab_options(table.width = pct(80))%>%
  tab_footnote(
    footnote = "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")

correlation_matrix

# create a covariance table ####
# Function for covariance matrix
covariance_table=cov(data[,-c(1,2)])%>%round(2)%>%as.data.frame()%>%gt( rownames_to_stub = TRUE) %>%tab_header(
  title = md("**Covariance Matrix**"),
  subtitle = "Package used : agricolae; PB-Perfect"
)%>%tab_source_note(
  source_note = "Source: https://medium.com/@victorallan"
)%>%
  tab_options(
    heading.subtitle.font.size = 12,
    heading.align = "left",
    table.border.top.color = "red",
    column_labels.border.bottom.color = "red",
    column_labels.border.bottom.width= px(3)
  )%>%opt_stylize(style = 6, color = "cyan")%>%
  tab_options(table.width = pct(80))

covariance_table

# extract eigenvalues from PCA ####
eigen_table=function(data){
  library(factoextra)
  library(FactoMineR)
  library(gt)
  data=data[,-c(1,2)]
  data=na.omit(data)
  
  pc=PCA(data,scale.unit = T, graph = F)
  x=get_eigenvalue(pc)
  x=t(x)
  rownames(x)=c("Eigen Value","Variance %", "Cumulative Variance %")
  colnames(x)=paste("PC",1:ncol(x),sep = "")
  x=round(x,2)
  
  y=eigen(cor(data))$vectors
  colnames(y)=paste("PC",1:ncol(x),sep = "")
  rownames(y)=colnames(data)
  
  x=rbind(x,y)
  
  x=as.data.frame(x)
  x=round(x,2)
  
  
  grp=c(rep("Eigen values and Variance %",3),rep("Eigenvectors",nrow(y)))
  
  x_out=cbind(grp,x)
  x=cbind(x,grp)
  
  
  x=x%>%gt(groupname_col = "grp", rownames_to_stub = TRUE)%>%
    tab_header(
      title = md("**Eigenvalues and Eigenvectors**"),
      subtitle = "Package used : agricolae; PB-Perfect")%>%
    tab_source_note(source_note = "Source: https://medium.com/@victorallan")%>%
    tab_options(
      heading.subtitle.font.size = 12,
      heading.align = "left",
      table.border.top.color = "red",
      column_labels.border.bottom.color = "red",
      column_labels.border.bottom.width= px(3)
    )%>%opt_stylize(style = 6, color = "cyan")%>%
    tab_options(table.width = pct(80))
  return(x)
  
}


eig_table=eigen_table(data)
eig_table

# create scree plot ####
pca_res=PCA(data[,-c(1,2)], graph = F)
screeplot=fviz_eig(pca_res, addlabels = TRUE,ncp = 200)
plot(screeplot)

ggsave("scree_plot_POC-C_8-10_10x10.pdf",
       width = 10, height = 10, units = "in",
       dpi = 300,
       device = cairo_pdf)

# calculate variable stats ####
var_stats=function(data){
  library(FactoMineR)
  library(dplyr)
  data=data[,-c(1,2)]
  
  
  res.pca=get_pca_var(PCA(data))
  res.pca=lapply(res.pca, as.data.frame)
  res.pca=lapply(res.pca, function(x){colnames(x)=c(paste("PC",1:ncol(x),sep = ""));x})
  res.pca=lapply(res.pca, round,2)
  
  cor_pca=res.pca$cor
  cos_pca=res.pca$cos2
  contrib_pca=res.pca$contrib
  
  grp=("Correation between PCs and Traits")
  cor_pca=cbind(cor_pca,grp)
  grp=("Quality of representation of Traits on PCs")
  cos_pca=cbind(cos_pca,grp)
  grp=c("Contributions(%) of the Traits to the Principal Components")
  contrib_pca=cbind(contrib_pca,grp)
  
  res=rbind(cor_pca, cos_pca, contrib_pca)
  x=res%>%gt(groupname_col = "grp", rownames_to_stub = TRUE)%>%
    tab_header(
      title = md("**PCA - Trait statistics**"),
      subtitle = "Package used : agricolae; PB-Perfect")%>%
    tab_source_note(source_note = "Source: https://medium.com/@victorallan")%>%
    tab_options(
      heading.subtitle.font.size = 12,
      heading.align = "left",
      table.border.top.color = "red",
      column_labels.border.bottom.color = "red",
      column_labels.border.bottom.width= px(3)
    )%>%opt_stylize(style = 6, color = "cyan")%>%
    tab_options(table.width = pct(80))
  
  return(x)
  
}

variable_stats = var_stats(data)
variable_stats

# visualize variable stats ####
varplot=function(data){
  library(FactoMineR)
  require(factoextra)
  library(dplyr)
  require(stringi)
  require(cowplot)
  require(gridExtra)
  
  
  data=data[,-c(1,2)]
  
  res=PCA(data,graph = F)
  var=get_pca_var(PCA(data, graph = F))
  res.pca=get_pca_var(PCA(data, graph = F))
  res.pca=lapply(res.pca, as.data.frame)
  res.pca=lapply(res.pca, function(x){colnames(x)=c(paste("PC",1:ncol(x),sep = ""));x})
  res.pca=lapply(res.pca, round,2)
  
  cor_pca=res.pca$cor
  cos_pca=res.pca$cos2
  contrib_pca=res.pca$contrib
  
  cos_contrib=fviz_pca_var(res, col.var = "contrib",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
  )
  
  cos_contrib$labels$x=gsub("Dim","PC",cos_contrib$labels$x)
  cos_contrib$labels$y=gsub("Dim","PC",cos_contrib$labels$y)
  
  
  cos_cor=fviz_pca_var(res, col.var = "cos2",
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                       repel = TRUE # Avoid text overlapping
  )
  
  cos_cor$labels$x=gsub("Dim","PC",cos_contrib$labels$x)
  cos_cor$labels$y=gsub("Dim","PC",cos_contrib$labels$y)
  
  
  
  temp1=fviz_cos2(res, choice = "var", axes = 1)
  temp2=fviz_cos2(res, choice = "var", axes = 2)
  temp3=fviz_cos2(res, choice = "var", axes = 3)
  
  # Contributions of variables to PC1
  temp4=fviz_contrib(res, choice = "var", axes = 1)
  # Contributions of variables to PC2
  temp5=fviz_contrib(res, choice = "var", axes = 2)
  temp6=fviz_contrib(res, choice = "var", axes = 3)
  
  temp1$labels$title="Cos2 of variables to PC-1"
  temp2$labels$title="Cos2 of variables to PC-2"
  temp3$labels$title="Cos2 of variables to PC-3"
  
  temp4$labels$title="Contribution of variables to PC-1"
  temp5$labels$title="Contribution of variables to PC-2"
  temp6$labels$title="Contribution of variables to PC-2"
  
  
  
  
  first_row=plot_grid(cos_contrib,cos_cor, nrow = 1)
  second_row=plot_grid(temp1,temp2,temp3,nrow = 1)
  third_row=plot_grid(temp4,temp5,temp6,nrow = 1)
  
  return(list(first_row,second_row,third_row,temp4,temp5))
  
}


variable_plot=varplot(data)

variable_plot[[1]]
variable_plot[[2]]
variable_plot[[3]]
variable_plot[[4]]
variable_plot[[5]]

# visualize individual contributions ####
ind_plot = fviz_contrib(pca_res, choice = "ind", axes = 1:2)
ind_plot$labels$title = "Contribution of Individuals to PC-1-2"
ind_plot

# PCA biplot ####
#Parameters you can adjust
pt_size=2
contrib_plot=T
ellipse=T

#To function for the biplot

biplot_allan=function(data,pca_res,var_plot,pt_size,contrib_plot=T,ellipse=T){
  require(ggpubr)
  grps=colnames(data)[2]
  grp_vec=as.factor(data[,2])
  data=data[,-c(1,2)]
  xplot <- var_plot[[4]]
  yplot <- var_plot[[5]]
  # Cleaning the plots
  
  xplot=xplot+xlab("")+theme_minimal()
  yplot <- yplot + ylab("Contribution (%)")+ xlab("")+rotate()
  
  
  
  y=fviz_pca_biplot(pca_res,
                    # Individuals
                    geom.ind = "point",
                    fill.ind = grp_vec, col.ind = "black",
                    pointshape = 21, pointsize = pt_size,
                    palette = "jco",
                    addEllipses = ellipse,
                    # Variables
                    alpha.var ="contrib", col.var = "contrib",
                    gradient.cols = "RdYlBu",
                    
                    legend.title = list(fill = grps, color = "Contrib",
                                        alpha = "Contrib")
  )
  
  
  
  y$labels$x=gsub("Dim","PC",y$labels$x)
  y$labels$y=gsub("Dim","PC",y$labels$y)
  
  z=ggarrange(xplot, NULL, y, yplot,
              ncol = 2, nrow = 2,  align = "hv",
              widths = c(2, 1), heights = c(1, 2),
              common.legend = TRUE)
  if(contrib_plot){
    return(z)
  }else{
    return(y)
  }
  
}



biplot=biplot_allan(data,pca_res,
                    pt_size=pt_size,
                    contrib_plot=contrib_plot,
                    ellipse=ellipse,
                    var_plot=variable_plot)

biplot

# option to save outputs ####
gtsave(correlation_matrix,"correlation_matrix_POC-C_8-10.docx") # For saving the correlation matrix
gtsave(covariance_table,"covariance_table_POC-C_8-10.docx") # For saving the covariance matrix
gtsave(eig_table,"eigen_table_POC-C_8-10.docx") # For saving the eigen table
gtsave(variable_stats,"variable_stats_POC-C_8-10.docx") # For saving the variable statistics

#You can also save the plot in different formats like
# "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only)
# by just changing the .pdf to a different extension

# You can also play with the Width and Height of the plot to get the optimum results.

ggsave(plot= screeplot,
       filename="scree_plot_POC-C_8-10.pdf",
       width = 10,
       height = 10) # For saving the Scree plot
ggsave(plot= variable_plot[[1]],
       filename="contrib_cos2_plots_POC-C_8-10.pdf",
       width = 16,
       height = 8) # For saving the variable plot
ggsave(plot= variable_plot[[2]],
       filename="cos2_barplots_POC-C_8-10.pdf",
       width = 16,
       height = 8) # For saving the variable plot
ggsave(plot= variable_plot[[3]],
       filename="contrib_barplos_POC-C_8-10.pdf",
       width = 16,
       height = 8) # For saving the variable plot
ggsave(plot= variable_plot[[4]],
       filename="contrib_barplot_PC1_POC-C_8-10.pdf",
       width = 10,
       height = 10) # For saving the variable plot
ggsave(plot= variable_plot[[5]],
       filename="contrib_barplot_PC2_POC-C_8-10.pdf",
       width = 10,
       height = 10) # For saving the variable plot
ggsave(plot= ind_plot,
       filename="individual_contrib_POC-C_8-10.pdf",
       width = 10,
       height = 10) # For saving the individual contribution plot
ggsave(plot= biplot,
       filename="pca_biplot_POC-C_8-10.pdf",
       width = 10,
       height = 10) # For saving the variable plot

# extract individual scores ####

# factoextra prcomp(...)  ####
data2 <- subset(envs_pocc8910, select=-c(sample, year, lat, lon, pop1, cohort_2))

str(data2)
# 'data.frame':	24 obs. of  8 variables:
# $ current_m    : num  0.0507 0.0131 -0.0155 -0.0155 -0.0155 ...
# $ current_z    : num  0.1245 0.0204 -0.0477 -0.0477 -0.0477 ...
# $ sst          : num  -1.6 -1.8 -1.84 -1.84 -1.84 ...
# $ sal          : num  33.9 34.3 34.1 34.1 34.1 ...
# $ sea_ice_thick: num  0.824 1.045 1.427 1.427 1.427 ...
# $ sea_ice_conc : num  0.724 0.925 0.949 0.949 0.949 ...
# $ mxl0.01      : num  52.9 59.6 132.8 132.8 132.8 ...
# $ mxl0.03      : num  68.2 63.5 155 155 155 ...

# conduct pca with prcomp with scaling
data2_pca <- prcomp(data2,  scale = TRUE)

var <- get_pca_var(data2_pca)
var
# Principal Component Analysis Results for variables
#   Name       Description                                    
# 1 "$coord"   "Coordinates for the variables"                
# 2 "$cor"     "Correlations between variables and dimensions"
# 3 "$cos2"    "Cos2 for the variables"                       
# 4 "$contrib" "contributions of the variables"  

# Coordinates of variables
head(var$coord)
#                 Dim.1       Dim.2       Dim.3       Dim.4       Dim.5       Dim.6        Dim.7        Dim.8
# current_m     -0.4094023 -0.84015325 -0.01012440  0.03828797 -0.30956993 -0.16959486  0.018669462  0.004393126
# current_z     -0.4671362 -0.70297992  0.01812055 -0.51384964  0.12847544  0.05901881  0.056939837 -0.001383790
# sst           -0.8950466 -0.31786390  0.10986123  0.11052769  0.16623415 -0.03414369 -0.211139861 -0.013736097
# sal           -0.2152437  0.10205939 -0.96854493 -0.04623657 -0.03406703  0.01846387 -0.038897857 -0.004736634
# sea_ice_thick  0.9085350  0.06616507  0.10548081 -0.28874815 -0.21116931  0.05062556 -0.168694363 -0.008451893
# sea_ice_conc   0.9298659 -0.08717284 -0.11535486 -0.09679399  0.21918442 -0.23860559 -0.007825476 -0.006234186

# Cos2 of variables
head(var$cos2)
#                 Dim.1       Dim.2        Dim.3       Dim.4       Dim.5        Dim.6        Dim.7        Dim.8
# current_m     0.16761023 0.705857490 0.0001025036 0.001465969 0.095833540 0.0287624155 3.485488e-04 1.929956e-05
# current_z     0.21821621 0.494180763 0.0003283542 0.264041452 0.016505938 0.0034832201 3.242145e-03 1.914874e-06
# sst           0.80110838 0.101037457 0.0120694901 0.012216370 0.027633791 0.0011657916 4.458004e-02 1.886803e-04
# sal           0.04632983 0.010416119 0.9380792751 0.002137820 0.001160563 0.0003409146 1.513043e-03 2.243570e-05
# sea_ice_thick 0.82543584 0.004377817 0.0111262023 0.083375492 0.044592479 0.0025629470 2.845779e-02 7.143449e-05
# sea_ice_conc  0.86465054 0.007599104 0.0133067438 0.009369077 0.048041810 0.0569326260 6.123807e-05 3.886507e-05

# Contribution of variables
head(var$contrib)
#                 Dim.1      Dim.2       Dim.3      Dim.4      Dim.5      Dim.6       Dim.7      Dim.8
# current_m      4.610195 28.2776750  0.01037234  0.3334663 39.8823586 26.2629751  0.43396508 0.18899316
# current_z      6.002135 19.7975983  0.03322618 60.0619278  6.8691579  3.1805299  4.03667363 0.01875163
# sst           22.034845  4.0477071  1.22131231  2.7788771 11.5001572  1.0644849 55.50494288 1.84767396
# sal            1.274323  0.4172848 94.92428889  0.4862934  0.4829831  0.3112893  1.88383367 0.21970418
# sea_ice_thick 22.703982  0.1753817  1.12586097 18.9655553 18.5577328  2.3402281 35.43172836 0.69953049
# sea_ice_conc  23.782600  0.3044311  1.34651007  2.1311989 19.9932164 51.9852076  0.07624523 0.38059067

ind <- get_pca_ind(data2_pca)
ind
# Principal Component Analysis Results for individuals
#   Name       Description                       
# 1 "$coord"   "Coordinates for the individuals" 
# 2 "$cos2"    "Cos2 for the individuals"        
# 3 "$contrib" "contributions of the individuals"

# Coordinates of individuals
head(ind$coord)
#     Dim.1      Dim.2      Dim.3      Dim.4       Dim.5       Dim.6       Dim.7        Dim.8
# 1 -1.9242806 -3.6901621  2.4768107 -0.7727372  0.67883974 -0.17745709 -0.28770189  0.022386473
# 2  0.5580213 -1.2326675 -1.3152217 -0.5039882 -0.07624246  0.21982871  0.37740223 -0.009793189
# 3  1.6484522  0.9207833  0.9570898 -0.6396569 -0.46472247  0.08676322 -0.05130875  0.055186578
# 4  1.6484522  0.9207833  0.9570898 -0.6396569 -0.46472247  0.08676322 -0.05130875  0.055186578
# 5  1.6484522  0.9207833  0.9570898 -0.6396569 -0.46472247  0.08676322 -0.05130875  0.055186578
# 6  1.6484522  0.9207833  0.9570898 -0.6396569 -0.46472247  0.08676322 -0.05130875  0.055186578

# Cos2 of individuals
head(ind$cos2)
#     Dim.1     Dim.2     Dim.3      Dim.4       Dim.5       Dim.6        Dim.7        Dim.8
# 1 0.15035480 0.5529315 0.2490956 0.02424622 0.018711776 0.001278696 0.0033609798 2.034941e-05
# 2 0.07762695 0.3787939 0.4312300 0.06332159 0.001449122 0.012047027 0.0355074900 2.390888e-05
# 3 0.53078352 0.1656073 0.1789246 0.07992063 0.042184417 0.001470402 0.0005142184 5.948831e-04
# 4 0.53078352 0.1656073 0.1789246 0.07992063 0.042184417 0.001470402 0.0005142184 5.948831e-04
# 5 0.53078352 0.1656073 0.1789246 0.07992063 0.042184417 0.001470402 0.0005142184 5.948831e-04
# 6 0.53078352 0.1656073 0.1789246 0.07992063 0.042184417 0.001470402 0.0005142184 5.948831e-04

# Contribution of individuals
head(ind$contrib)
#     Dim.1     Dim.2     Dim.3    Dim.4     Dim.5     Dim.6     Dim.7      Dim.8
# 1 4.2436962 22.730360 25.864984 5.659519 7.9907323 1.1981027 4.2940330 0.20448375
# 2 0.3568691  2.536339  7.293307 2.407447 0.1007966 1.8385530 7.3890547 0.03913233
# 3 3.1142983  1.415241  3.862175 3.878021 3.7448930 0.2864037 0.1365724 1.24266582
# 4 3.1142983  1.415241  3.862175 3.878021 3.7448930 0.2864037 0.1365724 1.24266582
# 5 3.1142983  1.415241  3.862175 3.878021 3.7448930 0.2864037 0.1365724 1.24266582
# 6 3.1142983  1.415241  3.862175 3.878021 3.7448930 0.2864037 0.1365724 1.24266582

# You can also use the function get_pca()
get_pca(data2_pca, "ind") # Results for individuals
# Principal Component Analysis Results for individuals
#   Name       Description                       
# 1 "$coord"   "Coordinates for the individuals" 
# 2 "$cos2"    "Cos2 for the individuals"        
# 3 "$contrib" "contributions of the individuals"

get_pca(data2_pca, "var") # Results for variable categories
# Principal Component Analysis Results for variables
#   Name       Description                                    
# 1 "$coord"   "Coordinates for the variables"                
# 2 "$cor"     "Correlations between variables and dimensions"
# 3 "$cos2"    "Cos2 for the variables"                       
# 4 "$contrib" "contributions of the variables"   

# FINAL SCORES TO USE FOR RDA ####

# get scores for all PCs generated with prcomp with scaling
data2_pca$x
# PC1        PC2        PC3         PC4         PC5         PC6         PC7          PC8
# 1  -1.9242806 -3.6901621  2.4768107 -0.77273720  0.67883974 -0.17745709 -0.28770189  0.022386473
# 2   0.5580213 -1.2326675 -1.3152217 -0.50398821 -0.07624246  0.21982871  0.37740223 -0.009793189
# 3   1.6484522  0.9207833  0.9570898 -0.63965692 -0.46472247  0.08676322 -0.05130875  0.055186578
# 4   1.6484522  0.9207833  0.9570898 -0.63965692 -0.46472247  0.08676322 -0.05130875  0.055186578
# 5   1.6484522  0.9207833  0.9570898 -0.63965692 -0.46472247  0.08676322 -0.05130875  0.055186578
# 6   1.6484522  0.9207833  0.9570898 -0.63965692 -0.46472247  0.08676322 -0.05130875  0.055186578
# 7   0.4206738 -0.8474015 -0.5121922  0.15234036  0.12701427 -0.05031290 -0.38779440 -0.109560203
# 8   0.4206738 -0.8474015 -0.5121922  0.15234036  0.12701427 -0.05031290 -0.38779440 -0.109560203
# 9   0.4206738 -0.8474015 -0.5121922  0.15234036  0.12701427 -0.05031290 -0.38779440 -0.109560203
# 10  1.5962896  0.8988142  0.2112992  0.68248711  0.26415277  0.46827903  0.01227464  0.011145715
# 11  1.4195217  0.7988403  0.1629379  0.74755965  0.39551858  0.43518276  0.15743163  0.027127451
# 12  0.5484575 -0.3579878 -0.4925126  1.05410376 -0.05414952 -0.28349234 -0.23740438  0.030174002
# 13  1.6652540  0.9913020  0.1276069  0.76958718  0.21807906  0.49291109 -0.09785641 -0.007526401
# 14 -0.1121485 -1.0424588 -0.3637642 -0.77414478  0.98722052  0.53337563  0.13499210 -0.084276621
# 15  0.5484575 -0.3579878 -0.4925126  1.05410376 -0.05414952 -0.28349234 -0.23740438  0.030174002
# 16  0.7153948 -1.4755195 -1.6242605 -0.45659641 -1.19213454 -0.19557281 -0.16084649 -0.022439293
# 17 -5.7724961 -0.9808033 -0.9855704  0.12307184 -0.40062607  0.56058662 -0.13001847  0.250107980
# 18 -0.8975032  4.1517982 -1.0623488 -0.69145252  1.00807996 -0.60400252 -0.24954232  0.160936100
# 19 -3.5205066  2.0934340  0.9369392  0.35603428 -0.34651952 -0.04359795  0.16064482 -0.171966691
# 20 -3.5205066  2.0934340  0.9369392  0.35603428 -0.34651952 -0.04359795  0.16064482 -0.171966691
# 21 -0.4865681  0.2070461 -1.1976513 -1.11695030  0.10942266 -0.29674494  0.28034947 -0.123287663
# 22  0.4948935 -0.6485154 -0.9441686 -0.07250404  0.19914099  0.01671699  0.62163906  0.006544674
# 23  0.4159448 -1.2947476  0.6668475  0.67349912  0.04386698 -0.49751854  0.43200690  0.080297224
# 24  0.4159448 -1.2947476  0.6668475  0.67349912  0.04386698 -0.49751854  0.43200690  0.080297224

# parameters principal_components(...) ####
# number of PCs not specified
data2_pca2 <- principal_components(data2)
data2_pca2
# Loadings from Principal Component Analysis (no rotation)

# Variable      |   PC1 |   PC2 | Complexity
# ------------------------------------------
# current_m     | -0.41 | -0.84 |       1.45
# current_z     | -0.47 | -0.70 |       1.74
# sst           | -0.90 | -0.32 |       1.25
# sal           | -0.22 |  0.10 |       1.43
# sea_ice_thick |  0.91 |  0.07 |       1.01
# sea_ice_conc  |  0.93 | -0.09 |       1.02
# mxl0.01       | -0.60 |  0.77 |       1.89
# mxl0.03       | -0.59 |  0.77 |       1.88
# 
# The 2 principal components accounted for 76.65% of the total variance of the original data (PC1 = 45.45%, PC2 = 31.20%).

pca2_scores <- get_scores(data2_pca2)
pca2_scores
# Component_1 Component_2
# 1     8.466503    30.31849
# 2     8.624807    30.77214
# 3     8.656482    71.92328
# 4     8.656482    71.92328
# 5     8.656482    71.92328
# 6     8.656482    71.92328
# 7     8.627867    37.37789
# 8     8.627867    37.37789
# 9     8.627867    37.37789
# 10    8.603977    47.14088
# 11    8.585268    46.35783
# 12    8.593889    38.91612
# 13    8.612398    46.72663
# 14    8.580943    41.17619
# 15    8.593889    38.91612
# 16    8.687157    31.81326
# 17    8.273970   109.52232
# 18    8.601398   165.19152
# 19    8.336503   151.29462
# 20    8.336503   151.29462
# 21    8.612487    85.28821
# 22    8.584961    40.66645
# 23    8.539264    32.03237
# 24    8.539264    32.03237

# all PCs
data2_pca3 <- principal_components(data2, n="all")
data2_pca3
# Loadings from Principal Component Analysis (no rotation)

# Variable      |   PC1 |   PC2 |   PC3 |   PC4 |   PC5 |   PC6 |       PC7 | Complexity
# --------------------------------------------------------------------------------------
# current_m     | -0.41 | -0.84 | -0.01 |  0.04 | -0.31 | -0.17 |      0.02 |       1.86
# current_z     | -0.47 | -0.70 |  0.02 | -0.51 |  0.13 |  0.06 |      0.06 |       2.76
# sst           | -0.90 | -0.32 |  0.11 |  0.11 |  0.17 | -0.03 |     -0.21 |       1.53
# sal           | -0.22 |  0.10 | -0.97 | -0.05 | -0.03 |  0.02 |     -0.04 |       1.13
# sea_ice_thick |  0.91 |  0.07 |  0.11 | -0.29 | -0.21 |  0.05 |     -0.17 |       1.45
# sea_ice_conc  |  0.93 | -0.09 | -0.12 | -0.10 |  0.22 | -0.24 | -7.83e-03 |       1.33
# mxl0.01       | -0.60 |  0.77 |  0.10 | -0.16 | -0.07 | -0.08 |      0.03 |       2.09
# mxl0.03       | -0.59 |  0.77 |  0.06 | -0.20 | -0.03 | -0.10 |     -0.03 |       2.10
# 
# The 7 principal components accounted for 99.87% of the total variance of the original data (PC1 = 45.45%, PC2 = 31.20%, PC3 = 12.35%, PC4 = 5.50%, PC5 = 3.00%, PC6 = 1.37%, PC7 = 1.00%).

pca3_scores <- get_scores(data2_pca3)
pca3_scores
# Component_1 Component_2 Component_3
# 1  -0.01864000    30.31849    33.92193
# 2   0.05720889    30.77214    34.32760
# 3   0.18017556    71.92328    34.08540
# 4   0.18017556    71.92328    34.08540
# 5   0.18017556    71.92328    34.08540
# 6   0.18017556    71.92328    34.08540
# 7   0.08831222    37.37789    34.24653
# 8   0.08831222    37.37789    34.24653
# 9   0.08831222    37.37789    34.24653
# 10  0.08713556    47.14088    34.15450
# 11  0.06125778    46.35783    34.15730
# 12  0.04654111    38.91612    34.23593
# 13  0.09476333    46.72663    34.16530
# 14  0.03081333    41.17619    34.23133
# 15  0.04654111    38.91612    34.23593
# 16  0.12532000    31.81326    34.37267
# 17 -0.42833937   109.52232    34.38090
# 18  0.01699778   165.19152    34.35460
# 19 -0.26845111   151.29462    34.15137
# 20 -0.26845111   151.29462    34.15137
# 21  0.03634889    85.28821    34.34090
# 22  0.01980333    40.66645    34.28043
# 23  0.02190778    32.03237    34.09133
# 24  0.02190778    32.03237    34.09133

# the scores are the same even if all PCs are used...

# princomp(...)  ####
# without scaling
data2_pca4 <- princomp(data2, cor = FALSE, scores = TRUE)
data2_pca4
# Call:
#   princomp(x = data2, cor = FALSE, scores = TRUE)
# 
# Standard deviations:
#   Comp.1       Comp.2       Comp.3       Comp.4       Comp.5       Comp.6       Comp.7       Comp.8 
# 114.18663768   9.71367909   0.30719411   0.11831234   0.09242296   0.04377803   0.02094705   0.01161439 
# 
# 8  variables and  24 observations.

data2_pca4$scores
#       Comp.1      Comp.2      Comp.3       Comp.4       Comp.5       Comp.6       Comp.7        Comp.8
# 1    93.21244  -6.4755549  0.40152564  0.276772111 -0.177225704  0.101615762  0.016550035 -0.0036196691
# 2    92.41094   1.6687145  0.06670873 -0.110171886  0.094743358  0.028945250 -0.032372574 -0.0018273975
# 3   -24.67905  -2.9210684 -0.45710719  0.152988648  0.032942235  0.002914857 -0.008781010  0.0036201965
# 4   -24.67905  -2.9210684 -0.45710719  0.152988648  0.032942235  0.002914857 -0.008781010  0.0036201965
# 5   -24.67905  -2.9210684 -0.45710719  0.152988648  0.032942235  0.002914857 -0.008781010  0.0036201965
# 6   -24.67905  -2.9210684 -0.45710719  0.152988648  0.032942235  0.002914857 -0.008781010  0.0036201965
# 7    73.68298   2.1779479  0.01182852 -0.053322235  0.032194313  0.017908676  0.035113568 -0.0050311901
# 8    73.68298   2.1779479  0.01182852 -0.053322235  0.032194313  0.017908676  0.035113568 -0.0050311901
# 9    73.68298   2.1779479  0.01182852 -0.053322235  0.032194313  0.017908676  0.035113568 -0.0050311901
# 10   45.86503   1.3508693 -0.11593132  0.019437921 -0.007480881 -0.066380676 -0.005135871 -0.0160635704
# 11   48.09210   1.3850183 -0.04567821 -0.008140086 -0.037715751 -0.068093179 -0.014353049 -0.0155852978
# 12   68.89165  -3.2791912  0.10688070 -0.070189730 -0.024305402 -0.051715660  0.026553378  0.0091055672
# 13   47.06497   1.7964944 -0.13351883  0.019592060  0.011425889 -0.070891655  0.003523223 -0.0172076815
# 14   62.84487   0.6910152  0.13241403 -0.053867341 -0.009675064  0.060604383 -0.017619109 -0.0298143540
# 15   68.89165  -3.2791912  0.10688070 -0.070189730 -0.024305402 -0.051715660  0.026553378  0.0091055672
# 16   89.35119   0.1540474 -0.10017628 -0.052502785  0.226768118  0.029936949  0.010174577  0.0221029325
# 17 -132.35746 -19.3367474  0.96735998  0.120655430  0.185041554 -0.030289795 -0.015949089  0.0003228935
# 18 -290.96474 -27.3059890 -0.29789749 -0.222152038 -0.127685195  0.012516994  0.010190622 -0.0016336258
# 19 -247.95844  21.9469061  0.19032050  0.040162911 -0.015431394 -0.005900876  0.007494939  0.0016537637
# 20 -247.95844  21.9469061  0.19032050  0.040162911 -0.015431394 -0.005900876  0.007494939  0.0016537637
# 21  -61.88237   5.9161258 -0.11846178 -0.192205291  0.028001292  0.087935825 -0.013016483  0.0028859973
# 22   64.50893   4.1708156  0.10265695 -0.155080663 -0.016349594  0.006957110 -0.041888461  0.0012113101
# 23   88.82746   1.9000953  0.16976969 -0.017135841 -0.159363155 -0.021504677 -0.019208558  0.0191612924
# 24   88.82746   1.9000953  0.16976969 -0.017135841 -0.159363155 -0.021504677 -0.019208558  0.0191612924

# with scaling
data2_pca5 <- princomp(data2, cor = TRUE, scores = TRUE)
data2_pca5
# Call:
#   princomp(x = data2, cor = TRUE, scores = TRUE)
# 
# Standard deviations:
#   Comp.1    Comp.2    Comp.3    Comp.4    Comp.5    Comp.6    Comp.7    Comp.8 
# 1.9067363 1.5799257 0.9941023 0.6630350 0.4901944 0.3309335 0.2834030 0.1010533 
# 
# 8  variables and  24 observations.

data2_pca5$scores
#       Comp.1     Comp.2     Comp.3      Comp.4      Comp.5      Comp.6      Comp.7       Comp.8
# 1   1.9656677  3.7695295 -2.5300816 -0.78935711 -0.69344012  0.18127381 -0.29388973  0.022867958
# 2  -0.5700231  1.2591796  1.3435093 -0.51482791  0.07788227 -0.22455675  0.38551934 -0.010003819
# 3  -1.6839069 -0.9405874 -0.9776747 -0.65341456  0.47471765 -0.08862931 -0.05241229  0.056373522
# 4  -1.6839069 -0.9405874 -0.9776747 -0.65341456  0.47471765 -0.08862931 -0.05241229  0.056373522
# 5  -1.6839069 -0.9405874 -0.9776747 -0.65341456  0.47471765 -0.08862931 -0.05241229  0.056373522
# 6  -1.6839069 -0.9405874 -0.9776747 -0.65341456  0.47471765 -0.08862931 -0.05241229  0.056373522
# 7  -0.4297215  0.8656273  0.5232083  0.15561687 -0.12974608  0.05139502 -0.39613502 -0.111916606
# 8  -0.4297215  0.8656273  0.5232083  0.15561687 -0.12974608  0.05139502 -0.39613502 -0.111916606
# 9  -0.4297215  0.8656273  0.5232083  0.15561687 -0.12974608  0.05139502 -0.39613502 -0.111916606
# 10 -1.6306223 -0.9181457 -0.2158438  0.69716593 -0.26983412 -0.47835070  0.01253864  0.011385436
# 11 -1.4500526 -0.8160216 -0.1664423  0.76363804 -0.40402533 -0.44454260  0.16081765  0.027710904
# 12 -0.5602537  0.3656873  0.5031055  1.07677525  0.05531416  0.28958965 -0.24251043  0.030822979
# 13 -1.7010700 -1.0126228 -0.1303514  0.78613934 -0.22276947 -0.50351254 -0.09996109 -0.007688278
# 14  0.1145605  1.0648799  0.3715880 -0.79079496 -1.00845349 -0.54484739  0.13789548 -0.086089228
# 15 -0.5602537  0.3656873  0.5031055  1.07677525  0.05531416  0.28958965 -0.24251043  0.030822979
# 16 -0.7307813  1.5072548  1.6591949 -0.46641681  1.21777477  0.19977916 -0.16430595 -0.022921913
# 17  5.8966500  1.0018983  1.0067679  0.12571884  0.40924267 -0.57264363 -0.13281489  0.255487262
# 18  0.9168066 -4.2410944  1.0851976 -0.70632417 -1.02976158  0.61699331 -0.25490944  0.164397487
# 19  3.5962251 -2.1384592 -0.9570908  0.36369180  0.35397241  0.04453564  0.16409995 -0.175665322
# 20  3.5962251 -2.1384592 -0.9570908  0.36369180  0.35397241  0.04453564  0.16409995 -0.175665322
# 21  0.4970331 -0.2114993  1.2234102 -1.14097348 -0.11177610  0.30312728  0.28637918 -0.125939314
# 22 -0.5055375  0.6624635  0.9644757 -0.07406345 -0.20342409 -0.01707653  0.63500917  0.006685436
# 23 -0.4248909  1.3225948 -0.6811900  0.68798463 -0.04481047  0.50821909  0.44129843  0.082024244
# 24 -0.4248909  1.3225948 -0.6811900  0.68798463 -0.04481047  0.50821909  0.44129843  0.082024244

# scores are quite different when scaling is done - but should scaling be done?