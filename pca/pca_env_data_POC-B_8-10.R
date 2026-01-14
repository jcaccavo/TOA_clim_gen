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
setwd("/Users/JMAC/Library/CloudStorage/Dropbox/Research/ClimGenAT/GEA/environmental_data/pca/POC-B_8-10")

# import data
envs_pocb8910 <- read.csv("pca_env_data_POC-B_8-10.csv")

str(envs_pocb8910)
# 'data.frame':	24 obs. of  14 variables:
# $ sample       : chr  "1-425" "168" "312" "313" ...
# $ year         : int  2005 2005 2000 2000 2000 2000 2005 2005 2005 2010 ...
# $ lat          : num  -64 -61.4 -74.6 -74.6 -74.6 ...
# $ lon          : num  178 -34.7 -27.8 -27.8 -27.8 ...
# $ pop1         : int  88 48 48 48 48 48 88 88 88 58 ...
# $ cohort_2     : chr  "post2000" "post2000" "post2000" "post2000" ...
# $ current_m    : num  -0.01235 0.0226 -0.00549 -0.00549 -0.00549 ...
# $ current_z    : num  0.1284 0.0412 -0.0672 -0.0672 -0.0672 ...
# $ sst          : num  -1.69 -1.84 -1.82 -1.82 -1.82 ...
# $ sal          : num  33.9 34.1 34.2 34.2 34.2 ...
# $ sea_ice_thick: num  0.549 1.298 1.121 1.121 1.121 ...
# $ sea_ice_conc : num  0.619 0.889 0.934 0.934 0.934 ...
# $ mxl0.01      : num  67.7 49.3 142.7 142.7 142.7 ...
# $ mxl0.03      : num  75.8 51 185.2 185.2 185.2 ...

# make year, pop1, and cohort_2 into factors
envs_pocb8910$year <- factor(envs_pocb8910$year, levels=unique(envs_pocb8910$year))
envs_pocb8910$pop1 <- factor(envs_pocb8910$pop1, levels=unique(envs_pocb8910$pop1))
envs_pocb8910$cohort_2 <- factor(envs_pocb8910$cohort_2, levels=unique(envs_pocb8910$cohort_2))

str(envs_pocb8910)
# 'data.frame':	24 obs. of  14 variables:
# $ sample       : chr  "1-425" "168" "312" "313" ...
# $ year         : Factor w/ 4 levels "2005","2000",..: 1 1 2 2 2 2 1 1 1 3 ...
# $ lat          : num  -64 -61.4 -74.6 -74.6 -74.6 ...
# $ lon          : num  178 -34.7 -27.8 -27.8 -27.8 ...
# $ pop1         : Factor w/ 3 levels "88","48","58": 1 2 2 2 2 2 1 1 1 3 ...
# $ cohort_2     : Factor w/ 2 levels "post2000","pre2000": 1 1 1 1 1 1 1 1 1 1 ...
# $ current_m    : num  -0.01235 0.0226 -0.00549 -0.00549 -0.00549 ...
# $ current_z    : num  0.1284 0.0412 -0.0672 -0.0672 -0.0672 ...
# $ sst          : num  -1.69 -1.84 -1.82 -1.82 -1.82 ...
# $ sal          : num  33.9 34.1 34.2 34.2 34.2 ...
# $ sea_ice_thick: num  0.549 1.298 1.121 1.121 1.121 ...
# $ sea_ice_conc : num  0.619 0.889 0.934 0.934 0.934 ...
# $ mxl0.01      : num  67.7 49.3 142.7 142.7 142.7 ...
# $ mxl0.03      : num  75.8 51 185.2 185.2 185.2 ...

# create correlation matrix ####

# for correlation matrix, first column needs to be sample IDs (all unique)
# second column needs to be groupings, and the rest of the columns need to be quantitative variables
# subset data frame to have this structure before creating correlation matrix
# data frame should be called "data"
data <- subset(envs_pocb8910, select=-c(year, lat, lon, cohort_2))

str(data)
# 'data.frame':	24 obs. of  10 variables:
# $ sample       : chr  "1-425" "168" "312" "313" ...
# $ pop1         : Factor w/ 3 levels "88","48","58": 1 2 2 2 2 2 1 1 1 3 ...
# $ current_m    : num  -0.01235 0.0226 -0.00549 -0.00549 -0.00549 ...
# $ current_z    : num  0.1284 0.0412 -0.0672 -0.0672 -0.0672 ...
# $ sst          : num  -1.69 -1.84 -1.82 -1.82 -1.82 ...
# $ sal          : num  33.9 34.1 34.2 34.2 34.2 ...
# $ sea_ice_thick: num  0.549 1.298 1.121 1.121 1.121 ...
# $ sea_ice_conc : num  0.619 0.889 0.934 0.934 0.934 ...
# $ mxl0.01      : num  67.7 49.3 142.7 142.7 142.7 ...
# $ mxl0.03      : num  75.8 51 185.2 185.2 185.2 ...

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

ggsave("scree_plot_POC-B_8-10_10x10.pdf",
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
gtsave(correlation_matrix,"correlation_matrix.docx") # For saving the correlation matrix
gtsave(covariance_table,"covariance_table.docx") # For saving the covariance matrix
gtsave(eig_table,"eigen_table.docx") # For saving the eigen table
gtsave(variable_stats,"variable_stats.docx") # For saving the variable statistics


#You can also save the plot in different formats like
# "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only)
# by just changing the .pdf to a different extension


# You can also play with the Width and Height of the plot to get the optimum results.

ggsave(plot= screeplot,
       filename="Scree_plot.pdf",
       width = 10,
       height = 10) # For saving the Scree plot



ggsave(plot= variable_plot[[1]],
       filename="varplot1.pdf",
       width = 16,
       height = 8) # For saving the variable plot
ggsave(plot= variable_plot[[2]],
       filename="varplot2.pdf",
       width = 16,
       height = 8) # For saving the variable plot
ggsave(plot= variable_plot[[3]],
       filename="varplot3.pdf",
       width = 16,
       height = 8) # For saving the variable plot
ggsave(plot= variable_plot[[4]],
       filename="varplot4.pdf",
       width = 10,
       height = 10) # For saving the variable plot
ggsave(plot= variable_plot[[5]],
       filename="varplot5.pdf",
       width = 10,
       height = 10) # For saving the variable plot



ggsave(plot= ind_plot,
       filename="ind_contrib.pdf",
       width = 10,
       height = 10) # For saving the individual contribution plot


ggsave(plot= biplot,
       filename="biplot.pdf",
       width = 10,
       height = 10) # For saving the variable plot

# extract individual scores ####

# factoextra prcomp(...)  ####
data2 <- subset(envs_pocb8910, select=-c(sample, year, lat, lon, pop1, cohort_2))

str(data2)
# 'data.frame':	24 obs. of  8 variables:
# $ current_m    : num  -0.01235 0.0226 -0.00549 -0.00549 -0.00549 ...
# $ current_z    : num  0.1284 0.0412 -0.0672 -0.0672 -0.0672 ...
# $ sst          : num  -1.69 -1.84 -1.82 -1.82 -1.82 ...
# $ sal          : num  33.9 34.1 34.2 34.2 34.2 ...
# $ sea_ice_thick: num  0.549 1.298 1.121 1.121 1.121 ...
# $ sea_ice_conc : num  0.619 0.889 0.934 0.934 0.934 ...
# $ mxl0.01      : num  67.7 49.3 142.7 142.7 142.7 ...
# $ mxl0.03      : num  75.8 51 185.2 185.2 185.2 ...

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
#                 Dim.1      Dim.2       Dim.3        Dim.4      Dim.5       Dim.6       Dim.7         Dim.8
# current_m      0.06250138 -0.2854882  0.90335666  0.234439310 -0.1868180  0.09067571  0.02125875 -8.198443e-05
# current_z      0.36248491 -0.4967157 -0.18576255  0.738917281  0.1749405 -0.09592659 -0.03956130  7.732165e-04
# sst           -0.79820363 -0.5260959  0.04643514 -0.162303922 -0.1222486 -0.05112709 -0.20006858  3.026572e-03
# sal           -0.63314036  0.4962973  0.40449781  0.003609679  0.3699120 -0.22848217  0.01229600  7.295683e-04
# sea_ice_thick  0.56915549  0.6268233  0.01234708  0.110938451 -0.4840846 -0.18377340 -0.05082480  1.004118e-03
# sea_ice_conc   0.65881939  0.6512689  0.13233829  0.042437520  0.2738775  0.15687432 -0.15123216 -1.358159e-03

# Cos2 of variables
head(var$cos2)
#                   Dim.1      Dim.2        Dim.3        Dim.4      Dim.5       Dim.6        Dim.7        Dim.8
# current_m     0.003906422 0.08150353 0.8160532554 5.496179e-02 0.03490098 0.008222084 0.0004519344 6.721447e-09
# current_z     0.131395309 0.24672644 0.0345077232 5.459987e-01 0.03060418 0.009201910 0.0015650963 5.978638e-07
# sst           0.637129037 0.27677688 0.0021562223 2.634256e-02 0.01494472 0.002613979 0.0400274366 9.160136e-06
# sal           0.400866721 0.24631105 0.1636184808 1.302978e-05 0.13683490 0.052204101 0.0001511917 5.322699e-07
# sea_ice_thick 0.323937966 0.39290747 0.0001524503 1.230734e-02 0.23433794 0.033772662 0.0025831599 1.008254e-06
# sea_ice_conc  0.434042982 0.42415121 0.0175134219 1.800943e-03 0.07500888 0.024609553 0.0228711676 1.844595e-06

# Contribution of variables
head(var$contrib)
#                   Dim.1     Dim.2      Dim.3        Dim.4     Dim.5     Dim.6      Dim.7       Dim.8
# current_m      0.1174511  4.044545 76.5809136  6.690935313  6.346950  5.552203  0.6668915 0.000110781
# current_z      3.9505504 12.243596  3.2383095 66.468765021  5.565551  6.213859  2.3095152 0.009853823
# sst           19.1560141 13.734824  0.2023464  3.206889477  2.717786  1.765166 59.0659986 0.150974774
# sal           12.0525170 12.222983 15.3544547  0.001586219 24.884238 35.252344  0.2231042 0.008772722
# sea_ice_thick  9.7395659 19.497709  0.0143064  1.498270245 42.615746 22.805976  3.8118085 0.016617752
# sea_ice_conc  13.0499993 21.048154  1.6435127  0.219243110 13.640810 16.618319 33.7495595 0.030402102

ind <- get_pca_ind(data2_pca)
ind
# Principal Component Analysis Results for individuals
#   Name       Description                       
# 1 "$coord"   "Coordinates for the individuals" 
# 2 "$cos2"    "Cos2 for the individuals"        
# 3 "$contrib" "contributions of the individuals"

# Coordinates of individuals
head(ind$coord)
#     Dim.1      Dim.2      Dim.3      Dim.4      Dim.5      Dim.6      Dim.7        Dim.8
# 1 1.2942393 -3.1135182 -1.8429308  1.1870577 -0.1932879  0.7677291 -0.2433881 -0.001555433
# 2 2.1731220 -0.2685826  0.3789990  0.8549136 -0.7109312 -0.4879533  0.0278441 -0.052719923
# 3 0.6763469  0.9124812 -0.2113816 -0.2476028 -0.5554805  0.2676188  0.1360421  0.083450156
# 4 0.6763469  0.9124812 -0.2113816 -0.2476028 -0.5554805  0.2676188  0.1360421  0.083450156
# 5 0.6763469  0.9124812 -0.2113816 -0.2476028 -0.5554805  0.2676188  0.1360421  0.083450156
# 6 0.6763469  0.9124812 -0.2113816 -0.2476028 -0.5554805  0.2676188  0.1360421  0.083450156

# Cos2 of individuals
head(ind$cos2)
#     Dim.1     Dim.2      Dim.3      Dim.4       Dim.5      Dim.6       Dim.7        Dim.8
# 1 0.09934756 0.5749510 0.20144017 0.08357409 0.002215834 0.03495780 0.003513392 1.434930e-07
# 2 0.73602239 0.0112429 0.02238715 0.11391143 0.078773097 0.03710901 0.000120834 4.331839e-04
# 3 0.25389551 0.4621297 0.02479994 0.03402729 0.171259088 0.03975111 0.010272189 3.865186e-03
# 4 0.25389551 0.4621297 0.02479994 0.03402729 0.171259088 0.03975111 0.010272189 3.865186e-03
# 5 0.25389551 0.4621297 0.02479994 0.03402729 0.171259088 0.03975111 0.010272189 3.865186e-03
# 6 0.25389551 0.4621297 0.02479994 0.03402729 0.171259088 0.03975111 0.010272189 3.865186e-03

# Contribution of individuals
head(ind$contrib)
#     Dim.1     Dim.2      Dim.3     Dim.4     Dim.5     Dim.6      Dim.7       Dim.8
# 1 2.0984358 20.044023 13.2803299 7.1475691 0.2830906 16.583958 3.64222812 0.001661476
# 2 5.9160894  0.149155  0.5616514 3.7073126 3.8297582  6.699289 0.04766887 1.908713658
# 3 0.5730671  1.721591  0.1747131 0.3109754 2.3380501  2.015140 1.13792978 4.782396153
# 4 0.5730671  1.721591  0.1747131 0.3109754 2.3380501  2.015140 1.13792978 4.782396153
# 5 0.5730671  1.721591  0.1747131 0.3109754 2.3380501  2.015140 1.13792978 4.782396153
# 6 0.5730671  1.721591  0.1747131 0.3109754 2.3380501  2.015140 1.13792978 4.782396153

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
#       PC1          PC2        PC3         PC4        PC5           PC6         PC7          PC8
# 1   1.29423931 -3.113518230 -1.8429308  1.18705768 -0.1932879  0.7677291379 -0.24338808 -0.001555433
# 2   2.17312201 -0.268582593  0.3789990  0.85491364 -0.7109312 -0.4879533190  0.02784410 -0.052719923
# 3   0.67634686  0.912481160 -0.2113816 -0.24760284 -0.5554805  0.2676187860  0.13604214  0.083450156
# 4   0.67634686  0.912481160 -0.2113816 -0.24760284 -0.5554805  0.2676187860  0.13604214  0.083450156
# 5   0.67634686  0.912481160 -0.2113816 -0.24760284 -0.5554805  0.2676187860  0.13604214  0.083450156
# 6   0.67634686  0.912481160 -0.2113816 -0.24760284 -0.5554805  0.2676187860  0.13604214  0.083450156
# 7   0.79793099 -0.123545934  0.4254513  0.05277647  0.8181943 -0.3745950549 -0.05225385 -0.049888744
# 8   0.79793099 -0.123545934  0.4254513  0.05277647  0.8181943 -0.3745950549 -0.05225385 -0.049888744
# 9   0.79793099 -0.123545934  0.4254513  0.05277647  0.8181943 -0.3745950549 -0.05225385 -0.049888744
# 10  0.64271035  0.890637200 -0.6463978 -0.58211329 -0.4762387  0.0728123142 -0.05444915 -0.057851637
# 11  0.63967306  0.880142773 -0.6242465 -0.57648516 -0.4928256  0.0699550880 -0.02467201 -0.065226447
# 12  0.46256521  0.087593131  0.1882574 -0.53607461  0.7716029 -0.0007439936  0.09522942 -0.066873198
# 13 -0.09481801  0.805641103 -0.1208155 -0.49144486  0.5221917 -0.0964492995  0.19789978 -0.025089775
# 14  1.60952376 -1.474104159 -1.4998891  0.91097195  0.5602140  0.0958539166 -0.34356199 -0.044809484
# 15  0.46256521  0.087593131  0.1882574 -0.53607461  0.7716029 -0.0007439936  0.09522942 -0.066873198
# 16  2.51575924  0.161691035  1.2167070  0.80467376 -1.9509128 -0.5520472287  0.01339734 -0.037048804
# 17 -1.62014464 -3.769645479 -0.4249305 -0.11164365  0.1998452 -0.4662139275  0.78904006  0.101428170
# 18 -0.38097431  2.033056281 -2.3776876 -1.62385207  0.4620604 -0.2397247779 -0.05908762  0.028968152
# 19 -4.53952513 -0.895469516  0.1499047 -0.55690640 -0.8618549 -0.0943214452 -0.27261627 -0.055645682
# 20 -4.53952513 -0.895469516  0.1499047 -0.55690640 -0.8618549 -0.0943214452 -0.27261627 -0.055645682
# 21 -3.48111247  2.924923497 -0.3030563  3.12807148  0.5318929  0.1401917505  0.25401698 -0.015754583
# 22 -0.28004158  0.003899379  0.9565658  0.35742738  0.6644365 -0.4748333393 -0.63028583  0.251212856
# 23  0.01840136 -0.368837439  2.0902654 -0.41976643  0.4156993  0.7070602914  0.02030654 -0.010324862
# 24  0.01840136 -0.368837439  2.0902654 -0.41976643  0.4156993  0.7070602914  0.02030654 -0.010324862

# parameters principal_components(...) ####
# number of PCs not specified
data2_pca2 <- principal_components(data2)
data2_pca2
# Loadings from Principal Component Analysis (no rotation)

# Variable      |   PC1 |   PC2 |   PC3 | Complexity
# --------------------------------------------------
# current_m     |  0.06 | -0.29 |  0.90 |       1.21
# current_z     |  0.36 | -0.50 | -0.19 |       2.15
# sst           | -0.80 | -0.53 |  0.05 |       1.74
# sal           | -0.63 |  0.50 |  0.40 |       2.65
# sea_ice_thick |  0.57 |  0.63 |  0.01 |       1.98
# sea_ice_conc  |  0.66 |  0.65 |  0.13 |       2.08
# mxl0.01       | -0.85 |  0.39 | -0.12 |       1.46
# mxl0.03       | -0.82 |  0.44 | -0.13 |       1.59
# 
# The 3 principal components accounted for 80.08% of the total variance of the original data (PC1 = 41.58%, PC2 = 25.19%, PC3 = 13.32%).

pca2_scores <- get_scores(data2_pca2)
pca2_scores
# Component_1 Component_2  Component_3
# 1     35.25784   0.3387844 -0.012347943
# 2     26.69766   0.6697000  0.022604267
# 3     72.24784   0.5266746 -0.005491928
# 4     72.24784   0.5266746 -0.005491928
# 5     72.24784   0.5266746 -0.005491928
# 6     72.24784   0.5266746 -0.005491928
# 7     37.04495   0.4149434  0.004484833
# 8     37.04495   0.4149434  0.004484833
# 9     37.04495   0.4149434  0.004484833
# 10    66.45435   0.5250909 -0.021735747
# 11    66.44099   0.5252153 -0.020877983
# 12    44.59424   0.3509650 -0.004510843
# 13    71.49434   0.3926625 -0.015647373
# 14    31.71893   0.4245110 -0.025646867
# 15    44.59424   0.3509650 -0.004510843
# 16    27.43696   0.8303593  0.050649900
# 17    45.75355   0.1041261  0.006463963
# 18    89.60828   0.4315962 -0.092040187
# 19   155.38221   0.1924249  0.005819276
# 20   155.38221   0.1924249  0.005819276
# 21   245.14263   0.4710700 -0.006936221
# 22    65.08722   0.4185264  0.017849621
# 23    49.57599   0.2757984  0.055709610
# 24    49.57599   0.2757984  0.055709610

# all PCs
data2_pca3 <- principal_components(data2, n="all")
data2_pca3
# Loadings from Principal Component Analysis (no rotation)

# Variable      |   PC1 |   PC2 |   PC3 |      PC4 |   PC5 |   PC6 |      PC7 | Complexity
# ----------------------------------------------------------------------------------------
# current_m     |  0.06 | -0.29 |  0.90 |     0.23 | -0.19 |  0.09 |     0.02 |       1.48
# current_z     |  0.36 | -0.50 | -0.19 |     0.74 |  0.17 | -0.10 |    -0.04 |       2.64
# sst           | -0.80 | -0.53 |  0.05 |    -0.16 | -0.12 | -0.05 |    -0.20 |       2.06
# sal           | -0.63 |  0.50 |  0.40 | 3.61e-03 |  0.37 | -0.23 |     0.01 |       3.71
# sea_ice_thick |  0.57 |  0.63 |  0.01 |     0.11 | -0.48 | -0.18 |    -0.05 |       3.17
# sea_ice_conc  |  0.66 |  0.65 |  0.13 |     0.04 |  0.27 |  0.16 |    -0.15 |       2.66
# mxl0.01       | -0.85 |  0.39 | -0.12 |     0.30 | -0.11 |  0.08 | 1.85e-03 |       1.81
# mxl0.03       | -0.82 |  0.44 | -0.13 |     0.30 | -0.10 |  0.10 |     0.01 |       1.98
# 
# The 7 principal components accounted for 99.92% of the total variance of the original data (PC1 = 41.58%, PC2 = 25.19%, PC3 = 13.32%, PC4 = 10.27%, PC5 = 6.87%, PC6 = 1.85%, PC7 = 0.85%).

pca3_scores <- get_scores(data2_pca3)
pca3_scores
#     Component_1 Component_2  Component_3  Component_4
# 1     35.25784   0.5492083 -0.012347943  0.128360433
# 2     26.69766   1.2981737  0.022604267  0.041226333
# 3     72.24784   1.1205430 -0.005491928 -0.067193883
# 4     72.24784   1.1205430 -0.005491928 -0.067193883
# 5     72.24784   1.1205430 -0.005491928 -0.067193883
# 6     72.24784   1.1205430 -0.005491928 -0.067193883
# 7     37.04495   0.8303190  0.004484833 -0.000432200
# 8     37.04495   0.8303190  0.004484833 -0.000432200
# 9     37.04495   0.8303190  0.004484833 -0.000432200
# 10    66.45435   1.1237817 -0.021735747 -0.073599793
# 11    66.44099   1.1242983 -0.020877983 -0.073867807
# 12    44.59424   0.7481088 -0.004510843 -0.046178753
# 13    71.49434   0.8522901 -0.015647373 -0.066965173
# 14    31.71893   0.7490793 -0.025646867  0.099942600
# 15    44.59424   0.7481088 -0.004510843 -0.046178753
# 16    27.43696   1.6550900  0.050649900  0.005628570
# 17    45.75355   0.1720212  0.006463963  0.036231006
# 18    89.60828   0.9913310 -0.092040187 -0.128138600
# 19   155.38221   0.4937593  0.005819276 -0.108909530
# 20   155.38221   0.4937593  0.005819276 -0.108909530
# 21   245.14263   0.9114724 -0.006936221  0.030667680
# 22    65.08722   0.8392055  0.017849621 -0.002152766
# 23    49.57599   0.6253310  0.055709610 -0.073734153
# 24    49.57599   0.6253310  0.055709610 -0.073734153

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
# 173.99781389  10.84201728   0.33618541   0.17178635   0.09727009   0.05453855   0.03741222   0.02311175 
# 
# 8  variables and  24 observations.

data2_pca4$scores
#       Comp.1       Comp.2       Comp.3      Comp.4      Comp.5       Comp.6       Comp.7        Comp.8
# 1   121.313496   3.80927088  0.397070001 -0.14197590 -0.33178955 -0.014967300  0.044263725 -7.650691e-05
# 2   152.039926   6.86809284 -0.418050424 -0.17750325  0.01789652  0.074399890  0.007476788  8.250651e-04
# 3   -10.057426 -14.64429979 -0.096437262 -0.11074506 -0.01383404 -0.036906151 -0.013189571  1.349797e-02
# 4   -10.057426 -14.64429979 -0.096437262 -0.11074506 -0.01383404 -0.036906151 -0.013189571  1.349797e-02
# 5   -10.057426 -14.64429979 -0.096437262 -0.11074506 -0.01383404 -0.036906151 -0.013189571  1.349797e-02
# 6   -10.057426 -14.64429979 -0.096437262 -0.11074506 -0.01383404 -0.036906151 -0.013189571  1.349797e-02
# 7   115.591379   6.22039609 -0.061629030  0.16695801  0.04661076  0.046677391 -0.001900850 -1.469488e-02
# 8   115.591379   6.22039609 -0.061629030  0.16695801  0.04661076  0.046677391 -0.001900850 -1.469488e-02
# 9   115.591379   6.22039609 -0.061629030  0.16695801  0.04661076  0.046677391 -0.001900850 -1.469488e-02
# 10   11.318838   4.46174263 -0.326608780 -0.05096439 -0.02053671 -0.055260253 -0.012173581 -6.793782e-03
# 11   11.400850   5.36974450 -0.335474640 -0.05269968 -0.02004914 -0.053472068 -0.015822823 -5.042608e-03
# 12   88.871608   6.61806456 -0.009684687  0.19853928  0.01432283 -0.006355990 -0.026522646 -9.059068e-04
# 13   -6.613996   0.31064898 -0.048994869  0.11966895  0.03482362 -0.003822559 -0.042465550 -3.586122e-03
# 14  134.213812   6.35998228  0.055230940  0.03952104 -0.20141493  0.022447252  0.032279805 -3.098734e-02
# 15   88.871608   6.61806456 -0.009684687  0.19853928  0.01432283 -0.006355990 -0.026522646 -9.059068e-04
# 16  149.330216   5.56420222 -0.687740196 -0.38843637  0.09503579  0.062935163  0.028088134  1.977998e-02
# 17   83.894870  -0.48186925  1.072898561 -0.26123601  0.04869399  0.089988998 -0.086834628  1.144948e-02
# 18  -71.109327  -8.34340125 -0.109109576  0.03443188  0.01060511 -0.089990313 -0.057791552 -5.686747e-02
# 19 -302.860247  18.47163577  0.306947541 -0.14937853  0.09368503 -0.053221906  0.036905624 -6.133830e-03
# 20 -302.860247  18.47163577  0.306947541 -0.14937853  0.09368503 -0.053221906  0.036905624 -6.133830e-03
# 21 -621.306382  -6.16227554 -0.315401419  0.16379620 -0.12761877  0.126378066 -0.015868382  1.326886e-02
# 22   14.930742 -28.06303556  0.354441455  0.05022014  0.13985667  0.047094169  0.099832837 -3.756844e-02
# 23   71.009901   0.02175374  0.168924690  0.25448103  0.02699277 -0.039491412  0.028355054  4.988555e-02
# 24   71.009901   0.02175374  0.168924690  0.25448103  0.02699277 -0.039491412  0.028355054  4.988555e-02

# with scaling
data2_pca5 <- princomp(data2, cor = TRUE, scores = TRUE)
data2_pca5
# Call:
#   princomp(x = data2, cor = TRUE, scores = TRUE)
# 
# Standard deviations:
#   Comp.1     Comp.2     Comp.3     Comp.4     Comp.5     Comp.6     Comp.7     Comp.8 
# 1.82373245 1.41955868 1.03228348 0.90633138 0.74154286 0.38482060 0.26032155 0.07789306 
# 
# 8  variables and  24 observations.

data2_pca5$scores
#       Comp.1       Comp.2     Comp.3      Comp.4     Comp.5        Comp.6      Comp.7       Comp.8
# 1   1.32207560  3.180483273 -1.8825683  1.21258872  0.1974451  0.7842413310 -0.24862283  0.001588887
# 2   2.21986116  0.274359223  0.3871504  0.87330098  0.7262218 -0.4984481394  0.02844297  0.053853815
# 3   0.69089362 -0.932106656 -0.2159280 -0.25292824  0.5674276  0.2733746872  0.13896812 -0.085244988
# 4   0.69089362 -0.932106656 -0.2159280 -0.25292824  0.5674276  0.2733746872  0.13896812 -0.085244988
# 5   0.69089362 -0.932106656 -0.2159280 -0.25292824  0.5674276  0.2733746872  0.13896812 -0.085244988
# 6   0.69089362 -0.932106656 -0.2159280 -0.25292824  0.5674276  0.2733746872  0.13896812 -0.085244988
# 7   0.81509276  0.126203139  0.4346019  0.05391157 -0.8357919 -0.3826517842 -0.05337772  0.050961743
# 8   0.81509276  0.126203139  0.4346019  0.05391157 -0.8357919 -0.3826517842 -0.05337772  0.050961743
# 9   0.81509276  0.126203139  0.4346019  0.05391157 -0.8357919 -0.3826517842 -0.05337772  0.050961743
# 10  0.65653366 -0.909792879 -0.6603004 -0.59463328  0.4864816  0.0743783496 -0.05562024  0.059095900
# 11  0.65343105 -0.899072740 -0.6376727 -0.58888410  0.5034253  0.0714596706 -0.02520265  0.066629327
# 12  0.47251399 -0.089477070  0.1923065 -0.54760442 -0.7881984 -0.0007599953  0.09727760  0.068311496
# 13 -0.09685734 -0.822968701 -0.1234140 -0.50201477 -0.5334229 -0.0985237153  0.20215618  0.025629402
# 14  1.64414113  1.505808951 -1.5321485  0.93056498 -0.5722630  0.0979155270 -0.35095127  0.045773239
# 15  0.47251399 -0.089477070  0.1923065 -0.54760442 -0.7881984 -0.0007599953  0.09727760  0.068311496
# 16  2.56986778 -0.165168660  1.2428757  0.82198055  1.9928727 -0.5639205704  0.01368549  0.037845643
# 17 -1.65499044  3.850722399 -0.4340698 -0.11404486 -0.2041434 -0.4762411806  0.80601061 -0.103609670
# 18 -0.38916824 -2.076782923 -2.4288265 -1.65877761 -0.4719983 -0.2448807394 -0.06035847 -0.029591194
# 19 -4.63716050  0.914729128  0.1531288 -0.56888425  0.8803916 -0.0963500955 -0.27847966  0.056842500
# 20 -4.63716050  0.914729128  0.1531288 -0.56888425  0.8803916 -0.0963500955 -0.27847966  0.056842500
# 21 -3.55598367 -2.987832274 -0.3095744  3.19534953 -0.5433328  0.1432069719  0.25948034  0.016093430
# 22 -0.28606467 -0.003983247  0.9771394  0.36511487 -0.6787270 -0.4850459774 -0.64384191 -0.256615901
# 23  0.01879713  0.376770334  2.1352224 -0.42879470 -0.4246401  0.7222676289  0.02074329  0.010546927
# 24  0.01879713  0.376770334  2.1352224 -0.42879470 -0.4246401  0.7222676289  0.02074329  0.010546927

# scores are quite different when scaling is done - but should scaling be done?