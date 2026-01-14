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
setwd(".../pca/POC-B_8-10")

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
