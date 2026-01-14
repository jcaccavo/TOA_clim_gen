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
setwd(".../pca/SG-B_8-10")

# import data
envs_sgb8910 <- read.csv("pca_env_data_SG-B_8-10.csv")

str(envs_sgb8910)
# 'data.frame':	24 obs. of  14 variables:
# $ sample       : chr  "1-425" "168" "312" "313" ...
# $ year         : int  2005 2005 2000 2000 2000 2000 2005 2005 2005 2010 ...
# $ lat          : num  -64.8 -60.7 -60.7 -60.7 -60.7 -60.7 -64.8 -64.8 -64.8 -62.2 ...
# $ lon          : num  167 -27.8 -27.8 -27.8 -27.8 -27.8 167 167 167 76.6 ...
# $ pop1         : int  88 48 48 48 48 48 88 88 88 58 ...
# $ cohort_2     : chr  "post2000" "post2000" "post2000" "post2000" ...
# $ current_m    : num  0.0265 0.0549 0.0498 0.0498 0.0498 ...
# $ current_z    : num  0.2177 0.0358 0.0385 0.0385 0.0385 ...
# $ sst          : num  -3.18 -1.82 -1.78 -1.78 -1.78 ...
# $ sal          : num  67.9 34.2 34.4 34.4 34.4 ...
# $ sea_ice_thick: num  1.017 0.899 0.762 0.762 0.762 ...
# $ sea_ice_conc : num  1.125 0.868 0.843 0.843 0.843 ...
# $ mxl0.01      : num  138.5 65.3 69.7 69.7 69.7 ...
# $ mxl0.03      : num  163.3 71.5 83.2 83.2 83.2 ...

# make year, pop1, and cohort_2 into factors
envs_sgb8910$year <- factor(envs_sgb8910$year, levels=unique(envs_sgb8910$year))
envs_sgb8910$pop1 <- factor(envs_sgb8910$pop1, levels=unique(envs_sgb8910$pop1))
envs_sgb8910$cohort_2 <- factor(envs_sgb8910$cohort_2, levels=unique(envs_sgb8910$cohort_2))

str(envs_sgb8910)
# 'data.frame':	24 obs. of  14 variables:
# $ sample       : chr  "1-425" "168" "312" "313" ...
# $ year         : Factor w/ 4 levels "2005","2000",..: 1 1 2 2 2 2 1 1 1 3 ...
# $ lat          : num  -64.8 -60.7 -60.7 -60.7 -60.7 -60.7 -64.8 -64.8 -64.8 -62.2 ...
# $ lon          : num  167 -27.8 -27.8 -27.8 -27.8 -27.8 167 167 167 76.6 ...
# $ pop1         : Factor w/ 3 levels "88","48","58": 1 2 2 2 2 2 1 1 1 3 ...
# $ cohort_2     : Factor w/ 2 levels "post2000","pre2000": 1 1 1 1 1 1 1 1 1 1 ...
# $ current_m    : num  0.0265 0.0549 0.0498 0.0498 0.0498 ...
# $ current_z    : num  0.2177 0.0358 0.0385 0.0385 0.0385 ...
# $ sst          : num  -3.18 -1.82 -1.78 -1.78 -1.78 ...
# $ sal          : num  67.9 34.2 34.4 34.4 34.4 ...
# $ sea_ice_thick: num  1.017 0.899 0.762 0.762 0.762 ...
# $ sea_ice_conc : num  1.125 0.868 0.843 0.843 0.843 ...
# $ mxl0.01      : num  138.5 65.3 69.7 69.7 69.7 ...
# $ mxl0.03      : num  163.3 71.5 83.2 83.2 83.2 ...

# create correlation matrix ####

# for correlation matrix, first column needs to be sample IDs (all unique)
# second column needs to be groupings, and the rest of the columns need to be quantitative variables
# subset data frame to have this structure before creating correlation matrix
# data frame should be called "data"
data <- subset(envs_sgb8910, select=-c(year, lat, lon, cohort_2))

str(data)
# 'data.frame':	24 obs. of  10 variables:
# $ sample       : chr  "1-425" "168" "312" "313" ...
# $ pop1         : Factor w/ 3 levels "88","48","58": 1 2 2 2 2 2 1 1 1 3 ...
# $ current_m    : num  0.0265 0.0549 0.0498 0.0498 0.0498 ...
# $ current_z    : num  0.2177 0.0358 0.0385 0.0385 0.0385 ...
# $ sst          : num  -3.18 -1.82 -1.78 -1.78 -1.78 ...
# $ sal          : num  67.9 34.2 34.4 34.4 34.4 ...
# $ sea_ice_thick: num  1.017 0.899 0.762 0.762 0.762 ...
# $ sea_ice_conc : num  1.125 0.868 0.843 0.843 0.843 ...
# $ mxl0.01      : num  138.5 65.3 69.7 69.7 69.7 ...
# $ mxl0.03      : num  163.3 71.5 83.2 83.2 83.2 ...

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

ggsave("scree_plot_SG-B_8-10_10x10.pdf",
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
gtsave(correlation_matrix,"Correlation_matrix.docx") # For saving the correlation matrix
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
       width = 10,
       height = 10) # For saving the variable plot
ggsave(plot= variable_plot[[2]],
       filename="varplot2.pdf",
       width = 10,
       height = 10) # For saving the variable plot
ggsave(plot= variable_plot[[3]],
       filename="varplot3.pdf",
       width = 10,
       height = 10) # For saving the variable plot
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
data2 <- subset(envs_sgb8910, select=-c(sample, year, lat, lon, pop1, cohort_2))

str(data2)
# 'data.frame':	24 obs. of  8 variables:
# $ current_m    : num  0.0265 0.0549 0.0498 0.0498 0.0498 ...
# $ current_z    : num  0.2177 0.0358 0.0385 0.0385 0.0385 ...
# $ sst          : num  -3.18 -1.82 -1.78 -1.78 -1.78 ...
# $ sal          : num  67.9 34.2 34.4 34.4 34.4 ...
# $ sea_ice_thick: num  1.017 0.899 0.762 0.762 0.762 ...
# $ sea_ice_conc : num  1.125 0.868 0.843 0.843 0.843 ...
# $ mxl0.01      : num  138.5 65.3 69.7 69.7 69.7 ...
# $ mxl0.03      : num  163.3 71.5 83.2 83.2 83.2 ...

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
#                 Dim.1      Dim.2        Dim.3       Dim.4       Dim.5        Dim.6         Dim.7         Dim.8
# current_m      0.698728048 -0.2331926 -0.674154507  0.05181766 -0.01490957 -0.002940179 -0.0001827478  1.197049e-17
# current_z     -0.890861391 -0.3560173 -0.248764362 -0.13037541  0.02093421  0.017259064  0.0002869563 -5.864581e-17
# sst            0.500587424  0.8531988 -0.004309695  0.13935094  0.03983300  0.020977540 -0.0002018655  1.619451e-17
# sal           -0.977196433 -0.1858706 -0.094180307 -0.03503090  0.02093903  0.001846036 -0.0005487393  1.086444e-16
# sea_ice_thick  0.001628094 -0.9831653  0.044838646  0.15956037  0.07668154 -0.005753226  0.0004906177 -1.326595e-17
# sea_ice_conc   0.016577322 -0.9788803  0.144929181  0.12676427 -0.06515818  0.014105058 -0.0005757823  7.755638e-18

# Cos2 of variables
head(var$cos2)
#                 Dim.1      Dim.2        Dim.3       Dim.4        Dim.5        Dim.6        Dim.7        Dim.8
# current_m     4.882209e-01 0.05437877 4.544843e-01 0.002685070 0.0002222953 8.644653e-06 3.339678e-08 1.432927e-34
# current_z     7.936340e-01 0.12674833 6.188371e-02 0.016997748 0.0004382412 2.978753e-04 8.234393e-08 3.439331e-33
# sst           2.505878e-01 0.72794821 1.857347e-05 0.019418685 0.0015866676 4.400572e-04 4.074968e-08 2.622622e-34
# sal           9.549129e-01 0.03454789 8.869930e-03 0.001227164 0.0004384431 3.407849e-06 3.011148e-07 1.180360e-32
# sea_ice_thick 2.650690e-06 0.96661393 2.010504e-03 0.025459513 0.0058800578 3.309961e-05 2.407057e-07 1.759854e-34
# sea_ice_conc  2.748076e-04 0.95820667 2.100447e-02 0.016069179 0.0042455885 1.989527e-04 3.315253e-07 6.014993e-35

# Contribution of variables
head(var$contrib)
#                 Dim.1     Dim.2        Dim.3      Dim.4     Dim.5      Dim.6     Dim.7      Dim.8
# current_m     1.196698e+01  1.690214 80.589786763  2.1633031  1.602015  0.8298699 0.3531195  0.8047159
# current_z     1.945308e+01  3.939622 10.973304927 13.6947181  3.158273 28.5954480 0.8706604 19.3148956
# sst           6.142256e+00 22.626263  0.003293473 15.6452152 11.434638 42.2446339 0.4308652  1.4728353
# sal           2.340625e+01  1.073826  1.572828344  0.9886993  3.159728  0.3271468 3.1838257 66.2876989
# sea_ice_thick 6.497210e-05 30.044529  0.356505387 20.5121796 42.375816  3.1774978 2.5450924  0.9883145
# sea_ice_conc  6.735918e-03 29.783212  3.724541329 12.9465906 30.596685 19.0990697 3.5053696  0.3377953

# FINAL SCORES TO USE FOR RDA ####

# get scores for all PCs generated with prcomp with scaling
data2_pca$x
# PC1        PC2         PC3        PC4          PC5           PC6           PC7           PC8
# 1  -1.315406 -2.2751804  0.03795461 -0.1897538  0.003365114 -0.0004608202 -7.751727e-05 -2.544622e-16
# 2   2.443966 -0.7385953 -0.42540457  0.4764610  0.211488656 -0.0584421625  1.555902e-03  3.851086e-16
# 3   2.257708 -0.3269004 -0.24368411  0.2746121 -0.036290868  0.0523632888 -8.082511e-04 -1.214306e-17
# 4   2.257708 -0.3269004 -0.24368411  0.2746121 -0.036290868  0.0523632888 -8.082511e-04 -1.214306e-17
# 5   2.257708 -0.3269004 -0.24368411  0.2746121 -0.036290868  0.0523632888 -8.082511e-04 -1.214306e-17
# 6   2.257708 -0.3269004 -0.24368411  0.2746121 -0.036290868  0.0523632888 -8.082511e-04 -1.214306e-17
# 7  -1.315406 -2.2751804  0.03795461 -0.1897538  0.003365114 -0.0004608202 -7.751727e-05 -2.544622e-16
# 8  -1.315406 -2.2751804  0.03795461 -0.1897538  0.003365114 -0.0004608202 -7.751727e-05 -2.544622e-16
# 9  -1.315406 -2.2751804  0.03795461 -0.1897538  0.003365114 -0.0004608202 -7.751727e-05 -2.544622e-16
# 10  1.925526  2.1398823  0.09471647 -0.7160777  0.111421763 -0.0070436087 -5.432196e-03 -4.857226e-17
# 11  1.925526  2.1398823  0.09471647 -0.7160777  0.111421763 -0.0070436087 -5.432196e-03 -4.857226e-17
# 12 -1.315406 -2.2751804  0.03795461 -0.1897538  0.003365114 -0.0004608202 -7.751727e-05 -2.544622e-16
# 13  1.862587  1.8481631  0.40816668 -0.6710264  0.018777927  0.0189914612  1.223763e-02  8.326673e-17
# 14 -1.315406 -2.2751804  0.03795461 -0.1897538  0.003365114 -0.0004608202 -7.751727e-05 -2.544622e-16
# 15 -1.315406 -2.2751804  0.03795461 -0.1897538  0.003365114 -0.0004608202 -7.751727e-05 -2.544622e-16
# 16  2.443966 -0.7385953 -0.42540457  0.4764610  0.211488656 -0.0584421625  1.555902e-03  3.851086e-16
# 17  2.276070  0.3936561 -1.20249689 -0.1349845 -0.428969836 -0.0703298345 -1.702671e-04  6.765422e-16
# 18  1.153771  0.3407687  3.24590443  0.3630048 -0.123021553 -0.0273952691 -1.083249e-03 -3.625572e-16
# 19 -2.309066  1.9747837 -0.18019065  0.1920113  0.001500049  0.0005796285  9.068270e-05 -2.324529e-16
# 20 -2.309066  1.9747837 -0.18019065  0.1920113  0.001500049  0.0005796285  9.068270e-05 -2.324529e-16
# 21 -2.309066  1.9747837 -0.18019065  0.1920113  0.001500049  0.0005796285  9.068270e-05 -2.324529e-16
# 22 -2.309066  1.9747837 -0.18019065  0.1920113  0.001500049  0.0005796285  9.068270e-05 -2.324529e-16
# 23 -2.309066  1.9747837 -0.18019065  0.1920113  0.001500049  0.0005796285  9.068270e-05 -2.324529e-16
# 24 -2.309066  1.9747837 -0.18019065  0.1920113  0.001500049  0.0005796285  9.068270e-05 -2.324529e-16
