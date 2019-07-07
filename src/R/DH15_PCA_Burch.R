# -------------------------------------
# Merritt Burch does things with data!
# Double Haploid PCA
# Summer 2014 Data
# 9/21/2016

# Information:
# http://gastonsanchez.com/how-to/2012/06/17/PCA-in-R/
# http://www.sthda.com/english/wiki/principal-component-analysis-in-r-prcomp-vs-princomp-r-software-and-data-mining
# http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization
# http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining
# --------------------------------------

# Refresh Directory
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/DH15")

# Import Data
df.15 <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH15_raw data_edited.csv")

# Load packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'stats', 'ggbiplot', 'factoextra', 'dplyr')
pack.man(packages)


# Remove extraneous rows
df.15 <- df.15[, -1]
tbl_df(df.15)

# Convert treatments to factors
df.15[, 1:6] <- lapply(df.15[, 1:6], factor)

# Extract critical variables for PCA
df.pca <- df.15[, c(7:11, 14:15)]


# ----------------------------
# Principal Component Analysis
# ----------------------------

# prcomp uses singular value decomposistion (SVD) --> has bettter numerical accuracy
# Center = If TRUE, the data will be centered and scaled before the analysis
# Scale =  a logical value indicating whether the variables should be scaled 
#           to have unit variance before the analysis takes place
dh.pca <- prcomp(df.pca,
                 center = TRUE,
                 scale = TRUE)

# Look at list of dh.pca options
str(dh.pca)
print(dh.pca)

# Variance retained by each PC
# Eigenvalues, variances, and cumulative variances
eig.val <- get_eigenvalue(dh.pca)
(dh.eigenval <- eig.val)

# Variable graphing names
# All variable: coordinates, corr b/w variables and axes, squared cosine and contributions
(var <- get_pca_var(dh.pca))


# -------------------------------------------------------------
# Scree plots
# These plots show the importance of princpal components (PCs) 
# -------------------------------------------------------------

# Percent Variance vs. Dimensions
fviz_screeplot(dh.pca, ncp=10) +
  theme(text = element_text(size = 30)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file="DH15_scree_variance.png", width = 10, height = 6)

# Eigenvalues vs. Dimensions
fviz_screeplot(dh.pca, ncp=10, choice="eigenvalue")
ggsave(file="DH15_scree_eigen.png")


# ---------------------
# Correlation Circles
# ---------------------

# Correlation Circle
fviz_pca_var(dh.pca)
ggsave(file="DH15_corr circle.png")


# Cos2 in correlation circles
# Cos2 : quality of representation for variables on the factor map
# cos2 calculated as the squared coordinates : var.cos2 = var.coord * var.coord
# The contribution of a variable to a given principal component is 
#   (in percentage) : (var.cos2 * 100) / (total cos2 of the component)
fviz_pca_var(dh.pca, col.var="contrib")+
  scale_color_gradient2(low="lightgreen", mid="coral", 
                        high="darkred", midpoint=50) 
ggsave(file="DH15_Cos2 contribution_corr circle.png", height = 10, width = 6)


# -----------------------------------------
# Contributions of variables to dimensions
# -----------------------------------------

# Contributions on axis 1
fviz_contrib(dh.pca, choice="var", axes = 1 )+
  labs(title = "Contributions to Dim 1") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 25))
ggsave(file="DH15_percent contrib_dim 1.png", width = 10, height = 6)

# Variable contributions on axes 1 + 2
fviz_contrib(dh.pca, choice="var", axes = 1:2)+
  labs(title = "Contributions to Dim 1+2")
ggsave(file="DH15_percent contrib_dim 12.png")

# Variable contributions on axes 1 : 5
fviz_contrib(dh.pca, choice="var", axes = 1:5)+
  labs(title = "Contributions to Dim 1:5")
ggsave(file="DH15_percent contrib_dim 15.png")


# -------------------------------
# Factor map/Graph of individuals
# -------------------------------

# Repel = TRUE to avoid overplotting
# cos2 = the quality of the individuals on the factor map
fviz_pca_ind(dh.pca, repel = TRUE, col.ind = "cos2")+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.45)
ggsave(file="DH15_percent contrib_individual FM_.png")

# Default individuals factor map
fviz_pca_ind(dh.pca)

# Individuals factor map with no names
fviz_pca_biplot(dh.pca, label ="var")

# Individuals and variables, biplot is heated
fviz_pca_biplot(dh.pca, label = "var", col.var="contrib", labelsize = 4.5, pointsize = 1)+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=30) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 25))
ggsave(file="DH15_heated biplot and individuals.png", width = 10, height = 6)

# Biplot with no heating, correct name and dot sizes
fviz_pca_biplot(dh.pca, axes = c(1, 2), geom = c("point"),
                label = "all", invisible = "none", labelsize = 4.5, pointsize = 1.2,
                col.ind = "black", col.ind.sup = "blue", col.var = "red") 

# Individual Factor map with heating
fviz_pca_ind(dh.pca, col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.45) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 25))
ggsave(file="DH15_heated individuals FM_.png", width = 10, height = 6)


# -------------------------------------------
# Look at outliers in individual factor maps
# -------------------------------------------

# Quadrant 1
q1 <- df.15[c(64,91,112,128,147,152,155,166,222,240), 4]

# Quadrant 2
q2 <- df.15[c(22,88,89), 4]

# Quadrant 3
q3 <- df.15[c(1,49,97,146,151,153,195),4]

# Quadrant 4
q4 <- df.15[c(9,15), 4]

# ----------------------------------
# Correlation Circles with variables
# ----------------------------------

# Block
fviz_pca_ind(dh.pca, geom = "point",
             habillage=df.15$Block, addEllipses=TRUE,
             ellipse.level= 0.95) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 25))
ggsave(file = "DH15_individual FM and corr circles by block.png", width = 10, height = 6)

# Lineage
fviz_pca_ind(dh.pca, geom = "point",
             habillage=df.15$lin, addEllipses=TRUE,
             ellipse.level= 0.95)  +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 25))
ggsave(file = "DH15_individual FM and corr circles by lin.png", width = 10, height = 6)


# ------------------------------------
# Correlation circles with biplots
# Similar to above, but with biplots
# ------------------------------------

# Block
fviz_pca_biplot(dh.pca,  label="var", habillage=df.15$Block,
                addEllipses=TRUE, ellipse.level=0.95)

# Family
fviz_pca_biplot(dh.pca,  label="var", habillage=df.15$Family,
                addEllipses=TRUE, ellipse.level=0.95)

# Lineage-Generation
fviz_pca_biplot(dh.pca,  label="var", habillage=df.15$lin,
                addEllipses=TRUE, ellipse.level=0.95)

