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
setwd("~/SDSU Masters/Vivek Data/Data Sets")

# Import Data
df.dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH raw data -2014_edited_V2.csv")


# Packages Required
library(stats) # comes with R, fasest and easiest way to do PCA
library(ggbiplot) # For PCA graphing
library(ggplot2) # For graphing
library(factoextra) # Best ggplot graphing and PCA helper ever


# Remove extraneous rows
df.dh <- df.dh[, c(-1, -2)]
tbl_df(df.dh)

# Reorder this dataframe (it's a mess!)
tmp <- c(1, 2, 3, 5, 4, 20, 19, 6:18)
df.dh <- df.dh[, tmp]

# Convert treatments to factors
df.dh[, 1:7] <- lapply(df.dh[, 1:7], factor)

# Extract critical variables for PCA
df.pca <- df.dh[, 8:20]


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
fviz_screeplot(dh.pca, ncp=10)

# Eigenvalues vs. Dimensions
fviz_screeplot(dh.pca, ncp=10, choice="eigenvalue")


# ---------------------
# Correlation Circles
# ---------------------

# Correlation Circle
fviz_pca_var(dh.pca)


# Cos2 in correlation circles
# Cos2 : quality of representation for variables on the factor map
# cos2 calculated as the squared coordinates : var.cos2 = var.coord * var.coord
# The contribution of a variable to a given principal component is 
#   (in percentage) : (var.cos2 * 100) / (total cos2 of the component)
fviz_pca_var(dh.pca, col.var="contrib")+
  scale_color_gradient2(low="cornflowerblue", mid="green3", 
                        high="red", midpoint=35) 


# -----------------------------------------
# Contributions of variables to dimensions
# -----------------------------------------

# Contributions on axis 1
fviz_contrib(dh.pca, choice="var", axes = 1 )+
  labs(title = "Contributions to Dim 1") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))


# Variable contributions on axes 1 + 2
fviz_contrib(dh.pca, choice="var", axes = 1:2)+
  labs(title = "Contributions to Dim 1+2")

# Variable contributions on axes 1 : 5
fviz_contrib(dh.pca, choice="var", axes = 1:5)+
  labs(title = "Contributions to Dim 1:5")


# -------------------------------
# Factor map/Graph of individuals
# -------------------------------

# Repel = TRUE to avoid overplotting
# cos2 = the quality of the individuals on the factor map
fviz_pca_ind(dh.pca, repel = TRUE, col.ind = "cos2")+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=0.45)

# Default individuals factor map
fviz_pca_ind(dh.pca)

# Individuals factor map with no names
fviz_pca_biplot(dh.pca, label ="var")

# Individuals and variables, biplot is heated
fviz_pca_biplot(dh.pca, label = "var", col.var="contrib", labelsize = 4.5, pointsize = 1)+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=30)
  
# Biplot with no heating, correct name and dot sizes
fviz_pca_biplot(dh.pca, axes = c(1, 2), geom = c("point"),
                label = "all", invisible = "none", labelsize = 4.5, pointsize = 1.2,
                col.ind = "black", col.ind.sup = "blue", col.var = "red") 

# Biplot with heating
fviz_pca_ind(dh.pca, col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
    high="red", midpoint=0.45) 


# -------------------------------------------
# Look at outliers in individual factor maps
# -------------------------------------------

# Quadrant 1
q1 <- dh[c(17,28,56,58,67,80,119,251,276,325,366,426),c(3:4,8,14,21:22)] 

# Quadrant 2
q2 <- dh[c(63,91,95,138,151,198,316,305,321,328,433,434,439),c(3:4,19:20,22)]

# Quadrant 3
q3 <- dh[c(55,57,107,169,165,184,190,191,292,194,246,262,263,276,377,367,369,378,384,399,404,422),c(3:4, 10:13, 22)]

# Quadrant 4
q4 <- dh[c(10,106,132,134,144,149,199,205,206,219,244,272,300,303,304,307,334,360,387,441), c(3:4,9,14:18,22)]


# ----------------------------------
# Correlation Circles with variables
# ----------------------------------

# Block
fviz_pca_ind(dh.pca, geom = "point",
                  habillage=dh$Block, addEllipses=TRUE,
                  ellipse.level= 0.95)

# Family
fviz_pca_ind(dh.pca, geom = "point",
                  habillage=dh$Family, addEllipses=TRUE,
                  ellipse.level= 0.95)

# Lineage Generation names
# Specify shapes is a hard one, hot mess graph
fviz_pca_ind(dh.pca, geom = "point",
                  habillage=dh$newtreat2, addEllipses=TRUE,
                  ellipse.level= 0.95)

# Year
fviz_pca_ind(dh.pca, geom = "point",
             habillage=dh$year, addEllipses=TRUE,
             ellipse.level= 0.95)


# ------------------------------------
# Correlation circles with biplots
# Similar to above, but with biplots
# ------------------------------------

# Block
fviz_pca_biplot(dh.pca,  label="var", habillage=dh$Block,
                addEllipses=TRUE, ellipse.level=0.95)

# Family
fviz_pca_biplot(dh.pca,  label="var", habillage=dh$Family,
                addEllipses=TRUE, ellipse.level=0.95)

# Lineage-Generation
fviz_pca_biplot(dh.pca,  label="var", habillage=dh$newtreat2,
                addEllipses=TRUE, ellipse.level=0.95)

# Year
fviz_pca_biplot(dh.pca,  label="var", habillage=dh$year,
                addEllipses=TRUE, ellipse.level=0.95)





