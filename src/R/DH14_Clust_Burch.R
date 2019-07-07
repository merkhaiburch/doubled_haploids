# ---------------------------------------------------------------------------
# Merritt Burch does things with data!
# Double Haploid Cluster Analysis
# Summer 2014 Data
# 9/21/2016
# For 3D Biplots
# http://planspace.org/2013/02/03/pca-3d-visualization-and-clustering-in-r/
# For clustering
# http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-
#     data-analyses-and-elegant-visualization
# ---------------------------------------------------------------------------

# Refresh environment
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/Data Sets")

# Import data
dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH raw data -2014_edited_V2.csv")

# Check Data
str(dh)
summary(dh)


# Packages Required
library(stats) # comes with R, fasest and easiest way to do PCA
library(ggbiplot) # For PCA graphing
library(ggplot2) # For graphing
library(factoextra) # Best ggplot graphing and PCA helper ever
library(rgl)

# Convert to factors (instead of integers)
# Hastags b/c I want to correlate everything
#dh$Family <- factor(dh$Family)
#dh$Block <- factor(dh$Block)
#dh$newtreat <- factor(dh$newtreat)
#dh$Treatment <- factor(dh$Treatment)
#dh$year <- factor(dh$year)




# Log transform data
#   Authors in this blog do not do this but it is standard practice elsewhere
#   [Row,Column]
#   All traits in clumns in 1-20, add column 5 for newtreatment names
subdh1 <- dh[,3:5]
subdh2 <- dh[, 8:21]
dh2 <- merge(subdh2, subdh2)
log.dh2 <- log(dh2)
#dh2.newtreat <- dh[, 5]


# prcomp uses singular value decomposistion (SVD)-->has bettter numerical accuracy
# Center and scaling data possible
# Center = If TRUE, the data will be centered and scaled before the analysis
# Scale =  a logical value indicating whether the variables should be scaled 
#           to have unit variance before the analysis takes place
# Run PCA
dh2.pca <- prcomp(log.dh2,
                  center = TRUE,
                  scale = TRUE) 


# Simple 3D plot
plot3d(dh2.pca$x[,1:5], col=dh$Family)
plot3d(dh2.pca$x[,1:5], col=dh$Block)
plot3d(dh2.pca$x[,1:5], col=dh$plantno)
plot3d(dh2.pca$x[,1:5], col=dh$year)
plot3d(dh2.pca$x[,1:5], col=dh$newtreat2)


# Biplot with 3D view, names are very small
text3d(dh2.pca$x[,1:5],texts=rownames(dh))
text3d(dh2.pca$rotation[,1:5], texts=rownames(dh2.pca$rotation), col="red")
coords <- NULL
for (i in 1:nrow(dh2.pca$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0),dh2.pca$rotation[i,1:3]))
}
lines3d(coords, col="red", lwd=4)


# Biplot with 3D view, only colmn names
text3d(dh2.pca$x[,1:5],texts=newtreat2(dh))
text3d(dh2.pca$rotation[,1:5], texts=rownames(dh2.pca$rotation), col="red")
coords <- NULL
for (i in 1:nrow(dh2.pca$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0),dh2.pca$rotation[i,1:3]))
}
lines3d(coords, col="red", lwd=4)



# ---------------------------------------
# Hericherial clustering


# To determine the optimal number of clusters
dh2_data <- scale(dh2)
fviz_nbclust(dh2_data, kmeans, method = "gap_stat")


# Compute hierarchical clustering and cut into 4 clusters
res <- hcut(dh2, k = 3, stand = TRUE)
# Visualize
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800"))















