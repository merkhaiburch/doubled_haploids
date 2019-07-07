# ---------------------------------------------------------------------------------
# Merritt Burch
# 11/29/2016
# Correlation matrix for DH data
# DH 15 Data
# Code source:
# http://stackoverflow.com/questions/19012529/correlation-corrplot-configuration
# ------------------------------------------------------------------------------


# Refresh environment
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/DH15")

# Import data
df.15 <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH15_raw data_edited.csv")

# Subset data
df.15 <- df.15[, -1]
df2 <- (df.15[, 7:15])

# Locad packages 
library(corrplot)
library(corrgram)


# ------------------------
# Plot Correlation matrix
# ------------------------

panel.shadeNtext <- function (x, y, corr = NULL, col.regions, ...) 
{
  corr <- cor(x, y, use = "pair")
  results <- cor.test(x, y, alternative = "two.sided")
  est <- results$p.value
  stars <- ifelse(est < 5e-4, "***", 
                  ifelse(est < 5e-3, "**", 
                         ifelse(est < 5e-2, "*", "")))
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, 
                                               length = ncol + 1), include.lowest = TRUE))
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], 
       border = NA)
  box(col = "lightgray")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- formatC(corr, digits = 2, format = "f")
  cex.cor <- 0.5/strwidth("-X.xx")
  fonts <- ifelse(stars != "", 2,1)
  # option 1: stars:
  text(0.5, 0.4, paste0(r,"\n", stars), cex = cex.cor)
  
  # option 2: bolding:
  # text(0.5, 0.5, r, cex = cex.cor, font=fonts)
}


# Plot matrix
pdf(file = "DH2017 Cor Matrix.pdf", width = 8, height = 5)
corrgram(pca.2017, type="data", lower.panel=panel.shadeNtext, 
         upper.panel=NULL)
dev.off()

pdf(file = "DH2014 Cor Matrix.pdf", width = 5, height = 6)

# Plot matrix w/upper panel of data scatterplots
corrgram(pca.2014, type="data", lower.panel=panel.shadeNtext, 
         upper.panel=panel.pts )


# Plot matrix w/upper panel of data eclipses w/confidence circles & distribution
corrgram(df15.pca, type="data", lower.panel=panel.shadeNtext, 
         upper.panel=panel.ellipse)






