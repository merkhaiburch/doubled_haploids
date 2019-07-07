# ---------------------------------------------------------------------------------
# Merritt Burch does things with data!
# 10/19/2016
# Correlation matrix for DH data
# DH 14 Data
# Code source:
# http://stackoverflow.com/questions/19012529/correlation-corrplot-configuration
# ---------------------------------------------------------------------------------


# Refresh environment
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/Data Sets")

# Import data
dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH raw data -2014_edited.csv")

# Subset data
dh2 <- (dh[, 7:20])

# Check data
str(dh)
summary(dh)

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
  cex.cor <- .8/strwidth("-X.xx")
  fonts <- ifelse(stars != "", 2,1)
  # option 1: stars:
  text(0.5, 0.4, paste0(r,"\n", stars), cex = cex.cor)
  # option 2: bolding:
  # text(0.5, 0.5, r, cex = cex.cor, font=fonts)
}


# Plot matrix
corrgram(dh2, type="data", lower.panel=panel.shadeNtext, 
         upper.panel=NULL)

# Plot matrix w/upper panel of data scatterplots
corrgram(dh2, type="data", lower.panel=panel.shadeNtext, 
         upper.panel=panel.pts )


# Plot matrix w/upper panel of data eclipses w/confidence circles & distribution
corrgram(dh2, type="data", lower.panel=panel.shadeNtext, 
         upper.panel=panel.ellipse)






