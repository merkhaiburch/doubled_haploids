# ------------------------------------
# Merritt Burch does things with data!
# Double Haploid Factor Analysis
# Summer 2014 Data
# 10/03/2016collected by Vivek
# ------------------------------------

# Clear Directory
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/Data Sets")

# Load Packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'dplyr', 'tidyr', 'psych')
pack.man(packages)

# Load datframe
df.dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH raw data -2014_edited_V2.csv")

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


# -----------------------------
# Pricipal Components Analysis
# -----------------------------

# Retaining 5 components, varimax = orthagonal rotation
ls.fit <- principal(df.pca, nfactors=5, rotate="varimax") 

# Print results
ls.fit

# Sort and extract factor loadings
print(ls.fit$loadings, digits=2, cutoff=.3, sort=T)

# See list
str(ls.fit)

# Assign 'scores' object to new data frame for easier manipulation
ls.fit.scores <- data.frame(ls.fit$scores)

# Merge treatments with scores
ls.fit.scores <- bind_cols(df.dh[, 1:7], ls.fit.scores)


#--------------
# Visualization
#--------------

# Make data long form
df.scores.gath <- gather(ls.fit.scores, component, value, RC1:RC5)

# Show distribution of each block for factor (RC) 1
ggplot(ls.fit.scores, aes(factor(newtreat2), RC1)) +
  geom_boxplot()

# Be 'lazy' and write a for loop to generate ALL DA GRAPHS
treatment <- colnames(ls.fit.scores[, 1:7])
comp      <- colnames(ls.fit.scores[, 8:12])


# For lineage-genrations
for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(newtreat, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggtitle(paste('Lineage.Generation vs.', comp[i]))
  print(plot)
}


# Year
for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(year, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggtitle(paste('Year vs.', comp[i]))
  print(plot)
}

# Plant number
for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(plantno, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggtitle(paste('Plant Number vs.', comp[i]))
  print(plot)
}

# Family
for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(x = Family, y = value)) +
    geom_boxplot(aes(fill = newtreat)) +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggtitle(paste('Family vs.', comp[i]))
  print(plot)
}


# Block
for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(Block, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggtitle(paste('Block vs.', comp[i]))
  print(plot)
}


# ------------------------------------------
# For fun: try with no rotation and compare
# ------------------------------------------

# Unrotated factors
ls.fit.unrotated <- principal(df.pca, nfactors=5, rotate="none") 

# Print results
ls.fit.unrotated

# Sort and extract factor loadings
print(ls.fit.unrotated$loadings, digits=2, cutoff=.3, sort=T)


# Compare

# Rotated
ls.fit <- principal(df.pca, nfactors=5, rotate="varimax") 

# Print results
ls.fit

# Sort and extract factor loadings
print(ls.fit$loadings, digits=2, cutoff=.3, sort=T)


# -------------------------
# For traditional FA plots
# -------------------------

# Unrotated - Factor 1 by Factor 2
load <- ls.fit.unrotated$loadings[,1:2]

# Set up plot
plot(load,type="n")

# Add variable names 
text(load,labels=names(df.pca),cex=.7)


# Rotated Compare

# Unrotated - Factor 1 by Factor 2
load <- ls.fit$loadings[,1:2]

# Set up plot
plot(load,type="n")

# Add variable names 
text(load,labels=names(df.pca),cex=.7)



