# ------------------------------------
# Merritt Burch does things with data!
# Double Haploid Factor Analysis
# Summer 2014 Data
# 10/03/2016 collected by Vivek
# ------------------------------------

# Clear Directory
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/DH_15")

# Load Packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'dplyr', 'tidyr', 'psych')
pack.man(packages)

# Import Data
df.15 <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH15_raw data_edited.csv")

# Remove column
df.15 <- df.15[-22, -1]

# Convert treatments to factors
df.15[, 1:6] <- lapply(df.15[, 1:6], factor)

# Extract critical variables for PCA
df.15.pca <- df.15[, c(7:11, 14:15)]

# -----------------------------
# Pricipal Components Analysis
# -----------------------------

# Retaining 5 components, varimax = orthagonal rotation
ls.fit <- principal(df.pca, nfactors = 4, rotate="varimax") 

# Print results
ls.fit

# Sort and extract factor loadings
print(ls.fit$loadings, digits=2, cutoff=.3, sort=T)

# See list
str(ls.fit)

# Assign 'scores' object to new data frame for easier manipulation
ls.fit.scores <- data.frame(ls.fit$scores)

# Merge treatments with scores
ls.fit.scores <- bind_cols(df.15[, 1:6], ls.fit.scores)

# Reorder data
temp <- c(1:8, 10, 9)
ls.fit.scores <- ls.fit.scores[, temp]

#--------------
# Visualization
#--------------

# Make data long form
df.scores.gath <- gather(ls.fit.scores, component, value, RC1:RC4)

# Be 'lazy' and write a for loop to generate ALL DA GRAPHS
treatment <- colnames(ls.fit.scores[, 1:6])
comp      <- colnames(ls.fit.scores[, 7:10])


# For lineage
for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(lin, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          text = element_text(size = 30)) +
    xlab("Lineage")
  ggtitle(paste('Lineage vs.', comp[i]))
  print(plot)
}


# Lineage.generation
for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(lineage.generation, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          text = element_text(size = 30)) +
    xlab("Lineage")
  ggtitle(paste('Lineage-Generation vs.', comp[i]))
  print(plot)
}

# block
for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(Block, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          text = element_text(size = 30))
  print(plot)
}

# -----------------------
# DO kruskal wallis test
# Lineage
# Lineage-Generation
# Block
# -----------------------

# Lineage
for (i in 7:10) {
  
  tmp <- kruskal.test(ls.fit.scores[[i]] ~ ls.fit.scores[[6]])
  print(tmp)
}

# Lineage-Generation
for (i in 7:10) {
  
  tmp <- kruskal.test(ls.fit.scores[[i]] ~ ls.fit.scores[[4]])
  print(tmp)
}

# Block
for (i in 7:10) {
  
  tmp <- kruskal.test(ls.fit.scores[[i]] ~ ls.fit.scores[[2]])
  print(tmp)
}



# ----------------------------------------
# Function for post hoc tests and graphing
# ----------------------------------------

# Lineage
ggLetters <- function(data, resp, group, file, width, height) {
  
  tmp.pwt <- stats::pairwise.wilcox.test(x        = data[[resp]],
                                         g        = data[[group]],
                                         p.adjust = 'none')
  
  tmp.pwt <- tmp.pwt$p.value
  tmp.pwt <- rcompanion::fullPTable(tmp.pwt)
  
  tmp.mcl <- multcompView::multcompLetters(tmp.pwt,
                                           compare   = '<',
                                           threshold = 0.05,
                                           Letters   = letters,
                                           reversed  = FALSE)
  
  tmp.let <- tmp.mcl$Letters
  tmp.mat <- matrix(unlist(tmp.let), ncol = 1, byrow = TRUE)
  tmp.fun <- function(x) {c(m = max(x))}
  tmp.frm <- paste(resp, '~', group)
  tmp.frm <- as.formula(tmp.frm)
  
  tmp.sum <- doBy::summaryBy(tmp.frm,
                             data = data,
                             FUN  = tmp.fun)
  
  tmp.sum$posth   <- tmp.mat
  tmp.sum$letters <- tmp.sum[2] + 0.5
  
  tmp.plot <- ggplot(ls.fit.scores, 
                     aes_string(x = group, y = resp)) +
    geom_boxplot() +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 0, hjust = 1)) +
    geom_text(data = tmp.sum, 
              aes(x = tmp.sum[1], y = letters, label = posth), 
              size = 8) +
    xlab("Lineage")
    ggsave(file = file, width = width, height = height)
  
  print(tmp.plot)
  
}


# Run function and produce graphs
ggLetters(data = ls.fit.scores, 
          resp = 'RC4', 
          group = 'lin', 
          file = 'DH15_lin_FA4_KW.WT_boxplot.png', 
          width = 10, 
          height = 6)


# Lineage-generation
ggLetters <- function(data, resp, group, file, width, height) {
  
  tmp.pwt <- stats::pairwise.wilcox.test(x        = data[[resp]],
                                         g        = data[[group]],
                                         p.adjust = 'none')
  
  tmp.pwt <- tmp.pwt$p.value
  tmp.pwt <- rcompanion::fullPTable(tmp.pwt)
  
  tmp.mcl <- multcompView::multcompLetters(tmp.pwt,
                                           compare   = '<',
                                           threshold = 0.05,
                                           Letters   = letters,
                                           reversed  = FALSE)
  
  tmp.let <- tmp.mcl$Letters
  tmp.mat <- matrix(unlist(tmp.let), ncol = 1, byrow = TRUE)
  tmp.fun <- function(x) {c(m = max(x))}
  tmp.frm <- paste(resp, '~', group)
  tmp.frm <- as.formula(tmp.frm)
  
  tmp.sum <- doBy::summaryBy(tmp.frm,
                             data = data,
                             FUN  = tmp.fun)
  
  tmp.sum$posth   <- tmp.mat
  tmp.sum$letters <- tmp.sum[2] + 0.5
  
  tmp.plot <- ggplot(ls.fit.scores, 
                     aes_string(x = group, y = resp)) +
    geom_boxplot() +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    geom_text(data = tmp.sum, 
              aes(x = tmp.sum[1], y = letters, label = posth), 
              size = 8) +
    xlab("Lineage - Generation")
  ggsave(file = file, width = width, height = height)
  
  print(tmp.plot)
  
}

ggLetters(data = ls.fit.scores, 
          resp = 'RC1', 
          group = 'lineage.generation', 
          file = 'DH15_lingen_FA1_KW.WT_boxplot.png', 
          width = 10, 
          height = 6)

ggLetters(data = ls.fit.scores, 
          resp = 'RC2', 
          group = 'lineage.generation', 
          file = 'DH15_lingen_FA2_KW.WT_boxplot.png', 
          width = 10, 
          height = 6)

ggLetters(data = ls.fit.scores, 
          resp = 'RC3', 
          group = 'lineage.generation', 
          file = 'DH15_lingen_FA3_KW.WT_boxplot.png', 
          width = 10, 
          height = 6)














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


