# ---------------------------------------
# Merritt Burch
# Double Haploid Factor Analysis
# Looking at rate of pollenation/silking
# Summer 2014 Data (By Vivek)
# 12/15/2016
# ---------------------------------------

# Set working directory
setwd("~/SDSU Masters/Vivek Data/DH14_Images/Rate of silk and pol")

# Clear Directory
rm(list=ls())

# Load Packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'dplyr', 'tidyr', 'psych', 'rcompanion',
              'multcompView', 'doBy', 'RCurl', 'stats', 'ggbiplot',
              'factoextra')
pack.man(packages)


# Load datframe
df.dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH14_BorisEdit.csv")

# Rename columns
names(df.dh)[names(df.dh) == "newtreat2"] <- "old.lingen"
names(df.dh)[names(df.dh) == "newtreat3"] <- "lineage.generation"
names(df.dh)[names(df.dh) == "newtreat4"] <- "lineage"

# Remove extraneous rows
df.dh <- df.dh[, -1]
tbl_df(df.dh)

# Convert treatments to factors
df.dh[, 1:7] <- lapply(df.dh[, 1:7], factor)


# ---------------------
# Add columns for rate
# ---------------------

df.dh$rate.silk <- 1/df.dh$dayssilk
df.dh$rate.pol  <- 1/df.dh$dayspol

# Extract critical variables for PCA
df.pca <- df.dh[, c(8:18, 21:22)]


# ----------
# Histograms
# ----------

par(mfrow = c(2,2))

# Normal
hist(df.dh$dayspol)
hist(df.dh$dayssilk)
hist(df.dh$rate.pol)
hist(df.dh$rate.silk)


# ----------------------------
# Pricipal Components Analysis
# ----------------------------

# PCA
dh.pca.rate <- prcomp(df.pca, center = TRUE, scale = TRUE)

# Variance retained by each PC
eig.val.rate <- get_eigenvalue(dh.pca.rate)
(dh.eigen.rate <- eig.val.rate)

# Plotting

# Scree plot with variance
fviz_screeplot(dh.pca.rate, ncp=10)
ggsave(file = "DH14_rate_PCA_scree plot with variance.png", width = 10.5, height = 6)


# Individuals and variables, biplot is heated
fviz_pca_biplot(dh.pca.rate, label = "var", col.var="contrib", labelsize = 4.5, pointsize = 1)+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=30)
ggsave(file = "DH14_rate_PCA_heated biplot with ind.png", width = 10.5, height = 6)

# Individual map with heating
fviz_pca_ind(dh.pca.rate, col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.45) 
ggsave(file = "DH14_rate_PCA_heated ind.png", width = 10.5, height = 6)

# Tryouts, block
fviz_pca_ind(dh.pca.rate, geom = "point",
             habillage=df.dh$Block, addEllipses=TRUE,
             ellipse.level= 0.95)
ggsave(file = "DH14_rate_PCA_blocks with ellipses.png", width = 10.5, height = 6)

# Tryouts year
fviz_pca_ind(dh.pca.rate, geom = "point",
             habillage=df.dh$year, addEllipses=TRUE,
             ellipse.level= 0.95)
ggsave(file = "DH14_rate_PCA_years with ellipses.png", width = 10.5, height = 6)

# ----------------
# Factor analysis
# ----------------

# Retaining 5 components, varimax = orthagonal rotation
ls.fit <- principal(df.pca, nfactors = 5, rotate = "varimax") 

# Print results
ls.fit

# Sort and extract factor loadings
print(ls.fit$loadings, digits = 2, cutoff = .3, sort = T)


# See list
str(ls.fit)

# Assign 'scores' object to new data frame for easier manipulation
ls.fit.scores <- data.frame(ls.fit$scores)

# Merge treatments with scores
ls.fit.scores <- bind_cols(df.dh[, 1:7], ls.fit.scores)


#----------------------
# Visualization set up
#----------------------

# Make data long form
df.scores.gath <- gather(ls.fit.scores, component, value, RC1:RC5)

# Be 'lazy' and write a for loop to generate ALL DA GRAPHS
treatment <- colnames(ls.fit.scores[, 1:7])
comp      <- colnames(ls.fit.scores[, 8:12])


# -----------------------------
# Histograms of factor loadings
# -----------------------------

for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(value)) +
    geom_histogram(binwidth = 0.5, colour = "white") +
    scale_y_continuous(paste('Factor ', comp[i]))
  print(plot)
}

# ---------------------------
# Boxplots of factor loadings
# ---------------------------

for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(comp[i],value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]))
  print(plot)
}


# ------------------------------------
# Kruskal Wallis of Factors by lineage
# ------------------------------------

# Loop
# Column 6 is lineage
for (i in 8:12) {
  
  tmp <- kruskal.test(ls.fit.scores[[i]] ~ ls.fit.scores[[6]])
  print(tmp)
  
}

# ------------------------------------------------
# Plot non significant lineages by factors (1 & 4)
# ------------------------------------------------

# Plot lineage by FA 1
ggplot(ls.fit.scores, aes(lineage, RC1)) +
  geom_boxplot() +
  theme(text = element_text(size = 30))
ggsave(file="DH14_rate_lineage_FA1_KW.WT_boxplot.png", width = 10.5, height = 6)

# Facets for lineage by block FA1
ggplot(ls.fit.scores, aes(lineage, RC1)) +
  geom_boxplot() +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ Block)
ggsave(file="DH14_rate_lineage by block_FA1_boxplot.png", width = 10.5, height = 6)

# Facets for lineage by year FA1
ggplot(ls.fit.scores, aes(lineage, RC1)) +
  geom_boxplot() +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ year)
ggsave(file="DH14_rate_lineage by year_FA1_boxplot.png", width = 10.5, height = 6)

# --------- FA4 -------------------------------------

# Plot FA4 lineage by block
ggplot(ls.fit.scores, aes(lineage, RC4)) +
  geom_boxplot() +
  theme(text = element_text(size = 30))
ggsave(file="DH14_rate_lineage_FA4_KW.WT_boxplot.png", width = 10.5, height = 6)

# Facets for lineage by block FA4
ggplot(ls.fit.scores, aes(lineage, RC4)) +
  geom_boxplot() +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ Block)
ggsave(file="DH14_rate_lineage by block_FA4_boxplot.png", width = 10.5, height = 6)

# Facets for lineage by year FA4
ggplot(ls.fit.scores, aes(lineage, RC4)) +
  geom_boxplot() +
  theme(text = element_text(size = 20)) +
  facet_grid(. ~ year)
ggsave(file="DH14_rate_lineage by year_FA4_boxplot.png", width = 10.5, height = 6)


#------------------------------------
# Post Hoc test, lineage vs factor 2
#------------------------------------

# Using full code for this example, future post hoc tests will
#   be done using the function
{ 
  m11 <- pairwise.wilcox.test(ls.fit.scores$RC2,
                            ls.fit.scores$lineage,
                            p.adjust.method="none")
  m11 <- m11$p.value
  m11 <- fullPTable(m11)
  
  n11 <- multcompLetters(m11,
                       compare = "<",
                       threshold = 0.05,
                       Letters = letters,
                       reversed = FALSE)
  (o11 <- n11$Letters)
  
  (output11 <- matrix(unlist(o11), ncol = 1, byrow = T))
  
  p11 <- summaryBy(RC2 ~ lineage, 
                 data = ls.fit.scores,
                 FUN = function(x) { c( m = max(x)) } )
  
  p11$posth <- output11
  
  # Create spot for letters
  p11$letters <- p11$RC2.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(lineage, RC2)) +
    geom_boxplot() +
    theme(text = element_text(size=30)) +
    geom_text(data = p11, aes(x = lineage, y = letters, label = posth), size = 8)
  ggsave(file="DH14_rate_lineage_FA2_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each block
  ggplot(ls.fit.scores, aes(lineage, RC2)) +
    geom_boxplot() +
    theme(text = element_text(size = 20)) +
    geom_text(data = p11, aes(x = lineage, y = letters, label = posth), size = 6) +
    facet_grid(. ~ Block)
  ggsave(file="DH14_rate_lineage by block_FA2_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each year
  ggplot(ls.fit.scores, aes(lineage, RC2)) +
    geom_boxplot() +
    theme(text = element_text(size = 20)) +
    geom_text(data = p11, aes(x = lineage, y = letters, label = posth), size = 6, angle = 35) +
    facet_grid(. ~ year)
  ggsave(file="DH14_rate_lineage by year_FA2_KW.WT_boxplot.png", width = 10.5, height = 6)
  
}




#------------------------------------
# Post Hoc test, lineage vs factor 3
#------------------------------------

{ 
  m12 <- pairwise.wilcox.test(ls.fit.scores$RC3,
                              ls.fit.scores$lineage,
                              p.adjust.method="none")
  m12 <- m12$p.value
  m12 <- fullPTable(m12)
  
  n12 <- multcompLetters(m12,
                         compare = "<",
                         threshold = 0.05,
                         Letters = letters,
                         reversed = FALSE)
  (o12 <- n12$Letters)
  
  (output12 <- matrix(unlist(o12), ncol = 1, byrow = T))
  
  p12 <- summaryBy(RC3 ~ lineage, 
                   data = ls.fit.scores,
                   FUN = function(x) { c( m = max(x)) } )
  
  p12$posth <- output12
  
  # Create spot for letters
  p12$letters <- p12$RC3.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(lineage, RC3)) +
    geom_boxplot() +
    theme(text = element_text(size=30)) +
    geom_text(data = p12, aes(x = lineage, y = letters, label = posth), size = 8)
  ggsave(file="DH14_rate_lineage_FA3_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each block
  ggplot(ls.fit.scores, aes(lineage, RC3)) +
    geom_boxplot() +
    theme(text = element_text(size = 20)) +
    geom_text(data = p12, aes(x = lineage, y = letters, label = posth), size = 6) +
    facet_grid(. ~ Block)
  ggsave(file="DH14_rate_lineage by block_FA3_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each year
  ggplot(ls.fit.scores, aes(lineage, RC3)) +
    geom_boxplot() +
    theme(text = element_text(size = 20)) +
    geom_text(data = p12, aes(x = lineage, y = letters, label = posth), size = 6, angle = 35) +
    facet_grid(. ~ year)
  ggsave(file="DH14_rate_lineage by year_FA3_KW.WT_boxplot.png", width = 10.5, height = 6)
  
}


#------------------------------------
# Post Hoc test, lineage vs factor 5
#------------------------------------

{ 
  m13 <- pairwise.wilcox.test(ls.fit.scores$RC5,
                              ls.fit.scores$lineage,
                              p.adjust.method="none")
  m13 <- m13$p.value
  m13 <- fullPTable(m13)
  
  n13 <- multcompLetters(m13,
                         compare = "<",
                         threshold = 0.05,
                         Letters = letters,
                         reversed = FALSE)
  (o13 <- n13$Letters)
  
  (output13 <- matrix(unlist(o13), ncol = 1, byrow = T))
  
  p13 <- summaryBy(RC5 ~ lineage, 
                   data = ls.fit.scores,
                   FUN = function(x) { c( m = max(x)) } )
  
  p13$posth <- output13
  
  # Create spot for letters
  p13$letters <- p13$RC5.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(lineage, RC5)) +
    geom_boxplot() +
    theme(text = element_text(size = 30)) +
    geom_text(data = p13, aes(x = lineage, y = letters, label = posth), size = 8)
  ggsave(file="DH14_rate_lineage_FA5_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each block
  ggplot(ls.fit.scores, aes(lineage, RC5)) +
    geom_boxplot() +
    theme(text = element_text(size = 20)) +
    geom_text(data = p13, aes(x = lineage, y = letters, label = posth), size = 6) +
    facet_grid(. ~ Block)
  ggsave(file="DH14_rate_lineage by block_FA5_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each year
  ggplot(ls.fit.scores, aes(lineage, RC5)) +
    geom_boxplot() +
    theme(text = element_text(size = 20)) +
    geom_text(data = p13, aes(x = lineage, y = letters, label = posth), size = 6, angle = 35) +
    facet_grid(. ~ year)
  ggsave(file="DH14_rate_lineage by year_FA5_KW.WT_boxplot.png", width = 10.5, height = 6)
  
}


# -----------------------------------------------
# Kruskal Wallis of Factors by lineage-genrations
# -----------------------------------------------

for (i in 8:12) {
  
  tmp <- kruskal.test(ls.fit.scores[[i]] ~ ls.fit.scores[[5]])
  print(tmp)
}

#-----------------------------------------------
# Post Hoc test, lineage generation vs factor 1
#-----------------------------------------------

 { 
    m <- pairwise.wilcox.test(ls.fit.scores$RC1,
                          ls.fit.scores$lineage.generation,
                           p.adjust.method="none")
    m <- m$p.value
    m <- fullPTable(m)

    n <- multcompLetters(m,
                compare = "<",
                threshold = 0.05,
                Letters = letters,
                reversed = FALSE)
    (o <- n$Letters)

    (output <- matrix(unlist(o), ncol = 1, byrow = T))
    
    p <- summaryBy(RC1 ~ lineage.generation, 
                    data = ls.fit.scores,
                    FUN = function(x) { c( m = max(x)) } )
    
    p$posth <- output
    
    # Create spot for letters
    p$letters <- p$RC1.m + 0.5
    
    # Plot
    ggplot(ls.fit.scores, aes(lineage.generation, RC1)) +
      geom_boxplot() +
      theme(text = element_text(size=20), 
            axis.text.x = element_text(angle = 35, hjust = 1)) +
      geom_text(data = p, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle  = 10)
    ggsave(file="DH14_rate_lin.gen_FA1_KW.WT_boxplot.png", width = 10.5, height = 6)
    
    # Facets for each block
    ggplot(ls.fit.scores, aes(lineage.generation, RC1)) +
      geom_boxplot() +
      theme(text = element_text(size = 15), 
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(data = p, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
      facet_grid(. ~ Block)
    ggsave(file="DH14_rate_lin.gen by block_FA1_KW.WT_boxplot.png", width = 10.5, height = 6)
    
    # Facets for each year
    ggplot(ls.fit.scores, aes(lineage.generation, RC1)) +
      geom_boxplot() +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(data = p, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
      facet_grid(. ~ year)
    ggsave(file="DH14_rate_lin.gen by year_FA1_KW.WT_boxplot.png", width = 10.5, height = 6)
    
}

#-----------------------------------------------
# Post Hoc test, lineage generation vs factor 2
#-----------------------------------------------

{ 
  m2 <- pairwise.wilcox.test(ls.fit.scores$RC2,
                            ls.fit.scores$lineage.generation,
                            p.adjust.method="none")
  m2 <- m2$p.value
  m2 <- fullPTable(m2)
  
  n2 <- multcompLetters(m2,
                       compare = "<",
                       threshold = 0.05,
                       Letters = letters,
                       reversed = FALSE)
  (o2 <- n2$Letters)
  
  (output2 <- matrix(unlist(o2), ncol = 1, byrow = T))
  
  p2 <- summaryBy(RC2 ~ lineage.generation, 
                  data = ls.fit.scores,
                  FUN = function(x) { c( m = max(x)) } )
  
  p2$posth <- output2
  
  # Create spot for letters
  p2$letters <- p2$RC2.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(lineage.generation, RC2)) +
    geom_boxplot() +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 35, hjust = 1)) +
    geom_text(data = p2, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 10)
  ggsave(file="DH14_rate_lin.gen_FA2_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each block
  ggplot(ls.fit.scores, aes(lineage.generation, RC2)) +
    geom_boxplot() +
    theme(text = element_text(size = 15), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = p2, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
    facet_grid(. ~ Block)
  ggsave(file="DH14_rate_lin.gen by block_FA2_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each year
  ggplot(ls.fit.scores, aes(lineage.generation, RC2)) +
    geom_boxplot() +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = p2, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
    facet_grid(. ~ year)
  ggsave(file="DH14_rate_lin.gen by year_FA2_KW.WT_boxplot.png", width = 10.5, height = 6)
  
}


#-----------------------------------------------
# Post Hoc test, lineage generation vs factor 3
#-----------------------------------------------

{ 
  m3 <- pairwise.wilcox.test(ls.fit.scores$RC3,
                             ls.fit.scores$lineage.generation,
                             p.adjust.method="none")
  m3 <- m3$p.value
  m3 <- fullPTable(m3)
  
  n3 <- multcompLetters(m3,
                        compare = "<",
                        threshold = 0.05,
                        Letters = letters,
                        reversed = FALSE)
  (o3 <- n3$Letters)
  
  (output3 <- matrix(unlist(o3), ncol = 1, byrow = T))
  
  p3 <- summaryBy(RC3 ~ lineage.generation, 
                  data = ls.fit.scores,
                  FUN = function(x) { c( m = max(x)) } )
  
  p3$posth <- output3
  
  # Create spot for letters
  p3$letters <- p3$RC3.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(lineage.generation, RC3)) +
    geom_boxplot() +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 35, hjust = 1)) +
    geom_text(data = p3, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 10)
  ggsave(file="DH14_rate_FA3_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each block
  ggplot(ls.fit.scores, aes(lineage.generation, RC3)) +
    geom_boxplot() +
    theme(text = element_text(size = 15), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = p3, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
    facet_grid(. ~ Block)
  ggsave(file="DH14_rate_lin.gen by block_FA3_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each year
  ggplot(ls.fit.scores, aes(lineage.generation, RC3)) +
    geom_boxplot() +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = p3, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
    facet_grid(. ~ year)
  ggsave(file="DH14_rate_lin.gen by year_FA3_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  
}


#-----------------------------------------------
# Post Hoc test, lineage generation vs factor 4
#-----------------------------------------------

{ 
  m4 <- pairwise.wilcox.test(ls.fit.scores$RC4,
                             ls.fit.scores$lineage.generation,
                             p.adjust.method="none")
  m4 <- m4$p.value
  m4 <- fullPTable(m4)
  
  n4 <- multcompLetters(m4,
                        compare = "<",
                        threshold = 0.05,
                        Letters = letters,
                        reversed = FALSE)
  (o4 <- n4$Letters)
  
  (output4 <- matrix(unlist(o4), ncol = 1, byrow = T))
  
  p4 <- summaryBy(RC4 ~ lineage.generation, 
                  data = ls.fit.scores,
                  FUN = function(x) { c( m = max(x)) } )
  
  p4$posth <- output4
  
  # Create spot for letters
  p4$letters <- p4$RC4.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(lineage.generation, RC4)) +
    geom_boxplot() +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 35, hjust = 1)) +
    geom_text(data = p4, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 25)
  ggsave(file="DH14_rate_FA4_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each block
  ggplot(ls.fit.scores, aes(lineage.generation, RC4)) +
    geom_boxplot() +
    theme(text = element_text(size = 15), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = p4, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
    facet_grid(. ~ Block)
  ggsave(file="DH14_rate_lin.gen by block_FA4_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each year
  ggplot(ls.fit.scores, aes(lineage.generation, RC4)) +
    geom_boxplot() +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = p4, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
    facet_grid(. ~ year)
  ggsave(file="DH14_rate_lin.gen by year_FA4_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  
}


#-----------------------------------------------
# Post Hoc test, lineage generation vs factor 5
#-----------------------------------------------

{ 
  m5 <- pairwise.wilcox.test(ls.fit.scores$RC5,
                             ls.fit.scores$lineage.generation,
                             p.adjust.method="none")
  m5 <- m5$p.value
  m5 <- fullPTable(m5)
  
  n5 <- multcompLetters(m5,
                        compare = "<",
                        threshold = 0.05,
                        Letters = letters,
                        reversed = FALSE)
  (o5 <- n5$Letters)
  
  (output5 <- matrix(unlist(o5), ncol = 1, byrow = T))
  
  p5 <- summaryBy(RC5 ~ lineage.generation, 
                  data = ls.fit.scores,
                  FUN = function(x) { c( m = max(x)) } )
  
  p5$posth <- output5
  
  # Create spot for letters
  p5$letters <- p5$RC5.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(lineage.generation, RC5)) +
    geom_boxplot() +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 35, hjust = 1)) +
    geom_text(data = p5, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 10)
  ggsave(file="DH14_rate_FA5_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each block
  ggplot(ls.fit.scores, aes(lineage.generation, RC5)) +
    geom_boxplot() +
    theme(text = element_text(size = 15), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = p5, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
    facet_grid(. ~ Block)
  ggsave(file="DH14_rate_lin.gen by block_FA5_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  # Facets for each year
  ggplot(ls.fit.scores, aes(lineage.generation, RC5)) +
    geom_boxplot() +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(data = p5, aes(x = lineage.generation, y = letters, label = posth), size = 6, angle = 80) +
    facet_grid(. ~ year)
  ggsave(file="DH14_rate_lin.gen by year_FA5_KW.WT_boxplot.png", width = 10.5, height = 6)
  
  
}


# -------------------------------------------
# Loop for Factor loadings and year produced
# -------------------------------------------

# ---------------------------------
# Kruskal Wallis of Factors by year
# ---------------------------------

for (i in 8:12) {
  
  tmp <- kruskal.test(ls.fit.scores[[i]] ~ ls.fit.scores[[7]])
  print(tmp)
  
}

# ------------
# Year graphs
# ------------

for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(year, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 25, hjust = 1)) 
  print(plot)
}

#---------------------------
# Factor 2 post hoc on year
# --------------------------

{ 
  m6 <- pairwise.wilcox.test(ls.fit.scores$RC2,
                             ls.fit.scores$year,
                             p.adjust.method="none")
  m6 <- m6$p.value
  m6 <- fullPTable(m6)
  
  n6 <- multcompLetters(m6,
                        compare = "<",
                        threshold = 0.05,
                        Letters = letters,
                        reversed = FALSE)
  (o6 <- n6$Letters)
  
  (output6 <- matrix(unlist(o6), ncol = 1, byrow = T))
  
  p6 <- summaryBy(RC2 ~ year, 
                  data = ls.fit.scores,
                  FUN = function(x) { c( m = max(x)) } )
  
  p6$posth <- output6
  
  # Create spot for letters
  p6$letters <- p6$RC2.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(year, RC2)) +
    geom_boxplot() +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 35, hjust = 1)) +
    geom_text(data = p6, aes(x = year, y = letters, label = posth), size = 8)
  ggsave(file="DH14_rate_year_FA2_KW.WT_boxplot.png", width = 10.5, height = 6)
  
}


#---------------------------
# Factor 4 post hoc on year
# --------------------------

{ 
  m7 <- pairwise.wilcox.test(ls.fit.scores$RC4,
                             ls.fit.scores$year,
                             p.adjust.method="none")
  m7 <- m7$p.value
  m7 <- fullPTable(m7)
  
  n7 <- multcompLetters(m7,
                        compare = "<",
                        threshold = 0.05,
                        Letters = letters,
                        reversed = FALSE)
  (o7 <- n7$Letters)
  
  (output7 <- matrix(unlist(o7), ncol = 1, byrow = T))
  
  p7 <- summaryBy(RC4 ~ year, 
                  data = ls.fit.scores,
                  FUN = function(x) { c( m = max(x)) } )
  
  p7$posth <- output7
  
  # Create spot for letters
  p7$letters <- p7$RC4.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(year, RC4)) +
    geom_boxplot() +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle = 35, hjust = 1)) +
    geom_text(data = p7, aes(x = year, y = letters, label = posth), size = 8)
  ggsave(file="DH14_rate_year_FA4_KW.WT_boxplot.png", width = 10.5, height = 6)
  
}

#







# -------------------------------------
# Kruskal-Wallis on blocks vs. factors
# -------------------------------------

for (i in 8:12) {
  
  tmp <- kruskal.test(ls.fit.scores[[i]] ~ ls.fit.scores[[2]])
  print(tmp)
  
}

# Block
for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(Block, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(text = element_text(size = 30))
  print(plot)
}

# --------------------------
# Post hoc for block in FA1
# --------------------------

{ 
  m8 <- pairwise.wilcox.test(ls.fit.scores$RC1,
                             ls.fit.scores$Block,
                             p.adjust.method="none")
  m8 <- m8$p.value
  m8 <- fullPTable(m8)
  
  n8 <- multcompLetters(m8,
                        compare = "<",
                        threshold = 0.05,
                        Letters = letters,
                        reversed = FALSE)
  (o8 <- n8$Letters)
  
  (output8 <- matrix(unlist(o8), ncol = 1, byrow = T))
  
  p8 <- summaryBy(RC1 ~ Block, 
                  data = ls.fit.scores,
                  FUN = function(x) { c( m = max(x)) } )
  
  p8$posth <- output8
  
  # Create spot for letters
  p8$letters <- p8$RC1.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(Block, RC1)) +
    geom_boxplot() +
    theme(text = element_text(size=30)) +
    geom_text(data = p8, aes(x = Block, y = letters, label = posth), size = 8)
  ggsave(file="DH14_rate_block_FA1_KW.WT_boxplot.png", width = 10.5, height = 6)
  
}


# --------------------------
# Post hoc for block in FA2
# --------------------------

{ 
  m9 <- pairwise.wilcox.test(ls.fit.scores$RC2,
                             ls.fit.scores$Block,
                             p.adjust.method="none")
  m9 <- m9$p.value
  m9 <- fullPTable(m9)
  
  n9 <- multcompLetters(m9,
                        compare = "<",
                        threshold = 0.05,
                        Letters = letters,
                        reversed = FALSE)
  (o9 <- n9$Letters)
  
  (output9 <- matrix(unlist(o9), ncol = 1, byrow = T))
  
  p9 <- summaryBy(RC2 ~ Block, 
                  data = ls.fit.scores,
                  FUN = function(x) { c( m = max(x)) } )
  
  p9$posth <- output9
  
  # Create spot for letters
  p9$letters <- p9$RC2.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(Block, RC2)) +
    geom_boxplot() +
    theme(text = element_text(size=30)) +
    geom_text(data = p9, aes(x = Block, y = letters, label = posth), size = 8)
  ggsave(file="DH14_rate_block_FA2_KW.WT_boxplot.png", width = 10.5, height = 6)
  
}


# --------------------------
# Post hoc for block in FA4
# --------------------------

{ 
  m10 <- pairwise.wilcox.test(ls.fit.scores$RC4,
                             ls.fit.scores$Block,
                             p.adjust.method="none")
  m10 <- m10$p.value
  m10 <- fullPTable(m10)
  
  n10 <- multcompLetters(m10,
                        compare = "<",
                        threshold = 0.05,
                        Letters = letters,
                        reversed = FALSE)
  (o10 <- n10$Letters)
  
  (output10 <- matrix(unlist(o10), ncol = 1, byrow = T))
  
  p10 <- summaryBy(RC4 ~ Block, 
                  data = ls.fit.scores,
                  FUN = function(x) { c( m = max(x)) } )
  
  p10$posth <- output10
  
  # Create spot for letters
  p10$letters <- p10$RC4.m + 0.5
  
  # Plot
  ggplot(ls.fit.scores, aes(Block, RC4)) +
    geom_boxplot() +
    theme(text = element_text(size=30)) +
    geom_text(data = p10, aes(x = Block, y = letters, label = posth), size = 8)
  ggsave(file="DH14_rate_block_FA4_KW.WT_boxplot.png", width = 10.5, height = 6)
  
}

