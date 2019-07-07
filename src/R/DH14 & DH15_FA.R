# ------------------------------------------
# Merritt Burch
# Double Haploid Factor Analysis
# Summer 2014  & 2015 data
# 2/1/17
# Redo of '14 analysis to only include
#   '15 traits, see if factor scores change
# ------------------------------------------

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


# Load dataframe
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
df.pca <- df.dh[, c(8, 10:13, 21:22)]


# ----------------------------
# Pricipal Components Analysis
# ----------------------------

# PCA
dh.pca.rate <- prcomp(df.pca, center = TRUE, scale = TRUE)

# Variance retained by each PC
eig.val.rate <- get_eigenvalue(dh.pca.rate)
(dh.eigen.rate <- eig.val.rate)
write.table(dh.eigen.rate, "dh.eigen.rate.txt", sep="\t") 

# Individuals and variables, biplot is heated
fviz_pca_biplot(dh.pca.rate, label = "var", col.var="contrib", labelsize = 4.5, pointsize = 1)+
  scale_color_gradient2(low="white", mid="blue",
                        high="red", midpoint=20)
# ggsave(file = "DH14_rate_PCA_heated biplot with ind.png", width = 10.5, height = 6)

# Individual map with heating
fviz_pca_ind(dh.pca.rate, col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.45) 
# ggsave(file = "DH14_rate_PCA_heated ind.png", width = 10.5, height = 6)


# ----------------
# Factor analysis
# ----------------

# Retaining 5 components, varimax = orthagonal rotation
ls.fit <- principal(df.pca, nfactors = 5, rotate = "varimax") 

# Sort and extract factor loadings
(r14 <- ls.fit$loadings[1:7 , ])
write.table(r14, "r14.txt", sep="\t") 

# -------------------------
# Plot rotated FA results
# -------------------------
load = ls.fit$loadings[,1:2]
plot(load, type="n") # set up plot 
text(load,labels=names(df.pca),cex=1.3) # add variable names

# --------------------------------
# Fancy barplot of factor loadings
# ---------------------------------

# Create table for factor loadings
r14.scores <- data.frame(r14)
tmp <- c(1,2,4,3,5)
r14.scores <- r14.scores[, tmp]
r14.scores <- setNames(cbind(rownames(r14.scores), r14.scores, row.names = NULL), 
         c("Trait", "RC1", "RC2", "RC3", "RC4", 'RC5'))
df2 <- gather(r14.scores, component, value, RC1:RC5)

#For each test, plot the loading as length and fill color of a bar
# note that the length will be the absolute value of the loading but the 
# fill color will be the signed value, more on this below
ggplot(df2, aes(Trait, abs(value), fill=value)) + 
  facet_wrap(~ component, nrow=1) + #place the factors in separate facets
  geom_bar(stat="identity") + #make the bars
  coord_flip() + #flip the axes so the test names can be horizontal  
  #define the fill color gradient: blue=positive, red=negative
  scale_fill_gradient2(name = "value", 
                       high = "blue", mid = "white", low = "red", 
                       midpoint=0) +
  ylab("Loading Strength") + #improve y-axis label
  theme_bw() + 
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file = "DH14_rotated data.png", width = 10, height = 6)



# --------------------------
# Do 2014 unrotated analysis
# --------------------------

# Retaining 5 components
ur14 <- principal(df.pca, nfactors = 5, rotate = "none") 

# -------------------------
# Plot rotated FA results
# -------------------------
load = ur14$loadings[,1:2]
plot(load, type="n") # set up plot 
text(load,labels=names(df.pca),cex=1.3) # add variable names


# Sort and extract factor loadings
(ur14 <- ur14$loadings[1:7 , ])

# Create table for factor loadings
ur14.scores <- data.frame(ur14)
tmp <- c(1,2,4,3,5)
ur14.scores <- ur14.scores[, tmp]
ur14.scores <- setNames(cbind(rownames(ur14.scores), ur14.scores, row.names = NULL), 
                       c("Trait", "RC1", "RC2", "RC3", "RC4", 'RC5'))
# Data in long form
ur14.scores <- gather(ur14.scores, component, value, RC1:RC5)

# Plot 
ggplot(ur14.scores, aes(Trait, abs(value), fill=value)) + 
  facet_wrap(~ component, nrow=1) + #place the factors in separate facets
  geom_bar(stat="identity") + #make the bars
  coord_flip() + #flip the axes so the test names can be horizontal  
  #define the fill color gradient: blue=positive, red=negative
  scale_fill_gradient2(name = "value", 
                       high = "blue", mid = "white", low = "red", 
                       midpoint=0) +
  ylab("Loading Strength") + #improve y-axis label
  theme_bw() + 
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file = "DH14_UNrotated data.png", width = 10, height = 6)






# --------------------
# Rotated 2015 portion
# --------------------

# Import Data
df.15 <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH15_raw data_edited.csv")

# Remove column
df.15 <- df.15[-22, -1]

# Convert treatments to factors
df.15[, 1:6] <- lapply(df.15[, 1:6], factor)

# Extract critical variables for PCA
df15.pca <- df.15[, c(7:11, 14:15)]

# PCA for eigenvalues
dh15.pca <- prcomp(df15.pca, center = TRUE, scale = TRUE)
eig.val.rate <- get_eigenvalue(dh15.pca)
(dh.eigen.rate <- eig.val.rate)
write.table(dh.eigen.rate, "temp.txt", sep="\t") 

# Retaining 5 components
r15 <- principal(df15.pca, nfactors = 5, rotate = "varimax") 

# -------------------------
# Plot rotated FA results
# -------------------------
load = r15$loadings[,1:2]
plot(load, type="n") # set up plot 
text(load,labels=names(df.pca),cex=.7) # add variable names


# Sort and extract factor loadings
(r15 <- r15$loadings[1:7 , 1:5])

# Create table for factor loadings
r15.scores <- data.frame(r15)
tmp <- c(1,5,3,4,2)
r15.scores <- r15.scores[, tmp]
r15.scores <- setNames(cbind(rownames(r15.scores), r15.scores, row.names = NULL), 
                        c("Trait", "RC1", "RC2", "RC3", "RC4", "RC5"))
# Data in long form
r15.scores <- gather(r15.scores, component, value, RC1:RC4)

# Plot 
ggplot(r15.scores, aes(Trait, abs(value), fill=value)) + 
  facet_wrap(~ component, nrow=1) + #place the factors in separate facets
  geom_bar(stat="identity") + #make the bars
  coord_flip() + #flip the axes so the test names can be horizontal  
  #define the fill color gradient: blue=positive, red=negative
  scale_fill_gradient2(name = "value", 
                       high = "blue", mid = "white", low = "red", 
                       midpoint=0) +
  ylab("Loading Strength") + #improve y-axis label
  theme_bw() + 
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file = "DH15_rotated data.png", width = 10, height = 6)


# -----------
#  Unrotated
# -----------

# Retaining 5 components
ur15 <- principal(df15.pca, nfactors = 5, rotate = "none") 

# -------------------------
# Plot rotated FA results
# -------------------------
load = ur14$loadings[,1:2]
plot(load, type="n") # set up plot 
text(load,labels=names(df.pca),cex=.7) # add variable names


# Sort and extract factor loadings
ur15 <- ur15$loadings[1:7 , ]

# Create table for factor loadings
ur15.scores <- data.frame(ur15)
tmp <- c(1,5,3,4,2)
ur15.scores <- ur15.scores[, tmp]
ur15.scores <- setNames(cbind(rownames(ur15.scores), ur15.scores, row.names = NULL), 
                       c("Trait", "RC1", "RC2", "RC3", "RC4", 'RC5'))
# Data in long form
ur15.scores <- gather(ur15.scores, component, value, RC1:RC5)

# Plot 
ggplot(ur15.scores, aes(Trait, abs(value), fill=value)) + 
  facet_wrap(~ component, nrow=1) + #place the factors in separate facets
  geom_bar(stat="identity") + #make the bars
  coord_flip() + #flip the axes so the test names can be horizontal  
  #define the fill color gradient: blue=positive, red=negative
  scale_fill_gradient2(name = "value", 
                       high = "blue", mid = "white", low = "red", 
                       midpoint=0) +
  ylab("Loading Strength") + #improve y-axis label
  theme_bw() + 
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file = "DH15_UNrotated data.png", width = 10, height = 6)

# Plot rotated FA results
load = ls.fit3$loadings[,1:2]
plot(load, type="n") # set up plot 
text(load,labels=names(df15.pca),cex=.7) # add variable names





