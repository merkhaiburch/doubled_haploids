#-----------------------------------------------------#
# Title:  Consultation - Factor Analysis (Merritt)    #
# Author: Brandon Monier (brandon.monier@sdstate.edu) #
# Date:   11.07.16                                    #
#-----------------------------------------------------#

#---------
# Preamble
#---------

# Load packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'dplyr', 'tidyr', 'psych')
pack.man(packages)

# Load datframe
df.corn <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH raw data -2014_edited_V2.csv")

# Remove extraneous rows
df.corn <- df.corn[, c(-1, -2)]
tbl_df(df.corn)

# Reorder this dataframe (it's a mess!)
tmp <- c(1, 2, 3, 5, 4, 20, 19, 6:18)
df.corn <- df.corn[, tmp]

# Convert treatments to factors
df.corn[, 1:7] <- lapply(df.corn[, 1:7], factor)

# Extract critical variables for PCA
df.pca <- df.corn[, 8:20]


#------------------------------
# Principal components analysis
#------------------------------

# Create PCA list
ls.fit <- principal(df.pca, nfactors = 5, rotate = 'varimax')

# Inspect PCA list
str(ls.fit)

# Assign 'scores' object to new data frame for easier manipulation
df.fit.scores <- data.frame(ls.fit$scores)

# Merge treatments with scores
df.fit.scores <- bind_cols(df.corn[, 1:7], df.fit.scores)


#----------------------
# Visualization example
#----------------------

# Make data long form
df.scores.gath <- gather(df.fit.scores, component, value, RC1:RC5)

# Show distribution of each block for factor (RC) 1
ggplot(df.fit.scores, aes(factor(newtreat2), RC1)) +
  geom_boxplot()

# Be 'lazy' and write a for loop to generate ALL DA GRAPHS
treatment <- colnames(df.fit.scores[, 1:7])
comp      <- colnames(df.fit.scores[, 8:12])

for (i in seq_along(comp)) {
  plot <-
    ggplot(subset(df.scores.gath, df.scores.gath$component == comp[i]),
           aes(newtreat2, value)) +
    geom_boxplot() +
    scale_y_continuous(paste('Factor ', comp[i]), limits = c(-4, 4)) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggtitle(paste('Lineage.Generation vs.', comp[i]))
  print(plot)
}
