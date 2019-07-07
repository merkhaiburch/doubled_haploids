# -----------------------------------------
# Merritt Burch
# 2/18/2017
# DH14 Analysis with 5 factors
# Redo of 14 data with less variables
# FA with kruskal test and appro. post hoc
# -----------------------------------------

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
df.pca <- df.dh[, c(8, 10:13, 21:22)]

# ----------------
# Factor analysis
# ----------------

# Retaining 5 components, varimax = orthagonal rotation
ls.fit <- principal(df.pca, nfactors = 5, rotate = "varimax") 

# Sort and extract factor loadings
(r14 <- ls.fit$loadings[1:7 , ])