# ----------------------------------------------------
# Merritt Does things with data!
# Double Haploid Linear mixed models & ANOVA Analysis
# Summer 2014 Data
# 8/30/2016
# ----------------------------------------------------

# Clear Directory
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/Data Sets")

# Import Data
df.dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH raw data -2014_edited_V2.csv")

# Remove extraneous rows
df.dh <- df.dh[, -2]
tbl_df(df.dh)

# Reorder this dataframe (it's a mess!)
tmp <- c(1, 2, 3, 5, 4, 20, 19, 6:18)
df.dh <- df.dh[, tmp]

# Convert treatments to factors
df.dh[, 1:7] <- lapply(df.dh[, 1:7], factor)



# ----- Correlation ---------------
# Correlation between Total number of kernels and average kernel wt
# ---------------------------------
plot(Average.kernel.wt~ Total.kernel.per.ear, data=df.dh)
abline(lm(Average.kernel.wt~ Total.kernel.per.ear, data=df.dh))

# Make Model for the data using lm function
kernel.reg <- lm(Average.kernel.wt~ Total.kernel.per.ear, data=df.dh)

# Model summary
summary(kernel.reg)

# Overall model results - ANOVA table
anova(kernel.reg)




# ------------
# Scatterplots
# ------------

plot(No.tassel ~ X.1, data=df.dh)
abline(lm(No.tassel ~ X.1, data=df.dh))

plot(Leafwidth.cm. ~ X.1, data=df.dh)
abline(lm(Leafwidth.cm. ~ X.1, data=df.dh))

plot(Leaflength.cm. ~ X, data=df.dh)
abline(lm(Leaflength.cm. ~ X, data=df.dh))

plot(Plant.Height.cm. ~ X, data=df.dh)
abline(lm(Plant.Height.cm. ~ X, data=df.dh))

plot(ear.position ~ X, data=df.dh)
abline(lm(ear.position ~ X, data=df.dh))

plot(no..nodes ~ X, data=df.dh)
abline(lm(no..nodes ~ X, data=df.dh))

plot(No..of.rows ~ X, data=df.dh)
abline(lm(No..of.rows ~ X, data=df.dh))

plot(Total.kernel.per.ear ~ X, data=df.dh)
abline(lm(Total.kernel.per.ear ~ X, data=df.dh))

plot(ear.length ~ X, data=df.dh)
abline(lm(ear.length ~ X, data=df.dh))

plot(Average.kernel.wt ~ X, data=df.dh)
abline(lm(Average.kernel.wt ~ X, data=df.dh))

plot(ear.cirum. ~ X.1, data=df.dh)
abline(lm(ear.cirum. ~ X, data=df.dh))

plot(dayspol ~ X, data=df.dh)
abline(lm(dayspol  ~ X, data=df.dh))

plot(dayssilk ~ X, data=df.dh)
abline(lm(dayssilk ~ X, data=df.dh))


