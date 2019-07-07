# -----------------------------------------------
# Merritt Burch
# Editing Data frames
# Adding columns of EV's to Double Haploid Data
# 9/9/2016
# Adding location grown and year ear produced
#------------------------------------------------

# Refresh environment
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data")

# Import Dataframe
dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH raw data -2014_edited_V2.csv")

# Check data
str(dh)

# convert treatments into factors
dh[, c(3:4, 6:7)] <- lapply(dh[, c(3:4, 6:7)], factor)


# --------------------
# Correct misspellings
# --------------------

names(dh)[names(dh)=="ear.postion"] <- "ear.position"
names(dh)[names(dh)=="eal.length"] <- "ear.length"


# --------------------
# Add column for year
# --------------------

# Create two columns of replicate data names to change
# Years ears were produced and location in which they were grown
dh$year <- dh$Treatment
dh$year <- as.character(dh$year)

# For Year (all 2009 years were produced in Hawaii)
dh$year[dh$year == "22-5"] <- "2009"
dh$year[dh$year == "22-6"] <- "2012"
dh$year[dh$year == "22-7"] <- "2010"
dh$year[dh$year == "23-5"] <- "2009"
dh$year[dh$year == "23-6"] <- "2012"
dh$year[dh$year == "24-5"] <- "2009"
dh$year[dh$year == "24-9"] <- "2011"
dh$year[dh$year == "24-10"] <- "2011"
dh$year[dh$year == "25-5"] <- "2009"
dh$year[dh$year == "25-8"] <- "2011"
dh$year[dh$year == "26-5"] <- "2009"
dh$year[dh$year == "26-6"] <- "2011"
dh$year[dh$year == "27-5"] <- "2009"
dh$year[dh$year == "27-6"] <- "2013"
dh$year[dh$year == "28-5"] <- "2009"
dh$year[dh$year == "28-6"] <- "2013"
dh$year[dh$year == "29-6"] <- "2011"
dh$year[dh$year == "29-9"] <- "2011"
dh$year[dh$year == "29-10"] <- "2011"
dh$year[dh$year == "30-5"] <- "2009"
dh$year[dh$year == "30-6"] <- "2011"
dh$year[dh$year == "31-5"] <- "2009"
dh$year[dh$year == "31-6"] <- "2013"


# ---------------
# Remove columns
# ---------------

# Remove column of hundred kernel weight
dh$hundred.grain.wt. <- NULL


# ----------------
# Remove outliers
# ----------------

# Remove two outrageous leaf width measurement (not 90cm!) (341,418)
# Remove three number of rows (ear) b/c Don checked them and they had more rows 
#     Also b/c they're abnormally small (346,398,533)
# Row 358 b/c of abnormally low total kernel number (could be bad pollination)
# Rows 199,358,500 b/c they have 0 as a weight, rest of the columns (other traits)
#       have data present but going to delete these rows to make data consistent 
#       (without any blanks)
dh <- dh[-c(199,341,346,358,558,398,418,500,533), ]


# Remove all rows with NA data, most to all of the NA rows have nothing 
#     in them so it isnt an issue
# This also removes many of the kenel weights that were 0 (maybe 0 by accident?)
dh <- dh[complete.cases(dh),]


# --------------
# Save dataframe
# --------------

# As new (CSV)
write.csv(dh, file = "DH raw data -2014_edited.csv")


# ------------------
# Add another column
# ------------------

# Work with latest updatad dataset
dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH raw data -2014_edited_V2.csv")

# Convert to factor
dh$Treatment <- factor(dh$Treatment)

# Create column
dh$newtreat2 <- dh$Treatment
dh$newtreat2 <- as.character(dh$newtreat2)

# For newtreat2 name reformatting
dh$newtreat2[dh$newtreat2 == "22-5"] <- "22.05"
dh$newtreat2[dh$newtreat2 == "22-6"] <- "22.06"
dh$newtreat2[dh$newtreat2 == "22-7"] <- "22.07"
dh$newtreat2[dh$newtreat2 == "23-5"] <- "23.05"
dh$newtreat2[dh$newtreat2 == "23-6"] <- "23.06"
dh$newtreat2[dh$newtreat2 == "24-5"] <- "24.05"
dh$newtreat2[dh$newtreat2 == "24-9"] <- "24.09"
dh$newtreat2[dh$newtreat2 == "24-10"] <- "24.10"
dh$newtreat2[dh$newtreat2 == "25-5"] <- "25.05"
dh$newtreat2[dh$newtreat2 == "25-8"] <- "25.08"
dh$newtreat2[dh$newtreat2 == "26-5"] <- "26.05"
dh$newtreat2[dh$newtreat2 == "26-6"] <- "26.05"
dh$newtreat2[dh$newtreat2 == "27-5"] <- "27.05"
dh$newtreat2[dh$newtreat2 == "27-6"] <- "27.06"
dh$newtreat2[dh$newtreat2 == "28-5"] <- "28.05"
dh$newtreat2[dh$newtreat2 == "28-6"] <- "28.06"
dh$newtreat2[dh$newtreat2 == "29-6"] <- "29.06"
dh$newtreat2[dh$newtreat2 == "29-9"] <- "29.09"
dh$newtreat2[dh$newtreat2 == "29-10"] <- "29.10"
dh$newtreat2[dh$newtreat2 == "30-5"] <- "30.05"
dh$newtreat2[dh$newtreat2 == "30-6"] <- "30.06"
dh$newtreat2[dh$newtreat2 == "31-5"] <- "31.05"
dh$newtreat2[dh$newtreat2 == "31-6"] <- "31.06"


# Save Dataset
# Save as new dataset (CSV)
write.csv(dh, file = "DH raw data -2014_edited_V2.csv")


