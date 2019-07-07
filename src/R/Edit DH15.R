# -----------------------------------------------------
# Merritt Burch
# Editing Data frames
# Adding columns of EV's to Double Haploid Data (2015)
# 11/09/2016
#------------------------------------------------------

# Refresh environment
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data")

# Import Dataframe
df.15 <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH raw data- 2015.csv")

# Check data
str(df.dh15)

# convert treatments into factors
df.15[, 1:3] <- lapply(df.15[, 1:3], as.character)


# ---------------
# Remove columns
# ---------------

df.15 <- df.15[, c(-12, -13)]


# ----------------------------------
# Create column for treatment names
# ----------------------------------

df.15$newtreat2 <- df.15$Trt

# For newtreat2 name reformatting
df.15$newtreat2[df.15$newtreat2 == "22-5"] <- "00.05"
df.15$newtreat2[df.15$newtreat2 == "22-6"] <- "00.06"
df.15$newtreat2[df.15$newtreat2 == "22-7"] <- "00.07"
df.15$newtreat2[df.15$newtreat2 == "23-5"] <- "01.05"
df.15$newtreat2[df.15$newtreat2 == "23-6"] <- "01.06"
df.15$newtreat2[df.15$newtreat2 == "23-7"] <- "01.07"
df.15$newtreat2[df.15$newtreat2 == "24-5"] <- "02.05"
df.15$newtreat2[df.15$newtreat2 == "24-9"] <- "02.09"
df.15$newtreat2[df.15$newtreat2 == "24-10"] <- "02.10"
df.15$newtreat2[df.15$newtreat2 == "24-11"] <- "02.11"
df.15$newtreat2[df.15$newtreat2 == "25-5"] <- "03.05"
df.15$newtreat2[df.15$newtreat2 == "25-6"] <- "03.06"
df.15$newtreat2[df.15$newtreat2 == "25-8"] <- "03.08"
df.15$newtreat2[df.15$newtreat2 == "26-5"] <- "04.05"
df.15$newtreat2[df.15$newtreat2 == "26-6"] <- "04.05"
df.15$newtreat2[df.15$newtreat2 == "27-5"] <- "05.05"
df.15$newtreat2[df.15$newtreat2 == "27-6"] <- "05.06"
df.15$newtreat2[df.15$newtreat2 == "28-5"] <- "06.05"
df.15$newtreat2[df.15$newtreat2 == "28-6"] <- "06.06"
df.15$newtreat2[df.15$newtreat2 == "29-6"] <- "07.06"
df.15$newtreat2[df.15$newtreat2 == "29-7"] <- "07.07"
df.15$newtreat2[df.15$newtreat2 == "29-9"] <- "07.09"
df.15$newtreat2[df.15$newtreat2 == "29-10"] <- "07.10"
df.15$newtreat2[df.15$newtreat2 == "29-11"] <- "07.11"
df.15$newtreat2[df.15$newtreat2 == "30-5"] <- "08.05"
df.15$newtreat2[df.15$newtreat2 == "30-6"] <- "08.06"
df.15$newtreat2[df.15$newtreat2 == "30-7"] <- "08.07"
df.15$newtreat2[df.15$newtreat2 == "31-5"] <- "09.05"
df.15$newtreat2[df.15$newtreat2 == "31-6"] <- "09.06"


# ------------------
# Reorder Dataframe
# ------------------

tmp <- c(1, 2, 3, 12, 4, 5:11)
df.15 <- df.15[, tmp]


# -------------------
# Remove missing data
# --------------------

# Removes 16 rows
df.15 <- df.15[complete.cases(df.15),]


# --------------
# Rename columns
# --------------

names(df.15)[names(df.15)=="Numtassel"] <- "No.tassel"
names(df.15)[names(df.15)=="Height"] <- "Plant.Height.cm."
names(df.15)[names(df.15)=="Leaflength"] <- "Leaflength.cm."
names(df.15)[names(df.15)=="Earpos"] <- "ear.position"
names(df.15)[names(df.15)=="Dayspoll"] <- "dayspol"
names(df.15)[names(df.15)=="Dayssilkk"] <- "dayssilk"
names(df.15)[names(df.15)=="Numnodes"] <- "no..nodes"

# -------------
# Save Dataset
# -------------
# Save as (CSV)
write.csv(df.15, file = "DH15_raw data_edited.csv")


