# -----------------------------------------------------
# Merritt Burch
# Editing Data frames
# Adding columns of EV's to Double Haploid Data (2015)
# 11/09/2016
#------------------------------------------------------

# Refresh environment
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/Data Sets")

# Import Dataframe
df.15 <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH15_raw data_edited.csv")

# Check data
str(df.15)

# convert treatments into factors
df.15[, 1:5] <- lapply(df.15[, 1:5], as.character)

# Remove columns
df.15 <- df.15[ , c(-1, -5)]


# ----------------------------------
# Create column for treatment names
# ----------------------------------

df.15$lineage.generation <- df.15$Trt

# For lineage.generation name reformatting
df.15$lineage.generation[df.15$lineage.generation == "22-5"] <- "01.05"
df.15$lineage.generation[df.15$lineage.generation == "22-6"] <- "01.06"
df.15$lineage.generation[df.15$lineage.generation == "22-7"] <- "01.07"
df.15$lineage.generation[df.15$lineage.generation == "23-5"] <- "02.05"
df.15$lineage.generation[df.15$lineage.generation == "23-6"] <- "02.06"
df.15$lineage.generation[df.15$lineage.generation == "23-7"] <- "02.07"
df.15$lineage.generation[df.15$lineage.generation == "24-5"] <- "03.05"
df.15$lineage.generation[df.15$lineage.generation == "24-9"] <- "03.09"
df.15$lineage.generation[df.15$lineage.generation == "24-10"] <- "03.10"
df.15$lineage.generation[df.15$lineage.generation == "24-11"] <- "03.11"
df.15$lineage.generation[df.15$lineage.generation == "25-5"] <- "04.05"
df.15$lineage.generation[df.15$lineage.generation == "25-6"] <- "04.06"
df.15$lineage.generation[df.15$lineage.generation == "25-8"] <- "04.08"
df.15$lineage.generation[df.15$lineage.generation == "26-5"] <- "05.05"
df.15$lineage.generation[df.15$lineage.generation == "26-6"] <- "05.05"
df.15$lineage.generation[df.15$lineage.generation == "27-5"] <- "06.05"
df.15$lineage.generation[df.15$lineage.generation == "27-6"] <- "06.06"
df.15$lineage.generation[df.15$lineage.generation == "28-5"] <- "07.05"
df.15$lineage.generation[df.15$lineage.generation == "28-6"] <- "07.06"
df.15$lineage.generation[df.15$lineage.generation == "29-6"] <- "08.06"
df.15$lineage.generation[df.15$lineage.generation == "29-7"] <- "08.07"
df.15$lineage.generation[df.15$lineage.generation == "29-9"] <- "08.09"
df.15$lineage.generation[df.15$lineage.generation == "29-10"] <- "08.10"
df.15$lineage.generation[df.15$lineage.generation == "29-11"] <- "08.11"
df.15$lineage.generation[df.15$lineage.generation == "30-5"] <- "09.05"
df.15$lineage.generation[df.15$lineage.generation == "30-6"] <- "09.06"
df.15$lineage.generation[df.15$lineage.generation == "30-7"] <- "09.07"
df.15$lineage.generation[df.15$lineage.generation == "31-5"] <- "09.05"
df.15$lineage.generation[df.15$lineage.generation == "31-6"] <- "09.06"


# -------------------------------------------
# Create column for modified old.lingen names
# -------------------------------------------

df.15$old.lingen <- df.15$Trt

# For old.lingen name reformatting
df.15$old.lingen[df.15$old.lingen == "22-5"] <- "22.05"
df.15$old.lingen[df.15$old.lingen == "22-6"] <- "22.06"
df.15$old.lingen[df.15$old.lingen == "22-7"] <- "22.07"
df.15$old.lingen[df.15$old.lingen == "23-5"] <- "23.05"
df.15$old.lingen[df.15$old.lingen == "23-6"] <- "23.06"
df.15$old.lingen[df.15$old.lingen == "23-7"] <- "23.07"
df.15$old.lingen[df.15$old.lingen == "24-5"] <- "24.05"
df.15$old.lingen[df.15$old.lingen == "24-9"] <- "24.09"
df.15$old.lingen[df.15$old.lingen == "24-10"] <- "24.10"
df.15$old.lingen[df.15$old.lingen == "24-11"] <- "24.11"
df.15$old.lingen[df.15$old.lingen == "25-5"] <- "25.05"
df.15$old.lingen[df.15$old.lingen == "25-6"] <- "25.06"
df.15$old.lingen[df.15$old.lingen == "25-8"] <- "25.08"
df.15$old.lingen[df.15$old.lingen == "26-5"] <- "26.05"
df.15$old.lingen[df.15$old.lingen == "26-6"] <- "26.05"
df.15$old.lingen[df.15$old.lingen == "27-5"] <- "27.05"
df.15$old.lingen[df.15$old.lingen == "27-6"] <- "27.06"
df.15$old.lingen[df.15$old.lingen == "28-5"] <- "28.05"
df.15$old.lingen[df.15$old.lingen == "28-6"] <- "28.06"
df.15$old.lingen[df.15$old.lingen == "29-6"] <- "29.06"
df.15$old.lingen[df.15$old.lingen == "29-7"] <- "29.07"
df.15$old.lingen[df.15$old.lingen == "29-9"] <- "29.09"
df.15$old.lingen[df.15$old.lingen == "29-10"] <- "29.10"
df.15$old.lingen[df.15$old.lingen == "29-11"] <- "29.11"
df.15$old.lingen[df.15$old.lingen == "30-5"] <- "30.05"
df.15$old.lingen[df.15$old.lingen == "30-6"] <- "30.06"
df.15$old.lingen[df.15$old.lingen == "30-7"] <- "30.07"
df.15$old.lingen[df.15$old.lingen == "31-5"] <- "31.05"
df.15$old.lingen[df.15$old.lingen == "31-6"] <- "31.06"


# ----------------------------------
# Create column for treatment names
# ----------------------------------

df.15$lin <- df.15$Trt

df.15$lin[df.15$lin == "22-5"] <- "01"
df.15$lin[df.15$lin == "22-6"] <- "01"
df.15$lin[df.15$lin == "22-7"] <- "01"
df.15$lin[df.15$lin == "23-5"] <- "02"
df.15$lin[df.15$lin == "23-6"] <- "02"
df.15$lin[df.15$lin == "23-7"] <- "02"
df.15$lin[df.15$lin == "24-5"] <- "03"
df.15$lin[df.15$lin == "24-9"] <- "03"
df.15$lin[df.15$lin == "24-10"] <- "03"
df.15$lin[df.15$lin == "24-11"] <- "03"
df.15$lin[df.15$lin == "25-5"] <- "04"
df.15$lin[df.15$lin == "25-6"] <- "04"
df.15$lin[df.15$lin == "25-8"] <- "04"
df.15$lin[df.15$lin == "26-5"] <- "05"
df.15$lin[df.15$lin == "26-6"] <- "05"
df.15$lin[df.15$lin == "27-5"] <- "06"
df.15$lin[df.15$lin == "27-6"] <- "06"
df.15$lin[df.15$lin == "28-5"] <- "07"
df.15$lin[df.15$lin == "28-6"] <- "07"
df.15$lin[df.15$lin == "29-6"] <- "08"
df.15$lin[df.15$lin == "29-7"] <- "08"
df.15$lin[df.15$lin == "29-9"] <- "08"
df.15$lin[df.15$lin == "29-10"] <- "08"
df.15$lin[df.15$lin == "29-11"] <- "08"
df.15$lin[df.15$lin == "30-5"] <- "09"
df.15$lin[df.15$lin == "30-6"] <- "09"
df.15$lin[df.15$lin == "30-7"] <- "09"
df.15$lin[df.15$lin == "31-5"] <- "09"
df.15$lin[df.15$lin == "31-6"] <- "09"

# -------------
# Add rate data
# -------------

df.15$ratepol  <- 1/df.15$dayspol
df.15$ratesilk <- 1/df.15$dayssilk


# ------------------
# Reorder Dataframe
# ------------------

tmp <- c(1:3, 6, 14:15, 7:13, 4:5)
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


