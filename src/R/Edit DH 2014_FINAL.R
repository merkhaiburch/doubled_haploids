# -----------------------------------------------------------------------------
# Editing of the DH14 data
# Merritt Burch
# 6/21/2017

# I edited this data before, as seen in the 'Edit DH 2014' Script
# I realize now that I need to edit it more, changes listed below
# Permanently add the rate of silk/pol
# Delete first columns (added by default by excel)
# Addlineage information (my environmental lineages)
# -----------------------------------------------------------------------------

# Load data
df.dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH14_BorisEdit.csv")

# Rename columns
names(df.dh)[names(df.dh) == "newtreat2"] <- "old.lingen"
names(df.dh)[names(df.dh) == "newtreat3"] <- "lineage.generation"
names(df.dh)[names(df.dh) == "newtreat4"] <- "lineage"

# Remove extra rows
df.dh <- df.dh[, c(-1,-8)]

# Convert treatments to factors
df.dh[, 1:7] <- lapply(df.dh[, 1:6], factor)

# Add columns for rate
df.dh$rate.silk <- 1/df.dh$dayssilk
df.dh$rate.pol  <- 1/df.dh$dayspol

# ------------------------
# Add lineage information
# ------------------------

# Create new column with all levels from lineage.generation
df.dh$env.lin <- df.dh$lineage.generation

# Turn factors into characters for easier conversion 
df.dh$env.lin <- as.character(df.dh$env.lin)

# Find and replace lineage.generation info in the new env.lin column
#   based on its environmental history/phylogeny
df.dh$env.lin[df.dh$env.lin == "1.05"] <- "HI09.SD14"
df.dh$env.lin[df.dh$env.lin == "2.05"] <- "HI09.SD14"
df.dh$env.lin[df.dh$env.lin == "2.06"] <- "HI09.SD12.SD14"
df.dh$env.lin[df.dh$env.lin == "3.05"] <- "HI09.SD14"
df.dh$env.lin[df.dh$env.lin == "3.09"] <- "HI09.SD09.HI10.SD10.HI11.SD14"
df.dh$env.lin[df.dh$env.lin == "3.1"] <- "HI09.SD09.HI10.SD10.HI11.SD11.SD14"
df.dh$env.lin[df.dh$env.lin == "4.05"] <- "HI09.SD14"
df.dh$env.lin[df.dh$env.lin == "4.08"] <- "HI09.SD09.HI10.SD11.SD14"
df.dh$env.lin[df.dh$env.lin == "5.05"] <- "HI09.SD14"
df.dh$env.lin[df.dh$env.lin == "6.05"] <- "HI09.SD14"
df.dh$env.lin[df.dh$env.lin == "6.06"] <- "HI09.SD13.SD14"
df.dh$env.lin[df.dh$env.lin == "7.05"] <- "HI09.SD14"
df.dh$env.lin[df.dh$env.lin == "8.06"] <- "HI09.SD11.SD14"
df.dh$env.lin[df.dh$env.lin == "8.09"] <- "HI09.SD09.HI10.SD10.HI11.SD14"
df.dh$env.lin[df.dh$env.lin == "8.1"] <- "HI09.SD09.HI10.SD10.HI11.SD11.SD14"
df.dh$env.lin[df.dh$env.lin == "9.05"] <- "HI09.SD14"
df.dh$env.lin[df.dh$env.lin == "9.06"] <- "HI09.SD11.SD14"
df.dh$env.lin[df.dh$env.lin == "10.05"] <- "HI09.SD14"
df.dh$env.lin[df.dh$env.lin == "10.06"] <- "HI09.SD13.SD14"

# Turn this column into factors again
df.dh$env.lin <- as.factor(df.dh$env.lin)

# Reorder this dataframe (it's a mess!)
tmp <- c(1:6, 22, 7:18, 21, 19,20)
df.dh <- df.dh[, tmp]

# Save this updated dataframe for future use
write.csv(df.dh, file = "DH14_FINAL.csv")
