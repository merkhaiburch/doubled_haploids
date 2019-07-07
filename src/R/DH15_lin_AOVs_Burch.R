# ----------------------------------------------------
# Merritt Does things with data!
# Double Haploid Linear mixed models & ANOVA Analysis
# Summer 2015 Data
# 11/09/2016
# ----------------------------------------------------

# Clear Directory
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/DH_15")

# Import Data
df.15 <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH15_raw data_edited.csv")

# Load packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'doBy', 'lmtest', 'lsmeans', 'agricolae', 'tidyr')
pack.man(packages)


# Remove column
df.15 <- df.15[, -1]

# Convert treatments to factors
df.15[, 1:6] <- lapply(df.15[, 1:6], factor)

# Summary statistics of data
stat <- lapply(df.15[,c(7:11, 14:15)] , 
               function(x) rbind( mean = mean(x),
                                  sd = sd(x),
                                  median = median(x),
                                  minimum = min(x),
                                  maximum = max(x),
                                  n = length(x) ) )
# Limit number of significant figures
summary <- data.frame(stat)
(summary <- signif(summary, digits = 3))

# transpose data
(summary <- t(summary))

# ----------------------- No.tassel --------------------------

# --------
# Boxplot
# --------

# All blocks
ggplot(aes(y = No.tassel, x = lin), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Number of Tassels") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_No.tassel_boxplot all blocks.png", width = 10, height = 6)


# By block
ggplot(aes(y = No.tassel, x = lin), data = df.15) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Number of Tassels") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_No.tassel_boxplot each block.png", width = 11, height = 6)


# --------------
# Density plots
# --------------

# All Blocks
ggplot(df.15, aes(x=No.tassel)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.15, aes(No.tassel)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# -------------------------
#  Summary Stats and Graphs
# -------------------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.No.tassel <- summaryBy(No.tassel ~ lin, data = df.15, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.No.tassel$se <- summary.No.tassel$No.tassel.sd/(sqrt(summary.No.tassel$No.tassel.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
No.tassel.2way <- aov(No.tassel ~ lin + Block + lin*Block, data = df.15)
anova(No.tassel.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(No.tassel.2way)

# Creates Post-hoc groupings using agricolate
#     Orders treatments in numberical order
#     Adds column of post-hoc groupings to summary table for barchart
(tpoho <- HSD.test(No.tassel.2way, "lin", group=TRUE))
nt <- data.frame(tpoho$groups)
nt <- nt[order(nt$trt),]
summary.No.tassel$poho <- nt$M


# --------------------------
# Linear mixed effect models
# --------------------------
rNo.tassel <- lmer(No.tassel ~ lin + (1|Block), data = df.15)
anova(rNo.tassel)
summary(rNo.tassel)

# To check unexplained variability
cNo.tassel <- lmer(No.tassel ~ (1| lin) + (1|Block) + (1|year), data = df.15)
summary(cNo.tassel)

# Post hoc test
lsm <- lsmeans(rNo.tassel, pairwise ~ lin, adjust="tukey")

# cld is a function used to extract and display information 
phnt <- cld(lsm, alpha=.05, Letters=letters) 
phnt2 <- subset(phnt, select=c(lin, .group))
phnt2 <- phnt2[order(phnt2$lin),]
summary.No.tassel$poho <- phnt2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.No.tassel$se.mean.p <- c(summary.No.tassel$No.tassel.mean + summary.No.tassel$se + 0.6)

# Make Graph
ggplot(summary.No.tassel, aes(lin, No.tassel.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = No.tassel.mean - se, ymax = No.tassel.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Number of Tassels") +
  geom_text(data = summary.No.tassel, aes(x = lin, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_No.tassel_barchart with post hoc.png", width = 10.5, height = 6)


# ------ Leaflength.cm. -------------------------------------------------------

# ---------
# Boxplots
# ---------

# Error in row 22 of data
# Remove for this treatment analysis
ll.temp <- df.15[-22,]

# All blocks
ggplot(aes(y = Leaflength.cm., x = lin), data = ll.temp) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Leaf Length (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_Leaflength.cm._boxplot all blocks.png", width = 10, height = 6)

# Each block  
ggplot(aes(y = Leaflength.cm., x = lin), data = ll.temp) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Leaf Length (cm)")  +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_Leaflength.cm._boxplot each block.png", width = 10.5, height = 6)


# ------------- 
# Density plots
# -------------

# All Blocks
ggplot(ll.temp, aes(x=Leaflength.cm.)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(ll.temp, aes(Leaflength.cm.)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

# Summary Stats
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.Leaflength.cm. <- summaryBy(Leaflength.cm. ~ lin, data = ll.temp, na.rm = TRUE,
                                    FUN = c(mean, sd, length2))
summary.Leaflength.cm.$se <- summary.Leaflength.cm.$Leaflength.cm..sd/(sqrt(summary.Leaflength.cm.$Leaflength.cm..length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
Leaflength.cm..2way <- aov(Leaflength.cm. ~ lin + Block + lin*Block, data = ll.temp)
anova(Leaflength.cm..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(Leaflength.cm..2way)

# To check unexplained variability
cLeaflength.cm. <- lmer(Leaflength.cm. ~ (1| lin) + (1|Block) + (1|year), data = ll.temp)
summary(cLeaflength.cm.)

# Linear Mixed Effect Model
rLeaflength.cm. <- lmer(Leaflength.cm. ~ lin + (1 | Block), data = ll.temp)
anova(rLeaflength.cm.)
summary(rLeaflength.cm.)


# ------------- 
# Post hoc test
# NONE!
# -------------


# -----------------
# Post-Hoc graphing
# -----------------


# Make the graph
ggplot(summary.Leaflength.cm., aes(lin, Leaflength.cm..mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = Leaflength.cm..mean - se, ymax = Leaflength.cm..mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Leaf Length (cm)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_Leaflength.cm._barchart with post hoc.png", width = 10, height = 6)



# ------ Plant.Height.cm. ------------------------------------------------------

# --------
# Boxplots
# --------

# All groups
ggplot(aes(y = Plant.Height.cm., x = lin), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Plant Height (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_Plant.Height.cm._boxplot all groups.png", width = 10, height = 6)


# Each group
ggplot(aes(y = Plant.Height.cm., x = lin), data = df.15) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Plant Height (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_Plant.Height.cm._boxplot each group.png", width = 10, height = 6)


# --------------
# Density plots
# --------------

# All Blocks
ggplot(df.15, aes(x=Plant.Height.cm.)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.15, aes(Plant.Height.cm.)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.Plant.Height.cm. <- summaryBy(Plant.Height.cm. ~ lin, data = df.15, na.rm = TRUE,
                                      FUN = c(mean, sd, length2))
summary.Plant.Height.cm.$se <- summary.Plant.Height.cm.$Plant.Height.cm..sd/(sqrt(summary.Plant.Height.cm.$Plant.Height.cm..length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
Plant.Height.cm..2way <- aov(Plant.Height.cm. ~ lin + Block + lin*Block, data = df.15)
anova(Plant.Height.cm..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(Plant.Height.cm..2way)

# To check unexplained variability
cPlant.Height.cm. <- lmer(Plant.Height.cm. ~ (1| lin) + (1|Block) + (1|year), data = df.15)
summary(cPlant.Height.cm.)

# Linear Mixed Effect Modeling
rPlant.Height.cm. <- lmer(Plant.Height.cm. ~ lin + (1 | Block), data = df.15)
anova(rPlant.Height.cm.)
summary(rPlant.Height.cm.)


# -------------
# Post hoc test
# -------------

lsm <- lsmeans(rPlant.Height.cm., pairwise~lin, adjust="tukey")

# cld is a function used to extract and display information 
ph.ph <- cld(lsm, alpha=.05, Letters=letters) 
ph.ph2 <- subset(ph.ph, select=c(lin, .group))
ph.ph2 <- ph.ph2[order(ph.ph2$lin),]
summary.Plant.Height.cm.$poho <- ph.ph2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.Plant.Height.cm.$se.mean.p <- c(summary.Plant.Height.cm.$Plant.Height.cm..mean + summary.Plant.Height.cm.$se + 7)

# Make the graph
ggplot(summary.Plant.Height.cm., aes(lin, Plant.Height.cm..mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = Plant.Height.cm..mean - se, ymax = Plant.Height.cm..mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Plant Height (cm)") +
  geom_text(data = summary.Plant.Height.cm., aes(x = lin, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_Plant.Height.cm._barchart with post hoc.png", width = 10.5, height = 6)


# ------ ear.position ----------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ear.position, x = lin), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Ear Position") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_ear.position_boxplot all blocks.png", width = 10, height = 6)

# Each Block
ggplot(aes(y = ear.position, x = lin), data = df.15) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Ear Position")  +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_ear.position_boxplot each block.png", width = 10, height = 6)


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.15, aes(x=ear.position)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.15, aes(ear.position)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# ---------------
#  Summary Stats
# ---------------
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.ear.position <- summaryBy(ear.position ~ lin, data = df.15, na.rm = TRUE,
                                  FUN = c(mean, sd, length2))
summary.ear.position$se <- summary.ear.position$ear.position.sd/(sqrt(summary.ear.position$ear.position.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ear.position.2way <- aov(ear.position ~ lin + Block + lin*Block, data = df.15)
anova(ear.position.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(ear.position.2way)

# To check unexplained variability
cear.position <- lmer(ear.position ~ (1|lin) + (1|Block) + (1|year), data = df.15)
summary(cear.position)

# Linear Mixed Effect Model
rear.position <- lmer(ear.position ~ lin + (1 | Block), data = df.15)
anova(rear.position)
summary(rear.position)


# -------------
# NO hoc test
# -------------

# Make the graph
ggplot(summary.ear.position, aes(lin, ear.position.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ear.position.mean - se, ymax = ear.position.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Ear Position") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_ear.position_barchart with NO post hoc.png", width = 10, height = 6)


# ------ no..nodes ------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = no..nodes, x = lin), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Number of Nodes") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_no..nodes_boxplot all blocks.png", width = 10, height = 6)

# Each block
ggplot(aes(y = no..nodes, x = lin), data = df.15) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Number of Nodes")  +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_no..nodes_boxplot each block.png", width = 10, height = 6)

# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.15, aes(x=no..nodes)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.15, aes(no..nodes)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)

# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.no..nodes <- summaryBy(no..nodes ~ lin, data = df.15, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.no..nodes$se <- summary.no..nodes$no..nodes.sd/(sqrt(summary.no..nodes$no..nodes.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
no..nodes.2way <- aov(no..nodes ~ lin + Block + lin*Block, data = df.15)
anova(no..nodes.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(no..nodes.2way)

# To check unexplained variability
cno..nodes <- lmer(no..nodes ~ (1| lin) + (1|Block) + (1| year), data = df.15)
summary(cno..nodes)

# Linear Mixed Effect Modeling
rno..nodes <- lmer(no..nodes ~ lin + (1 | Block), data = df.15)
anova(rno..nodes)
summary(rno..nodes)


# --------------
# NO Post hoc test
# --------------

# Make the graph
ggplot(summary.no..nodes, aes(lin, no..nodes.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = no..nodes.mean - se, ymax = no..nodes.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Number of Nodes") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_no..nodes_barchart with post hoc.png", width = 10, height = 6)


# ------ ratepol --------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ratepol, x = lin), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Rate of Pollen Shed") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_ratepol_boxplot all blocks.png", width = 10, height = 6)

# Each block
ggplot(aes(y = ratepol, x = lin), data = df.15) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Rate of Pollen Shed") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_ratepol_boxplot each block.png", width = 10, height = 6)


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.15, aes(x=ratepol)) + 
  geom_histogram(aes(y =..density..), binwidth = 0.0008, colour="black", fill="white") +
  geom_density(alpha = 0.5, fill = "#FF6666")
ggsave(file = "DH15_ratepol_histogram_all block.png", width = 10, height = 6)

# Each Block
ggplot(df.15, aes(ratepol)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.5)
ggsave(file = "DH15_ratepol_histogram_each block.png", width = 10, height = 6)


# --------------
# Summary Stats
# --------------

summary.ratepol <- summaryBy(ratepol ~ lin, data = df.15, na.rm = TRUE,
                             FUN = c(mean, sd, length2))
summary.ratepol$se <- summary.ratepol$ratepol.sd/(sqrt(summary.ratepol$ratepol.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ratepol.2way <- aov(ratepol ~ lin + Block + lin*Block, data = df.15)
anova(ratepol.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(ratepol.2way)

# To check unexplained variability
cratepol <- lmer(ratepol ~ (1| lin) + (1|Block) + (1|year), data = df.15)
summary(cratepol)

# Linear Mixed Effect Model
rratepol <- lmer(ratepol ~ lin + (1 | Block), data = df.15)
anova(rratepol)
summary(rratepol)


# --------------
# NO Post hoc test
# --------------

# Make the graph
ggplot(summary.ratepol, aes(lin, ratepol.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ratepol.mean - se, ymax = ratepol.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Rate of Pollen Shed") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_ratepol_barchart with post hoc.png", width = 10, height = 6)



# ------ ratesilk --------------------------------------------------------------

# --------
# Boxplots
# --------

# All groups
ggplot(aes(y = ratesilk, x = lin), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Days to Silk Emergence") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_ratesilk_boxplot all groups.png", width = 10, height = 6)

# Each group
ggplot(aes(y = ratesilk, x = lin), data = df.15) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Days to Silk Emergence") +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_ratepol_boxplot_each block.png", width = 10, height = 6)


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.15, aes(x=ratesilk)) + 
  geom_histogram(aes(y=..density..), binwidth=0.0008, colour="black", fill="white") +
  geom_density(alpha=0.5, fill="#FF6666")
ggsave(file = "DH15_ratesilk_hist_all groups.png", width = 10, height = 6)

# Each Block
ggplot(df.15, aes(ratesilk)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.5)
ggsave(file = "DH15_ratesilk_hist_each block.png", width = 10, height = 6)


# --------------
#  Summary Stats
# --------------

summary.ratesilk <- summaryBy(ratesilk ~ lin, data = df.15, na.rm = TRUE,
                              FUN = c(mean, sd, length2))
summary.ratesilk$se <- summary.ratesilk$ratesilk.sd/(sqrt(summary.ratesilk$ratesilk.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
ratesilk.2way <- aov(ratesilk ~ lin + Block + lin*Block, data = df.15)
anova(ratesilk.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(ratesilk.2way)

# To check unexplained variability
cratesilk <- lmer(ratesilk ~ (1| lin) + (1|Block) + (1|year), data = df.15)
summary(cratesilk)

# Linear Mixed Effect Model
rratesilk <- lmer(ratesilk ~ lin + (1 | Block), data = df.15)
anova(rratesilk)
summary(rratesilk)


# -------------
# NO Post hoc test
# --------------

# Make the graph
ggplot(summary.ratesilk, aes(lin, ratesilk.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ratesilk.mean - se, ymax = ratesilk.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Rate of Silk Emergence") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH15_lin_ratesilk_barchart with post hoc.png", width = 10, height = 6)




