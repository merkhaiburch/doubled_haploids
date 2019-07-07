----------------------------------------------------
  # Merritt Does things with data!
  # Double Haploid Linear mixed models & ANOVA Analysis
  # Summer 2015 Data
  # 11/09/2016
  # ----------------------------------------------------

# Clear Directory
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/Data Sets")

# Import Data
df.15 <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH15_raw data_edited.csv")

# Load packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'doBy', 'lsmeans', 'agricolae', 'lme4')
pack.man(packages)


# ---------------
# Remove columns
# ---------------

df.15 <- df.15[, -1]


# -----------------------------
# Convert treatments to factors
# -----------------------------

df.15[, 1:4] <- lapply(df.15[, 1:4], factor)


# --------------------------------------------------------------------------
# No.tassel
# --------------------------------------------------------------------------

# ---------
# Boxplots
# ---------

# All blocks
a <- ggplot(aes(y = No.tassel, x = newtreat22), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Tassels")

# By block
ggplot(aes(y = No.tassel, x = newtreat22), data = df.15) + 
  geom_boxplot( aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of tassels")


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
summary.No.tassel <- summaryBy(No.tassel ~ newtreat22, data = df.15, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.No.tassel$se <- summary.No.tassel$No.tassel.sd/(sqrt(summary.No.tassel$No.tassel.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
No.tassel.2way <- aov(No.tassel ~ newtreat22 + Block + newtreat22*Block, data = df.15)
anova(No.tassel.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(No.tassel.2way)

# Creates Post-hoc groupings using agricolate
#     Orders treatments in numberical order
#     Adds column of post-hoc groupings to summary table for barchart
(tpoho <- HSD.test(No.tassel.2way, "newtreat22", group=TRUE))
nt <- data.frame(tpoho$groups)
nt <- nt[order(nt$trt),]
summary.No.tassel$poho <- nt$M


# --------------------------
# Linear mixed effect models
# --------------------------
(rNo.tassel <- lmer(No.tassel ~ newtreat22 + (1 | Block), data = df.15))
anova(rNo.tassel)
summary(rNo.tassel)

# To check unexplained variability
cNo.tassel <- lmer(No.tassel ~ (1| newtreat22) + (1|Block), data = df.15)
summary(cNo.tassel)

# Post hoc test
lsm <- lsmeans(rNo.tassel, pairwise ~ newtreat22, adjust="tukey")

# cld is a function used to extract and display information 
phnt <- cld(lsm, alpha=.05, Letters=letters) 
phnt2 <- subset(phnt, select=c(newtreat22, .group))
phnt2 <- phnt2[order(phnt2$newtreat22),]
summary.No.tassel$poho <- phnt2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.No.tassel$se.mean.p <- c(summary.No.tassel$No.tassel.mean + summary.No.tassel$se + 0.6)

# Make Graph
ggplot(summary.No.tassel, aes(newtreat22, No.tassel.mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = No.tassel.mean - se, ymax = No.tassel.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Number of Tassels") +
  geom_text(data = summary.No.tassel, aes(x = newtreat22, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  guides(fill=FALSE)



# ------ Leaflength.cm. -------------------------------------------------------
# One outlier in row 22 that needs to be removed
# Not significant so no post hoc
# -----------------------------------------------------------------------------

df.15 <- df.15[-22, ]

# ---------
# Boxplots
# ---------

# All blocks
ggplot(aes(y = Leaflength.cm., x = newtreat2), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Leaf Length (cm)")

# Each block  
ggplot(aes(y = Leaflength.cm., x = newtreat2), data = df.15) + 
  geom_boxplot( aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Leaf Length (cm)")


# ------------- 
# Density plots
# -------------

# All Blocks
ggplot(df.15, aes(x=Leaflength.cm.)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.15, aes(Leaflength.cm.)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

# Summary Stats
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.Leaflength.cm. <- summaryBy(Leaflength.cm. ~ newtreat2, data = df.15, na.rm = TRUE,
                                    FUN = c(mean, sd, length2))
summary.Leaflength.cm.$se <- summary.Leaflength.cm.$Leaflength.cm..sd/(sqrt(summary.Leaflength.cm.$Leaflength.cm..length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
Leaflength.cm..2way <- aov(Leaflength.cm. ~ newtreat2 + Block + newtreat2*Block, data = df.15)
anova(Leaflength.cm..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Leaflength.cm..2way)

# To check unexplained variability
cLeaflength.cm. <- lmer(Leaflength.cm. ~ (1| newtreat2) + (1|Block), data = df.15)
summary(cLeaflength.cm.)

# Linear Mixed Effect Model
(rLeaflength.cm. <- lmer(Leaflength.cm. ~ newtreat2 + (1 | Block), data = df.15))
anova(rLeaflength.cm.)
summary(rLeaflength.cm.)


# -----------------
# NO Post-Hoc graphing
# -----------------

# Make the graph
ggplot(summary.Leaflength.cm., aes(newtreat2, Leaflength.cm..mean)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = Leaflength.cm..mean - se, ymax = Leaflength.cm..mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Leaf Length (cm)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))



# ------ Plant.Height.cm. ------------------------------------------------------

# ------------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All groups
ggplot(aes(y = Plant.Height.cm., x = newtreat2), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Plant Height (cm)")

# Each group
ggplot(aes(y = Plant.Height.cm., x = newtreat2), data = df.15) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Plant Height (cm)")


# --------------
# Density plots
# --------------

# All Blocks
ggplot(df.15, aes(x=Plant.Height.cm.)) + 
  geom_histogram(aes(y=..density..), binwidth=5, colour="black", fill="white") +
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
summary.Plant.Height.cm. <- summaryBy(Plant.Height.cm. ~ newtreat2, data = df.15, na.rm = TRUE,
                                      FUN = c(mean, sd, length2))
summary.Plant.Height.cm.$se <- summary.Plant.Height.cm.$Plant.Height.cm..sd/(sqrt(summary.Plant.Height.cm.$Plant.Height.cm..length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
Plant.Height.cm..2way <- aov(Plant.Height.cm. ~ newtreat2 + Block + newtreat2*Block, data = df.15)
anova(Plant.Height.cm..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Plant.Height.cm..2way)

# To check unexplained variability
cPlant.Height.cm. <- lmer(Plant.Height.cm. ~ (1| newtreat2) + (1|Block), data = df.15)
summary(cPlant.Height.cm.)

# Linear Mixed Effect Modeling
(rPlant.Height.cm. <- lmer(Plant.Height.cm. ~ newtreat2 + (1 | Block), data = df.15))
anova(rPlant.Height.cm.)
summary(rPlant.Height.cm.)


# -------------
# Post hoc test
# -------------

lsm <- lsmeans(rPlant.Height.cm., pairwise~newtreat2, adjust="tukey")

# cld is a function used to extract and display information 
ph.ph <- cld(lsm, alpha=.05, Letters=letters) 
ph.ph2 <- subset(ph.ph, select=c(newtreat2, .group))
ph.ph2 <- ph.ph2[order(ph.ph2$newtreat2),]
summary.Plant.Height.cm.$poho <- ph.ph2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.Plant.Height.cm.$se.mean.p <- c(summary.Plant.Height.cm.$Plant.Height.cm..mean + summary.Plant.Height.cm.$se + 7)

# Make the graph
ggplot(summary.Plant.Height.cm., aes(newtreat2, Plant.Height.cm..mean)) + 
  geom_bar(stat="identity", aes(fill = poho)) +
  geom_errorbar(aes(ymin = Plant.Height.cm..mean - se, ymax = Plant.Height.cm..mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Plant Height (cm)") +
  geom_text(data = summary.Plant.Height.cm., aes(x = newtreat2, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  guides(fill=FALSE)


# ------ ear.position ----------------------------------------------------------
# Bimodal distribution b/c so few intermediate values 
# Treatment and interaction significant in ANOVA
# High variability explaind by Block and interaction in Random blocking
#    Doing post-hoc based on random blocking with blocking variability removed
# ------------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ear.position, x = newtreat2), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Ear Position")

# Each Block
ggplot(aes(y = ear.position, x = newtreat2), data = df.15) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Ear Position")


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
summary.ear.position <- summaryBy(ear.position ~ newtreat2, data = df.15, na.rm = TRUE,
                                  FUN = c(mean, sd, length2))
summary.ear.position$se <- summary.ear.position$ear.position.sd/(sqrt(summary.ear.position$ear.position.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ear.position.2way <- aov(ear.position ~ newtreat2 + Block + newtreat2*Block, data = df.15)
anova(ear.position.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ear.position.2way)

# To check unexplained variability
cear.position <- lmer(ear.position ~ (1|newtreat2) + (1|Block), data = df.15)
summary(cear.position)

# Linear Mixed Effect Model
(rear.position <- lmer(ear.position ~ newtreat2 + (1 | Block), data = df.15))
anova(rear.position)
summary(rear.position)


# -------------
# Post hoc test
# -------------
lsm <- lsmeans(rear.position, pairwise~newtreat2, adjust="tukey")

# cld is a function used to extract and display information 
ph.ep <- cld(lsm, alpha=.05, Letters=letters) 
ph.ep2 <- subset(ph.ep, select=c(newtreat2, .group))
ph.ep2 <- ph.ep2[order(ph.ep2$newtreat2),]
summary.ear.position$poho <- ph.ep2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.ear.position$se.mean.p <- c(summary.ear.position$ear.position.mean + summary.ear.position$se + .3)

# Make the graph
ggplot(summary.ear.position, aes(newtreat2, ear.position.mean)) + 
  geom_bar(stat="identity", aes(fill = poho)) +
  geom_errorbar(aes(ymin = ear.position.mean - se, ymax = ear.position.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Ear Position") +
  geom_text(data = summary.ear.position, aes(x = newtreat2, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  guides(fill = FALSE)


# ------ no..nodes ------------------------------------------------------------
# Bimodal distribution b/c so few intermediate values 
# Treatment only significant in ANOVA
# Little variability explaind by Block and interaction in Random blocking
#    Doing post-hoc based on random blocking with blocking variability removed
# -----------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = no..nodes, x = newtreat2), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Nodes")

# Each block
ggplot(aes(y = no..nodes, x = newtreat2), data = df.15) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Nodes")


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
summary.no..nodes <- summaryBy(no..nodes ~ newtreat2, data = df.15, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.no..nodes$se <- summary.no..nodes$no..nodes.sd/(sqrt(summary.no..nodes$no..nodes.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
no..nodes.2way <- aov(no..nodes ~ newtreat2 + Block + newtreat2*Block, data = df.15)
anova(no..nodes.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(no..nodes.2way)

# To check unexplained variability
cno..nodes <- lmer(no..nodes ~ (1| newtreat2) + (1|Block), data = df.15)
summary(cno..nodes)

# Linear Mixed Effect Modeling
(rno..nodes <- lmer(no..nodes ~ newtreat2  , data = df.15))
anova(rno..nodes)
summary(rno..nodes)


# --------------
# Post hoc test
# --------------
lsm <- lsmeans(rno..nodes, pairwise~newtreat2, adjust="tukey")

# cld is a function used to extract and display information 
ph.nn <- cld(lsm, alpha=.05, Letters=letters)
ph.nn2 <- subset(ph.nn, select=c(newtreat2, .group))
ph.nn2 <- ph.nn2[order(ph.nn2$newtreat2),]
summary.no..nodes$poho <- ph.nn2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.no..nodes$se.mean.p <- c(summary.no..nodes$no..nodes.mean + summary.no..nodes$se + .7)

# Make the graph
ggplot(summary.no..nodes, aes(newtreat2, no..nodes.mean)) + 
  geom_bar(stat="identity", aes(fill = poho)) +
  geom_errorbar(aes(ymin = no..nodes.mean - se, ymax = no..nodes.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Number of Nodes") +
  geom_text(data = summary.no..nodes, aes(x = newtreat2, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  guides(FALSE)



# ------ dayspol --------------------------------------------------------------
# Kinda Normal distribution (three small values pull dist to left)
# Treatment/block/interaction significant in ANOVA
# Equal variability explaind by Block and interaction in Random blocking
#    Doing post-hoc based on random blocking with blocking variability removed
# -----------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = dayspol, x = newtreat2), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Days to Pollen Shed")

# Each block
ggplot(aes(y = dayspol, x = newtreat2), data = df.15) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Days to Pollen Shed")


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.15, aes(x=dayspol)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.15, aes(dayspol)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
# Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.dayspol <- summaryBy(dayspol ~ newtreat2, data = df.15, na.rm = TRUE,
                             FUN = c(mean, sd, length2))
summary.dayspol$se <- summary.dayspol$dayspol.sd/(sqrt(summary.dayspol$dayspol.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
dayspol.2way <- aov(dayspol ~ newtreat2 + Block + newtreat2*Block, data = df.15)
anova(dayspol.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(dayspol.2way)

# To check unexplained variability
cdayspol <- lmer(dayspol ~ (1| newtreat2) + (1|Block) + (1|year), data = df.15)
summary(cdayspol)

# Linear Mixed Effect Model
(rdayspol <- lmer(dayspol ~ newtreat2 + (1 | Block), data = df.15))
anova(rdayspol)
summary(rdayspol)


# --------------
# Post hoc test
# --------------

lsm <- lsmeans(rdayspol, pairwise~newtreat2, adjust="tukey")

# cld is a function used to extract and display information 
ph.dp <- cld(lsm, alpha=.05, Letters=letters) 
ph.dp2 <- subset(ph.dp, select=c(newtreat2, .group))
ph.dp2 <- ph.dp2[order(ph.dp2$newtreat2),]
summary.dayspol$poho <- ph.dp2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.dayspol$se.mean.p <- c(summary.dayspol$dayspol.mean + summary.dayspol$se + 3)

# Make the graph
ggplot(summary.dayspol, aes(newtreat2, dayspol.mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = dayspol.mean - se, ymax = dayspol.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Days to Pollen Shed") +
  geom_text(data = summary.dayspol, aes(x = newtreat2, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))



# ------ dayssilk --------------------------------------------------------------
# Kinda Normal distribution, strange shaped density plots
# Treatment/block/interaction significant in ANOVA
# Equal variability explaind by Block and interaction in Random blocking
#    Doing post-hoc based on random blocking with blocking variability removed
# ------------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All groups
ggplot(aes(y = dayssilk, x = newtreat2), data = df.15) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Days to Silk Emergence")

# Each group
ggplot(aes(y = dayssilk, x = newtreat2), data = df.15) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Days to Silk Emergence")


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.15, aes(x=dayssilk)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.15, aes(dayssilk)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.dayssilk <- summaryBy(dayssilk ~ newtreat2, data = df.15, na.rm = TRUE,
                              FUN = c(mean, sd, length2))
summary.dayssilk$se <- summary.dayssilk$dayssilk.sd/(sqrt(summary.dayssilk$dayssilk.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
dayssilk.2way <- aov(dayssilk ~ newtreat2 + Block + newtreat2*Block, data = df.15)
anova(dayssilk.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(dayssilk.2way)

# To check unexplained variability
cdayssilk <- lmer(dayssilk ~ (1| newtreat2) + (1|Block) + (1|year), data = df.15)
summary(cdayssilk)

# Linear Mixed Effect Model
(rdayssilk <- lmer(dayssilk ~ newtreat2 + (1 | Block), data = df.15))
anova(rdayssilk)
summary(rdayssilk)


# -------------
# Post hoc test
# --------------

lsm <- lsmeans(rdayssilk, pairwise~newtreat2, adjust="tukey")

# cld is a function used to extract and display information 
ph.ds <- cld(lsm, alpha=.05, Letters=letters) 
ph.ds2 <- subset(ph.ds, select=c(newtreat2, .group))
ph.ds2 <- ph.ds2[order(ph.ds2$newtreat2),]
summary.dayssilk$poho <- ph.ds2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.dayssilk$se.mean.p <- c(summary.dayssilk$dayssilk.mean + summary.dayssilk$se + 3)

# Make the graph
ggplot(summary.dayssilk, aes(newtreat2, dayssilk.mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = dayssilk.mean - se, ymax = dayssilk.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Days to Silk Emergence") +
  geom_text(data = summary.dayssilk, aes(x = newtreat2, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))





