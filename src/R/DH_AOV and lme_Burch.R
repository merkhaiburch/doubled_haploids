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

# Load packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'doBy', 'lsmeans', 'agricolae')
pack.man(packages)


# Remove extraneous rows
df.dh <- df.dh[, c(-1, -2)]
tbl_df(df.dh)

# Reorder this dataframe (it's a mess!)
tmp <- c(1, 2, 3, 5, 4, 20, 19, 6:18)
df.dh <- df.dh[, tmp]

# Convert treatments to factors
df.dh[, 1:7] <- lapply(df.dh[, 1:7], factor)


# ------ No.tassel -------------------------------------------
# Normal distribution
# Only treatment significant
# Post hoc graph created with matching random block and ANOVA
# -----------------------------------------------------------

# --------
# Boxplot
# --------

# All blocks
ggplot(aes(y = No.tassel, x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Tassels")

# By block
ggplot(aes(y = No.tassel, x = newtreat), data = df.dh) + 
  geom_boxplot( aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of tassels")

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

# --------------
# Density plots
# --------------

# All Blocks
ggplot(df.dh, aes(x=No.tassel)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(No.tassel)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# -------------------------
#  Summary Stats and Graphs
# -------------------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.No.tassel <- summaryBy(No.tassel ~ newtreat, data = df.dh, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.No.tassel$se <- summary.No.tassel$No.tassel.sd/(sqrt(summary.No.tassel$No.tassel.length))


# Barchart with error bars 
ggplot(summary.No.tassel, aes(newtreat, No.tassel.mean)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin = No.tassel.mean - se, ymax = No.tassel.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Number of Tassels")
  

# --------------------
# Statistical Testing 
# --------------------

# ANOVA
No.tassel.2way <- aov(No.tassel ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(No.tassel.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(No.tassel.2way)

# Creates Post-hoc groupings using agricolate
#     Orders treatments in numberical order
#     Adds column of post-hoc groupings to summary table for barchart
(tpoho <- HSD.test(No.tassel.2way, "newtreat", group=TRUE))
nt <- data.frame(tpoho$groups)
nt <- nt[order(nt$trt),]
summary.No.tassel$poho <- nt$M


# --------------------------
# Linear mixed effect models
# --------------------------
(rNo.tassel <- lmer(No.tassel ~ newtreat + (1 | Block), data = df.dh))
anova(rNo.tassel)
summary(rNo.tassel)

# To check unexplained variability
cNo.tassel <- lmer(No.tassel ~ (1| newtreat) + (1|Block) + (1|year), data = df.dh)
summary(cNo.tassel)

# Post hoc test
lsm <- lsmeans(rNo.tassel, pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
phnt <- cld(lsm, alpha=.05, Letters=letters) 
phnt2 <- subset(phnt, select=c(newtreat, .group))
phnt2 <- phnt2[order(phnt2$newtreat),]
summary.No.tassel$poho <- phnt2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.No.tassel$se.mean.p <- c(summary.No.tassel$No.tassel.mean + summary.No.tassel$se + 0.6)

# Make Graph
ggplot(summary.No.tassel, aes(newtreat, No.tassel.mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = No.tassel.mean - se, ymax = No.tassel.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Number of Tassels") +
  geom_text(data = summary.No.tassel, aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
   guides(fill=FALSE)


# ------ Leafwidth.cm. ---------------------------------------------------------
# Has two outrageous outliers
# Skewed heavily to the left but does not change
#       with log or square transformation
# Treatment and blocks significant, interaction borderline significant in ANOVA
# Random blocking = significance  in treatments
# doing post hoc
# ------------------------------------------------------------------------------

# --------
# Boxplots
#---------

# All blocks
ggplot(aes(y = Leafwidth.cm., x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Leaf Width (cm)")

# Each Block
ggplot(aes(y = Leafwidth.cm., x = newtreat), data = df.dh) + 
  geom_boxplot( aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Leaf Width (cm)")


# --------------
# Density plots
# --------------

# All Blocks
ggplot(df.dh, aes(x=Leafwidth.cm.)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(Leafwidth.cm.)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# -------------
# Summary Stats 
# -------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.Leafwidth.cm. <- summaryBy(Leafwidth.cm. ~ newtreat, data = df.dh, na.rm = TRUE,
                                   FUN = c(mean, sd, length2))
summary.Leafwidth.cm.$se <- summary.Leafwidth.cm.$Leafwidth.cm..sd/(sqrt(summary.Leafwidth.cm.$Leafwidth.cm..length))


# --------------------
# Statistical Testing
# --------------------

# ANOVA
Leafwidth.cm..2way <- aov(Leafwidth.cm. ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(Leafwidth.cm..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Leafwidth.cm..2way)

# To check unexplained variability
cLeafwidth.cm. <- lmer(Leafwidth.cm. ~ (1| newtreat) + (1|Block) + (1| newtreat:Block) , data = df.dh)
summary(cLeafwidth.cm.)

# Mixed Effect modeling
(rLeafwidth.cm. <- lmer(Leafwidth.cm. ~ newtreat + (1 | Block), data = df.dh))
anova(rLeafwidth.cm.)
summary(rLeafwidth.cm.)

# Post hoc test
lsm <- lsmeans(rLeafwidth.cm., pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.lw <- cld(lsm, alpha=.05, Letters=letters) 
ph.lw2 <- subset(ph.lw, select=c(newtreat, .group))
ph.lw2 <- ph.lw2[order(ph.lw2$newtreat),]
summary.Leafwidth.cm.$poho <- ph.lw2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.Leafwidth.cm.$se.mean.p <- c(summary.Leafwidth.cm.$Leafwidth.cm..mean + summary.Leafwidth.cm.$se + 0.5)

# Make the graph
ggplot(summary.Leafwidth.cm., aes(newtreat, Leafwidth.cm..mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = Leafwidth.cm..mean - se, ymax = Leafwidth.cm..mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Leaf Width (cm)") +
  geom_text(data = summary.Leafwidth.cm., aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


# ------ Leaflength.cm. -------------------------------------------------------
# Skewed slightly to the left
# Treatment, blocking and interaction significant in ANOVA
# High variability explaind by Block and interaction in Random blocking
#     Doing post-hoc based on random blocking with those variabilities removed
# -----------------------------------------------------------------------------
  
# ---------
# Boxplots
# ---------

# All blocks
ggplot(aes(y = Leaflength.cm., x = newtreat), data = df.dh) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Lineage - Generation")+
    ylab("Leaf Length (cm)")

# Each block  
ggplot(aes(y = Leaflength.cm., x = newtreat), data = df.dh) + 
    geom_boxplot( aes(fill = Block, stat="identity", position = "dodge")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Lineage - Generation")+
    ylab("Leaf Length (cm)")
  

# ------------- 
# Density plots
# -------------

# All Blocks
ggplot(df.dh, aes(x=Leaflength.cm.)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(Leaflength.cm.)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

# Summary Stats
length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)}
summary.Leaflength.cm. <- summaryBy(Leaflength.cm. ~ newtreat, data = df.dh, na.rm = TRUE,
                                      FUN = c(mean, sd, length2))
summary.Leaflength.cm.$se <- summary.Leaflength.cm.$Leaflength.cm..sd/(sqrt(summary.Leaflength.cm.$Leaflength.cm..length))
  

# --------------------
# Statistical Testing 
# --------------------
  
# ANOVA
Leaflength.cm..2way <- aov(Leaflength.cm. ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(Leaflength.cm..2way)
  
# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Leaflength.cm..2way)
  
# To check unexplained variability
cLeaflength.cm. <- lmer(Leaflength.cm. ~ (1| newtreat) + (1|Block) + (1|year), data = df.dh)
summary(cLeaflength.cm.)
  
# Linear Mixed Effect Model
(rLeaflength.cm. <- lmer(Leaflength.cm. ~ newtreat + (1 | Block), data = df.dh))
anova(rLeaflength.cm.)
summary(rLeaflength.cm.)
  

# ------------- 
# Post hoc test
# -------------

lsm <- lsmeans(rLeaflength.cm., pairwise~newtreat, adjust="tukey")
  
# cld is a function used to extract and display information 
phll <- cld(lsm, alpha=.05, Letters=letters) 
phll2 <- subset(phll, select=c(newtreat, .group))
phll2 <- phll2[order(phll2$newtreat),]
summary.Leaflength.cm.$poho <- phll2$.group
  
  
# -----------------
# Post-Hoc graphing
# -----------------
  
# Create location where post hoc groupings will be laid
summary.Leaflength.cm.$se.mean.p <- c(summary.Leaflength.cm.$Leaflength.cm..mean + summary.Leaflength.cm.$se + 5)
  
  # Make the graph
ggplot(summary.Leaflength.cm., aes(newtreat, Leaflength.cm..mean)) + 
    geom_bar(stat="identity",  aes(fill = poho)) +
    geom_errorbar(aes(ymin = Leaflength.cm..mean - se, ymax = Leaflength.cm..mean + se),
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("Lineage - Generation") + 
    ylab("Leaf Length (cm)") + 
    geom_text(data = summary.Leaflength.cm., aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  


# ------ Plant.Height.cm. ------------------------------------------------------
# Skewed slightly to the left b/c two small plant heights
#     Changed one value that said it was only 17.0, probably 
# Treatment, blocking and interaction significant in ANOVA
# High variability explaind by Block and interaction in Random blocking
#     Doing post-hoc based on random blocking with blocking variability removed
# ------------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All groups
ggplot(aes(y = Plant.Height.cm., x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Plant Height (cm)")

# Each group
ggplot(aes(y = Plant.Height.cm., x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Plant Height (cm)")


# --------------
# Density plots
# --------------

# All Blocks
ggplot(df.dh, aes(x=Plant.Height.cm.)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(Plant.Height.cm.)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.Plant.Height.cm. <- summaryBy(Plant.Height.cm. ~ newtreat, data = df.dh, na.rm = TRUE,
                                      FUN = c(mean, sd, length2))
summary.Plant.Height.cm.$se <- summary.Plant.Height.cm.$Plant.Height.cm..sd/(sqrt(summary.Plant.Height.cm.$Plant.Height.cm..length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
Plant.Height.cm..2way <- aov(Plant.Height.cm. ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(Plant.Height.cm..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Plant.Height.cm..2way)

# To check unexplained variability
cPlant.Height.cm. <- lmer(Plant.Height.cm. ~ (1| newtreat) + (1|Block) + (1|year), data = df.dh)
summary(cPlant.Height.cm.)

# Linear Mixed Effect Modeling
(rPlant.Height.cm. <- lmer(Plant.Height.cm. ~ newtreat + (1 | Block), data = df.dh))
anova(rPlant.Height.cm.)
summary(rPlant.Height.cm.)


# -------------
# Post hoc test
# -------------

lsm <- lsmeans(rPlant.Height.cm., pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.ph <- cld(lsm, alpha=.05, Letters=letters) 
ph.ph2 <- subset(ph.ph, select=c(newtreat, .group))
ph.ph2 <- ph.ph2[order(ph.ph2$newtreat),]
summary.Plant.Height.cm.$poho <- ph.ph2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.Plant.Height.cm.$se.mean.p <- c(summary.Plant.Height.cm.$Plant.Height.cm..mean + summary.Plant.Height.cm.$se + 7)

# Make the graph
ggplot(summary.Plant.Height.cm., aes(newtreat, Plant.Height.cm..mean)) + 
  geom_bar(stat="identity", aes(fill = poho)) +
  geom_errorbar(aes(ymin = Plant.Height.cm..mean - se, ymax = Plant.Height.cm..mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Plant Height (cm)") +
  geom_text(data = summary.Plant.Height.cm., aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


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
ggplot(aes(y = ear.position, x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Ear Position")

# Each Block
ggplot(aes(y = ear.position, x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Ear Position")


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.dh, aes(x=ear.position)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(ear.position)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# ---------------
#  Summary Stats
# ---------------
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.ear.position <- summaryBy(ear.position ~ newtreat, data = df.dh, na.rm = TRUE,
                                  FUN = c(mean, sd, length2))
summary.ear.position$se <- summary.ear.position$ear.position.sd/(sqrt(summary.ear.position$ear.position.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ear.position.2way <- aov(ear.position ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(ear.position.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ear.position.2way)

# To check unexplained variability
cear.position <- lmer(ear.position ~ (1|newtreat) + (1|Block) + (1|year), data = df.dh)
summary(cear.position)

# Linear Mixed Effect Model
(rear.position <- lmer(ear.position ~ newtreat + (1 | Block), data = df.dh))
anova(rear.position)
summary(rear.position)


# -------------
# Post hoc test
# -------------
lsm <- lsmeans(rear.position, pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.ep <- cld(lsm, alpha=.05, Letters=letters) 
ph.ep2 <- subset(ph.ep, select=c(newtreat, .group))
ph.ep2 <- ph.ep2[order(ph.ep2$newtreat),]
summary.ear.position$poho <- ph.ep2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.ear.position$se.mean.p <- c(summary.ear.position$ear.position.mean + summary.ear.position$se + .3)

# Make the graph
ggplot(summary.ear.position, aes(newtreat, ear.position.mean)) + 
  geom_bar(stat="identity", aes(fill = poho)) +
  geom_errorbar(aes(ymin = ear.position.mean - se, ymax = ear.position.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Ear Position") +
  geom_text(data = summary.ear.position, aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))


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
ggplot(aes(y = no..nodes, x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Nodes")

# Each block
ggplot(aes(y = no..nodes, x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Nodes")


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.dh, aes(x=no..nodes)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(no..nodes)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.no..nodes <- summaryBy(no..nodes ~ newtreat, data = df.dh, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.no..nodes$se <- summary.no..nodes$no..nodes.sd/(sqrt(summary.no..nodes$no..nodes.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
no..nodes.2way <- aov(no..nodes ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(no..nodes.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(no..nodes.2way)

# To check unexplained variability
cno..nodes <- lmer(no..nodes ~ (1| newtreat) + (1|Block) + (1| year), data = df.dh)
summary(cno..nodes)

# Linear Mixed Effect Modeling
(rno..nodes <- lmer(no..nodes ~ newtreat + (1 | Block), data = df.dh))
anova(rno..nodes)
summary(rno..nodes)


# --------------
# Post hoc test
# --------------
lsm <- lsmeans(rno..nodes, pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.nn <- cld(lsm, alpha=.05, Letters=letters)
ph.nn2 <- subset(ph.nn, select=c(newtreat, .group))
ph.nn2 <- ph.nn2[order(ph.nn2$newtreat),]
summary.no..nodes$poho <- ph.nn2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.no..nodes$se.mean.p <- c(summary.no..nodes$no..nodes.mean + summary.no..nodes$se + .7)

# Make the graph
ggplot(summary.no..nodes, aes(newtreat, no..nodes.mean)) + 
  geom_bar(stat="identity", aes(fill = poho)) +
  geom_errorbar(aes(ymin = no..nodes.mean - se, ymax = no..nodes.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Number of Nodes") +
  geom_text(data = summary.no..nodes, aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))



# ------ No..of.rows ------------------------------------------------------
# Trimodal distribution b/c so few intermediate values 
# Nothing significant in ANOVA
# Little variability explaind by treat/Block/interaction in Random blocking
# Not doing post-hoc analysis
# -------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = No..of.rows, x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Rows (per ear)")

# Each block
ggplot(aes(y = No..of.rows, x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Rows (per ear)")


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.dh, aes(x=No..of.rows)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(No..of.rows)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.No..of.rows <- summaryBy(No..of.rows ~ newtreat, data = df.dh, na.rm = TRUE,
                                 FUN = c(mean, sd, length2))
summary.No..of.rows$se <- summary.No..of.rows$No..of.rows.sd/(sqrt(summary.No..of.rows$No..of.rows.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
No..of.rows.2way <- aov(No..of.rows ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(No..of.rows.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(No..of.rows.2way)

# To check unexplained variability
cNo..of.rows <- lmer(No..of.rows ~ (1| newtreat) + (1|Block) + (1| year), data = df.dh)
summary(cNo..of.rows)

# Linear Mixed Effect Model
(rNo..of.rows <- lmer(No..of.rows ~ newtreat + (1 | Block), data = df.dh))
anova(rNo..of.rows)
summary(rNo..of.rows)


# --------------
# Post hoc test
# --------------
lsm <- lsmeans(rNo..of.rows, pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.nor <- cld(lsm, alpha=.05, Letters=letters) 
ph.nor2 <- subset(ph.nor, select=c(newtreat, .group))
ph.nor2 <- ph.nor2[order(ph.nor2$newtreat),]
summary.No..of.rows$poho <- ph.nor2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.No..of.rows$se.mean.p <- c(summary.No..of.rows$No..of.rows.mean + summary.No..of.rows$se + .7)

# Make the graph
ggplot(summary.No..of.rows, aes(newtreat, No..of.rows.mean)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = No..of.rows.mean - se, ymax = No..of.rows.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Number of Rows (per ear)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))



# ------ Total.kernel.per.ear -------------------------------------------------
# Normal distribution
# Treatment/interaction significant in ANOVA
# Some variability explaind by Block and interaction in Random blocking
#    Doing post-hoc based on random blocking with blocking variability removed
# -----------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = Total.kernel.per.ear, x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Total Number of Kernels (per ear)")

# Each block
ggplot(aes(y = Total.kernel.per.ear, x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Total Number of Kernels (per ear)")


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.dh, aes(x=Total.kernel.per.ear)) + 
  geom_histogram(aes(y=..density..), binwidth=25, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(Total.kernel.per.ear)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# Summary Stats
# Have to add b/c have missing data
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.Total.kernel.per.ear <- summaryBy(Total.kernel.per.ear ~ newtreat, data = df.dh, na.rm = TRUE,
                                          FUN = c(mean, sd, length2))
summary.Total.kernel.per.ear$se <- summary.Total.kernel.per.ear$Total.kernel.per.ear.sd/(sqrt(summary.Total.kernel.per.ear$Total.kernel.per.ear.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
Total.kernel.per.ear.2way <- aov(Total.kernel.per.ear ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(Total.kernel.per.ear.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Total.kernel.per.ear.2way)

# To check unexplained variability
cTotal.kernel.per.ear <- lmer(Total.kernel.per.ear ~ (1| newtreat) + (1|Block) + (1| year), data = df.dh)
summary(cTotal.kernel.per.ear)

# Linear Mized Effect Model
(rTotal.kernel.per.ear <- lmer(Total.kernel.per.ear ~ newtreat + (1 | Block), data = df.dh))
anova(rTotal.kernel.per.ear)
summary(rTotal.kernel.per.ear)


# -------------
# Post hoc test
# -------------

lsm <- lsmeans(rTotal.kernel.per.ear, pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.tkpe <- cld(lsm, alpha=.05, Letters=letters) 
ph.tkpe2 <- subset(ph.tkpe, select=c(newtreat, .group))
ph.tkpe2 <- ph.tkpe2[order(ph.tkpe2$newtreat),]
summary.Total.kernel.per.ear$poho <- ph.tkpe2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.Total.kernel.per.ear$se.mean.p <- c(summary.Total.kernel.per.ear$Total.kernel.per.ear.mean + summary.Total.kernel.per.ear$se + 15)

# Make the graph
ggplot(summary.Total.kernel.per.ear, aes(newtreat, Total.kernel.per.ear.mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = Total.kernel.per.ear.mean - se, ymax = Total.kernel.per.ear.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Total Number of Kernels (per ear)") +
  geom_text(data = summary.Total.kernel.per.ear, aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))



# ------ ear.length -----------------------------------------------------------
# Kinda Normal distribution
# Treatment significant (interaction almost) in ANOVA
# Some variability explaind by Block and interaction in Random blocking
#    Doing post-hoc based on random blocking with blocking variability removed
# -----------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ear.length, x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Ear Length (cm)")

# Each block
ggplot(aes(y = ear.length, x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Ear Length (cm)")


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.dh, aes(x=ear.length)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(ear.length)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.ear.length <- summaryBy(ear.length ~ newtreat, data = df.dh, na.rm = TRUE,
                                FUN = c(mean, sd, length2))
summary.ear.length$se <- summary.ear.length$ear.length.sd/(sqrt(summary.ear.length$ear.length.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ear.length.2way <- aov(ear.length ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(ear.length.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ear.length.2way)

# To check unexplained variability
cear.length <- lmer(ear.length ~ (1| newtreat) + (1|Block) + (1|year), data = df.dh)
summary(cear.length)

# Linear Mixed Effect Model
(rear.length <- lmer(ear.length ~ newtreat + (1 | Block), data = df.dh))
anova(rear.length)
summary(rear.length)


# --------------
# Post hoc test
# --------------

lsm <- lsmeans(rear.length, pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.el <- cld(lsm, alpha=.05, Letters=letters) 
ph.el2 <- subset(ph.el, select=c(newtreat, .group))
ph.el2 <- ph.el2[order(ph.el2$newtreat),]
summary.ear.length$poho <- ph.el2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.ear.length$se.mean.p <- c(summary.ear.length$ear.length.mean + summary.ear.length$se + 0.7)

# Make the graph
ggplot(summary.ear.length, aes(newtreat, ear.length.mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = ear.length.mean - se, ymax = ear.length.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Ear Length (cm)") +
  geom_text(data = summary.ear.length, aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))



# ------ Average.kernel.wt ----------------------------------------------------
# Kinda Normal distribution
# Treatment/block/interaction significant in ANOVA
# Some variability explaind by Block and interaction in Random blocking
#    Doing post-hoc based on random blocking with blocking variability removed
# Lots of 0's in the data, really skewing results, necessary?
# May be causing the post-hoc groupings to look funny
# -----------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = Average.kernel.wt, x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Average Kernel Weight (g)")

# Each block
ggplot(aes(y = Average.kernel.wt, x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Average Kernel Weight (g)")


# --------------
# Density plots
# ---------------

# All Blocks
ggplot(df.dh, aes(x=Average.kernel.wt)) + 
  geom_histogram(aes(y=..density..), binwidth=0.01, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(Average.kernel.wt)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# ---------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.Average.kernel.wt <- summaryBy(Average.kernel.wt ~ newtreat, data = df.dh, na.rm = TRUE,
                                       FUN = c(mean, sd, length2))
summary.Average.kernel.wt$se <- summary.Average.kernel.wt$Average.kernel.wt.sd/(sqrt(summary.Average.kernel.wt$Average.kernel.wt.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
Average.kernel.wt.2way <- aov(Average.kernel.wt ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(Average.kernel.wt.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Average.kernel.wt.2way)

# To check unexplained variability
cAverage.kernel.wt <- lmer(Average.kernel.wt ~ (1| newtreat) + (1|Block) + (1| year), data = df.dh)
summary(cAverage.kernel.wt)

# Linear Mixed Effect Model
(rAverage.kernel.wt <- lmer(Average.kernel.wt ~ newtreat + (1 | Block), data = df.dh))
anova(rAverage.kernel.wt)
summary(rAverage.kernel.wt)


# --------------
# Post hoc test
# --------------

lsm <- lsmeans(rAverage.kernel.wt, pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.akw <- cld(lsm, alpha=.05, Letters=letters) 
ph.akw2 <- subset(ph.akw, select=c(newtreat, .group))
ph.akw2 <- ph.akw2[order(ph.akw2$newtreat),]
summary.Average.kernel.wt$poho <- ph.akw2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.Average.kernel.wt$se.mean.p <- c(summary.Average.kernel.wt$Average.kernel.wt.mean + summary.Average.kernel.wt$se + 0.007)

# Make the graph
ggplot(summary.Average.kernel.wt, aes(newtreat, Average.kernel.wt.mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = Average.kernel.wt.mean - se, ymax = Average.kernel.wt.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Average Kernel Weight (g)") +
  geom_text(data = summary.Average.kernel.wt, aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))



# ------ ear.cirum. ------------------------------------------------------------
# Kinda Normal distribution
# Treatment/block/interaction significant in ANOVA
# Some variability explaind by Block and interaction in Random blocking
#    Doing post-hoc based on random blocking with blocking variability removed
# -------------------------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ear.cirum., x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Ear Circumference (cm)")

# Each block
ggplot(aes(y = ear.cirum., x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Ear Circumference (cm)")


# --------------
# Density plots
# --------------

# All Blocks
ggplot(df.dh, aes(x=ear.cirum.)) + 
  geom_histogram(aes(y=..density..), binwidth=0.3, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(ear.cirum.)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# ---------------
#  Summary Stats
# ---------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.ear.cirum. <- summaryBy(ear.cirum. ~ newtreat, data = df.dh, na.rm = TRUE,
                                FUN = c(mean, sd, length2))
summary.ear.cirum.$se <- summary.ear.cirum.$ear.cirum..sd/(sqrt(summary.ear.cirum.$ear.cirum..length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ear.cirum..2way <- aov(ear.cirum. ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(ear.cirum..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ear.cirum..2way)

# To check unexplained variability
cear.cirum. <- lmer(ear.cirum. ~ (1| newtreat) + (1|Block) + (1|year), data = df.dh)
summary(cear.cirum.)

# Linear Mixed Effect Model
(rear.cirum. <- lmer(ear.cirum. ~ newtreat + (1 | Block), data = df.dh))
anova(rear.cirum.)
summary(rear.cirum.)


# --------------
# Post hoc test
# --------------

lsm <- lsmeans(rear.cirum., pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.ec <- cld(lsm, alpha=.05, Letters=letters) 
ph.ec2 <- subset(ph.ec, select=c(newtreat, .group))
ph.ec2 <- ph.ec2[order(ph.ec2$newtreat),]
summary.ear.cirum.$poho <- ph.ec2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.ear.cirum.$se.mean.p <- c(summary.ear.cirum.$ear.cirum..mean + summary.ear.cirum.$se + 0.3)

# Make the graph
ggplot(summary.ear.cirum., aes(newtreat, ear.cirum..mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = ear.cirum..mean - se, ymax = ear.cirum..mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Ear Circumference (cm)") +
  geom_text(data = summary.ear.cirum., aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))



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
ggplot(aes(y = dayspol, x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Days to Pollen Shed")

# Each block
ggplot(aes(y = dayspol, x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Days to Pollen Shed")


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.dh, aes(x=dayspol)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(dayspol)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
# Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.dayspol <- summaryBy(dayspol ~ newtreat, data = df.dh, na.rm = TRUE,
                             FUN = c(mean, sd, length2))
summary.dayspol$se <- summary.dayspol$dayspol.sd/(sqrt(summary.dayspol$dayspol.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
dayspol.2way <- aov(dayspol ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(dayspol.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(dayspol.2way)

# To check unexplained variability
cdayspol <- lmer(dayspol ~ (1| newtreat) + (1|Block) + (1|year), data = df.dh)
summary(cdayspol)

# Linear Mixed Effect Model
(rdayspol <- lmer(dayspol ~ newtreat + (1 | Block), data = df.dh))
anova(rdayspol)
summary(rdayspol)


# --------------
# Post hoc test
# --------------

lsm <- lsmeans(rdayspol, pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.dp <- cld(lsm, alpha=.05, Letters=letters) 
ph.dp2 <- subset(ph.dp, select=c(newtreat, .group))
ph.dp2 <- ph.dp2[order(ph.dp2$newtreat),]
summary.dayspol$poho <- ph.dp2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.dayspol$se.mean.p <- c(summary.dayspol$dayspol.mean + summary.dayspol$se + 3)

# Make the graph
ggplot(summary.dayspol, aes(newtreat, dayspol.mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = dayspol.mean - se, ymax = dayspol.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Days to Pollen Shed") +
  geom_text(data = summary.dayspol, aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
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
ggplot(aes(y = dayssilk, x = newtreat), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Days to Silk Emergence")

# Each group
ggplot(aes(y = dayssilk, x = newtreat), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Days to Silk Emergence")


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.dh, aes(x=dayssilk)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

# Each Block
ggplot(df.dh, aes(dayssilk)) +
  geom_density(aes(group=Block, color=Block, fill=Block), alpha=0.4)


# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.dayssilk <- summaryBy(dayssilk ~ newtreat, data = df.dh, na.rm = TRUE,
                              FUN = c(mean, sd, length2))
summary.dayssilk$se <- summary.dayssilk$dayssilk.sd/(sqrt(summary.dayssilk$dayssilk.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
dayssilk.2way <- aov(dayssilk ~ newtreat + Block + newtreat*Block, data = df.dh)
anova(dayssilk.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(dayssilk.2way)

# To check unexplained variability
cdayssilk <- lmer(dayssilk ~ (1| newtreat) + (1|Block) + (1|year), data = df.dh)
summary(cdayssilk)

# Linear Mixed Effect Model
(rdayssilk <- lmer(dayssilk ~ newtreat + (1 | Block), data = df.dh))
anova(rdayssilk)
summary(rdayssilk)


# -------------
# Post hoc test
# --------------

lsm <- lsmeans(rdayssilk, pairwise~newtreat, adjust="tukey")

# cld is a function used to extract and display information 
ph.ds <- cld(lsm, alpha=.05, Letters=letters) 
ph.ds2 <- subset(ph.ds, select=c(newtreat, .group))
ph.ds2 <- ph.ds2[order(ph.ds2$newtreat),]
summary.dayssilk$poho <- ph.ds2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.dayssilk$se.mean.p <- c(summary.dayssilk$dayssilk.mean + summary.dayssilk$se + 3)

# Make the graph
ggplot(summary.dayssilk, aes(newtreat, dayssilk.mean)) + 
  geom_bar(stat="identity" , aes(fill = poho)) +
  geom_errorbar(aes(ymin = dayssilk.mean - se, ymax = dayssilk.mean + se),
                width=.2,
                position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Days to Silk Emergence") +
  geom_text(data = summary.dayssilk, aes(x = newtreat, y = se.mean.p, label = poho), size = 4) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))





