# ----------------------------------------------------
# Merritt Does things with data!
# Double Haploid Linear mixed models & ANOVA Analysis
# Summer 2014 Data
# 8/30/2016
# ----------------------------------------------------

# Clear Directory
rm(list=ls())

# Set working directory
setwd("~/SDSU Masters/Vivek Data/DH14_Images")

# Load packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'doBy', 'lmtest', 'lsmeans', 'agricolae', 'tidyr')
pack.man(packages)

# Load datframe
df.dh <- read.csv("~/SDSU Masters/Vivek Data/Data Sets/DH14_BorisEdit.csv")

# Rename columns
names(df.dh)[names(df.dh) == "newtreat2"] <- "old.lingen"
names(df.dh)[names(df.dh) == "newtreat3"] <- "lineage.generation"
names(df.dh)[names(df.dh) == "newtreat4"] <- "lineage"

# Remove extraneous rows
df.dh <- df.dh[, -1]

# Convert treatments to factors
df.dh[, 1:7] <- lapply(df.dh[, 1:7], factor)

# Add columns for rate
df.dh$rate.silk <- 1/df.dh$dayssilk
df.dh$rate.pol  <- 1/df.dh$dayspol

median(df.dh$[,8:22])

# ------ No.tassel -------------------------------------------

# --------
# Boxplot
# --------

# All blocks
ggplot(aes(y = No.tassel, x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Number of Tassels") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_No.tassel_boxplot all blocks.png", width = 10.5, height = 6)

# By block
ggplot(aes(y = No.tassel, x = lineage), data = df.dh) + 
  geom_boxplot( aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Number of tassels") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_No.tassel_boxplot each block.png", width = 10.5, height = 6)


# -------------------------
#  Summary Stats and Graphs
# -------------------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}

summary.No.tassel <- summaryBy(No.tassel ~ lineage, data = df.dh, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.No.tassel$se <- summary.No.tassel$No.tassel.sd/(sqrt(summary.No.tassel$No.tassel.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
No.tassel.2way <- aov(No.tassel ~ lineage + Block + lineage*Block, data = df.dh)
anova(No.tassel.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(No.tassel.2way)


# --------------------------
# Linear mixed effect models
# --------------------------
rNo.tassel <- lmer(No.tassel ~ lineage + (1 | Block), data = df.dh)
anova(rNo.tassel)
summary(rNo.tassel)

# To check unexplained variability
cNo.tassel <- lmer(No.tassel ~ (1| lineage) + (1|Block) + (1|year), data = df.dh)
summary(cNo.tassel)

# Post hoc test
lsm <- lsmeans(rNo.tassel, pairwise ~ lineage, adjust="tukey")

# cld is a function used to extract and display information 
phnt <- cld(lsm, alpha=.05, Letters=letters) 
phnt2 <- subset(phnt, select=c(lineage, .group))
phnt2 <- phnt2[order(phnt2$lineage),]
summary.No.tassel$poho <- phnt2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.No.tassel$se.mean.p <- c(summary.No.tassel$No.tassel.mean + summary.No.tassel$se + 0.6)

# Make Graph
ggplot(summary.No.tassel, aes(lineage, No.tassel.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = No.tassel.mean - se, ymax = No.tassel.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Number of Tassels") +
  geom_text(data = summary.No.tassel, aes(x = lineage, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_No.tassel_barchart with post hoc.png", width = 10.5, height = 6)

# ------ Leafwidth.cm. ---------------------------------------------------------

# --------
# Boxplots
#---------

# All blocks
ggplot(aes(y = Leafwidth.cm., x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Leaf Width (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Leafwidth.cm._boxplot all blocks.png", width = 10.5, height = 6)

# Each Block
ggplot(aes(y = Leafwidth.cm., x = lineage), data = df.dh) + 
  geom_boxplot( aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Leaf Width (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Leafwidth.cm._boxplot each block.png", width = 10.5, height = 6)


# -------------
# Summary Stats 
# -------------

summary.Leafwidth.cm. <- summaryBy(Leafwidth.cm. ~ lineage, data = df.dh, na.rm = TRUE,
                                   FUN = c(mean, sd, length2))
summary.Leafwidth.cm.$se <- summary.Leafwidth.cm.$Leafwidth.cm..sd/(sqrt(summary.Leafwidth.cm.$Leafwidth.cm..length))


# --------------------
# Statistical Testing
# --------------------

# ANOVA
Leafwidth.cm..2way <- aov(Leafwidth.cm. ~ lineage + Block + lineage*Block, data = df.dh)
anova(Leafwidth.cm..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Leafwidth.cm..2way)

# To check unexplained variability
cLeafwidth.cm. <- lmer(Leafwidth.cm. ~ (1| lineage) + (1|Block) + (1| lineage:Block) , data = df.dh)
summary(cLeafwidth.cm.)

# Mixed Effect modeling
rLeafwidth.cm. <- lmer(Leafwidth.cm. ~ lineage + (1 | Block), data = df.dh)
anova(rLeafwidth.cm.)
summary(rLeafwidth.cm.)

# Post hoc test
lsm <- lsmeans(rLeafwidth.cm., pairwise~lineage, adjust="tukey")

# cld is a function used to extract and display information 
ph.lw <- cld(lsm, alpha=.05, Letters=letters) 
ph.lw2 <- subset(ph.lw, select=c(lineage, .group))
ph.lw2 <- ph.lw2[order(ph.lw2$lineage),]
summary.Leafwidth.cm.$poho <- ph.lw2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.Leafwidth.cm.$se.mean.p <- c(summary.Leafwidth.cm.$Leafwidth.cm..mean + summary.Leafwidth.cm.$se + 0.5)

ggplot(summary.Leafwidth.cm., aes(lineage, Leafwidth.cm..mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = Leafwidth.cm..mean - se, ymax = Leafwidth.cm..mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Leaf Width (cm)") +
  geom_text(data = summary.Leafwidth.cm., aes(x = lineage, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Leafwidth.cm._barchart with post hoc.png", width = 10.5, height = 6)



# ------ Leaflength.cm. -------------------------------------------------------

# ---------
# Boxplots
# ---------

# All blocks
ggplot(aes(y = Leaflength.cm., x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Leaf Length (cm)")  +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Leaflength.cm._boxplot all blocks.png", width = 10.5, height = 6)

# Each block  
ggplot(aes(y = Leaflength.cm., x = lineage), data = df.dh) + 
  geom_boxplot( aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Leaf Length (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Leaflength.cm._boxplot each block.png", width = 10.5, height = 6)


# --------------
#  Summary Stats
# --------------

summary.Leaflength.cm. <- summaryBy(Leaflength.cm. ~ lineage, data = df.dh, na.rm = TRUE,
                                    FUN = c(mean, sd, length2))
summary.Leaflength.cm.$se <- summary.Leaflength.cm.$Leaflength.cm..sd/(sqrt(summary.Leaflength.cm.$Leaflength.cm..length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
Leaflength.cm..2way <- aov(Leaflength.cm. ~ lineage + Block + lineage*Block, data = df.dh)
anova(Leaflength.cm..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Leaflength.cm..2way)

# To check unexplained variability
cLeaflength.cm. <- lmer(Leaflength.cm. ~ (1| lineage) + (1|Block) + (1|year), data = df.dh)
summary(cLeaflength.cm.)

# Linear Mixed Effect Model
rLeaflength.cm. <- lmer(Leaflength.cm. ~ lineage + (1 | Block), data = df.dh)
anova(rLeaflength.cm.)
summary(rLeaflength.cm.)


# ------------- 
# Post hoc test
# -------------

# Not significant, not doing 

# ----------------------
# Non-Post-Hoc graphing
# ----------------------

ggplot(summary.Leaflength.cm., aes(lineage, Leaflength.cm..mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = Leaflength.cm..mean - se, ymax = Leaflength.cm..mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Leaf Length (cm)") +
   theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Leaflength.cm._barchart with NO post hoc.png", width = 10.5, height = 6)



# ------ Plant.Height.cm. ------------------------------------------------------

# --------
# Boxplots
# --------

# All groups
ggplot(aes(y = Plant.Height.cm., x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Plant Height (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Plant.Height.cm._boxplot all groups.png", width = 10, height = 6)

# Each group
ggplot(aes(y = Plant.Height.cm., x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block, stat="identity", position = "dodge")) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Plant Height (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Plant.Height.cm._boxplot each group.png", width = 10, height = 6)


# --------------
#  Summary Stats
# --------------

summary.Plant.Height.cm. <- summaryBy(Plant.Height.cm. ~ lineage, data = df.dh, na.rm = TRUE,
                                      FUN = c(mean, sd, length2))
summary.Plant.Height.cm.$se <- summary.Plant.Height.cm.$Plant.Height.cm..sd/(sqrt(summary.Plant.Height.cm.$Plant.Height.cm..length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
Plant.Height.cm..2way <- aov(Plant.Height.cm. ~ lineage + Block + lineage*Block, data = df.dh)
anova(Plant.Height.cm..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Plant.Height.cm..2way)

# To check unexplained variability
cPlant.Height.cm. <- lmer(Plant.Height.cm. ~ (1| lineage) + (1|Block) + (1|year), data = df.dh)
summary(cPlant.Height.cm.)

# Linear Mixed Effect Modeling
rPlant.Height.cm. <- lmer(Plant.Height.cm. ~ lineage + (1 | Block), data = df.dh)
anova(rPlant.Height.cm.)
summary(rPlant.Height.cm.)


# -----------------
# No-Post hoc test
# -----------------

# Not doing post hoc because test not significant

# --------
# Graphing
# --------
# Make the graph
ggplot(summary.Plant.Height.cm., aes(lineage, Plant.Height.cm..mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = Plant.Height.cm..mean - se, ymax = Plant.Height.cm..mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Plant Height (cm)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Plant.Height.cm._barchart with post hoc.png", width = 10.5, height = 6)


# ------ ear.position ----------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ear.position, x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Ear Position") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_ear.position_boxplot all blocks.png", width = 10, height = 6)


# Each Block
ggplot(aes(y = ear.position, x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Ear Position") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_ear.position_boxplot each block.png", width = 10, height = 6)


# ---------------
#  Summary Stats
# ---------------

summary.ear.position <- summaryBy(ear.position ~ lineage, data = df.dh, na.rm = TRUE,
                                  FUN = c(mean, sd, length2))
summary.ear.position$se <- summary.ear.position$ear.position.sd/(sqrt(summary.ear.position$ear.position.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ear.position.2way <- aov(ear.position ~ lineage + Block + lineage*Block, data = df.dh)
anova(ear.position.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ear.position.2way)

# To check unexplained variability
cear.position <- lmer(ear.position ~ (1|lineage) + (1|Block) + (1|year), data = df.dh)
summary(cear.position)

# Linear Mixed Effect Model
rear.position <- lmer(ear.position ~ lineage + (1 | Block), data = df.dh)
anova(rear.position)
summary(rear.position)


# -------------
# Post hoc test
# -------------
lsm <- lsmeans(rear.position, pairwise~lineage, adjust="tukey")

# cld is a function used to extract and display information 
ph.ep <- cld(lsm, alpha=.05, Letters=letters) 
ph.ep2 <- subset(ph.ep, select=c(lineage, .group))
ph.ep2 <- ph.ep2[order(ph.ep2$lineage),]
summary.ear.position$poho <- ph.ep2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.ear.position$se.mean.p <- c(summary.ear.position$ear.position.mean + summary.ear.position$se + .3)

# Make the graph
ggplot(summary.ear.position, aes(lineage, ear.position.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ear.position.mean - se, ymax = ear.position.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Ear Position") +
  geom_text(data = summary.ear.position, aes(x = lineage, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_ear.position_barchart with post hoc.png", width = 10.5, height = 6)


# ------ no..nodes ------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = no..nodes, x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Number of Nodes")  +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_no..nodes_boxplot all blocks.png", width = 10, height = 6)


# Each block
ggplot(aes(y = no..nodes, x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Number of Nodes") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_no..nodes_boxplot each block.png", width = 10, height = 6)


# --------------
#  Summary Stats
# --------------

summary.no..nodes <- summaryBy(no..nodes ~ lineage, data = df.dh, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.no..nodes$se <- summary.no..nodes$no..nodes.sd/(sqrt(summary.no..nodes$no..nodes.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
no..nodes.2way <- aov(no..nodes ~ lineage + Block + lineage*Block, data = df.dh)
anova(no..nodes.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(no..nodes.2way)

# To check unexplained variability
cno..nodes <- lmer(no..nodes ~ (1| lineage) + (1|Block) + (1| year), data = df.dh)
summary(cno..nodes)

# Linear Mixed Effect Modeling
rno..nodes <- lmer(no..nodes ~ lineage + (1 | Block), data = df.dh)
anova(rno..nodes)
summary(rno..nodes)


# --------------
# Post hoc test
# --------------

lsm <- lsmeans(rno..nodes, pairwise~lineage, adjust="tukey")

# cld is a function used to extract and display information 
ph.nn <- cld(lsm, alpha=.05, Letters=letters)
ph.nn2 <- subset(ph.nn, select=c(lineage, .group))
ph.nn2 <- ph.nn2[order(ph.nn2$lineage),]
summary.no..nodes$poho <- ph.nn2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.no..nodes$se.mean.p <- c(summary.no..nodes$no..nodes.mean + summary.no..nodes$se + .7)

# Make the graph
ggplot(summary.no..nodes, aes(lineage, no..nodes.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = no..nodes.mean - se, ymax = no..nodes.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Number of Nodes") +
  geom_text(data = summary.no..nodes, aes(x = lineage, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_no..nodes_barchart with post hoc.png", width = 10.5, height = 6)


# ------ No..of.rows ------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = No..of.rows, x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Number of Rows (per ear)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_No..of.rows_boxplot all blocks.png", width = 10, height = 6)


# Each block
ggplot(aes(y = No..of.rows, x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Number of Rows (per ear)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_No..of.rows_boxplot each block.png", width = 10, height = 6)


# --------------
#  Summary Stats
# --------------

summary.No..of.rows <- summaryBy(No..of.rows ~ lineage, data = df.dh, na.rm = TRUE,
                                 FUN = c(mean, sd, length2))
summary.No..of.rows$se <- summary.No..of.rows$No..of.rows.sd/(sqrt(summary.No..of.rows$No..of.rows.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
No..of.rows.2way <- aov(No..of.rows ~ lineage + Block + lineage*Block, data = df.dh)
anova(No..of.rows.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(No..of.rows.2way)

# To check unexplained variability
cNo..of.rows <- lmer(No..of.rows ~ (1| lineage) + (1|Block) + (1| year), data = df.dh)
summary(cNo..of.rows)

# Linear Mixed Effect Model
rNo..of.rows <- lmer(No..of.rows ~ lineage + (1 | Block), data = df.dh)
anova(rNo..of.rows)
summary(rNo..of.rows)


# ----------------
# No Post hoc test
# ----------------

# Not significant

# ---------
# Graphing
# ---------

# Make the graph
ggplot(summary.No..of.rows, aes(lineage, No..of.rows.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = No..of.rows.mean - se, ymax = No..of.rows.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Number of Rows (per ear)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_No..of.rows_barchart with NO post hoc.png", width = 10.5, height = 6)


# ------ Total.kernel.per.ear -------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = Total.kernel.per.ear, x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Total Number of Kernels (per ear)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Total.kernel.per.ear_boxplot all blocks.png", width = 10, height = 6)


# Each block
ggplot(aes(y = Total.kernel.per.ear, x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Total Number of Kernels (per ear)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Total.kernel.per.ear_boxplot each block.png", width = 10, height = 6)


# --------------
#  Summary Stats
# --------------

summary.Total.kernel.per.ear <- summaryBy(Total.kernel.per.ear ~ lineage, data = df.dh, na.rm = TRUE,
                                          FUN = c(mean, sd, length2))
summary.Total.kernel.per.ear$se <- summary.Total.kernel.per.ear$Total.kernel.per.ear.sd/(sqrt(summary.Total.kernel.per.ear$Total.kernel.per.ear.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
Total.kernel.per.ear.2way <- aov(Total.kernel.per.ear ~ lineage + Block + lineage*Block, data = df.dh)
anova(Total.kernel.per.ear.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Total.kernel.per.ear.2way)

# To check unexplained variability
cTotal.kernel.per.ear <- lmer(Total.kernel.per.ear ~ (1| lineage) + (1|Block) + (1| year), data = df.dh)
summary(cTotal.kernel.per.ear)

# Linear Mized Effect Model
rTotal.kernel.per.ear <- lmer(Total.kernel.per.ear ~ lineage + (1 | Block), data = df.dh)
anova(rTotal.kernel.per.ear)
summary(rTotal.kernel.per.ear)


# -------------
# Post hoc test
# -------------

lsm <- lsmeans(rTotal.kernel.per.ear, pairwise~lineage, adjust="tukey")

# cld is a function used to extract and display information 
ph.tkpe <- cld(lsm, alpha=.05, Letters=letters) 
ph.tkpe2 <- subset(ph.tkpe, select=c(lineage, .group))
ph.tkpe2 <- ph.tkpe2[order(ph.tkpe2$lineage),]
summary.Total.kernel.per.ear$poho <- ph.tkpe2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.Total.kernel.per.ear$se.mean.p <- c(summary.Total.kernel.per.ear$Total.kernel.per.ear.mean + summary.Total.kernel.per.ear$se + 15)

# Make the graph
ggplot(summary.Total.kernel.per.ear, aes(lineage, Total.kernel.per.ear.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = Total.kernel.per.ear.mean - se, ymax = Total.kernel.per.ear.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Total Kernels Per Ear") +
  geom_text(data = summary.Total.kernel.per.ear, aes(x = lineage, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Total.kernel.per.ear_barchart with post hoc.png", width = 10.5, height = 6)



# ------ ear.length -----------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ear.length, x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Ear Length (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_ear.length_boxplot all blocks.png", width = 10, height = 6)


# Each block
ggplot(aes(y = ear.length, x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Ear Length (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_ear.length_boxplot each block.png", width = 10, height = 6)


# --------------
#  Summary Stats
# --------------

summary.ear.length <- summaryBy(ear.length ~ lineage, data = df.dh, na.rm = TRUE,
                                FUN = c(mean, sd, length2))
summary.ear.length$se <- summary.ear.length$ear.length.sd/(sqrt(summary.ear.length$ear.length.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ear.length.2way <- aov(ear.length ~ lineage + Block + lineage*Block, data = df.dh)
anova(ear.length.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ear.length.2way)

# To check unexplained variability
cear.length <- lmer(ear.length ~ (1| lineage) + (1|Block) + (1|year), data = df.dh)
summary(cear.length)

# Linear Mixed Effect Model
rear.length <- lmer(ear.length ~ lineage + (1 | Block), data = df.dh)
anova(rear.length)
summary(rear.length)


# ----------------
# NO Post hoc test
# ----------------

# ---------
# Graphing
# ---------

ggplot(summary.ear.length, aes(lineage, ear.length.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ear.length.mean - se, ymax = ear.length.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Ear Length (cm)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_ear.length_barchart with NO post hoc.png", width = 10.5, height = 6)



# ------ Average.kernel.wt ----------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = Average.kernel.wt, x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Average Kernel Weight (g)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Average.kernel.wt_boxplot all blocks.png", width = 10, height = 6)


# Each block
ggplot(aes(y = Average.kernel.wt, x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Average Kernel Weight (g)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Average.kernel.wt_boxplot each block.png", width = 10, height = 6)


# ---------------
#  Summary Stats
# --------------

summary.Average.kernel.wt <- summaryBy(Average.kernel.wt ~ lineage, data = df.dh, na.rm = TRUE,
                                       FUN = c(mean, sd, length2))
summary.Average.kernel.wt$se <- summary.Average.kernel.wt$Average.kernel.wt.sd/(sqrt(summary.Average.kernel.wt$Average.kernel.wt.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
Average.kernel.wt.2way <- aov(Average.kernel.wt ~ lineage + Block + lineage*Block, data = df.dh)
anova(Average.kernel.wt.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(Average.kernel.wt.2way)

# To check unexplained variability
cAverage.kernel.wt <- lmer(Average.kernel.wt ~ (1| lineage) + (1|Block) + (1| year), data = df.dh)
summary(cAverage.kernel.wt)

# Linear Mixed Effect Model
rAverage.kernel.wt <- lmer(Average.kernel.wt ~ lineage + (1 | Block), data = df.dh)
anova(rAverage.kernel.wt)
summary(rAverage.kernel.wt)


# --------------
# Post hoc test
# --------------

lsm <- lsmeans(rAverage.kernel.wt, pairwise~lineage, adjust="tukey")

# cld is a function used to extract and display information 
ph.akw <- cld(lsm, alpha=.05, Letters=letters) 
ph.akw2 <- subset(ph.akw, select=c(lineage, .group))
ph.akw2 <- ph.akw2[order(ph.akw2$lineage),]
summary.Average.kernel.wt$poho <- ph.akw2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.Average.kernel.wt$se.mean.p <- c(summary.Average.kernel.wt$Average.kernel.wt.mean + summary.Average.kernel.wt$se + 0.007)

# Make the graph
ggplot(summary.Average.kernel.wt, aes(lineage, Average.kernel.wt.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = Average.kernel.wt.mean - se, ymax = Average.kernel.wt.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Average Kernel Weight (g)") +
  geom_text(data = summary.Average.kernel.wt, aes(x = lineage, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_Average.kernel.wt_barchart with post hoc.png", width = 10.5, height = 7)



# ------ ear.cirum. ------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ear.cirum., x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Ear Circumference (cm)") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_ear.cirum._boxplot all blocks.png", width = 10, height = 6)


# Each block
ggplot(aes(y = ear.cirum., x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Ear Circumference (cm)")  +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_ear.cirum._boxplot each block.png", width = 10, height = 6)


# ---------------
#  Summary Stats
# ---------------

summary.ear.cirum. <- summaryBy(ear.cirum. ~ lineage, data = df.dh, na.rm = TRUE,
                                FUN = c(mean, sd, length2))
summary.ear.cirum.$se <- summary.ear.cirum.$ear.cirum..sd/(sqrt(summary.ear.cirum.$ear.cirum..length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ear.cirum..2way <- aov(ear.cirum. ~ lineage + Block + lineage*Block, data = df.dh)
anova(ear.cirum..2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ear.cirum..2way)

# To check unexplained variability
cear.cirum. <- lmer(ear.cirum. ~ (1| lineage) + (1|Block) + (1|year), data = df.dh)
summary(cear.cirum.)

# Linear Mixed Effect Model
rear.cirum. <- lmer(ear.cirum. ~ lineage + (1 | Block), data = df.dh)
anova(rear.cirum.)
summary(rear.cirum.)


# --------------
# Post hoc test
# --------------

lsm <- lsmeans(rear.cirum., pairwise ~ lineage, adjust="tukey")

# cld is a function used to extract and display information 
ph.ec <- cld(lsm, alpha=.05, Letters=letters) 
ph.ec2 <- subset(ph.ec, select=c(lineage, .group))
ph.ec2 <- ph.ec2[order(ph.ec2$lineage),]
summary.ear.cirum.$poho <- ph.ec2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.ear.cirum.$se.mean.p <- c(summary.ear.cirum.$ear.cirum..mean + summary.ear.cirum.$se + 0.3)

# Make the graph
ggplot(summary.ear.cirum., aes(lineage, ear.cirum..mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ear.cirum..mean - se, ymax = ear.cirum..mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Ear Circumfernece (cm)") +
  geom_text(data = summary.ear.cirum., aes(x = lineage, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_ear.cirum._barchart with post hoc.png", width = 10, height = 6)


# ------ rate.pol --------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = rate.pol, x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Rate of Pollen Shed") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_rate.pol_boxplot all blocks.png", width = 10, height = 6)


# Each block
ggplot(aes(y = rate.pol, x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage")+
  ylab("Days to Pollen Shed") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_rate.pol_boxplot each block.png", width = 10, height = 6)


# --------------
# Summary Stats
# --------------

summary.rate.pol <- summaryBy(rate.pol ~ lineage, data = df.dh, na.rm = TRUE,
                              FUN = c(mean, sd, length2))
summary.rate.pol$se <- summary.rate.pol$rate.pol.sd/(sqrt(summary.rate.pol$rate.pol.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
rate.pol.2way <- aov(rate.pol ~ lineage + Block + lineage*Block, data = df.dh)
anova(rate.pol.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(rate.pol.2way)

# To check unexplained variability
crate.pol <- lmer(rate.pol ~ (1| lineage) + (1|Block) + (1|year), data = df.dh)
summary(crate.pol)

# Linear Mixed Effect Model
rrate.pol <- lmer(rate.pol ~ lineage + (1 | Block), data = df.dh)
anova(rrate.pol)
summary(rrate.pol)


# -----------------
# NO Post hoc test
# -----------------

# Graph
ggplot(summary.rate.pol, aes(lineage, rate.pol.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = rate.pol.mean - se, ymax = rate.pol.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Rate of Pollen Shed") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_rate.pol_barchart with NO post hoc.png", width = 10.5, height = 6)


# ------ rate.silk --------------------------------------------------------------

# --------
# Boxplots
# --------

# All groups
ggplot(aes(y = rate.silk, x = lineage), data = df.dh) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Days to Silk Emergence") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_rate.silk_boxplot all groups.png", width = 10, height = 6)

# Each group
ggplot(aes(y = rate.silk, x = lineage), data = df.dh) + 
  geom_boxplot(aes(fill = Block)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  xlab("Lineage") +
  ylab("Days to Silk Emergence") +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_rate.silk_boxplot each group.png", width = 10, height = 6)


# --------------
#  Summary Stats
# --------------

summary.rate.silk <- summaryBy(rate.silk ~ lineage, data = df.dh, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.rate.silk$se <- summary.rate.silk$rate.silk.sd/(sqrt(summary.rate.silk$rate.silk.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
rate.silk.2way <- aov(rate.silk ~ lineage + Block + lineage*Block, data = df.dh)
anova(rate.silk.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(rate.silk.2way)

# To check unexplained variability
crate.silk <- lmer(rate.silk ~ (1| lineage) + (1|Block) + (1|year), data = df.dh)
summary(crate.silk)

# Linear Mixed Effect Model
rrate.silk <- lmer(rate.silk ~ lineage + (1 | Block), data = df.dh)
anova(rrate.silk)
summary(rrate.silk)


# -----------------
# NO Post hoc test
# -----------------

# Graphing
ggplot(summary.rate.silk, aes(lineage, rate.silk.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = rate.silk.mean - se, ymax = rate.silk.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage") + 
  ylab("Rate of Silk Emergence") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH14_lin_rate.silk_barchart with NO post hoc.png", width = 10.5, height = 6)





