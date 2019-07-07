# ----------------------------------------------------
# Merritt Burch
# Merritt.b.burch@sdstate.edu
# Double Haploid Linear mixed models & ANOVA Analysis
# Summer 2016 Data
# 04/13/2017
# ----------------------------------------------------

# Clear Directory
rm(list=ls())

# Set working directory
setwd("C:/Users/Merritt/Box Sync/SDSU Masters V2/Thesis Projects/Doubled Haploid Project/Data Sets/All Years Final Formatted Data")

# Load packages
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c('ggplot2', 'doBy', 'lmtest', 'lsmeans', 'agricolae', 'tidyr')
pack.man(packages)

# Import Data
df.16 <- read.csv("DH 16 - copied from DH16 directory.csv")

# Convert treatments to factors
df.16[, 1:6] <- lapply(df.16[, 1:6], factor)


# --------------------------
# Summary statistics of data
# --------------------------
stat <- lapply(df.16[ , c(7:15)] , 
               function(x) rbind( mean = mean(x[!is.na(x)]),
                                  sd = sd(x[!is.na(x)]),
                                  median = median(x[!is.na(x)]),
                                  minimum = min(x[!is.na(x)]),
                                  maximum = max(x[!is.na(x)]),
                                  n = length(x[!is.na(x)])))

# Limit number of significant figures
summary <- data.frame(stat)
(summary <- signif(summary, digits = 3))

# transpose data
(summary <- t(summary))

# ----------------------- no.tassels --------------------------

# Create special subset where only tassel data is shown
tassel <- df.16[ which(df.16$no.tassels > 2) , ]

# --------
# Boxplot
# --------

# All blocks
ggplot(aes(y = no.tassels, x = lingen), data = tassel, na.rm = TRUE) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Tassel Branches") +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lingen_no.tassels_boxplot.png", width = 10, height = 6)


# --------------
# Density plots
# --------------

# All Blocks
ggplot(df.16, aes(x=no.tassels)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")
ggsave(file = "DH16_lingen_no.tassels_histogram.png", width = 10, height = 6)


# -------------------------
#  Summary Stats and Graphs
# -------------------------

# First option
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.no.tassels <- summaryBy(no.tassels ~ lingen, data = tassel, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.no.tassels$se <- summary.no.tassels$no.tassels.sd/(sqrt(summary.no.tassels$no.tassels.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
no.tassels.2way <- aov(no.tassels ~ lingen, data = df.16)
anova(no.tassels.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(no.tassels.2way)


# TEST NOT SIGNIFICANT


# ------------------
# Post-Hoc graphing
# ------------------

# Make Graph
ggplot(summary.no.tassels, aes(lingen, no.tassels.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = no.tassels.mean - se, ymax = no.tassels.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Number of Tassel Branches") +
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lin.gen_no.tassels_barchart with NO post hoc.png", width = 10, height = 6)



# ------ no.rows ------------------------------------------------------

# --------
# Boxplots
# --------

# All groups
pdf("temp.pdf", width = 10, height = 6)
temp <- ggplot(aes(y = no.rows, x = reorder(lingen)), data = df.16) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Number of Rows") +
  theme(text = element_text(size = 30))
print(temp)
dev.off()
ggsave(file = "DH16_lingen_no.rows_boxplot.png", width = 10, height = 6)



# --------------
# Density plots
# --------------

# All Blocks
ggplot(df.16, aes(x=no.rows)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")
ggsave(file = "DH16_lingen_no.rows_histogram.png", width = 10, height = 6)

# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.no.rows <- summaryBy(no.rows ~ lingen, data = df.16, na.rm = TRUE,
                                      FUN = c(mean, sd, length2))
summary.no.rows$se <- summary.no.rows$no.rows.sd/(sqrt(summary.no.rows$no.rows.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
no.rows.2way <- aov(no.rows ~ lingen, data = df.16)
anova(no.rows.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(no.rows.2way)


# TEST NOT SIGNIFICANT

# -----------------
# Post-Hoc graphing
# -----------------

# Make the graph
ggplot(summary.no.rows, aes(lingen, no.rows.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = no.rows.mean - se, ymax = no.rows.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Number of Rows") +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lin.gen_no.rows_barchart with NO post hoc.png", width = 10, height = 6)


# ------ ear.circum ----------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ear.circum, x = lingen), data = df.16) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Ear Circumference") +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lingen_ear.circum_boxplot.png", width = 10, height = 6)


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.16, aes(x=ear.circum)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")
ggsave(file = "DH16_lingen_ear.circum_histogram.png", width = 10, height = 6)


# ---------------
#  Summary Stats
# ---------------
length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.ear.circum <- summaryBy(ear.circum ~ lingen, data = df.16, na.rm = TRUE,
                                  FUN = c(mean, sd, length2))
summary.ear.circum$se <- summary.ear.circum$ear.circum.sd/(sqrt(summary.ear.circum$ear.circum.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ear.circum.2way <- aov(ear.circum ~ lingen, data = df.16)
anova(ear.circum.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
plot(ear.circum.2way)


# TEST NOT SIGNIFICANT

# -----------------
# Post-Hoc graphing
# -----------------

# Make the graph
ggplot(summary.ear.circum, aes(lingen, ear.circum.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ear.circum.mean - se, ymax = ear.circum.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Ear Circumference") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lin.gen_ear.circum_barchart with NO post hoc.png", width = 10, height = 6)


# ------ ear.length ------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ear.length, x = lingen), data = df.16) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation")+
  ylab("Kernel Length") +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lingen_ear.length_boxplot.png", width = 10, height = 6)


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.16, aes(x=ear.length)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")


# --------------
#  Summary Stats
# --------------

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)}
summary.ear.length <- summaryBy(ear.length ~ lingen, data = df.16, na.rm = TRUE,
                               FUN = c(mean, sd, length2))
summary.ear.length$se <- summary.ear.length$ear.length.sd/(sqrt(summary.ear.length$ear.length.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
ear.length.2way <- aov(ear.length ~ lingen, data = df.16)
anova(ear.length.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ear.length.2way)

# To check unexplained variability
cear.length <- lmer(ear.length ~ (1| lingen), data = df.16)
summary(cear.length)

# Linear Mixed Effect Modeling
rear.length <- lmer(ear.length ~ lingen, data = df.16)
anova(rear.length)
summary(rear.length)


# --------------
# Post hoc test
# --------------
lsm <- lsmeans(rear.length, pairwise~lingen, adjust="tukey")

# cld is a function used to extract and display information 
ph.nn <- cld(lsm, alpha=.05, Letters=letters)
ph.nn2 <- subset(ph.nn, select=c(lingen, .group))
ph.nn2 <- ph.nn2[order(ph.nn2$lingen),]
summary.ear.length$poho <- ph.nn2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.ear.length$se.mean.p <- c(summary.ear.length$ear.length.mean + summary.ear.length$se + .7)

# Make the graph
ggplot(summary.ear.length, aes(lingen, ear.length.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ear.length.mean - se, ymax = ear.length.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Ear Length (cm)") +
  geom_text(data = summary.ear.length, aes(x = lingen, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lin.gen_ear.length_barchart with post hoc.png", width = 10.5, height = 6)


# ------ ratepol --------------------------------------------------------------

# --------
# Boxplots
# --------

# All blocks
ggplot(aes(y = ratepol, x = lingen), data = df.16) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation") +
  ylab("Rate of Pollen Shed") +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lingen_ratepol_boxplot.png", width = 10, height = 6)


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.16, aes(x=ratepol)) + 
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")


# --------------
# Summary Stats
# --------------

summary.ratepol <- summaryBy(ratepol ~ lingen, data = df.16, na.rm = TRUE,
                             FUN = c(mean, sd, length2))
summary.ratepol$se <- summary.ratepol$ratepol.sd/(sqrt(summary.ratepol$ratepol.length))


# -------------------
# Statistical Testing 
# -------------------

# ANOVA
ratepol.2way <- aov(ratepol ~ lingen, data = df.16)
anova(ratepol.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ratepol.2way)

# To check unexplained variability
cratepol <- lmer(ratepol ~ (1| lingen), data = df.16)
summary(cratepol)

# Linear Mixed Effect Model
rratepol <- lmer(ratepol ~ lingen, data = df.16)
anova(rratepol)
summary(rratepol)


# --------------
# Post hoc test
# --------------

lsm <- lsmeans(rratepol, pairwise ~ lingen, adjust="tukey")

# cld is a function used to extract and display information 
ph.dp <- cld(lsm, alpha=.05, Letters=letters) 
ph.dp2 <- subset(ph.dp, select=c(lingen, .group))
ph.dp2 <- ph.dp2[order(ph.dp2$lingen),]
summary.ratepol$poho <- ph.dp2$.group


# ------------------
# Post-Hoc graphing
# ------------------

# Create location where post hoc groupings will be laid
summary.ratepol$se.mean.p <- c(summary.ratepol$ratepol.mean + summary.ratepol$se + 0.0005)

# Make the graph
ggplot(summary.ratepol, aes(lingen, ratepol.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ratepol.mean - se, ymax = ratepol.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Rate of Pollen Shed") +
  geom_text(data = summary.ratepol, aes(x = lingen, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lin.gen_ratepol_barchart with post hoc.png", width = 10.5, height = 6)



# ------ ratesilk --------------------------------------------------------------

# --------
# Boxplots
# --------

# All groups
ggplot(aes(y = ratesilk, x = lingen), data = df.16) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Lineage - Generation") +
  ylab("Days to Silk Emergence") +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lingen_ratesilk_boxplot all groups.png", width = 10, height = 6)


# -------------
# Density plots
# -------------

# All Blocks
ggplot(df.16, aes(x=ratesilk)) + 
  geom_histogram(aes(y=..density..), binwidth=0.0005, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")
ggsave(file = "DH16_rate silk_hist_all groups.png", width = 10.5, height = 6)


# --------------
#  Summary Stats
# --------------

summary.ratesilk <- summaryBy(ratesilk ~ lingen, data = df.16, na.rm = TRUE,
                              FUN = c(mean, sd, length2))
summary.ratesilk$se <- summary.ratesilk$ratesilk.sd/(sqrt(summary.ratesilk$ratesilk.length))


# --------------------
# Statistical Testing 
# --------------------

# ANOVA
ratesilk.2way <- aov(ratesilk ~ lingen, data = df.16)
anova(ratesilk.2way)

# Visualize residual plots of data
par(mfrow=c(2,2))
par(las=1)
par(mar= c(5,7,4,2))
plot(ratesilk.2way)

# To check unexplained variability
cratesilk <- lmer(ratesilk ~ (1| lingen) , data = df.16)
summary(cratesilk)

# Linear Mixed Effect Model
rratesilk <- lmer(ratesilk ~ lingen, data = df.16)
anova(rratesilk)
summary(rratesilk)


# -------------
# Post hoc test
# --------------

lsm <- lsmeans(rratesilk, pairwise ~ lingen, adjust="tukey")

# cld is a function used to extract and display information 
ph.ds <- cld(lsm, alpha=.05, Letters=letters) 
ph.ds2 <- subset(ph.ds, select=c(lingen, .group))
ph.ds2 <- ph.ds2[order(ph.ds2$lingen),]
summary.ratesilk$poho <- ph.ds2$.group


# -----------------
# Post-Hoc graphing
# -----------------

# Create location where post hoc groupings will be laid
summary.ratesilk$se.mean.p <- c(summary.ratesilk$ratesilk.mean + summary.ratesilk$se + 0.0005)

# Make the graph
ggplot(summary.ratesilk, aes(lingen, ratesilk.mean)) + 
  geom_bar(stat="identity", fill = "royalblue2") +
  geom_errorbar(aes(ymin = ratesilk.mean - se, ymax = ratesilk.mean + se),
                width=.2, position=position_dodge(.9)) +
  xlab("Lineage - Generation") + 
  ylab("Rate of Silk Emergence") +
  geom_text(data = summary.ratesilk, aes(x = lingen, y = se.mean.p, label = poho), size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill=FALSE) +
  theme(text = element_text(size = 30))
ggsave(file = "DH16_lin.gen_ratesilk_barchart with post hoc.png", width = 10, height = 6)




