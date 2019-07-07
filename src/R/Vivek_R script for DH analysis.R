# Merritt tries stuff --------------


rm(list=ls())


setwd("~/SDSU Masters/Vivek Data")


dhrevoutIQR <- read.csv("~/SDSU Masters/Vivek Data/dhdata14.csv")


head(dhrevoutIQR)
str(dhrevoutIQR)
summary(dhrevoutIQR)




dhrevoutIQR=read.csv("C:\\Users\\Mani\\Documents\\bioinfo\\dhrevoutIQR.csv",header=TRUE)
attach(dhrevoutIQR)
head(dhrevoutIQR)
fix(dhrevoutIQR)
par(mfrow=c(2,2))
summary(dhrevoutIQR)

## Numtassel
hist(Numtassel)
qqnorm(model1$residuals)
qqline(model1$residuals)
require(nortest)
ad.test(model1$residuals)
summary(a)
model1= aov(Numtassel~Trt+Block, data=dhrevoutIQR )
summary(model1)
attributes(model1)
model1$fitted.values
model1$residuals
plot(model1)
plot(model1)
library(lattice)
bwplot(Numtassel~Trt,main="Effect of lineage on No.of tassle branch",ylab="No. of tassel brach", xlab="lineage-generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
library(agricolae)
duncan = duncan.test(model1,Trt, main="Numtassel~Trt + block, data=a", console=TRUE)
duncan<- duncan.test(model1,trt=list("Trt"))
duncan

## Leafwidth
hist(Leafwidth)
qqnorm(Leafwidth)
qqline(Leafwidth)
require(nortest)
ad.test(Leafwidth)
model2= aov(Leafwidth~Trt+Block, data=dhrevoutIQR )
summary(model2)
plot(model2)
lmmodel2=lm(Leafwidth~Trt*Block,data=dhrevoutIQR)
summary(lmmodel2)
plot(lmmodel2)
duncan = duncan.test(model2,Trt, main="Leafwidth~Trt + block, data=a", console=TRUE)
duncan<- duncan.test(model2,trt=list("Trt"))
duncan
### Tkpe
hist(Tkpe)
qqnorm(Tkpe)
qqline(Tkpe)
require(nortest)
ad.test(Tkpe)
model3= aov(Tkpe~Trt+Block, data=dhrevoutIQR )
summary(model3)
plot(model3)
duncan = duncan.test(model3,Trt, main="Tkpe~Trt + block, data=a", console=TRUE)
duncan<- duncan.test(model3,trt=list("Trt"))
duncan
bwplot(Tkpe~Trt,main="Effect of lineage on Total kernel per ear",ylab="Total kernel per ear", xlab="lineage-generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

### Height
hist(Height)
qqnorm(Height)
qqline(Height)
require(nortest)
ad.test(Height)

model4= aov(Height~Trt+Block, data=dhrevoutIQR)
summary(model4)
plot(model4)
duncan = duncan.test(model4,Trt, main="Height~Trt + block, data=a", console=TRUE)
duncan<- duncan.test(model4,trt=list("Trt"))
duncan
bwplot(Height~Trt,main="Effect of lineage on Height",ylab="Height(cm)", xlab="lineage-generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

##Earcircum
hist(Earcircum)
qqnorm(Earcircum)
qqline(Earcircum)
require(nortest)
ad.test(Earcircum)
model5= aov(Earcircum~Trt+Block, data=dhrevoutIQR )
summary(model5)
plot(model5)
duncan = duncan.test(model5,Trt, main="Earcircum~Trt + block, data=a", console=TRUE)
duncan<- duncan.test(model5,trt=list("Trt"))
duncan

## Numrows
hist(Numrows)
qqnorm(Numrows)
qqline(Numrows)
require(nortest)
ad.test(Numrows)
model6= aov(Numrows~Trt+Block, data=dhrevoutIQR )
summary(model6)
plot(model6)
duncan = duncan.test(model6,Trt, main="Numrows~Trt + block, data=dhrevoutIQR", console=TRUE)
duncan<- duncan.test(model6,trt=list("Trt"))
duncan

##Dayspol
hist(Dayspol)
qqnorm(Dayspol)
qqline(Dayspol)
require(nortest)
ad.test(Dayspol)
model7= aov(mdayspol~Trt+Block, data=dhrevoutIQR )
summary(model7)
plot(model7)
duncan = duncan.test(model7,Trt, main="mdayspol~Trt + block, data=dhrevoutIQR", console=TRUE)
duncan<- duncan.test(model7,trt=list("Trt"))
duncan
bwplot(Dayspol~Trt,main="Effect of lineage on No. of days to pollen shed",ylab="No. of days to pollen shed", xlab="lineage-generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

##Dayssilk
hist(Dayssilk)
qqnorm(Dayssilk)
qqline(Dayssilk)
require(nortest)
ad.test(Dayssilk)
model8= aov(updaysilk~Trt+Block, data=dhrevoutIQR)
summary(model8)
plot(model8)
library(agricolae)
duncan = duncan.test(model8,Trt, main="updaysilk~Trt + block, data=dhrevoutIQR", console=TRUE)
duncan<- duncan.test(model8,trt=list("Trt"))
duncan
res=HSD.test(model8,"Trt",main="Tassel with diff gen")
res
bwplot(Dayssilk~Trt,main="Effect of lineage on No.of days to silk emergence",ylab="No. of days for silk emergence", xlab="lineage-generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

##Hgw
hist(Hgw)
qqnorm(Hgw)
qqline(Hgw)
require(nortest)
ad.test(Hgw)
model9= aov(Hgw~Trt+Block, data=dhrevoutIQR )
summary(model9)
plot(model9)
library(agricolae)
duncan = duncan.test(model9,Trt, main="Hgw~Trt + block, data=dhrevoutIQR", console=TRUE)
duncan<- duncan.test(model9,trt=list("Trt"))
duncan

##Leaflength
hist(Leaflength)
qqnorm(Leaflength)
qqline(Leaflength)
require(nortest)
ad.test(Leaflength)
model10= aov(Leaflength~Trt+Block, data=dhrevoutIQR)
summary(model10)
plot(model10)
library(agricolae)
duncan = duncan.test(model10,Trt, main="Leaflength~Trt + block, data=a", console=TRUE)
duncan<- duncan.test(model10,trt=list("Trt"))
duncan
bwplot(Leaflength~Trt,main="Effect of lineage on leaflength",ylab="Leaflength (cm)", xlab="lineage-generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
## Numnodes

bwplot(Numnodes~Trt,main="Effect of lineage on No.of nodes",ylab="No. of nodes", xlab="lineage-generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)


##

barplot(tapply(Dayspol,list(Trt),mean),beside=T,col=T)
tapply(Leaflength,list(Block),sd)
model=aov(a$Numtassel~a$Trt*a$Block,data=a,)
summary(model)
summary.lm(model)

### Yield
model201= aov(Yield~Trt+Block, data=dhrevoutIQR )
summary(model201)
plot(model201)
library(agricolae)
duncan = duncan.test(model201,Trt, main="Yield~Trt + block, data=dhrevoutIQR", console=TRUE)
duncan<- duncan.test(model201,trt=list("Trt"))
duncan

###Leafarea
model202= aov(Leafarea~Trt+Block, data=dhrevoutIQR )
summary(model202)
plot(model201)
library(agricolae)
duncan = duncan.test(model202,Trt, main="Leafarea~Trt + block, data=dhrevoutIQR", console=TRUE)
duncan<- duncan.test(model202,trt=list("Trt"))
duncan

### Regression analysis
lmmodel3=lm(Numtassel~Tkpe+Numrows+Earlength+Earcircum+Hgw+Akw+Dayspol+Dayssilk+Numnodes+Earpos+Height+Leafwidth+Leaflength,data=dhrevoutIQR)
summary(lmmodel3)
pairs(realdh)
cor(lmmodel3)
names(dhrevoutIQR)

### correlation in missing values
cor(na.omit(dhrevoutIQR[,-c(1,2,3)]))
cor(na.omit(dhrevoutIQR[,-c(1,2,3,4,6,7,10,12,13,14)]))

###
pairs(dhrevoutIQR,panel=panel.smooth)

####Data analysis for DH summer 2015 dataset
dhdata=read.csv("C:\\Users\\Mani\\Documents\\bioinfo\\dhdata15.csv",header=TRUE)
attach(dhdata)
dhdata
summary(dhdata)
fix(dhdata)

### tkpe
model1051= aov(avgtkpe~trt+block, data=dhdata)
summary(model1051)
par(mfrow=c(2,2))
plot(model1051)
library(agricolae)
duncan = duncan.test(model1051,trt, main="tkpe~trt + block, data=dhdata", console=TRUE)
duncan<- duncan.test(model1051,trt=list("trt"))
duncan

### dayspol
model51= aov(dayspol~trt+block, data=dhdata)
summary(model51)
par(mfrow=c(2,2))
plot(model51)
library(agricolae)
duncan = duncan.test(model51,trt, main="dayspol~trt + block, data=dhdata", console=TRUE)
duncan<- duncan.test(model51,trt=list("trt"))
duncan

### dayspol
model51= aov(dayspol~trt+block, data=dhdata)
summary(model51)
par(mfrow=c(2,2))
plot(model51)
library(agricolae)
duncan = duncan.test(model51,trt, main="dayspol~trt + block, data=dhdata", console=TRUE)
duncan<- duncan.test(model51,trt=list("trt"))
duncan

###dayssilk
model52= aov(dayssilk~trt+block, data=dhdata)
summary(model52)
par(mfrow=c(2,2))
plot(model52)
library(agricolae)
duncan = duncan.test(model52,trt, main="dayssilk~trt + block, data=dhdata", console=TRUE)
duncan<- duncan.test(model52,trt=list("trt"))
duncan

### nodes
model53= aov(nodes~trt+block, data=dhdata)
summary(model53)
par(mfrow=c(2,2))
plot(model53)
library(agricolae)
duncan = duncan.test(model53,trt, main="nodes~trt + block, data=dhdata", console=TRUE)
duncan<- duncan.test(model53,trt=list("trt"))
duncan

### tassel
model54= aov(tassel~trt+block, data=dhdata)
summary(model54)
par(mfrow=c(2,2))
plot(model54)
library(agricolae)
duncan = duncan.test(model54,trt, main="tassel~trt + block, data=dhdata", console=TRUE)
duncan<- duncan.test(model54,trt=list("trt"))
duncan
library(lattice)
bwplot(tassel~trt,main="Effect of lineage on Tassel branch",ylab="No. of tassel branches", xlab="lineage-generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

### height
model55= aov(height~trt+block, data=dhdata)
summary(model55)
par(mfrow=c(2,2))
plot(model55)
library(agricolae)
duncan = duncan.test(model55,trt, main="height~trt + block, data=dhdata", console=TRUE)
duncan<- duncan.test(model55,trt=list("trt"))
duncan

### leaflength
model56= aov(leaflength~trt+block, data=dhdata)
summary(model56)
par(mfrow=c(2,2))
plot(model56)
library(agricolae)
duncan = duncan.test(model56,trt, main="leaflength~trt + block, data=dhdata", console=TRUE)
duncan<- duncan.test(model56,trt=list("trt"))
duncan

##earpos
model57= aov(earpos~trt+block, data=dhdata)
summary(model57)
par(mfrow=c(2,2))
plot(model57)
library(agricolae)
duncan = duncan.test(model57,trt, main="earpos~trt + block, data=dhdata", console=TRUE)
duncan<- duncan.test(model57,trt=list("trt"))
duncan

### DH graph
dhgraph=read.csv("C:\\Users\\Mani\\Documents\\bioinfo\\dhgraph.csv",header=TRUE)
attach(dhgraph)
dhgraph
library(lattice)
fix(dhgraph)
bwplot(tassell~lineage)
bwplot(tassell~lineage,main="Effect of lineage on No.of tassle branch",ylab="tassell", xlab="lineage",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
model1=aov(tassell~lineage)
summary(model1)

### CSSA conference- no of tassel branch
csa014=read.csv("C:\\Users\\Mani\\Documents\\bioinfo\\cssa2014.csv",header=TRUE)
attach(csa014)
head(csa014)
fix(csa014)
library(lattice)
bwplot(Numtassel~Trt,main="Lineage-Generation effect on No. of Tassel Branch",ylab="No. of Tassel Branch", xlab="Lineage-Generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

## for 2015
csa015=read.csv("C:\\Users\\Mani\\Documents\\bioinfo\\cssa2015.csv",header=TRUE)
attach(csa015)
head(csa015)
fix(csa015)
library(lattice)
bwplot(tassel~trt,main="Lineage-Generation effect on No. of Tassel Branch",ylab="No. of Tassel Branch", xlab="lineage-generation",fill="red",col="blue",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)



