
### This script is the analysis of line transect data from Keo Seima Wildlife Sanctuary for the years 2010, 2011, 2013, 2014, 2016, 2018, and 2020. The data belong to WCS Cambodia. The anlaysis has been done by Matt Nuttall (m.n.nuttall1@stir.ac.uk / mattnuttall00@gmail.com) in collaboration with Olly Griffin from WCS (ogriffin@wcs.org)

# The metadata are stored publicly on the Global Biodiversity Information Facility and can be found here: https://doi.org/10.15468/37thhj. Due to the sensitive nature of the data (locations of endangered species), if you would like access to the raw data please get in touch with Olly Griffin or WCS Cambodia. 

# The analytical approach is based on an assessment of a set of different approaches (see Analytical_aproach_test_DS script). The approach that has been decided upon, in agreement with WCS, is as follows:

## Species:
  # BSD - Black-shanked douc
  # YCG - Southern yellow-cheeked crested gibbon
  # GSL - Germain's silver langur
  # LTM - Long-tailed macaque
  # PTM - Northern pig-tailed macaque
  # STM - Stump-tailed macaque
  # BTG - Banteng
  # GAU - Gaur
  # RMJ - Northern red muntjac
  # PIG - Wild pig
  # GPF - Green peafowl


# Black-shanked douc will analysed on an annual basis, as they have sufficient within-year observations. Previously, Red Muntjac was also analysed annually, but based on the latest data, I have decided to revert to pooled analysis for RMJ.

# ALL other species will have pooled detection functions (all years), but "year" will be tested as a continuous covariate

# All priority strata will be removed from all years, except for BSD, YCG, PTM (see explanation in the "Subset Data" sections for those species).



### Binning ###

# There was evidence of lumping at distance 0m for some species, and so binning of detection data was tested.

# species to be binned: BSD 10,11,13,16, RMJ, GSL, LTM, PIG, PTM

# For many of the above species, binning made virtually no difference to estimates, and for none of the species did the inferences regarding trends change. This is good, as it suggests that although the data are not perfect, they are fairly robust to the lumping issue. Below I will go through each species and explain whether or not we will use the binned results or the original. Just a note, in the Intro to DS book, it suggests that as a general indicator, if your estimates (from binned and unbinned) are within 1 SE of each other then you could argue the binning is not really doing much.

# BSD - Very little difference to estimates, and very little difference to precision (although some).  2013 and 2016 are the largest differences between estimates, and they are well within 1 SE.  For 2011 and 2013, the binned estimates are more precise. I can't see a reason why I can't use the binned estimates for those years, and the original estimtes for the other years. Each year is techincally a separate analysis with totally different data, and therefore we can say that the binning in 2011 and 2013 has improved the estimates for those years.  Obviously we can't do this for the other species where data are pooled across years.  
# Therefore the estimates we will use for BSD are 2010 (orig), 2011 (bin), 2013 (bin), 2014 (orig), 2016 (orig), 2018 (orig), 2020 (orig). I will combine the results table and save it in the results folder. I will keep the original results and the binned results saved in separate csv's fo reference. Now the results to use are BSD_results_final_combined.  

# RMJ - Here the differences in estimates between binned and unbinned are more obvious. Although the differences are still less than 1 SE, I feel like when you look at the plot what you are seeing is some systematic bias that is being accounted for (i.e. the differences in estimates are of a similar size and in a similar direction). Therefore I think there is good cause here to go with the binned results, as the binning is supposed to help when lumping is an issues, which we know it is. Therefore I think we can assume that the unbinned results are biased. 

# GSL - This looks similar to RMJ in that it appears there is some systematic bias that is being accounted for. Interestingly the effect disappears in later years. This makes sense as the lumping issue is not there in later years (the teams were trained specifically to not round to 0), and so in later years the bias is not there.  Because the data are pooled, we cannot pick and choose per year as we did with BSD, and therefore we need to use the binned results as it is dealing with the lumping (we hope!).  It's a bit unfortunate as we lose some precision, but such is life. 

# LTM - Some of the differences in estimates between binned and unbinned are greater than 1 SE, and I think the binned estimates paint a slightly different picture (although I don't think the trend analysis will change much). The binned results suggest a large drop in abundance between 2011 and 2013, and also suggest more of a steady increase from 2016 onwards. I think we need to stick with the binned results for this species. Again, we lose some precision in most years, which is a shame.

# PIG - in 2010 and 2014 the difference in estimates between binned and unbinned are greater than 1 SE, although the rest of the years are very similar. This is a tough one, as it seems a shame to lose precision in all years because only 2 years are biased!  But becuase the rest of the years are basically the same, It actually makes little difference if we go with the binned results, and then we are confident that the results are unbiased. we will lose some precision, but nothing catastrophic. I think 2010 is particularly important as the binned result suggest an even more dramatic loss of abundance between 2010 and 2011. This can't be caused by hunting alone....disease?

# PTM - Here I am less convinced this shows bias. The difference in estimates is smaller than 1 SE, and sometimes the binned estimates are higher and sometimes they are lower. This suggests to me that the differences are more like noise, rather than something systematic. Inferences regarding the trend of this species will not change, and the binning has not provided any increases in precision. Therefore I will stick with the original analysis for PTM.


#### Load libraries and data ####

library(Distance) # develpment version required - install from GitHub
library(dsm)
library(DSsim)
library(mads)
library(rmarkdown)
library(mrds) #  development version required - install from GitHub
library(knitr)
library(tidyverse)
library(patchwork)
library(devtools)
library(remotes)

load("./Output/Data/KSWS_MASTER.Rdata")

#### Prepare data ####


## centre and scale year (stratum) into a continuous variable in allData 
allData$year <- as.factor(allData$stratum)
allData$stratum <- as.vector(scale(allData$stratum, center = T, scale = T)) 
head(allData$stratum)

allData$obs.habitat <- as.factor(allData$obs.habitat)
allData$obs.observer <- as.factor(allData$obs.observer)

str(allData)

### Based on discussions with Olly, we have decided to remove T20 completely. This is becuase it falls outside the core zone, and therefore the survey area (i.e in region.table) is not correct, and thus density estimates are not correct. We are not sure why T20 falls outside the area.  We considered alternatives such as adding a buffer around the transect and increasing the effective survey area to include the buffer, but because we are working with so many species, each with very different home ranges, it was unclear how to do a buffer that would work for all species.  Fortunately, T20 only contributed about 5 observations over the entire study period, and so removing it is not a big issue. We remove the transect below (rather than in the raw data), so that it can easily be undone if future people decide it's not the right thing to do.

allData <- allData[allData$obs.transect != 20, ]
sample.table <- sample.table[sample.table$Sample.Label != 20,]
obs.table <- obs.table[obs.table$Sample.Label != 20,]

### After disucssion with Olly and Hannah O'Kelly, we are not sure where the total survey area of 1807km2 came from. According to GIS, the core area is 1879830403.780 m2 (1880 km2). Hannah thinks that perhaps it is a mistake, based on an older shapefile or similar.  Therefore we have decided to change the area in region.table to 1880. I am doing it here, rather than in the raw data, so that it can be easily undone.
new.area <- 1880000000
full.region.table$Area <- new.area


## check histograms for all species
par(mfrow=c(3,4))
hist(allData$distance[allData$species=="YCG"], main="YCG")
hist(allData$distance[allData$species=="GSL"], main="GSL")
hist(allData$distance[allData$species=="LTM"], main="LTM")
hist(allData$distance[allData$species=="PTM"], main="PTM")
hist(allData$distance[allData$species=="STM"], main="STM")
hist(allData$distance[allData$species=="RED"], main="RED")
hist(allData$distance[allData$species=="BAN"], main="BAN")
hist(allData$distance[allData$species=="GAU"], main="GAU")
hist(allData$distance[allData$species=="PIG"], main="PIG")
hist(allData$distance[allData$species=="GPF"], main="GPF")

hist(allData$distance[allData$species=="YCG"], main="YCG", breaks=40)
hist(allData$distance[allData$species=="GSL"], main="GSL", breaks=40)
hist(allData$distance[allData$species=="LTM"], main="LTM", breaks=40)
hist(allData$distance[allData$species=="PTM"], main="PTM", breaks=40)
hist(allData$distance[allData$species=="STM"], main="STM", breaks=40)
hist(allData$distance[allData$species=="RED"], main="RED", breaks=40)
hist(allData$distance[allData$species=="BAN"], main="BAN", breaks=40)
hist(allData$distance[allData$species=="GAU"], main="GAU", breaks=40)
hist(allData$distance[allData$species=="PIG"], main="PIG", breaks=40)
hist(allData$distance[allData$species=="GPF"], main="GPF", breaks=40)

hist(allData$distance[allData$species=="YCG"], main="YCG", breaks=60)
hist(allData$distance[allData$species=="GSL"], main="GSL", breaks=60)
hist(allData$distance[allData$species=="LTM"], main="LTM", breaks=60)
hist(allData$distance[allData$species=="PTM"], main="PTM", breaks=60)
hist(allData$distance[allData$species=="STM"], main="STM", breaks=60)
hist(allData$distance[allData$species=="RED"], main="RED", breaks=60)
hist(allData$distance[allData$species=="BAN"], main="BAN", breaks=60)
hist(allData$distance[allData$species=="GAU"], main="GAU", breaks=60)
hist(allData$distance[allData$species=="PIG"], main="PIG", breaks=60)
hist(allData$distance[allData$species=="GPF"], main="GPF", breaks=60)

nrow(allData[allData$distance==0,])
nrow(allData[allData$distance<1,])
nrow(allData[allData$distance<5,])
nrow(allData)
314/5052*100

#### Black-shanked douc ####################################
  ## Subset data ####

# subset BSD data
BSD.data <- allData[allData$species=="BSD",] 
head(BSD.data)

# Total number of groups from all years
length(BSD.data$distance) 
 
### 2020 ####

# number of groups in 2020 
length(BSD.data$distance[BSD.data$year==2020]) # 333

## Subset data for 2020 ####

BSD.2020.data <- as.data.frame(BSD.data[BSD.data$year==2020, ]) 
region.table.2020 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2020",])
sample.table.2020 <- as.data.frame(sample.table[sample.table$Region.Label =="2020",])
obs.table.2020 <- as.data.frame(obs.table[obs.table$Region.Label=="2020",])
        
  ## Exploratory plots & linear model ####

par(mfrow=c(1,2))

# distance histograms
hist(BSD.2020.data$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(BSD.2020.data$distance, main=NULL, xlab="Distance (m)", breaks=c(40))
# evidence of evasive movement.  For this species it is plausible, and probabyl likely, that p=1 until a some decent distance from the transect. The speices is noisy and live in large groups, that scatter and flee when disturbed. This histogram suggests p is high until around 40m. 

# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.BSD.20 <- 50 

# Count the number of observations discarded
nrow(BSD.2020.data[BSD.2020.data$distance>trunc.BSD.20,]) 

# total number of groups
length(BSD.2020.data$distance)

34/333*100 # = 10%. 

# There is one observation at >150m, which is messing with the plots below. It is beyond the trunction ditance so I will remove it
BSD.2020.data <- BSD.2020.data[BSD.2020.data$distance < 150,]

# Plot of distance against cluster size
par(mfrow=c(2,2))
plot(BSD.2020.data$size, BSD.2020.data$distance, main="a)", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.BSD20 <- lm(distance~size, data=BSD.2020.data)
lines(BSD.2020.data$size, as.vector(predict(lm.BSD20, BSD.2020.data)))
summary(lm.BSD20)
# No evidence of size bias

# Plot of Observer factor against distance
plot(BSD.2020.data$obs.observer,BSD.2020.data$distance, main="b)", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# Some variation between observers, but not massive.

# Plot of AMPM factor against distance
plot(as.factor(BSD.2020.data$obs.AMPM),BSD.2020.data$distance, main="c)", xlab="Time period",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# no evidence of a difference between morning and afternoon. More outliers in the afternoon

# Plot of habitat as factor against distance
plot(as.factor(BSD.2020.data$obs.habitat), BSD.2020.data$distance, main="d)", xlab="Habitat class",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))                      
# Evidence of difference here - distances are smallest in evergreen and semi-evergreen. Interestingly, distances are signifciantly larger in bamboo than any other habitat. 

# 1 = evergreen forest, 2 = semi-evergreen, 3 = mixed forest, 4 = deciduous forest, 5 = grass land, 6 = bamboo forest.  Larger distances in mixed forests

  ## Fit detection functions ####

    # Uniform ####

# uniform cosine
BSD.df.unif.cos.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                      sample.table=sample.table.2020, obs.table=obs.table.2020, 
                      truncation=trunc.BSD.20, key="unif", adjustment= "cos")
summary(BSD.df.unif.cos.20)
# cosine(1,2,3) selected

# uniform poly
BSD.df.unif.poly.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                      sample.table=sample.table.2020, obs.table=obs.table.2020, 
                      truncation=trunc.BSD.20, key="unif", adjustment= "poly")
summary(BSD.df.unif.poly.20)
# simple polynomial(2) selected

# compare uniform models
bsd.uni.comp.20 <- summarize_ds_models(BSD.df.unif.cos.20,BSD.df.unif.poly.20,
                                    output="plain")
bsd.uni.comp.20[ ,1:4]
bsd.uni.comp.20[ ,5:7]
# All CvM p values are above 0.05. dAIC very similar.

# plot the fits
par(mfrow=c(2,2))
plot(BSD.df.unif.cos.20, main = "BSD.df.unif.cos.20")

covar.fit <- ddf.gof(BSD.df.unif.cos.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

plot(BSD.df.unif.poly.20, main = "BSD.df.unif.poly.20")

covar.fit <- ddf.gof(BSD.df.unif.poly.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

## cos holds p steady for longer then drops off after 30m, whereas poly decreases more steadily. qq plot is better for cos, and I prefer the fit

# BSD.df.unif.cos.20 selected

    # Half-normal ####

# HN cosine
BSD.df.hn.cos.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.BSD.20, key="hn", adjustment="cos")
summary(BSD.df.hn.cos.20)
# cosine(2,3) selected

# HN hermite
BSD.df.hn.herm.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.BSD.20, key="hn", adjustment="herm")
summary(BSD.df.hn.herm.20)
# hermite(4) selected. No convergence for some models

# compare half normal models
bsd.hn.comp.20 <- summarize_ds_models(BSD.df.hn.cos.20,BSD.df.hn.herm.20,
                                    output="plain")
bsd.hn.comp.20[ ,1:4]
bsd.hn.comp.20[ ,5:7]
# CvM p values > 0.05. AIC similar

# Plot the fits
par(mfrow=c(2,2))

# hn cos
plot(BSD.df.hn.cos.20, main = "BSD.df.hn.cos.20")

covar.fit <- ddf.gof(BSD.df.hn.cos.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# hn herm
plot(BSD.df.hn.herm.20, main = "BSD.df.hn.herm.20")

covar.fit <- ddf.gof(BSD.df.hn.herm.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

## very similar fits. herm model decreases slightly eariler, but based on the data I think cos is the better fit

# BSD.df.hn.cos.20 selected

    # Hazard rate ####

# HR cosine
BSD.df.hr.cos.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.BSD.20, key="hr", adjustment="cos")
summary(BSD.df.hr.cos.20)
# no adjustment selected

# HR poly
BSD.df.hr.poly.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.BSD.20, key="hr", adjustment="poly")
summary(BSD.df.hr.poly.20)
# no adjustment selected


# Plot
par(mfrow=c(1,2))
plot(BSD.df.hr.cos.20, main = "BSD.df.hr.cos.20")
covar.fit <- ddf.gof(BSD.df.hr.cos.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM



    # Compare primary models ####

bsd.df.prim.comp.20 <- summarize_ds_models(BSD.df.unif.cos.20, BSD.df.hn.cos.20, BSD.df.hr.cos.20, 
                                           output = "plain")
bsd.df.prim.comp.20[ ,1:5]
bsd.df.prim.comp.20[ ,6:7]

# HR model has the most support (uni & hn dAIC > 2).  The wide shoulder of the HR fit does appear to fit the data, and in fact both the uni and hn models fit a similar shape.  I will just look at the HR and HN together

par(mfrow=c(2,2))

# hn cos
plot(BSD.df.hn.cos.20, main = "BSD.df.hn.cos.20")

covar.fit <- ddf.gof(BSD.df.hn.cos.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# hr cos
plot(BSD.df.hr.cos.20, main = "BSD.df.hr.cos.20")
covar.fit <- ddf.gof(BSD.df.hr.cos.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# To be fair they are very similar. The HR QQ plot looks better. 

# BSD.df.hr.cos.20 selected

    # Models with harsh truncation ####

# set trunc distance 
trunc.harsh <- 40

par(mfrow=c(1,2))

# HR cos
BSD.df.hr.cos.20.harsh <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.harsh, key="hr", adjustment="cos")
summary(BSD.df.hr.cos.20.harsh)
# No adjustment selected

# compare original with harsh truncation
ddf.gof(BSD.df.hr.cos.20.harsh$ddf, main = "BSD.df.hr.cos.20.harsh", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(BSD.df.hr.cos.20$ddf, main = "BSD.df.hr.cos.20",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# no real difference - the original model is taken forward

    # Models with covariates ####

## cluster size
BSD.df.hr.size.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.BSD.20, key="hr", formula = ~size)
summary(BSD.df.hn.size.20)

# plot
par(mfrow=c(1,2))
plot(BSD.df.hn.size.20, main = "BSD.df.hr.size.20")

covar.fit <- ddf.gof(BSD.df.hr.size.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# terrible fit

## Habitat
BSD.df.hr.covar.hab.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.BSD.20, key="hr", formula = ~obs.habitat)
summary(BSD.df.hr.covar.hab.20)

# plot
plot(BSD.df.hr.covar.hab.20, main = "BSD.df.hr.covar.hab.20")

covar.fit <- ddf.gof(BSD.df.hr.covar.hab.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# average fit - better than size

## AM/PM
BSD.df.hr.covar.AMPM.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.BSD.20, key="hr", formula = ~obs.AMPM)
summary(BSD.df.hr.covar.AMPM.20)

# plot
plot(BSD.df.hr.covar.AMPM.20, main = "BSD.df.hr.covar.AMPM.20")

covar.fit <- ddf.gof(BSD.df.hr.covar.AMPM.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Observer
BSD.df.hr.covar.obs.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.BSD.20, key="hr", formula = ~obs.observer)
summary(BSD.df.hn.covar.obs.20)

# plot
plot(BSD.df.hr.covar.obs.20, main = "BSD.df.hr.covar.obs.20")

covar.fit <- ddf.gof(BSD.df.hr.covar.obs.20$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# quite a lot of variation - not super qq plot


    # Compare covariate models ####

bsd.df.cov.comp.20 <- summarize_ds_models(BSD.df.hr.size.20, BSD.df.hr.covar.obs.20, BSD.df.hr.covar.hab.20,
                                       BSD.df.hr.covar.AMPM.20, output = "plain")

bsd.df.cov.comp.20[ ,1:5]
bsd.df.cov.comp.20[ ,6:7]

# model with observer has most support. I'd like to test adding habitat, seeing as it looked interesting in the plots in the previous section

BSD.df.hr.covar.obs.hab.20 <- ds(data=BSD.2020.data, region.table=region.table.2020, 
                    sample.table=sample.table.2020, obs.table=obs.table.2020,
                    truncation=trunc.BSD.20, key="hr", formula = ~obs.observer+obs.habitat)
summary(BSD.df.hr.covar.obs.hab.20)


# add to the summary
bsd.df.cov.comp.20 <- summarize_ds_models(BSD.df.hr.size.20, BSD.df.hr.covar.obs.20, BSD.df.hr.covar.hab.20,
                                       BSD.df.hr.covar.AMPM.20,BSD.df.hr.covar.obs.hab.20, output = "plain")

bsd.df.cov.comp.20[ ,1:5]
bsd.df.cov.comp.20[ ,6:7]

# model with observer and habitat has some support (dAIC=1.2). I'd like to see what difference there is in the estimates

summary(BSD.df.hr.covar.obs.20)
summary(BSD.df.hr.covar.obs.hab.20)
# CV's and SE's for the model with observer and habitat are wildly large

# BSD.df.hr.covar.obs.20 is selected

   # Compare primary and covariate models ####
  
bsd.df.final.compare.20 <- summarize_ds_models(BSD.df.hr.cos.20,BSD.df.hr.covar.obs.20, 
                                               output = "plain")
bsd.df.final.compare.20[ ,1:5]
bsd.df.final.compare.20[ ,6:7]

### the model with covariates is the best and final model

  ## Final results 2020 ####

# extract estimates and put into format to match the final results data.frame for all years so they can just be appended onto the bottom
BSD.grp.results.20 <- BSD.df.hr.covar.obs.20$dht$clusters$N
BSD.grp.results.20$Label <- "Grp"
BSD.grp.results.20$Species <- "BSD"
BSD.grp.results.20$Year <- 2020
BSD.grp.results.20$DetFun <- "annual"
BSD.grp.results.20$Key <- "Hr"
BSD.grp.results.20$Adjust <- NA
BSD.grp.results.20$Covar <- "observer"
BSD.grp.results.20 <- BSD.grp.results.20[,-7]
BSD.grp.results.20 <- BSD.grp.results.20 %>% 
                      select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)


BSD.ind.results.20 <- BSD.df.hr.covar.obs.20$dht$individuals$N
BSD.ind.results.20$Label <- "Ind"
BSD.ind.results.20$Species <- "BSD"
BSD.ind.results.20$Year <- 2020
BSD.ind.results.20$DetFun <- "annual"
BSD.ind.results.20$Key <- "Hr"
BSD.ind.results.20$Adjust <- NA
BSD.ind.results.20$Covar <- "observer"
BSD.ind.results.20 <- BSD.ind.results.20[,-7]
BSD.ind.results.20 <- BSD.ind.results.20 %>% 
                      select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)

# bind together
BSD.results.2020 <- rbind(BSD.grp.results.20,BSD.ind.results.20)
BSD.results.2020 <- BSD.results.2020 %>% rename(N = Estimate)



### Extract density estimates
BSD.grp.density.20 <- BSD.df.hr.covar.obs.20$dht$clusters$D
BSD.ind.density.20 <- BSD.df.hr.covar.obs.20$dht$individuals$D
BSD.density.20 <- rbind(BSD.grp.density.20,BSD.ind.density.20)
BSD.density.20 <- BSD.density.20 %>% rename(D = Estimate)
BSD.density.20 <- BSD.density.20[, -c(1,7)]

# merge abundance and density
BSD.results.2020 <- cbind(BSD.results.2020,BSD.density.20)


### 2018 ####

# Number of groups in 2018
length(BSD.data$distance[BSD.data$year==2018])  # 338

  ## Subset data for 2018 ####

BSD.2018.data <- as.data.frame(BSD.data[BSD.data$year==2018, ]) 
region.table.2018 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2018",])
sample.table.2018 <- as.data.frame(sample.table[sample.table$Region.Label =="2018",])
obs.table.2018 <- as.data.frame(obs.table[obs.table$Region.Label=="2018",])
        
  ## Exploratory plots & linear model ####

par(mfrow=c(1,2))

# Remove an outlier found in the plot to make it clearer
BSD.2018.data.forplot <- BSD.2018.data$distance[BSD.2018.data$distance<100]

# distance histograms
hist(BSD.2018.data.forplot, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(BSD.2018.data.forplot, main=NULL, xlab="Distance (m)", breaks=c(40))
# evidence of evasive movement.  For this species it is plausible, and probabyl likely, that p=1 until a distance of aroun 20m (where the histogram drops off). This is because the species is NOT cryptic - they are noisy and exist in large groups. Therefore it is very unlikely that they will go undetected so close to the line. This hump at around 10m will be evaasive movement.

# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.BSD.18 <- 60 

# Count the number of observations discarded
nrow(BSD.2018.data[BSD.2018.data$distance>trunc.BSD.18,]) 

# total number of groups
length(BSD.2018.data$distance)

14/338*100 # = 4%

# Plot of distance against cluster size
par(mfrow=c(2,2))
plot(BSD.2018.data$size, BSD.2018.data$distance, main="a)", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.BSD18 <- lm(distance~size, data=BSD.2018.data)
lines(BSD.2018.data$size, as.vector(predict(lm.BSD18, BSD.2018.data)))
summary(lm.BSD18)
# No evidence of size bias

# Plot of Observer factor against distance
plot(BSD.2018.data$obs.observer,BSD.2018.data$distance, main="b)", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# Not much variation between observers

# Plot of AMPM factor against distance
plot(BSD.2018.data$obs.AMPM ,BSD.2018.data$distance, main="c)", xlab="Time period",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# no evidence of a difference between morning and afternoon

# Plot of habitat as factor against distance
plot(as.factor(BSD.2018.data$obs.habitat), BSD.2018.data$distance, main="d)", xlab="Habitat class",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))                      
# Evidence of difference here - distances are smallest in evergreen, and get steadily larger through SE, MX, and DDF. Only 1 obs in grassland but large distance.  All makes sense - will need to test habitat in DF. 

# 1 = evergreen forest, 2 = semi-evergreen, 3 = mixed forest, 4 = deciduous forest, 5 = grass land, 6 = bamboo forest.  Larger distances in mixed forests

  ## Fit detection functions ####

    # Uniform ####

# uniform cosine
BSD.df.unif.cos.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                      sample.table=sample.table.2018, obs.table=obs.table.2018, 
                      truncation=trunc.BSD, key="unif", adjustment= "cos")
summary(BSD.df.unif.cos.18)
# cosine(1) selected

# uniform poly
BSD.df.unif.poly.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                      sample.table=sample.table.2018, obs.table=obs.table.2018, 
                      truncation=trunc.BSD, key="unif", adjustment= "poly")
summary(BSD.df.unif.poly.18)
# simple polynomial(2,4) selected

# compare uniform models
bsd.uni.comp.18 <- summarize_ds_models(BSD.df.unif.cos.18,BSD.df.unif.poly.18,
                                    output="plain")
bsd.uni.comp.18[ ,1:4]
bsd.uni.comp.18[ ,5:7]
# All CvM p values are above 0.05. dAIC very similar.

# plot the fits
par(mfrow=c(2,2))
plot(BSD.df.unif.cos.18, main = "BSD.df.unif.cos.18")

covar.fit <- ddf.gof(BSD.df.unif.cos.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

plot(BSD.df.unif.poly.18, main = "BSD.df.unif.poly.18")

covar.fit <- ddf.gof(BSD.df.unif.poly.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

## similar fits but cos drops off more quickly aroudn 25-30m, whereas poly is a steadier decrease.  Poly is more realistic to me

# BSD.df.unif.poly.18 selected

    # Half-normal ####

# HN cosine
BSD.df.hn.cos.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.BSD.18, key="hn", adjustment="cos")
summary(BSD.df.hn.cos.18)
# no adjustment selected

# HN hermite
BSD.df.hn.herm.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.BSD.18, key="hn", adjustment="herm")
summary(BSD.df.hn.herm.18)
# no adjustment selected


# Plot
par(mfrow=c(1,2))
plot(BSD.df.hn.cos.18, main = "BSD.df.hn.cos.18")

covar.fit <- ddf.gof(BSD.df.hn.cos.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

## I don't think the shoulder is wide enough here.  Detection drops off too quickly

    # Hazard rate ####

# HR cosine
BSD.df.hr.cos.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.BSD.18, key="hr", adjustment="cos")
summary(BSD.df.hr.cos.18)
# no adjustment selected

# HR poly
BSD.df.hr.poly.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.BSD.18, key="hr", adjustment="poly")
summary(BSD.df.hr.poly.18)
# no adjustment selected


# Plot
par(mfrow=c(1,2))
plot(BSD.df.hr.cos.18, main = "BSD.df.hr.cos.18")
covar.fit <- ddf.gof(BSD.df.hr.cos.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# This looks better to me. I think it is ecologially realistic that p = 1 up to at least 10m

    # Compare primary models ####

bsd.df.prim.comp.18 <- summarize_ds_models(BSD.df.unif.poly.18, BSD.df.hn.cos.18, BSD.df.hr.cos.18, 
                                           output = "plain")
bsd.df.prim.comp.18[ ,1:5]
bsd.df.prim.comp.18[ ,6:7]

# All have some suppport, although HR is dAIC just over 2.  I want to see HR and UNi together

par(mfrow=c(2,2))

# uni poly
plot(BSD.df.unif.poly.18, main = "BSD.df.unif.poly.18")

covar.fit <- ddf.gof(BSD.df.unif.poly.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# hr cos
plot(BSD.df.hr.cos.18, main = "BSD.df.hr.cos.18")
covar.fit <- ddf.gof(BSD.df.hr.cos.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# The models actually predict similar p at eaceh 10m interval, the main difference being that the HR model says p will be 1 for longer, whereas uni model says it starts to drop off immediately (albeit slowly). Based on my knowledge of this species, I think the HR model is the more realistic one.

# BSD.df.hr.cos.18 selected

    # Models with harsh truncation ####

# set trunc distance 
trunc.harsh <- 40

par(mfrow=c(1,2))

# HN cos
BSD.df.hr.cos.18.harsh <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.harsh, key="hr", adjustment="cos")
summary(BSD.df.hr.cos.18.harsh)
# No adjustment selected

# compare original with harsh truncation
ddf.gof(BSD.df.hr.cos.18.harsh$ddf, main = "BSD.df.hr.cos.18.harsh", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(BSD.df.hr.cos.18$ddf, main = "BSD.df.hr.cos.18",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# The original model is a better fit than the harsh truncation

    # Models with covariates ####

## cluster size
BSD.df.hn.size.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.BSD.18, key="hr", formula = ~size)
summary(BSD.df.hn.size.18)

# plot
par(mfrow=c(1,2))
plot(BSD.df.hn.size.18, main = "BSD.df.hn.size.18")

covar.fit <- ddf.gof(BSD.df.hn.size.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Habitat
BSD.df.hn.covar.hab.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.BSD.18, key="hr", formula = ~obs.habitat)
summary(BSD.df.hn.covar.hab.18)

# plot
plot(BSD.df.hn.covar.hab.18, main = "BSD.df.hn.covar.hab.18")

covar.fit <- ddf.gof(BSD.df.hn.covar.hab.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## AM/PM
BSD.df.hn.covar.AMPM.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.BSD.18, key="hr", formula = ~obs.AMPM)
summary(BSD.df.hn.covar.AMPM.18)

# plot
plot(BSD.df.hn.covar.AMPM.18, main = "BSD.df.hn.covar.AMPM.18")

covar.fit <- ddf.gof(BSD.df.hn.covar.AMPM.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Observer
BSD.df.hn.covar.obs.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.BSD.18, key="hr", formula = ~obs.observer)
summary(BSD.df.hn.covar.obs.18)

# plot
plot(BSD.df.hn.covar.obs.18, main = "BSD.df.hn.covar.obs.18")

covar.fit <- ddf.gof(BSD.df.hn.covar.obs.18$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


    # Compare covariate models ####

bsd.df.cov.comp.18 <- summarize_ds_models(BSD.df.hn.size.18, BSD.df.hn.covar.obs.18, BSD.df.hn.covar.hab.18,
                                       BSD.df.hn.covar.AMPM.18, output = "plain")

bsd.df.cov.comp.18[ ,1:5]
bsd.df.cov.comp.18[ ,6:7]

# model with observer has most support. I'd like to test adding habitat, seeing as it looked interesting in the plots in the previous section

BSD.df.hn.covar.obs.hab.18 <- ds(data=BSD.2018.data, region.table=region.table.2018, 
                    sample.table=sample.table.2018, obs.table=obs.table.2018,
                    truncation=trunc.BSD.18, key="hr", formula = ~obs.observer+obs.habitat)
summary(BSD.df.hn.covar.obs.hab.18)


# add to the summary
bsd.df.cov.comp.18 <- summarize_ds_models(BSD.df.hn.size.18, BSD.df.hn.covar.obs.18, BSD.df.hn.covar.hab.18,
                                       BSD.df.hn.covar.AMPM.18,BSD.df.hn.covar.obs.hab.18, output = "plain")

bsd.df.cov.comp.18[ ,1:5]
bsd.df.cov.comp.18[ ,6:7]

# model with observer and habitat has good support. I think this is plausible - both observer and habitat play an important part in observing primates

   # Compare primary and covariate models ####
  
bsd.df.final.compare.18 <- summarize_ds_models(BSD.df.hr.cos.18,BSD.df.hn.covar.obs.hab.18, output = "plain")
bsd.df.final.compare.18[ ,1:5]
bsd.df.final.compare.18[ ,6:7]

### the model with covariates is the best and final model

  ## Final results 2018 ####

BSD.grp.results.18 <- BSD.df.hn.covar.obs.hab.18$dht$clusters$N
BSD.grp.results.18$Label <- "Grp"
BSD.grp.results.18$Species <- "BSD"
BSD.grp.results.18$Year <- 2018
BSD.grp.results.18$DetFun <- "annual"
BSD.grp.results.18$Key <- "Hn"
BSD.grp.results.18$Adjust <- NA
BSD.grp.results.18$Covar <- "observer + habitat"
BSD.grp.results.18 <- BSD.grp.results.18[,-7]
BSD.grp.results.18 <- BSD.grp.results.18 %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)


BSD.ind.results.18 <- BSD.df.hn.covar.obs.hab.18$dht$individuals$N
BSD.ind.results.18$Label <- "Ind"
BSD.ind.results.18$Species <- "BSD"
BSD.ind.results.18$Year <- 2018
BSD.ind.results.18$DetFun <- "annual"
BSD.ind.results.18$Key <- "Hn"
BSD.ind.results.18$Adjust <- NA
BSD.ind.results.18$Covar <- "observer + habitat"
BSD.ind.results.18 <- BSD.ind.results.18[,-7]
BSD.ind.results.18 <- BSD.ind.results.18 %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)

# bind together
BSD.results.2018 <- rbind(BSD.grp.results.18,BSD.ind.results.18)
BSD.results.2018 <- BSD.results.2018 %>% rename(N = Estimate)


### extract density estimates
BSD.grp.density.18 <- BSD.df.hn.covar.obs.hab.18$dht$clusters$D
BSD.ind.density.18 <- BSD.df.hn.covar.obs.hab.18$dht$individuals$D
BSD.density.18 <- rbind(BSD.grp.density.18,BSD.ind.density.18)
BSD.density.18 <- BSD.density.18[, -c(1,7)]
BSD.density.18 <- BSD.density.18 %>% rename(D = Estimate)

# merge N and D
BSD.results.2018 <- cbind(BSD.results.2018,BSD.density.18)

### 2016 ####

### Previously, effort strata were ignored for this species, as they were for the other species. After discussions with Olly, we decided to check the effect of accounting for the effort strata for BSD in 2011, 2014, and 2016, as the estimates for these years were very high compared to other years, suggesting the strata were having an effect. I tested a couple of different approaches (see "Analytical_approach_test_CDS.R", lines 2190 onwards).  For BSD though, becasue each year is analysed separately, it is fairly simple (as no pooling is required). We think accounting for the strata here and in 2011 & 2014 is the optimal approach.  

# This means that the old data will be loaded and used, as it includes the strata. Using the old data doesn't make any other difference here, because BSD are analysed annually so the 2016 data doesn't change between the old and new data.


# number of groups in 2016
length(BSD.data$distance[BSD.data$year==2016]) #437

## Subset data for 2016 ####

# old data being loaded and used here (see explanation in section above)

# I still need to change the survey area size and remove T20 as I do for all the data at the top of this script.

# load old data (2010-2018) that have strata
load("./Output/Data/Archive/KSWS_MASTER.Rdata")

allData$obs.habitat <- as.factor(allData$obs.habitat)
allData$obs.observer <- as.factor(allData$obs.observer)


## Remove T20
allData <- allData[allData$obs.transect != 20, ]
sample.table <- sample.table[sample.table$Sample.Label != 20,]
obs.table <- obs.table[obs.table$Sample.Label != 20,]


# subset BSD data
BSD.data <- allData[allData$species=="BSD",] 


# subset data for 2016
BSD.2016.data <- as.data.frame(BSD.data[BSD.data$stratum=="2016_H" | BSD.data$stratum=="2016_L", ]) 

region.table.2016 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2016_H" 
                                                     | full.region.table$Region.Label=="2016_L",])
# change survey area
new.area <- 1880000000/2
region.table.2016$Area <- new.area


sample.table.2016 <- as.data.frame(sample.table[sample.table$Region.Label =="2016_H" | 
                                                  sample.table$Region.Label =="2016_L",])

obs.table.2016 <- as.data.frame(obs.table[obs.table$Region.Label=="2016_H" | obs.table$Region.Label=="2016_L",])


        
  ## Exploratory plots & linear model ####

par(mfrow=c(1,2))

# Remove an outlier found in the plot to make it clearer
BSD.2016.data.forplot <- BSD.2016.data$distance[BSD.2016.data$distance<100]

# distance histograms
hist(BSD.2016.data.forplot, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(BSD.2016.data.forplot, main=NULL, xlab="Distance (m)", breaks=c(40))
# evidence of evasive movement (worse than 2018 I think).  For this species it is plausible, and probabyl likely, that p=1 until a distance of aroun 20m (where the histogram drops off). This is because the species is NOT cryptic - they are noisy and exist in large groups. Therefore it is very unlikely that they will go undetected so close to the line. This hump at around 10m will be evaasive movement.

# Update - Olly uncovered issues with lumping around 0. This is due to rounding of angles, probably due to the emphasis during training of not missing animals on the line. This issues is only really present in 2016 and earlier, and this can be seen in the histos above. We will need to try binning the data

# test bin cutpoints
hist(BSD.2016.data.forplot, main=NULL, xlab="Distance (m)", breaks=c(0,10,20,30,40,50,60,80))
hist(BSD.2016.data.forplot, main=NULL, xlab="Distance (m)", breaks=c(0,5,10,15,20,30,40,50,60,80))

# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.BSD.16 <- 50 
trunc.BSD.16.bin <- 60

# Count the number of observations discarded
nrow(BSD.2016.data[BSD.2016.data$distance>trunc.BSD.16,]) 

# total number of groups
length(BSD.2016.data$distance)

27/437*100 # = 6.2%

# Plot of distance against cluster size
par(mfrow=c(2,2))
plot(BSD.2016.data$size, BSD.2016.data$distance, main="a)", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.BSD16 <- lm(distance~size, data=BSD.2016.data)
lines(BSD.2016.data$size, as.vector(predict(lm.BSD16, BSD.2016.data)))
summary(lm.BSD16)
# some evidence of size bias, but not huge

# Plot of Observer factor against distance
plot(BSD.2016.data$obs.observer,BSD.2016.data$distance, main="b)", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# Not much variation between observers

# Plot of AMPM factor against distance
plot(BSD.2016.data$obs.AMPM ,BSD.2016.data$distance, main="c)", xlab="Time period",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# slightly smaller distances in the afternoon. Perhaps because primates tend to be less active in the hot afternoon, and so might not be so easy to spot moving at greater distances

# Plot of habitat as factor against distance
plot(as.factor(BSD.2016.data$obs.habitat), BSD.2016.data$distance, main="d)", xlab="Habitat class",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))                      
# broadly similar pattern to 2018, although distances are larger in mixed forest than in DDF. 

# 1 = evergreen forest, 2 = semi-evergreen, 3 = mixed forest, 4 = deciduous forest, 5 = grass land, 6 = bamboo forest.  Larger distances in mixed forests

  ## Fit detection functions ####

    # Uniform ####

### Unbinned models

# uniform cosine
BSD.df.unif.cos.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                      sample.table=sample.table.2016, obs.table=obs.table.2016, 
                      truncation=trunc.BSD.16, key="unif", adjustment= "cos")
summary(BSD.df.unif.cos.16)
# cosine(1,2) selected

# uniform poly
BSD.df.unif.poly.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                      sample.table=sample.table.2016, obs.table=obs.table.2016, 
                      truncation=trunc.BSD.16, key="unif", adjustment= "poly")
summary(BSD.df.unif.poly.16)
# simple polynomial(2,4) selected

# compare uniform models
bsd.uni.comp.16 <- summarize_ds_models(BSD.df.unif.cos.16,BSD.df.unif.poly.16,
                                    output="plain")
bsd.uni.comp.16[ ,1:4]
bsd.uni.comp.16[ ,5:7]
# All CvM p values are above 0.05. dAIC lower for cos

# plot the fits
par(mfrow=c(2,2))
plot(BSD.df.unif.cos.16, main = "BSD.df.unif.cos.16")

covar.fit <- ddf.gof(BSD.df.unif.cos.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

plot(BSD.df.unif.poly.16, main = "BSD.df.unif.poly.16")

covar.fit <- ddf.gof(BSD.df.unif.poly.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

## similar fits but cos drops off slower, which fits the data and the reality better I think. Cosine qq plot also looks better

# BSD.df.unif.cos.16 selected


### Binned models

## c(0,10,20,30,40,50,60)

# uniform cosine
BSD.df.unif.cos.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                         sample.table=sample.table.2016, obs.table=obs.table.2016, 
                         truncation=trunc.BSD.16.bin, key="unif", adjustment= "cos",
                         cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.unif.cos.16.bin)
# cosine(1,2) selected

# uniform poly
BSD.df.unif.poly.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                          sample.table=sample.table.2016, obs.table=obs.table.2016, 
                          truncation=trunc.BSD.16.bin, key="unif", adjustment= "poly",
                          cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.unif.poly.16.bin)
# simple polynomial(2,4) selected

# compare binned uniform models
bsd.uni.comp.16.bin <- summarize_ds_models(BSD.df.unif.cos.16.bin,BSD.df.unif.poly.16.bin,
                                       output="plain")
bsd.uni.comp.16.bin[ ,1:4]
bsd.uni.comp.16.bin[ ,5:7]

# plot the fits
par(mfrow=c(2,2))

# cos
plot(BSD.df.unif.cos.16.bin, main = "BSD.df.unif.cos.16.bin")


# poly
plot(BSD.df.unif.poly.16.bin, main = "BSD.df.unif.poly.16.bin")



## c(0,5,10,15,20,30,40,50,60,80)

# uniform cosine
BSD.df.unif.cos.16.bin2 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                             sample.table=sample.table.2016, obs.table=obs.table.2016, 
                             truncation=trunc.BSD.16.bin, key="unif", adjustment= "cos",
                             cutpoints = c(0,5,10,15,20,30,40,50,60))
summary(BSD.df.unif.cos.16.bin2)
# cosine(1,2) selected

# uniform poly
BSD.df.unif.poly.16.bin2 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                              sample.table=sample.table.2016, obs.table=obs.table.2016, 
                              truncation=trunc.BSD.16.bin, key="unif", adjustment= "poly",
                              cutpoints = c(0,5,10,15,20,30,40,50,60))
summary(BSD.df.unif.poly.16.bin2)
# simple polynomial(2,4) selected

# compare binned uniform models
bsd.uni.comp.16.bin2 <- summarize_ds_models(BSD.df.unif.cos.16.bin2,BSD.df.unif.poly.16.bin2,
                                           output="plain")
bsd.uni.comp.16.bin2[ ,1:4]
bsd.uni.comp.16.bin2[ ,5:7]

# plot the fits
par(mfrow=c(2,2))

# cos
plot(BSD.df.unif.cos.16.bin2, main = "BSD.df.unif.cos.16.bin2")


# poly
plot(BSD.df.unif.poly.16.bin2, main = "BSD.df.unif.poly.16.bin2")


## c(0,5,10,15,20,25,30,35,40,45,50,55,60)

# uniform cosine
BSD.df.unif.cos.16.bin3 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                              sample.table=sample.table.2016, obs.table=obs.table.2016, 
                              truncation=trunc.BSD.16.bin, key="unif", adjustment= "cos",
                              cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60))
summary(BSD.df.unif.cos.16.bin3)
# cosine(1,2) selected

# plot compared with bin2
par(mfrow=c(1,2))
plot(BSD.df.unif.cos.16.bin2, main = "BSD.df.unif.cos.16.bin2")
plot(BSD.df.unif.cos.16.bin3, main = "BSD.df.unif.cos.16.bin3")

## compare the binned models estimates
bsd.16.nobin.comp <- data.frame(model = "nobin",
                                estimate = BSD.df.unif.cos.16$dht$clusters$N$Estimate,
                                cv = BSD.df.unif.cos.16$dht$clusters$N$cv,
                                se = BSD.df.unif.cos.16$dht$clusters$N$se,
                                lcl = BSD.df.unif.cos.16$dht$clusters$N$lcl,
                                ucl = BSD.df.unif.cos.16$dht$clusters$N$ucl)

bsd.16.bin.comp1 <- data.frame(model = "bin1",
                              estimate = BSD.df.unif.cos.16.bin$dht$clusters$N$Estimate,
                              cv = BSD.df.unif.cos.16.bin$dht$clusters$N$cv,
                              se = BSD.df.unif.cos.16.bin$dht$clusters$N$se,
                              lcl = BSD.df.unif.cos.16.bin$dht$clusters$N$lcl,
                              ucl = BSD.df.unif.cos.16.bin$dht$clusters$N$ucl)

bsd.16.bin.comp2 <- data.frame(model = "bin2",
                               estimate = BSD.df.unif.cos.16.bin2$dht$clusters$N$Estimate,
                               cv = BSD.df.unif.cos.16.bin2$dht$clusters$N$cv,
                               se = BSD.df.unif.cos.16.bin2$dht$clusters$N$se,
                               lcl = BSD.df.unif.cos.16.bin2$dht$clusters$N$lcl,
                               ucl = BSD.df.unif.cos.16.bin2$dht$clusters$N$ucl)

bsd.16.bin.comp3 <- data.frame(model = "bin3",
                               estimate = BSD.df.unif.cos.16.bin3$dht$clusters$N$Estimate,
                               cv = BSD.df.unif.cos.16.bin3$dht$clusters$N$cv,
                               se = BSD.df.unif.cos.16.bin3$dht$clusters$N$se,
                               lcl = BSD.df.unif.cos.16.bin3$dht$clusters$N$lcl,
                               ucl = BSD.df.unif.cos.16.bin3$dht$clusters$N$ucl)

bsd.16.bin.comp <- rbind(bsd.16.nobin.comp,bsd.16.bin.comp1,bsd.16.bin.comp2,bsd.16.bin.comp3)

ggplot(bsd.16.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,15000)

# the estimates are very very similar. The binned estimates are slightly lower than the unbinned estimate. I feel like the bin3 scheme retains more information close to the line and over the rest of the distances, and so I will continue with that

    # Half-normal ####

## unbinned 

# HN cosine
BSD.df.hn.cos.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016,
                    truncation=trunc.BSD.16, key="hn", adjustment="cos")
summary(BSD.df.hn.cos.16)
# cosine(2) selected

# HN hermite
BSD.df.hn.herm.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016,
                    truncation=trunc.BSD.16, key="hn", adjustment="herm")
summary(BSD.df.hn.herm.16)
# hermite(4) selected. Some no covergence

# compare HN models
bsd.hn.comp.16 <- summarize_ds_models(BSD.df.hn.cos.16,BSD.df.hn.herm.16,
                                    output="plain")
bsd.hn.comp.16[ ,1:4]
bsd.hn.comp.16[ ,5:7]
# very little difference in AIC

# Plot the fits
par(mfrow=c(2,2))

# HN cos
plot(BSD.df.hn.cos.16, main = "BSD.df.hn.cos.16")

covar.fit <- ddf.gof(BSD.df.hn.cos.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HN herm
plot(BSD.df.hn.herm.16, main = "BSD.df.hn.herm.16")

covar.fit <- ddf.gof(BSD.df.hn.herm.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# very similar. Hermite has a better qq plot. BSD.df.hn.herm.16 selected



### binned data 

# HN cosine
BSD.df.hn.cos.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                       sample.table=sample.table.2016, obs.table=obs.table.2016,
                       truncation=trunc.BSD.16.bin, key="hn", adjustment="cos",
                       cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60))
summary(BSD.df.hn.cos.16.bin)
# cosine(2)

# HN hermite
BSD.df.hn.herm.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                        sample.table=sample.table.2016, obs.table=obs.table.2016,
                        truncation=trunc.BSD.16.bin, key="hn", adjustment="herm",
                        cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60))
summary(BSD.df.hn.herm.16.bin)
# key function selected

# compare binned uniform models
bsd.hn.comp.16.bin <- summarize_ds_models(BSD.df.hn.cos.16.bin,BSD.df.hn.herm.16.bin,
                                            output="plain")
bsd.hn.comp.16.bin[ ,1:4]
bsd.hn.comp.16.bin[ ,5:7]
# cos selected

# plot bin versus unbinned
par(mfrow=c(1,2))
plot(BSD.df.hn.cos.16.bin, main="BSD.df.hn.cos.16.bin")
plot(BSD.df.hn.herm.16, main = "BSD.df.hn.herm.16")

# compare estimates
bsd.16.hn.cos.bin <- data.frame(model = "binCos",
                                 estimate = BSD.df.hn.cos.16.bin$dht$clusters$N$Estimate,
                                 cv = BSD.df.hn.cos.16.bin$dht$clusters$N$cv,
                                 se = BSD.df.hn.cos.16.bin$dht$clusters$N$se,
                                 lcl = BSD.df.hn.cos.16.bin$dht$clusters$N$lcl,
                                 ucl = BSD.df.hn.cos.16.bin$dht$clusters$N$ucl)

bsd.16.hn.herm.bin <- data.frame(model = "binerm",
                               estimate = BSD.df.hn.herm.16.bin$dht$clusters$N$Estimate,
                               cv = BSD.df.hn.herm.16.bin$dht$clusters$N$cv,
                               se = BSD.df.hn.herm.16.bin$dht$clusters$N$se,
                               lcl = BSD.df.hn.herm.16.bin$dht$clusters$N$lcl,
                               ucl = BSD.df.hn.herm.16.bin$dht$clusters$N$ucl)

bsd.16.hn.unbin <- data.frame(model = "unbin",
                            estimate = BSD.df.hn.herm.16$dht$clusters$N$Estimate,
                            cv = BSD.df.hn.herm.16$dht$clusters$N$cv,
                            se = BSD.df.hn.herm.16$dht$clusters$N$se,
                            lcl = BSD.df.hn.herm.16$dht$clusters$N$lcl,
                            ucl = BSD.df.hn.herm.16$dht$clusters$N$ucl)

bsd.16.hn.bin.comp <- rbind(bsd.16.hn.cos.bin,bsd.16.hn.herm.bin,bsd.16.hn.unbin)

ggplot(bsd.16.hn.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,15000)

# cos is selected for the binned models

    # Hazard rate ####

## unbinned

# HR cosine
BSD.df.hr.cos.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016,
                    truncation=trunc.BSD.16, key="hr", adjustment="cos")
summary(BSD.df.hr.cos.16)
# no adjustment selected. Some errors about not finding starting rates...

# HR poly
BSD.df.hr.poly.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016,
                    truncation=trunc.BSD.16, key="hr", adjustment="poly")
summary(BSD.df.hr.poly.16)
# no adjustment selected


# Plot
par(mfrow=c(1,2))
plot(BSD.df.hr.cos.16, main = "BSD.df.hr.cos.16")
covar.fit <- ddf.gof(BSD.df.hr.cos.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM



## binned 

# cos
BSD.df.hr.cos.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                       sample.table=sample.table.2016, obs.table=obs.table.2016,
                       truncation=trunc.BSD.16.bin, key="hr", adjustment="cos",
                       cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60))
summary(BSD.df.hr.cos.16.bin)
# key only

# poly
BSD.df.hr.poly.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                           sample.table=sample.table.2016, obs.table=obs.table.2016,
                           truncation=trunc.BSD.16.bin, key="hr", adjustment="poly",
                           cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60))
summary(BSD.df.hr.poly.16.bin)
# key only

# plot binned and unbinned
par(mfrow=c(1,2))
plot(BSD.df.hr.cos.16.bin, main="BSD.df.hr.cos.16.bin")
plot(BSD.df.hr.cos.16, main="BSD.df.hr.cos.16")
# binned looks way better to me!


# compare estimates binned and unbinned

bsd.16.hr.bin <- data.frame(model = "bin",
                                 estimate = BSD.df.hr.cos.16.bin$dht$clusters$N$Estimate,
                                 cv = BSD.df.hr.cos.16.bin$dht$clusters$N$cv,
                                 se = BSD.df.hr.cos.16.bin$dht$clusters$N$se,
                                 lcl = BSD.df.hr.cos.16.bin$dht$clusters$N$lcl,
                                 ucl = BSD.df.hr.cos.16.bin$dht$clusters$N$ucl)

bsd.16.hr.unbin <- data.frame(model = "unbin",
                              estimate = BSD.df.hr.cos.16$dht$clusters$N$Estimate,
                              cv = BSD.df.hr.cos.16$dht$clusters$N$cv,
                              se = BSD.df.hr.cos.16$dht$clusters$N$se,
                              lcl = BSD.df.hr.cos.16$dht$clusters$N$lcl,
                              ucl = BSD.df.hr.cos.16$dht$clusters$N$ucl)

bsd.16.hr.bin.comp <- rbind(bsd.16.hr.bin,bsd.16.hr.unbin)

ggplot(bsd.16.hr.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,15000)
# very little difference at all


    # Compare primary models ####

# unbinned models
bsd.df.prim.comp.16 <- summarize_ds_models(BSD.df.unif.cos.16, BSD.df.hn.herm.16, BSD.df.hr.cos.16, 
                                           output = "plain")
bsd.df.prim.comp.16[ ,1:5]
bsd.df.prim.comp.16[ ,6:7]

# HR and Unif have most support.

# plot them together
par(mfrow=c(2,2))

# uni cos
plot(BSD.df.unif.cos.16, main = "BSD.df.unif.cos.16")

covar.fit <- ddf.gof(BSD.df.unif.cos.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# hr cos
plot(BSD.df.hr.cos.16, main = "BSD.df.hr.cos.16")
covar.fit <- ddf.gof(BSD.df.hr.cos.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# The main difference here is that HR suggests that p=1 right up to 20m whereas Uni says it starts to drop off around 15m.  Not huge difference. HR qq plot is better. I think HR is a better fit to the data, and it is plausible

# BSD.df.hr.cos.16 selected



## binned models
bsd.df.prim.comp.16.bin <- summarize_ds_models(BSD.df.unif.cos.16.bin3, BSD.df.hn.cos.16.bin, 
                                               BSD.df.hr.cos.16.bin, 
                                           output = "plain")
bsd.df.prim.comp.16.bin[ ,1:5]
bsd.df.prim.comp.16.bin[ ,6:7]
# hr and hn models both have support (although hr is best)

# plot them together
par(mfrow=c(1,2))
plot(BSD.df.hn.cos.16.bin, main="BSD.df.hn.cos.16.bin")
plot(BSD.df.hr.cos.16.bin, main="BSD.df.hr.cos.16.bin")
# HR model looks the better fit to me, especially close to the line

# BSD.df.hr.cos.16.bin selected as the binned model


    # Models with harsh truncation ####

# set trunc distance 
trunc.harsh <- 40
trunc.harsh.bin <- 50


## unbinned 

par(mfrow=c(1,2))

# HN cos
BSD.df.hr.cos.16.harsh <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016,
                    truncation=trunc.harsh, key="hr", adjustment="cos")
summary(BSD.df.hr.cos.16.harsh)
# No adjustment selected

# compare original with harsh truncation
ddf.gof(BSD.df.hr.cos.16.harsh$ddf, main = "BSD.df.hr.cos.16.harsh", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(BSD.df.hr.cos.16$ddf, main = "BSD.df.hr.cos.16",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# Not much in it but the original model is a better fit than the harsh truncation



## binned

# HR key only
BSD.df.hr.16.bin.harsh <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                             sample.table=sample.table.2016, obs.table=obs.table.2016,
                             truncation=trunc.harsh.bin, key="hr", adjustment="cos",
                             cutpoints = c(0,5,10,15,20,25,30,35,40,45,50))
summary(BSD.df.hr.16.bin.harsh)

# plot fits together
plot(BSD.df.hr.16.bin.harsh, main="BSD.df.hr.16.bin.harsh")
plot(BSD.df.hr.cos.16.bin, main="BSD.df.hr.cos.16.bin")
# pretty similar

# compare estimates
bsd.16.harsh <- data.frame(model = "harsh",
                           estimate = BSD.df.hr.16.bin.harsh$dht$clusters$N$Estimate,
                           cv = BSD.df.hr.16.bin.harsh$dht$clusters$N$cv,
                           se = BSD.df.hr.16.bin.harsh$dht$clusters$N$se,
                           lcl = BSD.df.hr.16.bin.harsh$dht$clusters$N$lcl,
                           ucl = BSD.df.hr.16.bin.harsh$dht$clusters$N$ucl)

bsd.16.hr <- data.frame(model = "orig",
                           estimate = BSD.df.hr.cos.16.bin$dht$clusters$N$Estimate,
                           cv = BSD.df.hr.cos.16.bin$dht$clusters$N$cv,
                           se = BSD.df.hr.cos.16.bin$dht$clusters$N$se,
                           lcl = BSD.df.hr.cos.16.bin$dht$clusters$N$lcl,
                           ucl = BSD.df.hr.cos.16.bin$dht$clusters$N$ucl)

bsd.16.harshcomp <- rbind(bsd.16.harsh,bsd.16.hr)

ggplot(bsd.16.harshcomp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,15000)
# original is ever so slightly more precise, and so is selected


    # Models with covariates ####

### unbinned

## cluster size
BSD.df.hr.size.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016,
                    truncation=trunc.BSD.16, key="hr", formula = ~size)
summary(BSD.df.hr.size.16)

# plot
par(mfrow=c(1,2))
plot(BSD.df.hr.size.16, main = "BSD.df.hr.size.16")

covar.fit <- ddf.gof(BSD.df.hr.size.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Habitat
BSD.df.hr.covar.hab.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016,
                    truncation=trunc.BSD.16, key="hr", formula = ~obs.habitat)
summary(BSD.df.hr.covar.hab.16)

# plot
plot(BSD.df.hr.covar.hab.16, main = "BSD.df.hr.covar.hab.16")

covar.fit <- ddf.gof(BSD.df.hn.covar.hab.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## AM/PM
BSD.df.hr.covar.AMPM.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016,
                    truncation=trunc.BSD.16, key="hr", formula = ~obs.AMPM)
summary(BSD.df.hr.covar.AMPM.16)

# plot
plot(BSD.df.hr.covar.AMPM.16, main = "BSD.df.hr.covar.AMPM.16")

covar.fit <- ddf.gof(BSD.df.hr.covar.AMPM.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Observer
BSD.df.hr.covar.obs.16 <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016,
                    truncation=trunc.BSD.16, key="hr", formula = ~obs.observer)
summary(BSD.df.hr.covar.obs.16)

# plot
plot(BSD.df.hr.covar.obs.16, main = "BSD.df.hr.covar.obs.16")

covar.fit <- ddf.gof(BSD.df.hr.covar.obs.16$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



### binned

## cluster size
BSD.df.hr.size.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                            sample.table=sample.table.2016, obs.table=obs.table.2016,
                            truncation=trunc.BSD.16.bin, key="hr", adjustment="cos",
                            cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60),
                            formula = ~size)
summary(BSD.df.hn.size.16)

plot(BSD.df.hr.size.16.bin)

## habitat
BSD.df.hr.hab.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                            sample.table=sample.table.2016, obs.table=obs.table.2016,
                            truncation=trunc.BSD.16.bin, key="hr", adjustment="cos",
                            cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60),
                            formula = ~obs.habitat)
summary(BSD.df.hr.hab.16.bin)

plot(BSD.df.hr.hab.16.bin)

## AMPM
BSD.df.hr.ampm.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                           sample.table=sample.table.2016, obs.table=obs.table.2016,
                           truncation=trunc.BSD.16.bin, key="hr", adjustment="cos",
                           cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60),
                           formula = ~obs.AMPM)
summary(BSD.df.hr.ampm.16.bin)

plot(BSD.df.hr.ampm.16.bin)

## observer
BSD.df.hr.obs.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                            sample.table=sample.table.2016, obs.table=obs.table.2016,
                            truncation=trunc.BSD.16.bin, key="hr", adjustment="cos",
                            cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60),
                            formula = ~obs.observer)
summary(BSD.df.hr.obs.16.bin)

plot(BSD.df.hr.obs.16.bin)

## habitat and ampm
## AMPM
BSD.df.hr.hab.ampm.16.bin <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                            sample.table=sample.table.2016, obs.table=obs.table.2016,
                            truncation=trunc.BSD.16.bin, key="hr", adjustment="cos",
                            cutpoints = c(0,5,10,15,20,25,30,35,40,45,50,55,60),
                            formula = ~obs.AMPM+obs.habitat)
summary(BSD.df.hr.hab.ampm.16.bin)

plot(BSD.df.hr.hab.ampm.16.bin)



    # Compare covariate models ####

## unbinned 

bsd.df.cov.comp.16 <- summarize_ds_models(BSD.df.hn.size.16, BSD.df.hn.covar.obs.16, BSD.df.hn.covar.hab.16,
                                       BSD.df.hn.covar.AMPM.16, output = "plain")

bsd.df.cov.comp.16[ ,1:5]
bsd.df.cov.comp.16[ ,6:7]

# model with AMPM has most support. 


## binned
bsd.df.cov.comp.16.bin <- summarize_ds_models(BSD.df.hr.size.16.bin,BSD.df.hr.hab.16.bin, 
                                              BSD.df.hr.ampm.16.bin,BSD.df.hr.obs.16.bin,
                                              BSD.df.hr.hab.ampm.16.bin,
                                              output = "plain")

bsd.df.cov.comp.16.bin[ ,1:5]
bsd.df.cov.comp.16.bin[ ,6:7]
# models with habitat and AMPM have support (the other two do not). I retrospecitvely added another model with both ampm and habitat, and this model is the best by far

# check estiamtes are not crazy
summary(BSD.df.hr.hab.ampm.16.bin)
summary(BSD.df.hr.hab.16.bin)
summary(BSD.df.hr.ampm.16.bin)

bsd.16.cov.bin.comp <- data.frame(model = c("hab","ampm", "both"),
                                  estimate = c(BSD.df.hr.hab.16.bin$dht$clusters$N$Estimate,
                                               BSD.df.hr.ampm.16.bin$dht$clusters$N$Estimate,
                                               BSD.df.hr.hab.ampm.16.bin$dht$clusters$N$Estimate),
                                  cv = c(BSD.df.hr.hab.16.bin$dht$clusters$N$cv,
                                         BSD.df.hr.ampm.16.bin$dht$clusters$N$cv,
                                         BSD.df.hr.hab.ampm.16.bin$dht$clusters$N$cv),
                                  se = c(BSD.df.hr.hab.16.bin$dht$clusters$N$se,
                                         BSD.df.hr.ampm.16.bin$dht$clusters$N$se,
                                         BSD.df.hr.hab.ampm.16.bin$dht$clusters$N$se),
                                  lcl = c(BSD.df.hr.hab.16.bin$dht$clusters$N$lcl,
                                          BSD.df.hr.ampm.16.bin$dht$clusters$N$lcl,
                                          BSD.df.hr.hab.ampm.16.bin$dht$clusters$N$lcl),
                                  ucl = c(BSD.df.hr.hab.16.bin$dht$clusters$N$ucl,
                                          BSD.df.hr.ampm.16.bin$dht$clusters$N$ucl,
                                          BSD.df.hr.hab.ampm.16.bin$dht$clusters$N$ucl))
bsd.16.cov.bin.comp
# ampm is the only model with reasonable CV and SE.  The other two are really high

ggplot(bsd.16.cov.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,29000)

# AMPM model selected


   # Compare primary and covariate models ####
  
## unbinned

bsd.df.final.compare.16 <- summarize_ds_models(BSD.df.hr.cos.16,BSD.df.hn.covar.AMPM.16, output = "plain")
bsd.df.final.compare.16[ ,1:5]
bsd.df.final.compare.16[ ,6:7]

### the model with covariates is the best and final unbinned model


## binned

# plot together
par(mfrow=c(1,2))
plot(BSD.df.hr.cos.16.bin, main="BSD.df.hr.cos.16.bin")
plot(BSD.df.hr.ampm.16.bin, main="BSD.df.hr.ampm.16.bin")

bsd.16.final.bin.comp <- data.frame(model = c("ampm", "orig"),
                                  estimate = c(BSD.df.hr.cos.16.bin$dht$clusters$N$Estimate,
                                               BSD.df.hr.ampm.16.bin$dht$clusters$N$Estimate),
                                  cv = c(BSD.df.hr.cos.16.bin$dht$clusters$N$cv,
                                         BSD.df.hr.ampm.16.bin$dht$clusters$N$cv),
                                  se = c(BSD.df.hr.cos.16.bin$dht$clusters$N$se,
                                         BSD.df.hr.ampm.16.bin$dht$clusters$N$se),
                                  lcl = c(BSD.df.hr.cos.16.bin$dht$clusters$N$lcl,
                                          BSD.df.hr.ampm.16.bin$dht$clusters$N$lcl),
                                  ucl = c(BSD.df.hr.cos.16.bin$dht$clusters$N$ucl,
                                          BSD.df.hr.ampm.16.bin$dht$clusters$N$ucl))

ggplot(bsd.16.final.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,13000)

# models are virtually identical.  The ampm is maybe a tiny bit more precise. Techincally I would assume that for a primate habitat would affect the detection process, and therefore with no good argument against it, the covariate model is selected


  ## Final results 2016 ####

## Unbinned

# extract estimates
BSD.grp.results.16 <- BSD.df.hr.covar.AMPM.16$dht$clusters$N[3,]
BSD.ind.results.16 <- BSD.df.hr.covar.AMPM.16$dht$individuals$N[3,]

# bind results
BSD.results.2016 <- rbind(BSD.grp.results.16, BSD.ind.results.16)
BSD.results.2016[1,1] <- "Grp"
BSD.results.2016[2,1] <- "Ind"
BSD.results.2016$Year <- "2016"
BSD.results.2016$Species <- "BSD"
BSD.results.2016$DetFun <- "annual"
BSD.results.2016$Key <- "Hr"
BSD.results.2016$Adjust <- "NA"
BSD.results.2016$Covar <- "AM/PM"
BSD.results.2016 <- BSD.results.2016 %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)
BSD.results.2016 <- BSD.results.2016 %>% rename(N = Estimate)

### extract density
BSD.grp.density.16 <- BSD.df.hr.covar.AMPM.16$dht$clusters$D[3,]
BSD.ind.density.16 <- BSD.df.hr.covar.AMPM.16$dht$individuals$D[3,]
BSD.density.16 <- rbind(BSD.grp.density.16,BSD.ind.density.16)
BSD.density.16 <- BSD.density.16[, -c(1,7)]
BSD.density.16 <- BSD.density.16 %>% rename(D = Estimate)

# merge N and D
BSD.results.2016 <- cbind(BSD.results.2016,BSD.density.16)


## Binned

# extract estimates
BSD.grp.results.16.bin <- BSD.df.hr.ampm.16.bin$dht$clusters$N
BSD.ind.results.16.bin <- BSD.df.hr.ampm.16.bin$dht$individuals$N

# bind results
BSD.results.2016.bin <- rbind(BSD.grp.results.16.bin, BSD.ind.results.16.bin)
BSD.results.2016.bin$Label <- as.character(BSD.results.2016$Label)
BSD.results.2016.bin[1,1] <- "Grp"
BSD.results.2016.bin[2,1] <- "Ind"
BSD.results.2016.bin$Year <- "2016"
BSD.results.2016.bin$Species <- "BSD"
BSD.results.2016.bin$DetFun <- "annual"
BSD.results.2016.bin$Key <- "Hr"
BSD.results.2016.bin$Adjust <- "NA"
BSD.results.2016.bin$Covar <- "AM/PM"
BSD.results.2016.bin <- BSD.results.2016.bin %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)


### 2014 ####

### Previously, effort strata were ignored for this species, as they were for the other species. After discussions with Olly, we decided to check the effect of accounting for the effort strata for BSD in 2011, 2014, and 2016, as the estimates for these years were very high compared to other years, suggesting the strata were having an effect. I tested a couple of different approaches (see "Analytical_approach_test_CDS.R", lines 2190 onwards).  For BSD though, becasue each year is analysed separately, it is fairly simple (as no pooling is required). We think accounting for the strata here and in 2011 is the optimal approach.  

# This means that the old data will be loaded and used, as it includes the strata. Using the old data doesn't make any other difference here, because BSD are analysed annually so the 2014 data doesn't change between the old and new data.


# Number of groups in 2014
length(BSD.data$distance[BSD.data$year==2014]) # 535


  ## Subset data for 2014 ####

# old data being loaded and used here (see explanation in section above)

# I still need to change the survey area size and remove T20 as I do for all the data at the top of this script.

# load old data (2010-2018) that have strata
load("./Output/Data/Archive/KSWS_MASTER.Rdata")

allData$obs.habitat <- as.factor(allData$obs.habitat)
allData$obs.observer <- as.factor(allData$obs.observer)


## Remove T20
allData <- allData[allData$obs.transect != 20, ]
sample.table <- sample.table[sample.table$Sample.Label != 20,]
obs.table <- obs.table[obs.table$Sample.Label != 20,]


# subset BSD data
BSD.data <- allData[allData$species=="BSD",] 


# subset data for 2014
BSD.2014.data <- as.data.frame(BSD.data[BSD.data$stratum=="2014_H" | BSD.data$stratum=="2014_L", ]) 

region.table.2014 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2014_H" 
                                                     | full.region.table$Region.Label=="2014_L",])
# change survey area
new.area <- 1880000000/2
region.table.2014$Area <- new.area


sample.table.2014 <- as.data.frame(sample.table[sample.table$Region.Label =="2014_H" | 
                                                  sample.table$Region.Label =="2014_L",])

obs.table.2014 <- as.data.frame(obs.table[obs.table$Region.Label=="2014_H" | obs.table$Region.Label=="2014_L",])



  ## Exploratory plots & linear model ####

par(mfrow=c(1,2))

# Remove an outlier found in the plot to make it clearer
BSD.2014.data.forplot <- BSD.2014.data$distance[BSD.2014.data$distance<100]

# distance histograms
hist(BSD.2014.data.forplot, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(BSD.2014.data.forplot, main=NULL, xlab="Distance (m)", breaks=c(40))
# This histogram looks better than 2018 and 2016 - much less of a hump, suggesting evaisve movement isn't such an issue this year. Does show quite a steep drop off in detectability at around 40m. 


# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.BSD.14 <- 50 

# Count the number of observations discarded
nrow(BSD.2014.data[BSD.2014.data$distance>trunc.BSD.14,]) 

# number of groups
length(BSD.data$distance[BSD.data$stratum==4])

32/535*100 # 6%

# remove one outlier that is ruining the plots
BSD.2014.data.plot <- BSD.2014.data[BSD.2014.data$distance<100, ]

# Plot of distance against cluster size
par(mfrow=c(2,2))
plot(BSD.2014.data.plot$size, BSD.2014.data.plot$distance, main="a)", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.BSD14 <- lm(distance~size, data=BSD.2014.data.plot)
lines(BSD.2014.data.plot$size, as.vector(predict(lm.BSD14, BSD.2014.data.plot)))
summary(lm.BSD14)
# no strong evidence of size bias

# Plot of Observer factor against distance
plot(BSD.2014.data.plot$obs.observer,BSD.2014.data.plot$distance, main="b)", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# some variation between observers

# Plot of AMPM factor against distance
plot(BSD.2014.data.plot$obs.AMPM ,BSD.2014.data.plot$distance, main="c)", xlab="Time period",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# no bvious diferernece between AMPM

# Plot of habitat as factor against distance
plot(as.factor(BSD.2014.data.plot$obs.habitat), BSD.2014.data.plot$distance, 
     main="d)", xlab="Habitat class",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# similar pattern to 2018 - increasing distances as you move from evergreen through to open forest

  ## Fit detection functions ####
    # Uniform ####

# uniform cosine
BSD.df.unif.cos.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                      sample.table=sample.table.2014, obs.table=obs.table.2014, 
                      truncation=trunc.BSD.14, key="unif", adjustment= "cos")
summary(BSD.df.unif.cos.14)
# cosine(1) selected

# uniform poly
BSD.df.unif.poly.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                      sample.table=sample.table.2014, obs.table=obs.table.2014, 
                      truncation=trunc.BSD.14, key="unif", adjustment= "poly")
summary(BSD.df.unif.poly.14)
# simple polynomial(2,4) selected

# compare uniform models
bsd.uni.comp.14 <- summarize_ds_models(BSD.df.unif.cos.14,BSD.df.unif.poly.14,
                                    output="plain")
bsd.uni.comp.14[ ,1:4]
bsd.uni.comp.14[ ,5:7]
# All CvM p values are above 0.05. dAIC very similar

# plot the fits
par(mfrow=c(2,2))
plot(BSD.df.unif.cos.14, main = "BSD.df.unif.cos.14")

covar.fit <- ddf.gof(BSD.df.unif.cos.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

plot(BSD.df.unif.poly.14, main = "BSD.df.unif.poly.14")

covar.fit <- ddf.gof(BSD.df.unif.poly.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# cos model drops off quicker between 10 and 20m - is about 0.8 by 20m, whereas poly is still 0.9 at 20m. But cos drops slower between 0 and 10m, which is what I'm after.  The cos qq plot also looks better.

# BSD.df.unif.cos.14 selected


    # Half-normal ####

# HN cosine
BSD.df.hn.cos.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                    sample.table=sample.table.2014, obs.table=obs.table.2014,
                    truncation=trunc.BSD.14, key="hn", adjustment="cos")
summary(BSD.df.hn.cos)
# no adjustment selected

# HN hermite
BSD.df.hn.herm.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                    sample.table=sample.table.2014, obs.table=obs.table.2014,
                    truncation=trunc.BSD.14, key="hn", adjustment="herm")
summary(BSD.df.hn.herm.14)
# no adjustment selected


# Plot
par(mfrow=c(1,2))
plot(BSD.df.hn.cos.14, main = "BSD.df.hn.cos.14")

covar.fit <- ddf.gof(BSD.df.hn.cos.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# about 0.8 at 20m too. 

    # Hazard rate ####

# HR cosine
BSD.df.hr.cos.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                    sample.table=sample.table.2014, obs.table=obs.table.2014,
                    truncation=trunc.BSD.14, key="hr", adjustment="cos")
summary(BSD.df.hr.cos.14)
# cosine(2) selected

# HR poly
BSD.df.hr.poly.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                    sample.table=sample.table.2014, obs.table=obs.table.2014,
                    truncation=trunc.BSD.14, key="hr", adjustment="poly")
summary(BSD.df.hr.poly.14)
# no adjustment selected


bsd.df.hr.comp.14 <- summarize_ds_models(BSD.df.hr.cos.14, BSD.df.hr.poly.14, output = "plain")
bsd.df.hr.comp.14[ ,1:5]
bsd.df.hr.comp.14[ ,6:7]
# AIC is lower for HR cosine

# Plot both
par(mfrow=c(2,2))
plot(BSD.df.hr.cos.14, main = "BSD.df.hr.cos.14")
covar.fit <- ddf.gof(BSD.df.hr.cos.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM


plot(BSD.df.hr.poly.14, main = "BSD.df.hr.poly.14")
covar.fit <- ddf.gof(BSD.df.hr.poly.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly model suggests p = 1 right up to nearly 20m. This is too much for me

# BSD.df.hr.cos.14 selected

    # Compare primary models ####

bsd.df.prim.comp.14 <- summarize_ds_models(BSD.df.unif.cos.14, BSD.df.hn.cos.14, BSD.df.hr.cos.14, 
                                           output = "plain")
bsd.df.prim.comp.14[ ,1:5]
bsd.df.prim.comp.14[ ,6:7]

# both half normal and uniform models have some support (dAIC <2)

par(mfrow=c(2,2))

# HN
plot(BSD.df.hn.cos.14, main = "BSD.df.hn.cos.14")

covar.fit <- ddf.gof(BSD.df.hn.cos.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# uni
plot(BSD.df.unif.cos.14, main = "BSD.df.unif.cos.14")

covar.fit <- ddf.gof(BSD.df.unif.cos.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

## half normal model selected

    # Models with harsh truncation ####

# set trunc distance 
trunc.harsh.14 <- 40

par(mfrow=c(1,2))

# HN cos
BSD.df.hn.cos.harsh.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                    sample.table=sample.table.2014, obs.table=obs.table.2014,
                    truncation=trunc.harsh.14, key="hn", adjustment="cos")
summary(BSD.df.hn.cos.harsh.14)
# No adjustment selected

# compare original with harsh truncation
ddf.gof(BSD.df.hn.cos.harsh.14$ddf, main = "BSD.df.hn.cos.harsh.14", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(BSD.df.hn.cos.14$ddf, main = "BSD.df.hn.cos.14",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# The harsh model is no better, perhaps slightly worse fit

##  half normal model is still selected

    # Models with covariates ####

## cluster size
BSD.df.hn.size.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                    sample.table=sample.table.2014, obs.table=obs.table.2014,
                    truncation=trunc.BSD.14, key="hn", formula = ~size)
summary(BSD.df.hn.size.14)

# plot
par(mfrow=c(1,2))
plot(BSD.df.hn.size.14, main = "BSD.df.hn.size.14")

covar.fit <- ddf.gof(BSD.df.hn.size.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Habitat
BSD.df.hn.covar.hab.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                      sample.table=sample.table.2014, obs.table=obs.table.2014, 
                      truncation=trunc.BSD.14, key="hn", formula = ~obs.habitat)
summary(BSD.df.hn.covar.hab.14)

# plot
plot(BSD.df.hn.covar.hab.14, main = "BSD.df.hn.covar.hab.14")

covar.fit <- ddf.gof(BSD.df.hn.covar.hab.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## AM/PM
BSD.df.hn.covar.AMPM.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                      sample.table=sample.table.2014, obs.table=obs.table.2014, 
                      truncation=trunc.BSD.14, key="hn", formula = ~ obs.AMPM)
summary(BSD.df.hn.covar.AMPM.14)

# plot
plot(BSD.df.hn.covar.AMPM.14, main = "BSD.df.hn.covar.AMPM.14")

covar.fit <- ddf.gof(BSD.df.hn.covar.AMPM.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Observer
BSD.df.hn.covar.obs.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                      sample.table=sample.table.2014, obs.table=obs.table.2014, 
                      truncation=trunc.BSD.14, key="hn", formula = ~ obs.observer)
summary(BSD.df.hn.covar.obs.14)

# plot
plot(BSD.df.hn.covar.obs.14, main = "BSD.df.hn.covar.obs.14")

covar.fit <- ddf.gof(BSD.df.hn.covar.obs.14$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


    # Compare covariate models ####

bsd.df.cov.comp.14 <- summarize_ds_models(BSD.df.hn.size.14, BSD.df.hn.covar.obs.14, BSD.df.hn.covar.hab.14,
                                       BSD.df.hn.covar.AMPM.14, output = "plain")

bsd.df.cov.comp.14[ ,1:5]
bsd.df.cov.comp.14[ ,6:7]

# model with observer has overwhelming support.  I 'd like to test the inclusion of habitat though

BSD.df.hn.covar.obs.hab.14 <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                      sample.table=sample.table.2014, obs.table=obs.table.2014, 
                      truncation=trunc.BSD.14, key="hn", formula = ~obs.observer + obs.habitat)
summary(BSD.df.hn.covar.obs.hab.14)

# add to the summary
bsd.df.cov.comp.14 <- summarize_ds_models(BSD.df.hn.size.14, BSD.df.hn.covar.obs.14, BSD.df.hn.covar.hab.14,
                                       BSD.df.hn.covar.AMPM.14, BSD.df.hn.covar.obs.hab.14, 
                                       output = "plain")

bsd.df.cov.comp.14[ ,1:5]
bsd.df.cov.comp.14[ ,6:7]

# The model with observer and habitat has good support. Habitat is important in the observation process for primates. This model is selected

    # Compare primary and covariate models ####
  
bsd.df.final.compare.14 <- summarize_ds_models(BSD.df.hn.cos.14,BSD.df.hn.covar.obs.hab.14, 
                                               output = "plain")
bsd.df.final.compare.14[ ,1:5]
bsd.df.final.compare.14[ ,6:7]

### the model with covariates is the best and final model


  ## Final results 2014 ####

# extract estimates
BSD.grp.results.14 <- BSD.df.hn.covar.obs.hab.14$dht$clusters$N[3,]
BSD.ind.results.14 <- BSD.df.hn.covar.obs.hab.14$dht$individuals$N[3,]

# bind results
BSD.results.2014 <- rbind(BSD.grp.results.14, BSD.ind.results.14)
BSD.results.2014[1,1] <- "Grp"
BSD.results.2014[2,1] <- "Ind"
BSD.results.2014$Year <- "2014"
BSD.results.2014$Species <- "BSD"
BSD.results.2014$DetFun <- "annual"
BSD.results.2014$Key <- "Hn"
BSD.results.2014$Adjust <- "NA"
BSD.results.2014$Covar <- "Observer + Habitat"
BSD.results.2014 <- BSD.results.2014 %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)
BSD.results.2014 <- BSD.results.2014 %>% rename(N = Estimate)

### extract density 
BSD.grp.density.14 <- BSD.df.hn.covar.obs.hab.14$dht$clusters$D[3,]
BSD.ind.density.14 <- BSD.df.hn.covar.obs.hab.14$dht$individuals$D[3,]
BSD.density.14 <- rbind(BSD.grp.density.14,BSD.ind.density.14)
BSD.density.14 <- BSD.density.14[, -c(1,7)]
BSD.density.14 <- BSD.density.14 %>% rename(D = Estimate)

# merge N and D
BSD.results.2014 <- cbind(BSD.results.2014,BSD.density.14)


### 2013 ####

# Number of groups in 2013
length(BSD.data$distance[BSD.data$year==2013]) # 265

# roughly half the number of individuals in 2013 compared with 2014, and less than half the number of groups

  ## Subset data for 2013 ####

# If you are re-running all years make sure you re-load the new data (AllData), run the data preparation section, and re-subset BSD.data

BSD.2013.data <- as.data.frame(BSD.data[BSD.data$year==2013, ]) 
region.table.2013 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2013",])
sample.table.2013 <- as.data.frame(sample.table[sample.table$Region.Label =="2013",])
obs.table.2013 <- as.data.frame(obs.table[obs.table$Region.Label=="2013",])

  ## Exploratory plots and linear model ####

par(mfrow=c(1,2))

# distance histograms
hist(BSD.2013.data$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(BSD.2013.data$distance, main=NULL, xlab="Distance (m)", breaks=c(40))

# clumping around 10m again - evasive movement. Shoulder between 20-30m, which doesn't really reflect reality, and shouldn't be modelled

## Update - based on Olly finding issues with lumping on 0, we will be binning data for this year to remove the issue

# check bin widths
par(mfrow=c(1,2))
hist(BSD.2013.data$distance, breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80))
hist(BSD.2013.data$distance, breaks=c(0,10,20,30,40,50,60,70,80))
hist(BSD.2013.data$distance, breaks=c(0,7,14,21,28,35,42,49,56,63,70,80))
hist(BSD.2013.data$distance, breaks=c(0,8,16,24,32,40,48,56,64,72,80))

# I think the 4th one is best - it reduces the difference between the first and second bin


# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.BSD.13 <- 50 
trunc.BSD.13.bin <- 56


# Count the number of observations discarded
nrow(BSD.2013.data[BSD.2013.data$distance>trunc.BSD.13,]) 

length(BSD.data$distance[BSD.data$stratum==3])

17/265*100 # 6.4%

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(2,2))

plot(BSD.2013.data$size, BSD.2013.data$distance, main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.BSD13 <- lm(distance~size, data=BSD.2013.data)
lines(BSD.2013.data$size, as.vector(predict(lm.BSD13, BSD.2013.data)))
summary(lm.BSD13)
# no significant relationship between cluster size and distance

# Plot of Observer factor against distance
plot(BSD.2013.data$obs.observer,BSD.2013.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# some variation between observers

# Plot of AMPM factor against distance
plot(BSD.2013.data$obs.AMPM ,BSD.2013.data$distance, main="AM/PM", xlab="Time period",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# Slighter large mean distance in the mornings

# Plot of habitat as factor against distance
plot(as.factor(BSD.2013.data$obs.habitat), BSD.2013.data$distance, main="habitat", xlab="Habitat class",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# distances increase from EV to SEG to MX, but then drop again in DDF. MX has largest mean distance

# 1 = evergreen forest, 2 = semi-evergreen, 3 = mixed forest, 4 = deciduous forest, 5 = grass land, 6 = bamboo forest.  Larger distances in mixed forests


## plots of number of observations against covars

p1 <- ggplot(BSD.2013.data, aes(size))+geom_histogram(binwidth = 0.5)
p2 <- ggplot(BSD.2013.data, aes(obs.AMPM))+geom_histogram(stat="count")
p3 <- ggplot(BSD.2013.data, aes(obs.habitat))+geom_histogram(stat="count")
p4 <- ggplot(BSD.2013.data, aes(obs.observer))+geom_histogram(stat="count")
plot_grid(p1,p2,p3,p4)
# vast majority of observations are of groups size < 5, which of course is a massive underestimate for this species. This is why I don't think we should be reporting abundance of individuals. We should only report abundance of groups.  Not much difference in obs between morning and afternoon. Vast majority of observations in evergreen and semievergreen, as expected. 

  ## Fit detection functions ####
    # Uniform ####    

### unbinned 

# uniform cosine
BSD.df.unif.cos.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="unif", adjustment= "cos")
summary(BSD.df.unif.cos.13)
# cosine(1) selected

# uniform poly
BSD.df.unif.poly.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="unif", adjustment= "poly")
summary(BSD.df.unif.poly.13)
# simple polynomial(2) selected

# compare uniform models
bsd.uni.comp.13 <- summarize_ds_models(BSD.df.unif.cos.13,BSD.df.unif.poly.13,
                                    output="plain")
bsd.uni.comp.13[ ,1:4]
bsd.uni.comp.13[ ,5:7]
# All CvM p values are above 0.05. dAIC very similar.

# plot both fits
par(mfrow=c(2,2))

# Uni cos
plot(BSD.df.unif.cos.13, main = "BSD.df.unif.cos.13")

covar.fit <- ddf.gof(BSD.df.unif.cos.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# Uni poly
plot(BSD.df.unif.poly.13, main = "BSD.df.unif.poly.13")

covar.fit <- ddf.gof(BSD.df.unif.poly.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly model maintains higher p to a greater distance - doesn't get to 0.8 until >20m whereas cos gets to 0.8 around 15m. Poly p drops off quickly after 30m and has p of ~0 post 50m, whereas cos sees a flattening out of p at 0.2 between 40m and 50m. QQ plot for cos is better. Tough one, both have merits. I will select poly because it's more important to have a better model close to the line, and I believe that maintaining a high p up to 20m is ecologically more accurate.

# BSD.df.unif.poly.13 is selected



### Binned

# uniform cosine
BSD.df.unif.cos.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                         sample.table=sample.table.2013, obs.table=obs.table.2013, 
                         truncation=trunc.BSD.13.bin, key="unif", adjustment= "cos",
                         cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.unif.cos.13.bin)
# cosine(1)

# unform poly
BSD.df.unif.poly.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                             sample.table=sample.table.2013, obs.table=obs.table.2013, 
                             truncation=trunc.BSD.13.bin, key="unif", adjustment= "poly",
                             cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.unif.poly.13.bin)
# poly(2,4)

# compare uniform models
bsd.uni.comp.13.bin <- summarize_ds_models(BSD.df.unif.cos.13.bin,BSD.df.unif.poly.13.bin,
                                       output="plain")
bsd.uni.comp.13.bin[ ,1:4]
bsd.uni.comp.13.bin[ ,5:7]
# cos has most support

# plot the fits
plot(BSD.df.unif.cos.13.bin)
plot(BSD.df.unif.poly.13.bin)
# the cos fit is better closer to the line and so is selected


    # Half normal ####


### unbinned


# HN cosine
BSD.df.hn.cos.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="hn", adjustment= "cos")
summary(BSD.df.hn.cos.13)
# no adjustment selected 

# HN hermite
BSD.df.hn.herm.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="hn", adjustment= "herm")
summary(BSD.df.hn.cos.13)
# no adjustment selected 

## half normal with no adjustment is selected

# plot the fit
par(mfrow=c(1,2))
plot(BSD.df.hn.cos.13, main = "BSD.df.hn.cos.13")

covar.fit <- ddf.gof(BSD.df.hn.cos.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# perhaps a happy medium - p gets to 0.8 just before 20m. By 50m p= just under 0.2.



### binned

# HN cosine
BSD.df.hn.cos.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                       sample.table=sample.table.2013, obs.table=obs.table.2013, 
                       truncation=trunc.BSD.13.bin, key="hn", adjustment= "cos",
                       cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.hn.cos.13)
# no adjustment selected


# HN herm
BSD.df.hn.herm.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                           sample.table=sample.table.2013, obs.table=obs.table.2013, 
                           truncation=trunc.BSD.13.bin, key="hn", adjustment= "herm",
                           cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.hn.herm.13.bin)
# no adjustment selected

# plot the fit
plot(BSD.df.hn.cos.13.bin)
# looks pretty good



    # Hazard rate ####


### unbinned

# HR cosine
BSD.df.hr.cos.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="hr", adjustment= "cos")
summary(BSD.df.hr.cos.13)
# cosine(2) is selected 

# HR poly
BSD.df.hr.poly.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="hr", adjustment= "poly")
summary(BSD.df.hr.poly.13)
# poly(2) selected

# compare HR models
bsd.hr.comp.13 <- summarize_ds_models(BSD.df.hr.cos.13,BSD.df.hr.poly.13,
                                    output="plain")
bsd.hr.comp.13[ ,1:4]
bsd.hr.comp.13[ ,5:7]
# All CvM p values are above 0.05. dAIC all <2. Cosine is recommended, 


# plot the fits
par(mfrow=c(2,2))

# cos
plot(BSD.df.hr.cos.13, main = "BSD.df.hr.cos.13")

covar.fit <- ddf.gof(BSD.df.hr.cos.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(BSD.df.hr.poly.13, main = "BSD.df.hr.poly.13")

covar.fit <- ddf.gof(BSD.df.hr.poly.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# I don't like the look of cos, as it is modelling that hump between 20-30m which I said we didn't want to model

# BSD.df.hr.poly.13 is selected


### binned

# HR cosine
BSD.df.hr.cos.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                           sample.table=sample.table.2013, obs.table=obs.table.2013, 
                           truncation=trunc.BSD.13.bin, key="hr", adjustment= "cos",
                           cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.hr.cos.13.bin)
# cosine(2)


# HR poly
BSD.df.hr.poly.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                           sample.table=sample.table.2013, obs.table=obs.table.2013, 
                           truncation=trunc.BSD.13.bin, key="hr", adjustment= "poly",
                           cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.hr.poly.13.bin)
# poly(2)

# compare HR models
bsd.hr.comp.13.bin <- summarize_ds_models(BSD.df.hr.cos.13.bin,BSD.df.hr.poly.13.bin,
                                      output="plain")
bsd.hr.comp.13.bin[ ,1:4]
bsd.hr.comp.13.bin[ ,5:7]
# both models have similar support

# plot
plot(BSD.df.hr.cos.13.bin, main="BSD.df.hr.cos.13.bin")
plot(BSD.df.hr.poly.13.bin, main="BSD.df.hr.poly.13.bin")

# cos is overfitted.  Poly selected



    # Compare primary models ####

### unbinned

bsd.df.prim.comp.13 <- summarize_ds_models(BSD.df.unif.poly.13, BSD.df.hn.cos.13, BSD.df.hr.poly.13, 
                                           output = "plain")
bsd.df.prim.comp.13[ ,1:5]
bsd.df.prim.comp.13[ ,6:7]
# hazard rate has dAIC >2.

# plot uniform and HN together

par(mfrow=c(2,2))

# Uni poly
plot(BSD.df.unif.poly.13, main = "BSD.df.unif.poly.13")

covar.fit <- ddf.gof(BSD.df.unif.poly.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HN
plot(BSD.df.hn.cos.13, main = "BSD.df.hn.cos.13")

covar.fit <- ddf.gof(BSD.df.hn.cos.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# I prefer HN.  I think that the uniform model keep p too high for too long - it doen't get to 0.5 until nearly 40m.

# BSD.df.hn.cos.13 selected



### binned

bsd.df.prim.comp.13.bin <- summarize_ds_models(BSD.df.unif.cos.13.bin, BSD.df.hn.cos.13.bin, BSD.df.hr.poly.13.bin, 
                                           output = "plain")
bsd.df.prim.comp.13.bin[ ,1:5]
bsd.df.prim.comp.13.bin[ ,6:7]
# HR poly has least support. HN and UNI very similar

# plot them together
plot(BSD.df.unif.cos.13.bin, main="BSD.df.unif.cos.13.bin")
plot(BSD.df.hn.cos.13.bin, main="BSD.df.hn.cos.13.bin")
# basically identical


#  compare estimates
bsd.13.prim.bin.comp <- data.frame(model = c("uni","hn"),
                                  estimate = c(BSD.df.unif.cos.13.bin$dht$clusters$N$Estimate,
                                               BSD.df.hn.cos.13.bin$dht$clusters$N$Estimate),
                                  cv = c(BSD.df.unif.cos.13.bin$dht$clusters$N$cv,
                                         BSD.df.hn.cos.13.bin$dht$clusters$N$cv),
                                  se = c(BSD.df.unif.cos.13.bin$dht$clusters$N$se,
                                         BSD.df.hn.cos.13.bin$dht$clusters$N$se),
                                  lcl = c(BSD.df.unif.cos.13.bin$dht$clusters$N$lcl,
                                          BSD.df.hn.cos.13.bin$dht$clusters$N$lcl),
                                  ucl = c(BSD.df.unif.cos.13.bin$dht$clusters$N$ucl,
                                          BSD.df.hn.cos.13.bin$dht$clusters$N$ucl))
bsd.13.prim.bin.comp

ggplot(bsd.13.prim.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,9000)
# essentially identical.  I want to test covariates and therfore will go with HN


    # Model with harsh truncation ####

# set trunc distance
trunc.harsh.13 <- 35
trunc.harsh.13.bin <- 48


### unbinned

par(mfrow=c(1,2))

# HN cos
BSD.df.hn.cos.harsh.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.harsh.13, key="hn", adjustment= "cos")
summary(BSD.df.hn.cos.harsh.13)
# no adjustment selected

# compare original with harsh truncation
ddf.gof(BSD.df.hn.cos.harsh.13$ddf, main = "BSD.df.hn.cos.harsh.13", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(BSD.df.hn.cos.13$ddf, main = "BSD.df.hn.cos.13",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# Theoriginal model looks to be a better fit

# original selected



### binned

BSD.df.hn.cos.13.bin.harsh <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                           sample.table=sample.table.2013, obs.table=obs.table.2013, 
                           truncation=trunc.harsh.13.bin, key="hn",
                           cutpoints = c(0,8,16,24,32,40,48))

# plot together
plot(BSD.df.hn.cos.13.bin.harsh,main="BSD.df.hn.cos.13.bin.harsh")
plot(BSD.df.hn.cos.13.bin, main="BSD.df.hn.cos.13.bin")
# I think original looks to be a better fit

BSD.df.hn.cos.13.bin.harsh$dht$clusters$N$cv
BSD.df.hn.cos.13.bin$dht$clusters$N$cv
BSD.df.hn.cos.13.bin.harsh$dht$clusters$N$Estimate
BSD.df.hn.cos.13.bin$dht$clusters$N$Estimate
# very similar. CV for original is slightly smaller

# original selected

    # Models with covariates #####


### Unbinned

## cluster size
BSD.df.hn.size.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                    sample.table=sample.table.2013, obs.table=obs.table.2013,
                    truncation=trunc.BSD.13, key="hn", formula = ~size)
summary(BSD.df.hn.size.13)

# plot
par(mfrow=c(1,2))
plot(BSD.df.hn.size.13, main = "BSD.df.hn.size.13")

covar.fit <- ddf.gof(BSD.df.hn.size.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Habitat
BSD.df.hn.covar.hab.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="hn", formula = ~obs.habitat)
summary(BSD.df.hn.covar.hab.13)

# plot
plot(BSD.df.hn.covar.hab.13, main = "BSD.df.hn.covar.hab.13")

covar.fit <- ddf.gof(BSD.df.hn.covar.hab.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## AM/PM
BSD.df.hn.covar.AMPM.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="hn", formula = ~ obs.AMPM)
summary(BSD.df.hn.covar.AMPM.13)

# plot
plot(BSD.df.hn.covar.AMPM.13, main = "BSD.df.hn.covar.AMPM.13")

covar.fit <- ddf.gof(BSD.df.hn.covar.AMPM.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Observer
BSD.df.hn.covar.obs.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="hn", formula = ~ obs.observer)
summary(BSD.df.hn.covar.obs.13)

# plot
plot(BSD.df.hn.covar.obs.13, main = "BSD.df.hn.covar.obs.13")

covar.fit <- ddf.gof(BSD.df.hn.covar.obs.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## size and observer. First comparison identified size and observer as the two top models independently so I will try and combine the two
BSD.df.hn.covar.size.obs.13 <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                      sample.table=sample.table.2013, obs.table=obs.table.2013, 
                      truncation=trunc.BSD.13, key="hn", formula = ~ size + obs.observer)
summary(BSD.df.hn.covar.size.obs.13)

# plot
plot(BSD.df.hn.covar.size.obs.13, main = "BSD.df.hn.covar.size.obs.13")

covar.fit <- ddf.gof(BSD.df.hn.covar.size.obs.13$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



### binned

# size
BSD.df.hn.size.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                           sample.table=sample.table.2013, obs.table=obs.table.2013, 
                           truncation=trunc.BSD.13.bin, key="hn", formula = ~size,
                           cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.hn.size.13.bin)
plot(BSD.df.hn.size.13.bin)

# habitat
BSD.df.hn.hab.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                            sample.table=sample.table.2013, obs.table=obs.table.2013, 
                            truncation=trunc.BSD.13.bin, key="hn", formula = ~obs.habitat,
                            cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.hn.hab.13.bin)
plot(BSD.df.hn.hab.13.bin)

# AMPM
BSD.df.hn.ampm.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                           sample.table=sample.table.2013, obs.table=obs.table.2013, 
                           truncation=trunc.BSD.13.bin, key="hn", formula = ~obs.AMPM,
                           cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.hn.ampm.13.bin)
plot(BSD.df.hn.ampm.13.bin)

# observer 
BSD.df.hn.obs.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                            sample.table=sample.table.2013, obs.table=obs.table.2013, 
                            truncation=trunc.BSD.13.bin, key="hn", formula = ~obs.observer,
                            cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.hn.obs.13.bin)
plot(BSD.df.hn.obs.13.bin)

# ampm + size
BSD.df.hn.ampm.size.13.bin <- ds(data=BSD.2013.data, region.table=region.table.2013, 
                           sample.table=sample.table.2013, obs.table=obs.table.2013, 
                           truncation=trunc.BSD.13.bin, key="hn", formula = ~obs.AMPM + size,
                           cutpoints = c(0,8,16,24,32,40,48,56))
summary(BSD.df.hn.ampm.size.13.bin)
plot(BSD.df.hn.ampm.size.13.bin)



    # Compare covariate models ####

### unbinned

bsd.df.cov.comp.13 <- summarize_ds_models(BSD.df.hn.size.13, BSD.df.hn.covar.obs.13, 
                                          BSD.df.hn.covar.hab.13,BSD.df.hn.covar.AMPM.13,
                                          BSD.df.hn.covar.size.obs.13, output = "plain")

bsd.df.cov.comp.13[ ,1:5]
bsd.df.cov.comp.13[ ,6:7]

### BSD.df.hn.covar.size.obs is the best covariate model


### binned

bsd.df.cov.comp.13.bin <- summarize_ds_models(BSD.df.hn.size.13.bin,BSD.df.hn.hab.13.bin, 
                                              BSD.df.hn.ampm.13.bin,BSD.df.hn.obs.13.bin,
                                              BSD.df.hn.ampm.size.13.bin,
                                              output = "plain")

bsd.df.cov.comp.13.bin[ ,1:5]
bsd.df.cov.comp.13.bin[ ,6:7]
# all models except habitat have some support

# check estimates and CVs
bsd.13.cov.bin.comp <- data.frame(model = c("ampm","size", "ampm+size", "obs"),
                                  estimate = c(BSD.df.hn.ampm.13.bin$dht$clusters$N$Estimate,
                                               BSD.df.hn.size.13.bin$dht$clusters$N$Estimate,
                                               BSD.df.hn.ampm.size.13.bin$dht$clusters$N$Estimate,
                                               BSD.df.hn.obs.13.bin$dht$clusters$N$Estimate),
                                  cv = c(BSD.df.hn.ampm.13.bin$dht$clusters$N$cv,
                                         BSD.df.hn.size.13.bin$dht$clusters$N$cv,
                                         BSD.df.hn.ampm.size.13.bin$dht$clusters$N$cv,
                                         BSD.df.hn.obs.13.bin$dht$clusters$N$cv),
                                  se = c(BSD.df.hn.ampm.13.bin$dht$clusters$N$se,
                                         BSD.df.hn.size.13.bin$dht$clusters$N$se,
                                         BSD.df.hn.ampm.size.13.bin$dht$clusters$N$se,
                                         BSD.df.hn.obs.13.bin$dht$clusters$N$se),
                                  lcl = c(BSD.df.hn.ampm.13.bin$dht$clusters$N$lcl,
                                          BSD.df.hn.size.13.bin$dht$clusters$N$lcl,
                                          BSD.df.hn.ampm.size.13.bin$dht$clusters$N$lcl,
                                          BSD.df.hn.obs.13.bin$dht$clusters$N$lcl),
                                  ucl = c(BSD.df.hn.ampm.13.bin$dht$clusters$N$ucl,
                                          BSD.df.hn.size.13.bin$dht$clusters$N$ucl,
                                          BSD.df.hn.ampm.size.13.bin$dht$clusters$N$ucl,
                                          BSD.df.hn.obs.13.bin$dht$clusters$N$ucl))

bsd.13.cov.bin.comp

ggplot(bsd.13.cov.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl,ymax=ucl))+
  ylim(0,9200)
# basically no difference whatsoever.  obs model is slightly less precise.  I will go for the model that takes account of the most number of bias-inducing factors.

# BSD.df.hn.ampm.size.13.bin selected



    # Compare primary and covariate models ####
  
### unbinned

bsd.df.final.compare.13 <- summarize_ds_models(BSD.df.hn.cos.13,BSD.df.hn.covar.size.obs.13, 
                                               output = "plain")
bsd.df.final.compare.13[ ,1:3]
bsd.df.final.compare.13[ ,4:7]

### both models have support, and one could argue the use of the simpler model with no covariates. But Seeing as the model with covariate has strong support, I wouldn't want to exclude the covariates. We should have sufficient observations as well to deal with two covariates.

# BSD.df.hn.covar.size.obs.13 is the final model



### binned

# plot together
par(mfrow=c(1,2))
plot(BSD.df.hn.cos.13.bin, main="BSD.df.hn.cos.13.bin")
plot(BSD.df.hn.ampm.size.13.bin, main="BSD.df.hn.ampm.size.13.bin")

bsd.13.final.bin.comp <- data.frame(model = c("cov", "orig"),
                                    estimate = c(BSD.df.hn.ampm.size.13.bin$dht$clusters$N$Estimate,
                                                 BSD.df.hn.cos.13.bin$dht$clusters$N$Estimate),
                                    cv = c(BSD.df.hn.ampm.size.13.bin$dht$clusters$N$cv,
                                           BSD.df.hn.cos.13.bin$dht$clusters$N$cv),
                                    se = c(BSD.df.hn.ampm.size.13.bin$dht$clusters$N$se,
                                           BSD.df.hn.cos.13.bin$dht$clusters$N$se),
                                    lcl = c(BSD.df.hn.ampm.size.13.bin$dht$clusters$N$lcl,
                                            BSD.df.hn.cos.13.bin$dht$clusters$N$lcl),
                                    ucl = c(BSD.df.hn.ampm.size.13.bin$dht$clusters$N$ucl,
                                            BSD.df.hn.cos.13.bin$dht$clusters$N$ucl))

bsd.13.final.bin.comp

ggplot(bsd.13.final.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,9000)

# Basically identical.  Original model is a tiny wee bit more precise, and so is selected.


  ## Final results 2013 ####

### unbinned

# extract estimates
BSD.grp.results.13 <- BSD.df.hn.covar.size.obs.13$dht$clusters$N
BSD.ind.results.13 <- BSD.df.hn.covar.size.obs.13$dht$individuals$N

# bind results
BSD.results.2013 <- rbind(BSD.grp.results.13, BSD.ind.results.13)
BSD.results.2013$Label <- as.character(BSD.results.2013$Label)
BSD.results.2013[1,1] <- "Grp"
BSD.results.2013[2,1] <- "Ind"
BSD.results.2013$Year <- "2013"
BSD.results.2013 <- BSD.results.2013 %>% select(Year,Label,Estimate,se,cv,lcl,ucl)



### Binned
# extract estimates
BSD.grp.results.13.bin <- BSD.df.hn.cos.13.bin$dht$clusters$N
BSD.ind.results.13.bin <- BSD.df.hn.cos.13.bin$dht$individuals$N

# bind results
BSD.results.2013.bin <- rbind(BSD.grp.results.13.bin, BSD.ind.results.13.bin)
BSD.results.2013.bin$Label <- as.character(BSD.results.2013.bin$Label)
BSD.results.2013.bin[1,1] <- "Grp"
BSD.results.2013.bin[2,1] <- "Ind"
BSD.results.2013.bin$Year <- "2013"
BSD.results.2013.bin$Species <- "BSD"
BSD.results.2013.bin$DetFun <- "annual"
BSD.results.2013.bin$Key <- "Hn"
BSD.results.2013.bin$Adjust <- "NA"
BSD.results.2013.bin$Covar <- "NA"
BSD.results.2013.bin <- BSD.results.2013.bin %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)
BSD.results.2013.bin <- BSD.results.2013.bin %>% rename(N = Estimate)

### extract density
BSD.grp.density.13.bin <- BSD.df.hn.cos.13.bin$dht$clusters$D
BSD.ind.density.13.bin <- BSD.df.hn.cos.13.bin$dht$individuals$D
BSD.density.13.bin <- rbind(BSD.grp.density.13.bin,BSD.ind.density.13.bin)
BSD.density.13.bin <- BSD.density.13.bin[, -c(1,7)]
BSD.density.13.bin <- BSD.density.13.bin %>% rename(D = Estimate)

# merge N and D
BSD.results.2013.bin <- cbind(BSD.results.2013.bin,BSD.density.13.bin)



### 2011 ####

### Previously, effort strata were ignored for this species, as they were for the other species. After discussions with Olly, we decided to check the effect of accounting for the effort strata for BSD in 2011, 2014, and 2016, as the estimates for these years were very high compared to other years, suggesting the strata were having an effect. I tested a couple of different approaches (see "Analytical_approach_test_CDS.R", lines 2190 onwards).  For BSD though, becasue each year is analysed separately, it is fairly simple (as no pooling is required). We think accounting for the strata here and in 2014 is the optimal approach.  

# This means that the old data will be loaded and used, as it includes the strata. Using the old data doesn't make any other difference here, because BSD are analysed annually so the 2011 data doesn't change between the old and new data.



# Number of groups in 2011
length(BSD.data$distance[BSD.data$year==2011]) # 461

  ## subset data for 2011 ####

# old data being loaded and used here (see explanation in section above)

# I still need to change the survey area size and remove T20 as I do for all the data at the top of this script.

# load old data (2010-2018) that have strata
load("./Output/Data/Archive/KSWS_MASTER.Rdata")

allData$obs.habitat <- as.factor(allData$obs.habitat)
allData$obs.observer <- as.factor(allData$obs.observer)


## Remove T20
allData <- allData[allData$obs.transect != 20, ]
sample.table <- sample.table[sample.table$Sample.Label != 20,]
obs.table <- obs.table[obs.table$Sample.Label != 20,]


# subset BSD data
BSD.data <- allData[allData$species=="BSD",] 


# subset data for 2011
BSD.2011.data <- as.data.frame(BSD.data[BSD.data$stratum=="2011_H" | BSD.data$stratum=="2011_L", ]) 

region.table.2011 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2011_H" 
                                                     | full.region.table$Region.Label=="2011_L",])
# change survey area
new.area <- 1880000000/2
region.table.2011$Area <- new.area


sample.table.2011 <- as.data.frame(sample.table[sample.table$Region.Label =="2011_H" | 
                                                  sample.table$Region.Label =="2011_L",])

obs.table.2011 <- as.data.frame(obs.table[obs.table$Region.Label=="2011_H" | obs.table$Region.Label=="2011_L",])




  ## Exploratory plots and linear model ####

par(mfrow=c(1,2))

# Remove an outlier found in the plot to make it clearer
BSD.2011.data.forplot <- BSD.2011.data$distance[BSD.2011.data$distance<80]

# distance histograms
hist(BSD.2011.data.forplot, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(BSD.2011.data.forplot, main=NULL, xlab="Distance (m)", breaks=c(40))

# small clump jsut before 10m - evasive movement. Spike also just over 20m, which we shouln't model. Apart from that not bad data.


## Update - this year will be re-analysed with data binning

# test bin widths
par(mfrow=c(2,3))
hist(BSD.2011.data.forplot, breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70), freq = T)
hist(BSD.2011.data.forplot, breaks=c(0,10,20,30,40,50,60,70), freq = T)
hist(BSD.2011.data.forplot, breaks=c(0,7,14,21,28,35,42,49,56,63,70), freq = T)
hist(BSD.2011.data.forplot, breaks=c(0,8,16,24,32,40,48,56,64,70), freq = T)
hist(BSD.2011.data.forplot, breaks=c(0,4,8,12,16,20,25,30,35,40,45,50,60,70), freq = T)
hist(BSD.2011.data.forplot, breaks=c(0,2,4,6,9,12,15,18,21,24,27,30,33,36,39,42,45,50,60,70), freq = T)
hist(BSD.2011.data.forplot, breaks=c(0,1,4,6,8,10,13,16,19,21,24,26,28,30,32,34,36,39,42,45,50,60,70), freq = T)

# final choice
hist(BSD.2011.data.forplot[BSD.2011.data.forplot<55], breaks=c(0,3,7,12,17,22,27,32,37,42,47,55), freq = T)




# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.BSD.11 <- 50 
trunc.BSD.11.bin <- 55

# Count the number of observations discarded
nrow(BSD.2011.data[BSD.2011.data$distance>trunc.BSD.11,]) 

length(BSD.data$distance[BSD.data$stratum==2])

18/461*100 # 3.9%

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(2,2))
plot(BSD.2011.data$size, BSD.2011.data$distance, main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.BSD11 <- lm(distance~size, data=BSD.2011.data)
lines(BSD.2011.data$size, as.vector(predict(lm.BSD11, BSD.2011.data)))
summary(lm.BSD11)
# significant relationship between group size and distance.  Group size increases with increased distance

# Plot of Observer factor against distance
plot(BSD.2011.data$obs.observer,BSD.2011.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# some variation between observers

# time of day and habitat were not recorded in 2010 and 2011 so we can't plot them


## plots of number of observations against covars

p1 <- ggplot(BSD.2013.data, aes(size))+geom_histogram(binwidth = 0.5)
p2 <- ggplot(BSD.2013.data, aes(obs.observer))+geom_histogram(stat="count")
plot_grid(p1,p2)
# vast majority of observations are of groups size < 5, which of course is a massive underestimate for this species. This is why I don't think we should be reporting abundance of individuals. We should only report abundance of groups. 


  ## Fit a detection function ####
    # Uniform ####    

### unbinned

# uniform cosine
BSD.df.unif.cos.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                      sample.table=sample.table.2011, obs.table=obs.table.2011, 
                      truncation=trunc.BSD.11, key="unif", adjustment= "cos")
summary(BSD.df.unif.cos.11)
# cosine(1) selected

# uniform poly
BSD.df.unif.poly.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                      sample.table=sample.table.2011, obs.table=obs.table.2011, 
                      truncation=trunc.BSD.11, key="unif", adjustment= "poly")
summary(BSD.df.unif.poly.11)
# simple polynomial(2,4,6,8) selected

# compare uniform models
bsd.uni.comp.11 <- summarize_ds_models(BSD.df.unif.cos.11,BSD.df.unif.poly.11,
                                    output="plain")
bsd.uni.comp.11[ ,1:4]
bsd.uni.comp.11[ ,5:7]
# All CvM p values are above 0.05. dAIC similar

# plot all the fits
par(mfrow=c(2,2))

# uni cos
plot(BSD.df.unif.cos.11, main = "BSD.df.unif.cos.11")

covar.fit <- ddf.gof(BSD.df.unif.cos.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# uni poly
plot(BSD.df.unif.poly.11, main = "BSD.df.unif.poly.11")

covar.fit <- ddf.gof(BSD.df.unif.poly.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

## similar decrease in p clsoe to the line. poly appears to be trying to fit that clump between 20-30m which we don't want to do.  Cosine selected



### binned

# uniform cosine
BSD.df.unif.cos.11.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                         sample.table=sample.table.2011, obs.table=obs.table.2011, 
                         truncation=trunc.BSD.11.bin, key="unif", adjustment= "cos",
                         cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))
summary(BSD.df.unif.cos.11.bin)
# cosine(1) selected


# uniform poly
BSD.df.unif.poly.11.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                             sample.table=sample.table.2011, obs.table=obs.table.2011, 
                             truncation=trunc.BSD.11.bin, key="unif", adjustment= "poly",
                             cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))
summary(BSD.df.unif.poly.11.bin)
# poly(2,4) selected


# compare uniform models
bsd.uni.comp.11.bin <- summarize_ds_models(BSD.df.unif.cos.11.bin,BSD.df.unif.poly.11.bin,
                                       output="plain")
bsd.uni.comp.11.bin[ ,1:4]
bsd.uni.comp.11.bin[ ,5:7]
# both have some support

# plot both
plot(BSD.df.unif.cos.11.bin, main="cos")
plot(BSD.df.unif.poly.11.bin, main="poly")
# not great fits. but they both look basically the same so I will go with AIC

# BSD.df.unif.poly.11.bin selected


    # Half normal ####

### unbinned

# HN cosine
BSD.df.hn.cos.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                      sample.table=sample.table.2011, obs.table=obs.table.2011, 
                      truncation=trunc.BSD.11, key="hn", adjustment= "cos")
summary(BSD.df.hn.cos.11)
# no adjustment selected 

# HN hermite
BSD.df.hn.herm.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                      sample.table=sample.table.2011, obs.table=obs.table.2011, 
                      truncation=trunc.BSD.11, key="hn", adjustment= "herm")
summary(BSD.df.hn.cos.11)
# no adjustment selected 

## half normal with no adjustment is selected

# plot the fit
par(mfrow=c(1,2))
plot(BSD.df.hn.cos.11, main = "BSD.df.hn.cos.11")

covar.fit <- ddf.gof(BSD.df.hn.cos.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# I feel like this is a decent fit - it is not trying to fit the clumps, and the shoulder close to the line looks realistic 

# BSD.df.hn.cos.11 selected



### binned


# HN cosine
BSD.df.hn.cos.11.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                             sample.table=sample.table.2011, obs.table=obs.table.2011, 
                             truncation=trunc.BSD.11.bin, key="hn", adjustment= "cos",
                             cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))
summary(BSD.df.hn.cos.11.bin)
# no adjustment selected

# HN herm
BSD.df.hn.herm.11.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                           sample.table=sample.table.2011, obs.table=obs.table.2011, 
                           truncation=trunc.BSD.11.bin, key="hn", adjustment= "herm",
                           cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))
summary(BSD.df.hn.herm.11.bin)
# no adjustment selected

# plot
plot(BSD.df.hn.cos.11.bin)



    # Hazard rate ####

### unbinned

# HR cosine
BSD.df.hr.cos.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                      sample.table=sample.table.2011, obs.table=obs.table.2011, 
                      truncation=trunc.BSD.11, key="hr", adjustment= "cos")
summary(BSD.df.hr.cos.11)
# cosine(2) is selected 

# HR poly
BSD.df.hr.poly.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                      sample.table=sample.table.2011, obs.table=obs.table.2011, 
                      truncation=trunc.BSD.11, key="hr", adjustment= "poly")
summary(BSD.df.hr.poly.11)
# poly(2,4) selected

# compare HR models
bsd.hr.comp.11 <- summarize_ds_models(BSD.df.hr.cos.11,BSD.df.hr.poly.11,
                                    output="plain")
bsd.hr.comp.11[ ,1:4]
bsd.hr.comp.11[ ,5:7]
# All CvM p values are above 0.05. poly has dAIC > 2


# plot the fits
par(mfrow=c(2,2))

# cos
plot(BSD.df.hr.cos.11, main = "BSD.df.hr.cos.11")

covar.fit <- ddf.gof(BSD.df.hr.cos.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# poly
plot(BSD.df.hr.poly.11, main = "BSD.df.hr.poly.11")

covar.fit <- ddf.gof(BSD.df.hr.poly.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# cos trying to fit that clump between 20-30m, which is not great. However, poly decreases too fast - by 20m p = 0.6 which is not realistic. Therefore I wil select cos



### binned

# HR cosine
BSD.df.hr.cos.11.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                           sample.table=sample.table.2011, obs.table=obs.table.2011, 
                           truncation=trunc.BSD.11.bin, key="hr", adjustment= "cos",
                           cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))
summary(BSD.df.hr.cos.11.bin)
# cosine(2)


# HR poly
BSD.df.hr.poly.11.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                           sample.table=sample.table.2011, obs.table=obs.table.2011, 
                           truncation=trunc.BSD.11.bin, key="hr", adjustment= "poly",
                           cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))
summary(BSD.df.hr.cos.11.bin)
# poly(2,4)


# compare HR models
bsd.hr.comp.11.bin <- summarize_ds_models(BSD.df.hr.cos.11.bin,BSD.df.hr.poly.11.bin,
                                      output="plain")
bsd.hr.comp.11.bin[ ,1:4]
bsd.hr.comp.11.bin[ ,5:7]
# similar support

# plot
plot(BSD.df.hr.cos.11.bin, main="cos")
plot(BSD.df.hr.poly.11.bin, main="poly")
# cos is overfitted

# BSD.df.hr.poly.11.bin selected


    # Compare primary models ####

### unbinned


bsd.df.prim.comp.11 <- summarize_ds_models(BSD.df.unif.cos.11, BSD.df.hn.cos.11, BSD.df.hr.cos.11, 
                                           output = "plain")
bsd.df.prim.comp.11[ ,1:5]
bsd.df.prim.comp.11[ ,6:7]
# Uniform model has dAIC >2 and so is not considered. HN and HR both have support

# plot them together
par(mfrow=c(2,2))

# HN
plot(BSD.df.hn.cos.11, main = "BSD.df.hn.cos.11")

covar.fit <- ddf.gof(BSD.df.hn.cos.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# HR
plot(BSD.df.hr.cos.11, main = "BSD.df.hr.cos.11")

covar.fit <- ddf.gof(BSD.df.hr.cos.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# The HR model QQ plot looks better, but the model looks overfitted. Seeing as both models have support, I will select half normal


### binned

bsd.df.prim.comp.11.bin <- summarize_ds_models(BSD.df.unif.poly.11.bin, BSD.df.hn.cos.11.bin, 
                                               BSD.df.hr.poly.11.bin, output = "plain")
bsd.df.prim.comp.11.bin[ ,1:5]
bsd.df.prim.comp.11.bin[ ,6:7]
# All have some support

# plot together
par(mfrow=c(2,2))
plot(BSD.df.unif.cos.11.bin, main="unif")
plot(BSD.df.hn.cos.11.bin, main="hn")
plot(BSD.df.hr.poly.11.bin, main="hr")
# Pretty simlar.  HR models more of the shoulder


#  compare estimates
bsd.11.prim.bin.comp <- data.frame(model = c("uni","hn","hr"),
                                   estimate = c(BSD.df.unif.cos.11.bin$dht$clusters$N$Estimate,
                                                BSD.df.hn.cos.11.bin$dht$clusters$N$Estimate,
                                                BSD.df.hr.poly.11.bin$dht$clusters$N$Estimate),
                                   cv = c(BSD.df.unif.cos.11.bin$dht$clusters$N$cv,
                                          BSD.df.hn.cos.11.bin$dht$clusters$N$cv,
                                          BSD.df.hr.poly.11.bin$dht$clusters$N$cv),
                                   se = c(BSD.df.unif.cos.11.bin$dht$clusters$N$se,
                                          BSD.df.hn.cos.11.bin$dht$clusters$N$se,
                                          BSD.df.hr.poly.11.bin$dht$clusters$N$se),
                                   lcl = c(BSD.df.unif.cos.11.bin$dht$clusters$N$lcl,
                                           BSD.df.hn.cos.11.bin$dht$clusters$N$lcl,
                                           BSD.df.hr.poly.11.bin$dht$clusters$N$lcl),
                                   ucl = c(BSD.df.unif.cos.11.bin$dht$clusters$N$ucl,
                                           BSD.df.hn.cos.11.bin$dht$clusters$N$ucl,
                                           BSD.df.hr.poly.11.bin$dht$clusters$N$ucl))
bsd.11.prim.bin.comp

ggplot(bsd.11.prim.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,16200)

## HR is the least precise, and so HN is the best compromise


    # Models with harsh truncation ####

# set trunc distance to 40m
trunc.harsh.11 <- 40
trunc.harsh.11.bin <- 42


### unbinned

par(mfrow=c(1,2))

# HN cos
BSD.df.hn.cos.harsh.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                    sample.table=sample.table.2011, obs.table=obs.table.2011,
                    truncation=trunc.harsh.11, key="hn", adjustment="cos")
summary(BSD.df.hn.cos.harsh.11)
# Cosine(2,3) selected

# compare original with harsh truncation
ddf.gof(BSD.df.hn.cos.harsh.11$ddf, main = "BSD.df.hn.cos.harsh.11", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(BSD.df.hn.cos.11$ddf, main = "BSD.df.hn.cos.11",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# The harsh truncation appears to have improved the model fit. 

# plot both fits
par(mfrow=c(2,2))

# HN
plot(BSD.df.hn.cos.11, main = "BSD.df.hn.cos.11")

covar.fit <- ddf.gof(BSD.df.hn.cos.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# harsh
plot(BSD.df.hn.cos.harsh.11, main = "BSD.df.hn.cos.harsh.11")

covar.fit <- ddf.gof(BSD.df.hn.cos.harsh.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)


## No no no, the harsh truncation model is ovefitting the data, and grovelling to the clump

# original model selected



### binned

# HN cosine
BSD.df.hn.cos.11.harsh.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                           sample.table=sample.table.2011, obs.table=obs.table.2011, 
                           truncation=trunc.harsh.11.bin, key="hn", adjustment= "cos",
                           cutpoints = c(0,3,7,12,17,22,27,32,37,42))

# plot together
plot(BSD.df.hn.cos.11.harsh.bin, main="harsh")
plot(BSD.df.hn.cos.11.bin, main="orig")

#  compare estimates
bsd.11.harsh.bin.comp <- data.frame(model = c("harsh","orig"),
                                   estimate = c(BSD.df.hn.cos.11.harsh.bin$dht$clusters$N$Estimate,
                                                BSD.df.hn.cos.11.bin$dht$clusters$N$Estimate),
                                   cv = c(BSD.df.hn.cos.11.harsh.bin$dht$clusters$N$cv,
                                          BSD.df.hn.cos.11.bin$dht$clusters$N$cv),
                                   se = c(BSD.df.hn.cos.11.harsh.bin$dht$clusters$N$se,
                                          BSD.df.hn.cos.11.bin$dht$clusters$N$se),
                                   lcl = c(BSD.df.hn.cos.11.harsh.bin$dht$clusters$N$lcl,
                                           BSD.df.hn.cos.11.bin$dht$clusters$N$lcl),
                                   ucl = c(BSD.df.hn.cos.11.harsh.bin$dht$clusters$N$ucl,
                                           BSD.df.hn.cos.11.bin$dht$clusters$N$ucl))
bsd.11.harsh.bin.comp

ggplot(bsd.11.harsh.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,16200)
# very similar. orig is a tiny bit more precise and so is retained


    # Models with covariates ####

### unbinned

## cluster size
BSD.df.hn.size.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                    sample.table=sample.table.2011, obs.table=obs.table.2011,
                    truncation=trunc.BSD.11, key="hn", formula = ~size)
summary(BSD.df.hn.size.11)

# plot
par(mfrow=c(1,2))
plot(BSD.df.hn.size.11, main = "BSD.df.hn.size.11")

covar.fit <- ddf.gof(BSD.df.hn.size.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



## Observer
BSD.df.hn.covar.obs.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                      sample.table=sample.table.2011, obs.table=obs.table.2011, 
                      truncation=trunc.BSD.11, key="hn", formula = ~ obs.observer)
summary(BSD.df.hn.covar.obs.11)

# plot
plot(BSD.df.hn.covar.obs.11, main = "BSD.df.hn.covar.obs.11")

covar.fit <- ddf.gof(BSD.df.hn.covar.obs.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## obsever and size
BSD.df.hn.obs.size.11 <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                      sample.table=sample.table.2011, obs.table=obs.table.2011, 
                      truncation=trunc.BSD.11, key="hn", formula = ~ obs.observer + size)
summary(BSD.df.hn.obs.size.11)

# plot
plot(BSD.df.hn.obs.size.11, main = "BSD.df.hn.obs.size.11")

covar.fit <- ddf.gof(BSD.df.hn.obs.size.11$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



### binned

par(mfrow=c(1,1))

# size
BSD.df.hn.size.11.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                           sample.table=sample.table.2011, obs.table=obs.table.2011, 
                           truncation=trunc.BSD.11.bin, key="hn", formula=~size,
                           cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))
summary(BSD.df.hn.size.11.bin)
plot(BSD.df.hn.size.11.bin)


# observer
BSD.df.hn.obs.11.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                            sample.table=sample.table.2011, obs.table=obs.table.2011, 
                            truncation=trunc.BSD.11.bin, key="hn", formula=~obs.observer,
                            cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))
summary(BSD.df.hn.obs.11.bin)
plot(BSD.df.hn.obs.11.bin)


# size and observer
BSD.df.hn.obs.size.11.bin <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                           sample.table=sample.table.2011, obs.table=obs.table.2011, 
                           truncation=trunc.BSD.11.bin, key="hn", formula=~obs.observer+size,
                           cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))
summary(BSD.df.hn.obs.size.11.bin)
plot(BSD.df.hn.obs.size.11.bin)


    # Compare covariate models ####


### unbinned

bsd.df.cov.comp.11 <- summarize_ds_models(BSD.df.hn.size.11, BSD.df.hn.covar.obs.11,BSD.df.hn.obs.size.11, 
                                          output = "plain")

bsd.df.cov.comp.11[ ,1:5]
bsd.df.cov.comp.11[ ,6:7]

# model with observer has overwhelming support. Odd, seeing as the linear model for size bias suggested a relationship. I will test a model with both

# ok so the model with both is the better model. I want to inlcude size because of the linear model result.

# BSD.df.hn.obs.size.11 is selected



### binned
bsd.df.cov.comp.11.bin <- summarize_ds_models(BSD.df.hn.size.11.bin,BSD.df.hn.obs.11.bin, 
                                              BSD.df.hn.obs.size.11.bin, output = "plain")

bsd.df.cov.comp.11.bin[ ,1:5]
bsd.df.cov.comp.11.bin[ ,6:7]
# size is very poor model, but observer and observer+size have support. I am suspicious of the observer+size model based on the size model AIC

# check estimates
bsd.11.cov.bin.comp <- data.frame(model = c("size","obs","both"),
                                   estimate = c(BSD.df.hn.size.11.bin$dht$clusters$N$Estimate,
                                                BSD.df.hn.obs.11.bin$dht$clusters$N$Estimate,
                                                BSD.df.hn.obs.size.11.bin$dht$clusters$N$Estimate),
                                   cv = c(BSD.df.hn.size.11.bin$dht$clusters$N$cv,
                                          BSD.df.hn.obs.11.bin$dht$clusters$N$cv,
                                          BSD.df.hn.obs.size.11.bin$dht$clusters$N$cv),
                                   se = c(BSD.df.hn.size.11.bin$dht$clusters$N$se,
                                          BSD.df.hn.obs.11.bin$dht$clusters$N$se,
                                          BSD.df.hn.obs.size.11.bin$dht$clusters$N$se),
                                   lcl = c(BSD.df.hn.size.11.bin$dht$clusters$N$lcl,
                                           BSD.df.hn.obs.11.bin$dht$clusters$N$lcl,
                                           BSD.df.hn.obs.size.11.bin$dht$clusters$N$lcl),
                                   ucl = c(BSD.df.hn.size.11.bin$dht$clusters$N$ucl,
                                           BSD.df.hn.obs.11.bin$dht$clusters$N$ucl,
                                           BSD.df.hn.obs.size.11.bin$dht$clusters$N$ucl))
bsd.11.cov.bin.comp

ggplot(bsd.11.cov.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,16500)

# I don't really know why the size model has such a high dAIC...the CV and SE are lower than the other two, and the model is more precise.

# screw it I'm going with size


    # Compare primary and covariate models ####

### unbinned

bsd.df.final.compare.11 <- summarize_ds_models(BSD.df.hn.cos.11,BSD.df.hn.obs.size.11, output = "plain")
bsd.df.final.compare.11[ ,1:3]
bsd.df.final.compare.11[ ,4:7]

# Model with covariates has overwhelming support


### binned

bsd.df.final.compare.11.bin <- summarize_ds_models(BSD.df.hn.cos.11.bin,BSD.df.hn.size.11.bin, output = "plain")
bsd.df.final.compare.11.bin[ ,1:3]
bsd.df.final.compare.11.bin[ ,4:7]
# covariate model has most support

# check estimates
bsd.11.final.bin.comp <- data.frame(model = c("cov","orig"),
                                  estimate = c(BSD.df.hn.size.11.bin$dht$clusters$N$Estimate,
                                               BSD.df.hn.cos.11.bin$dht$clusters$N$Estimate),
                                  cv = c(BSD.df.hn.size.11.bin$dht$clusters$N$cv,
                                         BSD.df.hn.cos.11.bin$dht$clusters$N$cv),
                                  se = c(BSD.df.hn.size.11.bin$dht$clusters$N$se,
                                         BSD.df.hn.cos.11.bin$dht$clusters$N$se),
                                  lcl = c(BSD.df.hn.size.11.bin$dht$clusters$N$lcl,
                                          BSD.df.hn.cos.11.bin$dht$clusters$N$lcl),
                                  ucl = c(BSD.df.hn.size.11.bin$dht$clusters$N$ucl,
                                          BSD.df.hn.cos.11.bin$dht$clusters$N$ucl))
bsd.11.final.bin.comp

ggplot(bsd.11.final.bin.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,15000)

# both models produce very similar esimates, with virtually identical CV and SE.  The original model is slightly more precies, and so is selected


  ## Final results 2011 ####

### unbinned


# extract estimates
BSD.grp.results.11 <- BSD.df.hn.obs.size.11$dht$clusters$N
BSD.ind.results.11 <- BSD.df.hn.obs.size.11$dht$individuals$N

# bind results
BSD.results.2011 <- rbind(BSD.grp.results.11, BSD.ind.results.11)
BSD.results.2011$Label <- as.character(BSD.results.2011$Label)
BSD.results.2011[1,1] <- "Grp"
BSD.results.2011[2,1] <- "Ind"
BSD.results.2011$Year <- "2011"
BSD.results.2011 <- BSD.results.2011 %>% select(Year,Label,Estimate,se,cv,lcl,ucl)



### binned

# extract estimates
BSD.grp.results.11.bin <- BSD.df.hn.cos.11.bin$dht$clusters$N[3,]
BSD.ind.results.11.bin <- BSD.df.hn.cos.11.bin$dht$individuals$N[3,]

# bind results
BSD.results.2011.bin <- rbind(BSD.grp.results.11.bin, BSD.ind.results.11.bin)
BSD.results.2011.bin[1,1] <- "Grp"
BSD.results.2011.bin[2,1] <- "Ind"
BSD.results.2011.bin$Year <- "2011"
BSD.results.2011.bin$Species <- "BSD"
BSD.results.2011.bin$DetFun <- "annual"
BSD.results.2011.bin$Key <- "Hn"
BSD.results.2011.bin$Adjust <- "NA"
BSD.results.2011.bin$Covar <- "NA"
BSD.results.2011.bin <- BSD.results.2011.bin %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)
BSD.results.2011.bin <- BSD.results.2011.bin %>% rename(N = Estimate)


### extract density
BSD.grp.density.11.bin <- BSD.df.hn.cos.11.bin$dht$clusters$D[3,]
BSD.ind.density.11.bin <- BSD.df.hn.cos.11.bin$dht$individuals$D[3,]
BSD.density.11.bin <- rbind(BSD.grp.density.11.bin,BSD.ind.density.11.bin)
BSD.density.11.bin <- BSD.density.11.bin[, -c(1,7)]
BSD.density.11.bin <- BSD.density.11.bin %>% rename(D = Estimate)

# merge N and D
BSD.results.2011.bin <- cbind(BSD.results.2011.bin,BSD.density.11.bin)


### 2010 ####

# Number of groups in 2010
length(BSD.data$distance[BSD.data$year==2010]) # 333


  ## Subset data for 2010 ####

BSD.2010.data <- as.data.frame(BSD.data[BSD.data$year==2010, ]) 
region.table.2010 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2010",])
sample.table.2010 <- as.data.frame(sample.table[sample.table$Region.Label =="2010",])
obs.table.2010 <- as.data.frame(obs.table[obs.table$Region.Label=="2010",])

  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

# distance histograms
hist(BSD.2010.data$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(BSD.2010.data$distance, main=NULL, xlab="Distance (m)", breaks=c(40))

# Large number of observations at <5m. Dip in observations between 25-35m.  Long tail beyond 60m


### update - data for this year will be binned to get rid of the clumping at 0

# find appropriate bins
par(mfrow=c(2,2))
hist(BSD.2010.data$distance[BSD.2010.data$distance<60], 
     breaks=c(0,2,4,6,8,10,15,20,25,30,35,40,45,50,55,60), freq=T)

hist(BSD.2010.data$distance[BSD.2010.data$distance<60], 
     breaks=c(0,2,6,11,14,17,21,24,28,35,40,45,50,55,60), freq=T)

hist(BSD.2010.data$distance[BSD.2010.data$distance<60], 
     breaks=c(0,1,6,12,16,21,26,33,40,50,60), freq=T)

hist(BSD.2010.data$distance[BSD.2010.data$distance<60], 
     breaks=c(0,1,5,9,13,16,19,22,25,30,35,40,45,55,60), freq=T)

hist(BSD.2010.data$distance[BSD.2010.data$distance<60], 
     breaks=c(0,10,20,30,40,50,60), freq=T)

# this one is tricky. I think the final one is actually the best, despite losing the most amount of information




# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.BSD.10 <- 50 
trunc.BSD.10.bin <- 60

# Count the number of observations discarded
nrow(BSD.2010.data[BSD.2010.data$distance>trunc.BSD.10,]) 

length(BSD.data$distance[BSD.data$stratum==1])

17/333*100 # 5%

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(1,2))
plot(BSD.2010.data$size, BSD.2010.data$distance, main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.BSD10 <- lm(distance~size, data=BSD.2010.data)
lines(BSD.2010.data$size, as.vector(predict(lm.BSD10, BSD.2010.data)))
summary(lm.BSD10)
# no significant relationship between group size and distance.  

# Plot of Observer factor against distance
plot(BSD.2010.data$obs.observer,BSD.2010.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# time of day and habitat were not recorded in 2010 and 2011 so we can't plot them


## plots of number of observations against covars

p1 <- ggplot(BSD.2010.data, aes(size))+geom_histogram(binwidth = 0.5)
p2 <- ggplot(BSD.2010.data, aes(obs.observer))+geom_histogram(stat="count")
plot_grid(p1,p2)
# Majority of observations have group size of <5


  ## Fit a detection function ####
    # Uniform ####    

### unbinned

# uniform cosine
BSD.df.unif.cos.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                      sample.table=sample.table.2010, obs.table=obs.table.2010, 
                      truncation=trunc.BSD.10, key="unif", adjustment= "cos")
summary(BSD.df.unif.cos.10)
# cosine(1) selected


# uniform poly
BSD.df.unif.poly.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                      sample.table=sample.table.2010, obs.table=obs.table.2010, 
                      truncation=trunc.BSD.10, key="unif", adjustment= "poly")
summary(BSD.df.unif.poly.10)
# simple polynomial(2,4,6,8) selected

# compare uniform models
bsd.uni.comp.10 <- summarize_ds_models(BSD.df.unif.cos.10,BSD.df.unif.poly.10,
                                    output="plain")
bsd.uni.comp.10[ ,1:4]
bsd.uni.comp.10[ ,5:7]
# Cos is excluded as it has CvM p values of <0.05. 

# plot poly fit
par(mfrow=c(1,2))

plot(BSD.df.unif.poly.10, main = "BSD.df.unif.poly.10")

covar.fit <- ddf.gof(BSD.df.unif.poly.10$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# not great looking. I think p decreases too quickly and the hump around 40m is not right. Perhaps the harsh truncation model will be better


### binned


# uniform cosine
BSD.df.unif.cos.10.bin <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                         sample.table=sample.table.2010, obs.table=obs.table.2010, 
                         truncation=trunc.BSD.10.bin, key="unif", adjustment= "cos",
                         cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.unif.cos.10.bin)
# cosine(1,2) selected

# uniform poly
BSD.df.unif.poly.10.bin <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                             sample.table=sample.table.2010, obs.table=obs.table.2010, 
                             truncation=trunc.BSD.10.bin, key="unif", adjustment= "poly",
                             cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.unif.poly.10.bin)
# poly(2,4,6) selected

# compare uniform models
bsd.uni.comp.10.bin <- summarize_ds_models(BSD.df.unif.cos.10.bin,BSD.df.unif.poly.10.bin,
                                       output="plain")
bsd.uni.comp.10.bin[ ,1:4]
bsd.uni.comp.10.bin[ ,5:7]
# cos has more support

# plot both
par(mfrow=c(1,2))
plot(BSD.df.unif.cos.10.bin, main="cos")
plot(BSD.df.unif.poly.10.bin, main="poly")
# poly is doing some weird overfit at the tail end

# cos selected

    # Half normal ####


### unbinned

# HN cosine
BSD.df.hn.cos.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                      sample.table=sample.table.2010, obs.table=obs.table.2010, 
                      truncation=trunc.BSD.10, key="hn", adjustment= "cos")
summary(BSD.df.hn.cos.10)
# no adjustment selected 

# HN hermite
BSD.df.hn.herm.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                      sample.table=sample.table.2010, obs.table=obs.table.2010, 
                      truncation=trunc.BSD.10, key="hn", adjustment= "herm")
summary(BSD.df.hn.cos.10)
# no adjustment selected 


## half normal with no adjustment is selected

# plot the fit
par(mfrow=c(1,2))
plot(BSD.df.hn.cos.10, main = "BSD.df.hn.cos.10")

covar.fit <- ddf.gof(BSD.df.hn.cos.10$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# this looks better. p doesn't decrease as quickly, and the hump at 40m is not fitted



### binned

# HN cosine
BSD.df.hn.cos.10.bin <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                             sample.table=sample.table.2010, obs.table=obs.table.2010, 
                             truncation=trunc.BSD.10.bin, key="hn", adjustment= "cos",
                             cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.hn.cos.10.bin)
# cosine(2) selected


# HN hermite
BSD.df.hn.herm.10.bin <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                           sample.table=sample.table.2010, obs.table=obs.table.2010, 
                           truncation=trunc.BSD.10.bin, key="hn", adjustment= "herm",
                           cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.hn.herm.10.bin)
# no adjustment selected

# compare HN models
bsd.hn.comp.10.bin <- summarize_ds_models(BSD.df.hn.cos.10.bin,BSD.df.hn.herm.10.bin,
                                           output="plain")
bsd.hn.comp.10.bin[ ,1:4]
bsd.hn.comp.10.bin[ ,5:7]
# very similar AIC

# plot both
plot(BSD.df.hn.cos.10.bin, main="cos")
plot(BSD.df.hn.herm.10.bin, main="key only")

# check estimates
bsd.hn.comp <- data.frame(model=c("cos", "key only"),
                          estimate = c(BSD.df.hn.cos.10.bin$dht$clusters$N$Estimate,
                                       BSD.df.hn.herm.10.bin$dht$clusters$N$Estimate),
                          cv = c(BSD.df.hn.cos.10.bin$dht$clusters$N$cv,
                                 BSD.df.hn.herm.10.bin$dht$clusters$N$cv),
                          se = c(BSD.df.hn.cos.10.bin$dht$clusters$N$se,
                                 BSD.df.hn.herm.10.bin$dht$clusters$N$se),
                          lcl = c(BSD.df.hn.cos.10.bin$dht$clusters$N$lcl,
                                  BSD.df.hn.herm.10.bin$dht$clusters$N$lcl),
                          ucl = c(BSD.df.hn.cos.10.bin$dht$clusters$N$ucl,
                                  BSD.df.hn.herm.10.bin$dht$clusters$N$ucl))

bsd.hn.comp

ggplot(bsd.hn.comp, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,11200)
# key only model is more precise. All other things being equal, I will select that


    # Hazard rate ####

### unbinned

# HR cosine
BSD.df.hr.cos.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                      sample.table=sample.table.2010, obs.table=obs.table.2010, 
                      truncation=trunc.BSD.10, key="hr", adjustment= "cos")
summary(BSD.df.hr.cos.10)
# cosine(2) is selected 


# HR poly
BSD.df.hr.poly.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                      sample.table=sample.table.2010, obs.table=obs.table.2010, 
                      truncation=trunc.BSD.10, key="hr", adjustment= "poly")
summary(BSD.df.hr.poly.10)
# poly(2) selected. 

# compare HR models
bsd.hr.comp.10 <- summarize_ds_models(BSD.df.hr.cos.10,BSD.df.hr.poly.10,
                                    output="plain")
bsd.hr.comp.10[ ,1:4]
bsd.hr.comp.10[ ,5:7]
# All CvM p values are above 0.05. All models have some support based on AIC



# plot all the fits
par(mfrow=c(2,2))

# hr cos
plot(BSD.df.hr.cos.10, main = "BSD.df.hr.cos.10")

covar.fit <- ddf.gof(BSD.df.hr.cos.10$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)


# hr poly
plot(BSD.df.hr.poly.10, main = "BSD.df.hr.poly.10")

covar.fit <- ddf.gof(BSD.df.hr.poly.10$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

#  I prefer the cos fit - it doesn't get pulled down by the dip in obs around 10m

# BSD.df.hr.cos.10 is selected



### binned

# hr cos
BSD.df.hr.cos.10.bin <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                           sample.table=sample.table.2010, obs.table=obs.table.2010, 
                           truncation=trunc.BSD.10.bin, key="hr", adjustment= "cos",
                           cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.hr.cos.10.bin)
# no adjustment selected

# hr poly
BSD.df.hr.poly.10.bin <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                           sample.table=sample.table.2010, obs.table=obs.table.2010, 
                           truncation=trunc.BSD.10.bin, key="hr", adjustment= "poly",
                           cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.hr.poly.10.bin)
# no adjustment selected

# plot
plot(BSD.df.hr.cos.10.bin)
# looks like a decent fit to me, and ecologically plausible with such a noisy species that lives in large groups


    # Compare primary models ####

### unbinned

bsd.df.prim.comp.10 <- summarize_ds_models(BSD.df.unif.poly.10, BSD.df.hn.cos.10, BSD.df.hr.cos.10, 
                                           output = "plain")
bsd.df.prim.comp.10[ ,1:5]
bsd.df.prim.comp.10[ ,6:7]
# half normal has the most support, and is the best looking fit. Other two models AIC > 2


### binned

bsd.df.prim.comp.10.bin <- summarize_ds_models(BSD.df.unif.cos.10.bin, BSD.df.hn.herm.10.bin, 
                                           BSD.df.hr.cos.10.bin, 
                                           output = "plain")
bsd.df.prim.comp.10.bin[ ,1:5]
bsd.df.prim.comp.10.bin[ ,6:7]
# all models have some support

# plot
par(mfrow=c(2,2))
plot(BSD.df.unif.cos.10.bin, main="unif")
plot(BSD.df.hn.herm.10.bin, main="hn")
plot(BSD.df.hr.cos.10.bin, main="hr")
# uniform actually looks to be the best fit

# check estimates
bsd.prim.comp.bin <- data.frame(model=c("uni", "hn", "hr"),
                          estimate = c(BSD.df.unif.cos.10.bin$dht$clusters$N$Estimate,
                                       BSD.df.hn.herm.10.bin$dht$clusters$N$Estimate,
                                       BSD.df.hr.cos.10.bin$dht$clusters$N$Estimate),
                          cv = c(BSD.df.unif.cos.10.bin$dht$clusters$N$cv,
                                 BSD.df.hn.herm.10.bin$dht$clusters$N$cv,
                                 BSD.df.hr.cos.10.bin$dht$clusters$N$cv),
                          se = c(BSD.df.unif.cos.10.bin$dht$clusters$N$se,
                                 BSD.df.hn.herm.10.bin$dht$clusters$N$se,
                                 BSD.df.hr.cos.10.bin$dht$clusters$N$se),
                          lcl = c(BSD.df.unif.cos.10.bin$dht$clusters$N$lcl,
                                  BSD.df.hn.herm.10.bin$dht$clusters$N$lcl,
                                  BSD.df.hr.cos.10.bin$dht$clusters$N$lcl),
                          ucl = c(BSD.df.unif.cos.10.bin$dht$clusters$N$ucl,
                                  BSD.df.hn.herm.10.bin$dht$clusters$N$ucl,
                                  BSD.df.hr.cos.10.bin$dht$clusters$N$ucl))

bsd.prim.comp.bin

ggplot(bsd.prim.comp.bin, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,11000)

# not much in it.  I think the choice is between hn and uni, and seeing as I want to test covaras, I will go for hn

    # Models with harsh truncation ####

# set trunc distance 
trunc.harsh.10 <- 35
trunc.harsh.10.bin <- 50


### unbinned

par(mfrow=c(1,2))

# HN cos
BSD.df.hn.cos.harsh.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                    sample.table=sample.table.2010, obs.table=obs.table.2010,
                    truncation=trunc.harsh.10, key="hn", adjustment="cos")
summary(BSD.df.hn.cos.harsh.10)
# no adjustment selected

# compare original with harsh truncation
ddf.gof(BSD.df.hn.cos.harsh.10$ddf, main = "BSD.df.hn.cos.harsh.10", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(BSD.df.hn.cos.10$ddf, main = "BSD.df.hn.cos.10",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# original is better

# I will stick with original half normal model



### binned

BSD.df.hn.herm.10.bin.harsh <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                            sample.table=sample.table.2010, obs.table=obs.table.2010, 
                            truncation=trunc.harsh.10.bin, key="hn",
                            cutpoints = c(0,10,20,30,40,50))
summary(BSD.df.hn.herm.10.bin.harsh)

# plot both
par(mfrow=c(1,2))
plot(BSD.df.hn.herm.10.bin, main="orig")
plot(BSD.df.hn.herm.10.bin.harsh, main="harsh")
# harsh truncation looks to be a better fit close to the line

# check estimates
bsd.harsh.comp.bin <- data.frame(model=c("orig", "harsh"),
                                estimate = c(BSD.df.hn.herm.10.bin$dht$clusters$N$Estimate,
                                             BSD.df.hn.herm.10.bin.harsh$dht$clusters$N$Estimate),
                                cv = c(BSD.df.hn.herm.10.bin$dht$clusters$N$cv,
                                       BSD.df.hn.herm.10.bin.harsh$dht$clusters$N$cv),
                                se = c(BSD.df.hn.herm.10.bin$dht$clusters$N$se,
                                       BSD.df.hn.herm.10.bin.harsh$dht$clusters$N$se),
                                lcl = c(BSD.df.hn.herm.10.bin$dht$clusters$N$lcl,
                                        BSD.df.hn.herm.10.bin.harsh$dht$clusters$N$lcl),
                                ucl = c(BSD.df.hn.herm.10.bin$dht$clusters$N$ucl,
                                        BSD.df.hn.herm.10.bin.harsh$dht$clusters$N$ucl))

bsd.harsh.comp.bin

ggplot(bsd.harsh.comp.bin, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,11000)
# very similar.  Original is a little bit more precise and so is retained

    # Models with covariates ####


### unbinned


## cluster size
BSD.df.hn.size.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                    sample.table=sample.table.2010, obs.table=obs.table.2010,
                    truncation=trunc.BSD.10, key="hn", formula = ~size)
summary(BSD.df.hn.size.10)

# plot
par(mfrow=c(1,2))
plot(BSD.df.hn.size.10, main = "BSD.df.hn.size.10")

covar.fit <- ddf.gof(BSD.df.hn.size.10$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



## Observer
BSD.df.hn.covar.obs.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                      sample.table=sample.table.2010, obs.table=obs.table.2010, 
                      truncation=trunc.BSD.10, key="hn", formula = ~ obs.observer)
summary(BSD.df.hn.covar.obs.10)

# plot
plot(BSD.df.hn.covar.obs.10, main = "BSD.df.hn.covar.obs.10")

covar.fit <- ddf.gof(BSD.df.hn.covar.obs.10$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## Observer + size
BSD.df.hn.covar.obs.size.10 <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                      sample.table=sample.table.2010, obs.table=obs.table.2010, 
                      truncation=trunc.BSD.10, key="hn", formula = ~ obs.observer+size)
summary(BSD.df.hn.covar.obs.size.10)

# plot
plot(BSD.df.hn.covar.obs.size.10, main = "BSD.df.hn.covar.obs.size.10")

covar.fit <- ddf.gof(BSD.df.hn.covar.obs.size.10$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



### binned


# size
BSD.df.hn.size.10.bin <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                            sample.table=sample.table.2010, obs.table=obs.table.2010, 
                            truncation=trunc.BSD.10.bin, key="hn", formula=~size,
                            cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.hn.size.10.bin)
plot(BSD.df.hn.size.10.bin)

# observer
BSD.df.hn.observer.10.bin <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                            sample.table=sample.table.2010, obs.table=obs.table.2010, 
                            truncation=trunc.BSD.10.bin, key="hn", formula=~obs.observer,
                            cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.hn.observer.10.bin)
plot(BSD.df.hn.observer.10.bin)


# observer+size
BSD.df.hn.obs.size.10.bin <- ds(data=BSD.2010.data, region.table=region.table.2010, 
                                sample.table=sample.table.2010, obs.table=obs.table.2010, 
                                truncation=trunc.BSD.10.bin, key="hn", formula=~obs.observer+size,
                                cutpoints = c(0,10,20,30,40,50,60))
summary(BSD.df.hn.obs.size.10.bin)
plot(BSD.df.hn.obs.size.10.bin)



    # Compare covariate models ####


### unbinned

bsd.df.cov.comp.10 <- summarize_ds_models(BSD.df.hn.size.10, BSD.df.hn.covar.obs.10,
                                          BSD.df.hn.covar.obs.size.10,output = "plain")

bsd.df.cov.comp.10[ ,1:5]
bsd.df.cov.comp.10[ ,6:7]

# All models have some support. Based on previous years, observer seems to be important. And in this year there did not appera to be a size bias problem. Therefore I will select the model with observer only



### binned

bsd.df.cov.comp.10.bin <- summarize_ds_models(BSD.df.hn.size.10.bin,BSD.df.hn.observer.10.bin,
                                              BSD.df.hn.obs.size.10.bin,output = "plain")

bsd.df.cov.comp.10.bin[ ,1:5]
bsd.df.cov.comp.10.bin[ ,6:7]
# size has the most support

# check estimates
bsd.cov.comp.bin <- data.frame(model=c("size", "obs", "obs+size"),
                                estimate = c(BSD.df.hn.size.10.bin$dht$clusters$N$Estimate,
                                             BSD.df.hn.observer.10.bin$dht$clusters$N$Estimate,
                                             BSD.df.hn.obs.size.10.bin$dht$clusters$N$Estimate),
                                cv = c(BSD.df.hn.size.10.bin$dht$clusters$N$cv,
                                       BSD.df.hn.observer.10.bin$dht$clusters$N$cv,
                                       BSD.df.hn.obs.size.10.bin$dht$clusters$N$cv),
                                se = c(BSD.df.hn.size.10.bin$dht$clusters$N$se,
                                       BSD.df.hn.observer.10.bin$dht$clusters$N$se,
                                       BSD.df.hn.obs.size.10.bin$dht$clusters$N$se),
                                lcl = c(BSD.df.hn.size.10.bin$dht$clusters$N$lcl,
                                        BSD.df.hn.observer.10.bin$dht$clusters$N$lcl,
                                        BSD.df.hn.obs.size.10.bin$dht$clusters$N$lcl),
                                ucl = c(BSD.df.hn.size.10.bin$dht$clusters$N$ucl,
                                        BSD.df.hn.observer.10.bin$dht$clusters$N$ucl,
                                        BSD.df.hn.obs.size.10.bin$dht$clusters$N$ucl))

bsd.cov.comp.bin

ggplot(bsd.cov.comp.bin, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,11000)
# estimates are all very similar.  size model is the most precise and so is selected


    # Compare primary and covariate models ####

### unbinned

bsd.df.final.compare.10 <- summarize_ds_models(BSD.df.hn.cos.10,BSD.df.hn.covar.obs.10, output = "plain")
bsd.df.final.compare.10[ ,1:3]
bsd.df.final.compare.10[ ,4:7]

# both models have support. I will select the simplest model in this case, as neither appears better than the other.



### binned

bsd.df.final.compare.10.bin <- summarize_ds_models(BSD.df.hn.herm.10.bin, BSD.df.hn.size.10.bin, 
                                                   output = "plain")
bsd.df.final.compare.10.bin[ ,1:3]
bsd.df.final.compare.10.bin[ ,4:7]
# not really anything in it

# check estimates
bsd.final.comp.bin <- data.frame(model=c("orig", "size"),
                               estimate = c(BSD.df.hn.herm.10.bin$dht$clusters$N$Estimate,
                                            BSD.df.hn.size.10.bin$dht$clusters$N$Estimate),
                               cv = c(BSD.df.hn.herm.10.bin$dht$clusters$N$cv,
                                      BSD.df.hn.size.10.bin$dht$clusters$N$cv),
                               se = c(BSD.df.hn.herm.10.bin$dht$clusters$N$se,
                                      BSD.df.hn.size.10.bin$dht$clusters$N$se),
                               lcl = c(BSD.df.hn.herm.10.bin$dht$clusters$N$lcl,
                                       BSD.df.hn.size.10.bin$dht$clusters$N$lcl),
                               ucl = c(BSD.df.hn.herm.10.bin$dht$clusters$N$ucl,
                                       BSD.df.hn.size.10.bin$dht$clusters$N$ucl))

bsd.final.comp.bin

ggplot(bsd.final.comp.bin, aes(x=model, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  ylim(0,11000)
# virtually identical. I wiill go with AIC (and very slightly better precision) and go with original model


  ## Final results 2010 ####

### unbinned

BSD.grp.results.10 <- BSD.df.hn.cos.10$dht$clusters$N
BSD.ind.results.10 <- BSD.df.hn.cos.10$dht$individuals$N

# bind results
BSD.results.2010 <- rbind(BSD.grp.results.10, BSD.ind.results.10)
BSD.results.2010$Label <- as.character(BSD.results.2010$Label)
BSD.results.2010[1,1] <- "Grp"
BSD.results.2010[2,1] <- "Ind"
BSD.results.2010$Year <- "2010"
BSD.results.2010$Species <- "BSD"
BSD.results.2010$DetFun <- "annual"
BSD.results.2010$Key <- "Hn"
BSD.results.2010$Adjust <- "NA"
BSD.results.2010$Covar <- "NA"
BSD.results.2010 <- BSD.results.2010 %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)
BSD.results.2010 <- BSD.results.2010 %>% rename(N = Estimate)

### extract density
BSD.grp.density.10 <- BSD.df.hn.cos.10$dht$clusters$D
BSD.ind.density.10 <- BSD.df.hn.cos.10$dht$individuals$D
BSD.density.10 <- rbind(BSD.grp.density.10,BSD.ind.density.10)
BSD.density.10 <- BSD.density.10[, -c(1,7)]
BSD.density.10 <- BSD.density.10 %>% rename(D = Estimate)

# merge N and D
BSD.results.2010 <- cbind(BSD.results.2010,BSD.density.10)



### binned

# extract estimates
BSD.grp.results.10.bin <- BSD.df.hn.herm.10.bin$dht$clusters$N
BSD.ind.results.10.bin <- BSD.df.hn.herm.10.bin$dht$individuals$N

# bind results
BSD.results.2010.bin <- rbind(BSD.grp.results.10.bin, BSD.ind.results.10.bin)
BSD.results.2010.bin$Label <- as.character(BSD.results.2010.bin$Label)
BSD.results.2010.bin[1,1] <- "Grp"
BSD.results.2010.bin[2,1] <- "Ind"
BSD.results.2010.bin$Year <- "2010"
BSD.results.2010.bin$Species <- "BSD"
BSD.results.2010.bin$DetFun <- "annual"
BSD.results.2010.bin$Key <- "Hn"
BSD.results.2010.bin$Adjust <- "NA"
BSD.results.2010.bin$Covar <- "NA"
BSD.results.2010.bin <- BSD.results.2010.bin %>% 
                        select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)


#### BSD FINAL RESULTS ALL YEARS ####

### If you have re-run all of the above annual analyses, then do this:

BSD.results.final <- rbind(BSD.results.2010,BSD.results.2011.bin,BSD.results.2013.bin,BSD.results.2014,
                           BSD.results.2016,BSD.results.2018,BSD.results.2020)

colnames(BSD.results.final) <- c("Year","Species","DetFun","Key","Adjust","Covar","Label","N","n_se","n_cv",
                                 "n_lcl","n_ucl","D","d_se","d_cv","d_lcl","d_ucl")

BSD.results.final <- BSD.results.final %>% arrange(Label,Year)

BSD.results.final <- BSD.results.final %>% 
                      mutate(D = D*1000000) %>% 
                      mutate(d_se = d_se*1000000) %>% 
                      mutate(d_lcl = d_lcl*1000000) %>% 
                      mutate(d_ucl = d_ucl*1000000)

# save new results
write.csv(BSD.results.final, "Output/Results/BSD_results_final.csv")

### If you are just adding a new year, do the below (provided you have made the annual results table correctly (see the 2020 results section for example). But obviously adjusting the info/data/years
BSD_results <- read.csv("Output/Results/BSD_results_final.csv")
BSD_results <- BSD_results[,-1]

BSD_results <- rbind(BSD_results,BSD.results.2020)
BSD_results <- BSD_results %>% arrange(Label)

# save new results
write.csv(BSD.results, "Output/Results/BSD_results_final.csv")


## Updated results for 2011, 2014 & 2016 (02.09.20) which account for strata.

# load old results for all years
BSD.results <- read.csv("Output/Results/BSD_results_final.csv")
BSD.results <- BSD.results[ ,-1]

# remove old 2011 and 2014 results
BSD.results <- BSD.results[-c(2,4,5,9,11,12), ]

# add new results
BSD.results <- rbind(BSD.results,BSD.results.2011.bin,BSD.results.2014,BSD.results.2016)


# re-order
BSD.results <- BSD.results %>% arrange(Label,Year)



# plot
BSD_final_plot <- ggplot(BSD.results[BSD.results$Label=="Grp",], aes(x=as.factor(Year), y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  ylim(0,13100)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(BSD.results, "Output/Results/BSD_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/BSD_final_plot.png", BSD_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")



### BINNED

### collating results to replace 2010, 2011, 2013, 2016 with binned results

# load in original results
BSD.results.all.orig <- read.csv("Output/Results/BSD_results_final.csv")
BSD.results.all.orig <- BSD.results.all.orig[,-1]
BSD.results.all.orig$Analysis <- "Original"
BSD.results.all.bin <- BSD.results.all.orig 

# replace rows in new main results table with binned results
BSD.results.all.bin[1,] <- BSD.results.2010.bin[BSD.results.2010.bin$Label=="Grp",]
BSD.results.all.bin[8,] <- BSD.results.2010.bin[BSD.results.2010.bin$Label=="Ind",]
BSD.results.all.bin[2,] <- BSD.results.2011.bin[BSD.results.2011.bin$Label=="Grp",]
BSD.results.all.bin[9,] <- BSD.results.2011.bin[BSD.results.2011.bin$Label=="Ind",]
BSD.results.all.bin[3,] <- BSD.results.2013.bin[BSD.results.2013.bin$Label=="Grp",]
BSD.results.all.bin[10,] <- BSD.results.2013.bin[BSD.results.2013.bin$Label=="Ind",]
BSD.results.all.bin[5,] <- BSD.results.2016.bin[BSD.results.2016.bin$Label=="Grp",]
BSD.results.all.bin[12,] <- BSD.results.2016.bin[BSD.results.2016.bin$Label=="Ind",]
BSD.results.all.bin$Analysis <- "Binned"

## save the binned results
write.csv(BSD.results.all.bin, file="Output/Results/BSD_results_final_binned.csv")

BSD.results.all.bin <- read.csv("Output/Results/BSD_results_final_binned.csv")
BSD.results.all.bin <- BSD.results.all.bin[,-1]

## plot the original results with the binned results for comparison

# bind the two tables
BSD.results.all.comp <- rbind(BSD.results.all.orig,BSD.results.all.bin)

# plot together 
BSD_results_plot_comparison <- ggplot(BSD.results.all.comp[BSD.results.all.comp$Label=="Grp",], 
                               aes(x=Year, y=Estimate, group=Analysis, colour=Analysis))+
                          geom_point(position = position_dodge(width=0.5),shape=16, size=2)+
                          geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.5),width=0.3)+
                          ylim(0,16000)+
                          theme_bw()+
                          theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),panel.border = element_blank())+
                          theme(axis.line = element_line(color = 'black'))

ggsave("Output/Results/Plots/Point_estimates/BSD_comparison_binned.png", BSD_results_plot_comparison,
       dpi=300, width = 20, height = 20, units = "cm")

# plot binned results only
BSD_results_plot_binned <- ggplot(BSD.results.all.bin[BSD.results.all.bin$Label=="Grp",], 
                                  aes(x=Year, y=Estimate))+
                            geom_point(shape=16, size=2)+
                            geom_errorbar(aes(ymin=lcl, ymax=ucl, width=0.3))+
                            ylim(0,16000)+
                            theme_bw()+
                            theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                            theme(axis.line = element_line(color = 'black'))

ggsave("Output/Results/Plots/Point_estimates/BSD_final_plot_binned.png", BSD_results_plot_binned,
       dpi=300, width = 20, height = 20, units = "cm")


### Combining results to form the final set (see explanation under "binning" at the top of this script)

# extract 2011 and 2013 binned results
binned.grp <- BSD.results.all.bin[2:3,]
binned.ind <- BSD.results.all.bin[9:10,]

# slot them into a new dataframe with the rest of the original results
BSD.results.all.comb <- BSD.results.all.orig
BSD.results.all.comb[2:3,] <- binned.grp
BSD.results.all.comb[9:10,] <- binned.ind

# save the results
write.csv(BSD.results.all.comb, file="Output/Results/BSD_results_final_combined.csv")

# plot coloured by analysis type
ggplot(BSD.results.all.comb[BSD.results.all.comb$Label=="Grp",], aes(x=Year,y=Estimate, colour=Analysis))+
  geom_point(shape=16, size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl, width=0.3))+
  ylim(0,16000)+
  theme_bw()+
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank())+
  theme(axis.line = element_line(color = 'black'))

# final version to save
BSD_final_plot_combined <- ggplot(BSD.results.all.comb[BSD.results.all.comb$Label=="Grp",], 
                                  aes(x=as.factor(Year),y=Estimate))+
                          geom_point(shape=16, size=2)+
                          geom_errorbar(aes(ymin=lcl, ymax=ucl, width=0.3))+
                          ylim(0,16000)+
                          xlab("Year")+
                          ylab("Estimated group abundance")+
                          theme_bw()+
                          theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),panel.border = element_blank())+
                          theme(axis.line = element_line(color = 'black'))

ggsave("Output/Results/Plots/Point_estimates/BSD_final_plot_combined.png", BSD_final_plot_combined,
       dpi=300, width = 20, height = 20, units = "cm")


#### Yellow-cheeked crested gibbon ############################################
  ## Subset data ####

### Previously, effort strata were ignored for this species, as they were for the other species. After discussions with Olly, we decided to check the effect of accounting for the effort strata for YCG in 2011, 2014, and 2016, as the estimates for these years were very high compared to other years, suggesting the strata were having an effect. I tested a couple of different approaches (see "Analytical_approach_test_CDS.R", lines 2190 onwards).  For YCG, the approach we have settled on is to split the stratified and un-stratified years, and analyse them separately. This accounts for the effect of the strata, and removes any effects in the unstratified years.

# This means that the old data will be loaded and used for 2011, 2014, and 2016, but the new data will be used for all of the unstratified years.  

# note - I still need to remove T20, change the survey area, and add the stratum variable


# load old data (2010-2018) that have strata
load("./Output/Data/Archive/KSWS_MASTER.Rdata")

allData$obs.habitat <- as.factor(allData$obs.habitat)
allData$obs.observer <- as.factor(allData$obs.observer)

## Remove T20
allData <- allData[allData$obs.transect != 20, ]
sample.table <- sample.table[sample.table$Sample.Label != 20,]
obs.table <- obs.table[obs.table$Sample.Label != 20,]


# subset YCG data
YCG.data <- allData[allData$species=="YCG",] 


# subset data for 2011, 2014, 2016
YCG.data.strat <- as.data.frame(YCG.data[YCG.data$stratum=="2011_H" | YCG.data$stratum=="2011_L" 
                                         | YCG.data$stratum=="2014_H" | YCG.data$stratum=="2014_L"
                                         | YCG.data$stratum=="2016_H" | YCG.data$stratum=="2016_L", ]) 


region.table.strat <- as.data.frame(full.region.table[full.region.table$Region.Label=="2011_H"  |
                                                        full.region.table$Region.Label=="2011_L" |
                                                        full.region.table$Region.Label=="2014_H" |
                                                        full.region.table$Region.Label=="2014_L" |
                                                        full.region.table$Region.Label=="2016_H" |
                                                        full.region.table$Region.Label=="2016_L",])

# change survey area
new.area <- 1880000000/2
region.table.strat$Area <- new.area


# create region tables where other strata are set to 0
region.table.11 <- region.table.strat
region.table.11$Area[region.table.11$Region.Label != "2011_L" & region.table.11$Region.Label != "2011_H"] <- 0

region.table.14 <- region.table.strat
region.table.14$Area[region.table.14$Region.Label != "2014_L" & region.table.14$Region.Label != "2014_H"] <- 0

region.table.16 <- region.table.strat
region.table.16$Area[region.table.16$Region.Label != "2016_L" & region.table.16$Region.Label != "2016_H"] <- 0


sample.table.strat <- as.data.frame(sample.table[sample.table$Region.Label =="2011_H" | 
                                                 sample.table$Region.Label =="2011_L" |
                                                 sample.table$Region.Label =="2014_H" |
                                                 sample.table$Region.Label =="2014_L" |
                                                 sample.table$Region.Label =="2016_H" |
                                                 sample.table$Region.Label =="2016_L",])

obs.table.strat <- as.data.frame(obs.table[obs.table$Region.Label=="2011_H" | obs.table$Region.Label=="2011_L" |
                                             obs.table$Region.Label=="2014_H" | obs.table$Region.Label=="2014_L" |
                                             obs.table$Region.Label=="2016_H" | obs.table$Region.Label=="2016_L",])



### now using the updated data, I will get the 2010, 2013, 2018, 2020 data

# load the latest data
load("./Output/Data/KSWS_MASTER.Rdata")

## Remove T20
allData <- allData[allData$obs.transect != 20, ]
sample.table <- sample.table[sample.table$Sample.Label != 20,]
obs.table <- obs.table[obs.table$Sample.Label != 20,]


# subset YCG data
YCG.data <- as.data.frame(allData[allData$species=="YCG",]) 


YCG.data.nostrat <- YCG.data[YCG.data$stratum==2010 | YCG.data$stratum==2013 | YCG.data$stratum==2018 |
                               YCG.data$stratum==2020, ]

region.table.nostrat <- as.data.frame(full.region.table[full.region.table$Region.Label=="2010" | 
                                                          full.region.table$Region.Label=="2013" |
                                                          full.region.table$Region.Label=="2018" |
                                                          full.region.table$Region.Label=="2020", ])

# change survey area
new.area <- 1880000000
region.table.nostrat$Area <- new.area


sample.table.nostrat <- as.data.frame(
        sample.table[sample.table$Region.Label==2010 | sample.table$Region.Label==2013 |
        sample.table$Region.Label==2018 | sample.table$Region.Label==2020, ])

obs.table.nostrat <- as.data.frame(obs.table[obs.table$Region.Label==2010 | obs.table$Region.Label==2013 |
                                               obs.table$Region.Label==2018 | obs.table$Region.Label==2020, ])


# create scaled continuous stratum variable for the detection functions
YCG.data.nostrat$stratum <- as.vector(scale(YCG.data.nostrat$stratum, center = T, scale = T)) 
YCG.data.strat$year <- ifelse(YCG.data.strat$stratum=="2011_L"|YCG.data.strat$stratum=="2011_H",2011,
                          ifelse(YCG.data.strat$stratum=="2014_L"|YCG.data.strat$stratum=="2014_H",2014,
                              ifelse(YCG.data.strat$stratum=="2016_L"|YCG.data.strat$stratum=="2016_H",2016,NA)))
YCG.data.strat$stratum <- as.vector(scale(YCG.data.strat$year, center = T, scale = T)) 

  ## New analysis - splitting years by strata ####

# this analysis is from Sept 2020 and is aimed at accounting for strata in 2011, 2014, and 2016.  For the original analysis, go to the sections below starting at "Exploritory plots and linear models"

    # Stratified years ####

par(mfrow=c(1,2))
# histograms of the data
hist(YCG.data.strat$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(YCG.data.strat$distance, main=NULL, xlab="Distance (m)", breaks=c(40))

# set truncation distance
ycg.strat.trunc <- 60


### run DF models

## Primary models

# Uni cosine
YCG.df.uni.cos <- ds(data=YCG.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ycg.strat.trunc, key="uni", adjustment= "cos")

# Uni poly
YCG.df.uni.poly <- ds(data=YCG.data.strat, region.table=region.table.strat, 
                      sample.table=sample.table.strat, obs.table=obs.table.strat, 
                      truncation=ycg.strat.trunc, key="uni", adjustment= "poly")

# HN cosine
YCG.df.hn.cos <- ds(data=YCG.data.strat, region.table=region.table.strat, 
                    sample.table=sample.table.strat, obs.table=obs.table.strat, 
                    truncation=ycg.strat.trunc, key="hn", adjustment= "cos")


# HN hermite
YCG.df.hn.herm <- ds(data=YCG.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ycg.strat.trunc, key="hn", adjustment= "herm")


# HR cosine
YCG.df.hr.cos <- ds(data=YCG.data.strat, region.table=region.table.strat, 
                    sample.table=sample.table.strat, obs.table=obs.table.strat, 
                    truncation=ycg.strat.trunc, key="hr", adjustment= "cos")

# HR poly
YCG.df.hr.poly <- ds(data=YCG.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ycg.strat.trunc, key="hr", adjustment= "poly")

# compare models
ycg.df.prim.comp <- summarize_ds_models(YCG.df.uni.cos, YCG.df.uni.poly,
                                        YCG.df.hn.cos,  
                                        YCG.df.hr.cos,  
                                        output = "plain")
ycg.df.prim.comp[ ,1:5]
ycg.df.prim.comp[ ,6:7]
# all models have some support

# plot all
par(mfrow=c(2,2))

plot(YCG.df.uni.cos, main="YCG.df.uni.cos")
plot(YCG.df.uni.poly, main="YCG.df.uni.poly")
plot(YCG.df.hn.cos, main="YCG.df.hn.cos")
plot(YCG.df.hr.cos, main="YCG.df.hr.cos")

# hn cos selected


## covariate models

# cluster size
YCG.df.hn.size <- ds(data=YCG.data.strat, region.table=region.table.strat, 
                    sample.table=sample.table.strat, obs.table=obs.table.strat, 
                    truncation=ycg.strat.trunc, key="hn", formula=~size)

# observer
YCG.df.hn.observer <- ds(data=YCG.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ycg.strat.trunc, key="hn", formula=~obs.observer)

# stratum
YCG.df.hn.stratum <- ds(data=YCG.data.strat, region.table=region.table.strat, 
                         sample.table=sample.table.strat, obs.table=obs.table.strat, 
                         truncation=ycg.strat.trunc, key="hn", formula=~stratum)

# compare models
ycg.df.cov.comp <- summarize_ds_models(YCG.df.hn.size, YCG.df.hn.observer, YCG.df.hn.stratum,
                                        output = "plain")
ycg.df.cov.comp[ ,1:5]
ycg.df.cov.comp[ ,6:7]
# all models have some support


# extract and plot to see what the different models say
size <- data.frame(label = "size",
                   estimate = YCG.df.hn.size$dht$clusters$N$Estimate[c(1:6)],
                   lcl = YCG.df.hn.size$dht$clusters$N$lcl[c(1:6)],
                   ucl = YCG.df.hn.size$dht$clusters$N$ucl[c(1:6)])

observer <- data.frame(label = "strat",
                    estimate = YCG.df.hn.observer$dht$clusters$N$Estimate[c(1:6)],
                    lcl = YCG.df.hn.observer$dht$clusters$N$lcl[c(1:6)],
                    ucl = YCG.df.hn.observer$dht$clusters$N$ucl[c(1:6)])

strat <- data.frame(label = "strat.size",
                         estimate = YCG.df.hn.stratum$dht$clusters$N$Estimate[c(1:6)],
                         lcl = YCG.df.hn.stratum$dht$clusters$N$lcl[c(1:6)],
                         ucl = YCG.df.hn.stratum$dht$clusters$N$ucl[c(1:6)])

mod.check <- rbind(size,observer,strat)
mod.check$region <- rep(c("11_H","11_L","14_H","14_L","16_H","16_L"), times=3)

ggplot(mod.check, aes(x=region, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# size model has generally lower CVs and SEs

# size selected


## compare primary with covariate model
ycg.df.final.comp <- summarize_ds_models(YCG.df.hn.size, YCG.df.hn.cos,
                                       output = "plain")
ycg.df.final.comp[ ,1:5]
ycg.df.final.comp[ ,6:7]
# both models have support


# plot the estimates
prim <- data.frame(label = "primary",
                   estimate = YCG.df.hn.cos$dht$clusters$N$Estimate[c(1:6)],
                   lcl = YCG.df.hn.cos$dht$clusters$N$lcl[c(1:6)],
                   ucl = YCG.df.hn.cos$dht$clusters$N$ucl[c(1:6)])

final.mod.check <- rbind(size,prim)
final.mod.check$region <- rep(c("11_H","11_L","14_H","14_L","16_H","16_L"), times=2)

ggplot(final.mod.check, aes(x=region, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))
# very similar results. Not much difference in precision. I will select size model because group size is likely to be influential in the detection process in reality


## run final models to get estimates for each year (areas set to 0)

# 2011
YCG.11.hn.size <- ds(data=YCG.data.strat, region.table=region.table.11, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ycg.strat.trunc, key="hn", formula=~size)

# 2014
YCG.14.hn.size <- ds(data=YCG.data.strat, region.table=region.table.14, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ycg.strat.trunc, key="hn", formula=~size)

# 2016
YCG.16.hn.size <- ds(data=YCG.data.strat, region.table=region.table.16, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ycg.strat.trunc, key="hn", formula=~size)


    # Unstratified years ####

par(mfrow=c(1,2))
# histograms of the data
hist(YCG.data.nostrat$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(YCG.data.nostrat$distance, main=NULL, xlab="Distance (m)", breaks=c(40))

# set truncation distance
ycg.nostrat.trunc <- 60


### run DF models

## Primary models

# Uni cosine
YCG.df.uni.cos <- ds(data=YCG.data.nostrat, region.table=region.table.nostrat, 
                     sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                     truncation=ycg.nostrat.trunc, key="uni", adjustment= "cos")

# Uni poly
YCG.df.uni.poly <- ds(data=YCG.data.nostrat, region.table=region.table.nostrat, 
                      sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                      truncation=ycg.nostrat.trunc, key="uni", adjustment= "poly")

# HN cosine
YCG.df.hn.cos <- ds(data=YCG.data.nostrat, region.table=region.table.nostrat, 
                    sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                    truncation=ycg.nostrat.trunc, key="hn", adjustment= "cos")


# HN hermite
YCG.df.hn.herm <- ds(data=YCG.data.nostrat, region.table=region.table.nostrat, 
                     sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                     truncation=ycg.nostrat.trunc, key="hn", adjustment= "herm")


# HR cosine
YCG.df.hr.cos <- ds(data=YCG.data.nostrat, region.table=region.table.nostrat, 
                    sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                    truncation=ycg.nostrat.trunc, key="hr", adjustment= "cos")

# HR poly
YCG.df.hr.poly <- ds(data=YCG.data.nostrat, region.table=region.table.nostrat, 
                     sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                     truncation=ycg.nostrat.trunc, key="hr", adjustment= "poly")

# compare models
ycg.df.prim.comp <- summarize_ds_models(YCG.df.uni.cos, YCG.df.uni.poly,
                                        YCG.df.hn.cos,  
                                        YCG.df.hr.cos,  
                                        output = "plain")
ycg.df.prim.comp[ ,1:5]
ycg.df.prim.comp[ ,6:7]
# HN and Uni Cos have most support

# plot
plot(YCG.df.uni.cos, main="YCG.df.uni.cos")
plot(YCG.df.hn.cos, main="YCG.df.hn.cos")
# very very similar. I will go with HN as I want to test covariates


## covariate models 

# cluster size
YCG.df.hn.size <- ds(data=YCG.data.nostrat, region.table=region.table.nostrat, 
                     sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                     truncation=ycg.nostrat.trunc, key="hn", formula=~size)

# observer
YCG.df.hn.observer <- ds(data=YCG.data.nostrat, region.table=region.table.nostrat, 
                     sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                     truncation=ycg.nostrat.trunc, key="hn", formula=~obs.observer)

# stratum
YCG.df.hn.stratum <- ds(data=YCG.data.nostrat, region.table=region.table.nostrat, 
                         sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                         truncation=ycg.nostrat.trunc, key="hn", formula=~stratum)


# compare models
ycg.df.cov.comp <- summarize_ds_models(YCG.df.hn.size, YCG.df.hn.observer, YCG.df.hn.stratum,
                                        output = "plain")
ycg.df.cov.comp[ ,1:5]
ycg.df.cov.comp[ ,6:7]
# stratum model has overwhelming support


# compare primary and covariate models
ycg.df.final.comp <- summarize_ds_models(YCG.df.hn.stratum, YCG.df.hn.cos,
                                       output = "plain")
ycg.df.final.comp[ ,1:5]
ycg.df.final.comp[ ,6:7]

# stratum model has most support and is selected

# skip down to "YCG final results" section for results

#
  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

# distance histograms
hist(YCG.data$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(YCG.data$distance, main=NULL, xlab="Distance (m)", breaks=c(40))

# horrible spike at 40-50m

# Save the chosen truncation distance for later use. I'll start large, and then try the harsher truncation later
trunc.YCG <- 60 

# Count the number of observations discarded
nrow(YCG.data[YCG.data$distance>trunc.YCG,]) 

length(YCG.data$distance)

14/206*100 # 6%

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(1,2))
plot(YCG.data$size, YCG.data$distance, main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.YCG <- lm(distance~size, data=YCG.data)
lines(YCG.data$size, as.vector(predict(lm.YCG, YCG.data)))
summary(lm.YCG)
# no significant relationship between group size and distance.  

# Plot of Observer factor against distance
plot(YCG.data$obs.observer,YCG.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# a lot of variation between observers

# I can't include habitat or AMPM in the DF because those variables weren't recorded in earlier years


  ## Fit a detection function ####
    # Uniform ####

# Uni cosine
YCG.df.uni.cos <- ds(data=YCG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.YCG, key="uni", adjustment= "cos")
summary(YCG.df.uni.cos)
# cosine(1)

# Uni poly
YCG.df.uni.poly <- ds(data=YCG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.YCG, key="uni", adjustment= "poly")
summary(YCG.df.uni.poly)
# poly(2,4)

## compare uni models
ycg.uni.comp <- summarize_ds_models(YCG.df.uni.cos,YCG.df.uni.poly,
                                    output="plain")
ycg.uni.comp[ ,1:4]
ycg.uni.comp[ ,5:7]
# uni cos has most support (uni poly dAIC > 2)

# plot the fits
par(mfrow=c(2,2))

# cos
plot(YCG.df.uni.cos, main = "YCG.df.uni.cos")

covar.fit <- ddf.gof(YCG.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(YCG.df.uni.poly, main = "YCG.df.uni.poly")

covar.fit <- ddf.gof(YCG.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# very similar but cos qq plot is better

# YCG.df.uni.cos selected


    # Half normal ####

# HN cosine
YCG.df.hn.cos <- ds(data=YCG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.YCG, key="hn", adjustment= "cos")
summary(YCG.df.hn.cos)
# cosine(2)

# HN hermite
YCG.df.hn.herm <- ds(data=YCG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.YCG, key="hn", adjustment= "herm")
summary(YCG.df.hn.herm)
# no adjustment selected 

## compare uni models
ycg.hn.comp <- summarize_ds_models(YCG.df.hn.cos,YCG.df.hn.herm,
                                    output="plain")
ycg.hn.comp[ ,1:4]
ycg.hn.comp[ ,5:7]
# very similar

# plot the fits
par(mfrow=c(2,2))

# hn cos
plot(YCG.df.hn.cos, main = "YCG.df.hn.cos")

covar.fit <- ddf.gof(YCG.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# hn herm
plot(YCG.df.hn.herm, main = "YCG.df.hn.herm")

covar.fit <- ddf.gof(YCG.df.hn.herm$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# cos qq plot looks better. The cos fit drops quicker beyond 10m, whereas the herm fit is a slower decrease. There is not much in it though - p reaches 0.5 just before 30m for cos, and just after 30m for herm.  Based on the quicker decrease 9which I think is more realistic), and the qq plot, I will go with cos

# YCG.df.hn.cos selected

    # Hazard rate ####

# HR cosine
YCG.df.hr.cos <- ds(data=YCG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.YCG, key="hr", adjustment= "cos")
summary(YCG.df.hr.cos)
# no adjustment selected


# HR poly
YCG.df.hr.poly <- ds(data=YCG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.YCG, key="hr", adjustment= "poly")
summary(YCG.df.hr.poly)
# no adjustment selected 

# plot the fit
par(mfrow=c(1,2))

# hr cos
plot(YCG.df.hr.cos, main = "YCG.df.hr.cos")

covar.fit <- ddf.gof(YCG.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# shoulder up to just before 10m then steep drop, although p=0.5 roughly the same time as the HN model above. the HN and HR models both have p=~0.2 at 60m.



    # Compare primary models ####

ycg.df.prim.comp <- summarize_ds_models(YCG.df.uni.cos, YCG.df.hn.cos, YCG.df.hr.cos, 
                                           output = "plain")
ycg.df.prim.comp[ ,1:5]
ycg.df.prim.comp[ ,6:7]

# half normal and uniform have the most support (dAIC <2)

# plot the fits together
par(mfrow=c(2,2))

# Uni cos
plot(YCG.df.uni.cos, main = "YCG.df.uni.cos")

covar.fit <- ddf.gof(YCG.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HN cos
plot(YCG.df.hn.cos, main = "YCG.df.hn.cos")

covar.fit <- ddf.gof(YCG.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# basically the HN model has p dropping off quicker beyond ~7m, whereas the uni model has a flatter shoulder. p=0.5 before 30m for HN and ager 30m for UNi.  I think realistically the HN model is closer to reality. The QQ plot also looks better

# HN model selected

    # Models with harsh truncation ####

# first I will test harsher truncation, then I will test a more lenient truncation

# set trunc distance 
ycg.trunc.harsh <- 40

par(mfrow=c(1,2))

# HN cos
YCG.df.hn.cos.harsh <- ds(data=YCG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=ycg.trunc.harsh, key="hn", adjustment="cos")
summary(YCG.df.hn.cos.harsh)
# no adjustment selected

# compare original with harsh truncation
ddf.gof(YCG.df.hn.cos.harsh$ddf, main = "ycg.df.hn.cos.harsh", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(YCG.df.hn.cos$ddf, main = "ycg.df.hn.cos",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# original is better

# YCG.df.hn.cos selected

    # Models with covariates ####

## cluster size
YCG.df.hn.size <- ds(data=YCG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.YCG, key="hn", formula = ~size)
summary(YCG.df.hn.size)

# plot
par(mfrow=c(1,2))
plot(YCG.df.hn.size, main = "YCG.df.hn.size")

covar.fit <- ddf.gof(YCG.df.hn.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



## Observer
YCG.df.hn.obs <- ds(data=YCG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.YCG, key="hn", formula = ~obs.observer)
summary(YCG.df.hn.obs)

# plot
plot(YCG.df.hn.obs, main = "YCG.df.hn.obs")

covar.fit <- ddf.gof(YCG.df.hn.obs$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# lot of variation

## stratum (year)
YCG.df.hn.strat <- ds(data=YCG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.YCG, key="hn", formula = ~stratum)
summary(YCG.df.hn.strat)

# plot
plot(YCG.df.hn.strat, main = "YCG.df.hn.strat")

covar.fit <- ddf.gof(YCG.df.hn.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# looks like there are differences between years, but not sure if they're significant

    # Compare covariate models ####

ycg.df.cov.comp <- summarize_ds_models(YCG.df.hn.size, YCG.df.hn.obs,
                                          YCG.df.hn.strat,output = "plain")

ycg.df.cov.comp[ ,1:5]
ycg.df.cov.comp[ ,6:7]

# model with year as continuous covar is the best model, but the model with size has some support

summary(YCG.df.hn.strat)
# the coefficient for stratum is 0.13, which means for every subsequent year, average p increases by 0.13. This is of note I think.  This suggests that probability of detection is increasing by a non-trivial amount each year. This could be reflecting the skill of the observers increasing over time. I suppose in theory it could suggest habituation of gibbons to people/observers, but I doubt it. 

# I will just check a model with both stratum and size
YCG.df.hn.strat.size <- ds(data=YCG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.YCG, key="hn", formula = ~stratum+size)
summary(YCG.df.hn.strat.size)

# add to the comparison table
ycg.df.cov.comp <- summarize_ds_models(YCG.df.hn.size, YCG.df.hn.obs,YCG.df.hn.strat.size,
                                          YCG.df.hn.strat,output = "plain")

ycg.df.cov.comp[ ,1:5]
ycg.df.cov.comp[ ,6:7]

# all three models hav some support. I will have a closer look at the estimates
summary(YCG.df.hn.size) # cv between 0.2 - 0.35
summary(YCG.df.hn.strat) # cv beween 
summary(YCG.df.hn.strat.size)

# extract and plot to see what the different models say
size <- data.frame(label = "size",
                   estimate = YCG.df.hn.size$dht$clusters$N$Estimate[c(1:7)],
                   lcl = YCG.df.hn.size$dht$clusters$N$lcl[c(1:7)],
                   ucl = YCG.df.hn.size$dht$clusters$N$ucl[c(1:7)])

strat <- data.frame(label = "strat",
                   estimate = YCG.df.hn.strat$dht$clusters$N$Estimate[c(1:7)],
                   lcl = YCG.df.hn.strat$dht$clusters$N$lcl[c(1:7)],
                   ucl = YCG.df.hn.strat$dht$clusters$N$ucl[c(1:7)])

strat.size <- data.frame(label = "strat.size",
                   estimate = YCG.df.hn.strat.size$dht$clusters$N$Estimate[c(1:7)],
                   lcl = YCG.df.hn.strat.size$dht$clusters$N$lcl[c(1:7)],
                   ucl = YCG.df.hn.strat.size$dht$clusters$N$ucl[c(1:7)])

mod.check <- rbind(size,strat,strat.size)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=3)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# ok so our inferences regarding trend will not change based on which model I choose (which is good).  The size and size.strat models fluctuate slightly more, with each of those models being either the highest or lowest estiamte in each year.  The strat only model tends to sit in the middle fo the other two, which I think makes it a reasonable compromise. I also don't trust the size model as the linear model before suggested there wasn't really a relationship between distance and size. I also wonder whether having two covariates in for a species with so few observations is a bit dodge

# YCG.df.hn.strat selected

    # Compare primary and covariate models ####

ycg.df.final.compare <- summarize_ds_models(YCG.df.hn.cos,YCG.df.hn.strat, output = "plain")
ycg.df.final.compare[ ,1:3]
ycg.df.final.compare[ ,4:7]

# not much in it really. compare CV 

YCG.df.hn.cos$dht$clusters$N$cv
YCG.df.hn.strat$dht$clusters$N$cv

# YCG.df.hn.strat has lower CVs and so is selected

  ## YCG Final Results ####

## These first results are from the latest analysis (Sept 2020) which split the years nto stratified and unstratified.  Below these results I have kept in the code for the original results. 

## extract estimates

# 2011 abundance
YCG.grp.results.11 <- YCG.11.hn.size$dht$clusters$N[7,]
YCG.grp.results.11$Label <- "Grp"
YCG.grp.results.11$Year <- 2011

YCG.ind.results.11 <- YCG.11.hn.size$dht$individuals$N[7,]
YCG.ind.results.11$Label <- "Ind"
YCG.ind.results.11$Year <- 2011

# 2011 density
YCG.grp.density.11 <- YCG.11.hn.size$dht$clusters$D[7,]
YCG.grp.density.11$Label <- "Grp"
YCG.grp.density.11$Year <- "2011"

YCG.ind.density.11 <- YCG.11.hn.size$dht$individuals$D[7,]
YCG.ind.density.11$Label <- "Ind"
YCG.ind.density.11$Year <- "2011"

# 2014 abundance
YCG.grp.results.14 <- YCG.14.hn.size$dht$clusters$N[7,]
YCG.grp.results.14$Label <- "Grp"
YCG.grp.results.14$Year <- 2014

YCG.ind.results.14 <- YCG.14.hn.size$dht$individuals$N[7,]
YCG.ind.results.14$Label <- "Ind"
YCG.ind.results.14$Year <- 2014

# 2014 density
YCG.grp.density.14 <- YCG.14.hn.size$dht$clusters$D[7,]
YCG.grp.density.14$Label <- "Grp"
YCG.grp.density.14$Year <- "2014"

YCG.ind.density.14 <- YCG.14.hn.size$dht$individuals$D[7,]
YCG.ind.density.14$Label <- "Ind"
YCG.ind.density.14$Year <- "2014"


# 2016 abundance
YCG.grp.results.16 <- YCG.16.hn.size$dht$clusters$N[7,]
YCG.grp.results.16$Label <- "Grp"
YCG.grp.results.16$Year <- 2016

YCG.ind.results.16 <- YCG.16.hn.size$dht$individuals$N[7,]
YCG.ind.results.16$Label <- "Ind"
YCG.ind.results.16$Year <- 2016

# 2016 density
YCG.grp.density.16 <- YCG.16.hn.size$dht$clusters$D[7,]
YCG.grp.density.16$Label <- "Grp"
YCG.grp.density.16$Year <- "2016"

YCG.ind.density.16 <- YCG.16.hn.size$dht$individuals$D[7,]
YCG.ind.density.16$Label <- "Ind"
YCG.ind.density.16$Year <- "2016"


# no strata abundance
YCG.grp.results.nostrat <- YCG.df.hn.stratum$dht$clusters$N[1:4,]
YCG.grp.results.nostrat <- YCG.grp.results.nostrat %>% rename(Year = Label) %>% 
  mutate(Label = "Grp") %>% 
  select(Label,Estimate,se,cv,lcl,ucl,df,Year)

YCG.ind.results.nostrat <- YCG.df.hn.stratum$dht$individuals$N[1:4,]
YCG.ind.results.nostrat <- YCG.ind.results.nostrat %>% rename(Year = Label) %>% 
  mutate(Label = "Ind") %>% 
  select(Label,Estimate,se,cv,lcl,ucl,df,Year)

# no strata density
YCG.grp.density.nostrat <- YCG.df.hn.stratum$dht$clusters$D[1:4,]
YCG.grp.density.nostrat <- YCG.grp.density.nostrat %>% rename(Year = Label) %>% 
  mutate(Label = "Grp") %>% 
  select(Label,Estimate,se,cv,lcl,ucl,df,Year)

YCG.ind.density.nostrat <- YCG.df.hn.stratum$dht$individuals$D[1:4,]
YCG.ind.density.nostrat <- YCG.ind.density.nostrat %>% rename(Year = Label) %>% 
  mutate(Label = "Ind") %>% 
  select(Label,Estimate,se,cv,lcl,ucl,df,Year)

# bind all abundance together
YCG.results <- rbind(YCG.grp.results.11,YCG.ind.results.11,YCG.grp.results.14,YCG.ind.results.14,
                     YCG.grp.results.16,YCG.ind.results.16,YCG.grp.results.nostrat,
                     YCG.ind.results.nostrat)
YCG.results <- YCG.results[, -7]
colnames(YCG.results) <- c("Label","N","n_se","n_cv","n_lcl","n_ucl","Year")
YCG.results <- YCG.results %>% arrange(Label,Year)


# bind all density together
YCG.density <- rbind(YCG.grp.density.11,YCG.ind.density.11,YCG.grp.density.14,YCG.ind.density.14,
                     YCG.grp.density.16,YCG.ind.density.16,YCG.grp.density.nostrat,YCG.ind.density.nostrat)
YCG.density <- YCG.density[ ,-7]
colnames(YCG.density) <- c("Label","D","d_se","d_cv","d_lcl","d_ucl","Year")
YCG.density <- YCG.density %>% arrange(Label,Year)
YCG.density <- YCG.density[,-7]

# merge N and D
YCG.results.final <- cbind(YCG.results,YCG.density)
YCG.results.final <- YCG.results.final[,-8]

YCG.results.final$Species <- "YCG"

# add other fields and re-order
YCG.results.final <- YCG.results.final %>% 
                      mutate(DetFun = rep("pooled", times=14)) %>% 
                      mutate(Key = rep("Hn", times=14)) %>% 
                      mutate(Adjust = rep("NA", times=14)) %>% 
                      mutate(Covar = rep("year", times=14)) %>% 
                      select(Year,Species,DetFun,Key,Adjust,Covar,Label,N,n_se,n_cv,n_lcl,n_ucl,
                             D,d_se,d_cv,d_lcl,d_ucl) 

YCG.results.final$Covar[YCG.results.final$Year==2011 | YCG.results.final$Year==2014 | 
                    YCG.results.final$Year==2016] <- "Cluster size"

YCG.results.final <- YCG.results.final %>% 
                      mutate(D = D*1000000) %>% 
                      mutate(d_se = d_se*1000000) %>% 
                      mutate(d_lcl = d_lcl*1000000) %>% 
                      mutate(d_ucl = d_ucl*1000000)




YCG_final_plot <- ggplot(YCG.results[YCG.results$Label=="Grp",], aes(x=Year, y= Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  ylim(0,1510)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(YCG.results.final, "Output/Results/YCG_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/YCG_final_plot.png", YCG_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")





### Original results (all data pooled, no strata)

summary(YCG.df.hn.strat)

# extract estimates
YCG.grp.results <- YCG.df.hn.strat$dht$clusters$N[1:7, ]
YCG.ind.results <- YCG.df.hn.strat$dht$individuals$N[1:7, ]

# bind results
YCG.results <- rbind(YCG.grp.results, YCG.ind.results)
YCG.results <- YCG.results %>% rename(Year = Label) %>% 
                mutate(Label = rep(c("Grp", "Ind"), each=7)) %>% 
                mutate(Species = rep("YCG", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hn", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("year", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)



YCG_final_plot <- ggplot(YCG.results[YCG.results$Label=="Grp",], aes(x=Year, y= Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  ylim(0,1500)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(YCG.results, "Output/Results/YCG_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/YCG_final_plot.png", YCG_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")


#### Germain's silver langur #################################################
  ## Subset data ####

# subset GSL data
GSL.data <- allData[allData$species=="GSL",] 
GSL.data <- as.data.frame(GSL.data)
head(GSL.data)

# Total number of groups from all years
length(GSL.data$distance) 

# for GSL, I will pool all years for DF, but test year as a continous covariate in the DF model

# check 2020 is there
unique(GSL.data$year)

  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

# distance histograms
hist(GSL.data$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(GSL.data$distance, main=NULL, xlab="Distance (m)", breaks=c(40))

# clump at 20m - evasive movement. not a nice tail. Could truncate at 40m, 50m, or 60m. I will strat at 60


### update - need to bin data because of lumping at 0

# find appropriate bins
hist(GSL.data$distance[GSL.data$distance<60], breaks=c(0,5,13,18,22,28,33,43,50,60))
hist(GSL.data$distance[GSL.data$distance<60], breaks=c(0,10,19,30,50,60))
hist(GSL.data$distance[GSL.data$distance<60], breaks=c(0,10,30,50,60))


# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.GSL <- 60

# Count the number of observations discarded
nrow(GSL.data[GSL.data$distance>trunc.GSL,]) 

length(GSL.data$distance)

13/193*100 # 6.7%

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(1,2))
plot(GSL.data$size, GSL.data$distance, main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.GSL <- lm(distance~size, data=GSL.data)
lines(GSL.data$size, as.vector(predict(lm.GSL, GSL.data)))
summary(lm.GSL)
# significant relationship between group size and distance  

# Plot of Observer factor against distance
plot(GSL.data$obs.observer,GSL.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# some variation between observers

# I can't include habitat or AMPM in the DF because those variables weren't recorded in earlier years


  ## Fit a detection function ####
    # Uniform ####

## unbinned

# Uni cosine
GSL.df.uni.cos <- ds(data=GSL.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GSL, key="uni", adjustment= "cos")
summary(GSL.df.uni.cos)
# cosine(1)

# Uni poly
GSL.df.uni.poly <- ds(data=GSL.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.GSL, key="uni", adjustment= "poly")
summary(GSL.df.uni.poly)
# poly(2,4)

## compare uni models
gsl.uni.comp <- summarize_ds_models(GSL.df.uni.cos,GSL.df.uni.poly,
                                    output="plain")
gsl.uni.comp[ ,1:4]
gsl.uni.comp[ ,5:7]
# uni cos is the best model. poly has dAIC = 2.7

# plot the fit
par(mfrow=c(1,2))

# cos
plot(GSL.df.uni.cos, main = "GSL.df.uni.cos")

covar.fit <- ddf.gof(GSL.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM
# not a nice QQ plot. 

# GSL.df.uni.cos selected


### binned


# Uni cosine
GSL.df.uni.cos.bin <- ds(data=GSL.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                     truncation=trunc.GSL, key="uni", adjustment= "cos",
                     cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.uni.cos.bin)
# cosine(1)

# Uni poly
GSL.df.uni.poly.bin <- ds(data=GSL.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table, 
                         truncation=trunc.GSL, key="uni", adjustment= "poly",
                         cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.uni.poly.bin)
# poly(2,4)

## compare uni models
gsl.uni.comp.bin <- summarize_ds_models(GSL.df.uni.cos.bin,GSL.df.uni.poly.bin,
                                    output="plain")
gsl.uni.comp.bin[ ,1:4]
gsl.uni.comp.bin[ ,5:7]
# cos has most support

# plot
par(mfrow=c(1,2))
plot(GSL.df.uni.cos.bin, main="cos")
plot(GSL.df.uni.poly.bin, main="poly")

# GSL.df.uni.cos.bin selected


    # Half normal ####


## unbinned 


# HN cosine
GSL.df.hn.cos <- ds(data=GSL.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GSL, key="hn", adjustment= "cos")
summary(GSL.df.hn.cos)
# no adjustment selected 

# HN hermite
GSL.df.hn.herm <- ds(data=GSL.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GSL, key="hn", adjustment= "herm")
summary(GSL.df.hn.herm)
# no adjustment selected 


## half normal with no adjustment is selected

# plot the fit
par(mfrow=c(1,2))
plot(GSL.df.hn.cos, main = "GSL.df.hn.cos")

covar.fit <- ddf.gof(GSL.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# This looks sensible. The spike around 20m is not being modelled, and p at the various distances seem realistic to me.  p=0.5 around 30m, and close to 0 at 60m. Looks very similar to the UNi model.


### binned


# HN cosine
GSL.df.hn.cos.bin <- ds(data=GSL.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                    truncation=trunc.GSL, key="hn", adjustment= "cos",
                    cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.hn.cos.bin)
# key only

# HN herm
GSL.df.hn.herm.bin <- ds(data=GSL.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.GSL, key="hn", adjustment= "herm",
                        cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.hn.herm.bin)
# key only

plot(GSL.df.hn.cos.bin)



    # Hazard rate ####

### unbinned


# HR cosine
GSL.df.hr.cos <- ds(data=GSL.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GSL, key="hr", adjustment= "cos")
summary(GSL.df.hr.cos)
# no adjustment selected 


# HR poly
GSL.df.hr.poly <- ds(data=GSL.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GSL, key="hr", adjustment= "poly")
summary(GSL.df.hr.poly)
# no adjustment selected 


# plot the fit
par(mfrow=c(1,2))

# hr cos
plot(GSL.df.hr.cos, main = "GSL.df.hr.cos")

covar.fit <- ddf.gof(GSL.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# Hmmm, I suspected this is what the HR would do. It is taking the spike into account and giving it a shoulder of p=1 out to just below 20m.  This is not beyond the realms of possibility - GSL live in large groups and can be noisy.  But they can also freeze and hide in response to observers. 

# GSL.df.hr.cos selected


### binned


# HR cosine
GSL.df.hr.cos.bin <- ds(data=GSL.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                    truncation=trunc.GSL, key="hr", adjustment= "cos",
                    cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.hr.cos.bin)
# cosine(2)


# HR poly
GSL.df.hr.poly.bin <- ds(data=GSL.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.GSL, key="hr", adjustment= "poly",
                        cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.hr.poly.bin)
# poly(2)

## compare HR models
gsl.hr.comp.bin <- summarize_ds_models(GSL.df.hr.cos.bin,GSL.df.hr.poly.bin,
                                        output="plain")
gsl.hr.comp.bin[ ,1:4]
gsl.hr.comp.bin[ ,5:7]
# both models have support

par(mfrow=c(1,2))
plot(GSL.df.hr.cos.bin, main="cos")
plot(GSL.df.hr.poly.bin, main="poly")
# cos appears to overfit and wobbles around without really achieving much improvement I don't think

# GSL.df.hr.poly.bin selected


    # Compare primary models ####

### unbinned

gsl.df.prim.comp <- summarize_ds_models(GSL.df.uni.cos, GSL.df.hn.cos, GSL.df.hr.cos, 
                                           output = "plain")
gsl.df.prim.comp[ ,1:5]
gsl.df.prim.comp[ ,6:7]

# half normal has the most support, although uniform does too.  Some debate about which fit is the most ecologically accurate. The species lives in large groups, and are quite flighty. If I remember correctly, the group will panic when the observers are detected, and will flee quite noisily. Therefore the HR model is not necessarily wrong. Olly has checked with the monitoring team and they support what I remeber about the behaviour of the species. The HR model is potentially the most ecologically accurate

summary(GSL.df.hn.cos)
summary(GSL.df.hr.cos)

# extract and plot to see what the different models say
hn <- data.frame(label = "hn",
                   estimate = GSL.df.hn.cos$dht$clusters$N$Estimate[c(1:7)],
                   lcl = GSL.df.hn.cos$dht$clusters$N$lcl[c(1:7)],
                   ucl = GSL.df.hn.cos$dht$clusters$N$ucl[c(1:7)])

hr <- data.frame(label = "hr",
                   estimate = GSL.df.hr.cos$dht$clusters$N$Estimate[c(1:7)],
                   lcl = GSL.df.hr.cos$dht$clusters$N$lcl[c(1:7)],
                   ucl = GSL.df.hr.cos$dht$clusters$N$ucl[c(1:7)])


mod.check <- rbind(hn,hr)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=2)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))


# the models do not differ in the picture the paint, and our inference regarding trends will not change. HR model produces slightly smaller estimates than the HN model. CVs and SEs are smaller for the HR model, and I think that it is a better reflection of the observation process. 


### binned

gsl.df.prim.comp.bin <- summarize_ds_models(GSL.df.uni.cos.bin, GSL.df.hn.cos.bin, GSL.df.hr.poly.bin, 
                                        output = "plain")
gsl.df.prim.comp.bin[ ,1:5]
gsl.df.prim.comp.bin[ ,6:7]
# apparently HR has no support (which I'm surprised by)

# plot together
par(mfrow=c(2,2))
plot(GSL.df.uni.cos.bin, main="uni")
plot(GSL.df.hn.cos.bin, main="hn")
plot(GSL.df.hr.poly.bin, main="hr")
# I think HR looks the best to me

# check estimates
gsl.prim.bin.comp <- data.frame(model = rep(c("uni","hn", "hr"), each=7),
                                year = c(GSL.df.uni.cos.bin$dht$individuals$N$Label[1:7],
                                         GSL.df.hn.cos.bin$dht$individuals$N$Label[1:7],
                                         GSL.df.hr.poly.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(GSL.df.uni.cos.bin$dht$individuals$N$Estimate[1:7],
                                             GSL.df.hn.cos.bin$dht$individuals$N$Estimate[1:7],
                                             GSL.df.hr.poly.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(GSL.df.uni.cos.bin$dht$individuals$N$cv[1:7],
                                       GSL.df.hn.cos.bin$dht$individuals$N$cv[1:7],
                                       GSL.df.hr.poly.bin$dht$individuals$N$cv[1:7]),
                                se = c(GSL.df.uni.cos.bin$dht$individuals$N$se[1:7],
                                       GSL.df.hn.cos.bin$dht$individuals$N$se[1:7],
                                       GSL.df.hr.poly.bin$dht$individuals$N$se[1:7]),
                                lcl = c(GSL.df.uni.cos.bin$dht$individuals$N$lcl[1:7],
                                        GSL.df.hn.cos.bin$dht$individuals$N$lcl[1:7],
                                        GSL.df.hr.poly.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(GSL.df.uni.cos.bin$dht$individuals$N$ucl[1:7],
                                        GSL.df.hn.cos.bin$dht$individuals$N$ucl[1:7],
                                        GSL.df.hr.poly.bin$dht$individuals$N$ucl[1:7]))

gsl.prim.bin.comp

# plot cv
ggplot(gsl.prim.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# HR highest every year. Uni and Hn similar, but uni always lowest

# plot estimates
ggplot(gsl.prim.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# UNi the most precise by a fair way. HR not great precision. Overall though, trend is the same - only small differences in estimates. 

# GSL.df.uni.cos.bin selected


    # Models with harsh truncation ####

# set trunc distance 
gsl.trunc.harsh <- 45
gsl.trunc.harsh.bin <- 50


### unbinned

par(mfrow=c(1,2))

# HN cos
GSL.df.hr.cos.harsh <- ds(data=GSL.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=gsl.trunc.harsh, key="hr", adjustment="cos")
summary(GSL.df.hr.cos.harsh)
# cosine(2,3,4,5)

# compare original with harsh truncation
ddf.gof(GSL.df.hr.cos.harsh$ddf, main = "GSL.df.hr.cos.harsh", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(GSL.df.hr.cos$ddf, main = "GSL.df.hr.cos",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# harsh model much worse


### binned

# Uni cosine
GSL.df.uni.cos.bin.harsh <- ds(data=GSL.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table, 
                         truncation=gsl.trunc.harsh.bin, key="uni", adjustment= "cos",
                         cutpoints = c(0,10,19,30,50))

par(mfrow=c(1,2))
plot(GSL.df.uni.cos.bin.harsh, main="harsh")
plot(GSL.df.uni.cos.bin, main="orig")

# check estimates
gsl.harsh.bin.comp <- data.frame(model = rep(c("harsh","orig"), each=7),
                                year = c(GSL.df.uni.cos.bin.harsh$dht$individuals$N$Label[1:7],
                                         GSL.df.uni.cos.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(GSL.df.uni.cos.bin.harsh$dht$individuals$N$Estimate[1:7],
                                             GSL.df.uni.cos.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(GSL.df.uni.cos.bin.harsh$dht$individuals$N$cv[1:7],
                                       GSL.df.uni.cos.bin$dht$individuals$N$cv[1:7]),
                                se = c(GSL.df.uni.cos.bin.harsh$dht$individuals$N$se[1:7],
                                       GSL.df.uni.cos.bin$dht$individuals$N$se[1:7]),
                                lcl = c(GSL.df.uni.cos.bin.harsh$dht$individuals$N$lcl[1:7],
                                        GSL.df.uni.cos.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(GSL.df.uni.cos.bin.harsh$dht$individuals$N$ucl[1:7],
                                        GSL.df.uni.cos.bin$dht$individuals$N$ucl[1:7]))

gsl.harsh.bin.comp

# plot cv
ggplot(gsl.harsh.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# HR highest every year. Uni and Hn similar, but uni always lowest

# plot estimates
ggplot(gsl.harsh.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# estimates change very little, and the original is the more precise

# original selected


    # Models with covariates ####


### unbinned

## cluster size
GSL.df.hr.size <- ds(data=GSL.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.GSL, key="hr", formula = ~size)
summary(GSL.df.hr.size)

# plot
par(mfrow=c(1,2))
plot(GSL.df.hr.size, main = "GSL.df.hr.size")

covar.fit <- ddf.gof(GSL.df.hr.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# lot of variation



## Observer
GSL.df.hr.obs <- ds(data=GSL.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.GSL, key="hr", formula = ~obs.observer)
summary(GSL.df.hr.obs)

# plot
plot(GSL.df.hr.obs, main = "GSL.df.hr.obs")

covar.fit <- ddf.gof(GSL.df.hr.obs$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


## stratum (year)
GSL.df.hr.strat <- ds(data=GSL.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.GSL, key="hr", formula = ~stratum)
summary(GSL.df.hr.strat)

# plot
plot(GSL.df.hr.strat, main = "GSL.df.hr.strat")

covar.fit <- ddf.gof(GSL.df.hr.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



### binned


# HN size
GSL.df.hn.size.bin <- ds(data=GSL.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.GSL, key="hn", formula = ~size,
                        cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.hn.size.bin)
plot(GSL.df.hn.size.bin)

# HN observer
GSL.df.hn.obs.bin <- ds(data=GSL.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table, 
                         truncation=trunc.GSL, key="hn", formula = ~obs.observer,
                         cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.hn.obs.bin)
plot(GSL.df.hn.obs.bin)

# HN stratum
GSL.df.hn.strat.bin <- ds(data=GSL.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.GSL, key="hn", formula = ~stratum,
                        cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.hn.strat.bin)
plot(GSL.df.hn.strat.bin)


# HN stratum + size
GSL.df.hn.strat.size.bin <- ds(data=GSL.data, region.table=full.region.table, 
                          sample.table=sample.table, obs.table=obs.table, 
                          truncation=trunc.GSL, key="hn", formula = ~stratum+size,
                          cutpoints = c(0,10,19,30,50,60))
summary(GSL.df.hn.strat.bin)
plot(GSL.df.hn.strat.bin)


    # Compare covariate models ####

### unbinned

gsl.df.cov.comp <- summarize_ds_models(GSL.df.hr.size, GSL.df.hr.obs,
                                          GSL.df.hr.strat,output = "plain")

gsl.df.cov.comp[ ,1:5]
gsl.df.cov.comp[ ,6:7]

# models with size and stratum both have support.

# add model with both size and stratum
GSL.df.hr.strat.size <- ds(data=GSL.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.GSL, key="hr", formula = ~stratum+size)

# compare
gsl.df.cov.comp <- summarize_ds_models(GSL.df.hr.size, GSL.df.hr.obs,GSL.df.hr.strat.size,
                                       GSL.df.hr.strat,output = "plain")

gsl.df.cov.comp[ ,1:5]
gsl.df.cov.comp[ ,6:7]
# model with both is poor. So still decision between stratum and size. 

# extract and plot to see what the different models say
size <- data.frame(label = "size",
                   estimate = GSL.df.hr.size$dht$clusters$N$Estimate[c(1:7)],
                   lcl = GSL.df.hr.size$dht$clusters$N$lcl[c(1:7)],
                   ucl = GSL.df.hr.size$dht$clusters$N$ucl[c(1:7)])

strat <- data.frame(label = "strat",
                    estimate = GSL.df.hr.strat$dht$clusters$N$Estimate[c(1:7)],
                    lcl = GSL.df.hr.strat$dht$clusters$N$lcl[c(1:7)],
                    ucl = GSL.df.hr.strat$dht$clusters$N$ucl[c(1:7)])


mod.check <- rbind(size,strat)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=2)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# The overall inference regarding trend will not change based on the model chosen, which is good. The models agree relatively well.  The stratum model does suggest more of a downward trend, but the estimates from strat model also have quite high CV and confidence intervals, and thus precision is worse. Therefore the model with size is selected



### binned

gsl.df.cov.comp.bin <- summarize_ds_models(GSL.df.hn.size.bin, GSL.df.hn.strat.bin, 
                                           GSL.df.hn.strat.size.bin,output = "plain")
gsl.df.cov.comp.bin[ ,1:5]
gsl.df.cov.comp.bin[ ,6:7]
# strat+ size has most support

# check estimates
gsl.cov.bin.comp <- data.frame(model = rep(c("size", "strat", "strat+size"), each=7),
                                year = c(GSL.df.hn.size.bin$dht$individuals$N$Label[1:7],
                                         GSL.df.hn.strat.bin$dht$individuals$N$Label[1:7],
                                         GSL.df.hn.strat.size.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(GSL.df.hn.size.bin$dht$individuals$N$Estimate[1:7],
                                             GSL.df.hn.strat.bin$dht$individuals$N$Estimate[1:7],
                                             GSL.df.hn.strat.size.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(GSL.df.hn.size.bin$dht$individuals$N$cv[1:7],
                                       GSL.df.hn.strat.bin$dht$individuals$N$cv[1:7],
                                       GSL.df.hn.strat.size.bin$dht$individuals$N$cv[1:7]),
                                se = c(GSL.df.hn.size.bin$dht$individuals$N$se[1:7],
                                       GSL.df.hn.strat.bin$dht$individuals$N$se[1:7],
                                       GSL.df.hn.strat.size.bin$dht$individuals$N$se[1:7]),
                                lcl = c(GSL.df.hn.size.bin$dht$individuals$N$lcl[1:7],
                                        GSL.df.hn.strat.bin$dht$individuals$N$lcl[1:7],
                                        GSL.df.hn.strat.size.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(GSL.df.hn.size.bin$dht$individuals$N$ucl[1:7],
                                        GSL.df.hn.strat.bin$dht$individuals$N$ucl[1:7],
                                        GSL.df.hn.strat.size.bin$dht$individuals$N$ucl[1:7]))

gsl.cov.bin.comp

# plot cv
ggplot(gsl.cov.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# Hsome variation in CVs

# plot estimates
ggplot(gsl.cov.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))

# GSL.df.hn.strat.size.bin selected


    # Compare primary and covariate models ####


### unbinned

gsl.df.final.compare <- summarize_ds_models(GSL.df.hr.cos,GSL.df.hr.size, output = "plain")
gsl.df.final.compare[ ,1:3]
gsl.df.final.compare[ ,4:7]

# covariate model is the better model, and is selected


### binned


gsl.df.final.compare.bin <- summarize_ds_models(GSL.df.hn.strat.size.bin,GSL.df.uni.cos.bin, 
                                                output = "plain")
gsl.df.final.compare.bin[ ,1:3]
gsl.df.final.compare.bin[ ,4:7]
# original model has no support


# check estimates
gsl.final.bin.comp <- data.frame(model = rep(c("orig", "covar"), each=7),
                               year = c(GSL.df.uni.cos.bin$dht$individuals$N$Label[1:7],
                                        GSL.df.hn.strat.size.bin$dht$individuals$N$Label[1:7]),
                               estimate = c(GSL.df.uni.cos.bin$dht$individuals$N$Estimate[1:7],
                                            GSL.df.hn.strat.size.bin$dht$individuals$N$Estimate[1:7]),
                               cv = c(GSL.df.uni.cos.bin$dht$individuals$N$cv[1:7],
                                      GSL.df.hn.strat.size.bin$dht$individuals$N$cv[1:7]),
                               se = c(GSL.df.uni.cos.bin$dht$individuals$N$se[1:7],
                                      GSL.df.hn.strat.size.bin$dht$individuals$N$se[1:7]),
                               lcl = c(GSL.df.uni.cos.bin$dht$individuals$N$lcl[1:7],
                                       GSL.df.hn.strat.size.bin$dht$individuals$N$lcl[1:7]),
                               ucl = c(GSL.df.uni.cos.bin$dht$individuals$N$ucl[1:7],
                                       GSL.df.hn.strat.size.bin$dht$individuals$N$ucl[1:7]))

gsl.final.bin.comp

# plot cv
ggplot(gsl.final.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# some variation in CVs

# plot estimates
ggplot(gsl.final.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# covariate model is less precise in earlier years but more precise in later years. AIC suggests it is by far the better model

# covar model selected

  ## GSL Final Results ####

### unbinned

summary(GSL.df.hr.size)

# extract estimates
GSL.grp.results <- GSL.df.hr.size$dht$clusters$N[1:7, ]
GSL.ind.results <- GSL.df.hr.size$dht$individuals$N[1:7, ]

# bind results
GSL.results <- rbind(GSL.grp.results, GSL.ind.results)
GSL.results <- GSL.results %>% rename(Year = Label) %>% 
                mutate(Label = rep(c("Grp", "Ind"), each=7)) %>% 
                mutate(Species = rep("GSL", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hr", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("size", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)



GSL_final_plot <- ggplot(GSL.results[GSL.results$Label=="Grp",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  ylim(0,2100)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(GSL.results, "Output/Results/GSL_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/GSL_final_plot.png", GSL_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")



### binned

summary(GSL.df.hn.strat.size.bin)

# extract estimates
GSL.grp.abund.bin   <- GSL.df.hn.strat.size.bin$dht$clusters$N[1:7,1:6]
GSL.ind.abund.bin   <- GSL.df.hn.strat.size.bin$dht$individuals$N[1:7,1:6]
GSL.grp.density.bin <- GSL.df.hn.strat.size.bin$dht$clusters$D[1:7,1:6]
GSL.ind.density.bin <- GSL.df.hn.strat.size.bin$dht$individuals$D[1:7,1:6]

GSL.grp.abund.bin <- GSL.grp.abund.bin %>% 
                      rename(Year=Label,N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl) %>% 
                      mutate(Label="Grp")
GSL.ind.abund.bin <- GSL.ind.abund.bin %>% 
                      rename(Year=Label,N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl) %>% 
                      mutate(Label="Ind")
GSL.grp.density.bin <- GSL.grp.density.bin %>% 
                      rename(Year=Label,D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl) %>% 
                      mutate(Label="Grp")
GSL.ind.density.bin <- GSL.ind.density.bin %>% 
                      rename(Year=Label,D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl) %>% 
                      mutate(Label="Ind")


# merge 
GSL.N.bin <- rbind(GSL.grp.abund.bin,GSL.ind.abund.bin)
GSL.N.bin <- GSL.N.bin %>% arrange(Label,Year)
GSL.D.bin <- rbind(GSL.grp.density.bin,GSL.ind.density.bin)
GSL.D.bin <- GSL.D.bin %>% arrange(Label, Year)
GSL.D.bin <- GSL.D.bin[, -c(1,7)]
GSL.results.final <- cbind(GSL.N.bin,GSL.D.bin)


# bind results

GSL.results.final <- GSL.results.final %>% 
                      mutate(Species = rep("GSL", times=14)) %>% 
                      mutate(DetFun = rep("pooled", times=14)) %>% 
                      mutate(Key = rep("Hn", times=14)) %>% 
                      mutate(Adjust = rep("NA", times=14)) %>% 
                      mutate(Covar = rep("size + year", times=14)) %>% 
                      select(Year,Species,DetFun,Key,Adjust,Covar,Label,N,n_se,n_cv,n_lcl,n_ucl,
                             D,d_se,d_cv,d_lcl,d_ucl)

# convert density estimates from m2 to km2
GSL.results.final <- GSL.results.final %>% 
                      mutate(D = D*1000000) %>% 
                      mutate(d_se = d_se*1000000) %>% 
                      mutate(d_lcl = d_lcl*1000000) %>% 
                      mutate(d_ucl = d_ucl*1000000)


GSL_final_plot_binned <- ggplot(GSL.results.bin[GSL.results.bin$Label=="Grp",], aes(x=Year, y=Estimate))+
  geom_point(shape=16, size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
  ylim(0,2500)+
  theme_bw()+
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank())+
  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(GSL.results.final, "Output/Results/GSL_results_final_binned.csv")
ggsave("Output/Results/Plots/Point_estimates/GSL_final_plot_binned.png", GSL_final_plot_binned, 
       dpi=300, width = 20, height = 20, units = "cm")



### extract density
GSL.grp.density.bin <- GSL.df.hn.strat.size.bin$dht$clusters$D[1:7, ]
GSL.grp.density.bin <- GSL.grp.density.bin %>% rename(Year=Label) %>% 
                        mutate(Species="GSL") %>% 
                        select(Species,Year,Estimate,se,cv,lcl,ucl,df)

# compare binned and unbinned

# load unbinned
GSL.results <- read.csv("Output/Results/GSL_results_final.csv")
GSL.results <- GSL.results[,-1]
GSL.results$analysis <- "original"

GSL.results.bin$analysis <- "binned"

GSL.results.all <- rbind(GSL.results,GSL.results.bin)

GSL_analysis_comparison_plot <- ggplot(GSL.results.all[GSL.results.all$Label=="Grp",], 
                                       aes(x=Year, y=Estimate, colour=analysis))+
                    geom_point(shape=16, size=2, position = position_dodge(width=0.3))+
                    geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width=0.3))+
                    ylim(0,2400)+
                    theme_bw()+
                    theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),panel.border = element_blank())+
                    theme(axis.line = element_line(color = 'black'))

ggsave("Output/Results/Plots/Point_estimates/GSL_analysis_comparison_plot.png", 
       GSL_analysis_comparison_plot, dpi=300, width = 20, height = 20, units = "cm")


#### Long-tailed macaque #####################################################
  ## Subset data ####

# subset LTM data
LTM.data <- allData[allData$species=="LTM",] 
head(LTM.data)
LTM.data <- as.data.frame(LTM.data)

# Total number of groups from all years
length(LTM.data$distance) 

# for LTM, I will pool all years for DF, but test year as a continous covariate in the DF model

# check 2020 is there
unique(LTM.data$year)
str(LTM.data)

  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

# distance histograms
hist(LTM.data$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(LTM.data$distance, main=NULL, xlab="Distance (m)", breaks=c(40))

# not pretty data when more bins are used. 

### update - I will be binning data due to clumping at 0

# identify best bins
hist(LTM.data$distance[LTM.data$distance<50], breaks=c(0,10,15,20,25,35,50))
hist(LTM.data$distance[LTM.data$distance<50], breaks=c(0,9,13,20,25,35,50))
hist(LTM.data$distance[LTM.data$distance<50], breaks=c(0,10,20,27,35,50))


# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.LTM <- 50

# Count the number of observations discarded
nrow(LTM.data[LTM.data$distance>trunc.LTM,]) 

length(LTM.data$distance)

14/170*100 # 8.2%

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(1,2))
plot(LTM.data$size, LTM.data$distance, main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.LTM <- lm(distance~size, data=LTM.data)
lines(LTM.data$size, as.vector(predict(lm.LTM, LTM.data)))
summary(lm.LTM)
# no evidence of size bias 

# Plot of Observer factor against distance
plot(LTM.data$obs.observer,LTM.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# quite a lot of variation between observers

# I can't include habitat or AMPM in the DF because those variables weren't recorded in earlier years


  ## Fit a detection function ####
    # Uniform ####

### unbinned

# Uni cosine

LTM.df.uni.cos <- ds(data=LTM.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.LTM, key="uni", adjustment= "cos")
summary(LTM.df.uni.cos)
# cosine(1)

# Uni poly
LTM.df.uni.poly <- ds(data=LTM.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.LTM, key="uni", adjustment= "poly")
summary(LTM.df.uni.poly)
# poly(2)

## compare uni models
ltm.uni.comp <- summarize_ds_models(LTM.df.uni.cos, LTM.df.uni.poly,
                                    output="plain")
ltm.uni.comp[ ,1:4]
ltm.uni.comp[ ,5:7]
# poly is better but cos has dAIC < 2

# plot the fits
par(mfrow=c(2,2))

# cos
plot(LTM.df.uni.cos, main = "LTM.df.uni.cos")

covar.fit <- ddf.gof(LTM.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(LTM.df.uni.poly, main = "LTM.df.uni.poly")

covar.fit <- ddf.gof(LTM.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# The cos model QQ plot is nicer, but I prefer the fit of the poly model, as it suggest p decreases slower than the cos model. I think this is true ecolgically for similar reasons to the GSL. This species lives in large groups, and are noisy when moving / fleeing. Therefore I don't think p dops to 0.5 by 30m, as the cos model suggests.

# poly selected


### binned


# Uni cosine
LTM.df.uni.cos.bin <- ds(data=LTM.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                     truncation=trunc.LTM, key="uni", adjustment= "cos",
                     cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.uni.cos.bin)
# cosine(1)

# Uni poly
LTM.df.uni.poly.bin <- ds(data=LTM.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table, 
                         truncation=trunc.LTM, key="uni", adjustment= "poly",
                         cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.uni.poly.bin)
# poly(2)

## compare uni models
ltm.uni.comp.bin <- summarize_ds_models(LTM.df.uni.cos.bin, LTM.df.uni.poly.bin,
                                    output="plain")
ltm.uni.comp.bin[ ,1:4]
ltm.uni.comp.bin[ ,5:7]
# cos model has most support

# plot
par(mfrow=c(1,2))
plot(LTM.df.uni.cos.bin, main="cos")
plot(LTM.df.uni.poly.bin, main="poly")
# look virtually identical. cos selected


    # Half normal ####

### unbinned

# HN cosine
LTM.df.hn.cos <- ds(data=LTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,  
                      truncation=trunc.LTM, key="hn", adjustment= "cos")
summary(LTM.df.hn.cos)
# no adjustment selected 

# HN hermite
LTM.df.hn.herm <- ds(data=LTM.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.LTM, key="hn", adjustment= "herm")
summary(LTM.df.hn.herm)
# no adjustment selected 


## half normal with no adjustment is selected

# plot the fit
par(mfrow=c(1,2))
plot(LTM.df.hn.cos, main = "LTM.df.hn.cos")

covar.fit <- ddf.gof(LTM.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# This suggest p drops off quite quickly as well


### binned

# HN cosine
LTM.df.hn.cos.bin <- ds(data=LTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,  
                    truncation=trunc.LTM, key="hn", adjustment= "cos",
                    cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.hn.cos.bin)
# key only

# HN hermite
LTM.df.hn.herm.bin <- ds(data=LTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table,  
                        truncation=trunc.LTM, key="hn", adjustment= "herm",
                        cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.hn.herm.bin)
# key only

# plot
plot(LTM.df.hn.cos.bin, main="bin")



    # Hazard rate ####


### unbinned

# HR cosine
LTM.df.hr.cos <- ds(data=LTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.LTM, key="hr", adjustment= "cos")
summary(LTM.df.hr.cos)
# cosine(2)


# HR poly
LTM.df.hr.poly <- ds(data=LTM.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.LTM, key="hr", adjustment= "poly")
summary(LTM.df.hr.poly)
# no adjustment selected 


## compare HR models
ltm.hr.comp <- summarize_ds_models(LTM.df.hr.cos, LTM.df.hr.poly,
                                    output="plain")
ltm.hr.comp[ ,1:4]
ltm.hr.comp[ ,5:7]
# Cos is the preferred model, but not much in it - poly has dAIC=1.1

# plot the fits
par(mfrow=c(2,2))

# hr cos
plot(LTM.df.hr.cos, main = "LTM.df.hr.cos")

covar.fit <- ddf.gof(LTM.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# hr poly (no adjustment)
plot(LTM.df.hr.poly, main = "LTM.df.hr.poly")

covar.fit <- ddf.gof(LTM.df.hr.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# poly model (which is actually hr with no adjustment) looks much better to me. The cos model has p decreasing too fast, and is overfitting between 20-30m. The poly model has a small shoulder which I think is realistic. The poly model also has a nicer QQ plot

# poly selected


### binned

# HR cosine
LTM.df.hr.cos.bin <- ds(data=LTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                    truncation=trunc.LTM, key="hr", adjustment= "cos",
                    cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.hr.cos.bin)
# key only

# HR poly
LTM.df.hr.poly.bin <- ds(data=LTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.LTM, key="hr", adjustment= "poly",
                        cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.hr.poly.bin)
# key only

# plot
par(mfrow=c(1,1))
plot(LTM.df.hr.cos.bin)


    # Compare primary models ####

### unbinned

ltm.df.prim.comp <- summarize_ds_models(LTM.df.uni.poly, LTM.df.hn.cos,LTM.df.hr.poly, 
                                           output = "plain")
ltm.df.prim.comp[ ,1:5]
ltm.df.prim.comp[ ,6:7]

# HN and Uni have the most support.  HR dAIC<3. 

# plot all fits together
par(mfrow=c(3,2))

# uni poly
plot(LTM.df.uni.poly, main = "LTM.df.uni.poly")

covar.fit <- ddf.gof(LTM.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HN cos
plot(LTM.df.hn.cos, main = "LTM.df.hn.cos")

covar.fit <- ddf.gof(LTM.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HR poly
plot(LTM.df.hr.poly, main = "LTM.df.hr.poly")

covar.fit <- ddf.gof(LTM.df.hr.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# I prefer the UNi model as the other two suggest p reaches 0.5 by 30m, which I think is too soon for this species. Uni is more like 40m. The qq plot for Uni isn't the best, but it has support according to AIC. I will take the Uni model forward for teting harsh trunctation, but the HN model forward for testing covariates. I will compare the best covar model to the Uni model. 


### binned

ltm.df.prim.comp.bin <- summarize_ds_models(LTM.df.uni.cos.bin, LTM.df.hn.cos.bin,LTM.df.hr.cos.bin, 
                                        output = "plain")
ltm.df.prim.comp.bin[ ,1:5]
ltm.df.prim.comp.bin[ ,6:7]
# all models have some support

# plot
par(mfrow=c(2,2))
plot(LTM.df.uni.cos.bin, main="uni")
plot(LTM.df.hn.cos.bin, main="hn")
plot(LTM.df.hr.cos.bin, main="hr")

# check estimates
ltm.prim.bin.comp <- data.frame(model = rep(c("uni", "hn", "hr"), each=7),
                               year = c(LTM.df.uni.cos.bin$dht$individuals$N$Label[1:7],
                                        LTM.df.hn.cos.bin$dht$individuals$N$Label[1:7],
                                        LTM.df.hr.cos.bin$dht$individuals$N$Label[1:7]),
                               estimate = c(LTM.df.uni.cos.bin$dht$individuals$N$Estimate[1:7],
                                            LTM.df.hn.cos.bin$dht$individuals$N$Estimate[1:7],
                                            LTM.df.hr.cos.bin$dht$individuals$N$Estimate[1:7]),
                               cv = c(LTM.df.uni.cos.bin$dht$individuals$N$cv[1:7],
                                      LTM.df.hn.cos.bin$dht$individuals$N$cv[1:7],
                                      LTM.df.hr.cos.bin$dht$individuals$N$cv[1:7]),
                               se = c(LTM.df.uni.cos.bin$dht$individuals$N$se[1:7],
                                      LTM.df.hn.cos.bin$dht$individuals$N$se[1:7],
                                      LTM.df.hr.cos.bin$dht$individuals$N$se[1:7]),
                               lcl = c(LTM.df.uni.cos.bin$dht$individuals$N$lcl[1:7],
                                       LTM.df.hn.cos.bin$dht$individuals$N$lcl[1:7],
                                       LTM.df.hr.cos.bin$dht$individuals$N$lcl[1:7]),
                               ucl = c(LTM.df.uni.cos.bin$dht$individuals$N$ucl[1:7],
                                       LTM.df.hn.cos.bin$dht$individuals$N$ucl[1:7],
                                       LTM.df.hr.cos.bin$dht$individuals$N$ucl[1:7]))

ltm.prim.bin.comp

# plot cv
ggplot(ltm.prim.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# HR always the highest CV

# plot estimates
ggplot(ltm.prim.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# All very similar. Hn has consitently more precision

# LTM.df.hn.cos.bin selected


    # Models with harsh truncation ####


# set trunc distance 
ltm.trunc.harsh <- 40
ltm.trunc.harsh.bin <- 35


### unbinned

par(mfrow=c(1,2))

# HN cos
LTM.df.uni.poly.harsh <- ds(data=LTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=ltm.trunc.harsh, key="uni", adjustment="poly")
summary(LTM.df.uni.poly.harsh)
# poly(2) selected

# compare original with harsh truncation
ddf.gof(LTM.df.uni.poly.harsh$ddf, main = "LTM.df.uni.poly.harsh", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(LTM.df.uni.poly$ddf, main = "LTM.df.uni.poly",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# pretty much identical - harsh truncation doesn't appear to have improved the fit. 


# extract and plot to see what the different models say
harsh <- data.frame(label = "harsh",
                   estimate = LTM.df.uni.poly.harsh$dht$clusters$N$Estimate[c(1:7)],
                   lcl = LTM.df.uni.poly.harsh$dht$clusters$N$lcl[c(1:7)],
                   ucl = LTM.df.uni.poly.harsh$dht$clusters$N$ucl[c(1:7)])

orig <- data.frame(label = "orig",
                    estimate = LTM.df.uni.poly$dht$clusters$N$Estimate[c(1:7)],
                    lcl = LTM.df.uni.poly$dht$clusters$N$lcl[c(1:7)],
                    ucl = LTM.df.uni.poly$dht$clusters$N$ucl[c(1:7)])


mod.check <- rbind(harsh,orig)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=2)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# very similar estimates, with the same trend. Original model has slightly better precision and so is selected.


### binned

# HN cosine
LTM.df.hn.cos.harsh.bin <- ds(data=LTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table,  
                        truncation=ltm.trunc.harsh.bin, key="hn", adjustment= "cos",
                        cutpoints = c(0,10,20,27,35))
summary(LTM.df.hn.cos.harsh.bin)
# key only

# plot together
par(mfrow=c(1,2))
plot(LTM.df.hn.cos.bin, main="orig")
plot(LTM.df.hn.cos.harsh.bin, main="harsh")


# check estimates
ltm.harsh.bin.comp <- data.frame(model = rep(c("orig", "harsh"), each=7),
                                year = c(LTM.df.hn.cos.bin$dht$individuals$N$Label[1:7],
                                         LTM.df.hn.cos.harsh.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(LTM.df.hn.cos.bin$dht$individuals$N$Estimate[1:7],
                                             LTM.df.hn.cos.harsh.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(LTM.df.hn.cos.bin$dht$individuals$N$cv[1:7],
                                       LTM.df.hn.cos.harsh.bin$dht$individuals$N$cv[1:7]),
                                se = c(LTM.df.hn.cos.bin$dht$individuals$N$se[1:7],
                                       LTM.df.hn.cos.harsh.bin$dht$individuals$N$se[1:7]),
                                lcl = c(LTM.df.hn.cos.bin$dht$individuals$N$lcl[1:7],
                                        LTM.df.hn.cos.harsh.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(LTM.df.hn.cos.bin$dht$individuals$N$ucl[1:7],
                                        LTM.df.hn.cos.harsh.bin$dht$individuals$N$ucl[1:7]))

ltm.harsh.bin.comp

# plot cv
ggplot(ltm.harsh.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# harsh has much higher CV in most years

# plot estimates
ggplot(ltm.harsh.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# original model has greater precision and is selected

    # Models with covariates ####

### unbinned

## cluster size
LTM.df.hn.size <- ds(data=LTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.LTM, key="hn", formula = ~size)
summary(LTM.df.hn.size)

# plot
par(mfrow=c(1,2))
plot(LTM.df.hn.size, main = "LTM.df.hn.size")

covar.fit <- ddf.gof(LTM.df.hn.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# looks to be a very poor fit


## Observer
LTM.df.hn.obs <- ds(data=LTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.LTM, key="hn", formula = ~obs.observer)
summary(LTM.df.hn.obs)

# plot
plot(LTM.df.hn.obs, main = "LTM.df.hn.obs")

covar.fit <- ddf.gof(LTM.df.hn.obs$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# better


## stratum (year)
LTM.df.hn.strat <- ds(data=LTM.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.LTM, key="hn", formula = ~stratum)
summary(LTM.df.hn.strat)

# plot
plot(LTM.df.hn.strat, main = "LTM.df.hn.strat")

covar.fit <- ddf.gof(LTM.df.hn.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


### binned

# HN size
LTM.df.hn.size.bin <- ds(data=LTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table,  
                        truncation=trunc.LTM, key="hn", formula=~size,
                        cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.hn.size.bin)
plot(LTM.df.hn.size.bin)

# HN observer
LTM.df.hn.obs.bin <- ds(data=LTM.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table,  
                         truncation=trunc.LTM, key="hn", formula=~obs.observer,
                         cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.hn.obs.bin)
plot(LTM.df.hn.obs.bin)


# HN stratum
LTM.df.hn.strat.bin <- ds(data=LTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table,  
                        truncation=trunc.LTM, key="hn", formula=~stratum,
                        cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.hn.strat.bin)
plot(LTM.df.hn.strat.bin)

# HN stratum + size
LTM.df.hn.strat.size.bin <- ds(data=LTM.data, region.table=full.region.table, 
                          sample.table=sample.table, obs.table=obs.table,  
                          truncation=trunc.LTM, key="hn", formula=~stratum+size,
                          cutpoints = c(0,10,20,27,35,50))
summary(LTM.df.hn.strat.size.bin)
plot(LTM.df.hn.strat.size.bin)


    # Compare covariate models ####


### unbinned

ltm.df.cov.comp <- summarize_ds_models(LTM.df.hn.size, LTM.df.hn.obs,
                                          LTM.df.hn.strat,output = "plain")

ltm.df.cov.comp[ ,1:5]
ltm.df.cov.comp[ ,6:7]

# Models with stratum and observer have dAIC > 100 (!?), and so the model with size is selected


### binned

ltm.df.cov.comp.bin <- summarize_ds_models(LTM.df.hn.size.bin, LTM.df.hn.obs.bin,
                                           LTM.df.hn.strat.bin,LTM.df.hn.strat.size.bin,
                                           output = "plain")

ltm.df.cov.comp.bin[ ,1:5]
ltm.df.cov.comp.bin[ ,6:7]
# all models except observer have some support

# check estimates
ltm.cov.bin.comp <- data.frame(model = rep(c("size", "strat", "size+strat"), each=7),
                                year = c(LTM.df.hn.size.bin$dht$individuals$N$Label[1:7],
                                         LTM.df.hn.strat.bin$dht$individuals$N$Label[1:7],
                                         LTM.df.hn.strat.size.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(LTM.df.hn.size.bin$dht$individuals$N$Estimate[1:7],
                                             LTM.df.hn.strat.bin$dht$individuals$N$Estimate[1:7],
                                             LTM.df.hn.strat.size.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(LTM.df.hn.size.bin$dht$individuals$N$cv[1:7],
                                       LTM.df.hn.strat.bin$dht$individuals$N$cv[1:7],
                                       LTM.df.hn.strat.size.bin$dht$individuals$N$cv[1:7]),
                                se = c(LTM.df.hn.size.bin$dht$individuals$N$se[1:7],
                                       LTM.df.hn.strat.bin$dht$individuals$N$se[1:7],
                                       LTM.df.hn.strat.size.bin$dht$individuals$N$se[1:7]),
                                lcl = c(LTM.df.hn.size.bin$dht$individuals$N$lcl[1:7],
                                        LTM.df.hn.strat.bin$dht$individuals$N$lcl[1:7],
                                        LTM.df.hn.strat.size.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(LTM.df.hn.size.bin$dht$individuals$N$ucl[1:7],
                                        LTM.df.hn.strat.bin$dht$individuals$N$ucl[1:7],
                                        LTM.df.hn.strat.size.bin$dht$individuals$N$ucl[1:7]))

ltm.cov.bin.comp

# plot cv
ggplot(ltm.cov.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# size mostly has lowest CV, strat or size+strat have highest. 

# plot estimates
ggplot(ltm.cov.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# All pretty similar estimates. size has the best precision in earlier years but strat+size tends to be better in later years. The size+strat model tends to sit in between the other two estimates, and is quite a good compromise in terms of precision. 

# LTM.df.hn.strat.size.bin selected


    # Compare primary and covariate models ####

# unbinned


ltm.df.final.compare <- summarize_ds_models(LTM.df.uni.poly,LTM.df.hn.size, output = "plain")
ltm.df.final.compare[ ,1:3]
ltm.df.final.compare[ ,4:7]

# covariate model is the better model (by a loooong way), and is selected


# Just out of curiosity, extract and plot to see what the different models say
uni.poly <- data.frame(label = "uni.poly",
                    estimate = LTM.df.uni.poly$dht$clusters$N$Estimate[c(1:7)],
                    lcl = LTM.df.uni.poly$dht$clusters$N$lcl[c(1:7)],
                    ucl = LTM.df.uni.poly$dht$clusters$N$ucl[c(1:7)])

size <- data.frame(label = "size",
                   estimate = LTM.df.hn.size$dht$clusters$N$Estimate[c(1:7)],
                   lcl = LTM.df.hn.size$dht$clusters$N$lcl[c(1:7)],
                   ucl = LTM.df.hn.size$dht$clusters$N$ucl[c(1:7)])


mod.check <- rbind(uni.poly,size)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=2)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))


### binned

ltm.df.final.compare.bin <- summarize_ds_models(LTM.df.hn.cos.bin,LTM.df.hn.strat.size.bin, output = "plain")
ltm.df.final.compare.bin[ ,1:3]
ltm.df.final.compare.bin[ ,4:7]
# AIC prefers the original model over the covariate model

# check estimates
ltm.final.bin.comp <- data.frame(model = rep(c("orig", "covar"), each=7),
                                 year = c(LTM.df.hn.cos.bin$dht$individuals$N$Label[1:7],
                                          LTM.df.hn.strat.size.bin$dht$individuals$N$Label[1:7]),
                                 estimate = c(LTM.df.hn.cos.bin$dht$individuals$N$Estimate[1:7],
                                              LTM.df.hn.strat.size.bin$dht$individuals$N$Estimate[1:7]),
                                 cv = c(LTM.df.hn.cos.bin$dht$individuals$N$cv[1:7],
                                        LTM.df.hn.strat.size.bin$dht$individuals$N$cv[1:7]),
                                 se = c(LTM.df.hn.cos.bin$dht$individuals$N$se[1:7],
                                        LTM.df.hn.strat.size.bin$dht$individuals$N$se[1:7]),
                                 lcl = c(LTM.df.hn.cos.bin$dht$individuals$N$lcl[1:7],
                                         LTM.df.hn.strat.size.bin$dht$individuals$N$lcl[1:7]),
                                 ucl = c(LTM.df.hn.cos.bin$dht$individuals$N$ucl[1:7],
                                         LTM.df.hn.strat.size.bin$dht$individuals$N$ucl[1:7]))

ltm.final.bin.comp

# plot cv
ggplot(ltm.final.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# some variation but original has lower CV in more years

# plot estimates
ggplot(ltm.final.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# estimates are similar, excpet in 2018 where the difference is larger. The covar model is consistently estiamting a lower abundance than the original model. Because we know that macaques have been persecuted in certain years for the lab trade, and becasue they are partly terrestrial (and so vulnerable to snares), I am inclined to go with the covar model as it takes into account changes in p over time.  It is also produces more precise estimates. 

  ## LTM Final Results ####

### unbinned

summary(LTM.df.hn.size)

dht2.LTM.output <- dht2(LTM.df.hn.size, observations=obs.table, transects=sample.table, 
                        geo_strat=full.region.table, strat_formula = ~Region.Label, 
                        stratification="geographical")
print(dht2.LTM.output, report="Density")

# extract estimates
LTM.grp.results <- LTM.df.hn.size$dht$clusters$N[1:7, ]
LTM.ind.results <- LTM.df.hn.size$dht$individuals$N[1:7, ]

# bind results
LTM.results <- rbind(LTM.grp.results, LTM.ind.results)
LTM.results <- LTM.results %>% rename(Year = Label) %>% 
                mutate(Label = rep(c("Grp", "Ind"), each=7)) %>% 
                mutate(Species = rep("LTM", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hn", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("size", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)



LTM_final_plot <- ggplot(LTM.results[LTM.results$Label=="Grp",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(LTM.results, "Output/Results/LTM_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/LTM_final_plot.png", LTM_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")



### binned

summary(LTM.df.hn.strat.size.bin)

# extract estimates
LTM.grp.abund.bin   <- LTM.df.hn.strat.size.bin$dht$clusters$N[1:7,1:6]
LTM.ind.abund.bin   <- LTM.df.hn.strat.size.bin$dht$individuals$N[1:7,1:6]
LTM.grp.density.bin <- LTM.df.hn.strat.size.bin$dht$clusters$D[1:7,1:6]
LTM.ind.density.bin <- LTM.df.hn.strat.size.bin$dht$individuals$D[1:7,1:6]

LTM.grp.abund.bin <- LTM.grp.abund.bin %>% 
                      rename(Year=Label,N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl) %>% 
                      mutate(Label="Grp")
LTM.ind.abund.bin <- LTM.ind.abund.bin %>% 
                      rename(Year=Label,N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl) %>% 
                      mutate(Label="Ind")
LTM.grp.density.bin <- LTM.grp.density.bin %>% 
                      rename(Year=Label,D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl) %>% 
                      mutate(Label="Grp")
LTM.ind.density.bin <- LTM.ind.density.bin %>% 
                      rename(Year=Label,D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl) %>% 
                      mutate(Label="Ind")


# merge 
LTM.N.bin <- rbind(LTM.grp.abund.bin,LTM.ind.abund.bin)
LTM.N.bin <- LTM.N.bin %>% arrange(Label,Year)
LTM.D.bin <- rbind(LTM.grp.density.bin,LTM.ind.density.bin)
LTM.D.bin <- LTM.D.bin %>% arrange(Label, Year)
LTM.D.bin <- LTM.D.bin[, -c(1,7)]
LTM.results.final <- cbind(LTM.N.bin,LTM.D.bin)


# bind results

LTM.results.final <- LTM.results.final %>% 
                      mutate(Species = rep("LTM", times=14)) %>% 
                      mutate(DetFun = rep("pooled", times=14)) %>% 
                      mutate(Key = rep("Hn", times=14)) %>% 
                      mutate(Adjust = rep("NA", times=14)) %>% 
                      mutate(Covar = rep("size + year", times=14)) %>% 
                      select(Year,Species,DetFun,Key,Adjust,Covar,Label,N,n_se,n_cv,n_lcl,n_ucl,
                             D,d_se,d_cv,d_lcl,d_ucl)

# convert density estimates from m2 to km2
LTM.results.final <- LTM.results.final %>% 
                      mutate(D = D*1000000) %>% 
                      mutate(d_se = d_se*1000000) %>% 
                      mutate(d_lcl = d_lcl*1000000) %>% 
                      mutate(d_ucl = d_ucl*1000000)



LTM_final_plot_binned <- ggplot(LTM.results.bin[LTM.results.bin$Label=="Grp",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(LTM.results.final, "Output/Results/LTM_results_final_binned.csv")
ggsave("Output/Results/Plots/Point_estimates/LTM_final_plot_binned.png", LTM_final_plot_binned, 
       dpi=300, width = 20, height = 20, units = "cm")


### extract density
LTM.grp.density.bin <- LTM.df.hn.strat.size.bin$dht$clusters$D[1:7, ]
LTM.grp.density.bin <- LTM.grp.density.bin %>% rename(Year=Label) %>% 
                        mutate(Species="LTM") %>% 
                        select(Species,Year,Estimate,se,cv,lcl,ucl,df)

# compare original to binned

# load in original
LTM.results <- read.csv("Output/Results/LTM_results_final.csv")
LTM.results <- LTM.results[ ,-1]
LTM.results$analysis <- "original"

LTM.results.bin$analysis <- "binned"
LTM.results.comp <- rbind(LTM.results,LTM.results.bin)

LTM_analysis_comparison_plot <- ggplot(LTM.results.comp[LTM.results.comp$Label=="Grp",], 
                                       aes(x=Year, y=Estimate, colour=analysis))+
                  geom_point(shape=16, size=2, position=position_dodge(width=0.3))+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3, position = position_dodge(width=0.3))+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))
ggsave("Output/Results/Plots/Point_estimates/LTM_analysis_comparison_plot.png", 
       LTM_analysis_comparison_plot, dpi=300, width = 20, height = 20, units = "cm")



#### Pig-tailed macaque ######################################################
  ## Subset data ####

### Previously, effort strata were ignored for this species, as they were for the other species. After discussions with Olly, we decided to check the effect of accounting for the effort strata for PTM in 2011, 2014, and 2016, as the estimates for these years (well, 2011 and 2014 anyway) were very high compared to other years, suggesting the strata were having an effect. I tested a couple of different approaches (see "Analytical_approach_test_CDS.R", lines 2190 onwards).  For PTM, the approach we have settled on is to split the stratified and un-stratified years, and analyse them separately. This accounts for the effect of the strata, and removes any effects in the unstratified years.

# This means that the old data will be loaded and used for 2011, 2014, and 2016, but the new data will be used for all of the unstratified years.  

# note - I still need to remove T20, change the survey area, and add the stratum variable


# load old data (2010-2018) that have strata
load("./Output/Data/Archive/KSWS_MASTER.Rdata")

allData$obs.habitat <- as.factor(allData$obs.habitat)
allData$obs.observer <- as.factor(allData$obs.observer)

## Remove T20
allData <- allData[allData$obs.transect != 20, ]
sample.table <- sample.table[sample.table$Sample.Label != 20,]
obs.table <- obs.table[obs.table$Sample.Label != 20,]


# subset PTM data
PTM.data <- allData[allData$species=="PTM",] 


# subset data for 2011, 2014, 2016
PTM.data.strat <- as.data.frame(PTM.data[PTM.data$stratum=="2011_H" | PTM.data$stratum=="2011_L" 
                                         | PTM.data$stratum=="2014_H" | PTM.data$stratum=="2014_L"
                                         | PTM.data$stratum=="2016_H" | PTM.data$stratum=="2016_L", ]) 


region.table.strat <- as.data.frame(full.region.table[full.region.table$Region.Label=="2011_H"  |
                                                        full.region.table$Region.Label=="2011_L" |
                                                        full.region.table$Region.Label=="2014_H" |
                                                        full.region.table$Region.Label=="2014_L" |
                                                        full.region.table$Region.Label=="2016_H" |
                                                        full.region.table$Region.Label=="2016_L",])

# change survey area
new.area <- 1880000000/2
region.table.strat$Area <- new.area


# create region tables where other strata are set to 0
region.table.11 <- region.table.strat
region.table.11$Area[region.table.11$Region.Label != "2011_L" & region.table.11$Region.Label != "2011_H"] <- 0

region.table.14 <- region.table.strat
region.table.14$Area[region.table.14$Region.Label != "2014_L" & region.table.14$Region.Label != "2014_H"] <- 0

region.table.16 <- region.table.strat
region.table.16$Area[region.table.16$Region.Label != "2016_L" & region.table.16$Region.Label != "2016_H"] <- 0


sample.table.strat <- as.data.frame(sample.table[sample.table$Region.Label =="2011_H" | 
                                                   sample.table$Region.Label =="2011_L" |
                                                   sample.table$Region.Label =="2014_H" |
                                                   sample.table$Region.Label =="2014_L" |
                                                   sample.table$Region.Label =="2016_H" |
                                                   sample.table$Region.Label =="2016_L",])

obs.table.strat <- as.data.frame(obs.table[obs.table$Region.Label=="2011_H" | obs.table$Region.Label=="2011_L" |
                                             obs.table$Region.Label=="2014_H" | obs.table$Region.Label=="2014_L" |
                                             obs.table$Region.Label=="2016_H" | obs.table$Region.Label=="2016_L",])



### now using the updated data, I will get the 2010, 2013, 2018, 2020 data

# load the latest data
load("./Output/Data/KSWS_MASTER.Rdata")

## Remove T20
allData <- allData[allData$obs.transect != 20, ]
sample.table <- sample.table[sample.table$Sample.Label != 20,]
obs.table <- obs.table[obs.table$Sample.Label != 20,]


# subset PTM data
PTM.data <- as.data.frame(allData[allData$species=="PTM",]) 


PTM.data.nostrat <- PTM.data[PTM.data$stratum==2010 | PTM.data$stratum==2013 | PTM.data$stratum==2018 |
                               PTM.data$stratum==2020, ]

region.table.nostrat <- as.data.frame(full.region.table[full.region.table$Region.Label=="2010" | 
                                                          full.region.table$Region.Label=="2013" |
                                                          full.region.table$Region.Label=="2018" |
                                                          full.region.table$Region.Label=="2020", ])

# change survey area
new.area <- 1880000000
region.table.nostrat$Area <- new.area


sample.table.nostrat <- as.data.frame(
  sample.table[sample.table$Region.Label==2010 | sample.table$Region.Label==2013 |
                 sample.table$Region.Label==2018 | sample.table$Region.Label==2020, ])

obs.table.nostrat <- as.data.frame(obs.table[obs.table$Region.Label==2010 | obs.table$Region.Label==2013 |
                                               obs.table$Region.Label==2018 | obs.table$Region.Label==2020, ])


# create scaled continuous stratum variable for the detection functions
PTM.data.nostrat$stratum <- as.vector(scale(PTM.data.nostrat$stratum, center = T, scale = T)) 
PTM.data.strat$year <- ifelse(PTM.data.strat$stratum=="2011_L"|PTM.data.strat$stratum=="2011_H",2011,
                              ifelse(PTM.data.strat$stratum=="2014_L"|PTM.data.strat$stratum=="2014_H",2014,
                                     ifelse(PTM.data.strat$stratum=="2016_L"|PTM.data.strat$stratum=="2016_H",2016,NA)))
PTM.data.strat$stratum <- as.vector(scale(PTM.data.strat$year, center = T, scale = T)) 


  ## New analysis - splitting years by strata ####

# this analysis is from Sept 2020 and is aimed at accounting for strata in 2011, 2014, and 2016.  For the original analysis, go to the sections below starting at "Exploritory plots and linear models"

    # Stratified years ####

par(mfrow=c(1,2))
# histograms of the data
hist(PTM.data.strat$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(PTM.data.strat$distance, main=NULL, xlab="Distance (m)", breaks=c(40))

# set truncation distance
ptm.strat.trunc <- 60


### run DF models

## Primary models

# Uni cosine
PTM.df.uni.cos <- ds(data=PTM.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ptm.strat.trunc, key="uni", adjustment= "cos")

# Uni poly
PTM.df.uni.poly <- ds(data=PTM.data.strat, region.table=region.table.strat, 
                      sample.table=sample.table.strat, obs.table=obs.table.strat, 
                      truncation=ptm.strat.trunc, key="uni", adjustment= "poly")

# HN cosine
PTM.df.hn.cos <- ds(data=PTM.data.strat, region.table=region.table.strat, 
                    sample.table=sample.table.strat, obs.table=obs.table.strat, 
                    truncation=ptm.strat.trunc, key="hn", adjustment= "cos")


# HN hermite
PTM.df.hn.herm <- ds(data=PTM.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ptm.strat.trunc, key="hn", adjustment= "herm")


# HR cosine
PTM.df.hr.cos <- ds(data=PTM.data.strat, region.table=region.table.strat, 
                    sample.table=sample.table.strat, obs.table=obs.table.strat, 
                    truncation=ptm.strat.trunc, key="hr", adjustment= "cos")

# HR poly
PTM.df.hr.poly <- ds(data=PTM.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ptm.strat.trunc, key="hr", adjustment= "poly")

# compare models
ptm.df.prim.comp <- summarize_ds_models(PTM.df.uni.cos, PTM.df.uni.poly,
                                        PTM.df.hn.cos,  
                                        PTM.df.hr.cos,  
                                        output = "plain")
ptm.df.prim.comp[ ,1:5]
ptm.df.prim.comp[ ,6:7]

# HR and Uni cos have most support. Based on the data, the original analysis, and the species, I was expecting HR to be a solid choice

# plot them both just to check
plot(PTM.df.uni.cos, main="PTM.df.uni.cos")
plot(PTM.df.hr.cos, main="PTM.df.hr.cos")

# HR selected


### covariate models

# cluster size
PTM.df.hr.size <- ds(data=PTM.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ptm.strat.trunc, key="hr", formula=~size)

# observer
PTM.df.hr.observer <- ds(data=PTM.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ptm.strat.trunc, key="hr", formula=~obs.observer)

# stratum
PTM.df.hr.stratum <- ds(data=PTM.data.strat, region.table=region.table.strat, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ptm.strat.trunc, key="hr", formula=~stratum)

# compare models
ptm.df.cov.comp <- summarize_ds_models(PTM.df.hr.size, PTM.df.hr.observer, PTM.df.hr.stratum,
                                       output = "plain")
ptm.df.cov.comp[ ,1:5]
ptm.df.cov.comp[ ,6:7]

# size model has most support


## compare primary with covariate model
ptm.df.final.comp <- summarize_ds_models(PTM.df.hr.size, PTM.df.hr.cos,
                                         output = "plain")
ptm.df.final.comp[ ,1:5]
ptm.df.final.comp[ ,6:7]
# covariate model has most support. Cluster size is also a reasonable covariate for this species. 


## run final models to get estimates for each year (areas set to 0)

# 2011
PTM.11.hr.size <- ds(data=PTM.data.strat, region.table=region.table.11, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ptm.strat.trunc, key="hr", formula=~size)

# 2014
PTM.14.hr.size <- ds(data=PTM.data.strat, region.table=region.table.14, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ptm.strat.trunc, key="hr", formula=~size)

# 2016
PTM.16.hr.size <- ds(data=PTM.data.strat, region.table=region.table.16, 
                     sample.table=sample.table.strat, obs.table=obs.table.strat, 
                     truncation=ptm.strat.trunc, key="hr", formula=~size)


    # Unstratified years ####

par(mfrow=c(1,2))
# histograms of the data
hist(PTM.data.nostrat$distance, main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(PTM.data.nostrat$distance, main=NULL, xlab="Distance (m)", breaks=c(40))

# set truncation distance
ptm.nostrat.trunc <- 60


### run DF models

## Primary models

# Uni cosine
PTM.df.uni.cos <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                     sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                     truncation=ptm.nostrat.trunc, key="uni", adjustment= "cos")

# Uni poly
PTM.df.uni.poly <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                      sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                      truncation=ptm.nostrat.trunc, key="uni", adjustment= "poly")

# HN cosine
PTM.df.hn.cos <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                    sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                    truncation=ptm.nostrat.trunc, key="hn", adjustment= "cos")


# HN hermite
PTM.df.hn.herm <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                     sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                     truncation=ptm.nostrat.trunc, key="hn", adjustment= "herm")


# HR cosine
PTM.df.hr.cos <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                    sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                    truncation=ptm.nostrat.trunc, key="hr", adjustment= "cos")

# HR poly
PTM.df.hr.poly <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                     sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                     truncation=ptm.nostrat.trunc, key="hr", adjustment= "poly")

# compare models
PTM.df.prim.comp <- summarize_ds_models(PTM.df.uni.cos, PTM.df.uni.poly,
                                        PTM.df.hn.cos,  
                                        PTM.df.hr.cos, PTM.df.hr.poly,  
                                        output = "plain")
PTM.df.prim.comp[ ,1:5]
PTM.df.prim.comp[ ,6:7]
# Hn cos and Uni cos have most support

# plot
plot(PTM.df.uni.cos, main="PTM.df.uni.cos")
plot(PTM.df.hn.cos, main="PTM.df.hn.cos")
# very similar, barely any difference. I want to test covariates so HN selected


### models with covariates

# cluster size
PTM.df.hn.size <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                     sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                     truncation=ptm.nostrat.trunc, key="hn", formula=~size)

# observer
PTM.df.hn.observer <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                         sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                         truncation=ptm.nostrat.trunc, key="hn", formula=~obs.observer)

# stratum
PTM.df.hn.stratum <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                         sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                         truncation=ptm.nostrat.trunc, key="hn", formula=~stratum)

# size + stratum
PTM.df.hn.strat.size <- ds(data=PTM.data.nostrat, region.table=region.table.nostrat, 
                        sample.table=sample.table.nostrat, obs.table=obs.table.nostrat, 
                        truncation=ptm.nostrat.trunc, key="hn", formula=~stratum + size)

# compare models
ptm.df.cov.comp <- summarize_ds_models(PTM.df.hn.size, PTM.df.hn.observer, PTM.df.hn.stratum,
                                       PTM.df.hn.strat.size, output = "plain")
ptm.df.cov.comp[ ,1:5]
ptm.df.cov.comp[ ,6:7]
# model with stratum and size has significant support and is selected

# skip down to "PTM Final results" for the results


  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

# distance histograms
hist(PTM.data$distance, main="All years", xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(PTM.data$distance, main="All years", xlab="Distance (m)", breaks=c(40))

# interesting shape. plot all years separately to see if it is a consistent shape

par(mfrow=c(2,4))
hist(PTM.data$distance[PTM.data$year=="2010"], main="2010")
hist(PTM.data$distance[PTM.data$year=="2011"], main="2011")
hist(PTM.data$distance[PTM.data$year=="2013"], main="2013")
hist(PTM.data$distance[PTM.data$year=="2014"], main="2014")
hist(PTM.data$distance[PTM.data$year=="2016"], main="2016")
hist(PTM.data$distance[PTM.data$year=="2018"], main="2018")
hist(PTM.data$distance[PTM.data$year=="2020"], main="2020")
# Interesting shapes. The shape is consitent apart from 2011 and 2020 really.  It looks like detection is very likely between 0 and 35-40m, and then detection drops right off.  I wonder if this is a behavioural thng?  Perhaps they are able to stay silent and/or hide if they are far enough away, but if they are close to the observer they will flee loudly?
# After checking with the field teams, this shape is due to the species ecology and behaviour. The species lives in large, noisy groups, which are easy to detect. They call loudly in the forest, especially when foraging. However, they are only semi-arboreal, and will actually spend most of their time foraging on the ground. They wil also flee on the ground, rather than in the trees. And because they live in EG and SEG forest, if they are more than 30-40m away you will not be able to see them for a positive ID. 


### update - binning will be done for this species due to clumping around 0m

# find best bins
hist(PTM.data$distance[PTM.data$distance<60], breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60))
hist(PTM.data$distance[PTM.data$distance<60], breaks=c(0,5,10,15,20,25,33,45,60))
hist(PTM.data$distance[PTM.data$distance<60], breaks=c(0,10,20,30,45,60))
hist(PTM.data$distance[PTM.data$distance<60], breaks=c(0,7,15,25,33,39,52,60)) # try this one

hist(PTM.data$distance[PTM.data$distance<50], breaks=c(0,7,15,25,33,50))
hist(PTM.data$distance[PTM.data$distance<50], breaks=c(0,5,9,15,20,25,30,40,50))

hist(PTM.data$distance[PTM.data$distance<60], breaks=c(0,15,25,35,45,60)) # try this one with HR


# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.PTM <- 50
trunc.PTM.bin <- 60

# Count the number of observations discarded
nrow(PTM.data[PTM.data$distance>trunc.PTM,]) 

length(PTM.data$distance)

29/401*100 # 7.2%

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(1,2))
plot(PTM.data$size, PTM.data$distance, main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.PTM <- lm(distance~size, data=PTM.data)
lines(PTM.data$size, as.vector(predict(lm.PTM, PTM.data)))
summary(lm.LTM)
# some evidence of size bias although model is not strong 

# Plot of Observer factor against distance
plot(PTM.data$obs.observer,PTM.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# not as much variation between observers than some of the other species (apart from Samart).

# I can't include habitat or AMPM in the DF because those variables weren't recorded in earlier years


  ## Fit a detection function ####
    # Uniform ####

### unbinned


# Uni cosine

PTM.df.uni.cos <- ds(data=PTM.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.PTM, key="uni", adjustment= "cos")
summary(PTM.df.uni.cos)
# cosine(1)  

# Uni poly
PTM.df.uni.poly <- ds(data=PTM.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,  
                      truncation=trunc.PTM, key="uni", adjustment= "poly")
summary(PTM.df.uni.poly)
# poly(2,4)

## compare uni models
ptm.uni.comp <- summarize_ds_models(PTM.df.uni.cos, PTM.df.uni.poly,
                                    output="plain")
ptm.uni.comp[ ,1:4]
ptm.uni.comp[ ,5:7]
# cos is the preferred model but poly dAIC < 2

# plot the fits
par(mfrow=c(2,2))

# cos
plot(PTM.df.uni.cos, main = "PTM.df.uni.cos")

covar.fit <- ddf.gof(PTM.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(PTM.df.uni.poly, main = "PTM.df.uni.poly")

covar.fit <- ddf.gof(PTM.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# Virtually identical fits.

# cos selected


### binned


# Uni cosine
PTM.df.uni.cos.bin <- ds(data=PTM.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                     truncation=trunc.PTM.bin, key="uni", adjustment= "cos",
                     cutpoints = c(0,7,15,25,33,39,52,60))
summary(PTM.df.uni.cos.bin)
# cosine(1)

# Uni poly
PTM.df.uni.poly.bin <- ds(data=PTM.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table, 
                         truncation=trunc.PTM.bin, key="uni", adjustment= "poly",
                         cutpoints = c(0,7,15,25,33,39,52,60))
summary(PTM.df.uni.poly.bin)
# poly(2,4)

## compare uni models
ptm.uni.comp.bin <- summarize_ds_models(PTM.df.uni.cos.bin, PTM.df.uni.poly.bin,
                                    output="plain")
ptm.uni.comp.bin[ ,1:4]
ptm.uni.comp.bin[ ,5:7]
# cos perferred

# plot both
par(mfrow=c(1,2))
plot(PTM.df.uni.cos.bin, main="cos")
plot(PTM.df.uni.poly.bin, main="poly")
# don't look like great fits, and look pretty much identical so will go with AIC

# PTM.df.uni.cos.bin selected



    # Half normal ####


## unbinned

# HN cosine
PTM.df.hn.cos <- ds(data=PTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.PTM, key="hn", adjustment= "cos")
summary(PTM.df.hn.cos)
# no adjustment selected

# HN hermite
PTM.df.hn.herm <- ds(data=PTM.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.PTM, key="hn", adjustment= "herm")
summary(PTM.df.hn.herm)
# no adjustment selected


# plot the fit
par(mfrow=c(1,2))

# cos
plot(PTM.df.hn.cos, main = "PTM.df.hn.cos")

covar.fit <- ddf.gof(PTM.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM
# p=0.5 around the same distance as uni models (30m)

# cos is selected


### binned

# HN cosine
PTM.df.hn.cos.bin <- ds(data=PTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                    truncation=trunc.PTM.bin, key="hn", adjustment= "cos",
                    cutpoints = c(0,7,15,25,33,39,52,60))
summary(PTM.df.hn.cos.bin)
# no adjustment selected

# HN herm
PTM.df.hn.herm.bin <- ds(data=PTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.PTM.bin, key="hn", adjustment= "herm",
                        cutpoints = c(0,7,15,25,33,39,52,60))
summary(PTM.df.hn.herm.bin)
# no adjustment selected

plot(PTM.df.hn.cos.bin)
# also not a brilliant fit




    # Hazard rate ####


### unbinned

# HR cosine
PTM.df.hr.cos <- ds(data=PTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.PTM, key="hr", adjustment= "cos")
summary(PTM.df.hr.cos)
# no adjustment selected


# HR poly
PTM.df.hr.poly <- ds(data=PTM.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.PTM, key="hr", adjustment= "poly")
summary(PTM.df.hr.poly)
# no adjustment selected 


# plot the fit
par(mfrow=c(1,2))

# hr cos
plot(PTM.df.hr.cos, main = "PTM.df.hr.cos")

covar.fit <- ddf.gof(PTM.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# Probably a fairly realistic fit - p=1 out to 20m, then drops sharply off


### binned

# HR cosine
PTM.df.hr.cos.bin <- ds(data=PTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.PTM.bin, key="hr", adjustment= "cos",
                        cutpoints = c(0,7,15,25,33,39,52,60))
summary(PTM.df.hr.cos.bin)
# no adjustment selected

# HR poly
PTM.df.hr.poly.bin <- ds(data=PTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.PTM.bin, key="hr", adjustment= "poly",
                        cutpoints = c(0,7,15,25,33,39,52,60))
summary(PTM.df.hr.poly.bin)
# no adjustment selected

plot(PTM.df.hr.cos.bin)
# I wonder whether I re-jig the bins to fit a HR model

# try different cut points
# HR cosine
PTM.df.hr.cos.bin2 <- ds(data=PTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.PTM.bin, key="hr", adjustment= "cos",
                        cutpoints = c(0,15,25,35,45,60))
summary(PTM.df.hr.cos.bin)

plot(PTM.df.hr.cos.bin, main="1")
plot(PTM.df.hr.cos.bin2, main="2")

# check estimates between the different bins
ltm.hr.bin.comp <- data.frame(model = rep(c("1", "2"), each=7),
                                year = c(PTM.df.hr.cos.bin$dht$individuals$N$Label[1:7],
                                         PTM.df.hr.cos.bin2$dht$individuals$N$Label[1:7]),
                                estimate = c(PTM.df.hr.cos.bin$dht$individuals$N$Estimate[1:7],
                                             PTM.df.hr.cos.bin2$dht$individuals$N$Estimate[1:7]),
                                cv = c(PTM.df.hr.cos.bin$dht$individuals$N$cv[1:7],
                                       PTM.df.hr.cos.bin2$dht$individuals$N$cv[1:7]),
                                se = c(PTM.df.hr.cos.bin$dht$individuals$N$se[1:7],
                                       PTM.df.hr.cos.bin2$dht$individuals$N$se[1:7]),
                                lcl = c(PTM.df.hr.cos.bin$dht$individuals$N$lcl[1:7],
                                        PTM.df.hr.cos.bin2$dht$individuals$N$lcl[1:7]),
                                ucl = c(PTM.df.hr.cos.bin$dht$individuals$N$ucl[1:7],
                                        PTM.df.hr.cos.bin2$dht$individuals$N$ucl[1:7]))

ltm.hr.bin.comp

# plot cv
ggplot(ltm.hr.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# virtually identical CV

# plot estimates
ggplot(ltm.hr.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))

# they are exactly the same, so I will continue with the original bins


    # Compare primary models ####

### unbinned

ptm.df.prim.comp <- summarize_ds_models(PTM.df.uni.cos, PTM.df.hn.cos,PTM.df.hr.cos, 
                                           output = "plain")
ptm.df.prim.comp[ ,1:5]
ptm.df.prim.comp[ ,6:7]

# All models have good support (dAIC < 0.5). 

# plot all fits together
par(mfrow=c(3,2))

# uni poly
plot(PTM.df.uni.cos, main = "PTM.df.uni.cos")

covar.fit <- ddf.gof(PTM.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HN cos
plot(PTM.df.hn.cos, main = "PTM.df.hn.cos")

covar.fit <- ddf.gof(PTM.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HR cos
plot(PTM.df.hr.cos, main = "PTM.df.hr.cos")

covar.fit <- ddf.gof(PTM.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# I mean, HR really does fit the data better.  And this shape of data seems to be consistent across years, and therefore I am inclined to believe it. The HR QQ plot is also better than the other two. I selected HR previously (before including 2020), and I think I still agree with my past self. 

# PTM.df.hr.cos selected


### binned

ptm.df.prim.comp.bin <- summarize_ds_models(PTM.df.uni.cos.bin, PTM.df.hn.cos.bin,PTM.df.hr.cos.bin, 
                                        output = "plain")
ptm.df.prim.comp.bin[ ,1:5]
ptm.df.prim.comp.bin[ ,6:7]
# all models have some support

# plot together
par(mfrow=c(2,2))
plot(PTM.df.uni.cos.bin, main="uni")
plot(PTM.df.hn.cos.bin, main="hn")
plot(PTM.df.hr.cos.bin, main="hr")
# uni and hn are very similar. I think HR is the most realistic


# check estimates
ptm.prim.bin.comp <- data.frame(model = rep(c("uni", "hn", "hr"), each=7),
                                year = c(PTM.df.uni.cos.bin$dht$individuals$N$Label[1:7],
                                         PTM.df.hn.cos.bin$dht$individuals$N$Label[1:7],
                                         PTM.df.hr.cos.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(PTM.df.uni.cos.bin$dht$individuals$N$Estimate[1:7],
                                             PTM.df.hn.cos.bin$dht$individuals$N$Estimate[1:7],
                                             PTM.df.hr.cos.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(PTM.df.uni.cos.bin$dht$individuals$N$cv[1:7],
                                       PTM.df.hn.cos.bin$dht$individuals$N$cv[1:7],
                                       PTM.df.hr.cos.bin$dht$individuals$N$cv[1:7]),
                                se = c(PTM.df.uni.cos.bin$dht$individuals$N$se[1:7],
                                       PTM.df.hn.cos.bin$dht$individuals$N$se[1:7],
                                       PTM.df.hr.cos.bin$dht$individuals$N$se[1:7]),
                                lcl = c(PTM.df.uni.cos.bin$dht$individuals$N$lcl[1:7],
                                        PTM.df.hn.cos.bin$dht$individuals$N$lcl[1:7],
                                        PTM.df.hr.cos.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(PTM.df.uni.cos.bin$dht$individuals$N$ucl[1:7],
                                        PTM.df.hn.cos.bin$dht$individuals$N$ucl[1:7],
                                        PTM.df.hr.cos.bin$dht$individuals$N$ucl[1:7]))

ptm.prim.bin.comp

# plot cv
ggplot(ptm.prim.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# HR always the highest CV, uni always the lowest.  but they are all very similar

# plot estimates
ggplot(ptm.prim.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# Hr always predicting the lowest, hn always the highest. Fairly consistent variation, and the estimates are not wildly different. Hr is also the most precise. As mentioned above, I think HR reflects the ecology of the species the best

# PTM.df.hr.cos.bin selected


    # Models with harsh truncation ####


# set trunc distance 
ptm.trunc.harsh <- 35
ptm.trunc.harsh.bin <- 39


### unbinned

par(mfrow=c(1,2))

# HR cos harsh
PTM.df.hr.cos.harsh <- ds(data=PTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=ptm.trunc.harsh, key="hr", adjustment="cos")
summary(PTM.df.hr.cos.harsh)
# no adjustment selected

# compare original with harsh and lenient truncation
ddf.gof(PTM.df.hr.cos.harsh$ddf, main = "LTM.df.hn.cos.harsh", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(PTM.df.hr.cos$ddf, main = "PTM.df.hr.cos",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# Original is a better fit


# plot all the fits
par(mfrow=c(2,2))

# original
plot(PTM.df.hr.cos, main = "PTM.df.hr.cos")

covar.fit <- ddf.gof(PTM.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# harsh
plot(PTM.df.hr.cos.harsh, main = "PTM.df.hr.cos.harsh")

covar.fit <- ddf.gof(PTM.df.hr.cos.harsh$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

    
# PTM.df.hr.cos selected


### binned


PTM.df.hr.cos.bin.harsh <- ds(data=PTM.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=ptm.trunc.harsh.bin, key="hr", adjustment= "cos",
                        cutpoints = c(0,7,15,25,33,39))
summary(PTM.df.hr.cos.bin.harsh)

# plot together
par(mfrow=c(1,2))
plot(PTM.df.hr.cos.bin.harsh, main="harsh")
plot(PTM.df.hr.cos.bin, main="orig")

# compare estiamtes
ltm.harsh.bin.comp <- data.frame(model = rep(c("orig", "harsh"), each=7),
                              year = c(PTM.df.hr.cos.bin$dht$individuals$N$Label[1:7],
                                       PTM.df.hr.cos.bin.harsh$dht$individuals$N$Label[1:7]),
                              estimate = c(PTM.df.hr.cos.bin$dht$individuals$N$Estimate[1:7],
                                           PTM.df.hr.cos.bin.harsh$dht$individuals$N$Estimate[1:7]),
                              cv = c(PTM.df.hr.cos.bin$dht$individuals$N$cv[1:7],
                                     PTM.df.hr.cos.bin.harsh$dht$individuals$N$cv[1:7]),
                              se = c(PTM.df.hr.cos.bin$dht$individuals$N$se[1:7],
                                     PTM.df.hr.cos.bin.harsh$dht$individuals$N$se[1:7]),
                              lcl = c(PTM.df.hr.cos.bin$dht$individuals$N$lcl[1:7],
                                      PTM.df.hr.cos.bin.harsh$dht$individuals$N$lcl[1:7]),
                              ucl = c(PTM.df.hr.cos.bin$dht$individuals$N$ucl[1:7],
                                      PTM.df.hr.cos.bin.harsh$dht$individuals$N$ucl[1:7]))

ltm.harsh.bin.comp

# plot cv
ggplot(ltm.harsh.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# a lot of variation in cv. harsh trunc is lower in 2011 and 2013 but original is lower in other years

# plot estimates
ggplot(ltm.harsh.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# the model with harsh truncation appears to improve precision in most years. And the estimates are not very different, and inferences are the same

# PTM.df.hr.cos.bin.harsh selected


    # Models with covariates ####

### unbinned

## cluster size
PTM.df.hr.size <- ds(data=PTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.PTM, key="hr", formula = ~size)
summary(PTM.df.hr.size)

# plot
par(mfrow=c(1,2))
plot(PTM.df.hr.size, main = "PTM.df.hr.size")

covar.fit <- ddf.gof(PTM.df.hr.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



## Observer
PTM.df.hr.obs <- ds(data=PTM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.PTM, key="hr", formula = ~obs.observer)
summary(PTM.df.hr.obs)

# plot
plot(PTM.df.hr.obs, main = "PTM.df.hr.obs")

covar.fit <- ddf.gof(PTM.df.hr.obs$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# lot of variation

## stratum (year)
PTM.df.hr.strat <- ds(data=PTM.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.PTM, key="hr", formula = ~stratum)
summary(PTM.df.hr.strat)

# plot
plot(PTM.df.hr.strat, main = "PTM.df.hr.strat")

covar.fit <- ddf.gof(PTM.df.hr.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# less variation

## stratum (year) + size
PTM.df.hr.strat.size <- ds(data=PTM.data, region.table=full.region.table, 
                           sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.PTM, key="hr", formula = ~stratum+size)
summary(PTM.df.hr.strat.size)

# plot
plot(PTM.df.hr.strat.size, main = "PTM.df.hr.strat.size")

covar.fit <- ddf.gof(PTM.df.hr.strat.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)


### binned

# size
PTM.df.hr.size.bin.harsh <- ds(data=PTM.data, region.table=full.region.table, 
                              sample.table=sample.table, obs.table=obs.table, 
                              truncation=ptm.trunc.harsh.bin, key="hr", formula=~size,
                              cutpoints = c(0,7,15,25,33,39))
summary(PTM.df.hr.size.bin.harsh)

# observer
PTM.df.hr.obs.bin.harsh <- ds(data=PTM.data, region.table=full.region.table, 
                               sample.table=sample.table, obs.table=obs.table, 
                               truncation=ptm.trunc.harsh.bin, key="hr", formula=~obs.observer,
                               cutpoints = c(0,7,15,25,33,39))
summary(PTM.df.hr.obs.bin.harsh)

# stratum
PTM.df.hr.strat.bin.harsh <- ds(data=PTM.data, region.table=full.region.table, 
                              sample.table=sample.table, obs.table=obs.table, 
                              truncation=ptm.trunc.harsh.bin, key="hr", formula=~stratum,
                              cutpoints = c(0,7,15,25,33,39))
summary(PTM.df.hr.strat.bin.harsh)



    # Compare covariate models ####


### unbinned

ptm.df.cov.comp <- summarize_ds_models(PTM.df.hr.size, PTM.df.hr.obs,PTM.df.hr.strat.size,
                                          PTM.df.hr.strat,output = "plain")

ptm.df.cov.comp[ ,1:5]
ptm.df.cov.comp[ ,6:7]

summary(PTM.df.hr.strat)
summary(PTM.df.hr.strat.size)

# Model with size has the most support. obs model has dAIC > 7 so I will exclude. I added another model with size and stratum. it has some support.  

# extract and plot to see what the different models say
size <- data.frame(label = "size",
                   estimate = PTM.df.hr.size$dht$clusters$N$Estimate[c(1:7)],
                   lcl = PTM.df.hr.size$dht$clusters$N$lcl[c(1:7)],
                   ucl = PTM.df.hr.size$dht$clusters$N$ucl[c(1:7)])

strat <- data.frame(label = "strat",
                    estimate = PTM.df.hr.strat$dht$clusters$N$Estimate[c(1:7)],
                    lcl = PTM.df.hr.strat$dht$clusters$N$lcl[c(1:7)],
                    ucl = PTM.df.hr.strat$dht$clusters$N$ucl[c(1:7)])

strat.size <- data.frame(label = "strat.size",
                    estimate = PTM.df.hr.strat.size$dht$clusters$N$Estimate[c(1:7)],
                    lcl = PTM.df.hr.strat.size$dht$clusters$N$lcl[c(1:7)],
                    ucl = PTM.df.hr.strat.size$dht$clusters$N$ucl[c(1:7)])


mod.check <- rbind(size,strat,strat.size)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=3)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# estimates from the three models are all very similar, so the take home message is that it doesn' really matter!  Seeing as the models with the covariates have support and provide estiamtes with similar precision, I see no reason to exclude them - they are there to reduce bias so I'll keep them both.

# PTM.df.hr.strat.size selected


### binned

ptm.df.cov.comp.bin <- summarize_ds_models(PTM.df.hr.size.bin.harsh,PTM.df.hr.obs.bin.harsh,
                                           PTM.df.hr.strat.bin.harsh, output = "plain")

ptm.df.cov.comp.bin[ ,1:5]
ptm.df.cov.comp.bin[ ,6:7]
# strat and observer have no support

# PTM.df.hr.size.bin.harsh selected

    # Compare primary and covariate models ####


### unbinned

ptm.df.final.compare <- summarize_ds_models(PTM.df.hr.cos,PTM.df.hr.strat.size, output = "plain")
ptm.df.final.compare[ ,1:3]
ptm.df.final.compare[ ,4:7]

# Models have similar AIC support. Plot the estimates together

prim <- data.frame(label = "prim",
                    estimate = PTM.df.hr.cos$dht$clusters$N$Estimate[c(1:7)],
                    lcl = PTM.df.hr.cos$dht$clusters$N$lcl[c(1:7)],
                    ucl = PTM.df.hr.cos$dht$clusters$N$ucl[c(1:7)])

strat.size <- data.frame(label = "strat.size",
                         estimate = PTM.df.hr.strat.size$dht$clusters$N$Estimate[c(1:7)],
                         lcl = PTM.df.hr.strat.size$dht$clusters$N$lcl[c(1:7)],
                         ucl = PTM.df.hr.strat.size$dht$clusters$N$ucl[c(1:7)])


mod.check <- rbind(prim,strat.size)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=2)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# both models produce very similar estimates, with similar precision. As there is no good reason to exclude the covariates, I will use that model as the final model


### binned

ptm.df.final.compare.bin <- summarize_ds_models(PTM.df.hr.cos.bin.harsh,PTM.df.hr.size.bin.harsh, 
                                                output = "plain")
ptm.df.final.compare.bin[ ,1:3]
ptm.df.final.compare.bin[ ,4:7]
# covar model preferred

# check estimates 
ltm.final.bin.comp <- data.frame(model = rep(c("orig", "size"), each=7),
                              year = c(PTM.df.hr.cos.bin.harsh$dht$individuals$N$Label[1:7],
                                       PTM.df.hr.size.bin.harsh$dht$individuals$N$Label[1:7]),
                              estimate = c(PTM.df.hr.cos.bin.harsh$dht$individuals$N$Estimate[1:7],
                                           PTM.df.hr.size.bin.harsh$dht$individuals$N$Estimate[1:7]),
                              cv = c(PTM.df.hr.cos.bin.harsh$dht$individuals$N$cv[1:7],
                                     PTM.df.hr.size.bin.harsh$dht$individuals$N$cv[1:7]),
                              se = c(PTM.df.hr.cos.bin.harsh$dht$individuals$N$se[1:7],
                                     PTM.df.hr.size.bin.harsh$dht$individuals$N$se[1:7]),
                              lcl = c(PTM.df.hr.cos.bin.harsh$dht$individuals$N$lcl[1:7],
                                      PTM.df.hr.size.bin.harsh$dht$individuals$N$lcl[1:7]),
                              ucl = c(PTM.df.hr.cos.bin.harsh$dht$individuals$N$ucl[1:7],
                                      PTM.df.hr.size.bin.harsh$dht$individuals$N$ucl[1:7]))

ltm.final.bin.comp

# plot cv
ggplot(ltm.final.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# covar model has lwoer cv across the board

# plot estimates
ggplot(ltm.final.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# very similar estimates, but covar model is more precise


  ## PTM Final Results ####


## These first results are from the latest analysis (Sept 2020) which split the years nto stratified and unstratified.  Below these results I have kept in the code for the original results. 

## extract estimates

# 2011 abundance
PTM.grp.results.11 <- PTM.11.hr.size$dht$clusters$N[7,]
PTM.grp.results.11$Label <- "Grp"
PTM.grp.results.11$Year <- 2011

PTM.ind.results.11 <- PTM.11.hr.size$dht$individuals$N[7,]
PTM.ind.results.11$Label <- "Ind"
PTM.ind.results.11$Year <- 2011

# 2011 density
PTM.grp.density.11 <- PTM.11.hr.size$dht$clusters$D[7,]
PTM.grp.density.11$Year <- "2011"
PTM.grp.density.11 <- PTM.grp.density.11[ ,-1]
PTM.grp.density.11 <- PTM.grp.density.11 %>% select(Year,Estimate,se,cv,lcl,ucl,df)

PTM.ind.density.11 <- PTM.11.hr.size$dht$individuals$D[7,]
PTM.ind.density.11$Year <- "2011"
PTM.ind.density.11 <- PTM.ind.density.11[ ,-1]
PTM.ind.density.11 <- PTM.ind.density.11 %>% select(Year,Estimate,se,cv,lcl,ucl,df)


# 2014 abundance
PTM.grp.results.14 <- PTM.14.hr.size$dht$clusters$N[7,]
PTM.grp.results.14$Label <- "Grp"
PTM.grp.results.14$Year <- 2014

PTM.ind.results.14 <- PTM.14.hr.size$dht$individuals$N[7,]
PTM.ind.results.14$Label <- "Ind"
PTM.ind.results.14$Year <- 2014

# 2014 density
PTM.grp.density.14 <- PTM.14.hr.size$dht$clusters$D[7,]
PTM.grp.density.14$Year <- "2014"
PTM.grp.density.14 <- PTM.grp.density.14[ ,-1]
PTM.grp.density.14 <- PTM.grp.density.14 %>% select(Year,Estimate,se,cv,lcl,ucl,df)

PTM.ind.density.14 <- PTM.14.hr.size$dht$individuals$D[7,]
PTM.ind.density.14$Year <- "2014"
PTM.ind.density.14 <- PTM.ind.density.14[ ,-1]
PTM.ind.density.14 <- PTM.ind.density.14 %>% select(Year,Estimate,se,cv,lcl,ucl,df)

# 2016 abundance 
PTM.grp.results.16 <- PTM.16.hr.size$dht$clusters$N[7,]
PTM.grp.results.16$Label <- "Grp"
PTM.grp.results.16$Year <- 2016

PTM.ind.results.16 <- PTM.16.hr.size$dht$individuals$N[7,]
PTM.ind.results.16$Label <- "Ind"
PTM.ind.results.16$Year <- 2016

# 2016 density
PTM.grp.density.16 <- PTM.16.hr.size$dht$clusters$D[7,]
PTM.grp.density.16$Year <- "2016"
PTM.grp.density.16 <- PTM.grp.density.16[ ,-1]
PTM.grp.density.16 <- PTM.grp.density.16 %>% select(Year,Estimate,se,cv,lcl,ucl,df)

PTM.ind.density.16 <- PTM.16.hr.size$dht$individuals$D[7,]
PTM.ind.density.16$Year <- "2016"
PTM.ind.density.16 <- PTM.ind.density.16[ ,-1]
PTM.ind.density.16 <- PTM.ind.density.16 %>% select(Year,Estimate,se,cv,lcl,ucl,df)


# no strata years abundance
PTM.grp.results.nostrat <- PTM.df.hn.strat.size$dht$clusters$N[1:4,]
PTM.grp.results.nostrat <- PTM.grp.results.nostrat %>% rename(Year = Label) %>% 
  mutate(Label = "Grp") %>% 
  select(Label,Estimate,se,cv,lcl,ucl,df,Year)

PTM.ind.results.nostrat <- PTM.df.hn.strat.size$dht$individuals$N[1:4,]
PTM.ind.results.nostrat <- PTM.ind.results.nostrat %>% rename(Year = Label) %>% 
  mutate(Label = "Ind") %>% 
  select(Label,Estimate,se,cv,lcl,ucl,df,Year)

# no strata years density
PTM.grp.density.nostrat <- PTM.df.hn.strat.size$dht$clusters$D[1:4,]
PTM.grp.density.nostrat <- PTM.grp.density.nostrat %>% rename(Year=Label)
PTM.ind.density.nostrat <- PTM.df.hn.strat.size$dht$individuals$D[1:4,]
PTM.ind.density.nostrat <- PTM.ind.density.nostrat %>% rename(Year=Label)


# bind all abundance together
PTM.results <- rbind(PTM.grp.results.11,PTM.ind.results.11,PTM.grp.results.14,PTM.ind.results.14,
                     PTM.grp.results.16,PTM.ind.results.16,PTM.grp.results.nostrat,
                     PTM.ind.results.nostrat)

PTM.results <- PTM.results %>% arrange(Label,Year)
PTM.results <- PTM.results %>% rename(N = Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl)

# bind all density together
PTM.grp.density <- rbind(PTM.grp.density.11,PTM.grp.density.14,PTM.grp.density.16,PTM.grp.density.nostrat)
PTM.grp.density$Label <- "Grp"
PTM.grp.density <- PTM.grp.density %>% arrange(Year)
PTM.ind.density <- rbind(PTM.ind.density.11,PTM.ind.density.14,PTM.ind.density.16,PTM.ind.density.nostrat)
PTM.ind.density$Label <- "Ind"
PTM.ind.density <- PTM.ind.density %>% arrange(Year)
PTM.density <- rbind(PTM.grp.density,PTM.ind.density)
PTM.density <- PTM.density %>% rename(D = Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl)
PTM.density <- PTM.density[,-c(1,7)]

PTM.results <- cbind(PTM.results,PTM.density)
PTM.results <- PTM.results[,-c(7,14)]


# add other fields and re-order
PTM.results <- PTM.results %>% 
  mutate(Species = rep("PTM", times=14)) %>% 
  mutate(DetFun = rep("pooled", times=14)) %>% 
  mutate(Key = rep("Hn", times=14)) %>% 
  mutate(Adjust = rep("NA", times=14)) %>% 
  mutate(Covar = rep("Cluster size + Year", times=14)) %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,N,n_se,n_cv,n_lcl,n_ucl,D,d_se,d_cv,d_lcl,d_ucl) %>% 
  arrange(Label,Year)

PTM.results$Covar[PTM.results$Year==2011 | PTM.results$Year==2014 | 
                    PTM.results$Year==2016] <- "Cluster size"

PTM.results$Key[PTM.results$Year==2011 | PTM.results$Year==2014 | 
                    PTM.results$Year==2016] <- "Hr"

PTM.results <- PTM.results %>% mutate(D = D*1000000) %>% mutate(d_se = d_se*1000000) %>% 
                                mutate(d_lcl = d_lcl*1000000) %>% mutate(d_ucl = d_ucl*1000000)



PTM_final_plot <- ggplot(PTM.results[PTM.results$Label=="Grp",], aes(x=Year, y= Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  ylim(0,2500)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(PTM.results, "Output/Results/PTM_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/PTM_final_plot.png", PTM_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")





### unbinned

summary(PTM.df.hr.strat.size)

# extract estimates
PTM.grp.results <- PTM.df.hr.strat.size$dht$clusters$N[1:7, ]
PTM.ind.results <- PTM.df.hr.strat.size$dht$individuals$N[1:7, ]

# bind results
PTM.results <- rbind(PTM.grp.results, PTM.ind.results)
PTM.results <- PTM.results %>% rename(Year = Label) %>% 
                mutate(Label = rep(c("Grp", "Ind"), each=7)) %>% 
                mutate(Species = rep("PTM", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hr", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("stratum + size", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)



PTM_final_plot <- ggplot(PTM.results[PTM.results$Label=="Grp",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(PTM.results, "Output/Results/PTM_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/PTM_final_plot.png", PTM_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")



### binned


summary(PTM.df.hr.size.bin.harsh)

# extract estimates
PTM.grp.results.bin <- PTM.df.hr.size.bin.harsh$dht$clusters$N[1:7, ]
PTM.ind.results.bin <- PTM.df.hr.size.bin.harsh$dht$individuals$N[1:7, ]

# bind results
PTM.results.bin <- rbind(PTM.grp.results.bin, PTM.ind.results.bin)
PTM.results.bin <- PTM.results.bin %>% rename(Year = Label) %>% 
  mutate(Label = rep(c("Grp", "Ind"), each=7)) %>% 
  mutate(Species = rep("PTM", times=14)) %>% 
  mutate(DetFun = rep("pooled", times=14)) %>% 
  mutate(Key = rep("Hr", times=14)) %>% 
  mutate(Adjust = rep("NA", times=14)) %>% 
  mutate(Covar = rep("size", times=14)) %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)



PTM_final_plot_binned <- ggplot(PTM.results.bin[PTM.results.bin$Label=="Grp",], aes(x=Year, y=Estimate))+
  geom_point(shape=16, size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
  theme_bw()+
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank())+
  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(PTM.results.bin, "Output/Results/PTM_results_final_binned.csv")
ggsave("Output/Results/Plots/Point_estimates/PTM_final_plot_binned.png", PTM_final_plot_binned, 
       dpi=300, width = 20, height = 20, units = "cm")


# compare binned and original results
PTM.results <- read.csv("Output/Results/PTM_results_final.csv")
PTM.results <- PTM.results[,-1]
PTM.results$analysis <- "original"

PTM.results.bin$analysis <- "binned"

PTM.results.all <- rbind(PTM.results,PTM.results.bin)

PTM_analysis_comparison_plot <- ggplot(PTM.results.all[PTM.results.all$Label=="Grp",], 
                                       aes(x=Year, y=Estimate, colour=analysis))+
                    geom_point(shape=16, size=2, position = position_dodge(width=0.3))+
                    geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width=0.3))+
                    theme_bw()+
                    theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),panel.border = element_blank())+
                    theme(axis.line = element_line(color = 'black'))

ggsave("Output/Results/Plots/Point_estimates/PTM_analysis_comparison_plot.png", 
       PTM_analysis_comparison_plot, dpi=300, width = 20, height = 20, units = "cm")


#### Stump-tailed macaque #####################################################
  ## Subset data ####

# subset STM data
STM.data <- allData[allData$species=="STM",] 
STM.data <- as.data.frame(STM.data)
head(STM.data)

# Total number of groups from all years
length(STM.data$distance) 

# for STM, I will pool all years for DF, but test year as a continous covariate in the DF model

# check 2020 data is there (note there were no obs in 2018)
unique(STM.data$year)
str(STM.data)

  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

# distance histograms
hist(STM.data$distance, main="All years", xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(STM.data$distance, main="All years", xlab="Distance (m)", breaks=c(40))

# Henous data. I think this species is similar to PTM in terms of behaviour and ecology. Mostly ground dwelling, especially when fleeing. And they also prefer EG SEG forest, hence why the furthest observation is 50m. The data are very sparse, but I think a HR model with a decent shoulder and then rapid drop off will be appropriate here. If there were more observations, I think it would look very similar to PTM


# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.STM <- 35

# Count the number of observations discarded
nrow(STM.data[STM.data$distance>trunc.STM,]) 

length(STM.data$distance)

1/34*100 # 3%

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(1,2))
plot(STM.data$size, STM.data$distance, main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.STM <- lm(distance~size, data=STM.data)
lines(STM.data$size, as.vector(predict(lm.STM, STM.data)))
summary(lm.STM)
# some evidence of size bias although model is not strong - mostly caused by a couple of outlying points I think

# Plot of Observer factor against distance
plot(STM.data$obs.observer,STM.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# quite a lot of variation, but this is perhaps expected with such a small dataset

# I can't include habitat or AMPM in the DF because those variables weren't recorded in earlier years


  ## Fit a detection function ####
    # Uniform ####

# Uni cosine

STM.df.uni.cos <- ds(data=STM.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.STM, key="uni", adjustment= "cos")
summary(STM.df.uni.cos)
# cosine(1)  

# Uni poly
STM.df.uni.poly <- ds(STM.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.STM, key="uni", adjustment= "poly")
summary(STM.df.uni.poly)
# poly(2)

## compare uni models
stm.uni.comp <- summarize_ds_models(STM.df.uni.cos, STM.df.uni.poly,
                                    output="plain")
stm.uni.comp[ ,1:4]
stm.uni.comp[ ,5:7]
# no difference

# plot the fits
par(mfrow=c(2,2))

# cos
plot(STM.df.uni.cos, main = "STM.df.uni.cos")

covar.fit <- ddf.gof(STM.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(STM.df.uni.poly, main = "STM.df.uni.poly")

covar.fit <- ddf.gof(STM.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# cos drops off quicker - p=0.5 ~20m, whereas poly is closer to 30m. I think poly looks most realistics in terms of the detection process, although it does have p=0 at 35m which we know isn't true. Still, I'll take poly forward

# cos selected

    # Half normal ####

# HN cosine
STM.df.hn.cos <- ds(STM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.STM, key="hn", adjustment= "cos")
summary(STM.df.hn.cos)
# no adjustment selected

# HN hermite
STM.df.hn.herm <- ds(STM.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.STM, key="hn", adjustment= "herm")
summary(STM.df.hn.herm)
# no adjustment selected


# plot the fit
par(mfrow=c(1,2))

# cos
plot(STM.df.hn.cos, main = "STM.df.hn.cos")

covar.fit <- ddf.gof(STM.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM
# looks similar to UNi cos

# cos is selected

    # Hazard rate ####

# HR cosine
STM.df.hr.cos <- ds(STM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.STM, key="hr", adjustment= "cos")
summary(STM.df.hr.cos)
# no adjustment selected


# HR poly
STM.df.hr.poly <- ds(STM.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.STM, key="hr", adjustment= "poly")
summary(STM.df.hr.poly)
# no adjustment selected 


# plot the fit
par(mfrow=c(1,2))

# hr cos
plot(STM.df.hr.cos, main = "STM.df.hr.cos")

covar.fit <- ddf.gof(STM.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# This is more what I was thinking. 

    # Compare primary models ####

stm.df.prim.comp <- summarize_ds_models(STM.df.uni.poly, STM.df.hn.cos,STM.df.hr.cos, 
                                           output = "plain")
stm.df.prim.comp[ ,1:5]
stm.df.prim.comp[ ,6:7]

# HR has the least support (dAIC>2). Uni and Hn similar

# plot all fits together
par(mfrow=c(3,2))

# uni cos
plot(STM.df.uni.poly, main = "STM.df.uni.poly")

covar.fit <- ddf.gof(STM.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HN cos
plot(STM.df.hn.cos, main = "STM.df.hn.cos")

covar.fit <- ddf.gof(STM.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HR cos
plot(STM.df.hr.cos, main = "STM.df.hr.cos")

covar.fit <- ddf.gof(STM.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

STM.df.hn.cos$dht$clusters$N$cv
STM.df.uni.poly$dht$clusters$N$cv
STM.df.hr.cos$dht$clusters$N$cv

# I genuinely think HR is the most biologically realistic shape. There isn't much difference in inference between the models, and the HR model has higher precision than HN (similar precision to Uni).

# HR model selected

    # Models with harsh truncation ####

# because of the sparsity of data for this species, and the shape of the data, I do not think truncating further here is sensible or useful.



    # Models with covariates ####

## cluster size
STM.df.hr.size <- ds(data=STM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.STM, key="hr", formula = ~size)
summary(STM.df.hr.size)

# plot
par(mfrow=c(1,2))
plot(STM.df.hr.size, main = "STM.df.hr.size")

covar.fit <- ddf.gof(STM.df.hr.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



## Observer
STM.df.hr.obs <- ds(data=STM.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.STM, key="hr", formula = ~obs.observer)
summary(STM.df.hr.obs)

# plot
plot(STM.df.hr.obs, main = "STM.df.hr.obs")

covar.fit <- ddf.gof(STM.df.hr.obs$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# woah

## stratum (year)
STM.df.hr.strat <- ds(data=STM.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.STM, key="hr", formula = ~stratum)
summary(STM.df.hr.strat)

# plot
plot(STM.df.hr.strat, main = "STM.df.hr.strat")

covar.fit <- ddf.gof(STM.df.hr.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)
# woah


    # Compare covariate models ####

stm.df.cov.comp <- summarize_ds_models(STM.df.hr.size, STM.df.hr.obs,STM.df.hr.strat,output = "plain")

stm.df.cov.comp[ ,1:5]
stm.df.cov.comp[ ,6:7]

# Model with stratum has most support. The coefficient for stratum is large
# The other two models have dAIC >4. Because of the sparsity of data, I worry about trying to include a covariate... 

# stratum model selected


    # Compare primary and covariate models ####

stm.df.final.compare <- summarize_ds_models(STM.df.hr.cos,STM.df.hr.strat, output = "plain")
stm.df.final.compare[ ,1:3]
stm.df.final.compare[ ,4:7]

# covariate model is the better model by AIC, but is producing absolutley ludicrous estimates. There are too few data points to incldude a covariate.

# original HR model is selected

  ## STM Final Results ####

summary(STM.df.hr.cos)

# extract estimates
STM.grp.abund   <- STM.df.hr.cos$dht$clusters$N[1:7,1:6]
STM.ind.abund   <- STM.df.hr.cos$dht$individuals$N[1:7,1:6]
STM.grp.density <- STM.df.hr.cos$dht$clusters$D[1:7,1:6]
STM.ind.density <- STM.df.hr.cos$dht$individuals$D[1:7,1:6]

STM.grp.abund <- STM.grp.abund %>% rename(Year=Label) %>% 
                  mutate(Label="Grp")
STM.ind.abund <- STM.ind.abund %>% rename(Year=Label) %>% 
                  mutate(Label="Ind")
STM.results <- rbind(STM.grp.abund,STM.ind.abund)
STM.results <- STM.results %>% rename(N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl)

STM.grp.density <- STM.grp.density %>% rename(Year=Label) %>% 
                  mutate(Label="Grp")
STM.ind.density <- STM.ind.density %>% rename(Year=Label) %>% 
                  mutate(Label="Ind")
STM.density <- rbind(STM.grp.density,STM.ind.density)
STM.density <- STM.density %>% rename(D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl)
STM.density <- STM.density[,-c(1,7)]

# bind results
STM.results <- cbind(STM.results, STM.density)
STM.results <- STM.results %>% 
                mutate(Species = rep("STM", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hr", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("NA", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,N,n_se,n_cv,n_lcl,n_ucl,
                       D,d_se,d_cv,d_lcl,d_ucl)

STM.results <- STM.results %>% mutate(D = D*1000000, d_se=d_se*1000000, d_lcl=d_lcl*1000000,
                                      d_ucl=d_ucl*1000000)


STM_final_plot <- ggplot(STM.results[STM.results$Label=="Grp",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(STM.results, "Output/Results/STM_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/STM_final_plot.png", STM_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")


### extract density
STM.grp.density <- STM.df.hr.cos$dht$clusters$D[1:7, ]
STM.grp.density <- STM.grp.density %>% rename(Year=Label) %>% 
                    mutate(Species="STM") %>% 
                    select(Species,Year,Estimate,se,cv,lcl,ucl,df)


#### Banteng ##################################################################
  ## Subset data ####

### NO BANTENG OBSERVATIONS IN 2020

# Previously, I tested the inclusion of WWF banteng data from 2014 and 2016, to see what difference it made. I compared it to WCS-only data.  The inclusion of the extra data did not improve estimates (in fact made them worse). Therefore I will stick to just WCS data.  I have left the comparison code in below for reference, and so that it can be repeated in the future if required. 

## NOTE: Becuase I have not needed to re-run the analysis for 2020, I have not gone through and changed the code. Therefore for anyone re-doing this analysis, you will need to edit the Banteng code as you go. This mostly relates to the names of the input data. E.g. in the code below the data is called things like "region.table.nostrata" whereas it will need to be changed to "full.region.table". Check what the correct names are for the input data in the main "Load libraries & data" section at the top of the script. I have changed the code for the final BTG model (Uni.poly, WCS only) so that I could re-format the results to fit with the other species.


## Data from WCS only
# subset BTG data
BTG.data <- allData[allData$species=="BAN",] 
BTG.data <- as.data.frame(BTG.data)
head(BTG.data)

## Data from WWF 2014 and 2016
wwf_data <- read.csv("Input/BTG_extra_data.csv")

# merge it with our data
BTG.data.extra <- rbind(BTG.data, wwf_data)

# create new sample.table. The Sample.Labels and effort are meaningless here. The region label is just WWF, as that will keep the estiamtes separate from ours (the estimates for the WWF region will be meaningless)
sample.table.wwf <- data.frame(Sample.Label = c("wwf_14","wwf_16"),
                               Region.Label = c("wwf","wwf"),
                               Effort = c(40000,40000))

sample.table.extra <- rbind(sample.table.nostrata, sample.table.wwf) 

# create new obs.table
obs.table.wwf <- read.csv("Input/obs.table.wwf.csv")
obs.table.extra <- rbind(obs.table.nostrata, obs.table.wwf)

# create new region.table. Area here is meaningless, so I have just copied the WCS Area
region.table.wwf <- data.frame(Region.Label = "wwf",
                               Area = 1807000000) 

region.table.extra <- rbind(region.table.nostrata, region.table.wwf)


# Total number of groups from all years - WCS only
length(BTG.data$distance)  # 20

# total number of groups including WWF data
length(BTG.data.extra$distance) # 111

# I will run through the analysis for both datasets at once


  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

## WCS only data
# distance histograms
hist(BTG.data$distance[BTG.data$distance<80], main="All years", xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(BTG.data$distance[BTG.data$distance<80], main="All years", xlab="Distance (m)", breaks=c(40))

# pretty shitty data!  ONly 20 observations though, so to be expected

## WCS and WWF data
# distance histograms
hist(BTG.data.extra$distance[BTG.data.extra$distance<100], main="All years", xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(BTG.data.extra$distance[BTG.data.extra$distance<100], 
     main="All years", xlab="Distance (m)", breaks=c(40))

# still not great, but at least there's some shape to it!


# Save the chosen truncation distances for later use. Try harsher truncation to shrink CIs later
trunc.BTG <- 70
trunc.BTG.extra <- 75 

# Count the number of observations discarded WCS data only
nrow(BTG.data[BTG.data$distance>trunc.BTG,]) 
length(BTG.data$distance)

1/20*100 # 5%

# Count the number of observations discarded WCS & WWF data 
nrow(BTG.data.extra[BTG.data.extra$distance>trunc.BTG.extra,]) 
length(BTG.data.extra$distance)

11/111*100 # 10%


## Plots of covars against distance

# Plot of distance against cluster size WCS only
par(mfrow=c(1,2))

plot(BTG.data$size, BTG.data$distance, main="WCS", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.BTG <- lm(distance~size, data=BTG.data)
lines(BTG.data$size, as.vector(predict(lm.BTG, BTG.data)))
summary(lm.BTG)
# some evidence of size bias although model is not strong 

# Plot of distance against cluster size WCS & WWF
plot(BTG.data.extra$size, BTG.data.extra$distance, main="WCS & WWF", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.BTG.extra <- lm(distance~size, data=BTG.data.extra)
lines(BTG.data.extra$size, as.vector(predict(lm.BTG.extra, BTG.data.extra)))
summary(lm.BTG)

# The two models and plots show totally different slope directions.  The WCS only model is being influenced by one large group. the combined data show a negative slope.  Basically, not sure we can say much here...except that we should probably take the WCS-only data with a large pinch of salt

# Can't plot observer as a factor for the WWF data as I don't have that information. And even for the WCS data I am not sure there is any point for 20 observations


# I can't include habitat or AMPM in the DF because those variables weren't recorded in earlier years in WCS data, and I don't have this information in the WWF data


  ## Fit a detection function ####
    # Uniform ####


### WCS data only

# Uni cosine

BTG.df.uni.cos <- ds(data=BTG.data, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=trunc.BTG, key="uni", adjustment= "cos")
summary(BTG.df.uni.cos)
# cosine(1)  

# Uni poly
BTG.df.uni.poly <- ds(data=BTG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.BTG, key="uni", adjustment= "poly")
summary(BTG.df.uni.poly)
# poly(2)

## compare uni models
btg.uni.comp <- summarize_ds_models(BTG.df.uni.cos, BTG.df.uni.poly,
                                    output="plain")
btg.uni.comp[ ,1:4]
btg.uni.comp[ ,5:7]
# no difference

# plot the fits
par(mfrow=c(2,2))

# cos
plot(BTG.df.uni.cos, main = "BTG.df.uni.cos")
covar.fit <- ddf.gof(BTG.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(BTG.df.uni.poly, main = "BTG.df.uni.poly")
covar.fit <- ddf.gof(BTG.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# cos suggests a levelling off of p on 0.6 at 70m - not realistic. Poly has a much more realistic tail end 

# poly selected



### WCS and WWF data

# Uni cosine

BTG.df.uni.cos.extra <- ds(data=BTG.data.extra, region.table=region.table.extra, 
                      sample.table=sample.table.extra, obs.table=obs.table.extra, 
                      truncation=trunc.BTG.extra, key="uni", adjustment= "cos")
summary(BTG.df.uni.cos.extra)
# cosine(1)  

# Uni poly
BTG.df.uni.poly.extra <- ds(data=BTG.data.extra, region.table=region.table.extra, 
                      sample.table=sample.table.extra, obs.table=obs.table.extra, 
                      truncation=trunc.BTG.extra, key="uni", adjustment= "poly")
summary(BTG.df.uni.poly.extra)
# poly(2)

## compare uni models
btg.uni.comp.ex <- summarize_ds_models(BTG.df.uni.cos.extra, BTG.df.uni.poly.extra,
                                    output="plain")
btg.uni.comp.ex[ ,1:4]
btg.uni.comp.ex[ ,5:7]
# cos better, although poly dAIC < 2

# plot the fits
par(mfrow=c(2,2))

# cos
plot(BTG.df.uni.cos.extra, main = "BTG.df.uni.cos.extra")
covar.fit <- ddf.gof(BTG.df.uni.cos.extra$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(BTG.df.uni.poly.extra, main = "BTG.df.uni.poly.extra")
covar.fit <- ddf.gof(BTG.df.uni.poly.extra$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# main difference is that poly is a much gentler decrease. I think this is the more reasonable model - Banteng are large animals that live in open forest, therefore are relatively easy to see. Cos suggests p = 0.5 around 40m, which I think is too close

# poly selected

    # Half normal ####

### WCS only

# HN cosine
BTG.df.hn.cos <- ds(data=BTG.data, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=trunc.BTG, key="hn", adjustment= "cos")
summary(BTG.df.hn.cos)
# no adjustment selected

# HN hermite
BTG.df.hn.herm <- ds(data=BTG.data, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=trunc.BTG, key="hn", adjustment= "herm")
summary(BTG.df.hn.herm)
# no adjustment selected


# plot the fit
par(mfrow=c(1,2))

# cos
plot(BTG.df.hn.cos, main = "BTG.df.hn.cos")
covar.fit <- ddf.gof(BTG.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# cos is selected


### WCS and WWF

# HN cosine
BTG.df.hn.cos.extra <- ds(data=BTG.data.extra, region.table=region.table.extra, 
                      sample.table=sample.table.extra, obs.table=obs.table.extra, 
                      truncation=trunc.BTG.extra, key="hn", adjustment= "cos")
summary(BTG.df.hn.cos.extra)
# no adjustment selected

# HN hermite
BTG.df.hn.herm.extra <- ds(data=BTG.data.extra, region.table=region.table.extra, 
                      sample.table=sample.table.extra, obs.table=obs.table.extra, 
                      truncation=trunc.BTG.extra, key="hn", adjustment= "herm")
summary(BTG.df.hn.herm.extra)
# no adjustment selected


# plot the fit
par(mfrow=c(1,2))

# cos
plot(BTG.df.hn.cos.extra, main = "BTG.df.hn.cos.extra")
covar.fit <- ddf.gof(BTG.df.hn.cos.extra$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# probably too steep a drop off I reckon

    # Hazard rate ####

### WCS only

# HR cosine
BTG.df.hr.cos <- ds(data=BTG.data, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=trunc.BTG, key="hr", adjustment= "cos")
summary(BTG.df.hr.cos)
# no adjustment selected


# HR poly
BTG.df.hr.poly <- ds(data=BTG.data, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=trunc.BTG, key="hr", adjustment= "poly")
summary(BTG.df.hr.poly)
# no adjustment selected 


# plot the fit
par(mfrow=c(1,2))

# hr cos
plot(BTG.df.hr.cos, main = "BTG.df.hr.cos")
covar.fit <- ddf.gof(BTG.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# As I expected, massive shoulder - p=1 up to 45m.  Perhaps not entirley implausible?


### WCS and WWF

# HR cosine
BTG.df.hr.cos.extra <- ds(data=BTG.data.extra, region.table=region.table.extra, 
                      sample.table=sample.table.extra, obs.table=obs.table.extra, 
                      truncation=trunc.BTG.extra, key="hr", adjustment= "cos")
summary(BTG.df.hr.cos.extra)
# no adjustment selected


# HR poly
BTG.df.hr.poly.extra <- ds(data=BTG.data.extra, region.table=region.table.extra, 
                      sample.table=sample.table.extra, obs.table=obs.table.extra, 
                      truncation=trunc.BTG.extra, key="hr", adjustment= "poly")
summary(BTG.df.hr.poly.extra)
# no adjustment selected 


# plot the fit
par(mfrow=c(1,2))

# hr cos
plot(BTG.df.hr.cos.extra, main = "BTG.df.hr.cos.extra")
covar.fit <- ddf.gof(BTG.df.hr.cos.extra$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)


    # Compare primary models ####

### WCS only
stm.df.prim.comp <- summarize_ds_models(BTG.df.uni.poly, BTG.df.hn.cos,BTG.df.hr.cos, 
                                           output = "plain")
stm.df.prim.comp[ ,1:5]
stm.df.prim.comp[ ,6:7]

# ALl models have some support.   

# plot all fits together
par(mfrow=c(3,2))

# uni poly
plot(BTG.df.uni.poly, main = "BTG.df.uni.poly")

covar.fit <- ddf.gof(BTG.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HN cos
plot(BTG.df.hn.cos, main = "BTG.df.hn.cos")

covar.fit <- ddf.gof(BTG.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HR cos
plot(BTG.df.hr.cos, main = "BTG.df.hr.cos")

covar.fit <- ddf.gof(BTG.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# I think the uniform model has the most realistic fit here. I don't like how the HN fit flattens out at the tail, and the HR p=1 shulder is too large for me

# Uni model selected


### WCS and WWF

btg.df.prim.comp.ex <- summarize_ds_models(BTG.df.uni.poly.extra, BTG.df.hn.cos.extra,
                                           BTG.df.hr.cos.extra, 
                                           output = "plain")
btg.df.prim.comp.ex[ ,1:5]
btg.df.prim.comp.ex[ ,6:7]

# HR model dAIC =2.4  

# plot all fits together
par(mfrow=c(3,2))

# uni poly
plot(BTG.df.uni.poly.extra, main = "BTG.df.uni.poly.extra")

covar.fit <- ddf.gof(BTG.df.uni.poly.extra$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HN cos
plot(BTG.df.hn.cos.extra, main = "BTG.df.hn.cos.extra")

covar.fit <- ddf.gof(BTG.df.hn.cos.extra$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HR cos
plot(BTG.df.hr.cos.extra, main = "BTG.df.hr.cos.extra")

covar.fit <- ddf.gof(BTG.df.hr.cos.extra$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# I think HR drops off too quickly here.  It has p=0.5 at 30m which is too close. Uniform is still the best I think - slowest decrease in p

# uniform selected


    # Models with harsh truncation ####

# I will not try harsh truncation with the WCS data only, as there are only 20 obs. But I will try with the WWF data


# set trunc distance 
btg.trunc.harsh <- 50

par(mfrow=c(1,2))

# Uni poly harsh
BTG.df.unif.poly.harsh <- ds(data=BTG.data.extra, region.table=region.table.extra, 
                    sample.table=sample.table.extra, obs.table=obs.table.extra,
                    truncation=btg.trunc.harsh, key="uni", adjustment="poly")
summary(BTG.df.unif.poly.harsh)
# poly(2)

# compare original with harsh truncation
ddf.gof(BTG.df.unif.poly.harsh$ddf, main = "BTG.df.unif.poly.harsh", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(BTG.df.uni.poly.extra$ddf, main = "BTG.df.uni.poly.extra", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# hard to say 

summary(BTG.df.unif.poly.harsh) # CVs= 0.94, 0.46, 0.71, 0.41, 0.77, 1
summary(BTG.df.uni.poly.extra) # CVs = 0.94, 0.48, 0.71, 0.40, 0.77, 1
# NO real difference. I will stick with the original


    # Models with covariates ####

### I will only test covariates on the combined data, as there are not enough obs in the WCS-only data

## cluster size
BTG.df.hn.size <- ds(data=BTG.data.extra, region.table=region.table.extra, 
                    sample.table=sample.table.extra, obs.table=obs.table.extra,
                    truncation=trunc.BTG.extra, key="hn", formula = ~size)
summary(BTG.df.hn.size)

# plot
par(mfrow=c(1,2))
plot(BTG.df.hn.size, main = "BTG.df.hn.size")
covar.fit <- ddf.gof(BTG.df.hn.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)



## stratum (year). I know what year the WWF data are from so can test this
BTG.df.hn.strat <- ds(data=BTG.data.extra, region.table=region.table.extra, 
                    sample.table=sample.table.extra, obs.table=obs.table.extra,
                    truncation=trunc.BTG.extra, key="hn", formula = ~stratum)
summary(BTG.df.hn.strat)

# plot
plot(BTG.df.hn.strat, main = "BTG.df.hn.strat")
covar.fit <- ddf.gof(BTG.df.hn.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)


## stratum + size
BTG.df.hn.strat.size <- ds(data=BTG.data, region.table=region.table.extra, 
                    sample.table=sample.table.extra, obs.table=obs.table.extra,
                    truncation=trunc.BTG.extra, key="hn", formula = ~stratum)
summary(BTG.df.hn.strat)

# plot
plot(BTG.df.hn.strat, main = "BTG.df.hn.strat")
covar.fit <- ddf.gof(BTG.df.hn.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)


    # Compare covariate models ####

btg.df.cov.comp <- summarize_ds_models(BTG.df.hn.size,BTG.df.hn.strat,output = "plain")

btg.df.cov.comp[ ,1:5]
btg.df.cov.comp[ ,6:7]

# Both models have some support
summary(BTG.df.hn.size)
summary(BTG.df.hn.strat)

# the estimates for 2011 and 2013 are just nonsense. There is no way in hell there were over 2000 Banteng wandering around Seima in those years.



    # Compare primary and covariate models ####

# this is only a comparison for the combined data
btg.df.final.compare <- summarize_ds_models(BTG.df.uni.poly.extra,BTG.df.hn.strat, output = "plain")
btg.df.final.compare[ ,1:3]
btg.df.final.compare[ ,4:7]

# Original model has lower AIC
summary(BTG.df.uni.poly.extra)

# original model is selected

  ## BTG Final Results ####

    # Extract WCS results only ####

summary(BTG.df.uni.poly)

# extract estimates
BTG.grp.abund <- BTG.df.uni.poly$dht$clusters$N[1:7, ]
BTG.grp.abund <- BTG.grp.abund %>% rename(Year = Label) %>% mutate(Label="Grp")
BTG.ind.abund <- BTG.df.uni.poly$dht$individuals$N[1:7, ]
BTG.ind.abund <- BTG.ind.abund %>% rename(Year=Label) %>% mutate(Label="Ind")
BTG.abund <- rbind(BTG.grp.abund,BTG.ind.abund)
BTG.abund <- BTG.abund %>% rename(N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl)
BTG.abund <- BTG.abund[,-7]

### extract density
BTG.grp.density <- BTG.df.uni.poly$dht$clusters$D[1:7, ]
BTG.grp.density <- BTG.grp.density %>% rename(Year=Label) %>% mutate(Label="Grp")
BTG.ind.density <- BTG.df.uni.poly$dht$individuals$D[1:7, ]
BTG.ind.density <- BTG.ind.density %>% rename(Year=Label) %>% mutate(Label="Ind")
BTG.density <- rbind(BTG.grp.density,BTG.ind.density)
BTG.density <- BTG.density %>% rename(D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl)
BTG.density <- BTG.density[,-7]
BTG.density <- BTG.density[,-c(1,7)]


# bind results
BTG.results <- cbind(BTG.abund, BTG.density)
BTG.results <- BTG.results %>%  
                mutate(Species = rep("BTG", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Uni", times=14)) %>% 
                mutate(Adjust = rep("SimPoly", times=14)) %>% 
                mutate(Covar = rep("NA", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,
                       N,n_se,n_cv,n_lcl,n_ucl,
                       D,d_se,d_cv,d_lcl,d_ucl)

BTG.results <- BTG.results %>% mutate(D=D*1000000,d_se=d_se*1000000,d_lcl=d_lcl*1000000,d_ucl=d_ucl*1000000)





    # Extract combined results ####

# extract estimates
BTG.grp.results.comb <- BTG.df.uni.poly.extra$dht$clusters$N[1:6, ]
BTG.ind.results.comb <- BTG.df.uni.poly.extra$dht$individuals$N[1:6, ]

# bind results
BTG.results.comb <- rbind(BTG.grp.results.comb, BTG.ind.results.comb)
BTG.results.comb <- BTG.results.comb %>% rename(Year = Label) %>% 
                mutate(Label = rep(c("Grp_combined", "Ind_combined"), each=6)) %>% 
                mutate(Species = rep("BTG", times=12)) %>% 
                mutate(DetFun = rep("pooled", times=12)) %>% 
                mutate(Key = rep("Uni", times=12)) %>% 
                mutate(Adjust = rep("SimPoly", times=12)) %>% 
                mutate(Covar = rep("NA", times=12)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)


  # bind results & plot ####



BTG.results.all <- rbind(BTG.results, BTG.results.comb)

BTG_comparison_plot <- ggplot(BTG.results.all, aes(x=year, y= ind.abund, group=data, colour=data))+
                  geom_point(shape=16, size=3, position = position_dodge(width = 0.5))+
                  geom_errorbar(aes(ymin=ind.abund.lcl, ymax=ind.abund.ucl), width=0.4,
                                position = position_dodge(width = 0.5))+
                  ylim(0,6600)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

BTG_wcs_plot <- ggplot(BTG.results[BTG.results$Label=="Ind",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=3)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.4)+
                  ylim(0,2000)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(BTG.results, "Output/Results/BTG_results_final.csv")
ggsave("Output/Results/Plots/BTG_comparison_plot.png", BTG_comparison_plot, 
       dpi=300, width = 20, height = 20, units = "cm")
ggsave("Output/Results/Plots/BTG_wcs_plot.png", BTG_wcs_plot, 
       dpi=300, width = 20, height = 20, units = "cm")



#### Gaur ####################################################################
  ## Subset data ####

# subset GAU data
GAU.data <- allData[allData$species=="GAU",] 
GAU.data <- as.data.frame(GAU.data)
head(GAU.data)

# Total number of groups from all years
length(GAU.data$distance) 

# for GAU, I will pool all years for DF, and will NOT test year as a continous covariate in the DF model as there are too few observations

# check 2020 data is there
unique(GAU.data$year)


  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

# distance histograms
hist(GAU.data$distance[GAU.data$distance<100], main="All years", xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(GAU.data$distance[GAU.data$distance<100], main="All years", xlab="Distance (m)", breaks=c(40))

# Although these data are sparse, they make sense. Gaur are fucking massive, and you don't miss them when they leg it through the EG and SEG forest. There is no mistaking what they are. Therefore p is probably super high up until around 40m, after which you wouldn't be able to confidently ID them (at that distance the noise could be a group of pigs for example).


# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.GAU <- 40

# Count the number of observations discarded
nrow(GAU.data[GAU.data$distance>trunc.GAU,]) 

length(GAU.data$distance)

1/35*100 # 3%

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(1,2))
plot(GAU.data$size, GAU.data$distance, main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.GAU <- lm(distance~size, data=GAU.data)
lines(GAU.data$size, as.vector(predict(lm.GAU, GAU.data)))
summary(lm.GAU)
# no real evidence of size bias - looks to be an opposite effect if anything - larger groups closer to the line.

# Plot of Observer factor against distance
plot(GAU.data$obs.observer,GAU.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# quite a lot of variation, but this is perhaps expected with such a small dataset

# I can't include habitat or AMPM in the DF because those variables weren't recorded in earlier years


  ## Fit a detection function ####
    # Uniform ####

# Uni cosine

GAU.df.uni.cos <- ds(data=GAU.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GAU, key="uni", adjustment= "cos")
summary(GAU.df.uni.cos)
# cosine(1)  

# Uni poly
GAU.df.uni.poly <- ds(data=GAU.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GAU, key="uni", adjustment= "poly")
summary(GAU.df.uni.poly)
# poly(2,4)

## compare uni models
gau.uni.comp <- summarize_ds_models(GAU.df.uni.cos, GAU.df.uni.poly,
                                    output="plain")
gau.uni.comp[ ,1:4]
gau.uni.comp[ ,5:7]
# no big difference. Cos the preferred model but poly dAIC < 1.5

# plot the fits
par(mfrow=c(2,2))

# cos
plot(GAU.df.uni.cos, main = "GAU.df.uni.cos")

covar.fit <- ddf.gof(GAU.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(GAU.df.uni.poly, main = "GAU.df.uni.poly")

covar.fit <- ddf.gof(GAU.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# very similar - will go with AIC. I don't think there is enough of a shoulder here

# cos selected

    # Half normal ####

# HN cosine
GAU.df.hn.cos <- ds(data=GAU.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GAU, key="hn", adjustment= "cos")
summary(GAU.df.hn.cos)
# no adjustment selected

# HN hermite
GAU.df.hn.herm <- ds(data=GAU.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.GAU, key="hn", adjustment= "herm")
summary(STM.df.hn.herm)
# no adjustment selected


# plot the fit
par(mfrow=c(1,2))

# cos
plot(GAU.df.hn.cos, main = "GAU.df.hn.cos")

covar.fit <- ddf.gof(GAU.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# cos is selected. Looks very similar to UNi models

    # Hazard rate ####

# HR cosine
GAU.df.hr.cos <- ds(data=GAU.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GAU, key="hr", adjustment= "cos")
summary(GAU.df.hr.cos)
# no adjustment selected


# HR poly
GAU.df.hr.poly <- ds(data=GAU.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.GAU, key="hr", adjustment= "poly")
summary(GAU.df.hr.poly)
# no adjustment selected 


# plot the fit
par(mfrow=c(1,2))

# hr cos
plot(GAU.df.hr.cos, main = "GAU.df.hr.cos")

covar.fit <- ddf.gof(GAU.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# This is more like it

    # Compare primary models ####

gau.df.prim.comp <- summarize_ds_models(GAU.df.uni.cos, GAU.df.hn.cos,GAU.df.hr.cos, 
                                           output = "plain")
gau.df.prim.comp[ ,1:5]
gau.df.prim.comp[ ,6:7]

# All models have some support. Initially I had selected the HR model, but this model caused all sorts of issues in both the DSM and the trend analysis. Therefore I have reverted to the more stable HN

    # Models with harsh truncation ####

# because of the sparsity of data for this species, and the shape of the data, I do not think truncating further here is sensible or useful.



    # Models with covariates ####

# As with STM, I don't think we have enough data here to include covariates, but I'll have a look anyway

## cluster size
GAU.df.hr.size <- ds(data=GAU.data, region.table=region.table.nostrata, 
                    sample.table=sample.table.nostrata, obs.table=obs.table.nostrata,
                    truncation=trunc.GAU, key="hr", formula = ~size)
summary(GAU.df.hr.size)

# plot
par(mfrow=c(1,2))
plot(GAU.df.hr.size, main = "GAU.df.hr.size")

covar.fit <- ddf.gof(GAU.df.hr.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)



## Observer
GAU.df.hr.obs <- ds(data=GAU.data, region.table=region.table.nostrata, 
                    sample.table=sample.table.nostrata, obs.table=obs.table.nostrata,
                    truncation=trunc.GAU, key="hr", formula = ~obs.observer)
summary(GAU.df.hr.obs)
# model failed to converge


## stratum (year)
GAU.df.hr.strat <- ds(data=GAU.data, region.table=region.table.nostrata, 
                    sample.table=sample.table.nostrata, obs.table=obs.table.nostrata,
                    truncation=trunc.GAU, key="hr", formula = ~stratum)
summary(GAU.df.hr.strat)
# model failed to converge


    # Compare covariate models ####

# only 1 model
gau.df.cov.comp <- summarize_ds_models(STM.df.hr.size,output = "plain")

gau.df.cov.comp[ ,1:5]
gau.df.cov.comp[ ,6:7]

# CvM p value is 0.001 and so the model is rejected

    # Compare primary and covariate models ####

# no appropriate covariate models, and so the original model is selected

  ## GAU Final Results ####

summary(GAU.df.hn.cos)

# extract estimates
GAU.grp.abund <- GAU.df.hn.cos$dht$clusters$N[1:7, ]
GAU.grp.abund <- GAU.grp.abund %>% rename(Year = Label) %>% mutate(Label="Grp")
GAU.ind.abund <- GAU.df.hn.cos$dht$individuals$N[1:7, ]
GAU.ind.abund <- GAU.ind.abund %>% rename(Year=Label) %>% mutate(Label="Ind")
GAU.abund <- rbind(GAU.grp.abund,GAU.ind.abund)
GAU.abund <- GAU.abund %>% rename(N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl)
GAU.abund <- GAU.abund[,-7]

### extract density
GAU.grp.density <- GAU.df.hn.cos$dht$clusters$D[1:7, ]
GAU.grp.density <- GAU.grp.density %>% rename(Year=Label) %>% mutate(Label="Grp")
GAU.ind.density <- GAU.df.hn.cos$dht$individuals$D[1:7, ]
GAU.ind.density <- GAU.ind.density %>% rename(Year=Label) %>% mutate(Label="Ind")
GAU.density <- rbind(GAU.grp.density,GAU.ind.density)
GAU.density <- GAU.density %>% rename(D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl)
GAU.density <- GAU.density[,-7]
GAU.density <- GAU.density[,-c(1,7)]


# bind results
GAU.results <- cbind(GAU.abund, GAU.density)
GAU.results <- GAU.results %>%  
                mutate(Species = rep("GAU", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hn", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("NA", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,
                       N,n_se,n_cv,n_lcl,n_ucl,
                       D,d_se,d_cv,d_lcl,d_ucl)

GAU.results <- GAU.results %>% mutate(D=D*1000000,d_se=d_se*1000000,d_lcl=d_lcl*1000000,d_ucl=d_ucl*1000000)




GAU_final_plot <- ggplot(GAU.results[GAU.results$Label=="Ind",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  ylim(0,2700)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(GAU.results, "Output/Results/GAU_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/GAU_final_plot.png", GAU_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")


#### Wild pig ################################################################
  ## Subset data ####

# subset PIG data
PIG.data <- allData[allData$species=="PIG",] 
PIG.data <- as.data.frame(PIG.data)
head(PIG.data)

# Total number of groups from all years
length(PIG.data$distance) 

# for PIG, I will pool all years for DF, and will test year as a continous covariate in the DF model

# check 2020 data is there
unique(PIG.data$year)


  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

# distance histograms
hist(PIG.data$distance[PIG.data$distance<100], main="All years", xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(PIG.data$distance[PIG.data$distance<100], main="All years", xlab="Distance (m)", breaks=c(40))

# Some evidence of evasive movement between 0-10m. This is believable for pigs. steep drop off beyond 20m.  This also makes sense - beyond a certain distance in denser habitat then you will struggle to confidently ID the speceis as they leg it very quickly.  The long tail in the data will be from more open habitats when visibility is better


### update - binning will be done for PIG as there is clumping around 0

# identify best bins
hist(PIG.data$distance[PIG.data$distance<60], breaks=c(0,5,10,15,20,25,30,35,40,50,60))
hist(PIG.data$distance[PIG.data$distance<60], breaks=c(0,5,12,16,25,32,45,60))


# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.PIG <- 50
trunc.PIG.bin <- 60

# Count the number of observations discarded
nrow(PIG.data[PIG.data$distance>trunc.PIG,]) 

length(PIG.data$distance)

18/2013*100 # 1% - could be harsher - will test below

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(1,2))
plot(PIG.data$size, PIG.data$distance, 
     main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.PIG <- lm(distance~size, data=PIG.data)
lines(PIG.data$size, as.vector(predict(lm.PIG, PIG.data)))
summary(lm.PIG)
# no evidence ofpositive size bias - looks to be an opposite effect if anything - larger groups closer to the line. Which is not what I was expecting. Pigs live in large groups, and I wouls assume larger groups would be easier to spot. Although the lm could be being dragged down by that one observation of >25 individulas at close to 0m.  

# check the lm results when that outlier removed
lm.PIG2 <- lm(distance~size, data=PIG.data[PIG.data$size<25,])
lines(PIG.data$size, as.vector(predict(lm.PIG2, PIG.data)))
# still the same trend, just flatter slope. 

# Plot of Observer factor against distance
plot(PIG.data$obs.observer,PIG.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# Not as much variation as some other species. Quite a lot of large outlier distances

# I can't include habitat or AMPM in the DF because those variables weren't recorded in earlier years


  ## Fit a detection function ####
    # Uniform ####


### unbinned

# Uni cosine
PIG.df.uni.cos <- ds(data=PIG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.PIG, key="uni", adjustment= "cos")
summary(PIG.df.uni.cos)
# cosine(1)  

# Uni poly
PIG.df.uni.poly <- ds(data=PIG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.PIG, key="uni", adjustment= "poly")
summary(PIG.df.uni.poly)
# poly(2,4,6)

## compare uni models
pig.uni.comp <- summarize_ds_models(PIG.df.uni.cos, PIG.df.uni.poly,
                                    output="plain")
pig.uni.comp[ ,1:4]
pig.uni.comp[ ,5:7]
# similar AICS. cos preferred 

# plot the fits
par(mfrow=c(2,2))

# cos
plot(PIG.df.uni.cos, main = "PIG.df.uni.cos")

covar.fit <- ddf.gof(PIG.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(PIG.df.uni.poly, main = "PIG.df.uni.poly")

covar.fit <- ddf.gof(PIG.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# cos drops off slightly slower, and the fit includes more of the hump aournd 30m. But they are very similar. the QQ plot for cos is nicer

# cos selected


### binned

# Uni cosine
PIG.df.uni.cos.bin <- ds(data=PIG.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                     truncation=trunc.PIG.bin, key="uni", adjustment= "cos",
                     cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.uni.cos.bin)
# cosine(1,2) 

# Uni poly
PIG.df.uni.poly.bin <- ds(data=PIG.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table, 
                         truncation=trunc.PIG.bin, key="uni", adjustment= "poly",
                         cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.uni.poly.bin)
# poly(2,4,6)

## compare uni models
pig.uni.comp.bin <- summarize_ds_models(PIG.df.uni.cos.bin, PIG.df.uni.poly.bin,
                                    output="plain")
pig.uni.comp.bin[ ,1:4]
pig.uni.comp.bin[ ,5:7]
# cos preferred

# plot both
par(mfrow=c(1,2))
plot(PIG.df.uni.cos.bin, main="cos")
plot(PIG.df.uni.poly.bin, main="poly")

# PIG.df.uni.cos.bin selected

    # Half normal ####


### unbinned


# HN cosine
PIG.df.hn.cos <- ds(data=PIG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.PIG, key="hn", adjustment= "cos")
summary(PIG.df.hn.cos)
# no adjustment selected

# HN hermite
PIG.df.hn.herm <- ds(data=PIG.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.PIG, key="hn", adjustment= "herm")
summary(PIG.df.hn.herm)
# no adjustment selected


# plot the fit
par(mfrow=c(1,2))

# cos
plot(PIG.df.hn.cos, main = "PIG.df.hn.cos")

covar.fit <- ddf.gof(PIG.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# cos is selected. Looks similar to UNi models 


### binned

# Hn cosine
PIG.df.hn.cos.bin <- ds(data=PIG.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table, 
                         truncation=trunc.PIG.bin, key="hn", adjustment= "cos",
                         cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.hn.cos.bin)
# cosine(2,3)

# Hn herm
PIG.df.hn.herm.bin <- ds(data=PIG.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.PIG.bin, key="hn", adjustment= "herm",
                        cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.hn.herm.bin)
# key only

## compare hn models
pig.hn.comp.bin <- summarize_ds_models(PIG.df.hn.cos.bin, PIG.df.hn.herm.bin,
                                        output="plain")
pig.hn.comp.bin[ ,1:4]
pig.hn.comp.bin[ ,5:7]
# cos preferred

# plot both
plot(PIG.df.hn.cos.bin, main="cos")
plot(PIG.df.hn.herm.bin, main="herm")
# cos definitely the better fit, but it is because it is fitting the tail end shoulder. Overfitted? Potentially not, because the tail end shoulder probably reflects the detection process in more open habitat, whereas the steep drop from the shoulder close to the line reflects the process in dense forest.

# PIG.df.hn.cos.bin selected


    # Hazard rate ####

## unbinned

# HR cosine
PIG.df.hr.cos <- ds(data=PIG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.PIG, key="hr", adjustment= "cos")
summary(PIG.df.hr.cos)
# no adjustment selected


# HR poly
PIG.df.hr.poly <- ds(data=PIG.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.PIG, key="hr", adjustment= "poly")
summary(PIG.df.hr.poly)
# no adjustment selected 


# plot the fit
par(mfrow=c(1,2))

# hr cos
plot(PIG.df.hr.cos, main = "PIG.df.hr.cos")

covar.fit <- ddf.gof(PIG.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# wide shoulder - p=1 up to 15m.  Probably fairly realistic. Decent QQ plot. This feels right to me, based on the observation process


### binned

# HR cosine
PIG.df.hr.cos.bin <- ds(data=PIG.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.PIG.bin, key="hr", adjustment= "cos",
                        cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.hr.cos.bin)
# key only

# HR poly
PIG.df.hr.poly.bin <- ds(data=PIG.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.PIG.bin, key="hr", adjustment= "poly",
                        cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.hr.poly.bin)
# key only

plot(PIG.df.hr.cos.bin)


    # Compare primary models ####

### unbinned

pig.df.prim.comp <- summarize_ds_models(PIG.df.uni.poly, PIG.df.hn.cos,PIG.df.hr.cos, 
                                           output = "plain")
pig.df.prim.comp[ ,1:5]
pig.df.prim.comp[ ,6:7]

# Uni model has little support (dAIC > 3.5). HR model has the most support. It is the fit that I think is the most appropriate in reality. 

# HR model selected


### binned

pig.df.prim.comp.bin <- summarize_ds_models(PIG.df.uni.cos.bin, PIG.df.hn.cos.bin,PIG.df.hr.cos.bin, 
                                        output = "plain")
pig.df.prim.comp.bin[ ,1:5]
pig.df.prim.comp.bin[ ,6:7]
# all models have some support

# plot together
par(mfrow=c(2,2))
plot(PIG.df.uni.cos.bin, main="uni")
plot(PIG.df.hn.cos.bin, main="hn")
plot(PIG.df.hr.cos.bin, main="hr")
# the HR and HN models look to be the best fit to me. Uni model doesn't fit well close to the line. 


# check estimates
pig.prim.bin.comp <- data.frame(model = rep(c("hn", "hr"), each=7),
                                year = c(PIG.df.hn.cos.bin$dht$individuals$N$Label[1:7],
                                         PIG.df.hr.cos.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(PIG.df.hn.cos.bin$dht$individuals$N$Estimate[1:7],
                                             PIG.df.hr.cos.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(PIG.df.hn.cos.bin$dht$individuals$N$cv[1:7],
                                       PIG.df.hr.cos.bin$dht$individuals$N$cv[1:7]),
                                se = c(PIG.df.hn.cos.bin$dht$individuals$N$se[1:7],
                                       PIG.df.hr.cos.bin$dht$individuals$N$se[1:7]),
                                lcl = c(PIG.df.hn.cos.bin$dht$individuals$N$lcl[1:7],
                                        PIG.df.hr.cos.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(PIG.df.hn.cos.bin$dht$individuals$N$ucl[1:7],
                                        PIG.df.hr.cos.bin$dht$individuals$N$ucl[1:7]))

pig.prim.bin.comp

# plot cv
ggplot(pig.prim.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# HR always the lowest CV

# plot estimates
ggplot(pig.prim.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# HR the better model in terms of CV and precision


    # Models with harsh truncation ####

# set trunc distance 
pig.trunc.harsh <- 40
pig.trunc.harsh.bin <- 45

### unbinned

par(mfrow=c(1,2))

# HR cos harsh
PIG.df.hr.cos.harsh <- ds(data=PIG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=pig.trunc.harsh, key="hr", adjustment="cos")
summary(PIG.df.hr.cos.harsh)
# no adjustment selected


# compare original with harsh truncation
ddf.gof(PTM.df.hr.cos.harsh$ddf, main = "LTM.df.hr.cos.harsh", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(PTM.df.hr.cos$ddf, main = "PTM.df.hr.cos",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# original model is the better fit (I think). But check estimates below

# extract and plot to see what the different models say
harsh <- data.frame(label = "harsh",
                   estimate = PTM.df.hr.cos.harsh$dht$clusters$N$Estimate[c(1:7)],
                   lcl = PTM.df.hr.cos.harsh$dht$clusters$N$lcl[c(1:7)],
                   ucl = PTM.df.hr.cos.harsh$dht$clusters$N$ucl[c(1:7)])

orig <- data.frame(label = "orig",
                    estimate = PTM.df.hr.cos$dht$clusters$N$Estimate[c(1:7)],
                    lcl = PTM.df.hr.cos$dht$clusters$N$lcl[c(1:7)],
                    ucl = PTM.df.hr.cos$dht$clusters$N$ucl[c(1:7)])

mod.check <- rbind(harsh,orig)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=2)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# no obvious improvements in precision from the harsh truncation, and the trend inference does not change. therefore I will stick with the original



### binned

# HR cosine
PIG.df.hr.cos.bin.harsh <- ds(data=PIG.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=pig.trunc.harsh.bin, key="hr", adjustment= "cos",
                        cutpoints = c(0,5,12,16,25,32,45))
summary(PIG.df.hr.cos.bin.harsh)
# key only selected

# plot together
par(mfrow=c(1,2))
plot(PIG.df.hr.cos.bin.harsh, main="harsh")
plot(PIG.df.hr.cos.bin, main="orig")


# check estimates
pig.harsh.bin.comp <- data.frame(model = rep(c("harsh", "orig"), each=7),
                                year = c(PIG.df.hr.cos.bin.harsh$dht$individuals$N$Label[1:7],
                                         PIG.df.hr.cos.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(PIG.df.hr.cos.bin.harsh$dht$individuals$N$Estimate[1:7],
                                             PIG.df.hr.cos.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(PIG.df.hr.cos.bin.harsh$dht$individuals$N$cv[1:7],
                                       PIG.df.hr.cos.bin$dht$individuals$N$cv[1:7]),
                                se = c(PIG.df.hr.cos.bin.harsh$dht$individuals$N$se[1:7],
                                       PIG.df.hr.cos.bin$dht$individuals$N$se[1:7]),
                                lcl = c(PIG.df.hr.cos.bin.harsh$dht$individuals$N$lcl[1:7],
                                        PIG.df.hr.cos.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(PIG.df.hr.cos.bin.harsh$dht$individuals$N$ucl[1:7],
                                        PIG.df.hr.cos.bin$dht$individuals$N$ucl[1:7]))

pig.harsh.bin.comp

# plot cv
ggplot(pig.harsh.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# Harsh has worst CV every year

# plot estimates
ggplot(pig.harsh.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# original model is generally more precise. Estimates all very similar except 2010

# original selected

    # Models with covariates ####

### unbinned


## cluster size
PIG.df.hr.size <- ds(data=PIG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.PIG, key="hr", formula = ~size)
summary(PIG.df.hr.size)

# plot
par(mfrow=c(1,2))

plot(PIG.df.hr.size, main = "PIG.df.hr.size")
covar.fit <- ddf.gof(PIG.df.hr.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)



## Observer
PIG.df.hr.obs <- ds(data=PIG.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.PIG, key="hr", formula = ~obs.observer)
summary(PIG.df.hr.obs)

# plot
plot(PIG.df.hr.obs, main = "PIG.df.hr.obs")
covar.fit <- ddf.gof(PIG.df.hr.obs$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
# fair bit of variation


## stratum (year)
PIG.df.hr.strat <- ds(data=PIG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.PIG, key="hr", formula = ~stratum)
summary(PIG.df.hr.strat)


# plot
plot(PIG.df.hr.strat, main = "PIG.df.hr.strat")
covar.fit <- ddf.gof(PIG.df.hr.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)


## stratum and size
PIG.df.hr.strat.size <- ds(data=PIG.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.PIG, key="hr", formula = ~stratum+size)
summary(PIG.df.hr.strat.size)


# plot
plot(PIG.df.hr.strat.size, main = "PIG.df.hr.strat.size")
covar.fit <- ddf.gof(PIG.df.hr.strat.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)



### binned


# size
PIG.df.hr.size.bin <- ds(data=PIG.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.PIG.bin, key="hr", formula = ~size,
                        cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.hr.size.bin)
plot(PIG.df.hr.size.bin)


# observer
PIG.df.hr.obs.bin <- ds(data=PIG.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table, 
                         truncation=trunc.PIG.bin, key="hr", formula = ~obs.observer,
                         cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.hr.obs.bin)
plot(PIG.df.hr.obs.bin)

# stratum
PIG.df.hr.strat.bin <- ds(data=PIG.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.PIG.bin, key="hr", formula = ~stratum,
                        cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.hr.strat.bin)
plot(PIG.df.hr.strat.bin)

# stratum + size
PIG.df.hr.strat.size.bin <- ds(data=PIG.data, region.table=full.region.table, 
                          sample.table=sample.table, obs.table=obs.table, 
                          truncation=trunc.PIG.bin, key="hr", formula = ~stratum+size,
                          cutpoints = c(0,5,12,16,25,32,45,60))
summary(PIG.df.hr.strat.size.bin)
plot(PIG.df.hr.strat.size.bin)



    # Compare covariate models ####

### unbinned

pig.df.cov.comp <- summarize_ds_models(PIG.df.hr.size,PIG.df.hr.obs,PIG.df.hr.strat,
                                       PIG.df.hr.strat.size,output = "plain")

pig.df.cov.comp[ ,1:5]
pig.df.cov.comp[ ,6:7]

# model with observer has no support.  models with size, stratum, and size+stratum have similar support

# check 3 models
summary(PIG.df.hr.size)
summary(PIG.df.hr.strat)
summary(PIG.df.hr.strat.size)

# stratum model suggest increase in p over time (coeff = 0.1)
# stratum and size model suggest increase in p over time (coeff = 0.14), and decrease in p with increased group size (have I interpreted this correctly??)

# extract and plot to see what the different models say
size <- data.frame(label = "size",
                    estimate = PIG.df.hr.size$dht$clusters$N$Estimate[c(1:7)],
                    lcl = PIG.df.hr.size$dht$clusters$N$lcl[c(1:7)],
                    ucl = PIG.df.hr.size$dht$clusters$N$ucl[c(1:7)])

strat <- data.frame(label = "strat",
                   estimate = PIG.df.hr.strat$dht$clusters$N$Estimate[c(1:7)],
                   lcl = PIG.df.hr.strat$dht$clusters$N$lcl[c(1:7)],
                   ucl = PIG.df.hr.strat$dht$clusters$N$ucl[c(1:7)])

strat.size <- data.frame(label = "strat.size",
                    estimate = PIG.df.hr.strat.size$dht$clusters$N$Estimate[c(1:7)],
                    lcl = PIG.df.hr.strat.size$dht$clusters$N$lcl[c(1:7)],
                    ucl = PIG.df.hr.strat.size$dht$clusters$N$ucl[c(1:7)])

mod.check <- rbind(size,strat,strat.size)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=3)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# all three models agree on the braod estimates, and the trend inference will not change absed on which model is selected (which is good!). It looks like the model with both covars in is slightly less precise than the other two. 

# check QQ plots
par(mfrow=c(1,3))
ddf.gof(PIG.df.hr.size$ddf, main = "PIG.df.hr.size", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(PIG.df.hr.strat$ddf, main = "PIG.df.hr.strat",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(PIG.df.hr.strat.size$ddf, main = "PIG.df.hr.strat.size",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# basically identical.

## Not really sure the best way to select the final model, as there is so little in it. The strat model tends to be the estimate in the middle, and therefore the "compromise" between the three.  This is as good a reason as any as far as I can tell! I also think that as a heavily hunted species, it is important to account for changes in detection over the years.

# PIG.df.hr.strat selected


### binned

pig.df.cov.comp.bin <- summarize_ds_models(PIG.df.hr.size.bin,PIG.df.hr.obs.bin,PIG.df.hr.strat.bin,
                                           PIG.df.hr.strat.size.bin,  output = "plain")

pig.df.cov.comp.bin[ ,1:5]
pig.df.cov.comp.bin[ ,6:7]
# observer has no support. Other models all have some support

# check estimates
pig.cov.bin.comp <- data.frame(model = rep(c("size", "strat", "strat+size"), each=7),
                                year = c(PIG.df.hr.size.bin$dht$individuals$N$Label[1:7],
                                         PIG.df.hr.strat.bin$dht$individuals$N$Label[1:7],
                                         PIG.df.hr.strat.size.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(PIG.df.hr.size.bin$dht$individuals$N$Estimate[1:7],
                                             PIG.df.hr.strat.bin$dht$individuals$N$Estimate[1:7],
                                             PIG.df.hr.strat.size.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(PIG.df.hr.size.bin$dht$individuals$N$cv[1:7],
                                       PIG.df.hr.strat.bin$dht$individuals$N$cv[1:7],
                                       PIG.df.hr.strat.size.bin$dht$individuals$N$cv[1:7]),
                                se = c(PIG.df.hr.size.bin$dht$individuals$N$se[1:7],
                                       PIG.df.hr.strat.bin$dht$individuals$N$se[1:7],
                                       PIG.df.hr.strat.size.bin$dht$individuals$N$se[1:7]),
                                lcl = c(PIG.df.hr.size.bin$dht$individuals$N$lcl[1:7],
                                        PIG.df.hr.strat.bin$dht$individuals$N$lcl[1:7],
                                        PIG.df.hr.strat.size.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(PIG.df.hr.size.bin$dht$individuals$N$ucl[1:7],
                                        PIG.df.hr.strat.bin$dht$individuals$N$ucl[1:7],
                                        PIG.df.hr.strat.size.bin$dht$individuals$N$ucl[1:7]))

pig.cov.bin.comp

# plot cv
ggplot(pig.cov.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# strat+size always the highest CV.  strat most often lowest CV but sometimes size

# plot estimates
ggplot(pig.cov.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# Stratum model consistenlty the most precise

# strat model selected


    # Compare primary and covariate models ####

## unbinned

pig.df.comp.final <- summarize_ds_models(PIG.df.hr.strat, PIG.df.hr.cos,
                                       output = "plain")

pig.df.comp.final[ ,1:5]
pig.df.comp.final[ ,6:7]

# both models have support

# extract and plot to see what the different models say
orig <- data.frame(label = "orig",
                   estimate = PIG.df.hr.cos$dht$clusters$N$Estimate[c(1:7)],
                   lcl = PIG.df.hr.cos$dht$clusters$N$lcl[c(1:7)],
                   ucl = PIG.df.hr.cos$dht$clusters$N$ucl[c(1:7)])

strat <- data.frame(label = "strat",
                    estimate = PIG.df.hr.strat$dht$clusters$N$Estimate[c(1:7)],
                    lcl = PIG.df.hr.strat$dht$clusters$N$lcl[c(1:7)],
                    ucl = PIG.df.hr.strat$dht$clusters$N$ucl[c(1:7)])

mod.check <- rbind(orig,strat)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=2)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# very similar estimates, and the same trend. As I mentioned above, I think it is important for this species to account for changes in detection over time. Therefore strat model selected


### binned


pig.df.comp.final.bin <- summarize_ds_models(PIG.df.hr.strat.bin, PIG.df.hr.cos.bin,
                                         output = "plain")

pig.df.comp.final.bin[ ,1:5]
pig.df.comp.final.bin[ ,6:7]
# both models hav some support

# check estimates
pig.final.bin.comp <- data.frame(model = rep(c("orig", "strat"), each=7),
                               year = c(PIG.df.hr.cos.bin$dht$individuals$N$Label[1:7],
                                        PIG.df.hr.strat.bin$dht$individuals$N$Label[1:7]),
                               estimate = c(PIG.df.hr.cos.bin$dht$individuals$N$Estimate[1:7],
                                            PIG.df.hr.strat.bin$dht$individuals$N$Estimate[1:7]),
                               cv = c(PIG.df.hr.cos.bin$dht$individuals$N$cv[1:7],
                                      PIG.df.hr.strat.bin$dht$individuals$N$cv[1:7]),
                               se = c(PIG.df.hr.cos.bin$dht$individuals$N$se[1:7],
                                      PIG.df.hr.strat.bin$dht$individuals$N$se[1:7]),
                               lcl = c(PIG.df.hr.cos.bin$dht$individuals$N$lcl[1:7],
                                       PIG.df.hr.strat.bin$dht$individuals$N$lcl[1:7]),
                               ucl = c(PIG.df.hr.cos.bin$dht$individuals$N$ucl[1:7],
                                       PIG.df.hr.strat.bin$dht$individuals$N$ucl[1:7]))

pig.final.bin.comp

# plot cv
ggplot(pig.final.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# original model tends to have lower CV

# plot estimates
ggplot(pig.final.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))

# estimates are similar, and the original model tends to be more precise. however, I think the fluctuations in population size of this species, and hunting pressure, means that I am inclined to use the covariate model as I think detection will almost certainly vary across years, and so stratum I think is important.

# PIG.df.hr.strat.bin selected


  ## PIG Final Results ####

### unbinned

summary(PIG.df.hr.strat)

# extract estimates
PIG.grp.results <- PIG.df.hr.strat$dht$clusters$N[1:7, ]
PIG.ind.results <- PIG.df.hr.strat$dht$individuals$N[1:7, ]

# bind results
PIG.results <- rbind(PIG.grp.results, PIG.ind.results)
PIG.results <- PIG.results %>% rename(Year = Label) %>% 
                mutate(Label = rep(c("Grp", "Ind"), each=7)) %>% 
                mutate(Species = rep("PIG", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hr", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("stratum", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)



PIG_final_plot <- ggplot(PIG.results[PIG.results$Label=="Ind",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  ylim(0,5250)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(PIG.results, "Output/Results/PIG_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/PIG_final_plot.png", PIG_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")



### binned

summary(PIG.df.hr.strat.bin)

# extract estimates
PIG.grp.abund <- PIG.df.hr.strat.bin$dht$clusters$N[1:7, ]
PIG.grp.abund <- PIG.grp.abund %>% rename(Year = Label) %>% mutate(Label="Grp")
PIG.ind.abund <- PIG.df.hr.strat.bin$dht$individuals$N[1:7, ]
PIG.ind.abund <- PIG.ind.abund %>% rename(Year=Label) %>% mutate(Label="Ind")
PIG.abund <- rbind(PIG.grp.abund,PIG.ind.abund)
PIG.abund <- PIG.abund %>% rename(N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl)
PIG.abund <- PIG.abund[,-7]

### extract density
PIG.grp.density <- PIG.df.hr.strat.bin$dht$clusters$D[1:7, ]
PIG.grp.density <- PIG.grp.density %>% rename(Year=Label) %>% mutate(Label="Grp")
PIG.ind.density <- PIG.df.hr.strat.bin$dht$individuals$D[1:7, ]
PIG.ind.density <- PIG.ind.density %>% rename(Year=Label) %>% mutate(Label="Ind")
PIG.density <- rbind(PIG.grp.density,PIG.ind.density)
PIG.density <- PIG.density %>% rename(D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl)
PIG.density <- PIG.density[,-7]
PIG.density <- PIG.density[,-c(1,7)]


# bind results
PIG.results <- cbind(PIG.abund, PIG.density)
PIG.results <- PIG.results %>%  
                mutate(Species = rep("PIG", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hr", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("year", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,
                       N,n_se,n_cv,n_lcl,n_ucl,
                       D,d_se,d_cv,d_lcl,d_ucl)

PIG.results <- PIG.results %>% mutate(D=D*1000000,d_se=d_se*1000000,d_lcl=d_lcl*1000000,d_ucl=d_ucl*1000000)





PIG_final_plot_bin <- ggplot(PIG.results.bin[PIG.results.bin$Label=="Ind",], aes(x=Year, y=Estimate))+
                    geom_point(shape=16, size=2)+
                    geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                    ylim(0,6200)+
                    theme_bw()+
                    theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),panel.border = element_blank())+
                    theme(axis.line = element_line(color = 'black'))

# save results
write.csv(PIG.results, "Output/Results/PIG_results_final_binned.csv")
ggsave("Output/Results/Plots/Point_estimates/PIG_final_plot_binned.png", PIG_final_plot_bin, 
       dpi=300, width = 20, height = 20, units = "cm")


# compare with original resutls
PIG.results <- read.csv("Output/Results/PIG_results_final.csv")
PIG.results <- PIG.results[,-1]
PIG.results$analysis <- "original"

PIG.results.bin$analysis <- "binned"

PIG.results.all <- rbind(PIG.results,PIG.results.bin)

PIG_analysis_comparison_plot <- ggplot(PIG.results.all[PIG.results.all$Label=="Ind",], 
                                       aes(x=Year, y=Estimate, colour=analysis))+
                          geom_point(shape=16, size=2, position = position_dodge(width=0.3))+
                          geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))+
                          ylim(0,6100)+
                          theme_bw()+
                          theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),panel.border = element_blank())+
                          theme(axis.line = element_line(color = 'black'))

ggsave("Output/Results/Plots/Point_estimates/PIG_analysis_comparison_plot.png", 
       PIG_analysis_comparison_plot, dpi=300, width = 20, height = 20, units = "cm")


#### Green peafowl ###########################################################
  ## Subset data ####

# subset GPF data
GPF.data <- allData[allData$species=="GPF",] 
GPF.data <- as.data.frame(GPF.data)
head(GPF.data)

# Total number of groups from all years
length(GPF.data$distance) 

# for GPF, I will pool all years for DF, and will test year as a continous covariate in the DF model

# check 2020 data is there
unique(GPF.data$year)

  ## Exploratory plots and linear models ####

par(mfrow=c(1,2))

# distance histograms
hist(GPF.data$distance[GPF.data$distance<100], main="All years", xlab="Distance (m)")

# More bins to see better what is happening around 0
hist(GPF.data$distance[GPF.data$distance<100], main="All years", xlab="Distance (m)", breaks=c(40))

# I believe the larger distances - this species prefers open and mixed forest where visibility is better. interesting "step" looking shape - between 0-30, 30-60, and 60-90. I think for this species p is high clsoe to the line, as the data suggest. It is noisy and flighty, and is fairly easy to detect with confidence when it flees. I wonder whethe the next 'step' is because of all the observations in the open habitat, where visibility is greater. 


# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.GPF <- 80

# Count the number of observations discarded
nrow(GPF.data[GPF.data$distance>trunc.GPF,]) 

length(GPF.data$distance)

15/167*100 # 9% will defo test harsher truncation 

## Plots of covars against distance

# Plot of distance against cluster size
par(mfrow=c(1,2))

plot(GPF.data$size, GPF.data$distance, 
     main="size", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.GPF <- lm(distance~size, data=GPF.data)
lines(GPF.data$size, as.vector(predict(lm.GPF, GPF.data)))
summary(lm.GPF)
# some evidence of size bias but a weak model


# Plot of Observer factor against distance
plot(GPF.data$obs.observer,GPF.data$distance, main="observer", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# some variation in observers

# I can't include habitat or AMPM in the DF because those variables weren't recorded in earlier years


  ## Fit a detection function ####
    # Uniform ####

# Uni cosine

GPF.df.uni.cos <- ds(data=GPF.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GPF, key="uni", adjustment= "cos")
summary(GPF.df.uni.cos)
# cosine(1)  

# Uni poly
GPF.df.uni.poly <- ds(data=GPF.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GPF, key="uni", adjustment= "poly")
summary(GPF.df.uni.poly)
# poly(2,4)

## compare uni models
gpf.uni.comp <- summarize_ds_models(GPF.df.uni.cos, GPF.df.uni.poly,
                                    output="plain")
gpf.uni.comp[ ,1:4]
gpf.uni.comp[ ,5:7]
# cos is preferred (poly dAIC > 2)

# plot the fits
par(mfrow=c(2,2))

# cos
plot(GPF.df.uni.cos, main = "GPF.df.uni.cos")

covar.fit <- ddf.gof(GPF.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(GPF.df.uni.poly, main = "GPF.df.uni.poly")

covar.fit <- ddf.gof(GPF.df.uni.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# Look pretty identical to me. I will go with AIC

# cos selected

    # Half normal ####

# HN cosine
GPF.df.hn.cos <- ds(data=GPF.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.GPF, key="hn", adjustment= "cos")
summary(GPF.df.hn.cos)
# no adjustment selected

# HN hermite
GPF.df.hn.herm <- ds(data=GPF.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.GPF, key="hn", adjustment= "herm")
summary(GPF.df.hn.herm)
# no adjustment selected


# plot the fit
par(mfrow=c(1,2))

# cos
plot(GPF.df.hn.cos, main = "GPF.df.hn.cos")

covar.fit <- ddf.gof(GPF.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# cos is selected. Fairly similar to uni cos

    # Hazard rate ####

# HR cosine
GPF.df.hr.cos <- ds(data=GPF.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GPF, key="hr", adjustment= "cos")
summary(GPF.df.hr.cos)
# no adjustment selected


# HR poly
GPF.df.hr.poly <- ds(data=GPF.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.GPF, key="hr", adjustment= "poly")
summary(GPF.df.hr.poly)
# no adjustment selected 


# plot the fit
par(mfrow=c(1,2))

# hr cos
plot(GPF.df.hr.cos, main = "GPF.df.hr.cos")

covar.fit <- ddf.gof(GPF.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# wide shoulder - p=1 up to 10m, then steeper drop off than uni and hn.

    # Compare primary models ####

gpf.df.prim.comp <- summarize_ds_models(GPF.df.uni.cos, GPF.df.hn.cos,GPF.df.hr.cos, 
                                           output = "plain")
gpf.df.prim.comp[ ,1:5]
gpf.df.prim.comp[ ,6:7]

# All models have some support.  HR has the least with dAIC > 2 (only just).

# plot the fits together
par(mfrow=c(3,2))

# Uni cos
plot(GPF.df.uni.cos, main = "GPF.df.uni.cos")

covar.fit <- ddf.gof(GPF.df.uni.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# HN cos
plot(GPF.df.hn.cos, main = "GPF.df.hn.cos")

covar.fit <- ddf.gof(GPF.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# hr cos
plot(GPF.df.hr.cos, main = "GPF.df.hr.cos")

covar.fit <- ddf.gof(GPF.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)

# When plotted together, we can see that they are very similar. Uni and Hn are super similar. HR just assumes p=1 up to 20m, and then detection drops off quicker (p=0.5 ~ 40m, whereas the other two are ~45m). I think UNi or Hn are more realistic, and seeing as I want to test covars I will selected HN

# HN selected

    # Models with harsh truncation ####

# I would like to test two different harsh truncations 

# set trunc distance 
gpf.trunc.harsh1 <- 60
gpf.trunc.harsh2 <- 40

# HN cos harsh1
GPF.df.hn.cos.harsh1 <- ds(data=GPF.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=gpf.trunc.harsh1, key="hn", adjustment="cos")
summary(GPF.df.hn.cos.harsh1)
# no adjustment selected

# HN cos harsh2
GPF.df.hn.cos.harsh2 <- ds(data=GPF.data, region.table=full.region.table, 
                          sample.table=sample.table, obs.table=obs.table,
                          truncation=gpf.trunc.harsh2, key="hn", adjustment="cos")
summary(GPF.df.hn.cos.harsh)
# no adjustment selected


# compare original with harsh truncation
par(mfrow=c(1,3))

ddf.gof(GPF.df.hn.cos.harsh1$ddf, main = "GPF.df.hn.cos.harsh1", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(GPF.df.hn.cos.harsh2$ddf, main = "GPF.df.hn.cos.harsh2", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))

ddf.gof(GPF.df.hn.cos$ddf, main = "GPF.df.hn.cos",
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# not much in it. I think harsh2 is the worse

# extract and plot to see what the different models say
orig <- data.frame(label = "orig",
                   estimate = GPF.df.hn.cos$dht$clusters$N$Estimate[c(1:7)],
                   lcl = GPF.df.hn.cos$dht$clusters$N$lcl[c(1:7)],
                   ucl = GPF.df.hn.cos$dht$clusters$N$ucl[c(1:7)])

harsh1 <- data.frame(label = "harsh1",
                    estimate = GPF.df.hn.cos.harsh1$dht$clusters$N$Estimate[c(1:7)],
                    lcl = GPF.df.hn.cos.harsh1$dht$clusters$N$lcl[c(1:7)],
                    ucl = GPF.df.hn.cos.harsh1$dht$clusters$N$ucl[c(1:7)])

harsh2 <- data.frame(label = "harsh2",
                         estimate = GPF.df.hn.cos.harsh2$dht$clusters$N$Estimate[c(1:7)],
                         lcl = GPF.df.hn.cos.harsh2$dht$clusters$N$lcl[c(1:7)],
                         ucl = GPF.df.hn.cos.harsh2$dht$clusters$N$ucl[c(1:7)])

mod.check <- rbind(orig,harsh1,harsh2)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=3)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# Estimates are all similar. Harsh2 hs the most imprecise, whereas the original is the most precise. original selected


    # Models with covariates ####

## cluster size
GPF.df.hn.size <- ds(data=GPF.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.GPF, key="hn", formula = ~size)
summary(GPF.df.hn.size)

# plot
par(mfrow=c(1,2))

plot(GPF.df.hn.size, main = "GPF.df.hn.size")
covar.fit <- ddf.gof(GPF.df.hn.size$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
# the large distances potentially causing some issues here


## Observer
GPF.df.hn.obs <- ds(data=GPF.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.GPF, key="hn", formula = ~obs.observer)
summary(GPF.df.hn.obs)

# plot
plot(GPF.df.hn.obs, main = "GPF.df.hn.obs")
covar.fit <- ddf.gof(GPF.df.hn.obs$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
# interesting split


## stratum (year)
GPF.df.hn.strat <- ds(data=GPF.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                    truncation=trunc.GPF, key="hn", formula = ~stratum)
summary(GPF.df.hn.strat)


# plot
plot(GPF.df.hn.strat, main = "GPF.df.hn.strat")
covar.fit <- ddf.gof(GPF.df.hn.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
# doesn't look like much change over time


## stratum and size
GPF.df.hn.strat.size <- ds(data=GPF.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table,
                      truncation=trunc.GPF, key="hn", formula = ~stratum+size)
summary(GPF.df.hn.strat.size)


# plot
plot(GPF.df.hn.strat.size, main = "GPF.df.hn.strat.size")
covar.fit <- ddf.gof(GPF.df.hn.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)


    # Compare covariate models ####

gpf.df.cov.comp <- summarize_ds_models(GPF.df.hn.size,GPF.df.hn.obs,GPF.df.hn.strat,
                                       GPF.df.hn.strat.size, output = "plain")

gpf.df.cov.comp[ ,1:5]
gpf.df.cov.comp[ ,6:7]

# model with observer has no support. strat model and strat.size model have dAIC around 2.

# extract and plot to see what the different models say
size <- data.frame(label = "size",
                   estimate = GPF.df.hn.size$dht$clusters$N$Estimate[c(1:7)],
                   lcl = GPF.df.hn.size$dht$clusters$N$lcl[c(1:7)],
                   ucl = GPF.df.hn.size$dht$clusters$N$ucl[c(1:7)])

strat <- data.frame(label = "strat",
                     estimate = GPF.df.hn.strat$dht$clusters$N$Estimate[c(1:7)],
                     lcl = GPF.df.hn.strat$dht$clusters$N$lcl[c(1:7)],
                     ucl = GPF.df.hn.strat$dht$clusters$N$ucl[c(1:7)])

strat.size <- data.frame(label = "strat.size",
                     estimate = GPF.df.hn.strat.size$dht$clusters$N$Estimate[c(1:7)],
                     lcl = GPF.df.hn.strat.size$dht$clusters$N$lcl[c(1:7)],
                     ucl = GPF.df.hn.strat.size$dht$clusters$N$ucl[c(1:7)])

mod.check <- rbind(size,strat,strat.size)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=3)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))

# estimates from the three models are virtually identical, with very similar precision. as the covars are there to reduce bias, and there is no good reason to exclude them, I will use the model with both

# GPF.df.hn.strat.size selected

    # Compare primary and covariate models ####

gpf.df.comp.final <- summarize_ds_models(GPF.df.hn.strat.size, GPF.df.hn.cos,
                                       output = "plain")

gpf.df.comp.final[ ,1:5]
gpf.df.comp.final[ ,6:7]

# both models have support, although original model has more. 

# plot together
orig <- data.frame(label = "orig",
                    estimate = GPF.df.hn.cos$dht$clusters$N$Estimate[c(1:7)],
                    lcl = GPF.df.hn.cos$dht$clusters$N$lcl[c(1:7)],
                    ucl = GPF.df.hn.cos$dht$clusters$N$ucl[c(1:7)])

strat.size <- data.frame(label = "strat.size",
                         estimate = GPF.df.hn.strat.size$dht$clusters$N$Estimate[c(1:7)],
                         lcl = GPF.df.hn.strat.size$dht$clusters$N$lcl[c(1:7)],
                         ucl = GPF.df.hn.strat.size$dht$clusters$N$ucl[c(1:7)])

mod.check <- rbind(orig,strat.size)
mod.check$year <- rep(c("10","11","13","14","16","18","20"), times=2)

ggplot(mod.check, aes(x=year, y=estimate, group=label, colour=label))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3,position = position_dodge(width = 0.5))
# virtually identical. covar model provide more precise estiamtes most of the time, and so is selected

  ## GPF Final Results ####

summary(GPF.df.hn.strat.size)

# extract estimates
GPF.grp.abund <- GPF.df.hn.strat.size$dht$clusters$N[1:7, ]
GPF.grp.abund <- GPF.grp.abund %>% rename(Year = Label) %>% mutate(Label="Grp")
GPF.ind.abund <- GPF.df.hn.strat.size$dht$individuals$N[1:7, ]
GPF.ind.abund <- GPF.ind.abund %>% rename(Year=Label) %>% mutate(Label="Ind")
GPF.abund <- rbind(GPF.grp.abund,GPF.ind.abund)
GPF.abund <- GPF.abund %>% rename(N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl)
GPF.abund <- GPF.abund[,-7]

### extract density
GPF.grp.density <- GPF.df.hn.strat.size$dht$clusters$D[1:7, ]
GPF.grp.density <- GPF.grp.density %>% rename(Year=Label) %>% mutate(Label="Grp")
GPF.ind.density <- GPF.df.hn.strat.size$dht$individuals$D[1:7, ]
GPF.ind.density <- GPF.ind.density %>% rename(Year=Label) %>% mutate(Label="Ind")
GPF.density <- rbind(GPF.grp.density,GPF.ind.density)
GPF.density <- GPF.density %>% rename(D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl)
GPF.density <- GPF.density[,-7]
GPF.density <- GPF.density[,-c(1,7)]


# bind results
GPF.results <- cbind(GPF.abund, GPF.density)
GPF.results <- GPF.results %>%  
                mutate(Species = rep("GPF", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hn", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("year + size", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,
                       N,n_se,n_cv,n_lcl,n_ucl,
                       D,d_se,d_cv,d_lcl,d_ucl)

GPF.results <- GPF.results %>% mutate(D=D*1000000,d_se=d_se*1000000,d_lcl=d_lcl*1000000,d_ucl=d_ucl*1000000)




GPF_final_plot <- ggplot(GPF.results[GPF.results$Label=="Ind",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  ylim(0,2800)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(GPF.results, "Output/Results/GPF_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/GPF_final_plot.png", GPF_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")


#### Red muntjac #############################################################
  ## Subset data ####

### Previously I had analysed RMJ on an annual basis. Based on the continued decrease in annual observations, I am now going to analyse all years together. I have left the annual code below, so that in the future someone is welcome to have another crack at it.

# subset RMJ data
RMJ.data <- allData[allData$species=="RED",] 
RMJ.data <- as.data.frame(RMJ.data)
head(RMJ.data)

# Total number of groups from all years
length(RMJ.data$distance) 

# number of obs per year
length(RMJ.data$distance[RMJ.data$stratum==1]) # 169
length(RMJ.data$distance[RMJ.data$stratum==2]) # 175
length(RMJ.data$distance[RMJ.data$stratum==3]) # 167
length(RMJ.data$distance[RMJ.data$stratum==4]) # 181
length(RMJ.data$distance[RMJ.data$stratum==5]) # 93
length(RMJ.data$distance[RMJ.data$stratum==6]) # 58

# sufficient number of observations to do annual detection functions.  If 2018 data isn't great, we could always pool with one other year in get a better DF

# plot all histograms
par(mfrow=c(4,3))
hist(RMJ.data$distance, main="all yrs", xlab="Distance (m)")
hist(RMJ.data$distance, main="all yrs", xlab="Distance (m)", breaks=c(40))
hist(RMJ.2018.data$distance, main="2018", xlab="Distance (m)")
hist(RMJ.2018.data$distance, main="2018", xlab="Distance (m)", breaks=c(40))
hist(RMJ.2016.data$distance, main="2016", xlab="Distance (m)")
hist(RMJ.2016.data$distance, main="2016", xlab="Distance (m)", breaks=c(40))
hist(RMJ.2014.data$distance, main="2014", xlab="Distance (m)")
hist(RMJ.2014.data$distance, main="2014", xlab="Distance (m)", breaks=c(40))
hist(RMJ.2011.data$distance, main="2011", xlab="Distance (m)")
hist(RMJ.2011.data$distance, main="2011", xlab="Distance (m)", breaks=c(40))
hist(RMJ.2010.data$distance, main="2010", xlab="Distance (m)")
hist(RMJ.2010.data$distance, main="2010", xlab="Distance (m)", breaks=c(40))

#### POOLED ANALYSIS ####
  ## Exploratory plots & linear model ####

par(mfrow=c(1,2))

# distance histograms (remove outliers)
hist(RMJ.data$distance[RMJ.data$distance<100], main=NULL, xlab="Distance (m)")

# More bins to see better what is happening around 0 (remove outliers)
hist(RMJ.data$distance[RMJ.data$distance<100], main=NULL, xlab="Distance (m)", breaks=c(40))

# pretty nice data, apart from large spike right on the line

### UPDATE - due to lunping at 0, we are going to run analysis with bins

# identify appropritae bins (assuming truncation at 60)
hist(RMJ.data$distance[RMJ.data$distance<60], breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60))
hist(RMJ.data$distance[RMJ.data$distance<60], breaks=c(0,2,5,9,12,15,18,20,23,26,30,35,40,45,50,55,60), freq=T)

hist(RMJ.data$distance[RMJ.data$distance<60], breaks=c(0,1,4,7,10,12,15,18,21,24,27,30,33,36,39,43,48,53,58,60), 
     freq=T)

hist(RMJ.data$distance[RMJ.data$distance<60], breaks=c(0,1,4,7,10,12,15,18,21,24,28,31,35,38,41,44,48,50,53,
                                                       56,58,60), freq=T)

hist(RMJ.data$distance[RMJ.data$distance<60], breaks=c(0,0.5,4,7,10,12,15,18,21,24,28,31,35,38,41,44,48,50,53,
                                                       56,58,60), freq=T)

hist(RMJ.data$distance[RMJ.data$distance<60], breaks=c(0,7,14,21,28,35,42,49,56,60), freq=T)

# pretty challenging to get appropriate bins!  I think the last one jsut above is the best I can get it - it can't be perfect as it's real data! The first bin is still quite high, 

# Save the chosen truncation distance for later use. Try harsher truncation to shrink CIs later
trunc.RMJ <- 60

# Count the number of observations discarded
nrow(RMJ.data[RMJ.data$distance>trunc.RMJ,]) 

length(RMJ.data$distance)

67/878*100 # 7.6%

# Plot of distance against cluster size
par(mfrow=c(1,2))

plot(RMJ.data$size, RMJ.data$distance, main="a)", xlab="Group size",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))

# Fit a linear model
lm.RMJ <- lm(distance~size, data=RMJ.data)
lines(RMJ.data$size, as.vector(predict(lm.RMJ, RMJ.data)))
summary(lm.RMJ)

# testing for size bias isn't appropriate really as the species is (mostly) solitary


# Plot of Observer factor against distance
plot(RMJ.data$obs.observer,RMJ.data$distance, main="b)", xlab="Observer",
     ylab="Distance (m)", pch=19, cex=0.5, col=gray(0.7))
# some variation in observers

## now I am pooling data I can't plot AMPM or habitat as earlier years don't have them

  ## Fit a detection function ####
    # Uniform ####


### unbinned


# uniform cosine
RMJ.df.unif.cos <- ds(data=RMJ.data, region.table=full.region.table, 
                         sample.table=sample.table, obs.table=obs.table, 
                         truncation=trunc.RMJ, key="unif", adjustment= "cos")
summary(RMJ.df.unif.cos)
# cosine(1,1) selected


# uniform poly
RMJ.df.unif.poly <- ds(data=RMJ.data, region.table=full.region.table, 
                       sample.table=sample.table, obs.table=obs.table,  
                          truncation=trunc.RMJ, key="unif", adjustment= "poly")
summary(RMJ.df.unif.poly)
# simple polynomial(2,4,6) selected

# compare uniform models
rmj.uni.comp <- summarize_ds_models(RMJ.df.unif.cos,RMJ.df.unif.poly,
                                       output="plain")
rmj.uni.comp[ ,1:4]
rmj.uni.comp[ ,5:7]
#  poly dACI > 4


# plot the fits
par(mfrow=c(2,2))

# cos
plot(RMJ.df.unif.cos, main = "RMJ.df.unif.cos")

covar.fit <- ddf.gof(RMJ.df.unif.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# poly
plot(RMJ.df.unif.poly, main = "RMJ.df.unif.poly")

covar.fit <- ddf.gof(RMJ.df.unif.poly$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# look identical to me, so I will go with AIC
# cos selected


### binned

# uniform cosine
RMJ.df.unif.cos.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                      sample.table=sample.table, obs.table=obs.table, 
                      truncation=trunc.RMJ, key="unif", adjustment= "cos",
                      cutpoints = c(0,7,14,21,28,35,42,49,56,60))
summary(RMJ.df.unif.cos.bin)
# cosine(1,2,3) selected

# uniform poly
RMJ.df.unif.poly.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                          sample.table=sample.table, obs.table=obs.table, 
                          truncation=trunc.RMJ, key="unif", adjustment= "poly",
                          cutpoints = c(0,7,14,21,28,35,42,49,56,60))
summary(RMJ.df.unif.poly.bin)
# poly(2,4,6) selected


# compare uniform models
rmj.uni.comp.bin <- summarize_ds_models(RMJ.df.unif.cos.bin,RMJ.df.unif.poly.bin,
                                    output="plain")
rmj.uni.comp.bin[ ,1:4]
rmj.uni.comp.bin[ ,5:7]
# poly dAIC > 4

# plot them both
par(mfrow=c(1,2))
plot(RMJ.df.unif.cos.bin, main="cos")
plot(RMJ.df.unif.poly.bin, main="poly")

# RMJ.df.unif.cos.bin selected

    # Half-normal ####

### unbinned

# half-normal cosine
RMJ.df.hn.cos <- ds(data=RMJ.data, region.table=full.region.table, 
                       sample.table=sample.table, obs.table=obs.table, 
                       truncation=trunc.RMJ, key="hn", adjustment= "cos")
summary(RMJ.df.hn.cos)
# cosine(2)

# half-normal hermite 
RMJ.df.hn.herm <- ds(data=RMJ.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.RMJ, key="hn", adjustment= "herm")
summary(RMJ.df.hn.herm)
# no adjustment selected

# compare hn models
rmj.hn.comp <- summarize_ds_models(RMJ.df.hn.cos,RMJ.df.hn.herm,
                                    output="plain")
rmj.hn.comp[ ,1:4]
rmj.hn.comp[ ,5:7]
# herm model has CvM p value of 0.002 and so is rejected

# plot the fit
par(mfrow=c(1,2))

# cos
plot(RMJ.df.hn.cos, main = "RMJ.df.hn.cos")

covar.fit <- ddf.gof(RMJ.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# Wasn't keen on this fit at first glance. I thought p was dropping off too quickly, but actually, this is probably fairly true of this species. It is small and shy, and in thicker forest there is a good chance you won't get a look at it to be able to ID it.



### binned

# half-normal cosine
RMJ.df.hn.cos.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                    truncation=trunc.RMJ, key="hn", adjustment= "cos",
                    cutpoints = c(0,7,14,21,28,35,42,49,56,60))
summary(RMJ.df.hn.cos.bin)
# cosine(2)


# half-normal herminte
RMJ.df.hn.herm.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.RMJ, key="hn", adjustment= "herm",
                        cutpoints = c(0,7,14,21,28,35,42,49,56,60))
summary(RMJ.df.hn.cos.bin)
# no adjustment selected

# compare hn models
rmj.hn.comp.bin <- summarize_ds_models(RMJ.df.hn.cos.bin,RMJ.df.hn.herm.bin,
                                   output="plain")
rmj.hn.comp.bin[ ,1:4]
rmj.hn.comp.bin[ ,5:7]
# hermite dAIC = 14

# plot both
plot(RMJ.df.hn.cos.bin, main="cos")
plot(RMJ.df.hn.herm.bin, main="key")
# cos is doing a better job at fitting the weird shape in the data. Its the steep drop into a shoulder between 20m and 40m that is the problem.

# RMJ.df.hn.cos.bin selected


    # Hazard rate ####


### unbinned

# HR cosine
RMJ.df.hr.cos <- ds(data=RMJ.data, region.table=full.region.table, 
                       sample.table=sample.table, obs.table=obs.table, 
                       truncation=trunc.RMJ, key="hr", adjustment= "cos")
summary(RMJ.df.hr.cos)
# no adjustment selected

# HR poly
RMJ.df.hr.poly <- ds(data=RMJ.data, region.table=full.region.table, 
                     sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.RMJ, key="hr", adjustment= "poly")
summary(RMJ.df.hr.poly)
# no adjustment selected

# model with no adjustment selected

# Plot
par(mfrow=c(1,2))

plot(RMJ.df.hr.cos, main = "RMJ.df.hr.cos")
covar.fit <- ddf.gof(RMJ.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))

message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)

text(0.7, 0.1, message, cex=0.8)

covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# This feels like it fits the data better to me. Shoulder to just under 10m then steeper drop off


### binned

# HR cosine
RMJ.df.hr.cos.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                    truncation=trunc.RMJ, key="hr", adjustment= "cos",
                    cutpoints = c(0,7,14,21,28,35,42,49,56,60))
summary(RMJ.df.hr.cos.bin)
# key only

# HR poly
RMJ.df.hr.poly.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.RMJ, key="hr", adjustment= "poly",
                        cutpoints = c(0,7,14,21,28,35,42,49,56,60))
summary(RMJ.df.hr.poly.bin)
# key only

# plot
par(mfrow=c(1,1))
plot(RMJ.df.hr.cos.bin)
# not a bad looking fit to me


    # Compare primary models ####


### unbinned

rmj.df.prim.comp <- summarize_ds_models(RMJ.df.unif.cos, RMJ.df.hn.cos, 
                                           RMJ.df.hr.cos, 
                                           output = "plain")
rmj.df.prim.comp[ ,1:5]
rmj.df.prim.comp[ ,6:7]
# CvM p values all good. All models have support. That means this is essentially an ecolgical decision. DO we think that p is close to 1 up until 20m?  Or do we think the species is able to sneak away or hide that close to the line?  

# plot all fits together
par(mfrow=c(3,2))

# uni cos
plot(RMJ.df.unif.cos, main = "RMJ.df.unif.cos")
covar.fit <- ddf.gof(RMJ.df.unif.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# Hn cos
plot(RMJ.df.hn.cos, main = "RMJ.df.hn.cos")
covar.fit <- ddf.gof(RMJ.df.hn.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM

# Hr cos
plot(RMJ.df.hr.cos, main = "RMJ.df.hr.cos")
covar.fit <- ddf.gof(RMJ.df.hr.cos$ddf, lwd = 2, lty = 1, pch = ".", cex = 3,col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
covar.fit$dsgof$ks
covar.fit$dsgof$CvM


# so in fact, uni and hn aren't quite the same. Uniform maintains more of a shoulder - 0.8 is just before 20m, whereas HN 0.8 is aroun 10m.  Uni reacehs 0.5 just before 30m, HN and HR around 25m. the problem is that the different fits are more appropriate in the different habitat types. I think the uniform model is the best middle ground  

# uniform selected


### binned

rmj.df.prim.comp.bin <- summarize_ds_models(RMJ.df.unif.cos.bin, RMJ.df.hn.cos.bin, 
                                        RMJ.df.hr.cos.bin, 
                                        output = "plain")
rmj.df.prim.comp.bin[ ,1:5]
rmj.df.prim.comp.bin[ ,6:7]
# all models have some support

# plot together
par(mfrow=c(2,2))
plot(RMJ.df.unif.cos.bin, main="unif")
plot(RMJ.df.hn.cos.bin, main="hn")
plot(RMJ.df.hr.cos.bin, main="hr")

# check estimates
rmj.prim.bin.comp <- data.frame(model = rep(c("uni","hn", "hr"), each=7),
                                year = c(RMJ.df.unif.cos.bin$dht$individuals$N$Label[1:7],
                                         RMJ.df.hn.cos.bin$dht$individuals$N$Label[1:7],
                                         RMJ.df.hr.cos.bin$dht$individuals$N$Label[1:7]),
                                  estimate = c(RMJ.df.unif.cos.bin$dht$individuals$N$Estimate[1:7],
                                               RMJ.df.hn.cos.bin$dht$individuals$N$Estimate[1:7],
                                               RMJ.df.hr.cos.bin$dht$individuals$N$Estimate[1:7]),
                                  cv = c(RMJ.df.unif.cos.bin$dht$individuals$N$cv[1:7],
                                         RMJ.df.hn.cos.bin$dht$individuals$N$cv[1:7],
                                         RMJ.df.hr.cos.bin$dht$individuals$N$cv[1:7]),
                                  se = c(RMJ.df.unif.cos.bin$dht$individuals$N$se[1:7],
                                         RMJ.df.hn.cos.bin$dht$individuals$N$se[1:7],
                                         RMJ.df.hr.cos.bin$dht$individuals$N$se[1:7]),
                                  lcl = c(RMJ.df.unif.cos.bin$dht$individuals$N$lcl[1:7],
                                          RMJ.df.hn.cos.bin$dht$individuals$N$lcl[1:7],
                                          RMJ.df.hr.cos.bin$dht$individuals$N$lcl[1:7]),
                                  ucl = c(RMJ.df.unif.cos.bin$dht$individuals$N$ucl[1:7],
                                          RMJ.df.hn.cos.bin$dht$individuals$N$ucl[1:7],
                                          RMJ.df.hr.cos.bin$dht$individuals$N$ucl[1:7]))

rmj.prim.bin.comp

# plot cv
ggplot(rmj.prim.bin.comp, aes(x=year, y=cv, group=model, colour=model))+
  geom_point()
# HN model lowest CV every year

# plot estimates
ggplot(rmj.prim.bin.comp, aes(x=year, y=estimate, group=model, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# not much in it. Precision is similar. Uni and HN basically identical. Based on CV I will select HN

# RMJ.df.hn.cos.bin selected


    # Models with harsh truncation ####

# set trunc distance to 40m
trunc.harsh.rmj <- 40
trunc.harsh.rmj.bin <- 49


### unbinned


par(mfrow=c(1,2))

# uni cos
RMJ.df.uni.cos.harsh <- ds(data=RMJ.data, region.table=full.region.table, 
                              sample.table=sample.table, obs.table=obs.table, 
                              truncation=trunc.harsh.rmj, key="uni", adjustment= "cos")
summary(RMJ.df.uni.cos.harsh)
# cosine(1,2)

# compare original with harsh truncation 
ddf.gof(RMJ.df.uni.cos.harsh$ddf, main = "RMJ.df.uni.cos.harsh", 
        lwd = 2, lty = 1, pch = ".", 
        cex = 3, col = c(1,2))

ddf.gof(RMJ.df.unif.cos$ddf, main = "RMJ.df.unif.cos", 
        lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
# Not much in it but I think original is the better fit

# original selected


### binned

RMJ.df.hn.cos.harsh.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.harsh.rmj.bin, key="hn",adjustment = "cos",
                        cutpoints = c(0,7,14,21,28,35,42,49))
summary(RMJ.df.hn.cos.harsh.bin)
# cosine(2)

# plot together
par(mfrow=c(1,2))
plot(RMJ.df.hn.cos.harsh.bin, main="harsh")
plot(RMJ.df.hn.cos.bin, main="orig")

# compare estimates
rmj.harsh.bin.comp <- data.frame(model = rep(c("harsh","orig"), each=7),
                                year = c(RMJ.df.hn.cos.harsh.bin$dht$individuals$N$Label[1:7],
                                         RMJ.df.hn.cos.bin$dht$individuals$N$Label[1:7]),
                                estimate = c(RMJ.df.hn.cos.harsh.bin$dht$individuals$N$Estimate[1:7],
                                             RMJ.df.hn.cos.bin$dht$individuals$N$Estimate[1:7]),
                                cv = c(RMJ.df.hn.cos.harsh.bin$dht$individuals$N$cv[1:7],
                                       RMJ.df.hn.cos.bin$dht$individuals$N$cv[1:7]),
                                se = c(RMJ.df.hn.cos.harsh.bin$dht$individuals$N$se[1:7],
                                       RMJ.df.hn.cos.bin$dht$individuals$N$se[1:7]),
                                lcl = c(RMJ.df.hn.cos.harsh.bin$dht$individuals$N$lcl[1:7],
                                        RMJ.df.hn.cos.bin$dht$individuals$N$lcl[1:7]),
                                ucl = c(RMJ.df.hn.cos.harsh.bin$dht$individuals$N$ucl[1:7],
                                        RMJ.df.hn.cos.bin$dht$individuals$N$ucl[1:7]))
rmj.harsh.bin.comp

# plot cv
ggplot(rmj.harsh.bin.comp, aes(x=year, y=cv, colour=model))+
  geom_point()
# variation in which model is better

# plot estimates
ggplot(rmj.harsh.bin.comp, aes(x=year, y=estimate, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))
# not a simple decision. Each model performs better in different years. Thankfully the overall differences in results are small, and our inferences won't change. All other things being equal, I will stick with the model that uses the most data

# original selected


    # Models with covariates ####

### unbinned


# I will not test size here as the species is solitary. I also can't test habitat or AMPM as before, as the data are now pooled


## Observer
RMJ.df.hn.obs <- ds(data=RMJ.data, region.table=full.region.table, 
                       sample.table=sample.table, obs.table=obs.table, 
                       truncation=trunc.RMJ, key="hn", formula = ~obs.observer)
summary(RMJ.df.hn.obs)

# plot
plot(RMJ.df.hn.obs, main = "RMJ.df.hn.obs")
covar.fit <- ddf.gof(RMJ.df.hn.obs$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)
# not a great fit


## stratum
RMJ.df.hn.strat <- ds(data=RMJ.data, region.table=full.region.table, 
                    sample.table=sample.table, obs.table=obs.table, 
                    truncation=trunc.RMJ, key="hn", formula = ~stratum)
summary(RMJ.df.hn.strat)

# plot
plot(RMJ.df.hn.strat, main = "RMJ.df.hn.strat")
covar.fit <- ddf.gof(RMJ.df.hn.strat$ddf, lwd = 2, lty = 1, pch = ".", cex = 3, col = c(1,2))
message <- paste("CvM W=", round(covar.fit$dsgof$CvM$W,3), 
                 "(P=", round(covar.fit$dsgof$CvM$p,2),") ", "\nx2", 
                 round( covar.fit$chisquare$chi1$chisq, 1), " df", covar.fit$chisquare$chi1$df)
text(0.7, 0.1, message, cex=0.8)



### binned

# observer
RMJ.df.hn.obs.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.RMJ, key="hn", formula = ~obs.observer,
                        cutpoints = c(0,7,14,21,28,35,42,49,56,60))
summary(RMJ.df.hn.obs.bin)
plot(RMJ.df.hn.obs.bin)

# stratum
RMJ.df.hn.strat.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                        sample.table=sample.table, obs.table=obs.table, 
                        truncation=trunc.RMJ, key="hn", formula = ~stratum,
                        cutpoints = c(0,7,14,21,28,35,42,49,56,60))
summary(RMJ.df.hn.strat.bin)
plot(RMJ.df.hn.strat.bin)

# stat + obs
RMJ.df.hn.strat.obs.bin <- ds(data=RMJ.data, region.table=full.region.table, 
                          sample.table=sample.table, obs.table=obs.table, 
                          truncation=trunc.RMJ, key="hn", formula = ~stratum+obs.observer,
                          cutpoints = c(0,7,14,21,28,35,42,49,56,60))
summary(RMJ.df.hn.strat.obs.bin)
plot(RMJ.df.hn.strat.obs.bin)


    # Compare covariate models ####


## unbinned

rmj.df.cov.comp <- summarize_ds_models(RMJ.df.hn.obs,RMJ.df.hn.strat,
                                          output = "plain")

rmj.df.cov.comp[ ,1:5]
rmj.df.cov.comp[ ,6:7]
# observer model has no support (dAIC > 27). But both models have CvM p value < 0.05 so unlikely to be better than the uniform model

# RMJ.df.hn.strat selected 


### binned

rmj.df.cov.comp.bin <- summarize_ds_models(RMJ.df.hn.obs.bin,RMJ.df.hn.strat.bin,RMJ.df.hn.strat.obs.bin,
                                       output = "plain")

rmj.df.cov.comp.bin[ ,1:5]
rmj.df.cov.comp.bin[ ,6:7]
# stratum alone has no support. observer has most support, but model with both also does (just)

# check estimates from the two with support
rmj.cov.bin.comp <- data.frame(model = rep(c("obs","strat+obs"), each=7),
                                 year = c(RMJ.df.hn.obs.bin$dht$individuals$N$Label[1:7],
                                          RMJ.df.hn.strat.obs.bin$dht$individuals$N$Label[1:7]),
                                 estimate = c(RMJ.df.hn.obs.bin$dht$individuals$N$Estimate[1:7],
                                              RMJ.df.hn.strat.obs.bin$dht$individuals$N$Estimate[1:7]),
                                 cv = c(RMJ.df.hn.obs.bin$dht$individuals$N$cv[1:7],
                                        RMJ.df.hn.strat.obs.bin$dht$individuals$N$cv[1:7]),
                                 se = c(RMJ.df.hn.obs.bin$dht$individuals$N$se[1:7],
                                        RMJ.df.hn.strat.obs.bin$dht$individuals$N$se[1:7]),
                                 lcl = c(RMJ.df.hn.obs.bin$dht$individuals$N$lcl[1:7],
                                         RMJ.df.hn.strat.obs.bin$dht$individuals$N$lcl[1:7]),
                                 ucl = c(RMJ.df.hn.obs.bin$dht$individuals$N$ucl[1:7],
                                         RMJ.df.hn.strat.obs.bin$dht$individuals$N$ucl[1:7]))
rmj.cov.bin.comp

# plot cv
ggplot(rmj.cov.bin.comp, aes(x=year, y=cv, colour=model))+
  geom_point()
# variation in which model is better

# plot estimates
ggplot(rmj.cov.bin.comp, aes(x=year, y=estimate, colour=model))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))

# Something dodgy as shit going on with 2016 - wild cv, se, and estimates from both models. Neither selected.


    # Compare primary and covariate models ####


### unbinned

rmj.df.comp <- summarize_ds_models(RMJ.df.unif.cos, RMJ.df.hn.strat, 
                                      output = "plain")
rmj.df.comp[ ,1:3]
rmj.df.comp[ ,4:7]

# covar model has no support (large dAIC and small CvM p value). Uniform model selected


### binned


# I am not going to bother to compare primary and covar models - the covar models are wild (see above)

  ## RMJ Pooled Final Results ####


### unbinned

summary(RMJ.df.unif.cos)

# extract estimates
RMJ.grp.results <- RMJ.df.unif.cos$dht$clusters$N[1:7, ]
RMJ.ind.results <- RMJ.df.unif.cos$dht$individuals$N[1:7, ]

# bind results
RMJ.results <- rbind(RMJ.grp.results, RMJ.ind.results)
RMJ.results <- RMJ.results %>% rename(Year = Label) %>% 
  mutate(Label = rep(c("Grp", "Ind"), each=7)) %>% 
  mutate(Species = rep("RMJ", times=14)) %>% 
  mutate(DetFun = rep("pooled", times=14)) %>%
  mutate(Key = rep("Uni", times=14)) %>% 
  mutate(Adjust = rep("Cos", times=14)) %>% 
  mutate(Covar = "NA", times=14) %>% 
  select(Year,Species,DetFun,Key,Adjust,Covar,Label,Estimate,se,cv,lcl,ucl)



RMJ_final_plot <- ggplot(RMJ.results[RMJ.results$Label=="Ind",], aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                  ylim(0,5500)+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))

# save results
write.csv(RMJ.results, "Output/Results/RMJ_results_final.csv")
ggsave("Output/Results/Plots/Point_estimates/RMJ_final_plot.png", RMJ_final_plot, 
       dpi=300, width = 20, height = 20, units = "cm")


### binned

summary(RMJ.df.hn.cos.bin)

# extract estimates
RMJ.grp.abund <- RMJ.df.hn.cos.bin$dht$clusters$N[1:7, ]
RMJ.grp.abund <- RMJ.grp.abund %>% rename(Year = Label) %>% mutate(Label="Grp")
RMJ.ind.abund <- RMJ.df.hn.cos.bin$dht$individuals$N[1:7, ]
RMJ.ind.abund <- RMJ.ind.abund %>% rename(Year=Label) %>% mutate(Label="Ind")
RMJ.abund <- rbind(RMJ.grp.abund,RMJ.ind.abund)
RMJ.abund <- RMJ.abund %>% rename(N=Estimate,n_se=se,n_cv=cv,n_lcl=lcl,n_ucl=ucl)
RMJ.abund <- RMJ.abund[,-7]

### extract density
RMJ.grp.density <- RMJ.df.hn.cos.bin$dht$clusters$D[1:7, ]
RMJ.grp.density <- RMJ.grp.density %>% rename(Year=Label) %>% mutate(Label="Grp")
RMJ.ind.density <- RMJ.df.hn.cos.bin$dht$individuals$D[1:7, ]
RMJ.ind.density <- RMJ.ind.density %>% rename(Year=Label) %>% mutate(Label="Ind")
RMJ.density <- rbind(RMJ.grp.density,RMJ.ind.density)
RMJ.density <- RMJ.density %>% rename(D=Estimate,d_se=se,d_cv=cv,d_lcl=lcl,d_ucl=ucl)
RMJ.density <- RMJ.density[,-7]
RMJ.density <- RMJ.density[,-c(1,7)]


# bind results
RMJ.results <- cbind(RMJ.abund, RMJ.density)
RMJ.results <- RMJ.results %>%  
                mutate(Species = rep("RMJ", times=14)) %>% 
                mutate(DetFun = rep("pooled", times=14)) %>% 
                mutate(Key = rep("Hn", times=14)) %>% 
                mutate(Adjust = rep("NA", times=14)) %>% 
                mutate(Covar = rep("NA", times=14)) %>% 
                select(Year,Species,DetFun,Key,Adjust,Covar,Label,
                       N,n_se,n_cv,n_lcl,n_ucl,
                       D,d_se,d_cv,d_lcl,d_ucl)

RMJ.results <- RMJ.results %>% mutate(D=D*1000000,d_se=d_se*1000000,d_lcl=d_lcl*1000000,d_ucl=d_ucl*1000000)




RMJ_final_plot_bin <- ggplot(RMJ.results.bin[RMJ.results.bin$Label=="Ind",], aes(x=Year, y=Estimate))+
                      geom_point(shape=16, size=2)+
                      geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
                      ylim(0,6150)+
                      theme_bw()+
                      theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),panel.border = element_blank())+
                      theme(axis.line = element_line(color = 'black'))

# save results
write.csv(RMJ.results, "Output/Results/RMJ_results_final_binned.csv")
ggsave("Output/Results/Plots/Point_estimates/RMJ_final_plot_binned.png", RMJ_final_plot_bin, 
       dpi=300, width = 20, height = 20, units = "cm")



# plot comparison between binned and unbinned

# load in unbinned results
rmj.results.orig <- read.csv("Output/Results/RMJ_results_final.csv")
rmj.results.orig <- rmj.results.orig[ ,-1]
rmj.results.orig$analysis <- "orig"

rmj.results.bin <- read.csv("Output/Results/RMJ_results_final_binned.csv")
rmj.results.bin <- rmj.results.bin[,-1]
rmj.results.bin$analysis <- "binned"


RMJ.comp <- rbind(rmj.results.orig,rmj.results.bin)

RMJ_analysis_comparison_plot <- ggplot(
                              RMJ.comp[RMJ.comp$Label=="Ind",], aes(x=Year, y=Estimate, colour=analysis))+
                              geom_point(position = position_dodge(width=0.3))+
                              geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))+
                              theme_bw()+
                              theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),panel.border = element_blank())+
                              theme(axis.line = element_line(color = 'black'))
ggsave("Output/Results/Plots/Point_estimates/RMJ_analysis_comparison_plot.png", RMJ_analysis_comparison_plot, 
       dpi=300, width = 20, height = 20, units = "cm")



##### ALL SPECIES FINAL RESULTS ---------------------------------------------------------


## New results (Feb 21) include individual & group abundance and density for all species

# load all species results
BSD.results <- read.csv("Output/Results/BSD_results_final.csv")
YCG.results <- read.csv("Output/Results/YCG_results_final.csv")
GSL.results <- read.csv("Output/Results/GSL_results_final_binned.csv")
LTM.results <- read.csv("Output/Results/LTM_results_final_binned.csv")
PTM.results <- read.csv("Output/Results/PTM_results_final.csv")
STM.results <- read.csv("Output/Results/STM_results_final.csv")
BTG.results <- read.csv("Output/Results/BTG_results_final.csv")
GAU.results <- read.csv("Output/Results/GAU_results_final.csv")
PIG.results <- read.csv("Output/Results/PIG_results_final_binned.csv")
RMJ.results <- read.csv("Output/Results/RMJ_results_final_binned.csv")
GPF.results <- read.csv("Output/Results/GPF_results_final.csv")

# merge
Results_AllSpecies_2010_2020 <- rbind(BSD.results,YCG.results,GSL.results,LTM.results,PTM.results,
                                      STM.results,BTG.results,GAU.results,PIG.results,RMJ.results,
                                      GPF.results)
str(Results_AllSpecies_2010_2020)
Results_AllSpecies_2010_2020 <- Results_AllSpecies_2010_2020[,-1]

Results_AllSpecies_2010_2020$Year <- as.factor(Results_AllSpecies_2010_2020$Year)

# save
write.csv(Results_AllSpecies_2010_2020, file="Output/Results/Results_AllSpecies_2010_2020.csv")


### plotting all

# split prims & ungs
prim <- c("BSD","YCG","GSL","LTM","PTM","STM") 
prims <- Results_AllSpecies_2010_2020[Results_AllSpecies_2010_2020$Species %in% prim, ]
prims <- prims %>% filter(Label=="Grp")

ung <- c("BTG","GAU","PIG","RMJ","GPF")
ungs <- Results_AllSpecies_2010_2020[Results_AllSpecies_2010_2020$Species %in% ung, ]
ungs <- ungs %>% filter(Label=="Ind")


# plot
all_prims_plot <- ggplot(prims, aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
                  scale_y_continuous(limits=c(0,NA))+
                  facet_wrap(~prims$Species, scales = "free")+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))+
                  ylab("Estimated abundance of groups")

all_ungs_plot <- ggplot(ungs, aes(x=Year, y=Estimate))+
                  geom_point(shape=16, size=2)+
                  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
                  scale_y_continuous(limits=c(0,NA))+
                  facet_wrap(~ungs$Species, scales = "free")+
                  theme_bw()+
                  theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),panel.border = element_blank())+
                  theme(axis.line = element_line(color = 'black'))+
                  ylab("Estimated abundance of individuals")
  
ggsave("Output/Results/Plots/Point_estimates/all_prims_plot.png", all_prims_plot, 
       dpi=300, width = 20, height = 20, units = "cm")
ggsave("Output/Results/Plots/Point_estimates/all_ungs_plot.png", all_ungs_plot, 
       dpi=300, width = 20, height = 20, units = "cm")



### Density

# bind all density results together
all_spp_density <- rbind(BSD.density,YCG.grp.density,GSL.grp.density.bin,LTM.grp.density.bin,
                         PTM.grp.density,STM.grp.density,BTG.ind.density,GAU.ind.density,PIG.ind.density.bin,
                         GPF.ind.density,RMJ.ind.density.bin)

all_spp_density <- all_spp_density %>% arrange(Species,Year)

write.csv(all_spp_density, file="Output/Results/ALL_SPP_DENSITY.csv")
