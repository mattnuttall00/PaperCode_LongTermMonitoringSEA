### August 2020. This script is a comparison of different analytical approaches for the analysis of line transect distance sampling data from Keo Seima Wildlife Sanctuary between 2010 - 2020.  The purpose is to assess the differences in annual abundance estimates (and associated variances) produced when a) the geographical priority strata, and b) year as detection function covariate, are inlcuded or excluded. This script was written by Matt Nuttall (m.n.nuttall1@stir.ac.uk / mattnuttall00@gmail.com). The data belong to WCS. This was done in collaboration with WCS (Olly Griffin - ogriffin@wcs.org) 

### Background. The line transect surveys in Keo Seima have taken two different approaches over the years. In 2010, 2013, 2018, and 2020 effort across transects was approximately uniform (barring other aspects such as logistics). In 2011, 2014, and 2016 there was a different approach which allocated significantly higher effort towards "high priority" transects. These were the transects where anmial density is known to be higher, and where detections were higher. Because of the different approaches in different years, there were some analytical challenges and questions around whether or not the priority strata years needed to be stratified during analysis (because of differing levels of effort).

# Matt Nuttall had a conversation with Eric Rextad from St Andrews University about this, and his overall conclusion was that stratification probably wasn't necessary, but that in the end, the best way to decide would be to try both approaches and see what difference it made to the estimates. That is what was done below.

#### Load libraries and data ####

library(Distance)
library(dsm)
library(DSsim)
library(mads)
library(rmarkdown)
library(mrds)
library(knitr)
library(tidyverse)
library(patchwork)
library(devtools)
library(remotes)
install_github("DistanceDevelopment/mrds")
install_github("DistanceDevelopment/Distance")

load("./Output/Data/KSWS_MASTER.Rdata")


#### Yellow-cheeked crested gibbon ###########################
### test analysis with NO priority strata and POOLED detection function across years ####

# This analysis removes the high and low priority strata, and uses a single detection function for all years.

## remove strata

# all data
data.noStrata <- allData
data.noStrata$stratum <- as.character(data.noStrata$stratum)
data.noStrata <- data.noStrata %>% mutate(stratum = ifelse(stratum=="2011_H", "2011", stratum))
data.noStrata <- data.noStrata %>% mutate(stratum = ifelse(stratum=="2011_L", "2011", stratum))
data.noStrata <- data.noStrata %>% mutate(stratum = ifelse(stratum=="2014_H", "2014", stratum))
data.noStrata <- data.noStrata %>% mutate(stratum = ifelse(stratum=="2014_L", "2014", stratum))
data.noStrata <- data.noStrata %>% mutate(stratum = ifelse(stratum=="2016_H", "2016", stratum))
data.noStrata <- data.noStrata %>% mutate(stratum = ifelse(stratum=="2016_L", "2016", stratum))

data.noStrata$stratum <- as.factor(data.noStrata$stratum)
str(allData)

# sample table
sample.table.nostrata <- sample.table
sample.table.nostrata$Region.Label <- as.character(sample.table.nostrata$Region.Label)

sample.table.nostrata <- sample.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2011_H", "2011", Region.Label))
sample.table.nostrata <- sample.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2011_L", "2011", Region.Label))
sample.table.nostrata <- sample.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2014_H", "2014", Region.Label))
sample.table.nostrata <- sample.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2014_L", "2014", Region.Label))
sample.table.nostrata <- sample.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2016_H", "2016", Region.Label))
sample.table.nostrata <- sample.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2016_L", "2016", Region.Label))

# regional table
region.table.nostrata <- data.frame(Region.Label = c("2010","2011","2013","2014","2016","2018"),
                                    Area = c(full.region.table[1,2],
                                             sum(full.region.table[2,2],full.region.table[3,2]),
                                             full.region.table[4,2],
                                             sum(full.region.table[5,2],full.region.table[6,2]),
                                             sum(full.region.table[7,2],full.region.table[8,2]),
                                             full.region.table[9,2]))
# obs table
obs.table.nostrata <- obs.table
obs.table.nostrata$Region.Label <- as.character(obs.table.nostrata$Region.Label)

obs.table.nostrata <- obs.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2011_H", "2011", Region.Label))
obs.table.nostrata <- obs.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2011_L", "2011", Region.Label))
obs.table.nostrata <- obs.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2014_H", "2014", Region.Label))
obs.table.nostrata <- obs.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2014_L", "2014", Region.Label))
obs.table.nostrata <- obs.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2016_H", "2016", Region.Label))
obs.table.nostrata <- obs.table.nostrata %>% 
                          mutate(Region.Label = ifelse(Region.Label=="2016_L", "2016", Region.Label))



## subset data for YCG
YCG.data.nostrata <-data.noStrata[data.noStrata$species=="YCG",]

# plot histogram
hist(YCG.data.nostrata$distance)


## fit detection function

# uniform
ycg.nostrat.df.unif.co <- ds(data=YCG.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=60, key="unif", adjustment= "cos")

ycg.nostrat.df.unif.poly <- ds(data=YCG.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=60, key="unif", adjustment= "poly")

# half normal
ycg.nostrat.df.hn.co <- ds(data=YCG.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=60, key="hn", adjustment= "cos")

ycg.nostrat.df.hn.herm <- ds(data=YCG.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=60, key="hn", adjustment= "herm")

# hazard rate
ycg.nostrat.df.hr.co <- ds(data=YCG.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=60, key="hr", adjustment= "cos")

ycg.nostrat.df.hr.poly <- ds(data=YCG.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=60, key="hr", adjustment= "poly")

# compare models
ycg.nostrat.comp <- summarize_ds_models(ycg.nostrat.df.unif.co,ycg.nostrat.df.unif.poly,
                                        ycg.nostrat.df.hn.co,ycg.nostrat.df.hn.herm,
                                        ycg.nostrat.df.hr.co,ycg.nostrat.df.hr.poly,
                                             output = "plain")
ycg.nostrat.comp[ ,1:5]
ycg.nostrat.comp[ ,6:7]

par(mfrow=c(3,2))
plot(ycg.nostrat.df.unif.co)
plot(ycg.nostrat.df.unif.poly)
plot(ycg.nostrat.df.hn.co)
plot(ycg.nostrat.df.hn.herm)
plot(ycg.nostrat.df.hr.co)
plot(ycg.nostrat.df.hr.poly)
# hn co selected

summary(ycg.nostrat.df.hn.co)

## extract estimates, CIs, CVs

# 2010
YCG.nostrat.pooled.ind.abund.10 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2010"]

YCG.nostrat.pooled.LCL.10 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2010"]

YCG.nostrat.pooled.UCL.10 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2010"]

YCG.nostrat.pooled.CV.10 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$cv[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2010"]


# 2011
YCG.nostrat.pooled.ind.abund.11 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2011"]

YCG.nostrat.pooled.LCL.11 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2011"]

YCG.nostrat.pooled.UCL.11 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2011"]

YCG.nostrat.pooled.CV.11 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$cv[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2011"]

# 2013
YCG.nostrat.pooled.ind.abund.13 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2013"]

YCG.nostrat.pooled.LCL.13 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2013"]

YCG.nostrat.pooled.UCL.13 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2013"]

YCG.nostrat.pooled.CV.13 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$cv[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2013"]

# 2014
YCG.nostrat.pooled.ind.abund.14 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2014"]

YCG.nostrat.pooled.LCL.14 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2014"]

YCG.nostrat.pooled.UCL.14 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2014"]

YCG.nostrat.pooled.CV.14 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$cv[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2014"]

# 2016
YCG.nostrat.pooled.ind.abund.16 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2016"]

YCG.nostrat.pooled.LCL.16 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2016"]

YCG.nostrat.pooled.UCL.16 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2016"]

YCG.nostrat.pooled.CV.16 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$cv[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2016"]

# 2018
YCG.nostrat.pooled.ind.abund.18 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2018"]

YCG.nostrat.pooled.LCL.18 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2018"]

YCG.nostrat.pooled.UCL.18 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2018"]

YCG.nostrat.pooled.CV.18 <- 
  ycg.nostrat.df.hn.co$dht$individuals$N$cv[
   ycg.nostrat.df.hn.co$dht$individuals$N$Label == "2018"]


ycg_plot_dat_pool <- data.frame(DetFun = rep("pool", times=6),
                           year = c("2010","2011","2013","2014","2016","2018"),
                           abundance = c(YCG.nostrat.pooled.ind.abund.10,YCG.nostrat.pooled.ind.abund.11,
                                         YCG.nostrat.pooled.ind.abund.13,YCG.nostrat.pooled.ind.abund.14,
                                         YCG.nostrat.pooled.ind.abund.16,YCG.nostrat.pooled.ind.abund.18),
                           lcl = c(YCG.nostrat.pooled.LCL.10,YCG.nostrat.pooled.LCL.11,
                                   YCG.nostrat.pooled.LCL.13,YCG.nostrat.pooled.LCL.14,
                                   YCG.nostrat.pooled.LCL.16,YCG.nostrat.pooled.LCL.18),
                           ucl = c(YCG.nostrat.pooled.UCL.10,YCG.nostrat.pooled.UCL.11,
                                   YCG.nostrat.pooled.UCL.13,YCG.nostrat.pooled.UCL.14,
                                   YCG.nostrat.pooled.UCL.16,YCG.nostrat.pooled.UCL.18),
                           cv = c(YCG.nostrat.pooled.CV.10,YCG.nostrat.pooled.CV.11,
                                  YCG.nostrat.pooled.CV.13,YCG.nostrat.pooled.CV.14,
                                  YCG.nostrat.pooled.CV.16,YCG.nostrat.pooled.CV.18))



ggplot(ycg_plot_dat_pool, aes(x=year, y=abundance))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3)+
  ylim(0,4000)

### test analysis with no priority strata but year as a covariate in detection function ####


# half normal
ycg.nostrat.df.hn.yr <- ds(data=YCG.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=60, key="hn", formula = ~stratum)

# compare model with year as covariate and model with pooled DF
ycg.nostrat.comp.pool <- summarize_ds_models(ycg.nostrat.df.hn.co,ycg.nostrat.df.hn.yr,
                                             output = "plain")
ycg.nostrat.comp.pool[ ,1:5]
ycg.nostrat.comp.pool[ ,6:7]
# virtually no differece in AIC

# extract estimats from model with stratum as covariate

summary(ycg.nostrat.df.hn.yr)

# 2010
YCG.nostrat.yr.ind.abund.10 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2010"]

YCG.nostrat.yr.LCL.10 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2010"]

YCG.nostrat.yr.UCL.10 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2010"]

YCG.nostrat.yr.cv.10 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2010"]


# 2011
YCG.nostrat.yr.ind.abund.11 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2011"]

YCG.nostrat.yr.LCL.11 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2011"]

YCG.nostrat.yr.UCL.11 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2011"]

YCG.nostrat.yr.cv.11 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2011"]

# 2013
YCG.nostrat.yr.ind.abund.13 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2013"]

YCG.nostrat.yr.LCL.13 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2013"]

YCG.nostrat.yr.UCL.13 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2013"]

YCG.nostrat.yr.cv.13 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2013"]

# 2014
YCG.nostrat.yr.ind.abund.14 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2014"]

YCG.nostrat.yr.LCL.14 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2014"]

YCG.nostrat.yr.UCL.14 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2014"]

YCG.nostrat.yr.cv.14 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2014"]

# 2016
YCG.nostrat.yr.ind.abund.16 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2016"]

YCG.nostrat.yr.LCL.16 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2016"]

YCG.nostrat.yr.UCL.16 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2016"]

YCG.nostrat.yr.cv.16 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2016"]

# 2018
YCG.nostrat.yr.ind.abund.18 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2018"]

YCG.nostrat.yr.LCL.18 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2018"]

YCG.nostrat.yr.UCL.18 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2018"]

YCG.nostrat.yr.cv.18 <- 
  ycg.nostrat.df.hn.yr$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr$dht$individuals$N$Label == "2018"]


ycg_plot_dat_yr <- data.frame(DetFun = rep("yr", times=6),
                           year = c("2010","2011","2013","2014","2016","2018"),
                           abundance = c(YCG.nostrat.yr.ind.abund.10,YCG.nostrat.yr.ind.abund.11,
                                         YCG.nostrat.yr.ind.abund.13,YCG.nostrat.yr.ind.abund.14,
                                         YCG.nostrat.yr.ind.abund.16,YCG.nostrat.yr.ind.abund.18),
                           lcl = c(YCG.nostrat.yr.LCL.10,YCG.nostrat.yr.LCL.11,
                                   YCG.nostrat.yr.LCL.13,YCG.nostrat.yr.LCL.14,
                                   YCG.nostrat.yr.LCL.16,YCG.nostrat.yr.LCL.18),
                           ucl = c(YCG.nostrat.yr.UCL.10,YCG.nostrat.yr.UCL.11,
                                   YCG.nostrat.yr.UCL.13,YCG.nostrat.yr.UCL.14,
                                   YCG.nostrat.yr.UCL.16,YCG.nostrat.yr.UCL.18),
                           cv = c(YCG.nostrat.yr.cv.10,YCG.nostrat.yr.cv.11,
                                  YCG.nostrat.yr.cv.13,YCG.nostrat.yr.cv.14,
                                  YCG.nostrat.yr.cv.16,YCG.nostrat.yr.cv.18))

# merge the pooled and strtified dataframes
ycg_plot_dat_merge <- rbind(ycg_plot_dat_pool,ycg_plot_dat_yr)

# plot together
ggplot(ycg_plot_dat_merge, aes(x=year, y=abundance, group = DetFun, colour=DetFun))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3, position = position_dodge(width = 0.5))
 
 


### test analysis with priority strata in ALL years and pooled DF ####


# Here I am going to create priority strata for all years.  For 2011, 2014, and 2016 I will use the correct strata (i.e. the ones that were done).  For 2010 and 2013 I will use the strata from 2011. For 2018 I will use the strata from 2016. 

  ## create data with strata in ALL years ####

# need to change strata in allData, sample.table, region.table, and obs.table

    # allData ####

data_allstrata <- allData
str(data_allstrata)

# which transects are high / low priority in 2011?
H_pri_2011 <- unique(data_allstrata$obs.transect[data_allstrata$stratum=="2011_H"])
L_pri_2011 <- unique(data_allstrata$obs.transect[data_allstrata$stratum=="2011_L"])
sum(length(H_pri_2011), length(L_pri_2011)) # transect 23 missing
allData %>% filter(stratum=="2011_H"|stratum=="2011_L") %>% filter(obs.transect==23)
# allocate T23 to L strata in 2010 - very few obs

## change stratum

# 2010
data_allstrata$stratum <- as.character(data_allstrata$stratum)
data_allstrata <- data_allstrata %>% 
                  mutate(stratum = ifelse(
                    stratum=="2010" & obs.transect %in% H_pri_2011,
                    replace(stratum, stratum==stratum,"2010_H"),
                    stratum))
unique(data_allstrata$stratum)      


data_allstrata <- data_allstrata %>% 
                  mutate(stratum = ifelse(
                    stratum=="2010" & obs.transect %in% L_pri_2011,
                    replace(stratum, stratum==stratum,"2010_L"),
                    stratum))
unique(data_allstrata$stratum)      
data_allstrata %>% filter(stratum==2010)

data_allstrata <- data_allstrata %>% 
                  mutate(stratum = ifelse(
                    stratum=="2010" & obs.transect==23,
                    replace(stratum, stratum==stratum,"2010_L"),
                    stratum))
unique(data_allstrata$stratum)      

# 2013
data_allstrata <- data_allstrata %>% 
                  mutate(stratum = ifelse(
                    stratum=="2013" & obs.transect %in% H_pri_2011,
                    replace(stratum, stratum==stratum,"2013_H"),
                    stratum))

data_allstrata <- data_allstrata %>% 
                  mutate(stratum = ifelse(
                    stratum=="2013" & obs.transect %in% L_pri_2011,
                    replace(stratum, stratum==stratum,"2013_L"),
                    stratum))

data_allstrata <- data_allstrata %>% 
                  mutate(stratum = ifelse(
                    stratum=="2013" & obs.transect==23,
                    replace(stratum, stratum==stratum,"2013_L"),
                    stratum))

unique(data_allstrata$stratum)  
data_allstrata %>% filter(stratum==2013)

# 2018
# which transects are high / low priority in 2016?
H_pri_2016 <- unique(data_allstrata$obs.transect[data_allstrata$stratum=="2016_H"])
L_pri_2016 <- unique(data_allstrata$obs.transect[data_allstrata$stratum=="2016_L"])
sum(length(H_pri_2016), length(L_pri_2016)) # T20 and T30 missing
# T20 doesn't exist in 2018 either so thats fine, and only very few obs on T30 in 2018, so that will be allocated to low priority

data_allstrata <- data_allstrata %>% 
                  mutate(stratum = ifelse(
                    stratum=="2018" & obs.transect %in% H_pri_2016,
                    replace(stratum, stratum==stratum,"2018_H"),
                    stratum))

data_allstrata <- data_allstrata %>% 
                  mutate(stratum = ifelse(
                    stratum=="2018" & obs.transect %in% L_pri_2016,
                    replace(stratum, stratum==stratum,"2018_L"),
                    stratum))

data_allstrata <- data_allstrata %>% 
                  mutate(stratum = ifelse(
                    stratum=="2018" & obs.transect==30,
                    replace(stratum, stratum==stratum,"2018_L"),
                    stratum))

unique(data_allstrata$stratum) 

data_allstrata$stratum <- as.factor(data_allstrata$stratum)
levels(data_allstrata$stratum)

# subset data for UCG
YCG.data.allstrata <- data_allstrata[data_allstrata$species=="YCG", ]

    # sample.table #####

str(sample.table)
sample.table.allstrata <- sample.table
sample.table.allstrata$Region.Label <- as.character(sample.table.allstrata$Region.Label)
str(sample.table.allstrata)

# 2010
sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2010" & Sample.Label %in% H_pri_2011,
                            replace(Region.Label, Region.Label==Region.Label,"2010_H"),
                            Region.Label))

sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2010" & Sample.Label %in% L_pri_2011,
                            replace(Region.Label, Region.Label==Region.Label,"2010_L"),
                            Region.Label))

sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2010" & Sample.Label==23,
                            replace(Region.Label, Region.Label==Region.Label,"2010_L"),
                            Region.Label))

unique(sample.table.allstrata$Region.Label)  

# 2013
sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2013" & Sample.Label %in% H_pri_2011,
                            replace(Region.Label, Region.Label==Region.Label,"2013_H"),
                            Region.Label))

sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2013" & Sample.Label %in% L_pri_2011,
                            replace(Region.Label, Region.Label==Region.Label,"2013_L"),
                            Region.Label))

sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2013" & Sample.Label==23,
                            replace(Region.Label, Region.Label==Region.Label,"2013_L"),
                            Region.Label))

# 2018
sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2018" & Sample.Label %in% H_pri_2016,
                            replace(Region.Label, Region.Label==Region.Label,"2018_H"),
                            Region.Label))

sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2018" & Sample.Label %in% L_pri_2016,
                            replace(Region.Label, Region.Label==Region.Label,"2018_L"),
                            Region.Label))

sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2018" & Sample.Label==30,
                            replace(Region.Label, Region.Label==Region.Label,"2018_L"),
                            Region.Label))

sample.table.allstrata <- sample.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2018" & Sample.Label==20,
                            replace(Region.Label, Region.Label==Region.Label,"2018_L"),
                            Region.Label))

unique(sample.table.allstrata$Region.Label)  
sample.table.allstrata %>% filter(Region.Label=="2018")

    # Region.Table ####

str(full.region.table)

region.table.allstrata <- data.frame(Region.Label = c("2010_H","2010_L","2011_H","2011_L",
                                                      "2013_H","2013_L","2014_H","2014_L",
                                                      "2016_H","2016_L","2018_H","2018_L"),
                                     Area = c(903000000,903000000,903000000,903000000,903000000,903000000,
                                              903000000,903000000,858000000,949000000,858000000,949000000))

    # obs.table ####

str(obs.table)

obs.table.allstrata <- obs.table
obs.table.allstrata$Region.Label <- as.character(obs.table.allstrata$Region.Label)

# 2010
obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2010" & Sample.Label %in% H_pri_2011,
                            replace(Region.Label, Region.Label==Region.Label,"2010_H"),
                            Region.Label))

obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2010" & Sample.Label %in% L_pri_2011,
                            replace(Region.Label, Region.Label==Region.Label,"2010_L"),
                            Region.Label))

obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2010" & Sample.Label==23,
                            replace(Region.Label, Region.Label==Region.Label,"2010_L"),
                            Region.Label))

unique(obs.table.allstrata$Region.Label)

# 2013
obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2013" & Sample.Label %in% H_pri_2011,
                            replace(Region.Label, Region.Label==Region.Label,"2013_H"),
                            Region.Label))

obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2013" & Sample.Label %in% L_pri_2011,
                            replace(Region.Label, Region.Label==Region.Label,"2013_L"),
                            Region.Label))

obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2013" & Sample.Label==23,
                            replace(Region.Label, Region.Label==Region.Label,"2013_L"),
                            Region.Label))

# 2018
obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2018" & Sample.Label %in% H_pri_2016,
                            replace(Region.Label, Region.Label==Region.Label,"2018_H"),
                            Region.Label))

obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2018" & Sample.Label %in% L_pri_2016,
                            replace(Region.Label, Region.Label==Region.Label,"2018_L"),
                            Region.Label))

obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2018" & Sample.Label==30,
                            replace(Region.Label, Region.Label==Region.Label,"2018_L"),
                            Region.Label))

obs.table.allstrata <- obs.table.allstrata %>% 
                          mutate(Region.Label = ifelse(
                            Region.Label=="2018" & Sample.Label==20,
                            replace(Region.Label, Region.Label==Region.Label,"2018_L"),
                            Region.Label))

unique(obs.table.allstrata$Region.Label)

  ## fit detection function and get estimates ####
ycg.allstrat.df.hn.co <- ds(data=YCG.data.allstrata, region.table=region.table.allstrata, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", adjustment= "cos")

summary(ycg.allstrat.df.hn.co)

# the SE and CV's are not just the sums or means of the two strata, and so I have to go through to proces of setting all other strata areas to 0 for each year.

# set regions to 0 for 2011-2018
region.table.YCG.allstrata.pool.10 <- region.table.allstrata
region.table.YCG.allstrata.pool.10$Area[region.table.YCG.allstrata.pool.10$Region.Label != "2010_L" & 
                                          region.table.YCG.allstrata.pool.10$Region.Label != "2010_H"] <- 0

# re-run model
ycg.allstrat.df.hn.co.10 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.10, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", adjustment= "cos")
summary(ycg.allstrat.df.hn.co.10)

# set regions to 0 for all except 2011
region.table.YCG.allstrata.pool.11 <- region.table.allstrata
region.table.YCG.allstrata.pool.11$Area[region.table.YCG.allstrata.pool.11$Region.Label != "2011_L" & 
                                          region.table.YCG.allstrata.pool.11$Region.Label != "2011_H"] <- 0

# re-run model
ycg.allstrat.df.hn.co.11 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.11, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", adjustment= "cos")
summary(ycg.allstrat.df.hn.co.11)

# set regions to 0 for all except 2013
region.table.YCG.allstrata.pool.13 <- region.table.allstrata
region.table.YCG.allstrata.pool.13$Area[region.table.YCG.allstrata.pool.13$Region.Label != "2013_L" & 
                                          region.table.YCG.allstrata.pool.13$Region.Label != "2013_H"] <- 0

# re-run model
ycg.allstrat.df.hn.co.13 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.13, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", adjustment= "cos")
summary(ycg.allstrat.df.hn.co.13)


# set regions to 0 for all except 2014
region.table.YCG.allstrata.pool.14 <- region.table.allstrata
region.table.YCG.allstrata.pool.14$Area[region.table.YCG.allstrata.pool.14$Region.Label != "2014_L" & 
                                          region.table.YCG.allstrata.pool.14$Region.Label != "2014_H"] <- 0


# re-run model
ycg.allstrat.df.hn.co.14 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.14, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", adjustment= "cos")
summary(ycg.allstrat.df.hn.co.14)


# set regions to 0 for all except 2016
region.table.YCG.allstrata.pool.16 <- region.table.allstrata
region.table.YCG.allstrata.pool.16$Area[region.table.YCG.allstrata.pool.16$Region.Label != "2016_L" & 
                                          region.table.YCG.allstrata.pool.16$Region.Label != "2016_H"] <- 0


# re-run model
ycg.allstrat.df.hn.co.16 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.16, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", adjustment= "cos")
summary(ycg.allstrat.df.hn.co.16)


# set regions to 0 for all except 2018
region.table.YCG.allstrata.pool.18 <- region.table.allstrata
region.table.YCG.allstrata.pool.18$Area[region.table.YCG.allstrata.pool.18$Region.Label != "2018_L" & 
                                          region.table.YCG.allstrata.pool.18$Region.Label != "2018_H"] <- 0


# re-run model
ycg.allstrat.df.hn.co.18 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.18, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", adjustment= "cos")
summary(ycg.allstrat.df.hn.co.18)


# pull out estimates

# 2010
YCG.allstrat.pool.ind.abund.10 <- 
  ycg.allstrat.df.hn.co.10$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.co.10$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.LCL.10 <- 
  ycg.allstrat.df.hn.co.10$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.co.10$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.UCL.10 <- 
  ycg.allstrat.df.hn.co.10$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.co.10$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.cv.10 <- 
  ycg.allstrat.df.hn.co.10$dht$individuals$N$cv[
   ycg.allstrat.df.hn.co.10$dht$individuals$N$Label == "Total"]

# 2011
YCG.allstrat.pool.ind.abund.11 <- 
  ycg.allstrat.df.hn.co.11$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.co.11$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.LCL.11 <- 
  ycg.allstrat.df.hn.co.11$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.co.11$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.UCL.11 <- 
  ycg.allstrat.df.hn.co.11$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.co.11$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.cv.11 <- 
  ycg.allstrat.df.hn.co.11$dht$individuals$N$cv[
   ycg.allstrat.df.hn.co.11$dht$individuals$N$Label == "Total"]

# 2013
YCG.allstrat.pool.ind.abund.13 <- 
  ycg.allstrat.df.hn.co.13$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.co.13$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.LCL.13 <- 
  ycg.allstrat.df.hn.co.13$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.co.13$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.UCL.13 <- 
  ycg.allstrat.df.hn.co.13$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.co.13$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.cv.13 <- 
  ycg.allstrat.df.hn.co.13$dht$individuals$N$cv[
   ycg.allstrat.df.hn.co.13$dht$individuals$N$Label == "Total"]

# 2014
YCG.allstrat.pool.ind.abund.14 <- 
  ycg.allstrat.df.hn.co.14$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.co.14$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.LCL.14 <- 
  ycg.allstrat.df.hn.co.14$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.co.14$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.UCL.14 <- 
  ycg.allstrat.df.hn.co.14$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.co.14$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.cv.14 <- 
  ycg.allstrat.df.hn.co.14$dht$individuals$N$cv[
   ycg.allstrat.df.hn.co.14$dht$individuals$N$Label == "Total"]

# 2016
YCG.allstrat.pool.ind.abund.16 <- 
  ycg.allstrat.df.hn.co.16$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.co.16$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.LCL.16 <- 
  ycg.allstrat.df.hn.co.16$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.co.16$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.UCL.16 <- 
  ycg.allstrat.df.hn.co.16$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.co.16$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.cv.16 <- 
  ycg.allstrat.df.hn.co.16$dht$individuals$N$cv[
   ycg.allstrat.df.hn.co.16$dht$individuals$N$Label == "Total"]

# 2018
YCG.allstrat.pool.ind.abund.18 <- 
  ycg.allstrat.df.hn.co.18$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.co.18$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.LCL.18 <- 
  ycg.allstrat.df.hn.co.18$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.co.18$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.UCL.18 <- 
  ycg.allstrat.df.hn.co.18$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.co.18$dht$individuals$N$Label == "Total"]

YCG.allstrat.pool.cv.18 <- 
  ycg.allstrat.df.hn.co.18$dht$individuals$N$cv[
   ycg.allstrat.df.hn.co.18$dht$individuals$N$Label == "Total"]


# dataframe
ycg_plot_dat_allstrat_pool <- data.frame(DetFun = rep("stratPool", times=6),
                              year = c("2010","2011","2013","2014","2016","2018"),
                              abundance = c(YCG.allstrat.pool.ind.abund.10,YCG.allstrat.pool.ind.abund.11,
                                            YCG.allstrat.pool.ind.abund.13,YCG.allstrat.pool.ind.abund.14,
                                            YCG.allstrat.pool.ind.abund.16,YCG.allstrat.pool.ind.abund.18),
                              lcl = c(YCG.allstrat.pool.LCL.10,YCG.allstrat.pool.LCL.11,
                                      YCG.allstrat.pool.LCL.13,YCG.allstrat.pool.LCL.14,
                                      YCG.allstrat.pool.LCL.16,YCG.allstrat.pool.LCL.18),
                              ucl = c(YCG.allstrat.pool.UCL.10,YCG.allstrat.pool.UCL.11,
                                      YCG.allstrat.pool.UCL.13,YCG.allstrat.pool.UCL.14,
                                      YCG.allstrat.pool.UCL.16,YCG.allstrat.pool.UCL.18),
                             cv = c(YCG.allstrat.pool.cv.10,YCG.allstrat.pool.cv.11,
                                    YCG.allstrat.pool.cv.13,YCG.allstrat.pool.cv.14,
                                    YCG.allstrat.pool.cv.16,YCG.allstrat.pool.cv.18))

### test anlaysis with priority strata in ALL years and year as covariate in detection function ####

# fit detection function

ycg.allstrat.df.hn.yr <- ds(data=YCG.data.allstrata, region.table=region.table.allstrata, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", formula = ~stratum)

summary(ycg.allstrat.df.hn.yr)


# re-run model for 2010 estimates
ycg.allstrat.df.hn.yr.10 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.10, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", formula = ~stratum)
summary(ycg.allstrat.df.hn.yr.10)


# re-run model for 2011 estimates
ycg.allstrat.df.hn.yr.11 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.11, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", formula = ~stratum)
summary(ycg.allstrat.df.hn.yr.11)

# re-run model for 2013 estimates
ycg.allstrat.df.hn.yr.13 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.13, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", formula = ~stratum)
summary(ycg.allstrat.df.hn.yr.13)

# re-run model for 2014 estimates
ycg.allstrat.df.hn.yr.14 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.14, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", formula = ~stratum)
summary(ycg.allstrat.df.hn.yr.14)

# re-run model for 2016 estimates
ycg.allstrat.df.hn.yr.16 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.16, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", formula = ~stratum)
summary(ycg.allstrat.df.hn.yr.16)

# re-run model for 2018 estimates
ycg.allstrat.df.hn.yr.18 <- ds(data=YCG.data.allstrata, region.table=region.table.YCG.allstrata.pool.18, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=60, key="hn", formula = ~stratum)
summary(ycg.allstrat.df.hn.yr.18)


# pull out estimates

# 2010
YCG.allstrat.yr.ind.abund.10 <- 
  ycg.allstrat.df.hn.yr.10$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.yr.10$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.LCL.10 <- 
  ycg.allstrat.df.hn.yr.10$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.yr.10$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.UCL.10 <- 
  ycg.allstrat.df.hn.yr.10$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.yr.10$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.cv.10 <- 
  ycg.allstrat.df.hn.yr.10$dht$individuals$N$cv[
   ycg.allstrat.df.hn.yr.10$dht$individuals$N$Label == "Total"]

# 2011
YCG.allstrat.yr.ind.abund.11 <- 
  ycg.allstrat.df.hn.yr.11$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.yr.11$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.LCL.11 <- 
  ycg.allstrat.df.hn.yr.11$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.yr.11$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.UCL.11 <- 
  ycg.allstrat.df.hn.yr.11$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.yr.11$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.cv.11 <- 
  ycg.allstrat.df.hn.yr.11$dht$individuals$N$cv[
   ycg.allstrat.df.hn.yr.11$dht$individuals$N$Label == "Total"]

# 2013
YCG.allstrat.yr.ind.abund.13 <- 
  ycg.allstrat.df.hn.yr.13$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.yr.13$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.LCL.13 <- 
  ycg.allstrat.df.hn.yr.13$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.yr.13$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.UCL.13 <- 
  ycg.allstrat.df.hn.yr.13$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.yr.13$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.cv.13 <- 
  ycg.allstrat.df.hn.yr.13$dht$individuals$N$cv[
   ycg.allstrat.df.hn.yr.13$dht$individuals$N$Label == "Total"]

# 2014
YCG.allstrat.yr.ind.abund.14 <- 
  ycg.allstrat.df.hn.yr.14$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.yr.14$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.LCL.14 <- 
  ycg.allstrat.df.hn.yr.14$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.yr.14$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.UCL.14 <- 
  ycg.allstrat.df.hn.yr.14$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.yr.14$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.cv.14 <- 
  ycg.allstrat.df.hn.yr.14$dht$individuals$N$cv[
   ycg.allstrat.df.hn.yr.14$dht$individuals$N$Label == "Total"]

# 2016
YCG.allstrat.yr.ind.abund.16 <- 
  ycg.allstrat.df.hn.yr.16$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.yr.16$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.LCL.16 <- 
  ycg.allstrat.df.hn.yr.16$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.yr.16$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.UCL.16 <- 
  ycg.allstrat.df.hn.yr.16$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.yr.16$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.cv.16 <- 
  ycg.allstrat.df.hn.yr.16$dht$individuals$N$cv[
   ycg.allstrat.df.hn.yr.16$dht$individuals$N$Label == "Total"]

# 2018
YCG.allstrat.yr.ind.abund.18 <- 
  ycg.allstrat.df.hn.yr.18$dht$individuals$N$Estimate[
    ycg.allstrat.df.hn.yr.18$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.LCL.18 <- 
  ycg.allstrat.df.hn.yr.18$dht$individuals$N$lcl[
    ycg.allstrat.df.hn.yr.18$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.UCL.18 <- 
  ycg.allstrat.df.hn.yr.18$dht$individuals$N$ucl[
   ycg.allstrat.df.hn.yr.18$dht$individuals$N$Label == "Total"]

YCG.allstrat.yr.cv.18 <- 
  ycg.allstrat.df.hn.yr.18$dht$individuals$N$cv[
   ycg.allstrat.df.hn.yr.18$dht$individuals$N$Label == "Total"]


# dataframe
ycg_plot_dat_allstrat_yr <- data.frame(DetFun = rep("stratYear", times=6),
                              year = c("2010","2011","2013","2014","2016","2018"),
                              abundance = c(YCG.allstrat.yr.ind.abund.10,YCG.allstrat.yr.ind.abund.11,
                                            YCG.allstrat.yr.ind.abund.13,YCG.allstrat.yr.ind.abund.14,
                                            YCG.allstrat.yr.ind.abund.16,YCG.allstrat.yr.ind.abund.18),
                              lcl = c(YCG.allstrat.yr.LCL.10,YCG.allstrat.yr.LCL.11,
                                      YCG.allstrat.yr.LCL.13,YCG.allstrat.yr.LCL.14,
                                      YCG.allstrat.yr.LCL.16,YCG.allstrat.yr.LCL.18),
                              ucl = c(YCG.allstrat.yr.UCL.10,YCG.allstrat.yr.UCL.11,
                                      YCG.allstrat.yr.UCL.13,YCG.allstrat.yr.UCL.14,
                                      YCG.allstrat.yr.UCL.16,YCG.allstrat.yr.UCL.18),
                             cv = c(YCG.allstrat.yr.cv.10,YCG.allstrat.yr.cv.11,YCG.allstrat.yr.cv.13,
                                    YCG.allstrat.yr.cv.14,YCG.allstrat.yr.cv.16,YCG.allstrat.yr.cv.18))

### test analysis with no strata and year as a CONTINUOUS variable in the detection function ####

# make stratum continuous in the data
ycg.data.nostrata.cont <- YCG.data.nostrata
ycg.data.nostrata.cont$stratum <- as.character(ycg.data.nostrata.cont$stratum)
ycg.data.nostrata.cont$stratum <- as.numeric(ycg.data.nostrata.cont$stratum)
str(ycg.data.nostrata.cont)
unique(ycg.data.nostrata.cont$stratum)
ycg.data.nostrata.cont$stratum[ycg.data.nostrata.cont$stratum==2010] <- 1
ycg.data.nostrata.cont$stratum[ycg.data.nostrata.cont$stratum==2011] <- 2
ycg.data.nostrata.cont$stratum[ycg.data.nostrata.cont$stratum==2013] <- 3
ycg.data.nostrata.cont$stratum[ycg.data.nostrata.cont$stratum==2014] <- 4
ycg.data.nostrata.cont$stratum[ycg.data.nostrata.cont$stratum==2016] <- 5
ycg.data.nostrata.cont$stratum[ycg.data.nostrata.cont$stratum==2018] <- 6

# half normal
ycg.nostrat.df.hn.yr.cont <- ds(data=ycg.data.nostrata.cont, region.table=region.table.nostrata, 
                             sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                             truncation=60, key="hn", formula = ~stratum)

summary(ycg.nostrat.df.hn.yr.cont)
# the coefficient for the continuous year predictor is 0.1 (se=0.05), and so this suggests that time is not having a major impact on p

# extract estimates

# 2010
YCG.nostrat.yrCont.ind.abund.10 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2010"]

YCG.nostrat.yrCont.LCL.10 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2010"]

YCG.nostrat.yrCont.UCL.10 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2010"]

YCG.nostrat.yrCont.cv.10 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2010"]


# 2011
YCG.nostrat.yrCont.ind.abund.11 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2011"]

YCG.nostrat.yrCont.LCL.11 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2011"]

YCG.nostrat.yrCont.UCL.11 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2011"]

YCG.nostrat.yrCont.cv.11 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2011"]


# 2013
YCG.nostrat.yrCont.ind.abund.13 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2013"]

YCG.nostrat.yrCont.LCL.13 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2013"]

YCG.nostrat.yrCont.UCL.13 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2013"]

YCG.nostrat.yrCont.cv.13 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2013"]


# 2014
YCG.nostrat.yrCont.ind.abund.14 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2014"]

YCG.nostrat.yrCont.LCL.14 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2014"]

YCG.nostrat.yrCont.UCL.14 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2014"]

YCG.nostrat.yrCont.cv.14 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2014"]


# 2016
YCG.nostrat.yrCont.ind.abund.16 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2016"]

YCG.nostrat.yrCont.LCL.16 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2016"]

YCG.nostrat.yrCont.UCL.16 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2016"]

YCG.nostrat.yrCont.cv.16 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2016"]


# 2018
YCG.nostrat.yrCont.ind.abund.18 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2018"]

YCG.nostrat.yrCont.LCL.18 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2018"]

YCG.nostrat.yrCont.UCL.18 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2018"]

YCG.nostrat.yrCont.cv.18 <- 
  ycg.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   ycg.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2018"]



ycg_plot_dat_yr_cont <- data.frame(analysis = rep("NoStratYearCont", times=6),
                           year = c("2010","2011","2013","2014","2016","2018"),
                           abundance = c(YCG.nostrat.yrCont.ind.abund.10,YCG.nostrat.yrCont.ind.abund.11,
                                         YCG.nostrat.yrCont.ind.abund.13,YCG.nostrat.yrCont.ind.abund.14,
                                         YCG.nostrat.yrCont.ind.abund.16,YCG.nostrat.yrCont.ind.abund.18),
                           lcl = c(YCG.nostrat.yrCont.LCL.10,YCG.nostrat.yrCont.LCL.11,
                                   YCG.nostrat.yrCont.LCL.13,YCG.nostrat.yrCont.LCL.14,
                                   YCG.nostrat.yrCont.LCL.16,YCG.nostrat.yrCont.LCL.18),
                           ucl = c(YCG.nostrat.yrCont.UCL.10,YCG.nostrat.yrCont.UCL.11,
                                   YCG.nostrat.yrCont.UCL.13,YCG.nostrat.yrCont.UCL.14,
                                   YCG.nostrat.yrCont.UCL.16,YCG.nostrat.yrCont.UCL.18),
                           cv = c(YCG.nostrat.yrCont.cv.10,YCG.nostrat.yrCont.cv.11,
                                  YCG.nostrat.yrCont.cv.13,YCG.nostrat.yrCont.cv.14,
                                  YCG.nostrat.yrCont.cv.16,YCG.nostrat.yrCont.cv.18))


### plotting all analyses and model comparison ####

ycg_plot_dat_merge <- rbind(ycg_plot_dat_pool,ycg_plot_dat_yr)
ycg_plot_dat_merge <- rbind(ycg_plot_dat_merge,ycg_plot_dat_allstrat_pool)
ycg_plot_dat_merge <- rbind(ycg_plot_dat_merge,ycg_plot_dat_allstrat_yr)
ycg_plot_dat_merge$DetFun <- as.character(ycg_plot_dat_merge$DetFun)

ycg_plot_dat_merge <- ycg_plot_dat_merge %>% 
                      rename(analysis = DetFun) %>% 
                      mutate(analysis = ifelse(analysis=="pool",
                                               replace(analysis,analysis=="pool","NoStratPool"),
                                               analysis)) %>% 
                      mutate(analysis = ifelse(analysis=="yr",
                                               replace(analysis, analysis=="yr","NoStratYear"),
                                               analysis)) 

# add dataframe from nostrata and year as continous
ycg_plot_dat_merge <- rbind(ycg_plot_dat_merge,ycg_plot_dat_yr_cont)

# plot estimates
estplot <- ggplot(ycg_plot_dat_merge, aes(x=year,y=abundance, group=analysis, colour=analysis))+
           geom_point(position = position_dodge(width = 0.5),shape=16, size=2)+
           geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3, position = position_dodge(width = 0.5))+
           theme_bw()+
           theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),panel.border = element_blank())+
           theme(axis.line = element_line(color = 'black'))+
           stat_smooth(method='lm', formula = y~x, se=F)

# plot CV
cvplot <- ggplot(ycg_plot_dat_merge, aes(x=year, y=cv, group=analysis, colour=analysis))+
          geom_point(shape=16, size=5)+
          theme_bw()+
           theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),panel.border = element_blank())+
           theme(axis.line = element_line(color = 'black'))

plot_grid(estplot,cvplot)


## compare models with AIC

# add AIC values for stratified analysis

# allstrat pooled
allstrat.pool.AIC <- sum(AIC(ycg.allstrat.df.hn.co.10,ycg.allstrat.df.hn.co.11,
                              ycg.allstrat.df.hn.co.13,ycg.allstrat.df.hn.co.14,
                              ycg.allstrat.df.hn.co.16,ycg.allstrat.df.hn.co.18)$AIC)

# allstrat year
allstrat.yr.AIC <- sum(AIC(ycg.allstrat.df.hn.yr.10,ycg.allstrat.df.hn.yr.11,ycg.allstrat.df.hn.yr.13,
                           ycg.allstrat.df.hn.yr.14,ycg.allstrat.df.hn.yr.16,ycg.allstrat.df.hn.yr.18)$AIC)

# nostrat pooled
nostrat.pool.AIC <- AIC(ycg.nostrat.df.hn.co)$AIC

# nostrat year
nostrat.yr.AIC <- AIC(ycg.nostrat.df.hn.yr)$AIC

# nostrat year continous
nostrat.yrCont.AIC <- AIC(ycg.nostrat.df.hn.yr.cont)$AIC

comp.df <- data.frame(analysis = c("allstrat.pool.AIC","allstrat.yr.AIC","nostrat.pool.AIC",
                                   "nostrat.yr.AIC","nostrat.yrCont.AIC"),
                      aic = c(allstrat.pool.AIC,allstrat.yr.AIC,nostrat.pool.AIC,nostrat.yr.AIC,
                              nostrat.yrCont.AIC))

#### Green peafowl ###########################################
### subset data, plot histograms, select detection function ####

## subset data for GPF
GPF.data.nostrata <-data.noStrata[data.noStrata$species=="GPF",]

# plot histogram
hist(GPF.data.nostrata$distance)



### Test analysis with NO priority strata and pooled detection function across years ####
  ## select detection function ####

## fit detection function (no uniform as we will need covars)


# half normal
gpf.nostrat.df.hn.co <- ds(data=GPF.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=80, key="hn", adjustment= "cos")

gpf.nostrat.df.hn.herm <- ds(data=GPF.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=80, key="hn", adjustment= "herm")

# hazard rate
gpf.nostrat.df.hr.co <- ds(data=GPF.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=80, key="hr", adjustment= "cos")

gpf.nostrat.df.hr.poly <- ds(data=GPF.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=80, key="hr", adjustment= "poly")

# compare models
gpf.nostrat.comp <- summarize_ds_models(gpf.nostrat.df.hn.co,gpf.nostrat.df.hn.herm,
                                        gpf.nostrat.df.hr.co,gpf.nostrat.df.hr.poly,
                                             output = "plain")
gpf.nostrat.comp[ ,1:5]
gpf.nostrat.comp[ ,6:7]

par(mfrow=c(2,2))
plot(gpf.nostrat.df.hn.co)
plot(gpf.nostrat.df.hn.herm)
plot(gpf.nostrat.df.hr.co)
plot(gpf.nostrat.df.hr.poly)
# hn co selected


  ## extract estimates ####


## extract estimates, CIs, CVs

# 2010
GPF.nostrat.pooled.ind.abund.10 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2010"]

GPF.nostrat.pooled.LCL.10 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2010"]

GPF.nostrat.pooled.UCL.10 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2010"]

GPF.nostrat.pooled.CV.10 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$cv[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2010"]


# 2011
GPF.nostrat.pooled.ind.abund.11 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2011"]

GPF.nostrat.pooled.LCL.11 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2011"]

GPF.nostrat.pooled.UCL.11 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2011"]

GPF.nostrat.pooled.CV.11 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$cv[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2011"]

# 2013
GPF.nostrat.pooled.ind.abund.13 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2013"]

GPF.nostrat.pooled.LCL.13 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2013"]

GPF.nostrat.pooled.UCL.13 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2013"]

GPF.nostrat.pooled.CV.13 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$cv[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2013"]

# 2014
GPF.nostrat.pooled.ind.abund.14 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2014"]

GPF.nostrat.pooled.LCL.14 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2014"]

GPF.nostrat.pooled.UCL.14 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2014"]

GPF.nostrat.pooled.CV.14 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$cv[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2014"]

# 2016
GPF.nostrat.pooled.ind.abund.16 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2016"]

GPF.nostrat.pooled.LCL.16 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2016"]

GPF.nostrat.pooled.UCL.16 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2016"]

GPF.nostrat.pooled.CV.16 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$cv[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2016"]

# 2018
GPF.nostrat.pooled.ind.abund.18 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2018"]

GPF.nostrat.pooled.LCL.18 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2018"]

GPF.nostrat.pooled.UCL.18 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2018"]

GPF.nostrat.pooled.CV.18 <- 
  gpf.nostrat.df.hn.co$dht$individuals$N$cv[
   gpf.nostrat.df.hn.co$dht$individuals$N$Label == "2018"]


gpf_plot_dat_pool <- data.frame(DetFun = rep("pool", times=6),
                           year = c("2010","2011","2013","2014","2016","2018"),
                           abundance = c(GPF.nostrat.pooled.ind.abund.10,GPF.nostrat.pooled.ind.abund.11,
                                         GPF.nostrat.pooled.ind.abund.13,GPF.nostrat.pooled.ind.abund.14,
                                         GPF.nostrat.pooled.ind.abund.16,GPF.nostrat.pooled.ind.abund.18),
                           lcl = c(GPF.nostrat.pooled.LCL.10,GPF.nostrat.pooled.LCL.11,
                                   GPF.nostrat.pooled.LCL.13,GPF.nostrat.pooled.LCL.14,
                                   GPF.nostrat.pooled.LCL.16,GPF.nostrat.pooled.LCL.18),
                           ucl = c(GPF.nostrat.pooled.UCL.10,GPF.nostrat.pooled.UCL.11,
                                   GPF.nostrat.pooled.UCL.13,GPF.nostrat.pooled.UCL.14,
                                   GPF.nostrat.pooled.UCL.16,GPF.nostrat.pooled.UCL.18),
                           cv = c(GPF.nostrat.pooled.CV.10,GPF.nostrat.pooled.CV.11,
                                  GPF.nostrat.pooled.CV.13,GPF.nostrat.pooled.CV.14,
                                  GPF.nostrat.pooled.CV.16,GPF.nostrat.pooled.CV.18))
### test analysis with no priority strata but year as a covariate in detection function ####


# half normal
gpf.nostrat.df.hn.yr <- ds(data=GPF.data.nostrata, region.table=region.table.nostrata, 
                      sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                      truncation=80, key="hn", formula = ~stratum)

# compare model with year as covariate and model with pooled DF
gpf.nostrat.comp.pool <- summarize_ds_models(gpf.nostrat.df.hn.co,gpf.nostrat.df.hn.yr,
                                             output = "plain")
gpf.nostrat.comp.pool[ ,1:5]
gpf.nostrat.comp.pool[ ,6:7]
# model with pooled DF is much better

# extract estimats from model with stratum as covariate

summary(gpf.nostrat.df.hn.yr)

# 2010
GPF.nostrat.yr.ind.abund.10 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2010"]

GPF.nostrat.yr.LCL.10 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2010"]

GPF.nostrat.yr.UCL.10 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2010"]

GPF.nostrat.yr.cv.10 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2010"]


# 2011
GPF.nostrat.yr.ind.abund.11 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2011"]

GPF.nostrat.yr.LCL.11 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2011"]

GPF.nostrat.yr.UCL.11 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2011"]

GPF.nostrat.yr.cv.11 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2011"]

# 2013
GPF.nostrat.yr.ind.abund.13 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2013"]

GPF.nostrat.yr.LCL.13 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2013"]

GPF.nostrat.yr.UCL.13 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2013"]

GPF.nostrat.yr.cv.13 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2013"]

# 2014
GPF.nostrat.yr.ind.abund.14 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2014"]

GPF.nostrat.yr.LCL.14 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2014"]

GPF.nostrat.yr.UCL.14 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2014"]

GPF.nostrat.yr.cv.14 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2014"]

# 2016
GPF.nostrat.yr.ind.abund.16 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2016"]

GPF.nostrat.yr.LCL.16 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2016"]

GPF.nostrat.yr.UCL.16 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2016"]

GPF.nostrat.yr.cv.16 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2016"]

# 2018
GPF.nostrat.yr.ind.abund.18 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2018"]

GPF.nostrat.yr.LCL.18 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2018"]

GPF.nostrat.yr.UCL.18 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2018"]

GPF.nostrat.yr.cv.18 <- 
  gpf.nostrat.df.hn.yr$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr$dht$individuals$N$Label == "2018"]


gpf_plot_dat_yr <- data.frame(DetFun = rep("yr", times=6),
                           year = c("2010","2011","2013","2014","2016","2018"),
                           abundance = c(GPF.nostrat.yr.ind.abund.10,GPF.nostrat.yr.ind.abund.11,
                                         GPF.nostrat.yr.ind.abund.13,GPF.nostrat.yr.ind.abund.14,
                                         GPF.nostrat.yr.ind.abund.16,GPF.nostrat.yr.ind.abund.18),
                           lcl = c(GPF.nostrat.yr.LCL.10,GPF.nostrat.yr.LCL.11,
                                   GPF.nostrat.yr.LCL.13,GPF.nostrat.yr.LCL.14,
                                   GPF.nostrat.yr.LCL.16,GPF.nostrat.yr.LCL.18),
                           ucl = c(GPF.nostrat.yr.UCL.10,GPF.nostrat.yr.UCL.11,
                                   GPF.nostrat.yr.UCL.13,GPF.nostrat.yr.UCL.14,
                                   GPF.nostrat.yr.UCL.16,GPF.nostrat.yr.UCL.18),
                           cv = c(GPF.nostrat.yr.cv.10,GPF.nostrat.yr.cv.11,
                                  GPF.nostrat.yr.cv.13,GPF.nostrat.yr.cv.14,
                                  GPF.nostrat.yr.cv.16,GPF.nostrat.yr.cv.18))

# merge the pooled and strtified dataframes
gpf_plot_dat_merge <- rbind(gpf_plot_dat_pool,gpf_plot_dat_yr)

# plot together
ggplot(gpf_plot_dat_merge, aes(x=year, y=abundance, group = DetFun, colour=DetFun))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3, position = position_dodge(width = 0.5))
 
 


### test analysis with priority strata in ALL years and pooled DF ####


# Here I am going to create priority strata for all years.  For 2011, 2014, and 2016 I will use the correct strata (i.e. the ones that were done).  For 2010 and 2013 I will use the strata from 2011. For 2018 I will use the strata from 2016. 

# subset data for GPF
GPF.data.allstrata <- data_allstrata[data_allstrata$species=="GPF", ]

   ## fit detection function and get estimates ####

gpf.allstrat.df.hn.co <- ds(data=GPF.data.allstrata, region.table=region.table.allstrata, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", adjustment= "cos")

summary(gpf.allstrat.df.hn.co)

# the SE and CV's are not just the sums or means of the two strata, and so I have to go through to proces of setting all other strata areas to 0 for each year.

# set regions to 0 for 2011-2018
region.table.GPF.allstrata.pool.10 <- region.table.allstrata
region.table.GPF.allstrata.pool.10$Area[region.table.GPF.allstrata.pool.10$Region.Label != "2010_L" & 
                                          region.table.GPF.allstrata.pool.10$Region.Label != "2010_H"] <- 0

# re-run model
gpf.allstrat.df.hn.co.10 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.10, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", adjustment= "cos")
summary(gpf.allstrat.df.hn.co.10)

# set regions to 0 for all except 2011
region.table.GPF.allstrata.pool.11 <- region.table.allstrata
region.table.GPF.allstrata.pool.11$Area[region.table.GPF.allstrata.pool.11$Region.Label != "2011_L" & 
                                          region.table.GPF.allstrata.pool.11$Region.Label != "2011_H"] <- 0

# re-run model
gpf.allstrat.df.hn.co.11 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.11, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", adjustment= "cos")
summary(gpf.allstrat.df.hn.co.11)

# set regions to 0 for all except 2013
region.table.GPF.allstrata.pool.13 <- region.table.allstrata
region.table.GPF.allstrata.pool.13$Area[region.table.GPF.allstrata.pool.13$Region.Label != "2013_L" & 
                                          region.table.GPF.allstrata.pool.13$Region.Label != "2013_H"] <- 0

# re-run model
gpf.allstrat.df.hn.co.13 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.13, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", adjustment= "cos")
summary(gpf.allstrat.df.hn.co.13)


# set regions to 0 for all except 2014
region.table.GPF.allstrata.pool.14 <- region.table.allstrata
region.table.GPF.allstrata.pool.14$Area[region.table.GPF.allstrata.pool.14$Region.Label != "2014_L" & 
                                          region.table.GPF.allstrata.pool.14$Region.Label != "2014_H"] <- 0


# re-run model
gpf.allstrat.df.hn.co.14 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.14, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", adjustment= "cos")
summary(gpf.allstrat.df.hn.co.14)


# set regions to 0 for all except 2016
region.table.GPF.allstrata.pool.16 <- region.table.allstrata
region.table.GPF.allstrata.pool.16$Area[region.table.GPF.allstrata.pool.16$Region.Label != "2016_L" & 
                                          region.table.GPF.allstrata.pool.16$Region.Label != "2016_H"] <- 0


# re-run model
gpf.allstrat.df.hn.co.16 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.16, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", adjustment= "cos")
summary(gpf.allstrat.df.hn.co.16)


# set regions to 0 for all except 2018
region.table.GPF.allstrata.pool.18 <- region.table.allstrata
region.table.GPFGPF.allstrata.pool.18$Area[region.table.GPF.allstrata.pool.18$Region.Label != "2018_L" & 
                                          region.table.GPF.allstrata.pool.18$Region.Label != "2018_H"] <- 0


# re-run model
gpf.allstrat.df.hn.co.18 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.18, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", adjustment= "cos")
summary(gpf.allstrat.df.hn.co.18)


# pull out estimates

# 2010
GPF.allstrat.pool.ind.abund.10 <- 
  gpf.allstrat.df.hn.co.10$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.co.10$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.LCL.10 <- 
  gpf.allstrat.df.hn.co.10$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.co.10$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.UCL.10 <- 
  gpf.allstrat.df.hn.co.10$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.co.10$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.cv.10 <- 
  gpf.allstrat.df.hn.co.10$dht$individuals$N$cv[
   gpf.allstrat.df.hn.co.10$dht$individuals$N$Label == "Total"]

# 2011
GPF.allstrat.pool.ind.abund.11 <- 
  gpf.allstrat.df.hn.co.11$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.co.11$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.LCL.11 <- 
  gpf.allstrat.df.hn.co.11$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.co.11$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.UCL.11 <- 
  gpf.allstrat.df.hn.co.11$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.co.11$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.cv.11 <- 
  gpf.allstrat.df.hn.co.11$dht$individuals$N$cv[
   gpf.allstrat.df.hn.co.11$dht$individuals$N$Label == "Total"]

# 2013
GPF.allstrat.pool.ind.abund.13 <- 
  gpf.allstrat.df.hn.co.13$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.co.13$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.LCL.13 <- 
  gpf.allstrat.df.hn.co.13$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.co.13$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.UCL.13 <- 
  gpf.allstrat.df.hn.co.13$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.co.13$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.cv.13 <- 
  gpf.allstrat.df.hn.co.13$dht$individuals$N$cv[
   gpf.allstrat.df.hn.co.13$dht$individuals$N$Label == "Total"]

# 2014
GPF.allstrat.pool.ind.abund.14 <- 
  gpf.allstrat.df.hn.co.14$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.co.14$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.LCL.14 <- 
  gpf.allstrat.df.hn.co.14$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.co.14$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.UCL.14 <- 
  gpf.allstrat.df.hn.co.14$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.co.14$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.cv.14 <- 
  gpf.allstrat.df.hn.co.14$dht$individuals$N$cv[
   gpf.allstrat.df.hn.co.14$dht$individuals$N$Label == "Total"]

# 2016
GPF.allstrat.pool.ind.abund.16 <- 
  gpf.allstrat.df.hn.co.16$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.co.16$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.LCL.16 <- 
  gpf.allstrat.df.hn.co.16$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.co.16$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.UCL.16 <- 
  gpf.allstrat.df.hn.co.16$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.co.16$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.cv.16 <- 
  gpf.allstrat.df.hn.co.16$dht$individuals$N$cv[
   gpf.allstrat.df.hn.co.16$dht$individuals$N$Label == "Total"]

# 2018
GPF.allstrat.pool.ind.abund.18 <- 
  gpf.allstrat.df.hn.co.18$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.co.18$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.LCL.18 <- 
  gpf.allstrat.df.hn.co.18$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.co.18$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.UCL.18 <- 
  gpf.allstrat.df.hn.co.18$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.co.18$dht$individuals$N$Label == "Total"]

GPF.allstrat.pool.cv.18 <- 
  gpf.allstrat.df.hn.co.18$dht$individuals$N$cv[
   gpf.allstrat.df.hn.co.18$dht$individuals$N$Label == "Total"]


# dataframe
gpf_plot_dat_allstrat_pool <- data.frame(DetFun = rep("stratPool", times=6),
                              year = c("2010","2011","2013","2014","2016","2018"),
                              abundance = c(GPF.allstrat.pool.ind.abund.10,GPF.allstrat.pool.ind.abund.11,
                                            GPF.allstrat.pool.ind.abund.13,GPF.allstrat.pool.ind.abund.14,
                                            GPF.allstrat.pool.ind.abund.16,GPF.allstrat.pool.ind.abund.18),
                              lcl = c(GPF.allstrat.pool.LCL.10,GPF.allstrat.pool.LCL.11,
                                      GPF.allstrat.pool.LCL.13,GPF.allstrat.pool.LCL.14,
                                      GPF.allstrat.pool.LCL.16,GPF.allstrat.pool.LCL.18),
                              ucl = c(GPF.allstrat.pool.UCL.10,GPF.allstrat.pool.UCL.11,
                                      GPF.allstrat.pool.UCL.13,GPF.allstrat.pool.UCL.14,
                                      GPF.allstrat.pool.UCL.16,GPF.allstrat.pool.UCL.18),
                             cv = c(GPF.allstrat.pool.cv.10,GPF.allstrat.pool.cv.11,
                                    GPF.allstrat.pool.cv.13,GPF.allstrat.pool.cv.14,
                                    GPF.allstrat.pool.cv.16,GPF.allstrat.pool.cv.18))

### test anlaysis with priority strata in ALL years and year as covariate in detection function ####

# fit detection function

gpf.allstrat.df.hn.yr <- ds(data=GPF.data.allstrata, region.table=region.table.allstrata, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", formula = ~stratum)

summary(gpf.allstrat.df.hn.yr)


# re-run model for 2010 estimates
gpf.allstrat.df.hn.yr.10 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.10, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", formula = ~stratum)
summary(gpf.allstrat.df.hn.yr.10)


# re-run model for 2011 estimates
gpf.allstrat.df.hn.yr.11 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.11, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", formula = ~stratum)
summary(gpf.allstrat.df.hn.yr.11)

# re-run model for 2013 estimates
gpf.allstrat.df.hn.yr.13 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.13, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", formula = ~stratum)
summary(gpf.allstrat.df.hn.yr.13)

# re-run model for 2014 estimates
gpf.allstrat.df.hn.yr.14 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.14, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", formula = ~stratum)
summary(gpf.allstrat.df.hn.yr.14)

# re-run model for 2016 estimates
gpf.allstrat.df.hn.yr.16 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.16, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", formula = ~stratum)
summary(gpf.allstrat.df.hn.yr.16)

# re-run model for 2018 estimates
gpf.allstrat.df.hn.yr.18 <- ds(data=GPF.data.allstrata, region.table=region.table.GPF.allstrata.pool.18, 
                         sample.table=sample.table.allstrata, obs.table=obs.table.allstrata, 
                         truncation=80, key="hn", formula = ~stratum)
summary(gpf.allstrat.df.hn.yr.18)


# pull out estimates

# 2010
GPF.allstrat.yr.ind.abund.10 <- 
  gpf.allstrat.df.hn.yr.10$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.yr.10$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.LCL.10 <- 
  gpf.allstrat.df.hn.yr.10$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.yr.10$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.UCL.10 <- 
  gpf.allstrat.df.hn.yr.10$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.yr.10$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.cv.10 <- 
  gpf.allstrat.df.hn.yr.10$dht$individuals$N$cv[
   gpf.allstrat.df.hn.yr.10$dht$individuals$N$Label == "Total"]

# 2011
GPF.allstrat.yr.ind.abund.11 <- 
  gpf.allstrat.df.hn.yr.11$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.yr.11$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.LCL.11 <- 
  gpf.allstrat.df.hn.yr.11$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.yr.11$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.UCL.11 <- 
  gpf.allstrat.df.hn.yr.11$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.yr.11$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.cv.11 <- 
  gpf.allstrat.df.hn.yr.11$dht$individuals$N$cv[
   gpf.allstrat.df.hn.yr.11$dht$individuals$N$Label == "Total"]

# 2013
GPF.allstrat.yr.ind.abund.13 <- 
  gpf.allstrat.df.hn.yr.13$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.yr.13$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.LCL.13 <- 
  gpf.allstrat.df.hn.yr.13$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.yr.13$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.UCL.13 <- 
  gpf.allstrat.df.hn.yr.13$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.yr.13$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.cv.13 <- 
  gpf.allstrat.df.hn.yr.13$dht$individuals$N$cv[
   gpf.allstrat.df.hn.yr.13$dht$individuals$N$Label == "Total"]

# 2014
GPF.allstrat.yr.ind.abund.14 <- 
  gpf.allstrat.df.hn.yr.14$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.yr.14$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.LCL.14 <- 
  gpf.allstrat.df.hn.yr.14$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.yr.14$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.UCL.14 <- 
  gpf.allstrat.df.hn.yr.14$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.yr.14$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.cv.14 <- 
  gpf.allstrat.df.hn.yr.14$dht$individuals$N$cv[
   gpf.allstrat.df.hn.yr.14$dht$individuals$N$Label == "Total"]

# 2016
GPF.allstrat.yr.ind.abund.16 <- 
  gpf.allstrat.df.hn.yr.16$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.yr.16$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.LCL.16 <- 
  gpf.allstrat.df.hn.yr.16$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.yr.16$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.UCL.16 <- 
  gpf.allstrat.df.hn.yr.16$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.yr.16$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.cv.16 <- 
  gpf.allstrat.df.hn.yr.16$dht$individuals$N$cv[
   gpf.allstrat.df.hn.yr.16$dht$individuals$N$Label == "Total"]

# 2018
GPF.allstrat.yr.ind.abund.18 <- 
  gpf.allstrat.df.hn.yr.18$dht$individuals$N$Estimate[
    gpf.allstrat.df.hn.yr.18$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.LCL.18 <- 
  gpf.allstrat.df.hn.yr.18$dht$individuals$N$lcl[
    gpf.allstrat.df.hn.yr.18$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.UCL.18 <- 
  gpf.allstrat.df.hn.yr.18$dht$individuals$N$ucl[
   gpf.allstrat.df.hn.yr.18$dht$individuals$N$Label == "Total"]

GPF.allstrat.yr.cv.18 <- 
  gpf.allstrat.df.hn.yr.18$dht$individuals$N$cv[
   gpf.allstrat.df.hn.yr.18$dht$individuals$N$Label == "Total"]


# dataframe
gpf_plot_dat_allstrat_yr <- data.frame(DetFun = rep("stratYear", times=6),
                              year = c("2010","2011","2013","2014","2016","2018"),
                              abundance = c(GPF.allstrat.yr.ind.abund.10,GPF.allstrat.yr.ind.abund.11,
                                            GPF.allstrat.yr.ind.abund.13,GPF.allstrat.yr.ind.abund.14,
                                            GPF.allstrat.yr.ind.abund.16,GPF.allstrat.yr.ind.abund.18),
                              lcl = c(GPF.allstrat.yr.LCL.10,GPF.allstrat.yr.LCL.11,
                                      GPF.allstrat.yr.LCL.13,GPF.allstrat.yr.LCL.14,
                                      GPF.allstrat.yr.LCL.16,GPF.allstrat.yr.LCL.18),
                              ucl = c(GPF.allstrat.yr.UCL.10,GPF.allstrat.yr.UCL.11,
                                      GPF.allstrat.yr.UCL.13,GPF.allstrat.yr.UCL.14,
                                      GPF.allstrat.yr.UCL.16,GPF.allstrat.yr.UCL.18),
                             cv = c(GPF.allstrat.yr.cv.10,GPF.allstrat.yr.cv.11,GPF.allstrat.yr.cv.13,
                                    GPF.allstrat.yr.cv.14,GPF.allstrat.yr.cv.16,GPF.allstrat.yr.cv.18))

### test analysis with no strata and year as a CONTINUOUS variable in the detection function ####

# make stratum continuous in the data
gpf.data.nostrata.cont <- GPF.data.nostrata
gpf.data.nostrata.cont$stratum <- as.character(gpf.data.nostrata.cont$stratum)
gpf.data.nostrata.cont$stratum <- as.numeric(gpf.data.nostrata.cont$stratum)
str(gpf.data.nostrata.cont)
unique(gpf.data.nostrata.cont$stratum)
gpf.data.nostrata.cont$stratum[gpf.data.nostrata.cont$stratum==2010] <- 1
gpf.data.nostrata.cont$stratum[gpf.data.nostrata.cont$stratum==2011] <- 2
gpf.data.nostrata.cont$stratum[gpf.data.nostrata.cont$stratum==2013] <- 3
gpf.data.nostrata.cont$stratum[gpf.data.nostrata.cont$stratum==2014] <- 4
gpf.data.nostrata.cont$stratum[gpf.data.nostrata.cont$stratum==2016] <- 5
gpf.data.nostrata.cont$stratum[gpf.data.nostrata.cont$stratum==2018] <- 6

# half normal
gpf.nostrat.df.hn.yr.cont <- ds(data=gpf.data.nostrata.cont, region.table=region.table.nostrata, 
                             sample.table=sample.table.nostrata, obs.table=obs.table.nostrata, 
                             truncation=80, key="hn", formula = ~stratum)

summary(gpf.nostrat.df.hn.yr.cont)
# the coefficient for the continuous year predictor is 0.02 (se=0.05), and so this suggests that time is not having a major impact on p

# extract estimates

# 2010
GPF.nostrat.yrCont.ind.abund.10 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2010"]

GPF.nostrat.yrCont.LCL.10 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2010"]

GPF.nostrat.yrCont.UCL.10 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2010"]

GPF.nostrat.yrCont.cv.10 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2010"]


# 2011
GPF.nostrat.yrCont.ind.abund.11 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2011"]

GPF.nostrat.yrCont.LCL.11 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2011"]

GPF.nostrat.yrCont.UCL.11 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2011"]

GPF.nostrat.yrCont.cv.11 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2011"]


# 2013
GPF.nostrat.yrCont.ind.abund.13 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2013"]

GPF.nostrat.yrCont.LCL.13 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2013"]

GPF.nostrat.yrCont.UCL.13 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2013"]

GPF.nostrat.yrCont.cv.13 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2013"]


# 2014
GPF.nostrat.yrCont.ind.abund.14 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2014"]

GPF.nostrat.yrCont.LCL.14 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2014"]

GPF.nostrat.yrCont.UCL.14 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2014"]

GPF.nostrat.yrCont.cv.14 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2014"]


# 2016
GPF.nostrat.yrCont.ind.abund.16 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2016"]

GPF.nostrat.yrCont.LCL.16 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2016"]

GPF.nostrat.yrCont.UCL.16 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2016"]

GPF.nostrat.yrCont.cv.16 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2016"]


# 2018
GPF.nostrat.yrCont.ind.abund.18 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Estimate[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2018"]

GPF.nostrat.yrCont.LCL.18 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$lcl[
    gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2018"]

GPF.nostrat.yrCont.UCL.18 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$ucl[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2018"]

GPF.nostrat.yrCont.cv.18 <- 
  gpf.nostrat.df.hn.yr.cont$dht$individuals$N$cv[
   gpf.nostrat.df.hn.yr.cont$dht$individuals$N$Label == "2018"]



gpf_plot_dat_yr_cont <- data.frame(analysis = rep("NoStratYearCont", times=6),
                           year = c("2010","2011","2013","2014","2016","2018"),
                           abundance = c(GPF.nostrat.yrCont.ind.abund.10,GPF.nostrat.yrCont.ind.abund.11,
                                         GPF.nostrat.yrCont.ind.abund.13,GPF.nostrat.yrCont.ind.abund.14,
                                         GPF.nostrat.yrCont.ind.abund.16,GPF.nostrat.yrCont.ind.abund.18),
                           lcl = c(GPF.nostrat.yrCont.LCL.10,GPF.nostrat.yrCont.LCL.11,
                                   GPF.nostrat.yrCont.LCL.13,GPF.nostrat.yrCont.LCL.14,
                                   GPF.nostrat.yrCont.LCL.16,GPF.nostrat.yrCont.LCL.18),
                           ucl = c(GPF.nostrat.yrCont.UCL.10,GPF.nostrat.yrCont.UCL.11,
                                   GPF.nostrat.yrCont.UCL.13,GPF.nostrat.yrCont.UCL.14,
                                   GPF.nostrat.yrCont.UCL.16,GPF.nostrat.yrCont.UCL.18),
                           cv = c(GPF.nostrat.yrCont.cv.10,GPF.nostrat.yrCont.cv.11,
                                  GPF.nostrat.yrCont.cv.13,GPF.nostrat.yrCont.cv.14,
                                  GPF.nostrat.yrCont.cv.16,GPF.nostrat.yrCont.cv.18))


## plotting all analyses and model comparison ####

gpf_plot_dat_merge <- rbind(gpf_plot_dat_pool,gpf_plot_dat_yr)
gpf_plot_dat_merge <- rbind(gpf_plot_dat_merge,gpf_plot_dat_allstrat_pool)
gpf_plot_dat_merge <- rbind(gpf_plot_dat_merge,gpf_plot_dat_allstrat_yr)
gpf_plot_dat_merge$DetFun <- as.character(gpf_plot_dat_merge$DetFun)

gpf_plot_dat_merge <- gpf_plot_dat_merge %>% 
                      rename(analysis = DetFun) %>% 
                      mutate(analysis = ifelse(analysis=="pool",
                                               replace(analysis,analysis=="pool","NoStratPool"),
                                               analysis)) %>% 
                      mutate(analysis = ifelse(analysis=="yr",
                                               replace(analysis, analysis=="yr","NoStratYear"),
                                               analysis)) 

# add dataframe from nostrata and year as continous
gpf_plot_dat_merge <- rbind(gpf_plot_dat_merge,gpf_plot_dat_yr_cont)

# plot estimates
estplot.gpf <- ggplot(gpf_plot_dat_merge, aes(x=year,y=abundance, group=analysis, colour=analysis))+
           geom_point(position = position_dodge(width = 0.5),shape=16, size=2)+
           geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.3, position = position_dodge(width = 0.5))+
           theme_bw()+
           theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),panel.border = element_blank())+
           theme(axis.line = element_line(color = 'black'))+
           stat_smooth(method='lm', formula = y~x, se=F)

# plot CV
cvplot.gpf <- ggplot(gpf_plot_dat_merge, aes(x=year, y=cv, group=analysis, colour=analysis))+
          geom_point(shape=16, size=5)+
          theme_bw()+
           theme(plot.background = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),panel.border = element_blank())+
           theme(axis.line = element_line(color = 'black'))

plot_grid(estplot,cvplot)


## compare models with AIC

# add AIC values for stratified analysis

# allstrat pooled
allstrat.pool.AIC.gpf <- sum(AIC(gpf.allstrat.df.hn.co.10,gpf.allstrat.df.hn.co.11,
                              gpf.allstrat.df.hn.co.13,gpf.allstrat.df.hn.co.14,
                              gpf.allstrat.df.hn.co.16,gpf.allstrat.df.hn.co.18)$AIC)

# allstrat year
allstrat.yr.AIC.gpf <- sum(AIC(gpf.allstrat.df.hn.yr.10,gpf.allstrat.df.hn.yr.11,gpf.allstrat.df.hn.yr.13,
                           gpf.allstrat.df.hn.yr.14,gpf.allstrat.df.hn.yr.16,gpf.allstrat.df.hn.yr.18)$AIC)

# nostrat pooled
nostrat.pool.AIC.gpf <- AIC(gpf.nostrat.df.hn.co)$AIC

# nostrat year
nostrat.yr.AIC.gpf <- AIC(gpf.nostrat.df.hn.yr)$AIC

# nostrat year continous
nostrat.yrCont.AIC.gpf <- AIC(gpf.nostrat.df.hn.yr.cont)$AIC

comp.df.gpf <- data.frame(analysis = c("allstrat.pool.AIC","allstrat.yr.AIC","nostrat.pool.AIC",
                                   "nostrat.yr.AIC","nostrat.yrCont.AIC"),
                      aic = c(allstrat.pool.AIC.gpf,allstrat.yr.AIC.gpf,nostrat.pool.AIC.gpf,
                              nostrat.yr.AIC.gpf,nostrat.yrCont.AIC.gpf))

#### Black-shanked douc ######################################
## load data ####


# I want to check the influcence of the strata on BSD.  To do this, I need to use the old master data, because the strata have been removed in the new data. This won't make a difference to any of the other years as BSD is anaysed annually (i.e. the 2020 data does not affect the analysis of previous years)



# load old data (2010-2018) that have strata
load("./Output/Data/Archive/KSWS_MASTER.Rdata")

allData$obs.habitat <- as.factor(allData$obs.habitat)
allData$obs.observer <- as.factor(allData$obs.observer)

# subset BSD data
BSD.data <- allData[allData$species=="BSD",] 


# subset data for 2011
BSD.2011.data <- as.data.frame(BSD.data[BSD.data$stratum=="2011_H" | BSD.data$stratum=="2011_L", ]) 
region.table.2011 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2011_H" 
                                                     | full.region.table$Region.Label=="2011_L",])
sample.table.2011 <- as.data.frame(sample.table[sample.table$Region.Label =="2011_H" | 
                                                  sample.table$Region.Label =="2011_L",])
obs.table.2011 <- as.data.frame(obs.table[obs.table$Region.Label=="2011_H" | obs.table$Region.Label=="2011_L",])


# subset data for 2014
BSD.2014.data <- as.data.frame(BSD.data[BSD.data$stratum=="2014_H" | BSD.data$stratum=="2014_L", ]) 
region.table.2014 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2014_H" 
                                                     | full.region.table$Region.Label=="2014_L",])
sample.table.2014 <- as.data.frame(sample.table[sample.table$Region.Label =="2014_H" | 
                                                  sample.table$Region.Label =="2014_L",])
obs.table.2014 <- as.data.frame(obs.table[obs.table$Region.Label=="2014_H" | obs.table$Region.Label=="2014_L",])


# subset data for 2016
BSD.2016.data <- as.data.frame(BSD.data[BSD.data$stratum=="2016_H" | BSD.data$stratum=="2016_L", ]) 
region.table.2016 <- as.data.frame(full.region.table[full.region.table$Region.Label=="2016_H" 
                                                     | full.region.table$Region.Label=="2016_L",])
sample.table.2016 <- as.data.frame(sample.table[sample.table$Region.Label =="2016_H" | 
                                                  sample.table$Region.Label =="2016_L",])
obs.table.2016 <- as.data.frame(obs.table[obs.table$Region.Label=="2016_H" | obs.table$Region.Label=="2016_L",])



## now load the new data which does not have strata
load("./Output/Data/KSWS_MASTER.Rdata")

allData$obs.habitat <- as.factor(allData$obs.habitat)
allData$obs.observer <- as.factor(allData$obs.observer)
allData$year <- as.factor(allData$stratum)

# subset BSD data
BSD.data <- allData[allData$species=="BSD",] 

# subset 2011 data
BSD.2011.data.new <- as.data.frame(BSD.data[BSD.data$year==2011, ]) 
region.table.2011.new <- as.data.frame(full.region.table[full.region.table$Region.Label=="2011",])
sample.table.2011.new <- as.data.frame(sample.table[sample.table$Region.Label =="2011",])
obs.table.2011.new <- as.data.frame(obs.table[obs.table$Region.Label=="2011",])

# subset 2014 data
BSD.2014.data.new <- as.data.frame(BSD.data[BSD.data$year==2014, ]) 
region.table.2014.new <- as.data.frame(full.region.table[full.region.table$Region.Label=="2014",])
sample.table.2014.new <- as.data.frame(sample.table[sample.table$Region.Label =="2014",])
obs.table.2014.new <- as.data.frame(obs.table[obs.table$Region.Label=="2014",])

# subset 2016 data
BSD.2016.data.new <- as.data.frame(BSD.data[BSD.data$year==2016, ]) 
region.table.2016.new <- as.data.frame(full.region.table[full.region.table$Region.Label=="2016",])
sample.table.2016.new <- as.data.frame(sample.table[sample.table$Region.Label =="2016",])
obs.table.2016.new <- as.data.frame(obs.table[obs.table$Region.Label=="2016",])


## Detection function models ####

# set trunctation distance
BSD.trunc.11 <- 55
BSD.trunc.14 <- 50
BSD.trunc.16 <- 50

# model specification taken from CDS analysis


  ## 2011 detection function models ####

# I want to test the difference between strata and no strata


bsd.11.strata <- ds(data=BSD.2011.data, region.table=region.table.2011, 
                           sample.table=sample.table.2011, obs.table=obs.table.2011, 
                           truncation=BSD.trunc.11, key="hn", adjustment= "cos",
                           cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))

bsd.11.nostrata <- ds(data=BSD.2011.data.new, region.table=region.table.2011.new, 
                      sample.table=sample.table.2011.new, obs.table=obs.table.2011.new, 
                      truncation=BSD.trunc.11, key="hn", adjustment= "cos",
                      cutpoints = c(0,3,7,12,17,22,27,32,37,42,47,55))


  ## 2014 detection function models ####

bsd.14.strata <- ds(data=BSD.2014.data, region.table=region.table.2014, 
                    sample.table=sample.table.2014, obs.table=obs.table.2014, 
                    truncation=BSD.trunc.14, key="hn", formula = ~obs.observer + obs.habitat)

bsd.14.nostrata <- ds(data=BSD.2014.data.new, region.table=region.table.2014.new, 
                      sample.table=sample.table.2014.new, obs.table=obs.table.2014.new, 
                      truncation=BSD.trunc.14, key="hn", formula = ~obs.observer + obs.habitat)

  ## 2016 detection function models ####

bsd.16.strata <- ds(data=BSD.2016.data, region.table=region.table.2016, 
                    sample.table=sample.table.2016, obs.table=obs.table.2016, 
                    truncation=BSD.trunc.16, key="hr", formula = ~obs.AMPM)

bsd.16.nostrata <- ds(data=BSD.2016.data.new, region.table=region.table.2016.new, 
                      sample.table=sample.table.2016.new, obs.table=obs.table.2016.new, 
                      truncation=BSD.trunc.16, key="hr", formula = ~obs.AMPM)



## Plotting ####

# extract estimates

bsd.comp.results <- data.frame(year = rep(c("2011","2014","2016"), each=2),
                               data = rep(c("strata", "no strata"), times=3),
                               estimate = c(bsd.11.strata$dht$clusters$N[3,2],bsd.11.nostrata$dht$clusters$N$Estimate,
                                            bsd.14.strata$dht$clusters$N[3,2],bsd.14.nostrata$dht$clusters$N$Estimate,
                                            bsd.16.strata$dht$clusters$N[3,2],bsd.16.nostrata$dht$clusters$N$Estimate),
                               lcl = c(bsd.11.strata$dht$clusters$N[3,5],bsd.11.nostrata$dht$clusters$N$lcl,
                                       bsd.14.strata$dht$clusters$N[3,5],bsd.14.nostrata$dht$clusters$N$lcl,
                                       bsd.16.strata$dht$clusters$N[3,5],bsd.16.nostrata$dht$clusters$N$lcl),
                               ucl = c(bsd.11.strata$dht$clusters$N[3,6],bsd.11.nostrata$dht$clusters$N$ucl,
                                       bsd.14.strata$dht$clusters$N[3,6],bsd.14.nostrata$dht$clusters$N$ucl,
                                       bsd.16.strata$dht$clusters$N[3,6],bsd.16.nostrata$dht$clusters$N$ucl),
                               cv = c(bsd.11.strata$dht$clusters$N[3,4],bsd.11.nostrata$dht$clusters$N$cv,
                                      bsd.14.strata$dht$clusters$N[3,4],bsd.14.nostrata$dht$clusters$N$cv,
                                      bsd.16.strata$dht$clusters$N[3,4],bsd.16.nostrata$dht$clusters$N$cv))

# the estimates from the DF with strata generally have slightly lower CV


ggplot(bsd.comp.results, aes(x=year, y=estimate, group=data, colour=data))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl),position = position_dodge(width=0.3))+
  ylim(0,16000)


#
#### Yellow-cheecked crested gibbon ##########################
## load data ####

# I want to check the influcence of the strata on YCG.  


  ## Load data for splitting strata and nostrata ####


# load old data (2010-2018) that have strata
load("./Output/Data/Archive/KSWS_MASTER.Rdata")

allData$obs.habitat <- as.factor(allData$obs.habitat)
allData$obs.observer <- as.factor(allData$obs.observer)

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


# now using the updated data, I will get the 2010, 2013, 2018, 2020 data

# load the latest data
load("./Output/Data/KSWS_MASTER.Rdata")

# subset YCG data
YCG.data <- as.data.frame(allData[allData$species=="YCG",]) 


YCG.data.nostrat <- YCG.data[YCG.data$stratum==2010 | YCG.data$stratum==2013 | YCG.data$stratum==2018 |
                             YCG.data$stratum==2020, ]

region.table.nostrat <- as.data.frame(full.region.table[full.region.table$Region.Label=="2010" | 
                                          full.region.table$Region.Label=="2013" |
                                          full.region.table$Region.Label=="2018" |
                                          full.region.table$Region.Label=="2020", ])

sample.table.nostrat <- as.data.frame(sample.table[sample.table$Region.Label==2010 | sample.table$Region.Label==2013 |
                                     sample.table$Region.Label==2018 | sample.table$Region.Label==2020, ])

obs.table.nostrat <- as.data.frame(obs.table[obs.table$Region.Label==2010 | obs.table$Region.Label==2013 |
                               obs.table$Region.Label==2018 | obs.table$Region.Label==2020, ])


# create scaled continuous stratum variable for the detection functions
YCG.data.nostrat$stratum <- as.vector(scale(YCG.data.nostrat$stratum, center = T, scale = T)) 
YCG.data.strat$year <- ifelse(YCG.data.strat$stratum=="2011_L"|YCG.data.strat$stratum=="2011_H",2011,
                        ifelse(YCG.data.strat$stratum=="2014_L"|YCG.data.strat$stratum=="2014_H",2014,
                          ifelse(YCG.data.strat$stratum=="2016_L"|YCG.data.strat$stratum=="2016_H",2016,NA)))
YCG.data.strat$stratum <- as.vector(scale(YCG.data.strat$year, center = T, scale = T)) 


  ## load data for pooling all but with strata ####


# Here I am going to pool all years for the DF, but I am going to keep the strata in the years that have them. So to get the estimates for the stratified years, I will set all areas to 0 except for the year of interest.  For all non-stratified years, I will set 2011, 2014, 2016 areas to 0.  

#I think the easiest way to do this is to use the old data (which has the strata), but then to just rbind 2020 data on to it

# load the new data
load("./Output/Data/KSWS_MASTER.Rdata")

# save the 2020 region
region.table.20 <- full.region.table[7,1:2]

# save the 2020 sample.table
sample.table.20 <- sample.table[sample.table$Region.Label==2020,]

# save 2020 obs.table
obs.table.20 <- obs.table[obs.table$Region.Label==2020, ]

# save 2020 data
allData.20 <- allData[allData$stratum==2020, ]

## load old data (2010-2018) that have strata
load("./Output/Data/Archive/KSWS_MASTER.Rdata")

# add the 2020 data onto the different tables
full.region.table <- rbind(full.region.table, region.table.20)
sample.table.20$Region.Label <- as.factor(sample.table.20$Region.Label)
sample.table <- rbind(sample.table, sample.table.20)
obs.table.20$Region.Label <- as.factor(obs.table.20$Region.Label)
obs.table <- rbind(obs.table, obs.table.20)
allData.20$species <- as.factor(allData.20$species)
allData.20$stratum <- as.factor(allData.20$stratum)
allData.20$obs.time <- as.factor(allData.20$obs.time)
allData.20$obs.AMPM <- as.factor(allData.20$obs.AMPM)
allData.20$obs.observer <- as.factor(allData.20$obs.observer)
allData <- rbind(allData, allData.20)

# create region.tables for the different years. i.e. set areas to 0 for the years not of interest
region.table.11 <- full.region.table
region.table.11$Area[region.table.11$Region.Label != "2011_L" & region.table.11$Region.Label != "2011_H"] <- 0

region.table.14 <- full.region.table
region.table.14$Area[region.table.14$Region.Label != "2014_L" & region.table.14$Region.Label != "2014_H"] <- 0

region.table.16 <- full.region.table
region.table.16$Area[region.table.16$Region.Label != "2016_L" & region.table.16$Region.Label != "2016_H"] <- 0

region.table.nostrat <- full.region.table
region.table.nostrat$Area[region.table.nostrat$Region.Label != "2010" & region.table.nostrat$Region.Label != "2013" &
                          region.table.nostrat$Region.Label != "2018" & region.table.nostrat$Region.Label != "2020"]<-0

# subset data for YCG
YCG.data <- as.data.frame(allData[allData$species=="YCG",])

# create scaled continuous stratum variable for the detection functions
YCG.data$year <- ifelse(YCG.data$stratum=="2011_L"|YCG.data$stratum=="2011_H",2011,
                    ifelse(YCG.data$stratum=="2014_L"|YCG.data$stratum=="2014_H",2014,
                        ifelse(YCG.data$stratum=="2016_L"|YCG.data$stratum=="2016_H",2016,
                            ifelse(YCG.data$stratum=="2010", 2010,
                                ifelse(YCG.data$stratum=="2013",2013,
                                    ifelse(YCG.data$stratum=="2018", 2018,
                                        ifelse(YCG.data$stratum=="2020",2020,NA)))))))
YCG.data$stratum <- as.vector(scale(YCG.data$year, center = T, scale = T)) 

#
## Detection function models ####

# I will run two analyses - one that separates the stratified years and the unstratified years, both using pooling. I will need to run 3 models for the strata set because I will need to set the regions of the other strata to 0 to get a correct total estimate. The second analysis will pool all years (with strata being included) but I will set the areas to 0 for the different years to get estimates for the stratified and unstratified years separately.

# 1) split strata years and no strata years and analyse completely separately.
# 2) Pool all years, but for strata years set all other years area to 0, and then set strata years area to 0. 

# the best model from the CDS is HN strat, and so this will be used.

# set truncation
trunc.ycg <- 60

  ## 1) split years and analyse separately ####
    ## Strata ####

# 2011
ycg.strat.11 <- ds(YCG.data.strat, obs.table=obs.table.strat,sample.table=sample.table.strat,
                   region.table = region.table.11, truncation = trunc.ycg, key="hn",
                   formula=~stratum)

# 2014
ycg.strat.14 <- ds(YCG.data.strat, obs.table=obs.table.strat,sample.table=sample.table.strat,
                   region.table = region.table.14, truncation = trunc.ycg, key="hn",
                   formula=~stratum)

# 2016
ycg.strat.16 <- ds(YCG.data.strat, obs.table=obs.table.strat,sample.table=sample.table.strat,
                   region.table = region.table.16, truncation = trunc.ycg, key="hn",
                   formula=~stratum)

    ## No strata ####

ycg.nostrat <- ds(YCG.data.nostrat, obs.table = obs.table.nostrat, sample.table = sample.table.nostrat,
                  region.table = region.table.nostrat, truncation = trunc.ycg, key="hn",
                  formula=~stratum)


  ## 2) Pool all years ####

# Here I have pooled all years, but kept the strata, and just set appropriate areas to 0 to pull out estiamtes for stratified years and unstratified years.

# 2011
ycg.11 <- ds(YCG.data, obs.table=obs.table,sample.table=sample.table,
                   region.table = region.table.11, truncation = trunc.ycg, key="hn",
                   formula=~stratum)

# 2014
ycg.14 <- ds(YCG.data, obs.table=obs.table,sample.table=sample.table,
             region.table = region.table.14, truncation = trunc.ycg, key="hn",
             formula=~stratum)

# 2016
ycg.16 <- ds(YCG.data, obs.table=obs.table,sample.table=sample.table,
             region.table = region.table.16, truncation = trunc.ycg, key="hn",
             formula=~stratum)

# 2010, 2013, 2018, 2020
ycg.nostrat <- ds(YCG.data, obs.table=obs.table,sample.table=sample.table,
                  region.table = region.table.nostrat, truncation = trunc.ycg, key="hn",
                  formula=~stratum)


## Plotting ####
  ## results from analysis 1 ####

# extract estimates from above models
ycg.comp <- data.frame(Year = c("2010","2011","2013","2014","2016","2018","2020"),
                       Estimate = c(ycg.nostrat$dht$clusters$N$Estimate[1],
                                    ycg.strat.11$dht$clusters$N$Estimate[7],
                                    ycg.nostrat$dht$clusters$N$Estimate[3],
                                    ycg.strat.14$dht$clusters$N$Estimate[7],
                                    ycg.strat.16$dht$clusters$N$Estimate[7],
                                    ycg.nostrat$dht$clusters$N$Estimate[3],
                                    ycg.nostrat$dht$clusters$N$Estimate[4]),
                       lcl = c(ycg.nostrat$dht$clusters$N$lcl[1],
                                ycg.strat.11$dht$clusters$N$lcl[7],
                                ycg.nostrat$dht$clusters$N$lcl[3],
                                ycg.strat.14$dht$clusters$N$lcl[7],
                                ycg.strat.16$dht$clusters$N$lcl[7],
                                ycg.nostrat$dht$clusters$N$lcl[3],
                                ycg.nostrat$dht$clusters$N$lcl[4]),
                       ucl = c(ycg.nostrat$dht$clusters$N$ucl[1],
                               ycg.strat.11$dht$clusters$N$ucl[7],
                               ycg.nostrat$dht$clusters$N$ucl[3],
                               ycg.strat.14$dht$clusters$N$ucl[7],
                               ycg.strat.16$dht$clusters$N$ucl[7],
                               ycg.nostrat$dht$clusters$N$ucl[3],
                               ycg.nostrat$dht$clusters$N$ucl[4]),
                       cv = c(ycg.nostrat$dht$clusters$N$cv[1],
                              ycg.strat.11$dht$clusters$N$cv[7],
                              ycg.nostrat$dht$clusters$N$cv[3],
                              ycg.strat.14$dht$clusters$N$cv[7],
                              ycg.strat.16$dht$clusters$N$cv[7],
                              ycg.nostrat$dht$clusters$N$cv[3],
                              ycg.nostrat$dht$clusters$N$cv[4]),
                       Group = rep("strata", times=7))

ycg.comp

# load original YCG results
ycg.orig <- read.csv("Output/Results/YCG_results_final.csv", header=T)
ycg.orig$Year <- as.factor(ycg.orig$Year)
ycg.orig <- ycg.orig[ycg.orig$Label=="Grp",c(2,9,12:13,11)]
ycg.orig$Group <- "no strata"

ycg.comp <- rbind(ycg.comp,ycg.orig)


# plot together
ycg.plot.comp1 <- ggplot(ycg.comp, aes(x=Year, y=Estimate, group=Group, colour=Group))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))+
  ylim(0,1500)+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black"))

ggsave("Output/Results/Analysis_comparisons/YCG/ycg.plot.comp1.png", ycg.plot.comp1, height = 20, width = 20, 
       units = "cm", dpi=300)

# split by group
ggplot(ycg.comp, aes(x=Year, y=Estimate, colour=Group))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl))+
  facet_wrap(vars(Group), nrow=1)+
  ylim(0,1500)+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black"))

# predictably, the analysis with no strata (i.e. all data pooled) has produced more precise estimates. But when strata are taken into account, the estiamtes are lower in 2011, 2014, and 2016. 


  ## Results from analysis 2 ####

# extract results from models
ycg.comp2 <- data.frame(Year = c("2010","2011","2013","2014","2016","2018","2020"),
                       Estimate = c(ycg.nostrat$dht$clusters$N$Estimate[1],
                                    ycg.11$dht$clusters$N$Estimate[11],
                                    ycg.nostrat$dht$clusters$N$Estimate[4],
                                    ycg.14$dht$clusters$N$Estimate[11],
                                    ycg.16$dht$clusters$N$Estimate[11],
                                    ycg.nostrat$dht$clusters$N$Estimate[9],
                                    ycg.nostrat$dht$clusters$N$Estimate[10]),
                       lcl = c(ycg.nostrat$dht$clusters$N$lcl[1],
                               ycg.11$dht$clusters$N$lcl[11],
                               ycg.nostrat$dht$clusters$N$lcl[4],
                               ycg.14$dht$clusters$N$lcl[11],
                               ycg.16$dht$clusters$N$lcl[11],
                               ycg.nostrat$dht$clusters$N$lcl[9],
                               ycg.nostrat$dht$clusters$N$lcl[10]),
                       ucl = c(ycg.nostrat$dht$clusters$N$ucl[1],
                               ycg.11$dht$clusters$N$ucl[11],
                               ycg.nostrat$dht$clusters$N$ucl[4],
                               ycg.14$dht$clusters$N$ucl[11],
                               ycg.16$dht$clusters$N$ucl[11],
                               ycg.nostrat$dht$clusters$N$ucl[9],
                               ycg.nostrat$dht$clusters$N$ucl[10]),
                       cv = c(ycg.nostrat$dht$clusters$N$cv[1],
                              ycg.11$dht$clusters$N$cv[11],
                              ycg.nostrat$dht$clusters$N$cv[4],
                              ycg.14$dht$clusters$N$cv[11],
                              ycg.16$dht$clusters$N$cv[11],
                              ycg.nostrat$dht$clusters$N$cv[9],
                              ycg.nostrat$dht$clusters$N$cv[10]),
                       Group = rep("strata", times=7))

ycg.comp2 <- rbind(ycg.comp2,ycg.orig)


# plot together
ggplot(ycg.comp2, aes(x=Year, y=Estimate, group=Group, colour=Group))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))+
  ylim(0,1550)+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black"))


### look at the three analyses together

# first pull out the strata results from the first analysis
results1 <- ycg.comp[ycg.comp$Group=="strata",]
results1$Group <- "strata_split"

# re-name group in second analysis to differentiate from the first
ycg.comp2$Group <- ifelse(ycg.comp2$Group=="strata", "strata_pooled", ycg.comp2$Group)

# merge together
ycg.comp3 <- rbind(ycg.comp2, results1)

# plot all three
ycg.plot.comp.2 <- ggplot(ycg.comp3, aes(x=Year, y=Estimate, group=Group, colour=Group))+
  geom_point(position = position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position = position_dodge(width=0.3))+
  ylim(0,1550)+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line(colour = "black"))

ggsave("Output/Results/Analysis_comparisons/YCG/ycg.plot.comp2.png", ycg.plot.comp.2, height = 20, width = 20, 
       units = "cm", dpi=300)


#
##### CONCLUSIONS ############################################

#### conclusions for Gibbon & peafowl (first two analyses, up to line 2190)

# the positive overall conclusion is that the annual estimates from all of the different analytical approaches, for both species, show consistency in trends. In other words, the annual estimates from the different approaches are generally clumped together, and show the same overall pattern over the years. 

# For gibbons, there appears to be a wave-like pattern, which is hard to explain. I don't believe the populations actually fluctuate like that, and so this must be some artefact of the observation process, and thus the data. Perhaps there is an environmental factor that we are not accounting for that is affecting either behaviour, or observation. Perhaps something like how early/late the rains start, thus affecting foliage and therefore detection.  Something like that. It also could be an artefact of the priority strata, as there appears to be a pattern that the years with strata (2011, 2014, 2016) were generally higher estimates, but this doesn't hold true for 2016, which is a lower estimate. 

# The different approaches don't seem to reveal a strong pattern in terms of where the individual estimates sit within the group of estimates for that year. I.e. one approach isn't consistently high or consistently low relative to the other approaches. I think this is good, as it means that if we pick a particular approach, we won't necessarily be consistently over-estimating or under-estimating.

# when we fit linear models to the gibbon plot, we see that the approaches are split by pooling, rather than stratification.  i.e. the estimates noStratPool and stratPool show a similar trend, and the other approaches, which include stratfifed and non-stratified DF that ALL have year in as a covariate, are grouped. This suggests that including the stratification matters less than including year as a covariate. This is really good, as the stratification issue was one that has caused problems for years.  Turns out it's not as important as we thought.  This also means that for simplicity, we can continue with the analyses and ignore the stratification.  

# the CV plot for gibbons shows that the approach with no stratification and year in the DF as a continous covariate, produces the estimates with the lowest CV (4 out of 6 years). The next best is StratPool.

# the AIC comparison supports the above conclusion - the noStratyrCont approach as dAIC of 0.

# For Gren peafowl, the plot is particularly good for the years 2010 - 2016. All the estimates are tightly grouped together, showing strong agreement regardless of approach.  The 2018 estimate shows much more disagreement, although they all still agree that the population has increased. Interestingly, for GPF in 2018 it seems that stratification is the thing that is making the difference. The two estimates from the stratified approaches are the two that are significantly different from the rest. This can be seen further when you look at the plot with the linear models.

# The CV plot for peafowl is slightly more confusing. The approach with stratification and pooled DF has the most number of estimates with the lowest CV, whereas the approach with stratification and year as a covariate has the most number of estimates with the highest CV. But the AIC table shows that the approaches with stratification are significantly worse that the approaches with no stratification (~6000 comapared to ~1000). The appraoch with no stratification and year in the DF as a continuous covariate again has the dAIC of 0.

## Therefore the overall decision I have made is that:
# 1) I am removing the priority strata from all analyses
# 2) I will analyse BSD on an annual basis, as they have sufficient within-year observations. RMJ could probably be done annually although this is currently an objective assessment. A threshold should probably be agreed upon at some point
# 3) For all other species, I will use a POOLED DF with year as a continuous covariate. This will allow us to make infrences about changes in detectability over time.


### Conclusions for BSD & YCG second analyses - lines 2190 onwards

# In my discussion with Eric Rextad, he suggested ignoring the strata for the sake of simplicity.  This was fine, but for three species (BSD, YCG, PTM) the estiamtes were a bit odd, particulalry 2011 and 2014 which meant that there were quite unrealistic differences in estimates between years (e.g. between 2010 and 2011).  Eric did say that by ignoring the strata we would probably end up overestimating slightly in the years with strata.  For these species this was clearly happening quite badly. Therefore I tested two different approaches and compared the results with the original results

# You can find the plots in Outputs/Results/Analysis_comparison/YCG.  I didn't get round to saving the BSD plots, but the code above can easily be re-run. 

# The conclusion is not necessarily straightforward. Firstly, in 2018 and 2020 the different approaches don't seem to make much difference.  But they do in other years.  When strata are taken into account, but data from all years are pooled (strata_pooled), then the estimates for 2010 and 2013 are similar to the original analysis (no strata).  When the stratified years are removed altogeher, so that un-stratified years are analysed separately (strata_split) then 2010 and 2013 are quite a bit higher.  This suggests that the presence of the stratified years is pulling those estimates down for some reason.  When strata are taken into account in any way (strat_pooled and strat_split), the stratified years are always lower than when strata are ignored completely (no strata).  I would suspect that this is to do with weighting, but Olly asked this question on the google group and apparently the R package does not do weighting when combining estimates from different strata.  Therefore the only reason I can think is because when you split the H and L strata and estimate their abundance separately, the L strata is producing a very low density estimate (which is "correct"), and the H strata is producing a large estimate. These two estimates are then meaned (D) and summed (N). But when the strata are ignored, the data from the high density strata are swamping the data from the low strata, inflating the density within the low density strata.  Intuitively I think that analysing the strata separately makes more sense for these species, as that is exactly the purpose of having the strata - we think the that density is different between the strata.

# Therefore for YCG, PTM, and BSD I analysed 2011, 2014, 2016 separately from the other years to account for the strata.