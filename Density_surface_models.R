# This R script is for the density surface modelling of all key species in Keo Seima Wildlife Sanctuary in Cambodia. The analysis was done by Matt Nuttall (m.n.nuttal1@stir.ac.uk / mattnuttall00@gmail.com) and Olly Griffin (ogriffin@wcs.org).  The data belongs to the Wildlife Conservation Society Cambodia Program. For information on access to the data, see the blurb at the top of the conventional distance sampling (point estimate) script.  


# All of the detection functions are taken from the CDS analysis, as that analysis was where more effort and time was taken to determine the best DF model for each species. Therefore the DF model selection code under each of the species is obsolete, but has been left in for refrence purposes.

# The analysis was done for 2010-2018, but then in July 2020 the 2020 transect data was added, and the analysis re-done. 

# T20 was removed from all data (see explanation in the CDS script).

# The below script conducts the analysis for all years (i.e. uses all data from all years to fit models, and then predicts relative abundance for each year using the correct habitat/deforestation layer for that year). In the paper in Conservation Science and Practice however, we only present the results from 2020. 

# For guidance on running DSMs, see this really useful example:
# http://distancesampling.org/R/vignettes/mexico-analysis.html#fn2

# and this really useful presentation:
# https://bit.ly/2HwmCFy

# In the presentation, David Miller suggests using a thin plate regression spline with modified smoothing penalties (bs="ts"). He also suggests the following general strategy for selecting smooths terms, which I have broadly followed:

# 1) Build a model with the smooths you want
# 2) Make sure that the smooths are flexible enougth (k=...)
# 3) Remove smoothes that have been shrunk
# 4) Remove non-signifcant smooths


#### Load packages, functions, graphics settings ####

library("mgcv")
library("dsm")
library("tidyverse")
library("plyr")
library("Distance")
library("shapefiles")
library("rgdal")
library("maptools")
library("knitr")
library('patchwork')
library('broom')
library('viridis')
source("grid_plot_obj.R") # custom function used in plots later. See separate script

# Define some settings we'll use for ggplot2 later on
gg.opts <- theme(panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())



#### functions


## Plotting

# create abundance plotting spatial dataframe
abunPlotDF <- function(dat,predgrid){
  
row.names(dat) <- row.names(predgrid)
spdf_10.core <- SpatialPolygonsDataFrame(predgrid, dat)
spdf_10.core@data$id <- rownames(spdf_10.core@data)
spdf.points_10.core <- fortify(spdf_10.core) 
spdf.df_10.core <- join(spdf.points_10.core, spdf_10.core@data, by="id")
spdf.df_10.core$x <- spdf.df_10.core$long
spdf.df_10.core$y <- spdf.df_10.core$lat
return(spdf.df_10.core)
  
}

# greyscale continuous. The "species" term here just refers to the title of the plot so it can be any character string.
GSplotFun <- function(dat,datshp,fill, species){
  
  ggplot()+
    geom_polygon(aes_string(x="x",y="y",fill=fill, group="group"), 
    data=dat)+
    gg.opts+
    geom_path(aes(x=x, y=y),data=datshp)+
    labs(fill=fill)+
    scale_fill_gradient2(low = "white", high = "black")+
    coord_equal()+
    ggtitle(species)
}

# colour (viridis) continuous
CLplotFun <- function(dat,datshp,fill){
  
  ggplot()+
    geom_polygon(aes_string(x="x",y="y",fill=fill, group="group"), 
    data=dat) +
    gg.opts+
    geom_path(aes(x=x, y=y),data=datshp)+
    labs(fill=fill)+
    scale_fill_viridis()+
    coord_equal()
  
}

# greyscale, discrete bins. fill should be either"group2" or "CV", label should be either "Abundance" or "Variance", and legend should be either "Relative abundance" or "CV"
GSplotBin <- function(dat,fill,datshp,label,legend){
  
  ggplot()+
    geom_polygon(aes_string(x="x",y="y",fill=fill, group="group"), 
    data=dat) +
    gg.opts+
    geom_path(aes(x=x, y=y),data=survey.area.core)+
    labs(fill=label)+
    scale_fill_brewer(legend, palette = "Greys", direction=-1, 
                      guide="legend", drop=FALSE)+
    coord_equal()  
  
}
  
# save plot function
saveplot <- function(plot,filepath){
  ggsave(filepath, plot, width = 30, height = 20, units = "cm", dpi = 300)
}



## Variance functions

# estimate variance
varEstfun <- function(dat,model){
  
  # create block indexes of 50 entries 
chunk_start <- c(seq(1, nrow(dat)+1, by=50))
chunk_end <- c(chunk_start[2:length(chunk_start)]-1, nrow(dat))

# storage
vars <- rep(NA, nrow(dat))

# loop over the index
for(ii in seq_along(chunk_start)){

   # get the current chunk
   this_chunk <- dat[chunk_start[ii]:chunk_end[ii], ]
   # split
   pred.new <- split(this_chunk, 1:nrow(this_chunk))
   # estimate variance
   vargam <- dsm.var.gam(model, pred.new, off.set=40000, type.pred ="response")
   # extract the cell variances and predictions
   vars[chunk_start[ii]:chunk_end[ii]] <- vargam$pred.var
}

# Create new dataframe for plotting
df <- data.frame(id = 1:47801,
                         variance = vars)
return(df)

}

# create variance plotting spatial dataframe
varPlotDF <- function(dat,predgrid){
  
  row.names(dat) <- row.names(predgrid)
  var.spdf.core <- SpatialPolygonsDataFrame(predgrid, dat)
  var.spdf.core@data$id <- rownames(var.spdf.core@data)
  var.spdf.points <- fortify(var.spdf.core) 
  var.spdf.df <- join(var.spdf.points, var.spdf.core@data, by="id")
  var.spdf.df$x <- var.spdf.df$long
  var.spdf.df$y <- var.spdf.df$lat
  return(var.spdf.df)

}

# add CV to variance spdf (the column will be named "group2" just so that the above plotting function works for abundance and for variance)
CVaddFun <- function(Ndat, Vdat){

  Vdat$abundance <- Ndat$abundance
  Vdat$group2 <- sqrt(Vdat$variance)/Vdat$abundance
  return(Vdat)

}



## Binning functions

# put abundance into bins. "dat" needs to be the abundance spdf's
binFunAbun <- function(dat){
  
  Nquants25 <- quantile(dat$abundance, probs = 0.25)
  Nquants50 <- quantile(dat$abundance, probs = 0.5)
  Nquants75 <- quantile(dat$abundance, probs = 0.75)
  dat$group2 <- ifelse(dat$abundance < Nquants25, "Very low",
                       ifelse(dat$abundance < Nquants50, "Low", 
                              ifelse(dat$abundance < Nquants75, "Medium",
                                     ifelse(dat$abundance >= Nquants75, "High",NA)))) 
  dat$group2 <- factor(dat$group2, levels = c("High","Medium","Low","Very low"))
  return(dat)
  
}

# put variance into quartile bins. "dat" needs to be the variance spdf's that have had CV calculated already using CVaddFUn above
binFunVar <- function(dat){
  
  Nquants25 <- quantile(dat$group2, probs = 0.25)
  Nquants50 <- quantile(dat$group2, probs = 0.5)
  Nquants75 <- quantile(dat$group2, probs = 0.75)
  dat$group2 <- ifelse(dat$group2 < Nquants25, "Very low",
                       ifelse(dat$group2 < Nquants50, "Low", 
                              ifelse(dat$group2 < Nquants75, "Medium",
                                     ifelse(dat$group2 >= Nquants75, "High",NA)))) 
  dat$group2 <- factor(dat$group2, levels = c("High","Medium","Low","Very low"))
  return(dat)
  
}

# Put CV into discrete bins (currently 0-10%, 11-20%, 21-30%, 31-40%, 41-50%, 51-60%, >60%)
binFunVar2 <- function(dat){
  
  dat$CV <- ifelse(dat$CV < 0.11, "< 10%",
                       ifelse(dat$CV < 0.21, "11-20%", 
                              ifelse(dat$CV < 0.31, "21-30%",
                                     ifelse(dat$CV < 0.41, "31-40%",
                                            ifelse(dat$CV < 0.51, "41-50%",
                                                   ifelse(dat$CV < 0.61, "51-60%",
                                                          ifelse(dat$CV > 0.61, "> 60%", NA ))))))) 
  dat$CV <- factor(dat$CV, levels = c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
  return(dat)
  
}


#### Load preddata & segdata ####

# Load segdata (most recent segdata has unique transect IDs for each multi-day visit)
segdata <- read.csv("segdata.csv", header = TRUE)
segdata$Sample.Label <- as.factor(segdata$Sample.Label)
segdata$Transect.Label <- as.factor(segdata$Transect.Label)
segdata$Transect <- as.factor(segdata$Transect)
segdata <- segdata %>% select(-dstELC)
str(segdata)
head(segdata)


# Annual preddatas for the core area only
preddata10_core <- read.csv("Preddata/Annual_preddatas/core_only/preddata10_core.csv")
preddata11_core <- read.csv("Preddata/Annual_preddatas/core_only/preddata11_core.csv")
preddata13_core <- read.csv("Preddata/Annual_preddatas/core_only/preddata13_core.csv")
preddata14_core <- read.csv("Preddata/Annual_preddatas/core_only/preddata14_core.csv")
preddata16_core <- read.csv("Preddata/Annual_preddatas/core_only/preddata16_core.csv")
preddata18_core <- read.csv("Preddata/Annual_preddatas/core_only/preddata18_core.csv")
preddata20_core <- read.csv("Preddata/Annual_preddatas/core_only/preddata20_core.csv")


#### Load survey area + prediction grid shapefiles #### 

## Load survey area (core zone)
survey.area.core <- readOGR(dsn = "Spatial_data/ksws/spf_core.shp")

# Simplify the object (core zone)
survey.area.core <- data.frame(survey.area.core@polygons[[1]]@Polygons[[1]]@coords)
names(survey.area.core) <- c("x", "y")


# load core only prediction grid. 
pred.polys_200 <- readOGR(dsn = "Spatial_data/grid/grid_200_200_coreClip.shp")


#### Green Peafowl ############################################################
## Load data ####

# Observation data. Unique to species 
obsdata <-read.csv("Species_Data/GPF/R Data/obsdata.csv", header = TRUE)
obsdata$object <- as.factor(obsdata$object)
obsdata$Sample.Label <- as.factor(obsdata$Sample.Label)
str(obsdata)
head(obsdata)

# Transect data. Unique to species
distdata <- read.csv("Species_Data/GPF/R Data/distdata.csv", header = TRUE)
distdata$object <- as.factor(distdata$object)
distdata$NameObserver <- as.factor(distdata$NameObserver)
distdata$transect <- as.factor(distdata$transect)
distdata$date <- as.Date(distdata$date, format = "%d/%m/%Y")
distdata$stratum <- as.vector(scale(distdata$year, center = T, scale = T)) 
distdata$year <- as.factor(distdata$year)

str(distdata)
head(distdata)

## Plot the covariates across the grid, with group sizes #### 

# Warning - plots take a few minutes to run

## habitat
plot_obs_habitat <- ggplot() + grid_plot_obj(preddata200$habitat, "habitat", pred.polys200) + coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=distdata, colour="red", alpha=I(0.7))+
                    gg.opts
#ggsave("plot_obs_habitat.png", plot_obs_habitat, width = 20, height = 20, units = "cm", dpi = 300)

# dstWater
plot_obs_dstWater <- ggplot() + grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys200) + coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), data=distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("plot_obs_dstWater.png", plot_obs_dstWater, width = 20, height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_obs_dstStlmnt <- ggplot() + grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys200) + coord_equal()+
                      labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=distdata, colour="red", alpha=I(0.7))+
                      gg.opts
ggsave("plot_obs_dstStlmnt.png", plot_obs_dstStlmnt, width = 20, height = 20, units = "cm", dpi = 300)

# dstRoad
plot_obs_dstRoad <- ggplot() + grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys200) + coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=distdata, colour="red", alpha=I(0.7))+
                    gg.opts

ggsave("plot_obs_dstRoad.png", plot_obs_dstRoad, width = 20, height = 20, units = "cm", dpi = 300)

# dstBorder
plot_obs_dstBorder <- ggplot() + grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys200) + coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=distdata, colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("plot_obs_dstBorder.png", plot_obs_dstBorder, width = 20, height = 20, units = "cm", dpi = 300)

# dstStation
plot_obs_dstStation <- ggplot() + grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys200) + coord_equal()+
                       labs(fill="Distance to ranger station",x="x",y="y",size="Group size")+
                       geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                       geom_point(aes(x, y, size=size), data=distdata, colour="red", alpha=I(0.7))+
                       gg.opts

ggsave("plot_obs_dstStation.png", plot_obs_dstStation, width = 20, height = 20, units = "cm", dpi = 300)

# dstELC
plot_obs_dstELC <- ggplot() + grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys200) + coord_equal()+
                   labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                   geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                   geom_point(aes(x, y, size=size), data=distdata, colour="red", alpha=I(0.7))+
                   gg.opts

ggsave("plot_obs_dstELC.png", plot_obs_dstELC, width = 20, height = 20, units = "cm", dpi = 300)

# elevation
plot_obs_elev <- ggplot() + grid_plot_obj(preddata200$elevation, "elevation", pred.polys200) + coord_equal()+
                 labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                 geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                 geom_point(aes(x, y, size=size), data=distdata, colour="red", alpha=I(0.7))+
                 gg.opts

ggsave("plot_obs_elev.png", plot_obs_elev, width = 20, height = 20, units = "cm", dpi = 300)


## Exploratory plots & linear models ####


## Histograms

# Distance 
h1 <- ggplot(distdata, aes(distance))+ geom_histogram(binwidth = 1)
h2 <- ggplot(distdata, aes(distance))+ geom_histogram(binwidth = 5)
h3 <- ggplot(distdata, aes(distance))+ geom_histogram(binwidth = 10)
h4 <- ggplot(distdata, aes(distance))+ geom_histogram(binwidth = 15)
h5 <- ggplot(distdata, aes(distance))+ geom_histogram(binwidth = 20)
h6 <- ggplot(distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(h1,h2,h3,h4,h5,h6)

# cluster size, observer, habitat, year, month, transect
h7 <- ggplot(distdata, aes(size))+geom_histogram(binwidth = 0.5)
h8 <- ggplot(distdata, aes(NameObserver))+geom_histogram(stat="count")
h9 <- ggplot(distdata, aes(habitat))+geom_histogram(stat="count")
h10 <- ggplot(distdata, aes(year))+geom_histogram(stat="count")
h11 <- ggplot(distdata, aes(month))+geom_histogram(stat="count")
h12 <- ggplot(distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(h7,h8,h9,h10,h11,h12)

## Plots of distance against variables
plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

d1 <- ggplot(distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
d2 <- ggplot(distdata, aes(x=observer, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
d3 <- ggplot(distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")
d4 <- ggplot(distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
d5 <- ggplot(distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(d1,d2,d3,d4,d5)

## Plots of cluster size against variables
s1 <- ggplot(distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
s2 <- ggplot(distdata, aes(x=observer, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
s3 <- ggplot(distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")
s4 <- ggplot(distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
s5 <- ggplot(distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(s1,s2,s3,s4,s5)

## Linear models

# cluster size ~ distance
newdist <- data.frame(distance=seq(0,133,len=10))
lm1 <- lm(size~distance, data=distdata)
plot(distdata$size~distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# There is very little support for a relationship between cluster size and distance. The effect size is very small (0.008) and the model has a p-value of 0.19. This suggests that we will not need to worry about size bias in the detection function

## Estimating detection function ####

# The detection function model will now be taken from the CDS anaysis where more time and effort was taken in assessing the DF model fit. The DSMs and the CDS use the same data, and therefore I will no longer use the DF model selection below, but will use the model decided upon during the CDS analysis.  See the CDS script for the details of the model selection.   

# the DF model used in the CDS is hn.strat.size, and so I run that model below. 
gpfDF.hn.strat.size <- ds(distdata, truncation=50, key="hn", formula=~stratum+size)

# for old DF model selection see the DSM_obsolete_code script

 
## Fitting a spatial model #### 

# the DF model to be used is gpfDF.hn.strat.size

# Will build models according to parsimony e.g. starting complex then getting simpler until all remaining variables are important

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasi-possion, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated. This is because I don't think the group size estimates for this species are accurate enough, and the species is gregarious

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
obsdata <- obsdata %>% filter(distance <= 50)

# Comments about run time of models is only relevant when unique transect IDs are being used



   ## Quasipoisson response ####

# Below are the quasipoisson models 


# Start with a saturated model
gpfDSM.sat1 <- dsm(Nhat ~ s(dstWater, bs="ts")+ s(dstStlmnt, bs="ts")+ 
                          s(dstBorder, bs="ts")+ s(dstStation, bs="ts")+
                          s(elevation, bs="ts")+ 
                          habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = T)
summary(gpfDSM.sat1)
gam.check(gpfDSM.sat1)
par(mfrow=c(3,2))
plot(gpfDSM.sat1, scale=0)
# DE = 21.3. Only dstBorder and elevation sig

# I will increase k for the above model
gpfDSM.sat2 <- dsm(Nhat ~ s(dstWater, bs="ts",k=15)+ s(dstStlmnt, bs="ts",k=15)+ 
                          s(dstBorder, bs="ts",k=15)+ s(elevation, bs="ts",k=15)+
                          s(dstStation, bs="ts",k=15)+ habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.sat2)
gam.check(gpfDSM.sat2)
dev.off()
par(mfrow=c(3,2))
plot(gpfDSM.sat2,scale=0)
# R-sq is 0.027, DE is 25.4. dstStlmnt now the only sig term. Increasing k is not the right thing to do based on the smooth plots - overfitting

# anova for the two saturated models
anova(gpfDSM.sat1,gpfDSM.sat2, test = "Chisq")
# gpfDSM.sat2 is the better model, lower residual deviance, but this is to be expected as the model is overfit

# Reduce k from sat model
gpfDSM.sat3 <- dsm(Nhat ~ s(dstWater, bs="ts",k=4)+ s(dstStlmnt, bs="ts",k=4)+ 
                          s(dstBorder, bs="ts",k=4)+ s(elevation, bs="ts",k=4)+ 
                          s(dstStation, bs="ts",k=4)+ habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.sat3)
gam.check(gpfDSM.sat3)
dev.off()
par(mfrow=c(3,2))
plot(gpfDSM.sat3,scale=0)
# DE = 14.8, r2 = 0.01. dstWater, dstBorder, and dstStation sig. dstStlmnt looks like it could be linear.

# make dstStlmnt linear
gpfDSM.sat4 <- dsm(Nhat ~ s(dstWater, bs="ts",k=4)+ s(dstBorder, bs="ts",k=4)+
                      s(elevation, bs="ts",k=4)+ s(dstStation, bs="ts",k=4)+
                      habitat + dstStlmnt,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.sat4)
gam.check(gpfDSM.sat4)
dev.off()
par(mfrow=c(2,2))
plot(gpfDSM.sat4,scale=0)
# DE = 14.8, r2 = 0.01. dststlmtnt not sig as linear term. gam.check suggests k too low for dstborder.

# remove dstStlmnt and increase k
gpfDSM.sat5 <- dsm(Nhat ~ s(dstWater, bs="ts",k=5)+ s(dstBorder, bs="ts",k=5)+
                     s(elevation, bs="ts",k=5)+ s(dstStation, bs="ts",k=5)+
                     habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.sat5)
gam.check(gpfDSM.sat5)
dev.off()
par(mfrow=c(2,2))
plot(gpfDSM.sat5,scale=0)
# DE = 15.1, r2 = 0.009.  Plots look exactly the same, but sig has gone disappeard for dstWater and reduced for dstStation. gam.check still suggests k is too low, but if I increase it further I will get overfitting

# back to sat4 but remove elevation and dstStlmnt
gpfDSM.sat6 <- dsm(Nhat ~ s(dstWater, bs="ts",k=4)+ s(dstBorder, bs="ts",k=4)+
                    s(dstStation, bs="ts",k=4)+ habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.sat6)
gam.check(gpfDSM.sat6)
dev.off()
par(mfrow=c(2,2))
plot(gpfDSM.sat6,scale=0)
# DE = 14.2, r2= 0.009. all terms sig (dstWater only just). Not overfitted


## Best QP model is gpfDSM.sat6


   ## Tweedie response ####

# saturated model
gpfDSM.tw.sat1 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                             s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ 
                             s(elevation,bs="ts")+ habitat,
                  gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.tw.sat1)
gam.check(gpfDSM.tw.sat1)
par(mfrow=c(3,2))
plot(gpfDSM.tw.sat1, scale=0)
# AIC= 13970, DE=13.4. dstborder sig.

# increase k except dstBorder
gpfDSM.tw.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=10)+ s(dstStlmnt,bs="ts",k=10)+ 
                        s(dstBorder,bs="ts")+ s(dstStation,bs="ts",k=10)+ 
                        s(elevation,bs="ts",k=10)+ habitat,
                      gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                      family = tw(), engine = "gam",
                      segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.tw.sat2)
gam.check(gpfDSM.tw.sat2)
par(mfrow=c(3,2))
plot(gpfDSM.tw.sat2, scale=0)
# 13970, DE=13.4.  Increasing k made no difference.

# go back to sat1 and remove dstStation (lowest EDF)
gpfDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                        s(dstBorder,bs="ts")+  s(elevation,bs="ts")+ habitat,
                      gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                      family = tw(), engine = "gam",
                      segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.tw.3)
gam.check(gpfDSM.tw.3)
par(mfrow=c(3,2))
plot(gpfDSM.tw.3, scale=0)
# AIC= 13970, DE=13.4. dstborder sig.

# remove elevation
gpfDSM.tw.4 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                     s(dstBorder,bs="ts")+ habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.tw.4)
gam.check(gpfDSM.tw.4)
par(mfrow=c(2,2))
plot(gpfDSM.tw.4, scale=0)
# AIC= 13970, DE=13.4. dstborder sig.

# remove dstStlmnt
gpfDSM.tw.5 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstBorder,bs="ts")+ habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.tw.5)
gam.check(gpfDSM.tw.5)
par(mfrow=c(2,1))
plot(gpfDSM.tw.5, scale=0)
# AIC= 13969, DE=13.3. dstborder sig. AIC gone down a smidge.

# try dstWater as linear term
gpfDSM.tw.6 <- dsm(Nhat ~ s(dstBorder,bs="ts")+ habitat + dstWater,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.tw.6)
gam.check(gpfDSM.tw.6)
par(mfrow=c(2,1))
plot(gpfDSM.tw.6, scale=0)
# AIC= 13970, DE=13.4. Model slightly worse than above. dstborder sig. 


# remove dstWater
gpfDSM.tw.7 <- dsm(Nhat ~ s(dstBorder,bs="ts")+ habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.tw.7)
gam.check(gpfDSM.tw.7)
par(mfrow=c(2,1))
plot(gpfDSM.tw.7, scale=0)
# AIC= 13969, DE=13.2. Best model by AIC


# The best TW model is gpfDSM.tw.7


   ## Negative binomial response ####

# start with a fully saturated model with shrinkage 
gpfDSM.nb.sat1 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                             s(dstBorder,bs="ts")+ s(elevation,bs="ts")+
                             s(dstStation,bs="ts")+ 
                             habitat,
                  gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.nb.sat1)
dev.off()
par(mfrow=c(3,2))
gam.check(gpfDSM.nb.sat1)
plot(gpfDSM.nb.sat1,scale=0)
# AIC = 1387, DE = 24.3. dstborder and dstWater (just) sig. gam.check suggest k possibly too low

# as above but incrased k
gpfDSM.nb.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts", k=10)+ s(dstStlmnt,bs="ts", k=10)+ 
                        s(dstBorder,bs="ts", k=10)+ s(elevation,bs="ts", k=10)+
                        s(dstStation,bs="ts", k=10)+ 
                        habitat,
                      gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                      family = nb(), engine = "gam",
                      segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.nb.sat2)
dev.off()
par(mfrow=c(3,2))
gam.check(gpfDSM.nb.sat2)
plot(gpfDSM.nb.sat2,scale=0)
# AIC = 1387, DE = 24.3. no change

# remove elevation
gpfDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts", k=10)+ s(dstStlmnt,bs="ts", k=10)+ 
                        s(dstBorder,bs="ts", k=10)+ s(dstStation,bs="ts", k=10)+ 
                        habitat,
                      gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                      family = nb(), engine = "gam",
                      segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.nb.3)
dev.off()
par(mfrow=c(2,2))
gam.check(gpfDSM.nb.3)
plot(gpfDSM.nb.3,scale=0)
# AIC = 1387, DE = 24.3. no change


# remove dstStation
gpfDSM.nb.4 <- dsm(Nhat ~ s(dstWater,bs="ts", k=10)+ s(dstStlmnt,bs="ts", k=10)+ 
                     s(dstBorder,bs="ts", k=10)+  
                     habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.nb.4)
dev.off()
par(mfrow=c(2,2))
gam.check(gpfDSM.nb.4)
plot(gpfDSM.nb.4,scale=0)
# AIC = 1387, DE = 24.3. no change. dstWater potentially linear

# remove dstStlmnt
gpfDSM.nb.5 <- dsm(Nhat ~ s(dstWater,bs="ts", k=10)+  s(dstBorder,bs="ts", k=10)+  
                     habitat,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.nb.5)
dev.off()
par(mfrow=c(2,1))
gam.check(gpfDSM.nb.5)
plot(gpfDSM.nb.5,scale=0)
# AIC = 1387, DE = 24.3. no change. dstWater potentially linear

# dstWater as linear term
gpfDSM.nb.6 <- dsm(Nhat ~  s(dstBorder,bs="ts", k=10)+  
                     habitat + dstWater,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.nb.6)
dev.off()
par(mfrow=c(1,1))
gam.check(gpfDSM.nb.6)
plot(gpfDSM.nb.6,scale=0, xlab="Distance to border (m)")
# AIC = 1387, DE = 24.4. Very slight increase in DE. dstWater still jsut sig as linear term. gam.check suggests k may be too low

# increase k
gpfDSM.nb.7 <- dsm(Nhat ~  s(dstBorder,bs="ts", k=15)+  
                     habitat + dstWater,
                   gpfDF.hn.strat.size, segdata, obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group=T)
summary(gpfDSM.nb.7)
dev.off()
par(mfrow=c(1,1))
gam.check(gpfDSM.nb.7)
plot(gpfDSM.nb.7,scale=0)
# AIC = 1387, DE = 24.4. increasing k makes no difference.



### gpfDSM.nb.6 is the best NB model

   ## Model selection ####

# The best quasipoisson model is gpfDSM.sat6
# The best tweedie model is gpfDSM.tw.7
# The best negative binomial model is gpfDSM.nb.6

summary(gpfDSM.sat6) # DE = 14.2, 
summary(gpfDSM.tw.7) # DE = 13.2, AIC = 13969
summary(gpfDSM.nb.6) # DE = 24.4,   AIC = 1387
# the NB model explains the most deviance, and has much lower AIC than the tw() 

# Test the QP versus NB using anova
anova(gpfDSM.sat6,gpfDSM.nb.6, test = "Chisq")
# Both models have similar residual deviance

# Checking the Q-Q plots. 
par(mfrow=c(2,2))
gam.check(gpfDSM.sat6)
gam.check(gpfDSM.tw.7)
gam.check(gpfDSM.nb.6)
# NB qq plot is the best, and the response versus fitted values plot is the nicest

# gpfDSM.nb.6 is the best model.

## Autocorrelation ####

# Checking autocorrelation for  gpfDSM.nb.6

par(mfrow=c(1,1))
dsm.cor(gpfDSM.nb.6, max.lag=10, Segment.Label="Sample.Label")

# very slight autocorrelation at lag 5, but not enough to worry about


## Abundance estimation ####


# Annual predictions using negative binomial final model
gpf.Final.pred10.core <- predict(gpfDSM.nb.6, preddata10_core, off.set = 40000)
write.csv(gpf.Final.pred10.core, file="Results/GPF/core_only/gpf.pred10.core.csv")

gpf.Final.pred11.core <- predict(gpfDSM.nb.6, preddata11_core, off.set = 40000)
write.csv(gpf.Final.pred11.core, file="Results/GPF/core_only/gpf.pred11.core.csv")

gpf.Final.pred13.core <- predict(gpfDSM.nb.6, preddata13_core, off.set = 40000)
write.csv(gpf.Final.pred13.core, file="Results/GPF/core_only/gpf.pred13.core.csv")

gpf.Final.pred14.core <- predict(gpfDSM.nb.6, preddata14_core, off.set = 40000)
write.csv(gpf.Final.pred14.core, file="Results/GPF/core_only/gpf.pred14.core.csv")

gpf.Final.pred16.core <- predict(gpfDSM.nb.6, preddata16_core, off.set = 40000)
write.csv(gpf.Final.pred16.core, file="Results/GPF/core_only/gpf.pred16.core.csv")

gpf.Final.pred18.core <- predict(gpfDSM.nb.6, preddata18_core, off.set = 40000)
write.csv(gpf.Final.pred18.core, file="Results/GPF/core_only/gpf.pred18.core.csv")

gpf.Final.pred20.core <- predict(gpfDSM.nb.6, preddata20_core, off.set = 40000)
write.csv(gpf.Final.pred20.core, file="Results/GPF/core_only/gpf.pred20.core.csv")


# I want to make my own ggplot, as they will be easier to manipulate, rather than using their function.  First I put the model prediction results for each year into a dataframe
gpf.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = gpf.Final.pred10.core)

gpf.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = gpf.Final.pred11.core)

gpf.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = gpf.Final.pred13.core)

gpf.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = gpf.Final.pred14.core)

gpf.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = gpf.Final.pred16.core)

gpf.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = gpf.Final.pred18.core)

gpf.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = gpf.Final.pred20.core)



## This creates a dataframe that can be plotted as a map
gpf.spdf.df_10.core <- abunPlotDF(gpf.df.Final10.core, pred.polys_200)
gpf.spdf.df_11.core <- abunPlotDF(gpf.df.Final11.core, pred.polys_200)
gpf.spdf.df_13.core <- abunPlotDF(gpf.df.Final13.core, pred.polys_200)
gpf.spdf.df_14.core <- abunPlotDF(gpf.df.Final14.core, pred.polys_200)
gpf.spdf.df_16.core <- abunPlotDF(gpf.df.Final16.core, pred.polys_200)
gpf.spdf.df_18.core <- abunPlotDF(gpf.df.Final18.core, pred.polys_200)
gpf.spdf.df_20.core <- abunPlotDF(gpf.df.Final20.core, pred.polys_200)


# save SPDFs
write.csv(gpf.spdf.df_10.core,file="Results/GPF/Plots/core_only/spdf/gpf.spdf.df_10.core.csv")
write.csv(gpf.spdf.df_11.core,file="Results/GPF/Plots/core_only/spdf/gpf.spdf.df_11.core.csv")
write.csv(gpf.spdf.df_13.core,file="Results/GPF/Plots/core_only/spdf/gpf.spdf.df_13.core.csv")
write.csv(gpf.spdf.df_14.core,file="Results/GPF/Plots/core_only/spdf/gpf.spdf.df_14.core.csv")
write.csv(gpf.spdf.df_16.core,file="Results/GPF/Plots/core_only/spdf/gpf.spdf.df_16.core.csv")
write.csv(gpf.spdf.df_18.core,file="Results/GPF/Plots/core_only/spdf/gpf.spdf.df_18.core.csv")
write.csv(gpf.spdf.df_20.core,file="Results/GPF/Plots/core_only/spdf/gpf.spdf.df_20.core.csv")



    ## Plotting continuous ####

# Load spatial dataframes
spdf.df_10.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_10.core.csv")
spdf.df_11.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_11.core.csv")
spdf.df_13.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_13.core.csv")
spdf.df_14.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_14.core.csv")
spdf.df_16.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_16.core.csv")
spdf.df_18.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_18.core.csv")
spdf.df_20.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_20.core.csv")


# greyscale plots
GPF_plot_10_core <- GSplotFun(spdf.df_10.core, survey.area.core, "abundance", "2010")
GPF_plot_11_core <- GSplotFun(spdf.df_11.core, survey.area.core, "abundance", "2011")
GPF_plot_13_core <- GSplotFun(spdf.df_13.core, survey.area.core, "abundance", "2013")
GPF_plot_14_core <- GSplotFun(spdf.df_14.core, survey.area.core, "abundance", "2014")
GPF_plot_16_core <- GSplotFun(spdf.df_16.core, survey.area.core, "abundance", "2016")
GPF_plot_18_core <- GSplotFun(spdf.df_18.core, survey.area.core, "abundance", "2018")
GPF_plot_20_core <- GSplotFun(spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(GPF_plot_10_core,"Results/GPF/Plots/core_only/greyscale/GPF_plot_10_core.png")
saveplot(GPF_plot_11_core,"Results/GPF/Plots/core_only/greyscale/GPF_plot_11_core.png")
saveplot(GPF_plot_13_core,"Results/GPF/Plots/core_only/greyscale/GPF_plot_13_core.png")
saveplot(GPF_plot_14_core,"Results/GPF/Plots/core_only/greyscale/GPF_plot_14_core.png")
saveplot(GPF_plot_16_core,"Results/GPF/Plots/core_only/greyscale/GPF_plot_16_core.png")
saveplot(GPF_plot_18_core,"Results/GPF/Plots/core_only/greyscale/GPF_plot_18_core.png")
saveplot(GPF_plot_20_core,"Results/GPF/Plots/core_only/greyscale/GPF_plot_20_core.png")

# colour plots
GPF_plot_10_core_col <- CLplotFun(spdf.df_10.core, survey.area.core, "abundance")
GPF_plot_11_core_col <- CLplotFun(spdf.df_11.core, survey.area.core, "abundance")
GPF_plot_13_core_col <- CLplotFun(spdf.df_13.core, survey.area.core, "abundance")
GPF_plot_14_core_col <- CLplotFun(spdf.df_14.core, survey.area.core, "abundance")
GPF_plot_16_core_col <- CLplotFun(spdf.df_16.core, survey.area.core, "abundance")
GPF_plot_18_core_col <- CLplotFun(spdf.df_18.core, survey.area.core, "abundance")
GPF_plot_20_core_col <- CLplotFun(spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(GPF_plot_10_core_col,"Results/GPF/Plots/core_only/colour/GPF_plot_10_core_col.png")
saveplot(GPF_plot_11_core_col,"Results/GPF/Plots/core_only/colour/GPF_plot_11_core_col.png")
saveplot(GPF_plot_13_core_col,"Results/GPF/Plots/core_only/colour/GPF_plot_13_core_col.png")
saveplot(GPF_plot_14_core_col,"Results/GPF/Plots/core_only/colour/GPF_plot_14_core_col.png")
saveplot(GPF_plot_16_core_col,"Results/GPF/Plots/core_only/colour/GPF_plot_16_core_col.png")
saveplot(GPF_plot_18_core_col,"Results/GPF/Plots/core_only/colour/GPF_plot_18_core_col.png")
saveplot(GPF_plot_20_core_col,"Results/GPF/Plots/core_only/colour/GPF_plot_20_core_col.png")



## plot grids (abundance and variance)

# greyscale 
gpf_2yrs_gs <- 
  GPF_plot_10_core + GPF_varplot_final10.core.bw  +
  GPF_plot_20_core + GPF_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
gpf_2yrs_gs[[1]] <- gpf_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
gpf_2yrs_gs[[2]] <- gpf_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
gpf_2yrs_gs[[2]] <- gpf_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
gpf_2yrs_gs[[4]] <- gpf_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(gpf_2yrs_gs, "Results/GPF/Plots/core_only/greyscale/plot_grids/gpf_2yrs_gs.png")


    ## Plotting discrete bins ####
      # Add bins to SPDF - don't repeat ####

# this is the process of adding discrete bins to the abundance SPDFs. I will save them in a new folder so this only has to be done once

# Load original spatial dataframes
spdf.df_10.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_10.core.csv")
spdf.df_11.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_11.core.csv")
spdf.df_13.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_13.core.csv")
spdf.df_14.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_14.core.csv")
spdf.df_16.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_16.core.csv")
spdf.df_18.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_18.core.csv")
spdf.df_20.core <- read.csv("Results/GPF/Plots/core_only/spdf/spdf.df_20.core.csv")

# put spdf's into a list
dfs <- list(spdf.df_10.core,spdf.df_11.core,spdf.df_13.core,spdf.df_14.core,spdf.df_16.core,
            spdf.df_18.core,spdf.df_20.core)

# name the elements
names(dfs) <- c("spdf.df_10.core","spdf.df_11.core","spdf.df_13.core","spdf.df_14.core","spdf.df_16.core",
                "spdf.df_18.core","spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunAbun)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(spdf.df_10.core,file="Results/GPF/Plots/core_only/spdf/bins/spdf.df_10.core.csv")
write.csv(spdf.df_11.core,file="Results/GPF/Plots/core_only/spdf/bins/spdf.df_11.core.csv")
write.csv(spdf.df_13.core,file="Results/GPF/Plots/core_only/spdf/bins/spdf.df_13.core.csv")
write.csv(spdf.df_14.core,file="Results/GPF/Plots/core_only/spdf/bins/spdf.df_14.core.csv")
write.csv(spdf.df_16.core,file="Results/GPF/Plots/core_only/spdf/bins/spdf.df_16.core.csv")
write.csv(spdf.df_18.core,file="Results/GPF/Plots/core_only/spdf/bins/spdf.df_18.core.csv")
write.csv(spdf.df_20.core,file="Results/GPF/Plots/core_only/spdf/bins/spdf.df_20.core.csv")


      # Plotting ####

# Load spatial dataframes (with bins)
spdf.df_10.core <- read.csv("Results/GPF/Plots/core_only/spdf/bins/spdf.df_10.core.csv")
spdf.df_11.core <- read.csv("Results/GPF/Plots/core_only/spdf/bins/spdf.df_11.core.csv")
spdf.df_13.core <- read.csv("Results/GPF/Plots/core_only/spdf/bins/spdf.df_13.core.csv")
spdf.df_14.core <- read.csv("Results/GPF/Plots/core_only/spdf/bins/spdf.df_14.core.csv")
spdf.df_16.core <- read.csv("Results/GPF/Plots/core_only/spdf/bins/spdf.df_16.core.csv")
spdf.df_18.core <- read.csv("Results/GPF/Plots/core_only/spdf/bins/spdf.df_18.core.csv")
spdf.df_20.core <- read.csv("Results/GPF/Plots/core_only/spdf/bins/spdf.df_20.core.csv")

# change group2 (abundance) to factor and re-order. I've only done it for 2020 here but you can copy the code for the other years if you need to
spdf.df_20.core$group2 <- as.factor(spdf.df_20.core$group2)
spdf.df_20.core$group2 <- factor(spdf.df_20.core$group2, levels=c("High","Medium","Low","Very low"))

## plot greyscale
GPF_10_plot_bin_GS <- GSplotBin(spdf.df_10.core,survey.area.core,"Abundance","Relative abundance")
GPF_11_plot_bin_GS <- GSplotBin(spdf.df_11.core,survey.area.core,"Abundance","Relative abundance")
GPF_13_plot_bin_GS <- GSplotBin(spdf.df_13.core,survey.area.core,"Abundance","Relative abundance")
GPF_14_plot_bin_GS <- GSplotBin(spdf.df_14.core,survey.area.core,"Abundance","Relative abundance")
GPF_16_plot_bin_GS <- GSplotBin(spdf.df_16.core,survey.area.core,"Abundance","Relative abundance")
GPF_18_plot_bin_GS <- GSplotBin(spdf.df_18.core,survey.area.core,"Abundance","Relative abundance")
GPF_20_plot_bin_GS <- GSplotBin(spdf.df_20.core,"group2",survey.area.core,"Abundance","Relative abundance")


# save 
saveplot(GPF_10_plot_bin_GS,"Results/GPF/Plots/core_only/bins/GPF_10_plot_bin_GS.png")
saveplot(GPF_11_plot_bin_GS,"Results/GPF/Plots/core_only/bins/GPF_11_plot_bin_GS.png")
saveplot(GPF_13_plot_bin_GS,"Results/GPF/Plots/core_only/bins/GPF_13_plot_bin_GS.png")
saveplot(GPF_14_plot_bin_GS,"Results/GPF/Plots/core_only/bins/GPF_14_plot_bin_GS.png")
saveplot(GPF_16_plot_bin_GS,"Results/GPF/Plots/core_only/bins/GPF_16_plot_bin_GS.png")
saveplot(GPF_18_plot_bin_GS,"Results/GPF/Plots/core_only/bins/GPF_18_plot_bin_GS.png")
saveplot(GPF_20_plot_bin_GS,"Results/GPF/Plots/core_only/bins/GPF_20_plot_bin_GS.png")



## Variance estimation ####

## Becasue we have estimated Nhat we need to use the dsm.var.gam function to estiamte variance
## We can't etimate variance for each cell as there are to many cells (requires R to allocate a vector of size ~30Gb), so following advice fom David Miller we will split cells into chunks and estimate variance by chunk



# estimate variance
gpf.var.Final10.core <- varEstfun(preddata10_core, gpfDSM.nb.6)
gpf.var.Final11.core <- varEstfun(preddata11_core, gpfDSM.nb.6)
gpf.var.Final13.core <- varEstfun(preddata13_core, gpfDSM.nb.6)
gpf.var.Final14.core <- varEstfun(preddata14_core, gpfDSM.nb.6)
gpf.var.Final16.core <- varEstfun(preddata16_core, gpfDSM.nb.6)
gpf.var.Final18.core <- varEstfun(preddata18_core, gpfDSM.nb.6)
gpf.var.Final20.core <- varEstfun(preddata20_core, gpfDSM.nb.6)

# save variance estimates
write.csv(gpf.var.Final10.core, file="Results/GPF/core_only/gpf.var10.core.csv")
write.csv(gpf.var.Final11.core, file="Results/GPF/core_only/gpf.var11.core.csv")
write.csv(gpf.var.Final13.core, file="Results/GPF/core_only/gpf.var13.core.csv")
write.csv(gpf.var.Final14.core, file="Results/GPF/core_only/gpf.var14.core.csv")
write.csv(gpf.var.Final16.core, file="Results/GPF/core_only/gpf.var16.core.csv")
write.csv(gpf.var.Final18.core, file="Results/GPF/core_only/gpf.var18.core.csv")
write.csv(gpf.var.Final20.core, file="Results/GPF/core_only/gpf.var20.core.csv")

# create spdf's for plotting
gpf.var.spdf.df_10.core <- varPlotDF(gpf.var.Final10.core, pred.polys_200)
gpf.var.spdf.df_11.core <- varPlotDF(gpf.var.Final11.core, pred.polys_200)
gpf.var.spdf.df_13.core <- varPlotDF(gpf.var.Final13.core, pred.polys_200)
gpf.var.spdf.df_14.core <- varPlotDF(gpf.var.Final14.core, pred.polys_200)
gpf.var.spdf.df_16.core <- varPlotDF(gpf.var.Final16.core, pred.polys_200)
gpf.var.spdf.df_18.core <- varPlotDF(gpf.var.Final18.core, pred.polys_200)
gpf.var.spdf.df_20.core <- varPlotDF(gpf.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(gpf.var.spdf.df_10.core,
          file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_10.core.csv")
write.csv(gpf.var.spdf.df_11.core,
          file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_11.core.csv")
write.csv(gpf.var.spdf.df_13.core,
          file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_13.core.csv")
write.csv(gpf.var.spdf.df_14.core,
          file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_14.core.csv")
write.csv(gpf.var.spdf.df_16.core,
          file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_16.core.csv")
write.csv(gpf.var.spdf.df_18.core,
          file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_18.core.csv")
write.csv(gpf.var.spdf.df_20.core,
          file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_20.core.csv")


    ## Plotting variance ####
      # Calculate CV & add bins to SPDF - don't repeat ####

# Load spatial dataframes
gpf.var.spdf.df_10.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_10.core.csv")
gpf.var.spdf.df_11.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_11.core.csv")
gpf.var.spdf.df_13.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_13.core.csv")
gpf.var.spdf.df_14.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_14.core.csv")
gpf.var.spdf.df_16.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_16.core.csv")
gpf.var.spdf.df_18.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_18.core.csv")
gpf.var.spdf.df_20.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_20.core.csv")


# first need to calculate CV from the variance (note: need the abundance spdf's loaded)
gpf.var.spdf.df_10.core <- CVaddFun(spdf.df_10.core,gpf.var.spdf.df_10.core)
gpf.var.spdf.df_11.core <- CVaddFun(spdf.df_11.core,gpf.var.spdf.df_11.core)
gpf.var.spdf.df_13.core <- CVaddFun(spdf.df_13.core,gpf.var.spdf.df_13.core)
gpf.var.spdf.df_14.core <- CVaddFun(spdf.df_14.core,gpf.var.spdf.df_14.core)
gpf.var.spdf.df_16.core <- CVaddFun(spdf.df_16.core,gpf.var.spdf.df_16.core)
gpf.var.spdf.df_18.core <- CVaddFun(spdf.df_18.core,gpf.var.spdf.df_18.core)
gpf.var.spdf.df_20.core <- CVaddFun(spdf.df_20.core,gpf.var.spdf.df_20.core)


# Save the SPDF's with the CV value but no bins (for continuous plotting of the CV)
write.csv(gpf.var.spdf.df_10.core,file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_10.core.csv")
write.csv(gpf.var.spdf.df_11.core,file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_11.core.csv")
write.csv(gpf.var.spdf.df_13.core,file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_13.core.csv")
write.csv(gpf.var.spdf.df_14.core,file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_14.core.csv")
write.csv(gpf.var.spdf.df_16.core,file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_16.core.csv")
write.csv(gpf.var.spdf.df_18.core,file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_18.core.csv")
write.csv(gpf.var.spdf.df_20.core,file="Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_20.core.csv")


### add bins

## Quartiles

# put spdf's into a list
dfs <- list(gpf.var.spdf.df_10.core,gpf.var.spdf.df_11.core,gpf.var.spdf.df_13.core,gpf.var.spdf.df_14.core,
            gpf.var.spdf.df_16.core,gpf.var.spdf.df_18.core,gpf.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("gpf.var.spdf.df_10.core","gpf.var.spdf.df_11.core","gpf.var.spdf.df_13.core",
                "gpf.var.spdf.df_14.core","gpf.var.spdf.df_16.core","gpf.var.spdf.df_18.core",
                "gpf.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(gpf.var.spdf.df_10.core,file="Results/GPF/Plots/variance/core_only/spdf/bins/gpf.var.spdf.df_10.core.csv")
write.csv(gpf.var.spdf.df_11.core,file="Results/GPF/Plots/variance/core_only/spdf/bins/gpf.var.spdf.df_11.core.csv")
write.csv(gpf.var.spdf.df_13.core,file="Results/GPF/Plots/variance/core_only/spdf/bins/gpf.var.spdf.df_13.core.csv")
write.csv(gpf.var.spdf.df_14.core,file="Results/GPF/Plots/variance/core_only/spdf/bins/gpf.var.spdf.df_14.core.csv")
write.csv(gpf.var.spdf.df_16.core,file="Results/GPF/Plots/variance/core_only/spdf/bins/gpf.var.spdf.df_16.core.csv")
write.csv(gpf.var.spdf.df_18.core,file="Results/GPF/Plots/variance/core_only/spdf/bins/gpf.var.spdf.df_18.core.csv")
write.csv(gpf.var.spdf.df_20.core,file="Results/GPF/Plots/variance/core_only/spdf/bins/gpf.var.spdf.df_20.core.csv")




## custom bins

# put spdf's into a list
dfs <- list(gpf.var.spdf.df_10.core,gpf.var.spdf.df_11.core,gpf.var.spdf.df_13.core,gpf.var.spdf.df_14.core,
            gpf.var.spdf.df_16.core,gpf.var.spdf.df_18.core,gpf.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("gpf.var.spdf.df_10.core","gpf.var.spdf.df_11.core","gpf.var.spdf.df_13.core",
                "gpf.var.spdf.df_14.core","gpf.var.spdf.df_16.core","gpf.var.spdf.df_18.core",
                "gpf.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar2)

# split elements into original dataframes
list2env(dfs, globalenv())


# save the SPDFs with the custom bins
write.csv(gpf.var.spdf.df_10.core,
          file="Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_10.core.csv")
write.csv(gpf.var.spdf.df_11.core,
          file="Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_11.core.csv")
write.csv(gpf.var.spdf.df_13.core,
          file="Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_13.core.csv")
write.csv(gpf.var.spdf.df_14.core,
          file="Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_14.core.csv")
write.csv(gpf.var.spdf.df_16.core,
          file="Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_16.core.csv")
write.csv(gpf.var.spdf.df_18.core,
          file="Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_18.core.csv")
write.csv(gpf.var.spdf.df_20.core,
          file="Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_20.core.csv")


      # Discrete bins ####

# load spdfs (which has already had CV calculated and then put into bins)
gpf.var.spdf.df_10.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_10.core.csv")
gpf.var.spdf.df_11.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_11.core.csv")
gpf.var.spdf.df_13.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_13.core.csv")
gpf.var.spdf.df_14.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_14.core.csv")
gpf.var.spdf.df_16.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_16.core.csv")
gpf.var.spdf.df_18.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_18.core.csv")
gpf.var.spdf.df_20.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/bins/custom/gpf.var.spdf.df_20.core.csv")

# make CV into factor and re-order
gpf.var.spdf.df_10.core$CV <- as.factor(gpf.var.spdf.df_10.core$CV)
gpf.var.spdf.df_10.core$CV <- factor(gpf.var.spdf.df_10.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
gpf.var.spdf.df_11.core$CV <- as.factor(gpf.var.spdf.df_11.core$CV)
gpf.var.spdf.df_11.core$CV <- factor(gpf.var.spdf.df_11.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
gpf.var.spdf.df_13.core$CV <- as.factor(gpf.var.spdf.df_13.core$CV)
gpf.var.spdf.df_13.core$CV <- factor(gpf.var.spdf.df_13.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
gpf.var.spdf.df_14.core$CV <- as.factor(gpf.var.spdf.df_14.core$CV)
gpf.var.spdf.df_14.core$CV <- factor(gpf.var.spdf.df_14.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
gpf.var.spdf.df_16.core$CV <- as.factor(gpf.var.spdf.df_16.core$CV)
gpf.var.spdf.df_16.core$CV <- factor(gpf.var.spdf.df_16.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
gpf.var.spdf.df_18.core$CV <- as.factor(gpf.var.spdf.df_18.core$CV)
gpf.var.spdf.df_18.core$CV <- factor(gpf.var.spdf.df_18.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
gpf.var.spdf.df_20.core$CV <- as.factor(gpf.var.spdf.df_20.core$CV)
gpf.var.spdf.df_20.core$CV <- factor(gpf.var.spdf.df_20.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))

# plot CV in bins
GPF_10_plot_bin_GS_var <- GSplotBin(gpf.var.spdf.df_10.core,"CV",survey.area.core,"Variance","CV")
GPF_11_plot_bin_GS_var <- GSplotBin(gpf.var.spdf.df_11.core,"CV",survey.area.core,"Variance","CV")
GPF_13_plot_bin_GS_var <- GSplotBin(gpf.var.spdf.df_13.core,"CV",survey.area.core,"Variance","CV")
GPF_14_plot_bin_GS_var <- GSplotBin(gpf.var.spdf.df_14.core,"CV",survey.area.core,"Variance","CV")
GPF_16_plot_bin_GS_var <- GSplotBin(gpf.var.spdf.df_16.core,"CV",survey.area.core,"Variance","CV")
GPF_18_plot_bin_GS_var <- GSplotBin(gpf.var.spdf.df_18.core,"CV",survey.area.core,"Variance","CV")
GPF_20_plot_bin_GS_var <- GSplotBin(gpf.var.spdf.df_20.core,"CV",survey.area.core,"Variance","CV")

# save 
saveplot(GPF_10_plot_bin_GS_var,"Results/GPF/Plots/variance/core_only/greyscale/bins/GPF_10_plot_bin_GS_var.png")
saveplot(GPF_11_plot_bin_GS_var,"Results/GPF/Plots/variance/core_only/greyscale/bins/GPF_11_plot_bin_GS_var.png")
saveplot(GPF_13_plot_bin_GS_var,"Results/GPF/Plots/variance/core_only/greyscale/bins/GPF_13_plot_bin_GS_var.png")
saveplot(GPF_14_plot_bin_GS_var,"Results/GPF/Plots/variance/core_only/greyscale/bins/GPF_14_plot_bin_GS_var.png")
saveplot(GPF_16_plot_bin_GS_var,"Results/GPF/Plots/variance/core_only/greyscale/bins/GPF_16_plot_bin_GS_var.png")
saveplot(GPF_18_plot_bin_GS_var,"Results/GPF/Plots/variance/core_only/greyscale/bins/GPF_18_plot_bin_GS_var.png")
saveplot(GPF_20_plot_bin_GS_var,"Results/GPF/Plots/variance/core_only/greyscale/bins/GPF_20_plot_bin_GS_var.png")



      # Continuous ####

# Load spatial dataframes
gpf.var.spdf.df_10.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_10.core.csv")
gpf.var.spdf.df_11.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_11.core.csv")
gpf.var.spdf.df_13.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_13.core.csv")
gpf.var.spdf.df_14.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_14.core.csv")
gpf.var.spdf.df_16.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_16.core.csv")
gpf.var.spdf.df_18.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_18.core.csv")
gpf.var.spdf.df_20.core <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_20.core.csv")

# change group2 name to CV
gpf.var.spdf.df_10.core <- gpf.var.spdf.df_10.core %>% dplyr::rename(CV=group2)
gpf.var.spdf.df_11.core <- gpf.var.spdf.df_11.core %>% dplyr::rename(CV=group2)
gpf.var.spdf.df_13.core <- gpf.var.spdf.df_13.core %>% dplyr::rename(CV=group2)
gpf.var.spdf.df_14.core <- gpf.var.spdf.df_14.core %>% dplyr::rename(CV=group2)
gpf.var.spdf.df_16.core <- gpf.var.spdf.df_16.core %>% dplyr::rename(CV=group2)
gpf.var.spdf.df_18.core <- gpf.var.spdf.df_18.core %>% dplyr::rename(CV=group2)
gpf.var.spdf.df_20.core <- gpf.var.spdf.df_20.core %>% dplyr::rename(CV=group2)


# greyscale plots
GPF_varplot_final10.core.bw <- GSplotFun(gpf.var.spdf.df_10.core,survey.area.core,"CV", "")
GPF_varplot_final11.core.bw <- GSplotFun(gpf.var.spdf.df_11.core,survey.area.core,"group2", "2011")
GPF_varplot_final13.core.bw <- GSplotFun(gpf.var.spdf.df_13.core,survey.area.core,"group2", "2013")
GPF_varplot_final14.core.bw <- GSplotFun(gpf.var.spdf.df_14.core,survey.area.core,"group2", "2014")
GPF_varplot_final16.core.bw <- GSplotFun(gpf.var.spdf.df_16.core,survey.area.core,"group2", "2016")
GPF_varplot_final18.core.bw <- GSplotFun(gpf.var.spdf.df_18.core,survey.area.core,"group2", "2018")
GPF_varplot_final20.core.bw <- GSplotFun(gpf.var.spdf.df_20.core,survey.area.core,"group2", "2020")

# save greyscale
saveplot(GPF_varplot_final10.core.bw, "Results/GPF/Plots/variance/core_only/greyscale/2010_GPF_var.core.bw.png")
saveplot(GPF_varplot_final11.core.bw, "Results/GPF/Plots/variance/core_only/greyscale/2011_GPF_var.core.bw.png")
saveplot(GPF_varplot_final13.core.bw, "Results/GPF/Plots/variance/core_only/greyscale/2013_GPF_var.core.bw.png")
saveplot(GPF_varplot_final14.core.bw, "Results/GPF/Plots/variance/core_only/greyscale/2014_GPF_var.core.bw.png")
saveplot(GPF_varplot_final16.core.bw, "Results/GPF/Plots/variance/core_only/greyscale/2016_GPF_var.core.bw.png")
saveplot(GPF_varplot_final18.core.bw, "Results/GPF/Plots/variance/core_only/greyscale/2018_GPF_var.core.bw.png")
saveplot(GPF_varplot_final20.core.bw, "Results/GPF/Plots/variance/core_only/greyscale/2020_GPF_var.core.bw.png")

# colour plots
GPF_carplot_final10.core.col <- CLplotFun(gpf.var.spdf.df_10.core,survey.area.core,"variance")
GPF_carplot_final11.core.col <- CLplotFun(gpf.var.spdf.df_11.core,survey.area.core,"variance")
GPF_carplot_final13.core.col <- CLplotFun(gpf.var.spdf.df_13.core,survey.area.core,"variance")
GPF_carplot_final14.core.col <- CLplotFun(gpf.var.spdf.df_14.core,survey.area.core,"variance")
GPF_carplot_final16.core.col <- CLplotFun(gpf.var.spdf.df_16.core,survey.area.core,"variance")
GPF_carplot_final18.core.col <- CLplotFun(gpf.var.spdf.df_18.core,survey.area.core,"variance")
GPF_carplot_final20.core.col <- CLplotFun(gpf.var.spdf.df_20.core,survey.area.core,"variance")

# save colour
saveplot(GPF_carplot_final10.core.col, "Results/GPF/Plots/variance/core_only/colour/2010_GPF_var.core.col.png")
saveplot(GPF_carplot_final11.core.col, "Results/GPF/Plots/variance/core_only/colour/2011_GPF_var.core.col.png")
saveplot(GPF_carplot_final13.core.col, "Results/GPF/Plots/variance/core_only/colour/2013_GPF_var.core.col.png")
saveplot(GPF_carplot_final14.core.col, "Results/GPF/Plots/variance/core_only/colour/2014_GPF_var.core.col.png")
saveplot(GPF_carplot_final16.core.col, "Results/GPF/Plots/variance/core_only/colour/2016_GPF_var.core.col.png")
saveplot(GPF_carplot_final18.core.col, "Results/GPF/Plots/variance/core_only/colour/2018_GPF_var.core.col.png")
saveplot(GPF_carplot_final20.core.col, "Results/GPF/Plots/variance/core_only/colour/2020_GPF_var.core.col.png")



#### Black-shanked Douc #######################################################
## Load data ####

# Observation data. Unique to species 
bsd_obsdata <- read.csv("Species_Data/BSD/R Data/obsdata.csv", header = TRUE) 
bsd_obsdata$object <- as.factor(bsd_obsdata$object) 
bsd_obsdata$Sample.Label <- as.factor(bsd_obsdata$Sample.Label)
str(bsd_obsdata)
head(bsd_obsdata)

# Transect data. Unique to species
bsd_distdata <- read.csv("Species_Data/BSD/R Data/distdata.csv", header = TRUE)
bsd_distdata$object <- as.factor(bsd_distdata$object)
bsd_distdata$NameObserver <- as.factor(bsd_distdata$NameObserver)
bsd_distdata$transect <- as.factor(bsd_distdata$transect)
bsd_distdata$stratum <- as.vector(scale(bsd_distdata$year, center = T, scale = T)) 
bsd_distdata$year <- as.factor(bsd_distdata$year)
bsd_distdata$date <- as.Date(bsd_distdata$date, format = "%d/%m/%Y")
str(bsd_distdata)
head(bsd_distdata)

## Plot the covariates across the grid, with group sizes #### 

# Warning - plots take a few minutes to run

# habitat
plot_BSDobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata18_core$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=bsd_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/BSD/plot_BSDobs_habitat.png", plot_BSDobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_BSDobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=bsd_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/BSD/plot_BSDobs_dstWater.png", plot_BSDobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_BSDobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=bsd_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/BSD/plot_BSDobs_dstStlmnt.png", plot_BSDobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_BSDobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=bsd_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/BSD/plot_BSDobs_dstRoad.png", plot_BSDobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_BSDobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=bsd_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/BSD/plot_BSDobs_dstBorder.png", plot_BSDobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_BSDobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=bsd_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/BSD/plot_BSDobs_dstStation.png", plot_BSDobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_BSDobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=bsd_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/BSD/plot_BSDobs_dstELC.png", plot_BSDobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_BSDobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=bsd_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/BSD/plot_BSDobs_elev.png", plot_BSDobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)


## Exploratory plots & linear models ####

## Checking that there are no large gaps in the range of variables for douc observations

# subset segdata to get only the segments with BSD observations
bsd_varcheck <- segdata[match(bsd_obsdata$Sample.Label,segdata$Sample.Label), ]

# habitat - pretty low number of observations for the NF level (and DF & W but these are low in segdata too).  Is this a problem?
plot(segdata$habitat)
plot(bsd_varcheck$habitat)

# dstWater - looks good
hist(segdata$dstWater)
hist(bsd_varcheck$dstWater)

# dstStlmnt - bit of a gap between 8km and 10km. This is the same as GPF although for BSD there are observations beyond 10km.
hist(segdata$dstStlmnt)
hist(bsd_varcheck$dstStlmnt)

# dstRoad - quite a large gap here -there are no observations beyond 2km, and the variable range goes up to 4km (so the observations only cover 50% of the variable range). This could well be a problem.
hist(segdata$dstRoad)
hist(bsd_varcheck$dstRoad)

# dstBorder - small gap between 445km and 55km but not a problem I don't think
hist(segdata$dstBorder)
hist(bsd_varcheck$dstBorder)

# dstStation - fine
hist(segdata$dstStation)
hist(bsd_varcheck$dstStation)

# elevation - gap between 600m and 700m but not a problem I don't think
hist(segdata$elevation)
hist(bsd_varcheck$elevation)

## Histograms

# One massive outlier which I will remove for the histograms. Will also remove the 4 other observations >100m
bsd_distdata <- bsd_distdata %>% filter(distance<100)

# Distance 
bsd_h1 <- ggplot(bsd_distdata, aes(distance))+ geom_histogram(binwidth = 1)
bsd_h2 <- ggplot(bsd_distdata, aes(distance))+ geom_histogram(binwidth = 5)
bsd_h3 <- ggplot(bsd_distdata, aes(distance))+ geom_histogram(binwidth = 10)
bsd_h4 <- ggplot(bsd_distdata, aes(distance))+ geom_histogram(binwidth = 15)
bsd_h5 <- ggplot(bsd_distdata, aes(distance))+ geom_histogram(binwidth = 20)
bsd_h6 <- ggplot(bsd_distdata, aes(distance))+ geom_histogram(binwidth = 40)
bsd_h1+bsd_h2+bsd_h3+bsd_h4+bsd_h5+bsd_h6
# Evidence of evasive movement. Not surprising for this species. Also the lumping around 0m that we discovered during the CDS analysis is evident here. Binning the data will be necessary 

# find the correct bins. Data will be truncated at 50m (segment size)
hist(bsd_distdata$distance[bsd_distdata$distance<50], breaks=c(0,5,15,20,25,30,35,40,45,50))
hist(bsd_distdata$distance[bsd_distdata$distance<50], breaks=c(0,7,14,20,25,30,35,40,45,50)) # use


# cluster size, observer, habitat, year, month, transect
bsd_h7 <- ggplot(bsd_distdata, aes(size))+geom_histogram(binwidth = 0.5)
bsd_h8 <- ggplot(bsd_distdata, aes(NameObserver))+geom_histogram(stat="count")
bsd_h9 <- ggplot(bsd_distdata, aes(habitat))+geom_histogram(stat="count")
bsd_h10 <- ggplot(bsd_distdata, aes(year))+geom_histogram(stat="count")
bsd_h11 <- ggplot(bsd_distdata, aes(month))+geom_histogram(stat="count")
bsd_h12 <- ggplot(bsd_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(bsd_h7,bsd_h8,bsd_h9,bsd_h10,bsd_h11,bsd_h12)
# The most frequent group size is 2, which is clearly not the correct group size. It is extremely difficult to accurately observe all individuals in a group for this species. This suggests that we should be looking at density/abundance of groups. This will probably be the case for all gregarious species.

## Plots of distance against variables
plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}


bsd_d1 <- ggplot(bsd_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
bsd_d2 <- ggplot(bsd_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
bsd_d3 <- ggplot(bsd_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
bsd_d4 <- ggplot(bsd_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
bsd_d5 <- ggplot(bsd_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(bsd_d1,bsd_d2,bsd_d3,bsd_d4,bsd_d5)

# As we would expect, the mean distance of observations is higher in open forest and non-forest, when compared with dense forest. There is not a huge amount to take from the observer plot, firstly as there is not much difference between them and secondly because as Olly mentioned, the team members tend to stick to spatially clustered transects wich are likely to be similar habitats, and so it's not really a fair test.  The distance by month plot shows perhaps a slight overall increase in mean distance as the season progresses, but this is unlikely to be significant. The distance~group size plot looks interesting, but we will test that relationship with a model below


## Plots of cluster size against variables
bsd_s1 <- ggplot(bsd_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
bsd_s2 <- ggplot(bsd_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
bsd_s3 <- ggplot(bsd_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
bsd_s4 <- ggplot(bsd_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
bsd_s5 <- ggplot(bsd_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(bsd_s1,bsd_s2,bsd_s3,bsd_s4,bsd_s5)

# I am very sceptical of the group sizes that are recorded for primate species.  This is not a critisism, simply and observation.  It is extremely difficult to accurately count the number of doucs in a group in dense forest.  Therefore I am sceptical of these plots too.  Nevertheless, the group size~month plot shows that perhaps the mean group size increases throughout the season, which might make sense due to few leaves and thinner foliage as the dry season progresses.  It's good to see generally consistent group size estimates over the years (except 2016), as this suggests consistency in data collection from the monitoring team.   

## Linear models

# group size ~ distance
newdist <- data.frame(distance=seq(0,92,len=10))
lm1 <- lm(size~distance, data=bsd_distdata)
plot(bsd_distdata$size~bsd_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)

# The model suggests that there is a relationship between group size and distance, although based on the model coefficient, across our 50m truncation distance group size is only likely to increase by 0.5 individuals.  I don't think this represents a size bias issue that we need to worry about. I will test cluster size as a covariate in the DF models though

## Estimating detection function ####

# I will retain the below DF model selection for BSD because for the DSMs data from all years are used, rather than from individual years as in the CDS.  Initially I estiamted annual detection functions with the intention of modelling DSMs annually, but due to time constraints, and for the sake of consistency with the other species, I have used all data together.    


# All models will be truncated at 50m. 

# There is only 1 observation for NameObserver == 14 so I am going to remove it (it's messing with the models)
bsd_distdata <- bsd_distdata %>% filter(NameObserver != 14) 

# Subset the distdata into the different years
bsd_distdata10 <- bsd_distdata %>% filter(year == 2010)
bsd_distdata11 <- bsd_distdata %>% filter(year == 2011)
bsd_distdata13 <- bsd_distdata %>% filter(year == 2013)
bsd_distdata14 <- bsd_distdata %>% filter(year == 2014)
bsd_distdata16 <- bsd_distdata %>% filter(year == 2016)
bsd_distdata18 <- bsd_distdata %>% filter(year == 2018)


  ## All observations ####
    ## Models with NO covariates ####


# Model with no covariates, Uniform key and Cosine adjustment
bsdDF.un <- ds(bsd_distdata, truncation = 50, key = "unif", adjustment = "cos",
               cutpoints=c(0,7,14,20,25,30,35,40,45,50))
summary(bsdDF.un)
par(mfrow=c(1,2))
plot(bsdDF.un, showpoints=FALSE, pl.den=0, lwd=2)
ddf.gof(bsdDF.un$ddf)


# Model with no covariates, Uniform Key and simple polynomial adjustment 
bsdDF.un2 <- ds(bsd_distdata, truncation = 50, key = "unif", adjustment = "poly",
                cutpoints=c(0,7,14,20,25,30,35,40,45,50))
summary(bsdDF.un2)
par(mfrow=c(1,2))
plot(bsdDF.un2, showpoints=FALSE, pl.den=0, lwd=2)
ddf.gof(bsdDF.un2$ddf)
 


# Model with no covariates, half-normal key and no adjustment
bsdDF.hn <- ds(bsd_distdata, truncation = 50, key = "hn",
               cutpoints=c(0,7,14,20,25,30,35,40,45,50))
summary(bsdDF.hn)
par(mfrow=c(1,2))
plot(bsdDF.hn, showpoints=FALSE, pl.den=0, lwd=2)
ddf.gof(bsdDF.hn$ddf)



# Model with no covariates, half-normal key and cosine adjustment
bsdDF.hn2 <- ds(bsd_distdata, truncation = 50, key = "hn", adjustment = "cos",
                cutpoints=c(0,7,14,20,25,30,35,40,45,50))
summary(bsdDF.hn2)
par(mfrow=c(1,2))
plot(bsdDF.hn2, showpoints=FALSE, pl.den=0, lwd=2)
ddf.gof(bsdDF.hn2$ddf)


# Model with no covariates, half-normal key and hermite polynomial adjustment
bsdDF.hn3 <- ds(bsd_distdata, truncation = 50, key = "hn", adjustment = "herm",
                cutpoints=c(0,7,14,20,25,30,35,40,45,50))
summary(bsdDF.hn3)
par(mfrow=c(1,2))
plot(bsdDF.hn3, showpoints=FALSE, pl.den=0, lwd=2)
ddf.gof(bsdDF.hn3$ddf)
 



# Model with no covariates, hazard rate key, no adjustment
bsdDF.hr <- ds(bsd_distdata, truncation = 50, key = "hr",
               cutpoints=c(0,7,14,20,25,30,35,40,45,50))
summary(bsdDF.hr)
par(mfrow=c(1,2))
plot(bsdDF.hr, showpoints=FALSE, pl.den=0, lwd=2)
ddf.gof(bsdDF.hr$ddf)


# Model with no covariates, hazard rate key, cosine adjustment
bsdDF.hr2 <- ds(bsd_distdata, truncation = 50, key = "hr", adjustment = "cos",
                cutpoints=c(0,7,14,20,25,30,35,40,45,50))
summary(bsdDF.hr2)
par(mfrow=c(1,2))
plot(bsdDF.hr2, showpoints=FALSE, pl.den=0, lwd=2)
ddf.gof(bsdDF.hr2$ddf)
 


# Model with no covariates, hazard rate key, simple polynomial adjustment
bsdDF.hr4 <- ds(bsd_distdata, truncation = 50, key = "hr", adjustment = "poly",
                cutpoints=c(0,7,14,20,25,30,35,40,45,50))
summary(bsdDF.hr4)
par(mfrow=c(1,2))
plot(bsdDF.hr4, showpoints=FALSE, pl.den=0, lwd=2)
ddf.gof(bsdDF.hr4$ddf)


    ## Models WITH covariates ####
 
detfunc <- function(dat,key,covar) {
  
  name <- ds(data=dat, truncation = 50, key=key, formula = covar,
             cutpoints=c(0,7,14,20,25,30,35,40,45,50))
  
  par(mfrow=c(1,2))
  plot(name, showpoints=FALSE, pl.den=0, lwd=2)
  ddf.gof(name$ddf)
  summary(name)
}

bsdDF.cov.size <- detfunc(bsd_distdata,"hn",~size)
bsdDF.cov.obs <- detfunc(bsd_distdata[bsd_distdata$NameObserver != "14",],"hn",~NameObserver) 
# remove observer 14
bsdDF.cov.tran <- detfunc(bsd_distdata,"hn",~transect) # takes too long - too many levels
bsdDF.cov.strat <- detfunc(bsd_distdata,"hn",~stratum)
bsdDF.cov.month <- detfunc(bsd_distdata,"hn",~month)

# re-run model outside my function otherwise you get an error in the dsm() call
bsd_distdata <- bsd_distdata[bsd_distdata$NameObserver != "14",] # remove 14 from data
bsdDF.cov.obs <- ds(bsd_distdata, truncation=50, key="hn", formula=~NameObserver,
                    cutpoints=c(0,7,14,20,25,30,35,40,45,50))



    ## Results summaries ####

### no covars
BSD_DF_summary_key <- data.frame(name = c("bsdDF.un","bsdDF.un2","bsdDF.hn","bsdDF.hr",
                                      "bsdDF.hr2","bsdDF.hr4"),
                             p = c(0.67,0.65,0.64,0.65,0.65,0.66),
                             SE = c(0.02,0.01,0.01,0.02,0.02,0.02),
                             CV = c(0.03,0.02,0.02,0.03,0.03,0.03),
                             AIC = c(10323,10323,10323,10324,10324,10326))

BSD_DF_summary_key <- BSD_DF_summary_key %>% 
                  arrange(AIC,SE) 
BSD_DF_summary_key 
# all models except Hr4 have some support. 

# plot together
par(mfrow=c(3,2))
plot(bsdDF.hn, main="bsdDF.hn")
plot(bsdDF.un2, main="bsdDF.un2")
plot(bsdDF.un, main="bsdDF.un")
plot(bsdDF.hr, main="bsdDF.hr")
plot(bsdDF.hr2, main="bsdDF.hr2")

# all very simlar. I want to test covars so will exlude Unif models. HR models look to be slightly overfitting, and so I will select HN model.


### covars
BSD_DF_summary_cov <- data.frame(name = c("bsdDF.cov.size","bsdDF.cov.obs",
                                          "bsdDF.cov.strat","bsdDF.cov.month"),
                                 variables = c("size","observer","stratum","month"),
                                 p = c(0.64,0.63,0.64,0.64),
                                 SE = c(0.01,0.01,0.01,0.01),
                                 CV = c(0.02,0.02,0.02,0.02),
                                 AIC = c(10305,10276,10310,10324))

BSD_DF_summary_cov <- BSD_DF_summary_cov %>% 
  arrange(AIC,SE) 
# model with observer has overwhelming support

# model with covariate (observer) has significantly lower AIC than the HN-only model (dAIC == 47) and so is selected

# bsdDF.cov.obs

  
## Fitting a spatial model ####

# Best GLOBAL detection function model is bsdDF.cov.obs (observer)
# Best 2010 detection function model is bsdDF10.hn.cov7 (month)
# Best 2011 detection function model is bsdDF11.hn.cov4 (observer)
# Best 2013 detection function model is bsdDF13.hn.cov7 (month)
# Best 2014 detection function model is bsdDF14.hn.cov4 (observer)
# Best 2016 detection function model is bsdDF16.hr.cov5 (month + habitat) 
# Best 2018 detection function model is bsdDF18.un (no covariates)

# There are covariates in the DF so we have to use Nhat

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated. This is because I don't think the group size estimates for this species are accurate enough

# I am not including dstELC in the models for BSD

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
bsd_obsdata <- bsd_obsdata %>% filter(distance <= 50)


  ## Global (all years) ####
    ## Quasipoisson response ####

# Saturated model with all covariates. 
bsdDSM.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                         s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                         habitat,
                  bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.sat)
par(mfrow=c(3,2))
plot(bsdDSM.sat, scale = 0)
gam.check(bsdDSM.sat)
# DE = 27.1, R2=0.132. dstStation not sig, but could be linear. Plots looking overfitted. gam.check suggesting k too low but I disagree!

# Saturated model with reduced k 
bsdDSM.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=5)+ 
                    s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                    s(elevation,bs="ts",k=5)+ 
                    habitat,
                  bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.sat2)
par(mfrow=c(3,2))
plot(bsdDSM.sat2, scale = 0)
gam.check(bsdDSM.sat2)
# DE = 23.2, R2=0.107. All terms sig. These plots look much better - not overfitted.

# As above, but with k increased a bit 
bsdDSM.sat3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=7)+ s(dstStlmnt,bs="ts",k=7)+ 
                     s(dstBorder,bs="ts",k=7)+ s(dstStation,bs="ts",k=7)+ 
                     s(elevation,bs="ts",k=7)+ 
                     habitat,
                   bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.sat3)
par(mfrow=c(3,2))
plot(bsdDSM.sat3, scale = 0)
gam.check(bsdDSM.sat3)
# DE = 25.6, R2=0.127. All terms still sig. Some plots still fine, but dstBorder and dstStation are a bit too overfitted now. 



### the best QP model is bsdDSM.sat2

    ## Tweedie response ####

bsdDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(elevation,bs="ts")+
                            s(dstStation,bs="ts")+  
                            habitat,
                  bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.tw.sat)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.tw.sat)
dev.off()
par(mfrow=c(3,2))
plot(bsdDSM.tw.sat, scale=0)
# AIC = 20901, DE = 28. All terms, excpet dstWater are sig. All other plots look overfitted

# reduce k
bsdDSM.tw.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=5)+ 
                       s(dstBorder,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+
                       s(dstStation,bs="ts",k=5)+  
                       habitat,
                     bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                     family = tw(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.tw.sat2)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.tw.sat2)
dev.off()
par(mfrow=c(3,2))
plot(bsdDSM.tw.sat2, scale=0)
# AIC = 21468, DE = 23.1. dstWater now sig, although looks linear.  Plots look better

# dstWater as linear term
bsdDSM.tw.sat3 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=5)+ s(dstBorder,bs="ts",k=5)+ 
                        s(elevation,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+  
                        habitat + dstWater,
                      bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                      family = tw(), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.tw.sat3)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.tw.sat3)
dev.off()
par(mfrow=c(2,2))
plot(bsdDSM.tw.sat3, scale=0)
# AIC = 21469, DE = 23.1. AIC gone up a tiny bit, but dstWater sig as linear term.

# increase k a bit
bsdDSM.tw.sat4 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=7)+ s(dstBorder,bs="ts",k=7)+ 
                        s(elevation,bs="ts",k=7)+ s(dstStation,bs="ts",k=7)+  
                        habitat + dstWater,
                      bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                      family = tw(), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.tw.sat4)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.tw.sat4)
dev.off()
par(mfrow=c(2,2))
plot(bsdDSM.tw.sat4, scale=0)
# AIC = 21296, DE = 25.5. The dstStlmnt plot changes aren't too bad, but all the others start to look overfitted. Think I'll stick to the above




### bsdDSM.tw.sat3 is the best TW model

#
    ## Negative binomial response ####

bsdDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(elevation,bs="ts")+ 
                            s(dstStation,bs="ts")+ 
                            habitat,
                     bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.nb.sat)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.nb.sat)
dev.off()
par(mfrow=c(3,2))
plot(bsdDSM.nb.sat, scale=0)
# AIC = 14323, DE = 33. all terms sig. All plots look overfitted except dtWater. 


# reduce k for all except dstWater
bsdDSM.nb.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts",k=5)+ 
                       s(dstBorder,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                       s(dstStation,bs="ts",k=5)+ 
                       habitat,
                     bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                     family = nb(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.nb.sat2)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.nb.sat2)
dev.off()
par(mfrow=c(3,2))
plot(bsdDSM.nb.sat2, scale=0)
# AIC = 14679, DE = 27.4. dstStation no longer sig. plots look better

# increase k for dstStation
bsdDSM.nb.sat3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts",k=5)+ 
                        s(dstBorder,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                        s(dstStation,bs="ts",k=7)+ 
                        habitat,
                      bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                      family = nb(), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.nb.sat3)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.nb.sat3)
dev.off()
par(mfrow=c(3,2))
plot(bsdDSM.nb.sat3, scale=0)
# AIC = 14679, DE = 27.4. no difference, dstStation still not sig

# remove dstStation
bsdDSM.nb.4 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts",k=5)+ 
                        s(dstBorder,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                        habitat,
                      bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                      family = nb(), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.nb.4)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.nb.4)
dev.off()
par(mfrow=c(2,2))
plot(bsdDSM.nb.4, scale=0)
# AIC = 14688, DE = 27.2. AIC gone up but DE only down 0.2. Now dstWater not sig.

# set k for dstWater
bsdDSM.nb.5 <- dsm(Nhat ~ s(dstWater,bs="ts", k=10)+ s(dstStlmnt,bs="ts",k=5)+ 
                     s(dstBorder,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                     habitat,
                   bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.nb.5)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.nb.5)
dev.off()
par(mfrow=c(2,2))
plot(bsdDSM.nb.5, scale=0)
# AIC = 14688, DE = 27.2. dstWater still not sig

# remove dstWater
bsdDSM.nb.6 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=5)+ s(dstBorder,bs="ts",k=5)+ 
                     s(elevation,bs="ts",k=5)+ habitat,
                   bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.nb.6)
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.nb.6)
dev.off()
par(mfrow=c(2,2))
plot(bsdDSM.nb.6, scale=0)
# AIC = 14688, DE = 27.2. No change in aic or De. all terms sig


### bsdDSM.nb.6 is the best NB model

## Model selection ####

# Best QP model is bsdDSM.sat2
# Best TW model is bsdDSM.tw.sat3
# Best NB model is bsdDSM.nb.6

# NB distribution explains the most deviance
summary(bsdDSM.sat2) # DE = 23.2 %
summary(bsdDSM.tw.sat3) # DE = 23.1%
summary(bsdDSM.nb.6) # DE = 27.2%

# We can compare the TW and NB models with AIC
bsdDSM.tw.sat3$aic
bsdDSM.nb.6$aic
# Negative binomial model is better

# Compare NB mode with QP
anova(bsdDSM.sat2,bsdDSM.nb.6,test="Chisq")
# QP model has significantly less residual deviance

# Compare NB and QP models using Q-Q plot
dev.off()
par(mfrow=c(2,2))
gam.check(bsdDSM.sat2)
gam.check(bsdDSM.nb.6)
# NB qq plot looks better, although QP isn't horrendous. Response vs fitted plots for both look similar. 

### bsdDSM.sat2 is the best global model for BSD

## Autocorrelation ####


par(mfrow=c(1,1))
dsm.cor(bsdDSM.sat2, max.lag=10, Segment.Label="Sample.Label")
# yikes, some autocorrelation to deal with.

# after addding xy smooths
dsm.cor(bsdDSM.sat2.xy2, max.lag=10, Segment.Label="Sample.Label")


# add AR1 autocorrelation structure (always struggle with convergence)
segdata$sg.id <- as.numeric(sub("\\d+-","",segdata$Sample.Label))
segdata$tr.id <- as.numeric(segdata$Transect.Label)

bsdDSM.sat2.ar1 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=5)+ 
                     s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                     s(elevation,bs="ts",k=5)+ 
                     habitat,
                   bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gamm",
                   correlation=corAR1(form=~sg.id|tr.id),
                   segment.area = segdata$Sample.Area, group = TRUE)
# no convergence

# Add univeriate smooths of x and y
bsdDSM.sat2.xy <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=5)+ 
                         s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                         s(elevation,bs="ts",k=5)+ s(x,bs="ts")+s(y,bs="ts")+
                         habitat,
                       bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                       family = quasipoisson(link = "log"), engine = "gam",
                       segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.sat2.xy)
par(mfrow=c(3,3))
plot(bsdDSM.sat2.xy, scale = 0)
gam.check(bsdDSM.sat2.xy)
# DE = 28.5, R2=0.13. dstWater no longer sig

# remove dstWater
bsdDSM.sat2.xy2 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=5)+ s(dstBorder,bs="ts",k=5)+ 
                         s(dstStation,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                         s(x,bs="ts")+s(y,bs="ts")+ habitat,
                      bsdDF.cov.obs, segdata, bsd_obsdata, method = "REML",
                      family = quasipoisson(link = "log"), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(bsdDSM.sat2.xy2)
par(mfrow=c(3,2))
plot(bsdDSM.sat2.xy2, scale = 0, xlab="Distance to ranger station (m)")
gam.check(bsdDSM.sat2.xy2)
# DE = 28.5, R2=0.13.


# There was autocorrelation in the original best model (bsdDSM.sat2).  I tested adding a correlation structure to the model but there were issues with convergence. So I have added univariate smooths of x and y into the model, refit, and removed dstWater (as no longer sig).  There is still evidence of autocorrelation in the correlogram but according to David Miller if you have xy modelled then this sholdn't be too much of a problem (see https://groups.google.com/forum/#!searchin/distance-sampling/autocorrelation%7Csort:date/distance-sampling/BV2amUYMo6A/r3lWR9ztBQAJ).


## Abundance estimation ####

# bsdDSM.sat2.xy2 is the best model. 

# Predict over 2010 habitat
bsd.global.pred10.core <- predict(bsdDSM.sat2.xy2, preddata10_core, off.set = 40000)
write.csv(bsd.global.pred10.core, file="Results/BSD/core_only/bsd.pred10.core.csv")

# Predict over 2011 habitat
bsd.global.pred11.core <- predict(bsdDSM.sat2.xy2, preddata11_core, off.set = 40000)
write.csv(bsd.global.pred11.core, file="Results/BSD/core_only/bsd.pred11.core.csv")

# Predict over 2013 habitat
bsd.global.pred13.core <- predict(bsdDSM.sat2.xy2, preddata13_core, off.set = 40000)
write.csv(bsd.global.pred13.core, file="Results/BSD/core_only/bsd.pred13.core.csv")

# Predict over 2014 habitat
bsd.global.pred14.core <- predict(bsdDSM.sat2.xy2, preddata14_core, off.set = 40000)
write.csv(bsd.global.pred14.core, file="Results/BSD/core_only/bsd.pred14.core.csv")

# Predict over 2016 habitat
bsd.global.pred16.core <- predict(bsdDSM.sat2.xy2, preddata16_core, off.set = 40000)
write.csv(bsd.global.pred16.core, file="Results/BSD/core_only/bsd.pred16.core.csv")

# Predict over 2018 habitat
bsd.global.pred18.core <- predict(bsdDSM.sat2.xy2, preddata18_core, off.set = 40000)
write.csv(bsd.global.pred18.core, file="Results/BSD/core_only/bsd.pred18.core.csv")

# Predict over 2020 habitat
bsd.global.pred20.core <- predict(bsdDSM.sat2.xy2, preddata20_core, off.set = 40000)
write.csv(bsd.global.pred20.core, file="Results/BSD/core_only/bsd.pred20.core.csv")


# Create new dataframes
bsd.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = bsd.global.pred10.core)

bsd.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = bsd.global.pred11.core)

bsd.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = bsd.global.pred13.core)

bsd.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = bsd.global.pred14.core)

bsd.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = bsd.global.pred16.core)

bsd.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = bsd.global.pred18.core)

bsd.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = bsd.global.pred20.core)


## This creates a dataframe that can be plotted as a map
bsd.spdf.df_10.core <- abunPlotDF(bsd.df.Final10.core, pred.polys_200)
bsd.spdf.df_11.core <- abunPlotDF(bsd.df.Final11.core, pred.polys_200)
bsd.spdf.df_13.core <- abunPlotDF(bsd.df.Final13.core, pred.polys_200)
bsd.spdf.df_14.core <- abunPlotDF(bsd.df.Final14.core, pred.polys_200)
bsd.spdf.df_16.core <- abunPlotDF(bsd.df.Final16.core, pred.polys_200)
bsd.spdf.df_18.core <- abunPlotDF(bsd.df.Final18.core, pred.polys_200)
bsd.spdf.df_20.core <- abunPlotDF(bsd.df.Final20.core, pred.polys_200)

# save SPDFs
write.csv(bsd.spdf.df_10.core,file="Results/BSD/Plots/core_only/spdf/bsd.spdf.df_10.core.csv")
write.csv(bsd.spdf.df_11.core,file="Results/BSD/Plots/core_only/spdf/bsd.spdf.df_11.core.csv")
write.csv(bsd.spdf.df_13.core,file="Results/BSD/Plots/core_only/spdf/bsd.spdf.df_13.core.csv")
write.csv(bsd.spdf.df_14.core,file="Results/BSD/Plots/core_only/spdf/bsd.spdf.df_14.core.csv")
write.csv(bsd.spdf.df_16.core,file="Results/BSD/Plots/core_only/spdf/bsd.spdf.df_16.core.csv")
write.csv(bsd.spdf.df_18.core,file="Results/BSD/Plots/core_only/spdf/bsd.spdf.df_18.core.csv")
write.csv(bsd.spdf.df_20.core,file="Results/BSD/Plots/core_only/spdf/bsd.spdf.df_20.core.csv")




    ## Plotting continuous ####

# Load spatial dataframes
bsd.spdf.df_10.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_10.core.csv")
bsd.spdf.df_11.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_11.core.csv")
bsd.spdf.df_13.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_13.core.csv")
bsd.spdf.df_14.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_14.core.csv")
bsd.spdf.df_16.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_16.core.csv")
bsd.spdf.df_18.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_18.core.csv")
bsd.spdf.df_20.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_20.core.csv")


# greyscale plots
BSD_plot_10_core_gr <- GSplotFun(bsd.spdf.df_10.core, survey.area.core, "abundance","2010")
BSD_plot_11_core_gr <- GSplotFun(bsd.spdf.df_11.core, survey.area.core, "abundance")
BSD_plot_13_core_gr <- GSplotFun(bsd.spdf.df_13.core, survey.area.core, "abundance")
BSD_plot_14_core_gr <- GSplotFun(bsd.spdf.df_14.core, survey.area.core, "abundance")
BSD_plot_16_core_gr <- GSplotFun(bsd.spdf.df_16.core, survey.area.core, "abundance")
BSD_plot_18_core_gr <- GSplotFun(bsd.spdf.df_18.core, survey.area.core, "abundance")
BSD_plot_20_core_gr <- GSplotFun(bsd.spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(BSD_plot_10_core_gr,"Results/BSD/Plots/core_only/greyscale/BSD_plot_10_core_gr.png")
saveplot(BSD_plot_11_core_gr,"Results/BSD/Plots/core_only/greyscale/BSD_plot_11_core_gr.png")
saveplot(BSD_plot_13_core_gr,"Results/BSD/Plots/core_only/greyscale/BSD_plot_13_core_gr.png")
saveplot(BSD_plot_14_core_gr,"Results/BSD/Plots/core_only/greyscale/BSD_plot_14_core_gr.png")
saveplot(BSD_plot_16_core_gr,"Results/BSD/Plots/core_only/greyscale/BSD_plot_16_core_gr.png")
saveplot(BSD_plot_18_core_gr,"Results/BSD/Plots/core_only/greyscale/BSD_plot_18_core_gr.png")
saveplot(BSD_plot_20_core_gr,"Results/BSD/Plots/core_only/greyscale/BSD_plot_20_core_gr.png")

# colour plots
BSD_plot_10_core_col <- CLplotFun(bsd.spdf.df_10.core, survey.area.core, "abundance")
BSD_plot_11_core_col <- CLplotFun(bsd.spdf.df_11.core, survey.area.core, "abundance")
BSD_plot_13_core_col <- CLplotFun(bsd.spdf.df_13.core, survey.area.core, "abundance")
BSD_plot_14_core_col <- CLplotFun(bsd.spdf.df_14.core, survey.area.core, "abundance")
BSD_plot_16_core_col <- CLplotFun(bsd.spdf.df_16.core, survey.area.core, "abundance")
BSD_plot_18_core_col <- CLplotFun(bsd.spdf.df_18.core, survey.area.core, "abundance")
BSD_plot_20_core_col <- CLplotFun(bsd.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(BSD_plot_10_core_col,"Results/BSD/Plots/core_only/colour/BSD_plot_10_core_col.png")
saveplot(BSD_plot_11_core_col,"Results/BSD/Plots/core_only/colour/BSD_plot_11_core_col.png")
saveplot(BSD_plot_13_core_col,"Results/BSD/Plots/core_only/colour/BSD_plot_13_core_col.png")
saveplot(BSD_plot_14_core_col,"Results/BSD/Plots/core_only/colour/BSD_plot_14_core_col.png")
saveplot(BSD_plot_16_core_col,"Results/BSD/Plots/core_only/colour/BSD_plot_16_core_col.png")
saveplot(BSD_plot_18_core_col,"Results/BSD/Plots/core_only/colour/BSD_plot_18_core_col.png")
saveplot(BSD_plot_20_core_col,"Results/BSD/Plots/core_only/colour/BSD_plot_20_core_col.png")



## plot grids (abundance and variance)

# greyscale 
bsd_2yrs_gs <- 
  BSD_plot_10_core_gr + BSD_varplot_final10.core.bw  +
  BSD_plot_20_core_gr + BSD_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
bsd_2yrs_gs[[1]] <- bsd_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
bsd_2yrs_gs[[2]] <- bsd_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
bsd_2yrs_gs[[2]] <- bsd_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
bsd_2yrs_gs[[4]] <- bsd_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(bsd_2yrs_gs, "Results/BSD/Plots/core_only/greyscale/plot_grids/bsd_2yrs_gs.png")


    ## Plotting discrete bins ####
      # Add bins to SPDF - don't repeat ####

# this is the process of adding discrete bins to the abundance SPDFs. I will save them in a new folder so this only has to be done once

# Load original spatial dataframes
bsd.spdf.df_10.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_10.core.csv")
bsd.spdf.df_11.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_11.core.csv")
bsd.spdf.df_13.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_13.core.csv")
bsd.spdf.df_14.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_14.core.csv")
bsd.spdf.df_16.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_16.core.csv")
bsd.spdf.df_18.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_18.core.csv")
bsd.spdf.df_20.core <- read.csv("Results/BSD/Plots/core_only/spdf/bsd.spdf.df_20.core.csv")

# put spdf's into a list
dfs <- list(bsd.spdf.df_10.core,bsd.spdf.df_11.core,bsd.spdf.df_13.core,bsd.spdf.df_14.core,
            bsd.spdf.df_16.core,bsd.spdf.df_18.core,bsd.spdf.df_20.core)

# name the elements
names(dfs) <- c("bsd.spdf.df_10.core","bsd.spdf.df_11.core","bsd.spdf.df_13.core","bsd.spdf.df_14.core",
                "bsd.spdf.df_16.core","bsd.spdf.df_18.core","bsd.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunAbun)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(bsd.spdf.df_10.core,file="Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_10.core.csv")
write.csv(bsd.spdf.df_11.core,file="Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_11.core.csv")
write.csv(bsd.spdf.df_13.core,file="Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_13.core.csv")
write.csv(bsd.spdf.df_14.core,file="Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_14.core.csv")
write.csv(bsd.spdf.df_16.core,file="Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_16.core.csv")
write.csv(bsd.spdf.df_18.core,file="Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_18.core.csv")
write.csv(bsd.spdf.df_20.core,file="Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_20.core.csv")


      # Plotting ####

# Load spatial dataframes (with bins)
bsd.spdf.df_10.core <- read.csv("Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_10.core.csv")
bsd.spdf.df_11.core <- read.csv("Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_11.core.csv")
bsd.spdf.df_13.core <- read.csv("Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_13.core.csv")
bsd.spdf.df_14.core <- read.csv("Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_14.core.csv")
bsd.spdf.df_16.core <- read.csv("Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_16.core.csv")
bsd.spdf.df_18.core <- read.csv("Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_18.core.csv")
bsd.spdf.df_20.core <- read.csv("Results/BSD/Plots/core_only/spdf/bins/bsd.spdf.df_20.core.csv")

## plot greyscale
BSD_10_plot_bin_GS <- GSplotBin(bsd.spdf.df_10.core,survey.area.core,"Abundance","Relative abundance")
BSD_11_plot_bin_GS <- GSplotBin(bsd.spdf.df_11.core,survey.area.core,"Abundance","Relative abundance")
BSD_13_plot_bin_GS <- GSplotBin(bsd.spdf.df_13.core,survey.area.core,"Abundance","Relative abundance")
BSD_14_plot_bin_GS <- GSplotBin(bsd.spdf.df_14.core,survey.area.core,"Abundance","Relative abundance")
BSD_16_plot_bin_GS <- GSplotBin(bsd.spdf.df_16.core,survey.area.core,"Abundance","Relative abundance")
BSD_18_plot_bin_GS <- GSplotBin(bsd.spdf.df_18.core,survey.area.core,"Abundance","Relative abundance")
BSD_20_plot_bin_GS <- GSplotBin(bsd.spdf.df_20.core,"group2",survey.area.core,"Abundance","Relative abundance")


# save 
saveplot(BSD_10_plot_bin_GS,"Results/BSD/Plots/core_only/bins/BSD_10_plot_bin_GS.png")
saveplot(BSD_11_plot_bin_GS,"Results/BSD/Plots/core_only/bins/BSD_11_plot_bin_GS.png")
saveplot(BSD_13_plot_bin_GS,"Results/BSD/Plots/core_only/bins/BSD_13_plot_bin_GS.png")
saveplot(BSD_14_plot_bin_GS,"Results/BSD/Plots/core_only/bins/BSD_14_plot_bin_GS.png")
saveplot(BSD_16_plot_bin_GS,"Results/BSD/Plots/core_only/bins/BSD_16_plot_bin_GS.png")
saveplot(BSD_18_plot_bin_GS,"Results/BSD/Plots/core_only/bins/BSD_18_plot_bin_GS.png")
saveplot(BSD_20_plot_bin_GS,"Results/BSD/Plots/core_only/bins/BSD_20_plot_bin_GS.png")



## Variance estimation ####


# estimate variance
bsd.var.Final10.core <- varEstfun(preddata10_core, bsdDSM.sat2.xy2)
bsd.var.Final11.core <- varEstfun(preddata11_core, bsdDSM.sat2.xy2)
bsd.var.Final13.core <- varEstfun(preddata13_core, bsdDSM.sat2.xy2)
bsd.var.Final14.core <- varEstfun(preddata14_core, bsdDSM.sat2.xy2)
bsd.var.Final16.core <- varEstfun(preddata16_core, bsdDSM.sat2.xy2)
bsd.var.Final18.core <- varEstfun(preddata18_core, bsdDSM.sat2.xy2)
bsd.var.Final20.core <- varEstfun(preddata20_core, bsdDSM.sat2.xy2)

# save variance estimates
write.csv(bsd.var.Final10.core, file="Results/BSD/core_only/bsd.var10.core.csv")
write.csv(bsd.var.Final11.core, file="Results/BSD/core_only/bsd.var11.core.csv")
write.csv(bsd.var.Final13.core, file="Results/BSD/core_only/bsd.var13.core.csv")
write.csv(bsd.var.Final14.core, file="Results/BSD/core_only/bsd.var14.core.csv")
write.csv(bsd.var.Final16.core, file="Results/BSD/core_only/bsd.var16.core.csv")
write.csv(bsd.var.Final18.core, file="Results/BSD/core_only/bsd.var18.core.csv")
write.csv(bsd.var.Final20.core, file="Results/BSD/core_only/bsd.var20.core.csv")

# create spdf's for plotting
bsd.var.spdf.df_10.core <- varPlotDF(bsd.var.Final10.core, pred.polys_200)
bsd.var.spdf.df_11.core <- varPlotDF(bsd.var.Final11.core, pred.polys_200)
bsd.var.spdf.df_13.core <- varPlotDF(bsd.var.Final13.core, pred.polys_200)
bsd.var.spdf.df_14.core <- varPlotDF(bsd.var.Final14.core, pred.polys_200)
bsd.var.spdf.df_16.core <- varPlotDF(bsd.var.Final16.core, pred.polys_200)
bsd.var.spdf.df_18.core <- varPlotDF(bsd.var.Final18.core, pred.polys_200)
bsd.var.spdf.df_20.core <- varPlotDF(bsd.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(bsd.var.spdf.df_10.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_10.core.csv")
write.csv(bsd.var.spdf.df_11.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_11.core.csv")
write.csv(bsd.var.spdf.df_13.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_13.core.csv")
write.csv(bsd.var.spdf.df_14.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_14.core.csv")
write.csv(bsd.var.spdf.df_16.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_16.core.csv")
write.csv(bsd.var.spdf.df_18.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_18.core.csv")
write.csv(bsd.var.spdf.df_20.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_20.core.csv")


  
    ## Plotting variance ####

      # Calculate CV & add bins to SPDF - don't repeat ####

# Load spatial dataframes
bsd.var.spdf.df_10.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_10.core.csv")
bsd.var.spdf.df_11.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_11.core.csv")
bsd.var.spdf.df_13.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_13.core.csv")
bsd.var.spdf.df_14.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_14.core.csv")
bsd.var.spdf.df_16.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_16.core.csv")
bsd.var.spdf.df_18.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_18.core.csv")
bsd.var.spdf.df_20.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_20.core.csv")


# first need to calculate CV from the variance (note: need the abundance spdf's loaded)
bsd.var.spdf.df_10.core <- CVaddFun(bsd.spdf.df_10.core,bsd.var.spdf.df_10.core)
bsd.var.spdf.df_11.core <- CVaddFun(bsd.spdf.df_11.core,bsd.var.spdf.df_11.core)
bsd.var.spdf.df_13.core <- CVaddFun(bsd.spdf.df_13.core,bsd.var.spdf.df_13.core)
bsd.var.spdf.df_14.core <- CVaddFun(bsd.spdf.df_14.core,bsd.var.spdf.df_14.core)
bsd.var.spdf.df_16.core <- CVaddFun(bsd.spdf.df_16.core,bsd.var.spdf.df_16.core)
bsd.var.spdf.df_18.core <- CVaddFun(bsd.spdf.df_18.core,bsd.var.spdf.df_18.core)
bsd.var.spdf.df_20.core <- CVaddFun(bsd.spdf.df_20.core,bsd.var.spdf.df_20.core)


# Save the SPDF's with the CV value but no bins (for continuous plotting of the CV)
write.csv(bsd.var.spdf.df_10.core,file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_10.core.csv")
write.csv(bsd.var.spdf.df_11.core,file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_11.core.csv")
write.csv(bsd.var.spdf.df_13.core,file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_13.core.csv")
write.csv(bsd.var.spdf.df_14.core,file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_14.core.csv")
write.csv(bsd.var.spdf.df_16.core,file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_16.core.csv")
write.csv(bsd.var.spdf.df_18.core,file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_18.core.csv")
write.csv(bsd.var.spdf.df_20.core,file="Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_20.core.csv")


### add bins

## Quartiles

# put spdf's into a list
dfs <- list(bsd.var.spdf.df_10.core,bsd.var.spdf.df_11.core,bsd.var.spdf.df_13.core,bsd.var.spdf.df_14.core,
            bsd.var.spdf.df_16.core,bsd.var.spdf.df_18.core,bsd.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("bsd.var.spdf.df_10.core","bsd.var.spdf.df_11.core","bsd.var.spdf.df_13.core",
                "bsd.var.spdf.df_14.core","bsd.var.spdf.df_16.core","bsd.var.spdf.df_18.core",
                "bsd.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(bsd.var.spdf.df_10.core,file="Results/BSD/Plots/variance/core_only/spdf/bins/bsd.var.spdf.df_10.core.csv")
write.csv(bsd.var.spdf.df_11.core,file="Results/BSD/Plots/variance/core_only/spdf/bins/bsd.var.spdf.df_11.core.csv")
write.csv(bsd.var.spdf.df_13.core,file="Results/BSD/Plots/variance/core_only/spdf/bins/bsd.var.spdf.df_13.core.csv")
write.csv(bsd.var.spdf.df_14.core,file="Results/BSD/Plots/variance/core_only/spdf/bins/bsd.var.spdf.df_14.core.csv")
write.csv(bsd.var.spdf.df_16.core,file="Results/BSD/Plots/variance/core_only/spdf/bins/bsd.var.spdf.df_16.core.csv")
write.csv(bsd.var.spdf.df_18.core,file="Results/BSD/Plots/variance/core_only/spdf/bins/bsd.var.spdf.df_18.core.csv")
write.csv(bsd.var.spdf.df_20.core,file="Results/BSD/Plots/variance/core_only/spdf/bins/bsd.var.spdf.df_20.core.csv")




## custom bins

# change column name from group2 to CV
bsd.var.spdf.df_10.core <- bsd.var.spdf.df_10.core %>% dplyr::rename(CV=group2)
bsd.var.spdf.df_11.core <- bsd.var.spdf.df_11.core %>% dplyr::rename(CV=group2)
bsd.var.spdf.df_13.core <- bsd.var.spdf.df_13.core %>% dplyr::rename(CV=group2)
bsd.var.spdf.df_14.core <- bsd.var.spdf.df_14.core %>% dplyr::rename(CV=group2)
bsd.var.spdf.df_16.core <- bsd.var.spdf.df_16.core %>% dplyr::rename(CV=group2)
bsd.var.spdf.df_18.core <- bsd.var.spdf.df_18.core %>% dplyr::rename(CV=group2)
bsd.var.spdf.df_20.core <- bsd.var.spdf.df_20.core %>% dplyr::rename(CV=group2)


# put spdf's into a list
dfs <- list(bsd.var.spdf.df_10.core,bsd.var.spdf.df_11.core,bsd.var.spdf.df_13.core,bsd.var.spdf.df_14.core,
            bsd.var.spdf.df_16.core,bsd.var.spdf.df_18.core,bsd.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("bsd.var.spdf.df_10.core","bsd.var.spdf.df_11.core","bsd.var.spdf.df_13.core",
                "bsd.var.spdf.df_14.core","bsd.var.spdf.df_16.core","bsd.var.spdf.df_18.core",
                "bsd.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar2)

# split elements into original dataframes
list2env(dfs, globalenv())


# save the SPDFs with the custom bins
write.csv(bsd.var.spdf.df_10.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_10.core.csv")
write.csv(bsd.var.spdf.df_11.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_11.core.csv")
write.csv(bsd.var.spdf.df_13.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_13.core.csv")
write.csv(bsd.var.spdf.df_14.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_14.core.csv")
write.csv(bsd.var.spdf.df_16.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_16.core.csv")
write.csv(bsd.var.spdf.df_18.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_18.core.csv")
write.csv(bsd.var.spdf.df_20.core,
          file="Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_20.core.csv")


      # Discrete bins ####

# load spdfs (which has already had CV calculated and then put into bins)
bsd.var.spdf.df_10.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_10.core.csv")
bsd.var.spdf.df_11.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_11.core.csv")
bsd.var.spdf.df_13.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_13.core.csv")
bsd.var.spdf.df_14.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_14.core.csv")
bsd.var.spdf.df_16.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_16.core.csv")
bsd.var.spdf.df_18.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_18.core.csv")
bsd.var.spdf.df_20.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bins/custom/bsd.var.spdf.df_20.core.csv")

# make CV a factor and re-order. I have only done this for 2020 but you can copy the code for the other years if you need to
bsd.var.spdf.df_20.core$CV <- as.factor(bsd.var.spdf.df_20.core$CV)
bsd.var.spdf.df_20.core$CV <- factor(bsd.var.spdf.df_20.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))

# plot CV in bins
BSD_10_plot_bin_GS_var <- GSplotBin(bsd.var.spdf.df_10.core,survey.area.core,"Variance","CV")
BSD_11_plot_bin_GS_var <- GSplotBin(bsd.var.spdf.df_11.core,survey.area.core,"Variance","CV")
BSD_13_plot_bin_GS_var <- GSplotBin(bsd.var.spdf.df_13.core,survey.area.core,"Variance","CV")
BSD_14_plot_bin_GS_var <- GSplotBin(bsd.var.spdf.df_14.core,survey.area.core,"Variance","CV")
BSD_16_plot_bin_GS_var <- GSplotBin(bsd.var.spdf.df_16.core,survey.area.core,"Variance","CV")
BSD_18_plot_bin_GS_var <- GSplotBin(bsd.var.spdf.df_18.core,survey.area.core,"Variance","CV")
BSD_20_plot_bin_GS_var <- GSplotBin(bsd.var.spdf.df_20.core,"CV",survey.area.core,"Variance","CV")

# save 
saveplot(BSD_10_plot_bin_GS_var,"Results/BSD/Plots/variance/core_only/greyscale/bins/BSD_10_plot_bin_GS_var.png")
saveplot(BSD_11_plot_bin_GS_var,"Results/BSD/Plots/variance/core_only/greyscale/bins/BSD_11_plot_bin_GS_var.png")
saveplot(BSD_13_plot_bin_GS_var,"Results/BSD/Plots/variance/core_only/greyscale/bins/BSD_13_plot_bin_GS_var.png")
saveplot(BSD_14_plot_bin_GS_var,"Results/BSD/Plots/variance/core_only/greyscale/bins/BSD_14_plot_bin_GS_var.png")
saveplot(BSD_16_plot_bin_GS_var,"Results/BSD/Plots/variance/core_only/greyscale/bins/BSD_16_plot_bin_GS_var.png")
saveplot(BSD_18_plot_bin_GS_var,"Results/BSD/Plots/variance/core_only/greyscale/bins/BSD_18_plot_bin_GS_var.png")
saveplot(BSD_20_plot_bin_GS_var,"Results/BSD/Plots/variance/core_only/greyscale/bins/BSD_20_plot_bin_GS_var.png")



      # Continous ####


# Load spatial dataframes
bsd.var.spdf.df_10.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_10.core.csv")
bsd.var.spdf.df_11.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_11.core.csv")
bsd.var.spdf.df_13.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_13.core.csv")
bsd.var.spdf.df_14.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_14.core.csv")
bsd.var.spdf.df_16.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_16.core.csv")
bsd.var.spdf.df_18.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_18.core.csv")
bsd.var.spdf.df_20.core <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_20.core.csv")

# greyscale plots
BSD_varplot_final10.core.bw <- GSplotFun(bsd.var.spdf.df_10.core, survey.area.core, "group2", "2010")
BSD_varplot_final11.core.bw <- GSplotFun(bsd.var.spdf.df_11.core, survey.area.core, "variance")
BSD_varplot_final13.core.bw <- GSplotFun(bsd.var.spdf.df_13.core, survey.area.core, "variance")
BSD_varplot_final14.core.bw <- GSplotFun(bsd.var.spdf.df_14.core, survey.area.core, "variance")
BSD_varplot_final16.core.bw <- GSplotFun(bsd.var.spdf.df_16.core, survey.area.core, "variance")
BSD_varplot_final18.core.bw <- GSplotFun(bsd.var.spdf.df_18.core, survey.area.core, "variance")
BSD_varplot_final20.core.bw <- GSplotFun(bsd.var.spdf.df_20.core, survey.area.core, "variance", "2020")

# save greyscale
saveplot(BSD_varplot_final10.core.bw, "Results/BSD/Plots/variance/core_only/greyscale/2010_BSD_var.core.bw.png")
saveplot(BSD_varplot_final11.core.bw, "Results/BSD/Plots/variance/core_only/greyscale/2011_BSD_var.core.bw.png")
saveplot(BSD_varplot_final13.core.bw, "Results/BSD/Plots/variance/core_only/greyscale/2013_BSD_var.core.bw.png")
saveplot(BSD_varplot_final14.core.bw, "Results/BSD/Plots/variance/core_only/greyscale/2014_BSD_var.core.bw.png")
saveplot(BSD_varplot_final16.core.bw, "Results/BSD/Plots/variance/core_only/greyscale/2016_BSD_var.core.bw.png")
saveplot(BSD_varplot_final18.core.bw, "Results/BSD/Plots/variance/core_only/greyscale/2018_BSD_var.core.bw.png")
saveplot(BSD_varplot_final20.core.bw, "Results/BSD/Plots/variance/core_only/greyscale/2020_BSD_var.core.bw.png")

# colour plots
BSD_varplot_final10.core.col <- CLplotFun(bsd.var.spdf.df_10.core, survey.area.core, "variance")
BSD_varplot_final11.core.col <- CLplotFun(bsd.var.spdf.df_11.core, survey.area.core, "variance")
BSD_varplot_final13.core.col <- CLplotFun(bsd.var.spdf.df_13.core, survey.area.core, "variance")
BSD_varplot_final14.core.col <- CLplotFun(bsd.var.spdf.df_14.core, survey.area.core, "variance")
BSD_varplot_final16.core.col <- CLplotFun(bsd.var.spdf.df_16.core, survey.area.core, "variance")
BSD_varplot_final18.core.col <- CLplotFun(bsd.var.spdf.df_18.core, survey.area.core, "variance")
BSD_varplot_final20.core.col <- CLplotFun(bsd.var.spdf.df_20.core, survey.area.core, "variance")

# save colour
saveplot(BSD_varplot_final10.core.col, "Results/BSD/Plots/variance/core_only/colour/2010_BSD_var.core.col.png")
saveplot(BSD_varplot_final11.core.col, "Results/BSD/Plots/variance/core_only/colour/2011_BSD_var.core.col.png")
saveplot(BSD_varplot_final13.core.col, "Results/BSD/Plots/variance/core_only/colour/2013_BSD_var.core.col.png")
saveplot(BSD_varplot_final14.core.col, "Results/BSD/Plots/variance/core_only/colour/2014_BSD_var.core.col.png")
saveplot(BSD_varplot_final16.core.col, "Results/BSD/Plots/variance/core_only/colour/2016_BSD_var.core.col.png")
saveplot(BSD_varplot_final18.core.col, "Results/BSD/Plots/variance/core_only/colour/2018_BSD_var.core.col.png")
saveplot(BSD_varplot_final20.core.col, "Results/BSD/Plots/variance/core_only/colour/2020_BSD_var.core.col.png")


#### Red Muntjac ##############################################################
## Load data ####

# Observation data. Unique to species 
rmj_obsdata <- read.csv("Species_Data/RMJ/R Data/obsdata.csv", header = TRUE)
rmj_obsdata$object <- as.factor(rmj_obsdata$object)
rmj_obsdata$Sample.Label <- as.factor(rmj_obsdata$Sample.Label)
#rmj_obsdata$year <- as.factor(rmj_obsdata$year)
str(rmj_obsdata)
head(rmj_obsdata)

# Transect data. Unique to species
rmj_distdata <- read.csv("Species_Data/RMJ/R Data/distdata.csv", header = TRUE)
rmj_distdata$object <- as.factor(rmj_distdata$object)
rmj_distdata$NameObserver <- as.factor(rmj_distdata$NameObserver)
rmj_distdata$transect <- as.factor(rmj_distdata$transect)
rmj_distdata$year <- as.factor(rmj_distdata$year)
rmj_distdata$date <- as.Date(rmj_distdata$date, format = "%d/%m/%Y")
str(rmj_distdata)
head(rmj_distdata)

# Remove extreme outliers
rmj_distdata <- rmj_distdata %>% filter(distance < 100)

# check if there are any observations on transect 20
rmj_distdata[rmj_distdata$transect=="20",]
# no


## Plot the covariates across the grid, with group sizes ####

# Warning - plots take a few minutes to run

# habitat
plot_RMJobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata200$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=rmj_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/RMJ/plot_RMJobs_habitat.png", plot_RMJobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_RMJobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys_200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=rmj_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/RMJ/plot_RMJobs_dstWater.png", plot_RMJobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_RMJobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=rmj_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/RMJ/plot_RMJobs_dstStlmnt.png", plot_RMJobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_RMJobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=rmj_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/RMJ/plot_RMJobs_dstRoad.png", plot_RMJobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_RMJobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=rmj_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/RMJ/plot_RMJobs_dstBorder.png", plot_RMJobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_RMJobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=rmj_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/RMJ/plot_RMJobs_dstStation.png", plot_RMJobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_RMJobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=rmj_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/RMJ/plot_RMJobs_dstELC.png", plot_RMJobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_RMJobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys_200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=rmj_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/RMJ/plot_RMJobs_elev.png", plot_RMJobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)

## Exploratory plots & linear models ####

## Checking that there are no large gaps in the range of variables for muntjac observations

# subset segdata to get only the segments with RMJ observations
rmj_varcheck <- segdata[match(rmj_obsdata$Sample.Label,segdata$Sample.Label), ]

# habitat - fine
plot(segdata$habitat)
plot(rmj_varcheck$habitat)

# dstWater - fine
hist(segdata$dstWater)
hist(rmj_varcheck$dstWater)

# dstStlmnt - fine
hist(segdata$dstStlmnt)
hist(rmj_varcheck$dstStlmnt)

# dstRoad - fine
hist(segdata$dstRoad)
hist(rmj_varcheck$dstRoad)

# dstBorder - fine
hist(segdata$dstBorder)
hist(rmj_varcheck$dstBorder)

# dstStation - fine
hist(segdata$dstStation)
hist(rmj_varcheck$dstStation)

# elevation - fine
hist(segdata$elevation)
hist(rmj_varcheck$elevation)


## Histograms

# Distance 
rmj_h1 <- ggplot(rmj_distdata, aes(distance))+ geom_histogram(binwidth = 1)
rmj_h2 <- ggplot(rmj_distdata, aes(distance))+ geom_histogram(binwidth = 5)
rmj_h3 <- ggplot(rmj_distdata, aes(distance))+ geom_histogram(binwidth = 10)
rmj_h4 <- ggplot(rmj_distdata, aes(distance))+ geom_histogram(binwidth = 15)
rmj_h5 <- ggplot(rmj_distdata, aes(distance))+ geom_histogram(binwidth = 20)
rmj_h6 <- ggplot(rmj_distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(rmj_h1,rmj_h2,rmj_h3,rmj_h4,rmj_h5,rmj_h6)
# evidence of evasive movement, but not as bad as GPF or BSD. Histogram loks pretty good


# cluster size, observer, habitat, year, month, transect
rmj_h7 <- ggplot(rmj_distdata, aes(size))+geom_histogram(binwidth = 0.5)
rmj_h8 <- ggplot(rmj_distdata, aes(NameObserver))+geom_histogram(stat="count")
rmj_h9 <- ggplot(rmj_distdata, aes(habitat))+geom_histogram(stat="count")
rmj_h10 <- ggplot(rmj_distdata, aes(year))+geom_histogram(stat="count")
rmj_h11 <- ggplot(rmj_distdata, aes(month))+geom_histogram(stat="count")
rmj_h12 <- ggplot(rmj_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(rmj_h7,rmj_h8,rmj_h9,rmj_h10,rmj_h11,rmj_h12)
# Vast majority of observations are of single animals, and so estimating abundance should be for individuals not groups. Open habitat most common, but still plenty in Dense. Number of observations drop significantly in 2016 and 2018.  

## Plots of distance against variables

plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

rmj_d1 <- ggplot(rmj_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
rmj_d2 <- ggplot(rmj_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
rmj_d3 <- ggplot(rmj_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
rmj_d4 <- ggplot(rmj_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
rmj_d5 <- ggplot(rmj_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(rmj_d1,rmj_d2,rmj_d3,rmj_d4,rmj_d5)
# distances of observations slightly higher in open forest, as would be expected. Very little information about size bias (differences in distances based on group size) because most obs are either 1 or 2. 


## Plots of cluster size against variables
rmj_s1 <- ggplot(rmj_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
rmj_s2 <- ggplot(rmj_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
rmj_s3 <- ggplot(rmj_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
rmj_s4 <- ggplot(rmj_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
rmj_s5 <- ggplot(rmj_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(rmj_s1,rmj_s2,rmj_s3,rmj_s4,rmj_s5)
# All observations of group sizes of more than 2 are in open forest.  Makes sense, but there are so few of those observations we can't really draw any conclusions.   

## Linear models

# group size ~ distance
newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=rmj_distdata)
plot(rmj_distdata$size~rmj_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# No support for size bias


## Estimating detection function ####

### The below DF model selection is obsolete as I am now using the same DF model as decided upon in the CDS. I have kept the below in just for reference.  The CDS DF is Hn cos with bins (see binning explanation at top of CDS script)

rmjDF.hn.cos.bin <- ds(rmj_distdata, truncation=50, key="hn", adjustment ="cos",
                       cutpoints = c(0,7,14,21,28,35,42,50))



## Fitting a spatial model ####

# There are covariates in the DF so we have to use Nhat

# I am NOT setting group = TRUE which means abundance of individuals rather than groups will be estimated. This is because RMJ are not gregarious, and the vast majority of observations have a group size of 1. Therefore we will be estimating the abundance of individuals.

# I am not including dstELC in the models for RMJ, as I don't think it is an appropriate variable, especially now that we are going to be predicting into the buffer zone

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
rmj_obsdata <- rmj_obsdata %>% filter(distance <= 50)

  ## Quasipoisson response ####


# Saturated model with all covariates. 
rmjDSM.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(dstRoad,bs="ts")+ 
                         s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                         habitat,
                  rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(rmjDSM.sat)
par(mfrow=c(2,3))
plot(rmjDSM.sat, scale = 0)
gam.check(rmjDSM.sat)
# DE = 12.6. dstWater, dstRoad not significant.  All plots overfitted except dstBorder. 

# reduce k for all except dstBorder
rmjDSM.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts", k=5)+ s(dstStlmnt,bs="ts", k=5)+ 
                     s(dstRoad,bs="ts", k=5)+ s(dstBorder,bs="ts")+ 
                     s(dstStation,bs="ts", k=5)+ s(elevation,bs="ts", k=5)+ 
                    habitat,
                  rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(rmjDSM.sat2)
par(mfrow=c(2,3))
plot(rmjDSM.sat2, scale = 0)
gam.check(rmjDSM.sat2)
# DE = 11.1. All sig except dstRoad

# remove dstRoad
rmjDSM.sat3 <- dsm(Nhat ~ s(dstWater,bs="ts", k=5)+ s(dstStlmnt,bs="ts", k=5)+ 
                    s(dstBorder,bs="ts")+ s(dstStation,bs="ts", k=5)+ 
                    s(elevation,bs="ts", k=5)+ habitat,
                   rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area)
summary(rmjDSM.sat3)
par(mfrow=c(2,3))
plot(rmjDSM.sat3, scale = 0)
gam.check(rmjDSM.sat3)
# DE = 11. All terms sig.

# try increase k slightly
rmjDSM.sat4 <- dsm(Nhat ~ s(dstWater,bs="ts", k=7)+ s(dstStlmnt,bs="ts", k=7)+ 
                     s(dstBorder,bs="ts")+ s(dstStation,bs="ts", k=7)+ 
                     s(elevation,bs="ts", k=7)+ habitat,
                   rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area)
summary(rmjDSM.sat4)
par(mfrow=c(2,3))
plot(rmjDSM.sat4, scale = 0)
gam.check(rmjDSM.sat4)
# DE=11.5. dstWater no longer sig. Overfitted plots


# rmjDSM.sat3 is the selected QP model

  ## Tweedie response ####


rmjDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                             s(dstRoad,bs="ts")+ s(dstBorder,bs="ts")+
                             s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                             habitat,
                  rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(rmjDSM.tw.sat)
dev.off()
par(mfrow=c(3,2))
gam.check(rmjDSM.tw.sat)
dev.off()
par(mfrow=c(2,3))
plot(rmjDSM.tw.sat, scale=0)
# AIC = 15090, DE = 10.7.  dstRoad not sig.  dstStlmnt and maybe dstWater looking parametric. plots look good for dstBorder, dstStation, and elevation.

# increase k for dstWater, dstStlmnt, dstRoad
rmjDSM.tw.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+ 
                       s(dstRoad,bs="ts",k=15)+ s(dstBorder,bs="ts")+
                       s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                       habitat,
                     rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                     family = tw(), engine = "gam",
                     segment.area = segdata$Sample.Area)
summary(rmjDSM.tw.sat2)
dev.off()
par(mfrow=c(3,2))
gam.check(rmjDSM.tw.sat2)
dev.off()
par(mfrow=c(2,3))
plot(rmjDSM.tw.sat2, scale=0)
# AIC = 15090, DE = 10.7. very minor differences. 

# remove dstRoad
rmjDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+ 
                        s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ 
                        s(elevation,bs="ts")+ habitat,
                      rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                      family = tw(), engine = "gam",
                      segment.area = segdata$Sample.Area)
summary(rmjDSM.tw.3)
dev.off()
par(mfrow=c(3,2))
gam.check(rmjDSM.tw.3)
dev.off()
par(mfrow=c(2,3))
plot(rmjDSM.tw.3, scale=0)
# AIC = 15329, DE = 11.5. All terms except dstWater sig. dstStlmnt is overfit

# remove dstWater and reduce k for dstStlmnt
rmjDSM.tw.4 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=10)+ s(dstBorder,bs="ts")+ 
                     s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                     habitat,
                   rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area)
summary(rmjDSM.tw.4)
dev.off()
par(mfrow=c(3,2))
gam.check(rmjDSM.tw.4)
dev.off()
par(mfrow=c(2,2))
plot(rmjDSM.tw.4, scale=0)
# AIC = 15317, DE = 11.5. Improved AIC. All terms sig. Happy with the plots

# rmjDSM.tw.4 is the best tweedie model


  ## Negative binomial response ####

rmjDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                             s(dstRoad,bs="ts")+ s(dstBorder,bs="ts")+
                             s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                             habitat,
                  rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(rmjDSM.nb.sat)
dev.off()
par(mfrow=c(2,2))
gam.check(rmjDSM.nb.sat)
dev.off()
par(mfrow=c(2,3))
plot(rmjDSM.nb.sat, scale=0)
# AIC = 7385, DE = 15.1. dstWater, dstRoad, elevation not sig. dstStlmnt looks parametric

# increase k for dstWater, dstStlment, dstRoad, elevation
rmjDSM.nb.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=10)+ s(dstStlmnt,bs="ts",k=10)+ 
                       s(dstRoad,bs="ts",k=10)+ s(dstBorder,bs="ts")+
                       s(dstStation,bs="ts")+ s(elevation,bs="ts",k=10)+ 
                       habitat,
                     rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                     family = nb(), engine = "gam",
                     segment.area = segdata$Sample.Area)
summary(rmjDSM.nb.sat2)
dev.off()
par(mfrow=c(2,2))
gam.check(rmjDSM.nb.sat2)
dev.off()
par(mfrow=c(2,3))
plot(rmjDSM.nb.sat2, scale=0)
# AIC = 7385, DE = 15.1. No difference at all

# remove dstRoad
rmjDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=10)+ s(dstStlmnt,bs="ts",k=10)+ 
                        s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ 
                        s(elevation,bs="ts",k=10)+ habitat,
                      rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                      family = nb(), engine = "gam",
                      segment.area = segdata$Sample.Area)
summary(rmjDSM.nb.3)
dev.off()
par(mfrow=c(2,2))
gam.check(rmjDSM.nb.3)
dev.off()
par(mfrow=c(2,3))
plot(rmjDSM.nb.3, scale=0)
# AIC = 7385, DE = 15.1.

# remove elevation
rmjDSM.nb.4 <- dsm(Nhat ~ s(dstWater,bs="ts",k=10)+ s(dstStlmnt,bs="ts",k=10)+ 
                     s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ 
                     habitat,
                   rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area)
summary(rmjDSM.nb.4)
dev.off()
par(mfrow=c(2,2))
gam.check(rmjDSM.nb.4)
dev.off()
par(mfrow=c(2,2))
plot(rmjDSM.nb.4, scale=0)
# AIC = 7385, DE = 15.1.

# remove dsWater
rmjDSM.nb.5 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=10)+ s(dstBorder,bs="ts")+ 
                     s(dstStation,bs="ts")+ habitat,
                   rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area)
summary(rmjDSM.nb.5)
dev.off()
par(mfrow=c(2,2))
gam.check(rmjDSM.nb.5)
dev.off()
par(mfrow=c(2,2))
plot(rmjDSM.nb.5, scale=0)
# AIC = 7385, DE = 15.1. All terms sig. 

# increase k for dstStlment
rmjDSM.nb.6 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=15)+ s(dstBorder,bs="ts")+ 
                     s(dstStation,bs="ts")+ habitat,
                   rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area)
summary(rmjDSM.nb.6)
dev.off()
par(mfrow=c(2,2))
gam.check(rmjDSM.nb.6)
dev.off()
par(mfrow=c(2,2))
plot(rmjDSM.nb.6, scale=0, xlab="Distance to ranger station (m)")
# AIC = 7363, DE = 16.3. All terms sig. dstStlmnt no longer linear. 

# increase k for dstBorder and dstStation
rmjDSM.nb.7 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=15)+ s(dstBorder,bs="ts", k=15)+ 
                     s(dstStation,bs="ts", k=15)+ habitat,
                   rmjDF.hn.cos.bin, segdata, rmj_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area)
summary(rmjDSM.nb.7)
dev.off()
par(mfrow=c(2,2))
gam.check(rmjDSM.nb.7)
dev.off()
par(mfrow=c(2,2))
plot(rmjDSM.nb.7, scale=0)
# AIC = 7363, DE = 16.4. Tried various k settings and no real difference


### rmjDSM.nb.6 is the best model 

## Model selection ####

# The best quasipoisson model is rmjDSM.sat3
# The best tweedie model is rmjDSM.tw.4
# The best negative binomial model is rmjDSM.nb.6

summary(rmjDSM.sat3) # 11 %
summary(rmjDSM.tw.4) # 11.5 %
summary(rmjDSM.nb.6) # 16.3 %

# Negative binomial model has the highest DE. 
par(mfrow=c(2,2))
gam.check(rmjDSM.sat3)
gam.check(rmjDSM.tw.4)
gam.check(rmjDSM.nb.6)

# The Q-Q plot for the negative binomial model is much better than the other two.

rmjDSM.tw.4$aic
rmjDSM.nb.6$aic
# AIC for NB is better

anova(rmjDSM.sat3,rmjDSM.nb.6,test="Chisq")
# residual deviance is lower for NB model

# rmjDSM.nb.6 is the best global model for RMJ

## Autocorrelation ####

# rmjDSM.nb.6 is the best global model for RMJ
par(mfrow=c(1,1))
dsm.cor(rmjDSM.nb.6, max.lag=15, Segment.Label="Sample.Label")


# There is no significant autocorrelation


## Abundance estimation ####


# Predict over 2010 habitat
rmj.global.pred10.core <- predict(rmjDSM.nb.6, preddata10_core, off.set = 40000)
write.csv(rmj.global.pred10.core, file="Results/RMJ/core_only/rmj.pred10.core.csv")

# Predict over 2011 habitat
rmj.global.pred11.core <- predict(rmjDSM.nb.6, preddata11_core, off.set = 40000)
write.csv(rmj.global.pred11.core, file="Results/RMJ/core_only/rmj.pred11.core.csv")

# Predict over 2013 habitat
rmj.global.pred13.core <- predict(rmjDSM.nb.6, preddata13_core, off.set = 40000)
write.csv(rmj.global.pred13.core, file="Results/RMJ/core_only/rmj.pred13.core.csv")

# Predict over 2014 habitat
rmj.global.pred14.core <- predict(rmjDSM.nb.6, preddata14_core, off.set = 40000)
write.csv(rmj.global.pred14.core, file="Results/RMJ/core_only/rmj.pred14.core.csv")

# Predict over 2016 habitat
rmj.global.pred16.core <- predict(rmjDSM.nb.6, preddata16_core, off.set = 40000)
write.csv(rmj.global.pred16.core, file="Results/RMJ/core_only/rmj.pred16.core.csv")

# Predict over 2018 habitat
rmj.global.pred18.core <- predict(rmjDSM.nb.6, preddata18_core, off.set = 40000)
write.csv(rmj.global.pred18.core, file="Results/RMJ/core_only/rmj.pred18.core.csv")

# Predict over 2020 habitat
rmj.global.pred20.core <- predict(rmjDSM.nb.6, preddata20_core, off.set = 40000)
write.csv(rmj.global.pred20.core, file="Results/RMJ/core_only/rmj.pred20.core.csv")


# Create new dataframes for plotting
rmj.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = rmj.global.pred10.core)

rmj.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = rmj.global.pred11.core)

rmj.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = rmj.global.pred13.core)

rmj.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = rmj.global.pred14.core)

rmj.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = rmj.global.pred16.core)

rmj.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = rmj.global.pred18.core)

rmj.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = rmj.global.pred20.core)


## This creates a dataframe that can be plotted as a map
rmj.spdf.df_10.core <- abunPlotDF(rmj.df.Final10.core, pred.polys_200)
rmj.spdf.df_11.core <- abunPlotDF(rmj.df.Final11.core, pred.polys_200)
rmj.spdf.df_13.core <- abunPlotDF(rmj.df.Final13.core, pred.polys_200)
rmj.spdf.df_14.core <- abunPlotDF(rmj.df.Final14.core, pred.polys_200)
rmj.spdf.df_16.core <- abunPlotDF(rmj.df.Final16.core, pred.polys_200)
rmj.spdf.df_18.core <- abunPlotDF(rmj.df.Final18.core, pred.polys_200)
rmj.spdf.df_20.core <- abunPlotDF(rmj.df.Final20.core, pred.polys_200)

# save SPDFs
write.csv(rmj.spdf.df_10.core,file="Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_10.core.csv")
write.csv(rmj.spdf.df_11.core,file="Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_11.core.csv")
write.csv(rmj.spdf.df_13.core,file="Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_13.core.csv")
write.csv(rmj.spdf.df_14.core,file="Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_14.core.csv")
write.csv(rmj.spdf.df_16.core,file="Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_16.core.csv")
write.csv(rmj.spdf.df_18.core,file="Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_18.core.csv")
write.csv(rmj.spdf.df_20.core,file="Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_20.core.csv")




    ## Plotting continuous ####

## Load spatial dataframes
rmj.spdf.df_10.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_10.core.csv")
rmj.spdf.df_11.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_11.core.csv")
rmj.spdf.df_13.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_13.core.csv")
rmj.spdf.df_14.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_14.core.csv")
rmj.spdf.df_16.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_16.core.csv")
rmj.spdf.df_18.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_18.core.csv")
rmj.spdf.df_20.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_20.core.csv")

# greyscale plots
RMJ_plot_10_core_gr <- GSplotFun(rmj.spdf.df_10.core, survey.area.core, "abundance", "2010")
RMJ_plot_11_core_gr <- GSplotFun(rmj.spdf.df_11.core, survey.area.core, "abundance")
RMJ_plot_13_core_gr <- GSplotFun(rmj.spdf.df_13.core, survey.area.core, "abundance")
RMJ_plot_14_core_gr <- GSplotFun(rmj.spdf.df_14.core, survey.area.core, "abundance")
RMJ_plot_16_core_gr <- GSplotFun(rmj.spdf.df_16.core, survey.area.core, "abundance")
RMJ_plot_18_core_gr <- GSplotFun(rmj.spdf.df_18.core, survey.area.core, "abundance")
RMJ_plot_20_core_gr <- GSplotFun(rmj.spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(RMJ_plot_10_core_gr,"Results/RMJ/Plots/core_only/greyscale/RMJ_plot_10_core_gr.png")
saveplot(RMJ_plot_11_core_gr,"Results/RMJ/Plots/core_only/greyscale/RMJ_plot_11_core_gr.png")
saveplot(RMJ_plot_13_core_gr,"Results/RMJ/Plots/core_only/greyscale/RMJ_plot_13_core_gr.png")
saveplot(RMJ_plot_14_core_gr,"Results/RMJ/Plots/core_only/greyscale/RMJ_plot_14_core_gr.png")
saveplot(RMJ_plot_16_core_gr,"Results/RMJ/Plots/core_only/greyscale/RMJ_plot_16_core_gr.png")
saveplot(RMJ_plot_18_core_gr,"Results/RMJ/Plots/core_only/greyscale/RMJ_plot_18_core_gr.png")
saveplot(RMJ_plot_20_core_gr,"Results/RMJ/Plots/core_only/greyscale/RMJ_plot_20_core_gr.png")

# colour plots
RMJ_plot_10_core_col <- CLplotFun(rmj.spdf.df_10.core, survey.area.core, "abundance")
RMJ_plot_11_core_col <- CLplotFun(rmj.spdf.df_11.core, survey.area.core, "abundance")
RMJ_plot_13_core_col <- CLplotFun(rmj.spdf.df_13.core, survey.area.core, "abundance")
RMJ_plot_14_core_col <- CLplotFun(rmj.spdf.df_14.core, survey.area.core, "abundance")
RMJ_plot_16_core_col <- CLplotFun(rmj.spdf.df_16.core, survey.area.core, "abundance")
RMJ_plot_18_core_col <- CLplotFun(rmj.spdf.df_18.core, survey.area.core, "abundance")
RMJ_plot_20_core_col <- CLplotFun(rmj.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(RMJ_plot_10_core_col,"Results/RMJ/Plots/core_only/colour/RMJ_plot_10_core_col.png")
saveplot(RMJ_plot_11_core_col,"Results/RMJ/Plots/core_only/colour/RMJ_plot_11_core_col.png")
saveplot(RMJ_plot_13_core_col,"Results/RMJ/Plots/core_only/colour/RMJ_plot_13_core_col.png")
saveplot(RMJ_plot_14_core_col,"Results/RMJ/Plots/core_only/colour/RMJ_plot_14_core_col.png")
saveplot(RMJ_plot_16_core_col,"Results/RMJ/Plots/core_only/colour/RMJ_plot_16_core_col.png")
saveplot(RMJ_plot_18_core_col,"Results/RMJ/Plots/core_only/colour/RMJ_plot_18_core_col.png")
saveplot(RMJ_plot_20_core_col,"Results/RMJ/Plots/core_only/colour/RMJ_plot_20_core_col.png")



## plot grids (abundance and variance)

# greyscale 
rmj_2yrs_gs <- 
  RMJ_plot_10_core_gr + RMJ_varplot_final10.core.bw  +
  RMJ_plot_20_core_gr + RMJ_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
rmj_2yrs_gs[[1]] <- rmj_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
rmj_2yrs_gs[[2]] <- rmj_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
rmj_2yrs_gs[[2]] <- rmj_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
rmj_2yrs_gs[[4]] <- rmj_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(rmj_2yrs_gs, "Results/RMJ/Plots/core_only/greyscale/plot_grids/rmj_2yrs_gs.png")


    ## Plotting discrete bins ####
      # Add bins to SPDF - don't repeat ####

# this is the process of adding discrete bins to the abundance SPDFs. I will save them in a new folder so this only has to be done once

# Load original spatial dataframes
rmj.spdf.df_10.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_10.core.csv")
rmj.spdf.df_11.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_11.core.csv")
rmj.spdf.df_13.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_13.core.csv")
rmj.spdf.df_14.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_14.core.csv")
rmj.spdf.df_16.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_16.core.csv")
rmj.spdf.df_18.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_18.core.csv")
rmj.spdf.df_20.core <- read.csv("Results/RMJ/Plots/core_only/spdf/rmj.spdf.df_20.core.csv")

# put spdf's into a list
dfs <- list(rmj.spdf.df_10.core,rmj.spdf.df_11.core,rmj.spdf.df_13.core,rmj.spdf.df_14.core,
            rmj.spdf.df_16.core,rmj.spdf.df_18.core,rmj.spdf.df_20.core)

# name the elements
names(dfs) <- c("rmj.spdf.df_10.core","rmj.spdf.df_11.core","rmj.spdf.df_13.core","rmj.spdf.df_14.core",
                "rmj.spdf.df_16.core","rmj.spdf.df_18.core","rmj.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunAbun)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(rmj.spdf.df_10.core,file="Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_10.core.csv")
write.csv(rmj.spdf.df_11.core,file="Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_11.core.csv")
write.csv(rmj.spdf.df_13.core,file="Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_13.core.csv")
write.csv(rmj.spdf.df_14.core,file="Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_14.core.csv")
write.csv(rmj.spdf.df_16.core,file="Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_16.core.csv")
write.csv(rmj.spdf.df_18.core,file="Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_18.core.csv")
write.csv(rmj.spdf.df_20.core,file="Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_20.core.csv")


      # Plotting ####

# Load spatial dataframes (with bins)
rmj.spdf.df_10.core <- read.csv("Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_10.core.csv")
rmj.spdf.df_11.core <- read.csv("Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_11.core.csv")
rmj.spdf.df_13.core <- read.csv("Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_13.core.csv")
rmj.spdf.df_14.core <- read.csv("Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_14.core.csv")
rmj.spdf.df_16.core <- read.csv("Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_16.core.csv")
rmj.spdf.df_18.core <- read.csv("Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_18.core.csv")
rmj.spdf.df_20.core <- read.csv("Results/RMJ/Plots/core_only/spdf/bins/rmj.spdf.df_20.core.csv")

# change group2 (abundance) to factor and re-order. I've only done it for 2020 here but you can copy the code for the other years if you need to
rmj.spdf.df_20.core$group2 <- as.factor(rmj.spdf.df_20.core$group2)
rmj.spdf.df_20.core$group2 <- factor(rmj.spdf.df_20.core$group2, levels=c("High","Medium","Low","Very low"))

## plot greyscale
RMJ_10_plot_bin_GS <- GSplotBin(rmj.spdf.df_10.core,"group2",survey.area.core,"Abundance","Relative abundance")
RMJ_11_plot_bin_GS <- GSplotBin(rmj.spdf.df_11.core,"group2",survey.area.core,"Abundance","Relative abundance")
RMJ_13_plot_bin_GS <- GSplotBin(rmj.spdf.df_13.core,"group2",survey.area.core,"Abundance","Relative abundance")
RMJ_14_plot_bin_GS <- GSplotBin(rmj.spdf.df_14.core,"group2",survey.area.core,"Abundance","Relative abundance")
RMJ_16_plot_bin_GS <- GSplotBin(rmj.spdf.df_16.core,"group2",survey.area.core,"Abundance","Relative abundance")
RMJ_18_plot_bin_GS <- GSplotBin(rmj.spdf.df_18.core,"group2",survey.area.core,"Abundance","Relative abundance")
RMJ_20_plot_bin_GS <- GSplotBin(rmj.spdf.df_20.core,"group2",survey.area.core,"Abundance","Relative abundance")


# save 
saveplot(RMJ_10_plot_bin_GS,"Results/RMJ/Plots/core_only/bins/RMJ_10_plot_bin_GS.png")
saveplot(RMJ_11_plot_bin_GS,"Results/RMJ/Plots/core_only/bins/RMJ_11_plot_bin_GS.png")
saveplot(RMJ_13_plot_bin_GS,"Results/RMJ/Plots/core_only/bins/RMJ_13_plot_bin_GS.png")
saveplot(RMJ_14_plot_bin_GS,"Results/RMJ/Plots/core_only/bins/RMJ_14_plot_bin_GS.png")
saveplot(RMJ_16_plot_bin_GS,"Results/RMJ/Plots/core_only/bins/RMJ_16_plot_bin_GS.png")
saveplot(RMJ_18_plot_bin_GS,"Results/RMJ/Plots/core_only/bins/RMJ_18_plot_bin_GS.png")
saveplot(RMJ_20_plot_bin_GS,"Results/RMJ/Plots/core_only/bins/RMJ_20_plot_bin_GS.png")



## Variance estimation ####
 
# estimate variance
rmj.var.Final10.core <- varEstfun(preddata10_core, rmjDSM.nb.6)
rmj.var.Final11.core <- varEstfun(preddata11_core, rmjDSM.nb.6)
rmj.var.Final13.core <- varEstfun(preddata13_core, rmjDSM.nb.6)
rmj.var.Final14.core <- varEstfun(preddata14_core, rmjDSM.nb.6)
rmj.var.Final16.core <- varEstfun(preddata16_core, rmjDSM.nb.6)
rmj.var.Final18.core <- varEstfun(preddata18_core, rmjDSM.nb.6)
rmj.var.Final20.core <- varEstfun(preddata20_core, rmjDSM.nb.6)

# save variance estimates
write.csv(rmj.var.Final10.core, file="Results/RMJ/core_only/rmj.var10.core.csv")
write.csv(rmj.var.Final11.core, file="Results/RMJ/core_only/rmj.var11.core.csv")
write.csv(rmj.var.Final13.core, file="Results/RMJ/core_only/rmj.var13.core.csv")
write.csv(rmj.var.Final14.core, file="Results/RMJ/core_only/rmj.var14.core.csv")
write.csv(rmj.var.Final16.core, file="Results/RMJ/core_only/rmj.var16.core.csv")
write.csv(rmj.var.Final18.core, file="Results/RMJ/core_only/rmj.var18.core.csv")
write.csv(rmj.var.Final20.core, file="Results/RMJ/core_only/rmj.var20.core.csv")

# create spdf's for plotting
rmj.var.spdf.df_10.core <- varPlotDF(rmj.var.Final10.core, pred.polys_200)
rmj.var.spdf.df_11.core <- varPlotDF(rmj.var.Final11.core, pred.polys_200)
rmj.var.spdf.df_13.core <- varPlotDF(rmj.var.Final13.core, pred.polys_200)
rmj.var.spdf.df_14.core <- varPlotDF(rmj.var.Final14.core, pred.polys_200)
rmj.var.spdf.df_16.core <- varPlotDF(rmj.var.Final16.core, pred.polys_200)
rmj.var.spdf.df_18.core <- varPlotDF(rmj.var.Final18.core, pred.polys_200)
rmj.var.spdf.df_20.core <- varPlotDF(rmj.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(rmj.var.spdf.df_10.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_10.core.csv")
write.csv(rmj.var.spdf.df_11.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_11.core.csv")
write.csv(rmj.var.spdf.df_13.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_13.core.csv")
write.csv(rmj.var.spdf.df_14.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_14.core.csv")
write.csv(rmj.var.spdf.df_16.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_16.core.csv")
write.csv(rmj.var.spdf.df_18.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_18.core.csv")
write.csv(rmj.var.spdf.df_20.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_20.core.csv")


    
    ## Plotting variance ####
      # Calculate CV & add bins to SPDF - don't repeat ####

# Load spatial dataframes
rmj.var.spdf.df_10.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_10.core.csv")
rmj.var.spdf.df_11.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_11.core.csv")
rmj.var.spdf.df_13.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_13.core.csv")
rmj.var.spdf.df_14.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_14.core.csv")
rmj.var.spdf.df_16.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_16.core.csv")
rmj.var.spdf.df_18.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_18.core.csv")
rmj.var.spdf.df_20.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_20.core.csv")


# first need to calculate CV from the variance (note: need the abundance spdf's loaded)
rmj.var.spdf.df_10.core <- CVaddFun(rmj.spdf.df_10.core,rmj.var.spdf.df_10.core)
rmj.var.spdf.df_11.core <- CVaddFun(rmj.spdf.df_11.core,rmj.var.spdf.df_11.core)
rmj.var.spdf.df_13.core <- CVaddFun(rmj.spdf.df_13.core,rmj.var.spdf.df_13.core)
rmj.var.spdf.df_14.core <- CVaddFun(rmj.spdf.df_14.core,rmj.var.spdf.df_14.core)
rmj.var.spdf.df_16.core <- CVaddFun(rmj.spdf.df_16.core,rmj.var.spdf.df_16.core)
rmj.var.spdf.df_18.core <- CVaddFun(rmj.spdf.df_18.core,rmj.var.spdf.df_18.core)
rmj.var.spdf.df_20.core <- CVaddFun(rmj.spdf.df_20.core,rmj.var.spdf.df_20.core)


# Save the SPDF's with the CV value but no bins (for continuous plotting of the CV)
write.csv(rmj.var.spdf.df_10.core,file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_10.core.csv")
write.csv(rmj.var.spdf.df_11.core,file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_11.core.csv")
write.csv(rmj.var.spdf.df_13.core,file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_13.core.csv")
write.csv(rmj.var.spdf.df_14.core,file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_14.core.csv")
write.csv(rmj.var.spdf.df_16.core,file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_16.core.csv")
write.csv(rmj.var.spdf.df_18.core,file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_18.core.csv")
write.csv(rmj.var.spdf.df_20.core,file="Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_20.core.csv")


### add bins

## Quartiles

# put spdf's into a list
dfs <- list(rmj.var.spdf.df_10.core,rmj.var.spdf.df_11.core,rmj.var.spdf.df_13.core,rmj.var.spdf.df_14.core,
            rmj.var.spdf.df_16.core,rmj.var.spdf.df_18.core,rmj.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("rmj.var.spdf.df_10.core","rmj.var.spdf.df_11.core","rmj.var.spdf.df_13.core",
                "rmj.var.spdf.df_14.core","rmj.var.spdf.df_16.core","rmj.var.spdf.df_18.core",
                "rmj.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(rmj.var.spdf.df_10.core,file="Results/RMJ/Plots/variance/core_only/spdf/bins/rmj.var.spdf.df_10.core.csv")
write.csv(rmj.var.spdf.df_11.core,file="Results/RMJ/Plots/variance/core_only/spdf/bins/rmj.var.spdf.df_11.core.csv")
write.csv(rmj.var.spdf.df_13.core,file="Results/RMJ/Plots/variance/core_only/spdf/bins/rmj.var.spdf.df_13.core.csv")
write.csv(rmj.var.spdf.df_14.core,file="Results/RMJ/Plots/variance/core_only/spdf/bins/rmj.var.spdf.df_14.core.csv")
write.csv(rmj.var.spdf.df_16.core,file="Results/RMJ/Plots/variance/core_only/spdf/bins/rmj.var.spdf.df_16.core.csv")
write.csv(rmj.var.spdf.df_18.core,file="Results/RMJ/Plots/variance/core_only/spdf/bins/rmj.var.spdf.df_18.core.csv")
write.csv(rmj.var.spdf.df_20.core,file="Results/RMJ/Plots/variance/core_only/spdf/bins/rmj.var.spdf.df_20.core.csv")




## custom bins

# change column name from group2 to CV
rmj.var.spdf.df_10.core <- rmj.var.spdf.df_10.core %>% dplyr::rename(CV=group2)
rmj.var.spdf.df_11.core <- rmj.var.spdf.df_11.core %>% dplyr::rename(CV=group2)
rmj.var.spdf.df_13.core <- rmj.var.spdf.df_13.core %>% dplyr::rename(CV=group2)
rmj.var.spdf.df_14.core <- rmj.var.spdf.df_14.core %>% dplyr::rename(CV=group2)
rmj.var.spdf.df_16.core <- rmj.var.spdf.df_16.core %>% dplyr::rename(CV=group2)
rmj.var.spdf.df_18.core <- rmj.var.spdf.df_18.core %>% dplyr::rename(CV=group2)
rmj.var.spdf.df_20.core <- rmj.var.spdf.df_20.core %>% dplyr::rename(CV=group2)


# put spdf's into a list
dfs <- list(rmj.var.spdf.df_10.core,rmj.var.spdf.df_11.core,rmj.var.spdf.df_13.core,rmj.var.spdf.df_14.core,
            rmj.var.spdf.df_16.core,rmj.var.spdf.df_18.core,rmj.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("rmj.var.spdf.df_10.core","rmj.var.spdf.df_11.core","rmj.var.spdf.df_13.core",
                "rmj.var.spdf.df_14.core","rmj.var.spdf.df_16.core","rmj.var.spdf.df_18.core",
                "rmj.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar2)

# split elements into original dataframes
list2env(dfs, globalenv())


# save the SPDFs with the custom bins
write.csv(rmj.var.spdf.df_10.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_10.core.csv")
write.csv(rmj.var.spdf.df_11.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_11.core.csv")
write.csv(rmj.var.spdf.df_13.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_13.core.csv")
write.csv(rmj.var.spdf.df_14.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_14.core.csv")
write.csv(rmj.var.spdf.df_16.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_16.core.csv")
write.csv(rmj.var.spdf.df_18.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_18.core.csv")
write.csv(rmj.var.spdf.df_20.core,
          file="Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_20.core.csv")


      # Discrete bins ####

# load spdfs (which has already had CV calculated and then put into bins)
rmj.var.spdf.df_10.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_10.core.csv")
rmj.var.spdf.df_11.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_11.core.csv")
rmj.var.spdf.df_13.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_13.core.csv")
rmj.var.spdf.df_14.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_14.core.csv")
rmj.var.spdf.df_16.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_16.core.csv")
rmj.var.spdf.df_18.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_18.core.csv")
rmj.var.spdf.df_20.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/bins/custom/rmj.var.spdf.df_20.core.csv")

# make CV a factor and re-order. I have only done this for 2020 but you can copy the code for the other years if you need to
rmj.var.spdf.df_20.core$CV <- as.factor(rmj.var.spdf.df_20.core$CV)
rmj.var.spdf.df_20.core$CV <- factor(rmj.var.spdf.df_20.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))

# plot CV in bins
RMJ_10_plot_bin_GS_var <- GSplotBin(rmj.var.spdf.df_10.core,survey.area.core,"Variance","CV")
RMJ_11_plot_bin_GS_var <- GSplotBin(rmj.var.spdf.df_11.core,survey.area.core,"Variance","CV")
RMJ_13_plot_bin_GS_var <- GSplotBin(rmj.var.spdf.df_13.core,survey.area.core,"Variance","CV")
RMJ_14_plot_bin_GS_var <- GSplotBin(rmj.var.spdf.df_14.core,survey.area.core,"Variance","CV")
RMJ_16_plot_bin_GS_var <- GSplotBin(rmj.var.spdf.df_16.core,survey.area.core,"Variance","CV")
RMJ_18_plot_bin_GS_var <- GSplotBin(rmj.var.spdf.df_18.core,survey.area.core,"Variance","CV")
RMJ_20_plot_bin_GS_var <- GSplotBin(rmj.var.spdf.df_20.core,"CV",survey.area.core,"Variance","CV")

# save 
saveplot(RMJ_10_plot_bin_GS_var,"Results/RMJ/Plots/variance/core_only/greyscale/bins/RMJ_10_plot_bin_GS_var.png")
saveplot(RMJ_11_plot_bin_GS_var,"Results/RMJ/Plots/variance/core_only/greyscale/bins/RMJ_11_plot_bin_GS_var.png")
saveplot(RMJ_13_plot_bin_GS_var,"Results/RMJ/Plots/variance/core_only/greyscale/bins/RMJ_13_plot_bin_GS_var.png")
saveplot(RMJ_14_plot_bin_GS_var,"Results/RMJ/Plots/variance/core_only/greyscale/bins/RMJ_14_plot_bin_GS_var.png")
saveplot(RMJ_16_plot_bin_GS_var,"Results/RMJ/Plots/variance/core_only/greyscale/bins/RMJ_16_plot_bin_GS_var.png")
saveplot(RMJ_18_plot_bin_GS_var,"Results/RMJ/Plots/variance/core_only/greyscale/bins/RMJ_18_plot_bin_GS_var.png")
saveplot(RMJ_20_plot_bin_GS_var,"Results/RMJ/Plots/variance/core_only/greyscale/bins/RMJ_20_plot_bin_GS_var.png")



      # Continuous ####


## Load spatial dataframes
rmj.var.spdf.df_10.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_10.core.csv")
#rmj.var.spdf.df_11.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_11.core.csv")
#rmj.var.spdf.df_13.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_13.core.csv")
#rmj.var.spdf.df_14.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_14.core.csv")
#rmj.var.spdf.df_16.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_16.core.csv")
#rmj.var.spdf.df_18.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_18.core.csv")
rmj.var.spdf.df_20.core <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_20.core.csv")


# greyscale plots
RMJ_varplot_final10.core.bw <- GSplotFun(rmj.var.spdf.df_10.core, survey.area.core, "variance", "2010")
RMJ_varplot_final11.core.bw <- GSplotFun(rmj.var.spdf.df_11.core, survey.area.core, "variance")
RMJ_varplot_final13.core.bw <- GSplotFun(rmj.var.spdf.df_13.core, survey.area.core, "variance")
RMJ_varplot_final14.core.bw <- GSplotFun(rmj.var.spdf.df_14.core, survey.area.core, "variance")
RMJ_varplot_final16.core.bw <- GSplotFun(rmj.var.spdf.df_16.core, survey.area.core, "variance")
RMJ_varplot_final18.core.bw <- GSplotFun(rmj.var.spdf.df_18.core, survey.area.core, "variance")
RMJ_varplot_final20.core.bw <- GSplotFun(rmj.var.spdf.df_20.core, survey.area.core, "variance", "2020")

# save greyscale
saveplot(RMJ_varplot_final10.core.bw, "Results/RMJ/Plots/variance/core_only/greyscale/2010_RMJ_var.core.bw.png")
saveplot(RMJ_varplot_final11.core.bw, "Results/RMJ/Plots/variance/core_only/greyscale/2011_RMJ_var.core.bw.png")
saveplot(RMJ_varplot_final13.core.bw, "Results/RMJ/Plots/variance/core_only/greyscale/2013_RMJ_var.core.bw.png")
saveplot(RMJ_varplot_final14.core.bw, "Results/RMJ/Plots/variance/core_only/greyscale/2014_RMJ_var.core.bw.png")
saveplot(RMJ_varplot_final16.core.bw, "Results/RMJ/Plots/variance/core_only/greyscale/2016_RMJ_var.core.bw.png")
saveplot(RMJ_varplot_final18.core.bw, "Results/RMJ/Plots/variance/core_only/greyscale/2018_RMJ_var.core.bw.png")
saveplot(RMJ_varplot_final20.core.bw, "Results/RMJ/Plots/variance/core_only/greyscale/2020_RMJ_var.core.bw.png")

# colour plots
RMJ_varplot_final10.core.col <- CLplotFun(rmj.var.spdf.df_10.core, survey.area.core, "variance")
RMJ_varplot_final11.core.col <- CLplotFun(rmj.var.spdf.df_11.core, survey.area.core, "variance")
RMJ_varplot_final13.core.col <- CLplotFun(rmj.var.spdf.df_13.core, survey.area.core, "variance")
RMJ_varplot_final14.core.col <- CLplotFun(rmj.var.spdf.df_14.core, survey.area.core, "variance")
RMJ_varplot_final16.core.col <- CLplotFun(rmj.var.spdf.df_16.core, survey.area.core, "variance")
RMJ_varplot_final18.core.col <- CLplotFun(rmj.var.spdf.df_18.core, survey.area.core, "variance")
RMJ_varplot_final20.core.col <- CLplotFun(rmj.var.spdf.df_20.core, survey.area.core, "variance")

# save colour
saveplot(RMJ_varplot_final10.core.col, "Results/RMJ/Plots/variance/core_only/colour/2010_RMJ_var.core.col.png")
saveplot(RMJ_varplot_final11.core.col, "Results/RMJ/Plots/variance/core_only/colour/2011_RMJ_var.core.col.png")
saveplot(RMJ_varplot_final13.core.col, "Results/RMJ/Plots/variance/core_only/colour/2013_RMJ_var.core.col.png")
saveplot(RMJ_varplot_final14.core.col, "Results/RMJ/Plots/variance/core_only/colour/2014_RMJ_var.core.col.png")
saveplot(RMJ_varplot_final16.core.col, "Results/RMJ/Plots/variance/core_only/colour/2016_RMJ_var.core.col.png")
saveplot(RMJ_varplot_final18.core.col, "Results/RMJ/Plots/variance/core_only/colour/2018_RMJ_var.core.col.png")
saveplot(RMJ_varplot_final20.core.col, "Results/RMJ/Plots/variance/core_only/colour/2020_RMJ_var.core.col.png")


#### Banteng ##################################################################
## Load data ####

### There were no banteng observations in 2020, and therefore I have not re-done the modelling for this species (as without any new observations, the models will be exactly the same).  I will just predict abundance over the new 2020 preddata

# Observation data. Unique to species 
btg_obsdata <- read.csv("Species_Data/BTG/R Data/obsdata.csv", header = TRUE)
btg_obsdata$object <- as.factor(btg_obsdata$object)
btg_obsdata$Sample.Label <- as.factor(btg_obsdata$Sample.Label)
#btg_obsdata$year <- as.factor(btg_obsdata$year)
str(btg_obsdata)
head(btg_obsdata)

# Transect data. Unique to species
btg_distdata <- read.csv("Species_Data/BTG/R Data/distdata.csv", header = TRUE)
btg_distdata$object <- as.factor(btg_distdata$object)
btg_distdata$NameObserver <- as.factor(btg_distdata$NameObserver)
btg_distdata$transect <- as.factor(btg_distdata$transect)
btg_distdata$year <- as.factor(btg_distdata$year)
btg_distdata$date <- as.Date(btg_distdata$date, format = "%d/%m/%Y")
str(btg_distdata)
head(rmj_distdata)

## Plot the covariates across the grid, with group sizes ####

# Warning - plots take a few minutes to run

# habitat
plot_BTGobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata200$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=btg_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/BTG/plot_BTGobs_habitat.png", plot_BTGobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_BTGobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys_200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=btg_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/BTG/plot_BTGobs_dstWater.png", plot_BTGobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_BTGobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=btg_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/BTG/plot_BTGobs_dstStlmnt.png", plot_BTGobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_BTGobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=btg_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/BTG/plot_BTGobs_dstRoad.png", plot_BTGobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_BTGobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=btg_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/BTG/plot_BTGobs_dstBorder.png", plot_BTGobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_BTGobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=btg_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/BTG/plot_BTGobs_dstStation.png", plot_BTGobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_BTGobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=btg_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/BTG/plot_BTGobs_dstELC.png", plot_BTGobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_BTGobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys_200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=btg_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/BTG/plot_RMJobs_elev.png", plot_BTGobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)


## Exploratory plots & linear models ####

## Checking that there are no large gaps in the range of variables for banteng observations.  Bit of a silly exercise because there are only 20 observations so of course there are going to be large gaps...

# subset segdata to get only the segments with BTG observations
btg_varcheck <- segdata[match(btg_obsdata$Sample.Label,segdata$Sample.Label), ]

# habitat - no observations in dense forest
plot(segdata$habitat)
plot(btg_varcheck$habitat)

# dstWater - actually not too bad. Gap between 1000 and 1200, and aobut 1700 and 2500
hist(segdata$dstWater)
hist(btg_varcheck$dstWater)

# dstStlmnt - again surprisingly not horrible. Gap 6000 and 7000, and then no observations beyond 8000 which isnt great as the variable goes up to 12000.  This is the same as some of the above species
hist(segdata$dstStlmnt)
hist(btg_varcheck$dstStlmnt)

# dstRoad - no observations beyond 1600, whereas the variable goes up to 4000
hist(segdata$dstRoad)
hist(btg_varcheck$dstRoad)

# dstBorder - gap between 0-15km
hist(segdata$dstBorder)
hist(btg_varcheck$dstBorder)

# dstStation - gap between 0-5km
hist(segdata$dstStation)
hist(btg_varcheck$dstStation)

# elevation - no obs beyond 350m whereas variable goes up to 750m
hist(segdata$elevation)
hist(btg_varcheck$elevation)


## Histograms

# Distance 
btg_h1 <- ggplot(btg_distdata, aes(distance))+ geom_histogram(binwidth = 1)
btg_h2 <- ggplot(btg_distdata, aes(distance))+ geom_histogram(binwidth = 5)
btg_h3 <- ggplot(btg_distdata, aes(distance))+ geom_histogram(binwidth = 10)
btg_h4 <- ggplot(btg_distdata, aes(distance))+ geom_histogram(binwidth = 15)
btg_h5 <- ggplot(btg_distdata, aes(distance))+ geom_histogram(binwidth = 20)
btg_h6 <- ggplot(btg_distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(btg_h1,btg_h2,btg_h3,btg_h4,btg_h5,btg_h6)
# Histogram not looking good...too few observations really 


# cluster size, observer, habitat, year, month, transect
btg_h7 <- ggplot(btg_distdata, aes(size))+geom_histogram(binwidth = 0.5)
btg_h8 <- ggplot(btg_distdata, aes(NameObserver))+geom_histogram(stat="count")
btg_h9 <- ggplot(btg_distdata, aes(habitat))+geom_histogram(stat="count")
btg_h10 <- ggplot(btg_distdata, aes(year))+geom_histogram(stat="count")
btg_h11 <- ggplot(btg_distdata, aes(month))+geom_histogram(stat="count")
btg_h12 <- ggplot(btg_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(btg_h7,btg_h8,btg_h9,btg_h10,btg_h11,btg_h12)
# Most observations in open forest. Most observations are of single animals but there are some of larger groups. This species is known to exist in herds and so abundance of groups should be estimated.

## Plots of distance against variables

plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

btg_d1 <- ggplot(btg_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
btg_d2 <- ggplot(btg_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
btg_d3 <- ggplot(btg_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
btg_d4 <- ggplot(btg_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
btg_d5 <- ggplot(btg_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(btg_d1,btg_d2,btg_d3,btg_d4,btg_d5)
# No difference in distances between habitats. Doesn't look like there is a relationship between distance and group size, but really there are too few observations


## Plots of cluster size against variables
btg_s1 <- ggplot(btg_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
btg_s2 <- ggplot(btg_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
btg_s3 <- ggplot(btg_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
btg_s4 <- ggplot(btg_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
btg_s5 <- ggplot(btg_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(btg_s1,btg_s2,btg_s3,btg_s4,btg_s5)
# Groups size by year is an interesting plot.  Shows what looks like a decrease in group size from 2010 to 2018.  Again, too few obs to say there is a relationship but it would tie in with the assumtion that bantenge have been heavily persecuted over that time period. 

## Linear models

# group size ~ distance
newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=btg_distdata)
plot(btg_distdata$size~btg_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# It looks like there is the beginning of a positive relationship between distance and group size, but the model doesn't provide much support (distance coefficient: 0.02, p=0.5, Rsq=0.02).  Perhaps with more data that relationship would be stronger. 

## Estimating the detection function ####

# the below DF is the same as for the CDS. The original model selection from DSM (in obsolete code script) in fact came to the same conclusion

# Model with no covariates, Uniform Key and simple polynomial adjustment 
btgDF.un2 <- ds(btg_distdata, truncation = 65, key = "unif", adjustment = "poly")

  
## Fitting a spatial model ####


# The best detection function model for BTG is btgDF.un2

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated. This is because despite most of the observations being of single individuals, the species is gregarious (or generally should be anyway). 

# I am not including dstELC in the models for BTG, as I don't think it is an appropriate variable, especially now that we are going to be predicting into the buffer zone

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (65m)
btg_obsdata <- btg_obsdata %>% filter(distance <= 65)

  ## Quasipoisson response ####

# Saturated model with all covariates. 
btgDSM.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(dstRoad,bs="ts")+ 
                         s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                         habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.sat)
par(mfrow=c(2,3))
plot(btgDSM.sat, scale = 0)
gam.check(btgDSM.sat)
# DE = 65.3. All terms sig. gam.check suggests k is fine. Plots aren't great. dstBorder and elevation look vaguely sensible. And perhaps dstStation although just not enough data.  I think with this few observations we need to have as few variables as possible

# remove dstWater
btgDSM.qp.2 <- dsm(Nhat ~ s(dstStlmnt,bs="ts")+ s(dstRoad,bs="ts")+ 
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                          habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.qp.2)
par(mfrow=c(2,3))
plot(btgDSM.qp.2, scale = 0)
gam.check(btgDSM.qp.2)
# DE = 61.9.  All terms sig. elevation could use an increase in k now

# compare models
anova(btgDSM.sat,btgDSM.qp.2,test="Chisq")
# saturated model is better but I think its got way too many vars for so few obs

# remove dstRoad and increase k for elevation
btgDSM.qp.3 <- dsm(Nhat ~ s(dstStlmnt,bs="ts")+ s(dstStation,bs="ts")+
                          s(dstBorder,bs="ts")+ s(elevation,bs="ts",k=15)+ 
                          habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.qp.3)
dev.off()
par(mfrow=c(2,2))
plot(btgDSM.qp.3, scale = 0)
gam.check(btgDSM.qp.3)
# DE = 57.4. All terms sig. gam.check now suggests k needs to be increased for all.

# Compare models
anova(btgDSM.qp.2,btgDSM.qp.3,test="Chisq")
# qp.3 is worse than qp.2.

# Try a very basic model with just dstBorder, elevation, and habitat. These are the vars that I think in reality are probably quite good predictors. I will increase k for all
btgDSM.qp.4 <- dsm(Nhat ~ s(dstBorder,bs="ts",k=15)+ s(elevation,bs="ts",k=15)+ 
                          habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.qp.4)
dev.off()
par(mfrow=c(2,2))
plot(btgDSM.qp.4, scale = 0)
gam.check(btgDSM.qp.4)
# DE = 30.1. None of the terms are sig.

# as qp.3 but with k increased for all
btgDSM.qp.5 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=15)+ s(dstStation,bs="ts",k=15)+
                          s(dstBorder,bs="ts",k=15)+ s(elevation,bs="ts",k=15)+ 
                          habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.qp.5)
dev.off()
par(mfrow=c(2,2))
plot(btgDSM.qp.5, scale = 0)
gam.check(btgDSM.qp.5)
# DE = 61.3. All terms sig. 

# Compare models
anova(btgDSM.qp.5,btgDSM.qp.3,test="Chisq")
# qp.5 is better

# Could try predictions with qp.5

  ## Tweedie response ####

btgDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                             s(dstRoad,bs="ts")+ s(dstBorder,bs="ts")+
                             s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                             habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.tw.sat)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.tw.sat)
dev.off()
par(mfrow=c(3,2))
plot(btgDSM.tw.sat, scale=0)
# AIC = 12408.17, DE = 19. Only dstWater, dstStlmnt, and maybe dstRoad sig. dstWater and dstStlmnt could be parametric. gam.check doesn't suggest k should be increased.

# Remove dstborder, dstStation, and elevation
btgDSM.tw.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                          s(dstRoad,bs="ts")+ 
                          habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.tw.2)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.tw.2)
dev.off()
par(mfrow=c(3,2))
plot(btgDSM.tw.2, scale=0)
# AIC = 12396.35, DE = 41.4.  dstRoad no longer sig. gam.check suggests k is fine. dstWater and dstStlmnt could be parametric still

# Remove dstRoad
btgDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                          habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.tw.3)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.tw.3)
dev.off()
par(mfrow=c(3,2))
plot(btgDSM.tw.3, scale=0)
# AIC = 12405.41, DE= 12.8 Not sure why there is such a big difference in AIC and DE between tw.2 and tw.3.

# As tw.2 but with dstWater as para term
btgDSM.tw.4 <- dsm(Nhat ~ s(dstStlmnt,bs="ts")+ s(dstRoad,bs="ts")+
                          habitat + dstWater,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.tw.4)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.tw.4)
dev.off()
par(mfrow=c(1,2))
plot(btgDSM.tw.4, scale=0)
# AIC = 12396.11, DE = 41.3. dstWater sig as a para term. dsstRoad no longer sig. Plots look sensible

# Remove dstRoad
btgDSM.tw.5 <- dsm(Nhat ~ s(dstStlmnt,bs="ts")+ habitat + dstWater,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.tw.5)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.tw.5)
dev.off()
par(mfrow=c(1,2))
plot(btgDSM.tw.5, scale=0)
# AIC = 12406, DE = 8.98. Strange. dstRoad is not sig but seems to have a huge impact on the model. I will try it as a para term

# as tw.4 but with dstRoad as a para term
btgDSM.tw.6 <- dsm(Nhat ~ s(dstStlmnt,bs="ts")+
                          habitat + dstWater + dstRoad,
                  btgDF.un2, segdata, btg_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.tw.6)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.tw.6)
dev.off()
par(mfrow=c(1,2))
plot(btgDSM.tw.6, scale=0)
# AIC = 12408.3, DE = 8.25. Same thing has happened with dstRoad. It's not sig as a para term.  Perhaps I just need to keep it in! 

# Try predictions with btgDSM.tw.4

  ## Negative binomial response ####

btgDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                             s(dstRoad,bs="ts")+ s(dstBorder,bs="ts")+
                             s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                             habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML", group=TRUE,
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.nb.sat)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.nb.sat)
dev.off()
par(mfrow=c(3,2))
plot(btgDSM.nb.sat, scale=0)
# AIC = 269.02, DE = 26.5. Only dstWater and dstStlmnt sig. Both could be para terms. gam.check not suggesting an increase in k

# Remove all except dstWater and dstStlmnt
btgDSM.nb.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                           habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML", group=TRUE,
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.nb.2)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.nb.2)
dev.off()
par(mfrow=c(2,1))
plot(btgDSM.nb.2, scale=0)
# AIC = 269.03, DE = 26.5. both terms still sig, but defo should check as para terms. Plots look sensible

# Check incresing k
btgDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+ 
                           habitat,
                  btgDF.un2, segdata, btg_obsdata, method = "REML", group=TRUE,
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.nb.3)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.nb.3)
dev.off()
par(mfrow=c(1,2))
plot(btgDSM.nb.3, scale=0)
# AIC = 269.03, DE = 26.5. increasing k has done nowt. Probably because DEF is close to 0.  

# dstWater as a paremtric term
btgDSM.nb.4 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=15)+ 
                           habitat + dstWater,
                  btgDF.un2, segdata, btg_obsdata, method = "REML", group=TRUE,
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.nb.4)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.nb.4)
dev.off()
par(mfrow=c(1,1))
plot(btgDSM.nb.4, scale=0)
# AIC = 269.12, DE = 26.6.  dstWater sig as para term.

# Try increase k for dstStlmnt
btgDSM.nb.5 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=25)+ 
                           habitat + dstWater,
                  btgDF.un2, segdata, btg_obsdata, method = "REML", group=TRUE,
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.nb.5)
dev.off()
par(mfrow=c(2,2))
gam.check(btgDSM.nb.5)
dev.off()
par(mfrow=c(1,1))
plot(btgDSM.nb.5, scale=0)
# AIC = 269.12, DE = 26.6. 


## best GAM model is probably btgDSM.nb.4 

  ## GLM ####

# The only significant terms had EDF of <1 (but close to 1), and so using GAMs and smoothing the terms doesn't seem appropriate.

# glm with no smooth terms. 
# poisson because we're looking at abundance (not density) and so the response is count data
btgDSM.glm.1 <- dsm(Nhat ~ habitat + dstStlmnt + dstWater, 
                  btgDF.un2, segdata, btg_obsdata, method = "REML", group=TRUE,
                  family = poisson(link = "log"), engine = "glm",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.glm.1)
par(mfrow=c(2,2))
plot(btgDSM.glm.1)

# See what happens when I remove dstStlmnt
btgDSM.glm.2 <- dsm(Nhat ~ habitat + dstWater, 
                  btgDF.un2, segdata, btg_obsdata, method = "REML", group=TRUE,
                  family = poisson(link = "log"), engine = "glm",
                  segment.area = segdata$Sample.Area)
summary(btgDSM.glm.2)
par(mfrow=c(2,2))
plot(btgDSM.glm.2)

# compare models
anova(btgDSM.glm1,btgDSM.glm.2, test="Chisq")

# model with both terms is better

## Plot models and effects

# create new data for dstWater
newdat.btg.dstWater <- crossing(habitat=c("D","DEF","NF","O","W"), 
                        dstWater=seq(0,3874,length=100),
                        dstStlmnt=rep(4729,100))

# create new data for dstStlmnt
newdat.btg.dstStlmnt <- crossing(habitat=c("D","DEF","NF","O","W"), 
                        dstStlmnt=seq(0,12000,length=100),
                        dstWater=rep(636,length=100))

# predict over new data
pred.btg.dstWater <- predict(btgDSM.glm.1, newdata = newdat.btg.dstWater, 
                             type = "response", se = T, off.set = 40000)

pred.btg.dstStlmnt <- predict(btgDSM.glm.1, newdata = newdat.btg.dstStlmnt, 
                              type = "response", se = T, off.set = 40000)

# bind new data and prediction together
btg.mod.dstWater <- bind_cols(newdat.btg.dstWater, as.tibble(pred.btg.dstWater))
btg.mod.dstWater$abundance <- btg.mod.dstWater$fit
btg.mod.dstStlmnt <- bind_cols(newdat.btg.dstStlmnt, as.tibble(pred.btg.dstStlmnt))
btg.mod.dstStlmnt$abundance <- btg.mod.dstStlmnt$fit

# scatter plots
plot(btg.mod.dstWater$dstWater, pred.btg.dstWater$fit)
plot(btg.mod.dstStlmnt$dstStlmnt, pred.btg.dstStlmnt$fit)

# calculate CIs
btg.mod.dstWater$upr <- btg.mod.dstWater$fit+(1.96*btg.mod.dstWater$se.fit) 
btg.mod.dstWater$lwr <- btg.mod.dstWater$fit-(1.96*btg.mod.dstWater$se.fit)
btg.mod.dstStlmnt$upr <- btg.mod.dstStlmnt$fit+(1.96*btg.mod.dstStlmnt$se.fit)
btg.mod.dstStlmnt$lwr <- btg.mod.dstStlmnt$fit-(1.96*btg.mod.dstStlmnt$se.fit)

# plot dstWater
ggplot(data=btg.mod.dstWater,aes(x=dstWater, y=abundance, group=habitat, colour=habitat))+
  geom_ribbon(data=btg.mod.dstWater,aes(ymin=lwr, ymax=upr),fill=NA, linetype = 2, size=0.75)+
  scale_color_viridis(discrete=TRUE)+
  geom_line(data=btg.mod.dstWater,aes(x=dstWater, y=abundance),size=1)

# dstStlmnt
ggplot(data=btg.mod.dstStlmnt,aes(x=dstStlmnt, y=abundance, group=habitat, colour=habitat))+
  geom_ribbon(data=btg.mod.dstStlmnt,aes(ymin=lwr, ymax=upr),fill=NA, linetype = 3, size=0.75)+
  scale_color_viridis(discrete=TRUE)+
  geom_line(data=btg.mod.dstStlmnt,aes(x=dstStlmnt, y=abundance),size=1)


## Model selection ####

# below was the model selection for the GAMs, but I have now decided that the GLM btgDSM.glm.1 is the best

# best QP model is btgDSM.qp.5
# best TW model is btgDSM.tw.4
# best NB model is btgDSM.nb.4

summary(btgDSM.qp.5) 
summary(btgDSM.tw.4)
summary(btgDSM.nb.4)

# All have decent DE.  NB model has better AIC than the TW model

# compare QP with NB
anova(btgDSM.qp.5,btgDSM.nb.4,test="Chisq")
# NB has much less residual deviance

# check Q-Q plots
par(mfrow=c(2,2))
gam.check(btgDSM.qp.5)
gam.check(btgDSM.nb.4)


## Autocorrelation ####

# autocorrelation for GLM
dsm.cor(btgDSM.glm.1, max.lag=15,Segment.Label="Sample.Label", resid.type = "pearson")

# No autocorrelation

## Abundance estimation ####


# predict over 2010 habitat using btgDSM.glm.1 (glm)
btg.glm.pred10.core <- predict(btgDSM.glm.1, preddata10_core, off.set=40000)
#write.csv(btg.glm.pred10.core,file="Results/BTG/core_only/btg.pred10.core.csv")

btg.glm.pred11.core <- predict(btgDSM.glm.1, preddata11_core, off.set=40000)
#write.csv(btg.glm.pred11.core,file="Results/BTG/core_only/btg.pred11.core.csv")

btg.glm.pred13.core <- predict(btgDSM.glm.1, preddata13_core, off.set=40000)
#write.csv(btg.glm.pred13.core,file="Results/BTG/core_only/btg.pred13.core.csv")

btg.glm.pred14.core <- predict(btgDSM.glm.1, preddata14_core, off.set=40000)
#write.csv(btg.glm.pred14.core,file="Results/BTG/core_only/btg.pred14.core.csv")

btg.glm.pred16.core <- predict(btgDSM.glm.1, preddata16_core, off.set=40000)
#write.csv(btg.glm.pred16.core,file="Results/BTG/core_only/btg.pred16.core.csv")

btg.glm.pred18.core <- predict(btgDSM.glm.1, preddata18_core, off.set=40000)
#write.csv(btg.glm.pred18.core,file="Results/BTG/core_only/btg.pred18.core.csv")

btg.glm.pred20.core <- predict(btgDSM.glm.1, preddata20_core, off.set=40000)
#write.csv(btg.glm.pred20.core,file="Results/BTG/core_only/btg.pred20.core.csv")



# create dataframe for plotting
btg.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = btg.glm.pred10.core)

btg.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = btg.glm.pred11.core)

btg.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = btg.glm.pred13.core)

btg.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = btg.glm.pred14.core)

btg.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = btg.glm.pred16.core)

btg.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = btg.glm.pred18.core)

btg.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = btg.glm.pred20.core)


## This creates a dataframe that can be plotted as a map
btg.spdf.df_10.core <- abunPlotDF(btg.df.Final10.core, pred.polys_200)
btg.spdf.df_11.core <- abunPlotDF(btg.df.Final11.core, pred.polys_200)
btg.spdf.df_13.core <- abunPlotDF(btg.df.Final13.core, pred.polys_200)
btg.spdf.df_14.core <- abunPlotDF(btg.df.Final14.core, pred.polys_200)
btg.spdf.df_16.core <- abunPlotDF(btg.df.Final16.core, pred.polys_200)
btg.spdf.df_18.core <- abunPlotDF(btg.df.Final18.core, pred.polys_200)
btg.spdf.df_20.core <- abunPlotDF(btg.df.Final20.core, pred.polys_200)

# save SPDFs
#write.csv(btg.spdf.df_10.core,file="Results/BTG/Plots/core_only/spdf/btg.spdf.df_10.core.csv")
#write.csv(btg.spdf.df_11.core,file="Results/BTG/Plots/core_only/spdf/btg.spdf.df_11.core.csv")
#write.csv(btg.spdf.df_13.core,file="Results/BTG/Plots/core_only/spdf/btg.spdf.df_13.core.csv")
#write.csv(btg.spdf.df_14.core,file="Results/BTG/Plots/core_only/spdf/btg.spdf.df_14.core.csv")
#write.csv(btg.spdf.df_16.core,file="Results/BTG/Plots/core_only/spdf/btg.spdf.df_16.core.csv")
#write.csv(btg.spdf.df_18.core,file="Results/BTG/Plots/core_only/spdf/btg.spdf.df_18.core.csv")
#write.csv(btg.spdf.df_20.core,file="Results/BTG/Plots/core_only/spdf/btg.spdf.df_20.core.csv")

    ## Plotting ####

# load spatial dataframes
btg.spdf.df_10.core <- read.csv("Results/BTG/Plots/core_only/spdf/btg.spdf.df_10.core.csv")
btg.spdf.df_11.core <- read.csv("Results/BTG/Plots/core_only/spdf/btg.spdf.df_11.core.csv")
btg.spdf.df_13.core <- read.csv("Results/BTG/Plots/core_only/spdf/btg.spdf.df_13.core.csv")
btg.spdf.df_14.core <- read.csv("Results/BTG/Plots/core_only/spdf/btg.spdf.df_14.core.csv")
btg.spdf.df_16.core <- read.csv("Results/BTG/Plots/core_only/spdf/btg.spdf.df_16.core.csv")
btg.spdf.df_18.core <- read.csv("Results/BTG/Plots/core_only/spdf/btg.spdf.df_18.core.csv")

# greyscale plots
BTG_plot_10_core_gr <- GSplotFun(btg.spdf.df_10.core, survey.area.core, "abundance")
BTG_plot_11_core_gr <- GSplotFun(btg.spdf.df_11.core, survey.area.core, "abundance")
BTG_plot_13_core_gr <- GSplotFun(btg.spdf.df_13.core, survey.area.core, "abundance")
BTG_plot_14_core_gr <- GSplotFun(btg.spdf.df_14.core, survey.area.core, "abundance")
BTG_plot_16_core_gr <- GSplotFun(btg.spdf.df_16.core, survey.area.core, "abundance")
BTG_plot_18_core_gr <- GSplotFun(btg.spdf.df_18.core, survey.area.core, "abundance")
BTG_plot_20_core_gr <- GSplotFun(btg.spdf.df_20.core, survey.area.core, "abundance")

# save greyscale
saveplot(BTG_plot_10_core_gr,"Results/BTG/Plots/core_only/greyscale/BTG_plot_10_core_gr.png")
saveplot(BTG_plot_11_core_gr,"Results/BTG/Plots/core_only/greyscale/BTG_plot_11_core_gr.png")
saveplot(BTG_plot_13_core_gr,"Results/BTG/Plots/core_only/greyscale/BTG_plot_13_core_gr.png")
saveplot(BTG_plot_14_core_gr,"Results/BTG/Plots/core_only/greyscale/BTG_plot_14_core_gr.png")
saveplot(BTG_plot_16_core_gr,"Results/BTG/Plots/core_only/greyscale/BTG_plot_16_core_gr.png")
saveplot(BTG_plot_18_core_gr,"Results/BTG/Plots/core_only/greyscale/BTG_plot_18_core_gr.png")
saveplot(BTG_plot_20_core_gr,"Results/BTG/Plots/core_only/greyscale/BTG_plot_20_core_gr.png")

# colour plots
BTG_plot_10_core_col <- CLplotFun(btg.spdf.df_10.core, survey.area.core, "abundance")
BTG_plot_11_core_col <- CLplotFun(btg.spdf.df_11.core, survey.area.core, "abundance")
BTG_plot_13_core_col <- CLplotFun(btg.spdf.df_13.core, survey.area.core, "abundance")
BTG_plot_14_core_col <- CLplotFun(btg.spdf.df_14.core, survey.area.core, "abundance")
BTG_plot_16_core_col <- CLplotFun(btg.spdf.df_16.core, survey.area.core, "abundance")
BTG_plot_18_core_col <- CLplotFun(btg.spdf.df_18.core, survey.area.core, "abundance")
BTG_plot_20_core_col <- CLplotFun(btg.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(BTG_plot_10_core_col,"Results/BTG/Plots/core_only/colour/BTG_plot_10_core_col.png")
saveplot(BTG_plot_11_core_col,"Results/BTG/Plots/core_only/colour/BTG_plot_11_core_col.png")
saveplot(BTG_plot_13_core_col,"Results/BTG/Plots/core_only/colour/BTG_plot_13_core_col.png")
saveplot(BTG_plot_14_core_col,"Results/BTG/Plots/core_only/colour/BTG_plot_14_core_col.png")
saveplot(BTG_plot_16_core_col,"Results/BTG/Plots/core_only/colour/BTG_plot_16_core_col.png")
saveplot(BTG_plot_18_core_col,"Results/BTG/Plots/core_only/colour/BTG_plot_18_core_col.png")
saveplot(BTG_plot_20_core_col,"Results/BTG/Plots/core_only/colour/BTG_plot_20_core_col.png")


## Variance estimation ####
  ## Core only ####
    
## I am not sure how to get the variance estimates when a GLM instead of a GAM is used.  I am probably not going to incldue the BTG DSM anyway, so not a big issue at the moment. 



    ## Plotting variance ####

# NA

#### Yellow-cheeked crested gibbon ############################################
## Load data ####

# Observation data. Unique to species 
ycg_obsdata <- read.csv("Species_Data/YCG/R Data/obsdata.csv", header = TRUE)
ycg_obsdata$object <- as.factor(ycg_obsdata$object)
ycg_obsdata$Sample.Label <- as.factor(ycg_obsdata$Sample.Label)
str(ycg_obsdata)
head(ycg_obsdata)

# Transect data. Unique to species. YCG needs the scaled and centered continous variable stratum, as this is used in the DF model (taken from CDS analysis)
ycg_distdata <- read.csv("Species_Data/YCG/R Data/distdata.csv", header = TRUE)
ycg_distdata$object <- as.factor(ycg_distdata$object)
ycg_distdata$NameObserver <- as.factor(ycg_distdata$NameObserver)
ycg_distdata$transect <- as.factor(ycg_distdata$transect)
ycg_distdata$stratum <- as.vector(scale(ycg_distdata$year, center=T, scale=T))
ycg_distdata$year <- as.factor(ycg_distdata$year)
ycg_distdata$date <- as.Date(ycg_distdata$date, format = "%d/%m/%Y")
str(ycg_distdata)
head(ycg_distdata)

# check for any observations on T20
ycg_obsdata[ycg_obsdata$Sample.Label=="20",]
ycg_distdata[ycg_distdata$transect=="20",]

## Plot the covariates across the grid, with group sizes ####

# Warning - plots take a few minutes to run

# habitat
plot_YCGobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata200$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=ycg_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/YCG/plot_YCGobs_habitat.png", plot_YCGobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_YCGobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys_200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=ycg_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/YCG/plot_YCGobs_dstWater.png", plot_YCGobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_YCGobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=ycg_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/YCG/plot_YCGobs_dstStlmnt.png", plot_YCGobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_YCGobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=ycg_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/YCG/plot_YCGobs_dstRoad.png", plot_YCGobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_YCGobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=ycg_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/YCG/plot_YCGobs_dstBorder.png", plot_YCGobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_YCGobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=ycg_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/YCG/plot_YCGobs_dstStation.png", plot_YCGobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_YCGobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=ycg_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/YCG/plot_YCGobs_dstELC.png", plot_YCGobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_YCGobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys_200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=ycg_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/YCG/plot_YCGobs_elev.png", plot_YCGobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)

## Exploratory plots & linear models ####

## Checking that there are no large gaps in the range of variables for gibbon observations.  

# subset segdata to get only the segments with YCG observations
ycg_varcheck <- segdata[match(ycg_obsdata$Sample.Label,segdata$Sample.Label), ]

# habitat - fine
plot(segdata$habitat)
plot(ycg_varcheck$habitat)

# dstWater - fine
hist(segdata$dstWater)
hist(ycg_varcheck$dstWater)

# dstStlmnt - gap between 8km and 10km
hist(segdata$dstStlmnt)
hist(ycg_varcheck$dstStlmnt)

# dstRoad - no observations between 1750 and 4000
hist(segdata$dstRoad)
hist(ycg_varcheck$dstRoad)

# dstBorder - no observations between 45000 and 55000. Don't think that's a massive problem
hist(segdata$dstBorder)
hist(ycg_varcheck$dstBorder)

# dstStation - fine
hist(segdata$dstStation)
hist(ycg_varcheck$dstStation)

# elevation - no observations between 600-750
hist(segdata$elevation)
hist(ycg_varcheck$elevation)


## Histograms

# Distance 
ycg_h1 <- ggplot(ycg_distdata, aes(distance))+ geom_histogram(binwidth = 1)
ycg_h2 <- ggplot(ycg_distdata, aes(distance))+ geom_histogram(binwidth = 5)
ycg_h3 <- ggplot(ycg_distdata, aes(distance))+ geom_histogram(binwidth = 10)
ycg_h4 <- ggplot(ycg_distdata, aes(distance))+ geom_histogram(binwidth = 15)
ycg_h5 <- ggplot(ycg_distdata, aes(distance))+ geom_histogram(binwidth = 20)
ycg_h6 <- ggplot(ycg_distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(ycg_h1,ycg_h2,ycg_h3,ycg_h4,ycg_h5,ycg_h6)
# Evidence of evasive movement - not surprising for this species. Quite a sharp drop off around 20-25m, so hazard rate may be a good DF model. Spike of observations around 45-50m.  Might cause some issues. 


# cluster size, observer, habitat, year, month, transect
ycg_h7 <- ggplot(ycg_distdata, aes(size))+geom_histogram(binwidth = 0.5)
ycg_h8 <- ggplot(ycg_distdata, aes(NameObserver))+geom_histogram(stat="count")
ycg_h9 <- ggplot(ycg_distdata, aes(habitat))+geom_histogram(stat="count")
ycg_h10 <- ggplot(ycg_distdata, aes(year))+geom_histogram(stat="count")
ycg_h11 <- ggplot(ycg_distdata, aes(month))+geom_histogram(stat="count")
ycg_h12 <- ggplot(ycg_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(ycg_h7,ycg_h8,ycg_h9,ycg_h10,ycg_h11,ycg_h12)
# cluster size actually looks quite realistic for this species as they live in small family groups. I will estimate abundance of groups, not individuals. Most observations in dense forest, but still a few in open habitat. These obs may be in riverine habitat, as this species is unlikely to actually live in open forest. 

## Plots of distance against variables

plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

ycg_d1 <- ggplot(ycg_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
ycg_d2 <- ggplot(ycg_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
ycg_d3 <- ggplot(ycg_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
ycg_d4 <- ggplot(ycg_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
ycg_d5 <- ggplot(ycg_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(ycg_d1,ycg_d2,ycg_d3,ycg_d4,ycg_d5)
# ditances are greater in NF and O habitat which makes sense. Interestingly distances are much greater in December compared to other months. December is when there is still a lot of foliage and vegetation, which is not what I would expect.  Group size ~ distance looks like the opposite relationship you would expect with size bias. 


## Plots of cluster size against variables
ycg_s1 <- ggplot(ycg_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
ycg_s2 <- ggplot(ycg_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
ycg_s3 <- ggplot(ycg_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
ycg_s4 <- ggplot(ycg_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
ycg_s5 <- ggplot(ycg_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(ycg_s1,ycg_s2,ycg_s3,ycg_s4,ycg_s5)
# No apparent relationship between group size and year

## Linear models

# group size ~ distance
newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=ycg_distdata)
plot(ycg_distdata$size~ycg_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# no strong support for size bias. There is a small negative relationship, but not worth worrying about 


## Estimating the detection function ####

# Below is the DF to be used (from CDS)

ycgDF.hn.cov5 <- ds(ycg_distdata, truncation = 50, formula = ~stratum)


   
## Fitting a spatial model ####


# The best detection function model for YCG is ycgDF.hr.cov4

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated.  
# I am not including dstELC in the models for YCG, as I don't think it is an appropriate variable, especially now that we are going to be predicting into the buffer zone

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
ycg_obsdata <- ycg_obsdata %>% filter(distance <= 50)


  ## Quasipoisson response ####

ycgDSM.qp.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+  
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ 
                            s(elevation,bs="ts")+ habitat,
                  ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.qp.sat)
par(mfrow=c(2,3))
plot(ycgDSM.qp.sat, scale = 0)
gam.check(ycgDSM.qp.sat)
# DE = 14.3, R2=0.01. All terms sig. All plots overfitted

# reduce k
ycgDSM.qp.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=5)+  
                       s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                       s(elevation,bs="ts",k=5)+ habitat,
                     ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                     family = quasipoisson(link = "log"), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.qp.sat2)
par(mfrow=c(2,3))
plot(ycgDSM.qp.sat2, scale = 0)
gam.check(ycgDSM.qp.sat2)
# DE = 9.3, R2=0.006. elevation no longer sig. dstWater, dstStlmnt still look overfitted. dstStation potentially linear

# reduce k for dstWater & dstStlmnt
ycgDSM.qp.sat3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=4)+ s(dstStlmnt,bs="ts",k=4)+  
                        s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                        s(elevation,bs="ts",k=5)+ habitat,
                      ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                      family = quasipoisson(link = "log"), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.qp.sat3)
par(mfrow=c(2,3))
plot(ycgDSM.qp.sat3, scale = 0)
gam.check(ycgDSM.qp.sat3)
# DE = 8.4, R2=0.006. dstStlmnt no longer sig. dstWater now looking linear

# remove elevation
ycgDSM.qp.4 <- dsm(Nhat ~ s(dstWater,bs="ts",k=4)+ s(dstStlmnt,bs="ts",k=4)+  
                        s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                         habitat,
                      ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                      family = quasipoisson(link = "log"), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.qp.4)
par(mfrow=c(2,2))
plot(ycgDSM.qp.4, scale = 0)
gam.check(ycgDSM.qp.4)
# DE = 8.11, R2=0.006. dstStlmnt still not sig. dstwater & dstStation may be linear

# remove dstStlmnt
ycgDSM.qp.5 <- dsm(Nhat ~ s(dstWater,bs="ts",k=4)+ s(dstBorder,bs="ts",k=5)+ 
                     s(dstStation,bs="ts",k=5)+ habitat,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.qp.5)
par(mfrow=c(2,2))
plot(ycgDSM.qp.5, scale = 0)
gam.check(ycgDSM.qp.5)
# DE = 7.85, R2=0.006. All terms sig

# increase k for dstWater and dstStation
ycgDSM.qp.6 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstBorder,bs="ts",k=4)+ 
                     s(dstStation,bs="ts",k=5)+ habitat,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.qp.6)
par(mfrow=c(2,2))
plot(ycgDSM.qp.6, scale = 0)
gam.check(ycgDSM.qp.6)
# DE = 8.08, R2=0.006. now only dststation looks linear

# dstStation as parametric term
ycgDSM.qp.7 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstBorder,bs="ts",k=4)+ 
                          habitat + dstStation,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.qp.7)
par(mfrow=c(2,2))
plot(ycgDSM.qp.7, scale = 0)
gam.check(ycgDSM.qp.7)
# DE = 8.08, R2=0.006. term still sig, no other difference


# ycgDSM.qp.7 is the best QP model


  ## Tweedie response ####

ycgDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                            habitat,
                  ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.tw.sat)
dev.off()
par(mfrow=c(2,3))
plot(ycgDSM.tw.sat, scale = 0)
gam.check(ycgDSM.tw.sat)
# AIC = 14096, DE = 8.53. dstStlmnt & elevation not sig. dstWater, dstStation potentially linear

# increase k for all except dstBorder
ycgDSM.tw.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+ 
                       s(dstBorder,bs="ts")+ s(dstStation,bs="ts",k=15)+ s(elevation,bs="ts",k=15)+ 
                       habitat,
                     ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                     family = tw(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.tw.sat2)
dev.off()
par(mfrow=c(2,3))
plot(ycgDSM.tw.sat2, scale = 0)
gam.check(ycgDSM.tw.sat2)
# AIC = 14103, DE = 9.97. sig gone down for dstWater and dstStation but up for elevation

# reduce k for dstWater and dstStation and remove dtStlmnt
ycgDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=10)+  s(dstBorder,bs="ts")+ 
                        s(dstStation,bs="ts",k=10)+ s(elevation,bs="ts",k=15)+ 
                        habitat,
                      ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                      family = tw(), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.tw.3)
dev.off()
par(mfrow=c(2,2))
plot(ycgDSM.tw.3, scale = 0)
gam.check(ycgDSM.tw.3)
# AIC = 14096, DE = 8.53. elevation not sig. This suggest dstWater and dstStation are linear relatinships

# increase k for elevaiton
ycgDSM.tw.4 <- dsm(Nhat ~ s(dstWater,bs="ts",k=10)+  s(dstBorder,bs="ts")+ 
                     s(dstStation,bs="ts",k=10)+ s(elevation,bs="ts",k=20)+ 
                     habitat,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.tw.4)
dev.off()
par(mfrow=c(2,2))
plot(ycgDSM.tw.4, scale = 0)
gam.check(ycgDSM.tw.4)
# AIC = 14104, DE = 10.1. increasing k for elevaiton makes it sig, but removes importance of dstStation and dstWater. the model is worse according to AIC, but has more DE. out of elevaiton, dstStation, dstWater, I think dstWater is ecologically the most irrelevant and so I will test its removal

# remove dstWater
ycgDSM.tw.5 <- dsm(Nhat ~ s(dstBorder,bs="ts")+ s(dstStation,bs="ts",k=10)+ 
                     s(elevation,bs="ts",k=15)+ habitat,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.tw.5)
dev.off()
par(mfrow=c(2,2))
plot(ycgDSM.tw.5, scale = 0)
gam.check(ycgDSM.tw.5)
# AIC = 14094, DE = 8.53. 

# same as tw3 but with dstWater and dstStation as linear terms and elevation removed
ycgDSM.tw.6 <- dsm(Nhat ~  s(dstBorder,bs="ts")+ 
                     habitat + dstWater + dstStation,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.tw.6)
dev.off()
par(mfrow=c(2,2))
plot(ycgDSM.tw.6, scale = 0)
gam.check(ycgDSM.tw.6)
# AIC = 14094, DE = 8.53. Better model than tw3 based on AIC. Both term sig as linear.


### ycgDSM.tw.6 is the best tweedie model 

  ## Negative binomial response ####

ycgDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                            habitat,
                     ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.nb.sat)
dev.off()
par(mfrow=c(2,3))
plot(ycgDSM.nb.sat, scale = 0)
gam.check(ycgDSM.nb.sat)
# AIC = 2156, DE = 15.5. dstStlmnt and elevation not sig. dstWater & dstStation linear

# remove dstStlmnt
ycgDSM.nb.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstBorder,bs="ts")+ 
                          s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                          habitat,
                     ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                     family = nb(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.nb.2)
dev.off()
par(mfrow=c(2,2))
plot(ycgDSM.nb.2, scale = 0)
gam.check(ycgDSM.nb.2)
# AIC = 2156, DE = 15.5. no change

# increase k for elevation
ycgDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstBorder,bs="ts")+ 
                     s(dstStation,bs="ts")+ s(elevation,bs="ts",k=15)+ 
                     habitat,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.nb.3)
dev.off()
par(mfrow=c(2,2))
plot(ycgDSM.nb.3, scale = 0)
gam.check(ycgDSM.nb.3)
# AIC = 2156, DE = 15.5. no difference

# remove elevation
ycgDSM.nb.4 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstBorder,bs="ts")+ 
                     s(dstStation,bs="ts")+ 
                     habitat,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.nb.4)
dev.off()
par(mfrow=c(2,2))
plot(ycgDSM.nb.4, scale = 0)
gam.check(ycgDSM.nb.4)
# AIC = 2156, DE = 15.5. No change. dstWater and dtStation linear

# make dstWater & dstStation linear terms
ycgDSM.nb.5 <- dsm(Nhat ~ s(dstBorder,bs="ts")+ 
                          habitat + dstWater + dstStation,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.nb.5)
dev.off()
par(mfrow=c(2,2))
plot(ycgDSM.nb.5, scale = 0)
gam.check(ycgDSM.nb.5)
# AIC = 2154, DE = 6. Model ran with warning. AIC gone down but DE dropped massively

# make dstWater smooth term
ycgDSM.nb.6 <- dsm(Nhat ~ s(dstBorder,bs="ts")+ s(dstWater,bs="ts")+
                     habitat + dstStation,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.nb.6)
dev.off()
par(mfrow=c(2,1))
plot(ycgDSM.nb.6, scale = 0)
gam.check(ycgDSM.nb.6)
# AIC = 2154, DE = 15.5. Better. No warning, AIC as above but DE back up to what it was.



### ycgDSM.nb.6 is the best NB model


## Model selection ####

# Best QP model is ycgDSM.qp.7
# Best TW model is ycgDSM.tw.6
# Best NB model is ycgDSM.nb.6

summary(ycgDSM.qp.7) # DE = 8.08
summary(ycgDSM.tw.6) # DE = 8.53
summary(ycgDSM.nb.6) # DE = 15.5


ycgDSM.tw.6$aic
ycgDSM.nb.6$aic
# The NB model is better than the TW model

# Compare NB and QP models
anova(ycgDSM.qp.7,ycgDSM.nb.6, test="Chisq")
# the NB model has less residual deviance

# Compare diagnostic plots
par(mfrow=c(2,2))
gam.check(ycgDSM.qp.7)
gam.check(ycgDSM.nb.6)
# NB model much better

# ycgDSM.nb.6 is the best model

## Autocorrelation ####

par(mfrow=c(1,1))
dsm.cor(ycgDSM.nb.6, max.lag=15, Segment.Label="Sample.Label")
dsm.cor(ycgDSM.nb.6x.y, max.lag=15, Segment.Label="Sample.Label")
# adding univariate smooths improves the model (AIC and DE), and should account for the autocorrelation.

# some very minor correlation at lag 3, 10, 15.

# add bivariate smooth
ycgDSM.nb.6xy <- dsm(Nhat ~ s(dstBorder,bs="ts")+ s(dstWater,bs="ts")+
                            s(x,y,bs="ts")+
                            habitat + dstStation,
                   ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.nb.6xy)
dev.off()
par(mfrow=c(2,1))
plot(ycgDSM.nb.6xy, scale = 0)
gam.check(ycgDSM.nb.6xy)
# AIC = 2154, DE = 15.8. Correlation increases

# add uivariate smooths
ycgDSM.nb.6x.y <- dsm(Nhat ~ s(dstBorder,bs="ts")+ s(dstWater,bs="ts")+
                       s(x,bs="ts")+ s(y,bs="ts")+
                       habitat + dstStation,
                     ycgDF.hn.cov5, segdata, ycg_obsdata, method = "REML",
                     family = nb(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(ycgDSM.nb.6x.y)
dev.off()
par(mfrow=c(2,1))
plot(ycgDSM.nb.6x.y, scale = 0)
gam.check(ycgDSM.nb.6x.y)
# AIC = 2140, DE = 18.5. 


# ycgDSM.nb.6x.y is the final model

## Abundance estimation ####

# Predict over 2010 habitat
ycg.global.pred10.core <- predict(ycgDSM.nb.6x.y, preddata10_core, off.set = 40000)
write.csv(ycg.global.pred10.core, file="Results/YCG/core_only/ycg.pred10.core.csv")

# Predict over 2011 habitat
ycg.global.pred11.core <- predict(ycgDSM.nb.6x.y, preddata11_core, off.set = 40000)
write.csv(ycg.global.pred11.core, file="Results/YCG/core_only/ycg.pred11.core.csv")

# Predict over 2013 habitat
ycg.global.pred13.core <- predict(ycgDSM.nb.6x.y, preddata13_core, off.set = 40000)
write.csv(ycg.global.pred13.core, file="Results/YCG/core_only/ycg.pred13.core.csv")

# Predict over 2014 habitat
ycg.global.pred14.core <- predict(ycgDSM.nb.6x.y, preddata14_core, off.set = 40000)
write.csv(ycg.global.pred14.core, file="Results/YCG/core_only/ycg.pred14.core.csv")

# Predict over 2016 habitat
ycg.global.pred16.core <- predict(ycgDSM.nb.6x.y, preddata16_core, off.set = 40000)
write.csv(ycg.global.pred16.core, file="Results/YCG/core_only/ycg.pred16.core.csv")

# Predict over 2018 habitat
ycg.global.pred18.core <- predict(ycgDSM.nb.6x.y, preddata18_core, off.set = 40000)
write.csv(ycg.global.pred18.core, file="Results/YCG/core_only/ycg.pred18.core.csv")

# Predict over 2020 habitat
ycg.global.pred20.core <- predict(ycgDSM.nb.6x.y, preddata20_core, off.set = 40000)
write.csv(ycg.global.pred20.core, file="Results/YCG/core_only/ycg.pred20.core.csv")


# Create new dataframes for plotting
ycg.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = ycg.global.pred10.core)

ycg.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = ycg.global.pred11.core)

ycg.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = ycg.global.pred13.core)

ycg.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = ycg.global.pred14.core)

ycg.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = ycg.global.pred16.core)

ycg.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = ycg.global.pred18.core)

ycg.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = ycg.global.pred20.core)


## This creates a dataframe that can be plotted as a map
ycg.spdf.df_10.core <- abunPlotDF(ycg.df.Final10.core, pred.polys_200)
ycg.spdf.df_11.core <- abunPlotDF(ycg.df.Final11.core, pred.polys_200)
ycg.spdf.df_13.core <- abunPlotDF(ycg.df.Final13.core, pred.polys_200)
ycg.spdf.df_14.core <- abunPlotDF(ycg.df.Final14.core, pred.polys_200)
ycg.spdf.df_16.core <- abunPlotDF(ycg.df.Final16.core, pred.polys_200)
ycg.spdf.df_18.core <- abunPlotDF(ycg.df.Final18.core, pred.polys_200)
ycg.spdf.df_20.core <- abunPlotDF(ycg.df.Final20.core, pred.polys_200)

# save SPDFs
write.csv(ycg.spdf.df_10.core,file="Results/YCG/Plots/core_only/spdf/ycg.spdf.df_10.core.csv")
write.csv(ycg.spdf.df_11.core,file="Results/YCG/Plots/core_only/spdf/ycg.spdf.df_11.core.csv")
write.csv(ycg.spdf.df_13.core,file="Results/YCG/Plots/core_only/spdf/ycg.spdf.df_13.core.csv")
write.csv(ycg.spdf.df_14.core,file="Results/YCG/Plots/core_only/spdf/ycg.spdf.df_14.core.csv")
write.csv(ycg.spdf.df_16.core,file="Results/YCG/Plots/core_only/spdf/ycg.spdf.df_16.core.csv")
write.csv(ycg.spdf.df_18.core,file="Results/YCG/Plots/core_only/spdf/ycg.spdf.df_18.core.csv")
write.csv(ycg.spdf.df_20.core,file="Results/YCG/Plots/core_only/spdf/ycg.spdf.df_20.core.csv")



    ## Plotting continuous ####

# Load spatial dataframes
ycg.spdf.df_10.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_10.core.csv")
#ycg.spdf.df_11.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_11.core.csv")
#ycg.spdf.df_13.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_13.core.csv")
#ycg.spdf.df_14.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_14.core.csv")
#ycg.spdf.df_16.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_16.core.csv")
#ycg.spdf.df_18.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_18.core.csv")
ycg.spdf.df_20.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_20.core.csv")

# greyscale plots
YCG_plot_10_core_gr <- GSplotFun(ycg.spdf.df_10.core, survey.area.core, "abundance", "2010")
YCG_plot_11_core_gr <- GSplotFun(ycg.spdf.df_11.core, survey.area.core, "abundance")
YCG_plot_13_core_gr <- GSplotFun(ycg.spdf.df_13.core, survey.area.core, "abundance")
YCG_plot_14_core_gr <- GSplotFun(ycg.spdf.df_14.core, survey.area.core, "abundance")
YCG_plot_16_core_gr <- GSplotFun(ycg.spdf.df_16.core, survey.area.core, "abundance")
YCG_plot_18_core_gr <- GSplotFun(ycg.spdf.df_18.core, survey.area.core, "abundance")
YCG_plot_20_core_gr <- GSplotFun(ycg.spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(YCG_plot_10_core_gr,"Results/YCG/Plots/core_only/greyscale/YCG_plot_10_core_gr.png")
saveplot(YCG_plot_11_core_gr,"Results/YCG/Plots/core_only/greyscale/YCG_plot_11_core_gr.png")
saveplot(YCG_plot_13_core_gr,"Results/YCG/Plots/core_only/greyscale/YCG_plot_13_core_gr.png")
saveplot(YCG_plot_14_core_gr,"Results/YCG/Plots/core_only/greyscale/YCG_plot_14_core_gr.png")
saveplot(YCG_plot_16_core_gr,"Results/YCG/Plots/core_only/greyscale/YCG_plot_16_core_gr.png")
saveplot(YCG_plot_18_core_gr,"Results/YCG/Plots/core_only/greyscale/YCG_plot_18_core_gr.png")
saveplot(YCG_plot_20_core_gr,"Results/YCG/Plots/core_only/greyscale/YCG_plot_20_core_gr.png")

# colour plots
YCG_plot_10_core_col <- CLplotFun(ycg.spdf.df_10.core, survey.area.core, "abundance")
YCG_plot_11_core_col <- CLplotFun(ycg.spdf.df_11.core, survey.area.core, "abundance")
YCG_plot_13_core_col <- CLplotFun(ycg.spdf.df_13.core, survey.area.core, "abundance")
YCG_plot_14_core_col <- CLplotFun(ycg.spdf.df_14.core, survey.area.core, "abundance")
YCG_plot_16_core_col <- CLplotFun(ycg.spdf.df_16.core, survey.area.core, "abundance")
YCG_plot_18_core_col <- CLplotFun(ycg.spdf.df_18.core, survey.area.core, "abundance")
YCG_plot_20_core_col <- CLplotFun(ycg.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(YCG_plot_10_core_col,"Results/YCG/Plots/core_only/colour/YCG_plot_10_core_col.png")
saveplot(YCG_plot_11_core_col,"Results/YCG/Plots/core_only/colour/YCG_plot_11_core_col.png")
saveplot(YCG_plot_13_core_col,"Results/YCG/Plots/core_only/colour/YCG_plot_13_core_col.png")
saveplot(YCG_plot_14_core_col,"Results/YCG/Plots/core_only/colour/YCG_plot_14_core_col.png")
saveplot(YCG_plot_16_core_col,"Results/YCG/Plots/core_only/colour/YCG_plot_16_core_col.png")
saveplot(YCG_plot_18_core_col,"Results/YCG/Plots/core_only/colour/YCG_plot_18_core_col.png")
saveplot(YCG_plot_20_core_col,"Results/YCG/Plots/core_only/colour/YCG_plot_20_core_col.png")



## plot grids (abundance and variance)

# greyscale 
ycg_2yrs_gs <- 
  YCG_plot_10_core_gr + YCG_varplot_final10.core.bw  +
  YCG_plot_20_core_gr + YCG_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
ycg_2yrs_gs[[1]] <- ycg_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
ycg_2yrs_gs[[2]] <- ycg_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
ycg_2yrs_gs[[2]] <- ycg_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
ycg_2yrs_gs[[4]] <- ycg_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(ycg_2yrs_gs, "Results/YCG/Plots/core_only/greyscale/plot_grids/ycg_2yrs_gs.png")


    ## Plotting discrete bins ####
      # Add bins to SPDF - don't repeat ####

# this is the process of adding discrete bins to the abundance SPDFs. I will save them in a new folder so this only has to be done once

# Load original spatial dataframes
ycg.spdf.df_10.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_10.core.csv")
ycg.spdf.df_11.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_11.core.csv")
ycg.spdf.df_13.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_13.core.csv")
ycg.spdf.df_14.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_14.core.csv")
ycg.spdf.df_16.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_16.core.csv")
ycg.spdf.df_18.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_18.core.csv")
ycg.spdf.df_20.core <- read.csv("Results/YCG/Plots/core_only/spdf/ycg.spdf.df_20.core.csv")

# put spdf's into a list
dfs <- list(ycg.spdf.df_10.core,ycg.spdf.df_11.core,ycg.spdf.df_13.core,ycg.spdf.df_14.core,
            ycg.spdf.df_16.core,ycg.spdf.df_18.core,ycg.spdf.df_20.core)

# name the elements
names(dfs) <- c("ycg.spdf.df_10.core","ycg.spdf.df_11.core","ycg.spdf.df_13.core","ycg.spdf.df_14.core",
                "ycg.spdf.df_16.core","ycg.spdf.df_18.core","ycg.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunAbun)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(ycg.spdf.df_10.core,file="Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_10.core.csv")
write.csv(ycg.spdf.df_11.core,file="Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_11.core.csv")
write.csv(ycg.spdf.df_13.core,file="Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_13.core.csv")
write.csv(ycg.spdf.df_14.core,file="Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_14.core.csv")
write.csv(ycg.spdf.df_16.core,file="Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_16.core.csv")
write.csv(ycg.spdf.df_18.core,file="Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_18.core.csv")
write.csv(ycg.spdf.df_20.core,file="Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_20.core.csv")


      # Plotting ####

# Load spatial dataframes (with bins)
ycg.spdf.df_10.core <- read.csv("Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_10.core.csv")
ycg.spdf.df_11.core <- read.csv("Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_11.core.csv")
ycg.spdf.df_13.core <- read.csv("Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_13.core.csv")
ycg.spdf.df_14.core <- read.csv("Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_14.core.csv")
ycg.spdf.df_16.core <- read.csv("Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_16.core.csv")
ycg.spdf.df_18.core <- read.csv("Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_18.core.csv")
ycg.spdf.df_20.core <- read.csv("Results/YCG/Plots/core_only/spdf/bins/ycg.spdf.df_20.core.csv")

## plot greyscale
YCG_10_plot_bin_GS <- GSplotBin(ycg.spdf.df_10.core,"group2",survey.area.core,"Abundance","Relative abundance")
YCG_11_plot_bin_GS <- GSplotBin(ycg.spdf.df_11.core,"group2",survey.area.core,"Abundance","Relative abundance")
YCG_13_plot_bin_GS <- GSplotBin(ycg.spdf.df_13.core,"group2",survey.area.core,"Abundance","Relative abundance")
YCG_14_plot_bin_GS <- GSplotBin(ycg.spdf.df_14.core,"group2",survey.area.core,"Abundance","Relative abundance")
YCG_16_plot_bin_GS <- GSplotBin(ycg.spdf.df_16.core,"group2",survey.area.core,"Abundance","Relative abundance")
YCG_18_plot_bin_GS <- GSplotBin(ycg.spdf.df_18.core,"group2",survey.area.core,"Abundance","Relative abundance")
YCG_20_plot_bin_GS <- GSplotBin(ycg.spdf.df_20.core,"group2",survey.area.core,"Abundance","Relative abundance")


# save 
saveplot(YCG_10_plot_bin_GS,"Results/YCG/Plots/core_only/bins/YCG_10_plot_bin_GS.png")
saveplot(YCG_11_plot_bin_GS,"Results/YCG/Plots/core_only/bins/YCG_11_plot_bin_GS.png")
saveplot(YCG_13_plot_bin_GS,"Results/YCG/Plots/core_only/bins/YCG_13_plot_bin_GS.png")
saveplot(YCG_14_plot_bin_GS,"Results/YCG/Plots/core_only/bins/YCG_14_plot_bin_GS.png")
saveplot(YCG_16_plot_bin_GS,"Results/YCG/Plots/core_only/bins/YCG_16_plot_bin_GS.png")
saveplot(YCG_18_plot_bin_GS,"Results/YCG/Plots/core_only/bins/YCG_18_plot_bin_GS.png")
saveplot(YCG_20_plot_bin_GS,"Results/YCG/Plots/core_only/bins/YCG_20_plot_bin_GS.png")



## Variance estimation ####
  

# estimate variance
ycg.var.Final10.core <- varEstfun(preddata10_core, ycgDSM.nb.6x.y)
ycg.var.Final11.core <- varEstfun(preddata11_core, ycgDSM.nb.6x.y)
ycg.var.Final13.core <- varEstfun(preddata13_core, ycgDSM.nb.6x.y)
ycg.var.Final14.core <- varEstfun(preddata14_core, ycgDSM.nb.6x.y)
ycg.var.Final16.core <- varEstfun(preddata16_core, ycgDSM.nb.6x.y)
ycg.var.Final18.core <- varEstfun(preddata18_core, ycgDSM.nb.6x.y)
ycg.var.Final20.core <- varEstfun(preddata20_core, ycgDSM.nb.6x.y)

# save variance estimates
write.csv(ycg.var.Final10.core, file="Results/YCG/core_only/ycg.var10.core.csv")
write.csv(ycg.var.Final11.core, file="Results/YCG/core_only/ycg.var11.core.csv")
write.csv(ycg.var.Final13.core, file="Results/YCG/core_only/ycg.var13.core.csv")
write.csv(ycg.var.Final14.core, file="Results/YCG/core_only/ycg.var14.core.csv")
write.csv(ycg.var.Final16.core, file="Results/YCG/core_only/ycg.var16.core.csv")
write.csv(ycg.var.Final18.core, file="Results/YCG/core_only/ycg.var18.core.csv")
write.csv(ycg.var.Final20.core, file="Results/YCG/core_only/ycg.var20.core.csv")

# create spdf's for plotting
ycg.var.spdf.df_10.core <- varPlotDF(ycg.var.Final10.core, pred.polys_200)
ycg.var.spdf.df_11.core <- varPlotDF(ycg.var.Final11.core, pred.polys_200)
ycg.var.spdf.df_13.core <- varPlotDF(ycg.var.Final13.core, pred.polys_200)
ycg.var.spdf.df_14.core <- varPlotDF(ycg.var.Final14.core, pred.polys_200)
ycg.var.spdf.df_16.core <- varPlotDF(ycg.var.Final16.core, pred.polys_200)
ycg.var.spdf.df_18.core <- varPlotDF(ycg.var.Final18.core, pred.polys_200)
ycg.var.spdf.df_20.core <- varPlotDF(ycg.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(ycg.var.spdf.df_10.core,
          file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_10.core.csv")
write.csv(ycg.var.spdf.df_11.core,
          file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_11.core.csv")
write.csv(ycg.var.spdf.df_13.core,
          file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_13.core.csv")
write.csv(ycg.var.spdf.df_14.core,
          file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_14.core.csv")
write.csv(ycg.var.spdf.df_16.core,
          file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_16.core.csv")
write.csv(ycg.var.spdf.df_18.core,
          file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_18.core.csv")
write.csv(ycg.var.spdf.df_20.core,
          file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_20.core.csv")



    # Calculate CV & add bins to SPDF - don't repeat ####

# Load spatial dataframes
ycg.var.spdf.df_10.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_10.core.csv")
ycg.var.spdf.df_11.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_11.core.csv")
ycg.var.spdf.df_13.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_13.core.csv")
ycg.var.spdf.df_14.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_14.core.csv")
ycg.var.spdf.df_16.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_16.core.csv")
ycg.var.spdf.df_18.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_18.core.csv")
ycg.var.spdf.df_20.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_20.core.csv")


# first need to calculate CV from the variance (note: need the abundance spdf's loaded)
ycg.var.spdf.df_10.core <- CVaddFun(ycg.spdf.df_10.core,ycg.var.spdf.df_10.core)
ycg.var.spdf.df_11.core <- CVaddFun(ycg.spdf.df_11.core,ycg.var.spdf.df_11.core)
ycg.var.spdf.df_13.core <- CVaddFun(ycg.spdf.df_13.core,ycg.var.spdf.df_13.core)
ycg.var.spdf.df_14.core <- CVaddFun(ycg.spdf.df_14.core,ycg.var.spdf.df_14.core)
ycg.var.spdf.df_16.core <- CVaddFun(ycg.spdf.df_16.core,ycg.var.spdf.df_16.core)
ycg.var.spdf.df_18.core <- CVaddFun(ycg.spdf.df_18.core,ycg.var.spdf.df_18.core)
ycg.var.spdf.df_20.core <- CVaddFun(ycg.spdf.df_20.core,ycg.var.spdf.df_20.core)


# Save the SPDF's with the CV value but no bins (for continuous plotting of the CV)
write.csv(ycg.var.spdf.df_10.core,file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_10.core.csv")
write.csv(ycg.var.spdf.df_11.core,file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_11.core.csv")
write.csv(ycg.var.spdf.df_13.core,file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_13.core.csv")
write.csv(ycg.var.spdf.df_14.core,file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_14.core.csv")
write.csv(ycg.var.spdf.df_16.core,file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_16.core.csv")
write.csv(ycg.var.spdf.df_18.core,file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_18.core.csv")
write.csv(ycg.var.spdf.df_20.core,file="Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_20.core.csv")


### add bins

## Quartiles

# put spdf's into a list
dfs <- list(ycg.var.spdf.df_10.core,ycg.var.spdf.df_11.core,ycg.var.spdf.df_13.core,ycg.var.spdf.df_14.core,
            ycg.var.spdf.df_16.core,ycg.var.spdf.df_18.core,ycg.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("ycg.var.spdf.df_10.core","ycg.var.spdf.df_11.core","ycg.var.spdf.df_13.core",
                "ycg.var.spdf.df_14.core","ycg.var.spdf.df_16.core","ycg.var.spdf.df_18.core",
                "ycg.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(ycg.var.spdf.df_10.core,file="Results/YCG/Plots/variance/core_only/spdf/bins/ycg.var.spdf.df_10.core.csv")
write.csv(ycg.var.spdf.df_11.core,file="Results/YCG/Plots/variance/core_only/spdf/bins/ycg.var.spdf.df_11.core.csv")
write.csv(ycg.var.spdf.df_13.core,file="Results/YCG/Plots/variance/core_only/spdf/bins/ycg.var.spdf.df_13.core.csv")
write.csv(ycg.var.spdf.df_14.core,file="Results/YCG/Plots/variance/core_only/spdf/bins/ycg.var.spdf.df_14.core.csv")
write.csv(ycg.var.spdf.df_16.core,file="Results/YCG/Plots/variance/core_only/spdf/bins/ycg.var.spdf.df_16.core.csv")
write.csv(ycg.var.spdf.df_18.core,file="Results/YCG/Plots/variance/core_only/spdf/bins/ycg.var.spdf.df_18.core.csv")
write.csv(ycg.var.spdf.df_20.core,file="Results/YCG/Plots/variance/core_only/spdf/bins/ycg.var.spdf.df_20.core.csv")




## custom bins

# change column name from group2 to CV
ycg.var.spdf.df_10.core <- ycg.var.spdf.df_10.core %>% dplyr::rename(CV=group2)
ycg.var.spdf.df_11.core <- ycg.var.spdf.df_11.core %>% dplyr::rename(CV=group2)
ycg.var.spdf.df_13.core <- ycg.var.spdf.df_13.core %>% dplyr::rename(CV=group2)
ycg.var.spdf.df_14.core <- ycg.var.spdf.df_14.core %>% dplyr::rename(CV=group2)
ycg.var.spdf.df_16.core <- ycg.var.spdf.df_16.core %>% dplyr::rename(CV=group2)
ycg.var.spdf.df_18.core <- ycg.var.spdf.df_18.core %>% dplyr::rename(CV=group2)
ycg.var.spdf.df_20.core <- ycg.var.spdf.df_20.core %>% dplyr::rename(CV=group2)


# put spdf's into a list
dfs <- list(ycg.var.spdf.df_10.core,ycg.var.spdf.df_11.core,ycg.var.spdf.df_13.core,ycg.var.spdf.df_14.core,
            ycg.var.spdf.df_16.core,ycg.var.spdf.df_18.core,ycg.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("ycg.var.spdf.df_10.core","ycg.var.spdf.df_11.core","ycg.var.spdf.df_13.core",
                "ycg.var.spdf.df_14.core","ycg.var.spdf.df_16.core","ycg.var.spdf.df_18.core",
                "ycg.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar2)

# split elements into original dataframes
list2env(dfs, globalenv())


# save the SPDFs with the custom bins
write.csv(ycg.var.spdf.df_10.core,
          file="Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_10.core.csv")
write.csv(ycg.var.spdf.df_11.core,
          file="Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_11.core.csv")
write.csv(ycg.var.spdf.df_13.core,
          file="Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_13.core.csv")
write.csv(ycg.var.spdf.df_14.core,
          file="Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_14.core.csv")
write.csv(ycg.var.spdf.df_16.core,
          file="Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_16.core.csv")
write.csv(ycg.var.spdf.df_18.core,
          file="Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_18.core.csv")
write.csv(ycg.var.spdf.df_20.core,
          file="Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_20.core.csv")


    # Discrete bins ####

# load spdfs (which has already had CV calculated and then put into bins)
ycg.var.spdf.df_10.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_10.core.csv")
ycg.var.spdf.df_11.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_11.core.csv")
ycg.var.spdf.df_13.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_13.core.csv")
ycg.var.spdf.df_14.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_14.core.csv")
ycg.var.spdf.df_16.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_16.core.csv")
ycg.var.spdf.df_18.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_18.core.csv")
ycg.var.spdf.df_20.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/bins/custom/ycg.var.spdf.df_20.core.csv")

# make CV a factor and re-order. I have only done this for 2020 but you can copy the code for the other years if you need to
ycg.var.spdf.df_20.core$CV <- as.factor(ycg.var.spdf.df_20.core$CV)
ycg.var.spdf.df_20.core$CV <- factor(ycg.var.spdf.df_20.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))

# plot CV in bins
YCG_10_plot_bin_GS_var <- GSplotBin(ycg.var.spdf.df_10.core,"CV",survey.area.core,"Variance","CV")
YCG_11_plot_bin_GS_var <- GSplotBin(ycg.var.spdf.df_11.core,"CV",survey.area.core,"Variance","CV")
YCG_13_plot_bin_GS_var <- GSplotBin(ycg.var.spdf.df_13.core,"CV",survey.area.core,"Variance","CV")
YCG_14_plot_bin_GS_var <- GSplotBin(ycg.var.spdf.df_14.core,"CV",survey.area.core,"Variance","CV")
YCG_16_plot_bin_GS_var <- GSplotBin(ycg.var.spdf.df_16.core,"CV",survey.area.core,"Variance","CV")
YCG_18_plot_bin_GS_var <- GSplotBin(ycg.var.spdf.df_18.core,"CV",survey.area.core,"Variance","CV")
YCG_20_plot_bin_GS_var <- GSplotBin(ycg.var.spdf.df_20.core,"CV",survey.area.core,"Variance","CV")

# save 
saveplot(YCG_10_plot_bin_GS_var,"Results/YCG/Plots/variance/core_only/greyscale/bins/YCG_10_plot_bin_GS_var.png")
saveplot(YCG_11_plot_bin_GS_var,"Results/YCG/Plots/variance/core_only/greyscale/bins/YCG_11_plot_bin_GS_var.png")
saveplot(YCG_13_plot_bin_GS_var,"Results/YCG/Plots/variance/core_only/greyscale/bins/YCG_13_plot_bin_GS_var.png")
saveplot(YCG_14_plot_bin_GS_var,"Results/YCG/Plots/variance/core_only/greyscale/bins/YCG_14_plot_bin_GS_var.png")
saveplot(YCG_16_plot_bin_GS_var,"Results/YCG/Plots/variance/core_only/greyscale/bins/YCG_16_plot_bin_GS_var.png")
saveplot(YCG_18_plot_bin_GS_var,"Results/YCG/Plots/variance/core_only/greyscale/bins/YCG_18_plot_bin_GS_var.png")
saveplot(YCG_20_plot_bin_GS_var,"Results/YCG/Plots/variance/core_only/greyscale/bins/YCG_20_plot_bin_GS_var.png")



    # Plotting continuous ####

# Load spatial dataframes
ycg.var.spdf.df_10.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_10.core.csv")
ycg.var.spdf.df_11.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_11.core.csv")
ycg.var.spdf.df_13.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_13.core.csv")
ycg.var.spdf.df_14.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_14.core.csv")
ycg.var.spdf.df_16.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_16.core.csv")
ycg.var.spdf.df_18.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_18.core.csv")
ycg.var.spdf.df_20.core <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_20.core.csv")


# greyscale plots
YCG_varplot_final10.core.bw <- GSplotFun(ycg.var.spdf.df_10.core, survey.area.core, "variance", "2010")
YCG_varplot_final11.core.bw <- GSplotFun(ycg.var.spdf.df_11.core, survey.area.core, "variance")
YCG_varplot_final13.core.bw <- GSplotFun(ycg.var.spdf.df_13.core, survey.area.core, "variance")
YCG_varplot_final14.core.bw <- GSplotFun(ycg.var.spdf.df_14.core, survey.area.core, "variance")
YCG_varplot_final16.core.bw <- GSplotFun(ycg.var.spdf.df_16.core, survey.area.core, "variance")
YCG_varplot_final18.core.bw <- GSplotFun(ycg.var.spdf.df_18.core, survey.area.core, "variance")
YCG_varplot_final20.core.bw <- GSplotFun(ycg.var.spdf.df_20.core, survey.area.core, "variance", "2020")

# save greyscale
saveplot(YCG_varplot_final10.core.bw, "Results/YCG/Plots/variance/core_only/greyscale/2010_YCG_var.core.bw.png")
saveplot(YCG_varplot_final11.core.bw, "Results/YCG/Plots/variance/core_only/greyscale/2011_YCG_var.core.bw.png")
saveplot(YCG_varplot_final13.core.bw, "Results/YCG/Plots/variance/core_only/greyscale/2013_YCG_var.core.bw.png")
saveplot(YCG_varplot_final14.core.bw, "Results/YCG/Plots/variance/core_only/greyscale/2014_YCG_var.core.bw.png")
saveplot(YCG_varplot_final16.core.bw, "Results/YCG/Plots/variance/core_only/greyscale/2016_YCG_var.core.bw.png")
saveplot(YCG_varplot_final18.core.bw, "Results/YCG/Plots/variance/core_only/greyscale/2018_YCG_var.core.bw.png")
saveplot(YCG_varplot_final20.core.bw, "Results/YCG/Plots/variance/core_only/greyscale/2020_YCG_var.core.bw.png")

# colour plots
YCG_varplot_final10.core.col <- CLplotFun(ycg.var.spdf.df_10.core, survey.area.core, "variance")
YCG_varplot_final11.core.col <- CLplotFun(ycg.var.spdf.df_11.core, survey.area.core, "variance")
YCG_varplot_final13.core.col <- CLplotFun(ycg.var.spdf.df_13.core, survey.area.core, "variance")
YCG_varplot_final14.core.col <- CLplotFun(ycg.var.spdf.df_14.core, survey.area.core, "variance")
YCG_varplot_final16.core.col <- CLplotFun(ycg.var.spdf.df_16.core, survey.area.core, "variance")
YCG_varplot_final18.core.col <- CLplotFun(ycg.var.spdf.df_18.core, survey.area.core, "variance")
YCG_varplot_final20.core.col <- CLplotFun(ycg.var.spdf.df_20.core, survey.area.core, "variance")

# save colour
saveplot(YCG_varplot_final10.core.col, "Results/YCG/Plots/variance/core_only/colour/2010_YCG_var.core.col.png")
saveplot(YCG_varplot_final11.core.col, "Results/YCG/Plots/variance/core_only/colour/2011_YCG_var.core.col.png")
saveplot(YCG_varplot_final13.core.col, "Results/YCG/Plots/variance/core_only/colour/2013_YCG_var.core.col.png")
saveplot(YCG_varplot_final14.core.col, "Results/YCG/Plots/variance/core_only/colour/2014_YCG_var.core.col.png")
saveplot(YCG_varplot_final16.core.col, "Results/YCG/Plots/variance/core_only/colour/2016_YCG_var.core.col.png")
saveplot(YCG_varplot_final18.core.col, "Results/YCG/Plots/variance/core_only/colour/2018_YCG_var.core.col.png")
saveplot(YCG_varplot_final20.core.col, "Results/YCG/Plots/variance/core_only/colour/2020_YCG_var.core.col.png")


#### Silver langur ###########################################################
## Load data ####

# Observation data. Unique to species 
gsl_obsdata <- read.csv("Species_Data/GSL/R Data/obsdata.csv", header = TRUE)
gsl_obsdata$object <- as.factor(gsl_obsdata$object)
gsl_obsdata$Sample.Label <- as.factor(gsl_obsdata$Sample.Label)
str(gsl_obsdata)
head(gsl_obsdata)

# Transect data. Unique to species. DF model to be used includes strat and size, and so need to create continuous stratum covar
gsl_distdata <- read.csv("Species_Data/GSL/R Data/distdata.csv", header = TRUE)
gsl_distdata$object <- as.factor(gsl_distdata$object)
gsl_distdata$NameObserver <- as.factor(gsl_distdata$NameObserver)
gsl_distdata$transect <- as.factor(gsl_distdata$transect)
gsl_distdata$stratum <- as.vector(scale(gsl_distdata$year, scale=T, center=T))
gsl_distdata$year <- as.factor(gsl_distdata$year)
gsl_distdata$date <- as.Date(gsl_distdata$date, format = "%d/%m/%Y")
str(gsl_distdata)
head(gsl_distdata)

# check for any observations on T20
gsl_distdata[gsl_distdata$transect=="20",]
# none


## Plot the covariates across the grid, with group sizes ####

# Warning - plots take a few minutes to run

# habitat
plot_GSLobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata200$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=gsl_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/GSL/plot_GSLobs_habitat.png", plot_GSLobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_GSLobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys_200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=gsl_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/GSL/plot_GSLobs_dstWater.png", plot_GSLobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_GSLobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=gsl_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/GSL/plot_GSLobs_dstStlmnt.png", plot_GSLobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_GSLobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=gsl_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/GSL/plot_GSLobs_dstRoad.png", plot_GSLobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_GSLobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=gsl_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/GSL/plot_GSLobs_dstBorder.png", plot_GSLobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_GSLobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=gsl_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/GSL/plot_GSLobs_dstStation.png", plot_GSLobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_GSLobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=gsl_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/GSL/plot_GSLobs_dstELC.png", plot_GSLobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_GSLobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys_200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=gsl_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/GSL/plot_GSLobs_elev.png", plot_GSLobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)

## Exploratory plots & linear models ####

## Checking that there are no large gaps in the range of variables for langur observations.  

# subset segdata to get only the segments with GSL observations
gsl_varcheck <- segdata[match(gsl_obsdata$Sample.Label,segdata$Sample.Label), ]

# habitat - fine
plot(segdata$habitat)
plot(gsl_varcheck$habitat)

# dstWater - large gap between 1000 - 1750. This may be a problem
hist(segdata$dstWater)
hist(gsl_varcheck$dstWater)

# dstStlmnt - gaps between 6000-8000 and 9000 and 13000. Could be a problem
hist(segdata$dstStlmnt)
hist(gsl_varcheck$dstStlmnt)

# dstRoad - fine
hist(segdata$dstRoad)
hist(gsl_varcheck$dstRoad)

# dstBorder - fine
hist(segdata$dstBorder)
hist(gsl_varcheck$dstBorder)

# dstStation - fine
hist(segdata$dstStation)
hist(gsl_varcheck$dstStation)

# elevation - no observations beyond 310 but variable goes up to 750
hist(segdata$elevation)
hist(gsl_varcheck$elevation)

## Histograms

# Distance 
gsl_h1 <- ggplot(gsl_distdata, aes(distance))+ geom_histogram(binwidth = 1)
gsl_h2 <- ggplot(gsl_distdata, aes(distance))+ geom_histogram(binwidth = 5)
gsl_h3 <- ggplot(gsl_distdata, aes(distance))+ geom_histogram(binwidth = 10)
gsl_h4 <- ggplot(gsl_distdata, aes(distance))+ geom_histogram(binwidth = 15)
gsl_h5 <- ggplot(gsl_distdata, aes(distance))+ geom_histogram(binwidth = 20)
gsl_h6 <- ggplot(gsl_distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(gsl_h1,gsl_h2,gsl_h3,gsl_h4,gsl_h5,gsl_h6)
# Slight evidence of evasive movement.  Interesting sharp drop in observations around 30m. Perhaps this is because it is a riverine species, so the canopy is very thick? 


# cluster size, observer, habitat, year, month, transect
gsl_h7 <- ggplot(gsl_distdata, aes(size))+geom_histogram(binwidth = 0.5)
gsl_h8 <- ggplot(gsl_distdata, aes(NameObserver))+geom_histogram(stat="count")
gsl_h9 <- ggplot(gsl_distdata, aes(habitat))+geom_histogram(stat="count")
gsl_h10 <- ggplot(gsl_distdata, aes(year))+geom_histogram(stat="count")
gsl_h11 <- ggplot(gsl_distdata, aes(month))+geom_histogram(stat="count")
gsl_h12 <- ggplot(gsl_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(gsl_h7,gsl_h8,gsl_h9,gsl_h10,gsl_h11,gsl_h12)
# Most observations are of group sizes of 1:7.  This is an underestimate based on work by Jessica Moody.  I will be estimating abundance of group sizes.  Most observations are in dense forest but there are a lot in open forest.  I guess most of these open forest observations are actually in riverine dense forest in wider open forest. Looking at the observations on GIS they are all very close to rivers.

## Plots of distance against variables

plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

gsl_d1 <- ggplot(gsl_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
gsl_d2 <- ggplot(gsl_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
gsl_d3 <- ggplot(gsl_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
gsl_d4 <- ggplot(gsl_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
gsl_d5 <- ggplot(gsl_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(gsl_d1,gsl_d2,gsl_d3,gsl_d4,gsl_d5)
# distances do not seem to change with habitat. This again could be because they are mostly in riverine habitat which is as thick as dense forest. Distance by month seems to increase slightly between Jan-Apr.  It then decreases in May and June.  There doesn't seem to be a relationship between group size and distance


## Plots of cluster size against variables
gsl_s1 <- ggplot(gsl_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
gsl_s2 <- ggplot(gsl_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
gsl_s3 <- ggplot(gsl_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
gsl_s4 <- ggplot(gsl_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
gsl_s5 <- ggplot(gsl_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(gsl_s1,gsl_s2,gsl_s3,gsl_s4,gsl_s5)
# No real relationships except for group size across years.  It looks like group size is incraesing with year. This could be because the team's observation skills improve over time?

## Linear models

# group size ~ distance
newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=gsl_distdata)
plot(gsl_distdata$size~gsl_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)

#  There is a small trend showing an increase group size with increased distance.  But when truncated at 50m this is probably not important.

# test for size bias again but with the 50m truncation distance

gsl_distdata_trunc <- gsl_distdata %>% filter(distance <= 50)
newdist <- data.frame(distance=seq(0,50,len=10))
lm1 <- lm(size~distance, data=gsl_distdata_trunc)
plot(gsl_distdata_trunc$size~gsl_distdata_trunc$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# The effect size is larger when using the truncated distances and so cluster size will be tested in the DF models

## Estimating the detection function ####


### I am taking th DF model from the CDS.  The DF model is hn.size.strat, with cutpoints.

gslDF.hn.size.strat <- ds(gsl_distdata, truncation=50, key="hn", formula=~stratum+size,
                          cutpoints = c(0,10,19,30,50))



 
## Fitting a spatial model ####

# The best detection function model for GSL is gslDF.hn.cov8

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated.  

# I am not including dstELC in the models for GSL, as I don't think it is an appropriate variable, especially now that we are going to be predicting into the buffer zone

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
gsl_obsdata <- gsl_obsdata %>% filter(distance <= 50)


  ## Quasipoisson response ####

gslDSM.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                         s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                         habitat,
                  gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.sat)
dev.off()
par(mfrow=c(2,3))
plot(gslDSM.sat, scale = 0)
gam.check(gslDSM.sat)
# DE = 48.1. All terms sig. All plots look overfitted. 

# reduce k for all
gslDSM.qp.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=5)+ 
                    s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                    s(elevation,bs="ts",k=5)+ 
                    habitat,
                  gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.qp.sat2)
dev.off()
par(mfrow=c(2,3))
plot(gslDSM.qp.sat2, scale = 0)
gam.check(gslDSM.qp.sat2)
# DE = 34. All terms sig. Plots look better. Perhaps still a wee bit overfitted for some

# reduce k more for dstStlmnt and dstStation
gslDSM.qp.sat3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=4)+ 
                        s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=4)+ 
                        s(elevation,bs="ts",k=5)+ 
                        habitat,
                      gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                      family = quasipoisson(link = "log"), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.qp.sat3)
dev.off()
par(mfrow=c(2,3))
plot(gslDSM.qp.sat3, scale = 0)
gam.check(gslDSM.qp.sat3)
# DE = 29.9. dstStlmnt and dstborder no longer sig. dstStation looks better.

# increase k again for dstSltmnt
gslDSM.qp.sat4 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=5)+ 
                        s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=4)+ 
                        s(elevation,bs="ts",k=5)+ 
                        habitat,
                      gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                      family = quasipoisson(link = "log"), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.qp.sat4)
dev.off()
par(mfrow=c(2,3))
plot(gslDSM.qp.sat4, scale = 0)
gam.check(gslDSM.qp.sat4)
# DE = 30. dstSltmnt and dstBorder still not sig. plots didn't really change. 


### gslDSM.qp.sat2 is the best QP model


  ## Tweedie response ####

gslDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ 
                            habitat,
                  gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.tw.sat)
dev.off()
par(mfrow=c(2,3))
plot(gslDSM.tw.sat, scale = 0)
gam.check(gslDSM.tw.sat)
# AIC = 14128, DE = 35.6. dstStlmnt & dstBorder not sig. dstStation possibly linear

# remove dstStlmnt
gslDSM.tw.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts")+ 
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ 
                          habitat,
                     gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                     family = tw(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.tw.2)
dev.off()
par(mfrow=c(2,2))
plot(gslDSM.tw.2, scale = 0)
gam.check(gslDSM.tw.2)
# AIC = 14122, DE = 44.4. All terms sig, DE gone up and AIC gone down. dstBorder overfitted. dstStation looks linear

# reduce k for dstborder
gslDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts")+ 
                          s(dstBorder,bs="ts",k=10)+ s(dstStation,bs="ts")+ 
                          habitat,
                   gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.tw.3)
dev.off()
par(mfrow=c(2,2))
plot(gslDSM.tw.3, scale = 0)
gam.check(gslDSM.tw.3)
#  Varied k for dstBorder from 5 to 10.  It is non-sig flat line up until k=10, after which it is sig and very wiggly.  I am not convinced that is a true ecological relationship - it's just overfitting the points once it gets enough freedom.

# remove dstBorder
gslDSM.tw.4 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts")+ 
                          s(dstStation,bs="ts")+ habitat,
                   gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.tw.4)
dev.off()
par(mfrow=c(2,2))
plot(gslDSM.tw.4, scale = 0)
gam.check(gslDSM.tw.4)
# AIC = 14128, DE = 35.5. ran with warning. All terms sig, but dstStation linear

# dstStation as linear term
gslDSM.tw.5 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts")+ 
                          habitat + dstStation,
                   gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.tw.5)
dev.off()
par(mfrow=c(2,1))
plot(gslDSM.tw.5, scale = 0)
gam.check(gslDSM.tw.5)
# AIC = 14127, DE = 35.5. Still runs with warning. Better model by AIC. Still not sure that elevation needs 5 EDF - its basically a flat line

# reduce k for elevation
gslDSM.tw.6 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts",k=4)+ 
                     habitat + dstStation,
                   gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.tw.6)
dev.off()
par(mfrow=c(2,1))
plot(gslDSM.tw.6, scale = 0)
gam.check(gslDSM.tw.6)
# AIC = 14128, DE = 33.5. ran without warning. Elevation plot looks better. All terms sig



### Best tw model is gslDSM.tw.6  


  ## Negative binomial response ####

gslDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ 
                            habitat,
                     gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                     family = nb(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.nb.sat)
dev.off()
par(mfrow=c(2,3))
plot(gslDSM.nb.sat, scale = 0)
gam.check(gslDSM.nb.sat)
# AIC = 1487, DE = 65.6.All term sig. dstStlmnt and elevation may be linear. All others look overfitted

# reduce k for dstWater, dstBorder, dstStation
gslDSM.nb.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+ 
                       s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                       habitat,
                     gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                     family = nb(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.nb.sat2)
dev.off()
par(mfrow=c(2,3))
plot(gslDSM.nb.sat2, scale = 0)
gam.check(gslDSM.nb.sat2)
# AIC = 1594, DE = 54.2. dstStlmnt no longer sig. elevation and dstBorder linear. Plots for dstWater and dstStation look better

# remove dstStlmnt
gslDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+  s(elevation,bs="ts")+ 
                        s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                        habitat,
                      gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                      family = nb(), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.nb.3)
dev.off()
par(mfrow=c(2,2))
plot(gslDSM.nb.3, scale = 0)
gam.check(gslDSM.nb.3)
# AIC = 1594, DE = 54.2. All term sig. elevation and dstBorder linear.

# dstBorder linear
gslDSM.nb.4 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+  s(elevation,bs="ts")+ 
                          s(dstStation,bs="ts",k=5)+ 
                          habitat + dstBorder,
                   gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.nb.4)
dev.off()
par(mfrow=c(2,2))
plot(gslDSM.nb.4, scale = 0)
gam.check(gslDSM.nb.4)
# AIC = 1592, DE = 54.2. AIC down - better model.

# elevation linear
gslDSM.nb.5 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                          habitat + dstBorder + elevation,
                   gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.nb.5)
dev.off()
par(mfrow=c(1,2))
plot(gslDSM.nb.5, scale = 0)
gam.check(gslDSM.nb.5)
# AIC = 1592, DE = 54.1.



### gslDSM.nb.5 is the best model

  
## Model selection ####

# best QP is gslDSM.qp.sat2
# best TW is gslDSM.tw.6 
# best NB is gslDSM.nb.5

summary(gslDSM.qp.sat2) # DE = 34
summary(gslDSM.tw.6) # DE = 33.5
summary(gslDSM.nb.5) # DE = 54.1

gslDSM.tw.6$aic
gslDSM.nb.5$aic

# NB model is better than TW model (AIC)

# Compare NB and QP
anova(gslDSM.nb.5,gslDSM.qp.sat2,test="Chisq")
# NB has less residual deviance

# Compare diagnostic plots
par(mfrow=c(2,2))
gam.check(gslDSM.qp.sat2)
gam.check(gslDSM.nb.5)
# NB plots are much better

### gslDSM.nb.5 is the best model for GSL


## Autocorrelation ####

par(mfrow=c(1,1))
dsm.cor(gslDSM.nb.5, max.lag=15, Segment.Label="Sample.Label")
# significant correlation at lag 1  


# Add autoregressive structure

# First recode sample labels and transect labels as numeric
segdata$sg.id <- as.numeric(sub("\\d+-","",segdata$Sample.Label))
segdata$tr.id <- as.numeric(segdata$Transect.Label)

gslDSM.nb.5.cor <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                              habitat + elevation + dstBorder,
                   gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML",
                   engine = "gamm", correlation=corAR1(form=~sg.id|tr.id), niterPQL=50,
                   segment.area = segdata$Sample.Area, group = TRUE)
# no convergence even after niterPQL increased to 50

# add univariate smooths of location into the model
gslDSM.nb.5.xy <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+
                             s(x,bs="ts") + s(y,bs="ts")+
                             habitat + elevation + dstBorder,
                   gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML", 
                   family = nb(),engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.nb.5.xy)
dev.off()
par(mfrow=c(2,2))
plot(gslDSM.nb.5.xy, scale = 0)
gam.check(gslDSM.nb.5.xy)
# AIC = 1497, DE = 62.5. x not sig

# try add bivariate smooths
gslDSM.nb.5.xy2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+
                        s(x,y,bs="ts") +
                        habitat + elevation + dstBorder,
                      gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML", 
                      family = nb(),engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.nb.5.xy2)
par(mfrow=c(2,2))
plot(gslDSM.nb.5.xy2, scale = 0)
gam.check(gslDSM.nb.5.xy2)
# AIC = 1473, DE = 66.2. Much better model than above. dstBorder no longer sig as smooth term

# re-add dstBorder as smooth term
gslDSM.nb.5.xy3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+
                         s(x,y,bs="ts") + s(dstBorder,bs="ts",k=10)+
                         habitat + elevation,
                       gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML", 
                       family = nb(),engine = "gam",
                       segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.nb.5.xy3)
par(mfrow=c(2,2))
plot(gslDSM.nb.5.xy3, scale = 0)
gam.check(gslDSM.nb.5.xy3)
# AIC = 1474, DE = 66.1.

# remove dstBorder
gslDSM.nb.5.xy4 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+
                         s(x,y,bs="ts") + 
                         habitat + elevation,
                       gslDF.hn.size.strat, segdata, gsl_obsdata, method = "REML", 
                       family = nb(),engine = "gam",
                       segment.area = segdata$Sample.Area, group = TRUE)
summary(gslDSM.nb.5.xy4)
par(mfrow=c(2,2))
plot(gslDSM.nb.5.xy4, scale = 0)
gam.check(gslDSM.nb.5.xy4)
# AIC = 1474, DE = 66.1. No difference with dstBorder removed

par(mfrow=c(1,1))
dsm.cor(gslDSM.nb.5.xy4, max.lag=15, Segment.Label="Sample.Label")
 


### gslDSM.nb.5.xy4 will be the final model. Autocorrelation has been reduced and should be accounted for with the xy smooths
 
#
## Abundance estimation ####

# Predict over 2010 habitat
gsl.global.pred10.core <- predict(gslDSM.nb.5.xy4, preddata10_core, off.set = 40000)
write.csv(gsl.global.pred10.core, file="Results/GSL/core_only/gsl.pred10.core.csv")

# Predict over 2011 habitat
gsl.global.pred11.core <- predict(gslDSM.nb.5.xy4, preddata11_core, off.set = 40000)
write.csv(gsl.global.pred11.core, file="Results/GSL/core_only/gsl.pred11.core.csv")

# Predict over 2013 habitat
gsl.global.pred13.core <- predict(gslDSM.nb.5.xy4, preddata13_core, off.set = 40000)
write.csv(gsl.global.pred13.core, file="Results/GSL/core_only/gsl.pred13.core.csv")

# Predict over 2014 habitat
gsl.global.pred14.core <- predict(gslDSM.nb.5.xy4, preddata14_core, off.set = 40000)
write.csv(gsl.global.pred14.core, file="Results/GSL/core_only/gsl.pred14.core.csv")

# Predict over 2016 habitat
gsl.global.pred16.core <- predict(gslDSM.nb.5.xy4, preddata16_core, off.set = 40000)
write.csv(gsl.global.pred16.core, file="Results/GSL/core_only/gsl.pred16.core.csv")

# Predict over 2018 habitat
gsl.global.pred18.core <- predict(gslDSM.nb.5.xy4, preddata18_core, off.set = 40000)
write.csv(gsl.global.pred18.core, file="Results/GSL/core_only/gsl.pred18.core.csv")

# Predict over 2020 habitat
gsl.global.pred20.core <- predict(gslDSM.nb.5.xy4, preddata20_core, off.set = 40000)
write.csv(gsl.global.pred20.core, file="Results/GSL/core_only/gsl.pred20.core.csv")



# Create new dataframes for plotting
gsl.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = gsl.global.pred10.core)

gsl.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = gsl.global.pred11.core)

gsl.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = gsl.global.pred13.core)

gsl.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = gsl.global.pred14.core)

gsl.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = gsl.global.pred16.core)

gsl.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = gsl.global.pred18.core)

gsl.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = gsl.global.pred20.core)


## This creates a dataframe that can be plotted as a map
gsl.spdf.df_10.core <- abunPlotDF(gsl.df.Final10.core, pred.polys_200)
gsl.spdf.df_11.core <- abunPlotDF(gsl.df.Final11.core, pred.polys_200)
gsl.spdf.df_13.core <- abunPlotDF(gsl.df.Final13.core, pred.polys_200)
gsl.spdf.df_14.core <- abunPlotDF(gsl.df.Final14.core, pred.polys_200)
gsl.spdf.df_16.core <- abunPlotDF(gsl.df.Final16.core, pred.polys_200)
gsl.spdf.df_18.core <- abunPlotDF(gsl.df.Final18.core, pred.polys_200)
gsl.spdf.df_20.core <- abunPlotDF(gsl.df.Final20.core, pred.polys_200)

# save SPDFs
write.csv(gsl.spdf.df_10.core,file="Results/GSL/Plots/core_only/spdf/gsl.spdf.df_10.core.csv")
write.csv(gsl.spdf.df_11.core,file="Results/GSL/Plots/core_only/spdf/gsl.spdf.df_11.core.csv")
write.csv(gsl.spdf.df_13.core,file="Results/GSL/Plots/core_only/spdf/gsl.spdf.df_13.core.csv")
write.csv(gsl.spdf.df_14.core,file="Results/GSL/Plots/core_only/spdf/gsl.spdf.df_14.core.csv")
write.csv(gsl.spdf.df_16.core,file="Results/GSL/Plots/core_only/spdf/gsl.spdf.df_16.core.csv")
write.csv(gsl.spdf.df_18.core,file="Results/GSL/Plots/core_only/spdf/gsl.spdf.df_18.core.csv")
write.csv(gsl.spdf.df_20.core,file="Results/GSL/Plots/core_only/spdf/gsl.spdf.df_20.core.csv")



    ## Plotting continuous ####

# Load spatial dataframes
gsl.spdf.df_10.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_10.core.csv")
#gsl.spdf.df_11.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_11.core.csv") 
#gsl.spdf.df_13.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_13.core.csv") 
#gsl.spdf.df_14.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_14.core.csv") 
#gsl.spdf.df_16.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_16.core.csv") 
#gsl.spdf.df_18.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_18.core.csv") 
gsl.spdf.df_20.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_20.core.csv") 

# greyscale plots
GSL_plot_10_core_gr <- GSplotFun(gsl.spdf.df_10.core, survey.area.core, "abundance", "2010")
GSL_plot_11_core_gr <- GSplotFun(gsl.spdf.df_11.core, survey.area.core, "abundance")
GSL_plot_13_core_gr <- GSplotFun(gsl.spdf.df_13.core, survey.area.core, "abundance")
GSL_plot_14_core_gr <- GSplotFun(gsl.spdf.df_14.core, survey.area.core, "abundance")
GSL_plot_16_core_gr <- GSplotFun(gsl.spdf.df_16.core, survey.area.core, "abundance")
GSL_plot_18_core_gr <- GSplotFun(gsl.spdf.df_18.core, survey.area.core, "abundance")
GSL_plot_20_core_gr <- GSplotFun(gsl.spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(GSL_plot_10_core_gr,"Results/GSL/Plots/core_only/greyscale/GSL_plot_10_core_gr.png")
saveplot(GSL_plot_11_core_gr,"Results/GSL/Plots/core_only/greyscale/GSL_plot_11_core_gr.png")
saveplot(GSL_plot_13_core_gr,"Results/GSL/Plots/core_only/greyscale/GSL_plot_13_core_gr.png")
saveplot(GSL_plot_14_core_gr,"Results/GSL/Plots/core_only/greyscale/GSL_plot_14_core_gr.png")
saveplot(GSL_plot_16_core_gr,"Results/GSL/Plots/core_only/greyscale/GSL_plot_16_core_gr.png")
saveplot(GSL_plot_18_core_gr,"Results/GSL/Plots/core_only/greyscale/GSL_plot_18_core_gr.png")
saveplot(GSL_plot_20_core_gr,"Results/GSL/Plots/core_only/greyscale/GSL_plot_20_core_gr.png")

# colour plots
GSL_plot_10_core_col <- CLplotFun(gsl.spdf.df_10.core, survey.area.core, "abundance")
GSL_plot_11_core_col <- CLplotFun(gsl.spdf.df_11.core, survey.area.core, "abundance")
GSL_plot_13_core_col <- CLplotFun(gsl.spdf.df_13.core, survey.area.core, "abundance")
GSL_plot_14_core_col <- CLplotFun(gsl.spdf.df_14.core, survey.area.core, "abundance")
GSL_plot_16_core_col <- CLplotFun(gsl.spdf.df_16.core, survey.area.core, "abundance")
GSL_plot_18_core_col <- CLplotFun(gsl.spdf.df_18.core, survey.area.core, "abundance")
GSL_plot_20_core_col <- CLplotFun(gsl.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(GSL_plot_10_core_col,"Results/GSL/Plots/core_only/colour/GSL_plot_10_core_col.png")
saveplot(GSL_plot_11_core_col,"Results/GSL/Plots/core_only/colour/GSL_plot_11_core_col.png")
saveplot(GSL_plot_13_core_col,"Results/GSL/Plots/core_only/colour/GSL_plot_13_core_col.png")
saveplot(GSL_plot_14_core_col,"Results/GSL/Plots/core_only/colour/GSL_plot_14_core_col.png")
saveplot(GSL_plot_16_core_col,"Results/GSL/Plots/core_only/colour/GSL_plot_16_core_col.png")
saveplot(GSL_plot_18_core_col,"Results/GSL/Plots/core_only/colour/GSL_plot_18_core_col.png")
saveplot(GSL_plot_20_core_col,"Results/GSL/Plots/core_only/colour/GSL_plot_20_core_col.png")




## plot grids (abundance and variance)

# greyscale 
gsl_2yrs_gs <- 
  GSL_plot_10_core_gr + GSL_varplot_final10.core.bw  +
  GSL_plot_20_core_gr + GSL_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
gsl_2yrs_gs[[1]] <- gsl_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
gsl_2yrs_gs[[2]] <- gsl_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
gsl_2yrs_gs[[2]] <- gsl_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
gsl_2yrs_gs[[4]] <- gsl_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(gsl_2yrs_gs, "Results/GSL/Plots/core_only/greyscale/plot_grids/gsl_2yrs_gs.png")


## changing scale so that only 0 values are white, and all on-zero values are non-white

# first change all values in the data that are lower than the 1st quarter value (0.00008) to 0
gsl.spdf.df_10.core2 <- gsl.spdf.df_10.core
gsl.spdf.df_10.core2$abundance[gsl.spdf.df_10.core$abundance<1] <- NA

ggplot()+
  geom_polygon(aes_string(x="x",y="y",fill="abundance", group="group"), 
               data=gsl.spdf.df_10.core2)+
  gg.opts+
  geom_path(aes(x=x, y=y),data=survey.area.core)+
  labs(fill="abundance")+
  scale_fill_gradient(low = "grey", high = "black", na.value = "white")+
  coord_equal()




    ## Plotting discrtete ####
      # Add bins to SPDF - don't repeat ####

# this is the process of adding discrete bins to the abundance SPDFs. I will save them in a new folder so this only has to be done once

# Load original spatial dataframes
gsl.spdf.df_10.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_10.core.csv")
gsl.spdf.df_11.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_11.core.csv")
gsl.spdf.df_13.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_13.core.csv")
gsl.spdf.df_14.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_14.core.csv")
gsl.spdf.df_16.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_16.core.csv")
gsl.spdf.df_18.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_18.core.csv")
gsl.spdf.df_20.core <- read.csv("Results/GSL/Plots/core_only/spdf/gsl.spdf.df_20.core.csv")

# put spdf's into a list
dfs <- list(gsl.spdf.df_10.core,gsl.spdf.df_11.core,gsl.spdf.df_13.core,gsl.spdf.df_14.core,
            gsl.spdf.df_16.core,gsl.spdf.df_18.core,gsl.spdf.df_20.core)

# name the elements
names(dfs) <- c("gsl.spdf.df_10.core","gsl.spdf.df_11.core","gsl.spdf.df_13.core","gsl.spdf.df_14.core",
                "gsl.spdf.df_16.core","gsl.spdf.df_18.core","gsl.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunAbun)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(gsl.spdf.df_10.core,file="Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_10.core.csv")
write.csv(gsl.spdf.df_11.core,file="Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_11.core.csv")
write.csv(gsl.spdf.df_13.core,file="Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_13.core.csv")
write.csv(gsl.spdf.df_14.core,file="Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_14.core.csv")
write.csv(gsl.spdf.df_16.core,file="Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_16.core.csv")
write.csv(gsl.spdf.df_18.core,file="Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_18.core.csv")
write.csv(gsl.spdf.df_20.core,file="Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_20.core.csv")


      # Plotting ####

# Load spatial dataframes (with bins)
gsl.spdf.df_10.core <- read.csv("Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_10.core.csv")
gsl.spdf.df_11.core <- read.csv("Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_11.core.csv")
gsl.spdf.df_13.core <- read.csv("Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_13.core.csv")
gsl.spdf.df_14.core <- read.csv("Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_14.core.csv")
gsl.spdf.df_16.core <- read.csv("Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_16.core.csv")
gsl.spdf.df_18.core <- read.csv("Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_18.core.csv")
gsl.spdf.df_20.core <- read.csv("Results/GSL/Plots/core_only/spdf/bins/gsl.spdf.df_20.core.csv")

## plot greyscale
GSL_10_plot_bin_GS <- GSplotBin(gsl.spdf.df_10.core,"group2",survey.area.core,"Abundance","Relative abundance")
GSL_11_plot_bin_GS <- GSplotBin(gsl.spdf.df_11.core,"group2",survey.area.core,"Abundance","Relative abundance")
GSL_13_plot_bin_GS <- GSplotBin(gsl.spdf.df_13.core,"group2",survey.area.core,"Abundance","Relative abundance")
GSL_14_plot_bin_GS <- GSplotBin(gsl.spdf.df_14.core,"group2",survey.area.core,"Abundance","Relative abundance")
GSL_16_plot_bin_GS <- GSplotBin(gsl.spdf.df_16.core,"group2",survey.area.core,"Abundance","Relative abundance")
GSL_18_plot_bin_GS <- GSplotBin(gsl.spdf.df_18.core,"group2",survey.area.core,"Abundance","Relative abundance")
GSL_20_plot_bin_GS <- GSplotBin(gsl.spdf.df_20.core,"group2",survey.area.core,"Abundance","Relative abundance")


# save 
saveplot(GSL_10_plot_bin_GS,"Results/GSL/Plots/core_only/bins/GSL_10_plot_bin_GS.png")
saveplot(GSL_11_plot_bin_GS,"Results/GSL/Plots/core_only/bins/GSL_11_plot_bin_GS.png")
saveplot(GSL_13_plot_bin_GS,"Results/GSL/Plots/core_only/bins/GSL_13_plot_bin_GS.png")
saveplot(GSL_14_plot_bin_GS,"Results/GSL/Plots/core_only/bins/GSL_14_plot_bin_GS.png")
saveplot(GSL_16_plot_bin_GS,"Results/GSL/Plots/core_only/bins/GSL_16_plot_bin_GS.png")
saveplot(GSL_18_plot_bin_GS,"Results/GSL/Plots/core_only/bins/GSL_18_plot_bin_GS.png")
saveplot(GSL_20_plot_bin_GS,"Results/GSL/Plots/core_only/bins/GSL_20_plot_bin_GS.png")



## Variance estimation ####

# estimate variance
gsl.var.Final10.core <- varEstfun(preddata10_core, gslDSM.nb.5.xy4)
gsl.var.Final11.core <- varEstfun(preddata11_core, gslDSM.nb.5.xy4)
gsl.var.Final13.core <- varEstfun(preddata13_core, gslDSM.nb.5.xy4)
gsl.var.Final14.core <- varEstfun(preddata14_core, gslDSM.nb.5.xy4)
gsl.var.Final16.core <- varEstfun(preddata16_core, gslDSM.nb.5.xy4)
gsl.var.Final18.core <- varEstfun(preddata18_core, gslDSM.nb.5.xy4)
gsl.var.Final20.core <- varEstfun(preddata20_core, gslDSM.nb.5.xy4)

# save variance estimates
write.csv(gsl.var.Final10.core, file="Results/GSL/core_only/gsl.var10.core.csv")
write.csv(gsl.var.Final11.core, file="Results/GSL/core_only/gsl.var11.core.csv")
write.csv(gsl.var.Final13.core, file="Results/GSL/core_only/gsl.var13.core.csv")
write.csv(gsl.var.Final14.core, file="Results/GSL/core_only/gsl.var14.core.csv")
write.csv(gsl.var.Final16.core, file="Results/GSL/core_only/gsl.var16.core.csv")
write.csv(gsl.var.Final18.core, file="Results/GSL/core_only/gsl.var18.core.csv")
write.csv(gsl.var.Final20.core, file="Results/GSL/core_only/gsl.var20.core.csv")

# create spdf's for plotting
gsl.var.spdf.df_10.core <- varPlotDF(gsl.var.Final10.core, pred.polys_200)
gsl.var.spdf.df_11.core <- varPlotDF(gsl.var.Final11.core, pred.polys_200)
gsl.var.spdf.df_13.core <- varPlotDF(gsl.var.Final13.core, pred.polys_200)
gsl.var.spdf.df_14.core <- varPlotDF(gsl.var.Final14.core, pred.polys_200)
gsl.var.spdf.df_16.core <- varPlotDF(gsl.var.Final16.core, pred.polys_200)
gsl.var.spdf.df_18.core <- varPlotDF(gsl.var.Final18.core, pred.polys_200)
gsl.var.spdf.df_20.core <- varPlotDF(gsl.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(gsl.var.spdf.df_10.core,
          file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_10.core.csv")
write.csv(gsl.var.spdf.df_11.core,
          file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_11.core.csv")
write.csv(gsl.var.spdf.df_13.core,
          file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_13.core.csv")
write.csv(gsl.var.spdf.df_14.core,
          file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_14.core.csv")
write.csv(gsl.var.spdf.df_16.core,
          file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_16.core.csv")
write.csv(gsl.var.spdf.df_18.core,
          file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_18.core.csv")
write.csv(gsl.var.spdf.df_20.core,
          file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_20.core.csv")
   


    # Calculate CV & add bins to SPDF - don't repeat ####

# Load spatial dataframes
gsl.var.spdf.df_10.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_10.core.csv")
gsl.var.spdf.df_11.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_11.core.csv")
gsl.var.spdf.df_13.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_13.core.csv")
gsl.var.spdf.df_14.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_14.core.csv")
gsl.var.spdf.df_16.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_16.core.csv")
gsl.var.spdf.df_18.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_18.core.csv")
gsl.var.spdf.df_20.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_20.core.csv")


# first need to calculate CV from the variance (note: need the abundance spdf's loaded)
gsl.var.spdf.df_10.core <- CVaddFun(gsl.spdf.df_10.core,gsl.var.spdf.df_10.core)
gsl.var.spdf.df_11.core <- CVaddFun(gsl.spdf.df_11.core,gsl.var.spdf.df_11.core)
gsl.var.spdf.df_13.core <- CVaddFun(gsl.spdf.df_13.core,gsl.var.spdf.df_13.core)
gsl.var.spdf.df_14.core <- CVaddFun(gsl.spdf.df_14.core,gsl.var.spdf.df_14.core)
gsl.var.spdf.df_16.core <- CVaddFun(gsl.spdf.df_16.core,gsl.var.spdf.df_16.core)
gsl.var.spdf.df_18.core <- CVaddFun(gsl.spdf.df_18.core,gsl.var.spdf.df_18.core)
gsl.var.spdf.df_20.core <- CVaddFun(gsl.spdf.df_20.core,gsl.var.spdf.df_20.core)


# Save the SPDF's with the CV value but no bins (for continuous plotting of the CV)
write.csv(gsl.var.spdf.df_10.core,file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_10.core.csv")
write.csv(gsl.var.spdf.df_11.core,file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_11.core.csv")
write.csv(gsl.var.spdf.df_13.core,file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_13.core.csv")
write.csv(gsl.var.spdf.df_14.core,file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_14.core.csv")
write.csv(gsl.var.spdf.df_16.core,file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_16.core.csv")
write.csv(gsl.var.spdf.df_18.core,file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_18.core.csv")
write.csv(gsl.var.spdf.df_20.core,file="Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_20.core.csv")


### add bins

## Quartiles

# put spdf's into a list
dfs <- list(gsl.var.spdf.df_10.core,gsl.var.spdf.df_11.core,gsl.var.spdf.df_13.core,gsl.var.spdf.df_14.core,
            gsl.var.spdf.df_16.core,gsl.var.spdf.df_18.core,gsl.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("gsl.var.spdf.df_10.core","gsl.var.spdf.df_11.core","gsl.var.spdf.df_13.core",
                "gsl.var.spdf.df_14.core","gsl.var.spdf.df_16.core","gsl.var.spdf.df_18.core",
                "gsl.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(gsl.var.spdf.df_10.core,file="Results/GSL/Plots/variance/core_only/spdf/bins/gsl.var.spdf.df_10.core.csv")
write.csv(gsl.var.spdf.df_11.core,file="Results/GSL/Plots/variance/core_only/spdf/bins/gsl.var.spdf.df_11.core.csv")
write.csv(gsl.var.spdf.df_13.core,file="Results/GSL/Plots/variance/core_only/spdf/bins/gsl.var.spdf.df_13.core.csv")
write.csv(gsl.var.spdf.df_14.core,file="Results/GSL/Plots/variance/core_only/spdf/bins/gsl.var.spdf.df_14.core.csv")
write.csv(gsl.var.spdf.df_16.core,file="Results/GSL/Plots/variance/core_only/spdf/bins/gsl.var.spdf.df_16.core.csv")
write.csv(gsl.var.spdf.df_18.core,file="Results/GSL/Plots/variance/core_only/spdf/bins/gsl.var.spdf.df_18.core.csv")
write.csv(gsl.var.spdf.df_20.core,file="Results/GSL/Plots/variance/core_only/spdf/bins/gsl.var.spdf.df_20.core.csv")




## custom bins

# change column name from group2 to CV
gsl.var.spdf.df_10.core <- gsl.var.spdf.df_10.core %>% dplyr::rename(CV=group2)
gsl.var.spdf.df_11.core <- gsl.var.spdf.df_11.core %>% dplyr::rename(CV=group2)
gsl.var.spdf.df_13.core <- gsl.var.spdf.df_13.core %>% dplyr::rename(CV=group2)
gsl.var.spdf.df_14.core <- gsl.var.spdf.df_14.core %>% dplyr::rename(CV=group2)
gsl.var.spdf.df_16.core <- gsl.var.spdf.df_16.core %>% dplyr::rename(CV=group2)
gsl.var.spdf.df_18.core <- gsl.var.spdf.df_18.core %>% dplyr::rename(CV=group2)
gsl.var.spdf.df_20.core <- gsl.var.spdf.df_20.core %>% dplyr::rename(CV=group2)


# put spdf's into a list
dfs <- list(gsl.var.spdf.df_10.core,gsl.var.spdf.df_11.core,gsl.var.spdf.df_13.core,gsl.var.spdf.df_14.core,
            gsl.var.spdf.df_16.core,gsl.var.spdf.df_18.core,gsl.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("gsl.var.spdf.df_10.core","gsl.var.spdf.df_11.core","gsl.var.spdf.df_13.core",
                "gsl.var.spdf.df_14.core","gsl.var.spdf.df_16.core","gsl.var.spdf.df_18.core",
                "gsl.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar2)

# split elements into original dataframes
list2env(dfs, globalenv())


# save the SPDFs with the custom bins
write.csv(gsl.var.spdf.df_10.core,
          file="Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_10.core.csv")
write.csv(gsl.var.spdf.df_11.core,
          file="Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_11.core.csv")
write.csv(gsl.var.spdf.df_13.core,
          file="Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_13.core.csv")
write.csv(gsl.var.spdf.df_14.core,
          file="Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_14.core.csv")
write.csv(gsl.var.spdf.df_16.core,
          file="Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_16.core.csv")
write.csv(gsl.var.spdf.df_18.core,
          file="Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_18.core.csv")
write.csv(gsl.var.spdf.df_20.core,
          file="Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_20.core.csv")


    # Discrete bins ####

# load spdfs (which has already had CV calculated and then put into bins)
gsl.var.spdf.df_10.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_10.core.csv")
gsl.var.spdf.df_11.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_11.core.csv")
gsl.var.spdf.df_13.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_13.core.csv")
gsl.var.spdf.df_14.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_14.core.csv")
gsl.var.spdf.df_16.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_16.core.csv")
gsl.var.spdf.df_18.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_18.core.csv")
gsl.var.spdf.df_20.core <- read.csv("Results/GSL/Plots/variance/core_only/spdf/bins/custom/gsl.var.spdf.df_20.core.csv")

# make CV a factor and re-order. I have only done this for 2020 but you can copy the code for the other years if you need to
gsl.var.spdf.df_20.core$CV <- as.factor(gsl.var.spdf.df_20.core$CV)
gsl.var.spdf.df_20.core$CV <- factor(gsl.var.spdf.df_20.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))


# plot CV in bins
GSL_10_plot_bin_GS_var <- GSplotBin(gsl.var.spdf.df_10.core,"CV",survey.area.core,"Variance","CV")
GSL_11_plot_bin_GS_var <- GSplotBin(gsl.var.spdf.df_11.core,"CV",survey.area.core,"Variance","CV")
GSL_13_plot_bin_GS_var <- GSplotBin(gsl.var.spdf.df_13.core,"CV",survey.area.core,"Variance","CV")
GSL_14_plot_bin_GS_var <- GSplotBin(gsl.var.spdf.df_14.core,"CV",survey.area.core,"Variance","CV")
GSL_16_plot_bin_GS_var <- GSplotBin(gsl.var.spdf.df_16.core,"CV",survey.area.core,"Variance","CV")
GSL_18_plot_bin_GS_var <- GSplotBin(gsl.var.spdf.df_18.core,"CV",survey.area.core,"Variance","CV")
GSL_20_plot_bin_GS_var <- GSplotBin(gsl.var.spdf.df_20.core,"CV",survey.area.core,"Variance","CV")

# save 
saveplot(GSL_10_plot_bin_GS_var,"Results/GSL/Plots/variance/core_only/greyscale/bins/GSL_10_plot_bin_GS_var.png")
saveplot(GSL_11_plot_bin_GS_var,"Results/GSL/Plots/variance/core_only/greyscale/bins/GSL_11_plot_bin_GS_var.png")
saveplot(GSL_13_plot_bin_GS_var,"Results/GSL/Plots/variance/core_only/greyscale/bins/GSL_13_plot_bin_GS_var.png")
saveplot(GSL_14_plot_bin_GS_var,"Results/GSL/Plots/variance/core_only/greyscale/bins/GSL_14_plot_bin_GS_var.png")
saveplot(GSL_16_plot_bin_GS_var,"Results/GSL/Plots/variance/core_only/greyscale/bins/GSL_16_plot_bin_GS_var.png")
saveplot(GSL_18_plot_bin_GS_var,"Results/GSL/Plots/variance/core_only/greyscale/bins/GSL_18_plot_bin_GS_var.png")
saveplot(GSL_20_plot_bin_GS_var,"Results/GSL/Plots/variance/core_only/greyscale/bins/GSL_20_plot_bin_GS_var.png")



    # Plotting continous ####

# Load spatial dataframes
gsl.var.spdf.df_10.core<-read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_10.core.csv")
#gsl.var.spdf.df_11.core<-read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_11.core.csv")
#gsl.var.spdf.df_13.core<-read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_13.core.csv")
#gsl.var.spdf.df_14.core<-read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_14.core.csv")
#gsl.var.spdf.df_16.core<-read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_16.core.csv")
#gsl.var.spdf.df_18.core<-read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_18.core.csv")
gsl.var.spdf.df_20.core<-read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_20.core.csv")

# greyscale plots
GSL_varplot_final10.core.bw <- GSplotFun(gsl.var.spdf.df_10.core, survey.area.core, "variance", "2010")
GSL_varplot_final11.core.bw <- GSplotFun(gsl.var.spdf.df_11.core, survey.area.core, "variance")
GSL_varplot_final13.core.bw <- GSplotFun(gsl.var.spdf.df_13.core, survey.area.core, "variance")
GSL_varplot_final14.core.bw <- GSplotFun(gsl.var.spdf.df_14.core, survey.area.core, "variance")
GSL_varplot_final16.core.bw <- GSplotFun(gsl.var.spdf.df_16.core, survey.area.core, "variance")
GSL_varplot_final18.core.bw <- GSplotFun(gsl.var.spdf.df_18.core, survey.area.core, "variance")
GSL_varplot_final20.core.bw <- GSplotFun(gsl.var.spdf.df_20.core, survey.area.core, "variance", "2020")

# save greyscale
saveplot(GSL_varplot_final10.core.bw, "Results/GSL/Plots/variance/core_only/greyscale/2010_GSL_var.core.bw.png")
saveplot(GSL_varplot_final11.core.bw, "Results/GSL/Plots/variance/core_only/greyscale/2011_GSL_var.core.bw.png")
saveplot(GSL_varplot_final13.core.bw, "Results/GSL/Plots/variance/core_only/greyscale/2013_GSL_var.core.bw.png")
saveplot(GSL_varplot_final14.core.bw, "Results/GSL/Plots/variance/core_only/greyscale/2014_GSL_var.core.bw.png")
saveplot(GSL_varplot_final16.core.bw, "Results/GSL/Plots/variance/core_only/greyscale/2016_GSL_var.core.bw.png")
saveplot(GSL_varplot_final18.core.bw, "Results/GSL/Plots/variance/core_only/greyscale/2018_GSL_var.core.bw.png")
saveplot(GSL_varplot_final20.core.bw, "Results/GSL/Plots/variance/core_only/greyscale/2020_GSL_var.core.bw.png")

# colour plots
GSL_varplot_final10.core.col <- CLplotFun(gsl.var.spdf.df_10.core, survey.area.core, "variance")
GSL_varplot_final11.core.col <- CLplotFun(gsl.var.spdf.df_11.core, survey.area.core, "variance")
GSL_varplot_final13.core.col <- CLplotFun(gsl.var.spdf.df_13.core, survey.area.core, "variance")
GSL_varplot_final14.core.col <- CLplotFun(gsl.var.spdf.df_14.core, survey.area.core, "variance")
GSL_varplot_final16.core.col <- CLplotFun(gsl.var.spdf.df_16.core, survey.area.core, "variance")
GSL_varplot_final18.core.col <- CLplotFun(gsl.var.spdf.df_18.core, survey.area.core, "variance")
GSL_varplot_final20.core.col <- CLplotFun(gsl.var.spdf.df_20.core, survey.area.core, "variance")

# save colour
saveplot(GSL_varplot_final10.core.col, "Results/GSL/Plots/variance/core_only/colour/2010_GSL_var.core.col.png")
saveplot(GSL_varplot_final11.core.col, "Results/GSL/Plots/variance/core_only/colour/2011_GSL_var.core.col.png")
saveplot(GSL_varplot_final13.core.col, "Results/GSL/Plots/variance/core_only/colour/2013_GSL_var.core.col.png")
saveplot(GSL_varplot_final14.core.col, "Results/GSL/Plots/variance/core_only/colour/2014_GSL_var.core.col.png")
saveplot(GSL_varplot_final16.core.col, "Results/GSL/Plots/variance/core_only/colour/2016_GSL_var.core.col.png")
saveplot(GSL_varplot_final18.core.col, "Results/GSL/Plots/variance/core_only/colour/2018_GSL_var.core.col.png")
saveplot(GSL_varplot_final20.core.col, "Results/GSL/Plots/variance/core_only/colour/2020_GSL_var.core.col.png")



#### Wild pig ################################################################
## Load data ####

# Observation data. Unique to species 
pig_obsdata <- read.csv("Species_Data/PIG/R Data/obsdata.csv", header = TRUE)
pig_obsdata$object <- as.factor(pig_obsdata$object)
pig_obsdata$Sample.Label <- as.factor(pig_obsdata$Sample.Label)
str(pig_obsdata)
head(pig_obsdata)

# Transect data. Unique to species. CDS DF is Hr Cos so stratum not needed
pig_distdata <- read.csv("Species_Data/PIG/R Data/distdata.csv", header = TRUE)
pig_distdata$object <- as.factor(pig_distdata$object)
pig_distdata$NameObserver <- as.factor(pig_distdata$NameObserver)
pig_distdata$transect <- as.factor(pig_distdata$transect)
pig_distdata$year <- as.factor(pig_distdata$year)
pig_distdata$date <- as.Date(pig_distdata$date, format = "%d/%m/%Y")
str(pig_distdata)
head(pig_distdata)

# check there is no T20
pig_distdata[pig_distdata$transect=="20",]

## Plot the covariates across the grid with group sizes ####

# Warning - plots take a few minutes to run

# habitat
plot_PIGobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata200$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=pig_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/PIG/plot_PIGobs_habitat.png", plot_PIGobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_PIGobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys_200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=pig_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/PIG/plot_PIGobs_dstWater.png", plot_PIGobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_PIGobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=pig_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/PIG/plot_PIGobs_dstStlmnt.png", plot_PIGobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_PIGobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=pig_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/PIG/plot_PIGobs_dstRoad.png", plot_PIGobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_PIGobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=pig_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/PIG/plot_PIGobs_dstBorder.png", plot_PIGobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_PIGobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=pig_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/PIG/plot_PIGobs_dstStation.png", plot_PIGobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_PIGobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=pig_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/PIG/plot_PIGobs_dstELC.png", plot_PIGobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_PIGobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys_200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=pig_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/PIG/plot_PIGobs_elev.png", plot_PIGobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)

## Exploratory plots & linear models ####

## Checking that there are no large gaps in the range of variables for pig observations.  

# subset segdata to get only the segments with PIG observations
pig_varcheck <- segdata[match(pig_obsdata$Sample.Label,segdata$Sample.Label), ]

# habitat - fine
plot(segdata$habitat)
plot(pig_varcheck$habitat)

# dstWater - fine
hist(segdata$dstWater)
hist(pig_varcheck$dstWater)

# dstStlmnt - fine
hist(segdata$dstStlmnt)
hist(pig_varcheck$dstStlmnt)

# dstRoad - fine
hist(segdata$dstRoad)
hist(pig_varcheck$dstRoad)

# dstBorder - fine
hist(segdata$dstBorder)
hist(pig_varcheck$dstBorder)

# dstStation - fine
hist(segdata$dstStation)
hist(pig_varcheck$dstStation)

# elevation - small gap between 650-750 but shouldn't be a problem
hist(segdata$elevation)
hist(pig_varcheck$elevation)

## Histograms

# Distance 
pig_h1 <- ggplot(pig_distdata, aes(distance))+ geom_histogram(binwidth = 1)
pig_h2 <- ggplot(pig_distdata, aes(distance))+ geom_histogram(binwidth = 5)
pig_h3 <- ggplot(pig_distdata, aes(distance))+ geom_histogram(binwidth = 10)
pig_h4 <- ggplot(pig_distdata, aes(distance))+ geom_histogram(binwidth = 15)
pig_h5 <- ggplot(pig_distdata, aes(distance))+ geom_histogram(binwidth = 20)
pig_h6 <- ggplot(pig_distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(pig_h1,pig_h2,pig_h3,pig_h4,pig_h5,pig_h6)
# evidence of evasive movement. Spike at around 10m and then sharp drop at 25m. 


# cluster size, observer, habitat, year, month, transect
pig_h7 <- ggplot(pig_distdata, aes(size))+geom_histogram(binwidth = 0.5)
pig_h8 <- ggplot(pig_distdata, aes(NameObserver))+geom_histogram(stat="count")
pig_h9 <- ggplot(pig_distdata, aes(habitat))+geom_histogram(stat="count")
pig_h10 <- ggplot(pig_distdata, aes(year))+geom_histogram(stat="count")
pig_h11 <- ggplot(pig_distdata, aes(month))+geom_histogram(stat="count")
pig_h12 <- ggplot(pig_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(pig_h7,pig_h8,pig_h9,pig_h10,pig_h11,pig_h12)
# The vast majority of observations are of single animals, but there are a lot of observations where groups size > 1.  The species do tend to group, especially natal groups.  I think we need to estimate abudnace of groups, not individuals. The species generalism in terms of habitat prefrence is clear here.  Slightly more observations in dense forest but plenty in open forest too. Over the years, there has been a decrease in the number of obervations of this species which could reflect the hunting pressure.  

## Plots of distance against variables

plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

pig_d1 <- ggplot(pig_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
pig_d2 <- ggplot(pig_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
pig_d3 <- ggplot(pig_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
pig_d4 <- ggplot(pig_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
pig_d5 <- ggplot(pig_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(pig_d1,pig_d2,pig_d3,pig_d4,pig_d5)
# Distances are slightly greater in open habitat which makes sense. This is also reflected in the distances by transect - the transects in teh dense part of KSWS have smaller distances. 


## Plots of cluster size against variables
pig_s1 <- ggplot(pig_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
pig_s2 <- ggplot(pig_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
pig_s3 <- ggplot(pig_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
pig_s4 <- ggplot(pig_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
pig_s5 <- ggplot(pig_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(pig_s1,pig_s2,pig_s3,pig_s4,pig_s5)
# These plots don't reveal much

## Linear models

# group size ~ distance
newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=pig_distdata)
plot(pig_distdata$size~pig_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# No signficant relationship between group size and distance

## Estimating the detection function ####


### Taking DF mdoel from CDS so below model selection is obsolete. I'm leaving it in for reference.
# the PIG DF model is hr with no adjustment

pigDF.hr <- ds(pig_distdata, truncation=50, key="hr", adjustment=NULL)


    
## Fitting a spatial model ####

# The best detection function model for PIG is pigDF.hn 

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated.  

# I am not including dstELC in the models for PIG, as I don't think it is an appropriate variable, especially now that we are going to be predicting into the buffer zone

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
pig_obsdata <- pig_obsdata %>% filter(distance <= 50)

  ## Quasipoisson response ####

pigDSM.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+ 
                         s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                         habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.sat)
dev.off()
par(mfrow=c(2,3))
plot(pigDSM.sat, scale = 0)
gam.check(pigDSM.sat)
# DE = 5.5. dstWater & dstStation not sig. dstStlmt, elevation and dstBorder look overfitted

# reduce k for the above
pigDSM.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                         s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts")+  
                         habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.sat2)
dev.off()
par(mfrow=c(2,3))
plot(pigDSM.sat2, scale = 0)
gam.check(pigDSM.sat2)
# DE = 3. only dstBorder sig now. will try to remove a couple of terms

# remove dstWater
pigDSM.qp.3 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                          s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts")+  
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.qp.3)
dev.off()
par(mfrow=c(2,2))
plot(pigDSM.qp.3, scale = 0)
gam.check(pigDSM.qp.3)
# DE = 2.4. dstStation nearly sig but overfitted.

# remove elevation
pigDSM.qp.4 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=5)+ 
                          s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts")+  
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.qp.4)
dev.off()
par(mfrow=c(2,2))
plot(pigDSM.qp.4, scale = 0)
gam.check(pigDSM.qp.4)
# DE = 2.4. I don't really buy the dstStation relationship - looks like it's just fitting the points but that there is no real trend. Converely, I think dstStlmnt looks like there is an actual trend

# reduce k for dstStation and increase it for dstStlmnt
pigDSM.qp.5 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=6)+ 
                          s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+  
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.qp.5)
dev.off()
par(mfrow=c(2,2))
plot(pigDSM.qp.5, scale = 0)
gam.check(pigDSM.qp.5)
# DE = 1.58. dstStlmnt now sig and dstStation not. Very low DE

# incrase k for dstStation a tiny bit
pigDSM.qp.6 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=6)+ 
                          s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=7)+  
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.qp.6)
dev.off()
par(mfrow=c(2,2))
plot(pigDSM.qp.6, scale = 0)
gam.check(pigDSM.qp.6)
# DE = 1.27. Tried several different k settings - flogging dead horse with dstStation

# remove dstStation
pigDSM.qp.7 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=6)+ s(dstBorder,bs="ts",k=5)+   
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.qp.7)
par(mfrow=c(2,1))
plot(pigDSM.qp.7, scale = 0)
gam.check(pigDSM.qp.7)
# DE = 1.21. dstStlmnt is basically linear

# try dstStlmnt as linear term
pigDSM.qp.8 <- dsm(Nhat ~ s(dstBorder,bs="ts",k=5)+   
                          habitat + dstStlmnt,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.qp.8)
par(mfrow=c(2,1))
plot(pigDSM.qp.8, scale = 0)
gam.check(pigDSM.qp.8)
# DE = 1.2.

anova(pigDSM.qp.7,pigDSM.qp.8, test="Chisq")
# the models with dstStlmnt as linear and smooth terms are basically the same, except DE is tiny bit higher when it is a smooth term


### Not great models, but the best compromise is pigDSM.qp.7 


  ## Tweedie response ####

pigDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                            habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.tw.sat)
par(mfrow=c(2,3))
plot(pigDSM.tw.sat, scale = 0)
gam.check(pigDSM.tw.sat)
# AIC = 14063, DE = 0.56. No sig terms... gam.check suggests k is fine.

# increase k across the board
pigDSM.tw.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+ 
                             s(dstBorder,bs="ts",k=15)+ s(dstStation,bs="ts",k=15)+ 
                             s(elevation,bs="ts",k=15)+ habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.tw.sat2)
par(mfrow=c(2,3))
plot(pigDSM.tw.sat2, scale = 0)
gam.check(pigDSM.tw.sat2)
# AIC = 14063, DE = 0.56. nothing

# remove elevaiton
pigDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+ 
                             s(dstBorder,bs="ts",k=15)+ s(dstStation,bs="ts",k=15)+ 
                             habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.tw.3)
par(mfrow=c(2,2))
plot(pigDSM.tw.3, scale = 0)
gam.check(pigDSM.tw.3)
# AIC = 14063, DE = 0.56.

# remove dstStation
pigDSM.tw.4 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+ 
                             s(dstBorder,bs="ts",k=15)+ habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.tw.4)
par(mfrow=c(2,2))
plot(pigDSM.tw.4, scale = 0)
gam.check(pigDSM.tw.4)
# AIC = 14063, DE = 0.56.

# remove dstWater
pigDSM.tw.5 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=15)+ s(dstBorder,bs="ts",k=15)+ 
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.tw.5)
par(mfrow=c(2,1))
plot(pigDSM.tw.5, scale = 0)
gam.check(pigDSM.tw.5)
# AIC = 14063, DE = 0.56.

# remove dstStlmnt
pigDSM.tw.6 <- dsm(Nhat ~ s(dstBorder,bs="ts",k=15)+ 
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.tw.6)
par(mfrow=c(2,1))
plot(pigDSM.tw.6, scale = 0)
gam.check(pigDSM.tw.6)
# AIC = 14064, DE = 0.52.



### pigDSM.tw.6 is the best TW model, despite being shite

  ## Negative binomial response ####

pigDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+  
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+ s(elevation,bs="ts")+ 
                            habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.nb.sat)
par(mfrow=c(2,3))
plot(pigDSM.nb.sat, scale = 0)
gam.check(pigDSM.nb.sat)
# AIC = 2642, DE = 0.684. No terms are sig.  gam.check suggests k is fine

# increase k
pigDSM.nb.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+  
                             s(dstBorder,bs="ts",k=15)+ s(dstStation,bs="ts",k=15)+ 
                             s(elevation,bs="ts",k=15)+ 
                             habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.nb.sat2)
par(mfrow=c(2,3))
plot(pigDSM.nb.sat2, scale = 0)
gam.check(pigDSM.nb.sat2)
# AIC = 2642, DE = 0.684. nada

# remove dstStation
pigDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+  
                             s(dstBorder,bs="ts",k=15)+ s(elevation,bs="ts",k=15)+ 
                             habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.nb.3)
par(mfrow=c(2,3))
plot(pigDSM.nb.3, scale = 0)
gam.check(pigDSM.nb.3)
# AIC = 2642, DE = 0.685.

# remove elevation
pigDSM.nb.4 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+  
                          s(dstBorder,bs="ts",k=15)+
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.nb.4)
par(mfrow=c(2,3))
plot(pigDSM.nb.4, scale = 0)
gam.check(pigDSM.nb.4)
# AIC = 2642, DE = 0.685.

# remove dstWater
pigDSM.nb.5 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=15)+ s(dstBorder,bs="ts",k=15)+
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.nb.5)
par(mfrow=c(2,3))
plot(pigDSM.nb.5, scale = 0)
gam.check(pigDSM.nb.5)
# AIC = 2642, DE = 0.685.

# remove dstStlmnt
pigDSM.nb.6 <- dsm(Nhat ~ s(dstBorder,bs="ts",k=15)+
                          habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.nb.6)
par(mfrow=c(2,3))
plot(pigDSM.nb.6, scale = 0)
gam.check(pigDSM.nb.6)
# AIC = 2642, DE = 0.685.

### No NB models have any support


  ## GLM ####

# Poisson distribution as we're looking at abundance not density

# Saturated GLM
pigDSM.glm.1 <- dsm(Nhat ~ habitat + dstBorder + dstStlmnt + dstWater + dstStation + elevation,
                   pigDF.hr, segdata, pig_obsdata, method = "REML",
                   family = poisson(link = "log"), engine = "glm",
                   segment.area = segdata$Sample.Area, group=TRUE)
summary(pigDSM.glm.1)
# dstBorder is sig

# remove dstStation
pigDSM.glm.2 <- dsm(Nhat ~ habitat + dstBorder + dstStlmnt + dstWater + elevation,
                   pigDF.hr, segdata, pig_obsdata, method = "REML",
                   family = poisson(link = "log"), engine = "glm",
                   segment.area = segdata$Sample.Area, group=TRUE)
summary(pigDSM.glm.2)

# remove elevation
pigDSM.glm.3 <- dsm(Nhat ~ habitat + dstBorder + dstStlmnt + dstWater,
                   pigDF.hr, segdata, pig_obsdata, method = "REML",
                   family = poisson(link = "log"), engine = "glm",
                   segment.area = segdata$Sample.Area, group=TRUE)
summary(pigDSM.glm.3)

# remove dstWater
pigDSM.glm.4 <- dsm(Nhat ~ habitat + dstBorder + dstStlmnt,
                   pigDF.hr, segdata, pig_obsdata, method = "REML",
                   family = poisson(link = "log"), engine = "glm",
                   segment.area = segdata$Sample.Area, group=TRUE)
summary(pigDSM.glm.4)
# dstBorder sig, dstStlmnt p=0.08

# remove dstStlmnt
pigDSM.glm.5 <- dsm(Nhat ~ habitat + dstBorder,
                   pigDF.hr, segdata, pig_obsdata, method = "REML",
                   family = poisson(link = "log"), engine = "glm",
                   segment.area = segdata$Sample.Area, group=TRUE)
summary(pigDSM.glm.5)

# compare glm.4 and glm.5
anova(pigDSM.glm.4,pigDSM.glm.5, test="Chisq")
# close, but the model with dstStlmnt will be of more use

# pigDSM.glm.4


## Model selection ####

# pigDSM.qp.7 has dstBorder and dstStlmnt as sig terms, and pigDSM.tw.6 has dstBorder as the only (just) sig smooth. pigDSM.glm.4 has dstBorder and dstStlmnt as sig terms

summary(pigDSM.qp.7) # DE = 1.21
summary(pigDSM.tw.6) # DE = 0.52
summary(pigDSM.glm.4)

anova(pigDSM.qp.7,pigDSM.tw.6,test="Chisq")
# pigDSM.qp.7 is significantly better than tw.6

anova(pigDSM.qp.7,pigDSM.glm.4,test="Chisq")
# qp.7 has less residual deviance but not necessarily significant.

par(mfrow=c(2,2))
gam.check(pigDSM.qp.7)
plot(pigDSM.glm.4)
# both pretty hideous.

### pigDSM.qp.7 is the best PIG model

## Autocorrelation ####

par(mfrow=c(1,1))
dsm.cor(pigDSM.qp.7, max.lag=15, Segment.Label="Sample.Label")

# Some evidence of autocorrelation at lag 1, but nowhere else. 

# add univariate smooths
pigDSM.qp.7xy <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=6)+ s(dstBorder,bs="ts",k=5)+
                            s(x,bs="ts") + s(y,bs="ts")+
                            habitat,
                  pigDF.hr, segdata, pig_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(pigDSM.qp.7xy)
par(mfrow=c(2,2))
plot(pigDSM.qp.7xy, scale = 0)
gam.check(pigDSM.qp.7xy)
# DE = 4.39.

# pigDSM.qp.7xy final model


## Abundance estimation ####

# Predict over 2010 habitat 
pig.global.pred10.core <- predict(pigDSM.qp.7xy, preddata10_core, off.set = 40000)
write.csv(pig.global.pred10.core, file="Results/PIG/core_only/pig.pred10.core.csv")

pig.global.pred11.core <- predict(pigDSM.qp.7xy, preddata11_core, off.set = 40000)
write.csv(pig.global.pred11.core, file="Results/PIG/core_only/pig.pred11.core.csv")

pig.global.pred13.core <- predict(pigDSM.qp.7xy, preddata13_core, off.set = 40000)
write.csv(pig.global.pred13.core, file="Results/PIG/core_only/pig.pred13.core.csv")

pig.global.pred14.core <- predict(pigDSM.qp.7xy, preddata14_core, off.set = 40000)
write.csv(pig.global.pred14.core, file="Results/PIG/core_only/pig.pred14.core.csv")

pig.global.pred16.core <- predict(pigDSM.qp.7xy, preddata16_core, off.set = 40000)
write.csv(pig.global.pred16.core, file="Results/PIG/core_only/pig.pred16.core.csv")

pig.global.pred18.core <- predict(pigDSM.qp.7xy, preddata18_core, off.set = 40000)
write.csv(pig.global.pred18.core, file="Results/PIG/core_only/pig.pred18.core.csv")

pig.global.pred20.core <- predict(pigDSM.qp.7xy, preddata20_core, off.set = 40000)
write.csv(pig.global.pred20.core, file="Results/PIG/core_only/pig.pred20.core.csv")

# dataframes for plotting
pig.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = pig.global.pred10.core)

pig.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = pig.global.pred11.core)

pig.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = pig.global.pred13.core)

pig.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = pig.global.pred14.core)

pig.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = pig.global.pred16.core)

pig.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = pig.global.pred18.core)

pig.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = pig.global.pred20.core)


## This creates a dataframe that can be plotted as a map
pig.spdf.df_10.core <- abunPlotDF(pig.df.Final10.core, pred.polys_200)
pig.spdf.df_11.core <- abunPlotDF(pig.df.Final11.core, pred.polys_200)
pig.spdf.df_13.core <- abunPlotDF(pig.df.Final13.core, pred.polys_200)
pig.spdf.df_14.core <- abunPlotDF(pig.df.Final14.core, pred.polys_200)
pig.spdf.df_16.core <- abunPlotDF(pig.df.Final16.core, pred.polys_200)
pig.spdf.df_18.core <- abunPlotDF(pig.df.Final18.core, pred.polys_200)
pig.spdf.df_20.core <- abunPlotDF(pig.df.Final20.core, pred.polys_200)

# save SPDFs
write.csv(pig.spdf.df_10.core,file="Results/PIG/Plots/core_only/spdf/pig.spdf.df_10.core.csv")
write.csv(pig.spdf.df_11.core,file="Results/PIG/Plots/core_only/spdf/pig.spdf.df_11.core.csv")
write.csv(pig.spdf.df_13.core,file="Results/PIG/Plots/core_only/spdf/pig.spdf.df_13.core.csv")
write.csv(pig.spdf.df_14.core,file="Results/PIG/Plots/core_only/spdf/pig.spdf.df_14.core.csv")
write.csv(pig.spdf.df_16.core,file="Results/PIG/Plots/core_only/spdf/pig.spdf.df_16.core.csv")
write.csv(pig.spdf.df_18.core,file="Results/PIG/Plots/core_only/spdf/pig.spdf.df_18.core.csv")
write.csv(pig.spdf.df_20.core,file="Results/PIG/Plots/core_only/spdf/pig.spdf.df_20.core.csv")


    ## Plotting continuous ####

# Load spatial dataframes
pig.spdf.df_10.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_10.core.csv")
#pig.spdf.df_11.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_11.core.csv") 
#pig.spdf.df_13.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_13.core.csv") 
#pig.spdf.df_14.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_14.core.csv") 
#pig.spdf.df_16.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_16.core.csv") 
#pig.spdf.df_18.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_18.core.csv") 
pig.spdf.df_20.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_20.core.csv")

# greyscale plots
PIG_plot_10_core_gr <- GSplotFun(pig.spdf.df_10.core, survey.area.core, "abundance", "2010")
PIG_plot_11_core_gr <- GSplotFun(pig.spdf.df_11.core, survey.area.core, "abundance")
PIG_plot_13_core_gr <- GSplotFun(pig.spdf.df_13.core, survey.area.core, "abundance")
PIG_plot_14_core_gr <- GSplotFun(pig.spdf.df_14.core, survey.area.core, "abundance")
PIG_plot_16_core_gr <- GSplotFun(pig.spdf.df_16.core, survey.area.core, "abundance")
PIG_plot_18_core_gr <- GSplotFun(pig.spdf.df_18.core, survey.area.core, "abundance")
PIG_plot_20_core_gr <- GSplotFun(pig.spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(PIG_plot_10_core_gr,"Results/PIG/Plots/core_only/greyscale/PIG_plot_10_core_gr.png")
saveplot(PIG_plot_11_core_gr,"Results/PIG/Plots/core_only/greyscale/PIG_plot_11_core_gr.png")
saveplot(PIG_plot_13_core_gr,"Results/PIG/Plots/core_only/greyscale/PIG_plot_13_core_gr.png")
saveplot(PIG_plot_14_core_gr,"Results/PIG/Plots/core_only/greyscale/PIG_plot_14_core_gr.png")
saveplot(PIG_plot_16_core_gr,"Results/PIG/Plots/core_only/greyscale/PIG_plot_16_core_gr.png")
saveplot(PIG_plot_18_core_gr,"Results/PIG/Plots/core_only/greyscale/PIG_plot_18_core_gr.png")
saveplot(PIG_plot_20_core_gr,"Results/PIG/Plots/core_only/greyscale/PIG_plot_20_core_gr.png")

# colour plots
PIG_plot_10_core_col <- CLplotFun(pig.spdf.df_10.core, survey.area.core, "abundance")
PIG_plot_11_core_col <- CLplotFun(pig.spdf.df_11.core, survey.area.core, "abundance")
PIG_plot_13_core_col <- CLplotFun(pig.spdf.df_13.core, survey.area.core, "abundance")
PIG_plot_14_core_col <- CLplotFun(pig.spdf.df_14.core, survey.area.core, "abundance")
PIG_plot_16_core_col <- CLplotFun(pig.spdf.df_16.core, survey.area.core, "abundance")
PIG_plot_18_core_col <- CLplotFun(pig.spdf.df_18.core, survey.area.core, "abundance")
PIG_plot_20_core_col <- CLplotFun(pig.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(PIG_plot_10_core_col,"Results/PIG/Plots/core_only/colour/PIG_plot_10_core_col.png")
saveplot(PIG_plot_11_core_col,"Results/PIG/Plots/core_only/colour/PIG_plot_11_core_col.png")
saveplot(PIG_plot_13_core_col,"Results/PIG/Plots/core_only/colour/PIG_plot_13_core_col.png")
saveplot(PIG_plot_14_core_col,"Results/PIG/Plots/core_only/colour/PIG_plot_14_core_col.png")
saveplot(PIG_plot_16_core_col,"Results/PIG/Plots/core_only/colour/PIG_plot_16_core_col.png")
saveplot(PIG_plot_18_core_col,"Results/PIG/Plots/core_only/colour/PIG_plot_18_core_col.png")
saveplot(PIG_plot_20_core_col,"Results/PIG/Plots/core_only/colour/PIG_plot_20_core_col.png")



## plot grids (abundance and variance)

# greyscale 
pig_2yrs_gs <- 
  PIG_plot_10_core_gr + PIG_varplot_final10.core.bw  +
  PIG_plot_20_core_gr + PIG_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
pig_2yrs_gs[[1]] <- pig_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
pig_2yrs_gs[[2]] <- pig_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
pig_2yrs_gs[[2]] <- pig_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
pig_2yrs_gs[[4]] <- pig_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(pig_2yrs_gs, "Results/PIG/Plots/core_only/greyscale/plot_grids/pig_2yrs_gs.png")



    ## Plotting discrete ####
      # Add bins to SPDF - don't repeat ####

# this is the process of adding discrete bins to the abundance SPDFs. I will save them in a new folder so this only has to be done once

# Load original spatial dataframes
pig.spdf.df_10.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_10.core.csv")
pig.spdf.df_11.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_11.core.csv")
pig.spdf.df_13.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_13.core.csv")
pig.spdf.df_14.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_14.core.csv")
pig.spdf.df_16.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_16.core.csv")
pig.spdf.df_18.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_18.core.csv")
pig.spdf.df_20.core <- read.csv("Results/PIG/Plots/core_only/spdf/pig.spdf.df_20.core.csv")

# put spdf's into a list
dfs <- list(pig.spdf.df_10.core,pig.spdf.df_11.core,pig.spdf.df_13.core,pig.spdf.df_14.core,
            pig.spdf.df_16.core,pig.spdf.df_18.core,pig.spdf.df_20.core)

# name the elements
names(dfs) <- c("pig.spdf.df_10.core","pig.spdf.df_11.core","pig.spdf.df_13.core","pig.spdf.df_14.core",
                "pig.spdf.df_16.core","pig.spdf.df_18.core","pig.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunAbun)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(pig.spdf.df_10.core,file="Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_10.core.csv")
write.csv(pig.spdf.df_11.core,file="Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_11.core.csv")
write.csv(pig.spdf.df_13.core,file="Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_13.core.csv")
write.csv(pig.spdf.df_14.core,file="Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_14.core.csv")
write.csv(pig.spdf.df_16.core,file="Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_16.core.csv")
write.csv(pig.spdf.df_18.core,file="Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_18.core.csv")
write.csv(pig.spdf.df_20.core,file="Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_20.core.csv")


      # Plotting ####

# Load spatial dataframes (with bins)
pig.spdf.df_10.core <- read.csv("Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_10.core.csv")
pig.spdf.df_11.core <- read.csv("Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_11.core.csv")
pig.spdf.df_13.core <- read.csv("Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_13.core.csv")
pig.spdf.df_14.core <- read.csv("Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_14.core.csv")
pig.spdf.df_16.core <- read.csv("Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_16.core.csv")
pig.spdf.df_18.core <- read.csv("Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_18.core.csv")
pig.spdf.df_20.core <- read.csv("Results/PIG/Plots/core_only/spdf/bins/pig.spdf.df_20.core.csv")

## plot greyscale
PIG_10_plot_bin_GS <- GSplotBin(pig.spdf.df_10.core,"group2",survey.area.core,"Abundance","Relative abundance")
PIG_11_plot_bin_GS <- GSplotBin(pig.spdf.df_11.core,"group2",survey.area.core,"Abundance","Relative abundance")
PIG_13_plot_bin_GS <- GSplotBin(pig.spdf.df_13.core,"group2",survey.area.core,"Abundance","Relative abundance")
PIG_14_plot_bin_GS <- GSplotBin(pig.spdf.df_14.core,"group2",survey.area.core,"Abundance","Relative abundance")
PIG_16_plot_bin_GS <- GSplotBin(pig.spdf.df_16.core,"group2",survey.area.core,"Abundance","Relative abundance")
PIG_18_plot_bin_GS <- GSplotBin(pig.spdf.df_18.core,"group2",survey.area.core,"Abundance","Relative abundance")
PIG_20_plot_bin_GS <- GSplotBin(pig.spdf.df_20.core,"group2",survey.area.core,"Abundance","Relative abundance")


# save 
saveplot(PIG_10_plot_bin_GS,"Results/PIG/Plots/core_only/bins/PIG_10_plot_bin_GS.png")
saveplot(PIG_11_plot_bin_GS,"Results/PIG/Plots/core_only/bins/PIG_11_plot_bin_GS.png")
saveplot(PIG_13_plot_bin_GS,"Results/PIG/Plots/core_only/bins/PIG_13_plot_bin_GS.png")
saveplot(PIG_14_plot_bin_GS,"Results/PIG/Plots/core_only/bins/PIG_14_plot_bin_GS.png")
saveplot(PIG_16_plot_bin_GS,"Results/PIG/Plots/core_only/bins/PIG_16_plot_bin_GS.png")
saveplot(PIG_18_plot_bin_GS,"Results/PIG/Plots/core_only/bins/PIG_18_plot_bin_GS.png")
saveplot(PIG_20_plot_bin_GS,"Results/PIG/Plots/core_only/bins/PIG_20_plot_bin_GS.png")



## Variance estimation ####

# estimate variance
pig.var.Final10.core <- varEstfun(preddata10_core, pigDSM.qp.7xy)
pig.var.Final11.core <- varEstfun(preddata11_core, pigDSM.qp.7xy)
pig.var.Final13.core <- varEstfun(preddata13_core, pigDSM.qp.7xy)
pig.var.Final14.core <- varEstfun(preddata14_core, pigDSM.qp.7xy)
pig.var.Final16.core <- varEstfun(preddata16_core, pigDSM.qp.7xy)
pig.var.Final18.core <- varEstfun(preddata18_core, pigDSM.qp.7xy)
pig.var.Final20.core <- varEstfun(preddata20_core, pigDSM.qp.7xy)

# save variance estimates
write.csv(pig.var.Final10.core, file="Results/PIG/core_only/pig.var10.core.csv")
write.csv(pig.var.Final11.core, file="Results/PIG/core_only/pig.var11.core.csv")
write.csv(pig.var.Final13.core, file="Results/PIG/core_only/pig.var13.core.csv")
write.csv(pig.var.Final14.core, file="Results/PIG/core_only/pig.var14.core.csv")
write.csv(pig.var.Final16.core, file="Results/PIG/core_only/pig.var16.core.csv")
write.csv(pig.var.Final18.core, file="Results/PIG/core_only/pig.var18.core.csv")
write.csv(pig.var.Final20.core, file="Results/PIG/core_only/pig.var20.core.csv")

# create spdf's for plotting
pig.var.spdf.df_10.core <- varPlotDF(pig.var.Final10.core, pred.polys_200)
pig.var.spdf.df_11.core <- varPlotDF(pig.var.Final11.core, pred.polys_200)
pig.var.spdf.df_13.core <- varPlotDF(pig.var.Final13.core, pred.polys_200)
pig.var.spdf.df_14.core <- varPlotDF(pig.var.Final14.core, pred.polys_200)
pig.var.spdf.df_16.core <- varPlotDF(pig.var.Final16.core, pred.polys_200)
pig.var.spdf.df_18.core <- varPlotDF(pig.var.Final18.core, pred.polys_200)
pig.var.spdf.df_20.core <- varPlotDF(pig.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(pig.var.spdf.df_10.core,
          file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_10.core.csv")
write.csv(pig.var.spdf.df_11.core,
          file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_11.core.csv")
write.csv(pig.var.spdf.df_13.core,
          file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_13.core.csv")
write.csv(pig.var.spdf.df_14.core,
          file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_14.core.csv")
write.csv(pig.var.spdf.df_16.core,
          file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_16.core.csv")
write.csv(pig.var.spdf.df_18.core,
          file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_18.core.csv")
write.csv(pig.var.spdf.df_20.core,
          file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_20.core.csv") 



    # Calculate CV & add bins to SPDF - don't repeat ####

# Load spatial dataframes
pig.var.spdf.df_10.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_10.core.csv")
pig.var.spdf.df_11.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_11.core.csv")
pig.var.spdf.df_13.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_13.core.csv")
pig.var.spdf.df_14.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_14.core.csv")
pig.var.spdf.df_16.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_16.core.csv")
pig.var.spdf.df_18.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_18.core.csv")
pig.var.spdf.df_20.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_20.core.csv")


# first need to calculate CV from the variance (note: need the abundance spdf's loaded)
pig.var.spdf.df_10.core <- CVaddFun(pig.spdf.df_10.core,pig.var.spdf.df_10.core)
pig.var.spdf.df_11.core <- CVaddFun(pig.spdf.df_11.core,pig.var.spdf.df_11.core)
pig.var.spdf.df_13.core <- CVaddFun(pig.spdf.df_13.core,pig.var.spdf.df_13.core)
pig.var.spdf.df_14.core <- CVaddFun(pig.spdf.df_14.core,pig.var.spdf.df_14.core)
pig.var.spdf.df_16.core <- CVaddFun(pig.spdf.df_16.core,pig.var.spdf.df_16.core)
pig.var.spdf.df_18.core <- CVaddFun(pig.spdf.df_18.core,pig.var.spdf.df_18.core)
pig.var.spdf.df_20.core <- CVaddFun(pig.spdf.df_20.core,pig.var.spdf.df_20.core)


# Save the SPDF's with the CV value but no bins (for continuous plotting of the CV)
write.csv(pig.var.spdf.df_10.core,file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_10.core.csv")
write.csv(pig.var.spdf.df_11.core,file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_11.core.csv")
write.csv(pig.var.spdf.df_13.core,file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_13.core.csv")
write.csv(pig.var.spdf.df_14.core,file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_14.core.csv")
write.csv(pig.var.spdf.df_16.core,file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_16.core.csv")
write.csv(pig.var.spdf.df_18.core,file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_18.core.csv")
write.csv(pig.var.spdf.df_20.core,file="Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_20.core.csv")


### add bins

## Quartiles

# put spdf's into a list
dfs <- list(pig.var.spdf.df_10.core,pig.var.spdf.df_11.core,pig.var.spdf.df_13.core,pig.var.spdf.df_14.core,
            pig.var.spdf.df_16.core,pig.var.spdf.df_18.core,pig.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("pig.var.spdf.df_10.core","pig.var.spdf.df_11.core","pig.var.spdf.df_13.core",
                "pig.var.spdf.df_14.core","pig.var.spdf.df_16.core","pig.var.spdf.df_18.core",
                "pig.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(pig.var.spdf.df_10.core,file="Results/PIG/Plots/variance/core_only/spdf/bins/pig.var.spdf.df_10.core.csv")
write.csv(pig.var.spdf.df_11.core,file="Results/PIG/Plots/variance/core_only/spdf/bins/pig.var.spdf.df_11.core.csv")
write.csv(pig.var.spdf.df_13.core,file="Results/PIG/Plots/variance/core_only/spdf/bins/pig.var.spdf.df_13.core.csv")
write.csv(pig.var.spdf.df_14.core,file="Results/PIG/Plots/variance/core_only/spdf/bins/pig.var.spdf.df_14.core.csv")
write.csv(pig.var.spdf.df_16.core,file="Results/PIG/Plots/variance/core_only/spdf/bins/pig.var.spdf.df_16.core.csv")
write.csv(pig.var.spdf.df_18.core,file="Results/PIG/Plots/variance/core_only/spdf/bins/pig.var.spdf.df_18.core.csv")
write.csv(pig.var.spdf.df_20.core,file="Results/PIG/Plots/variance/core_only/spdf/bins/pig.var.spdf.df_20.core.csv")




## custom bins

# change column name from group2 to CV
pig.var.spdf.df_10.core <- pig.var.spdf.df_10.core %>% dplyr::rename(CV=group2)
pig.var.spdf.df_11.core <- pig.var.spdf.df_11.core %>% dplyr::rename(CV=group2)
pig.var.spdf.df_13.core <- pig.var.spdf.df_13.core %>% dplyr::rename(CV=group2)
pig.var.spdf.df_14.core <- pig.var.spdf.df_14.core %>% dplyr::rename(CV=group2)
pig.var.spdf.df_16.core <- pig.var.spdf.df_16.core %>% dplyr::rename(CV=group2)
pig.var.spdf.df_18.core <- pig.var.spdf.df_18.core %>% dplyr::rename(CV=group2)
pig.var.spdf.df_20.core <- pig.var.spdf.df_20.core %>% dplyr::rename(CV=group2)


# put spdf's into a list
dfs <- list(pig.var.spdf.df_10.core,pig.var.spdf.df_11.core,pig.var.spdf.df_13.core,pig.var.spdf.df_14.core,
            pig.var.spdf.df_16.core,pig.var.spdf.df_18.core,pig.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("pig.var.spdf.df_10.core","pig.var.spdf.df_11.core","pig.var.spdf.df_13.core",
                "pig.var.spdf.df_14.core","pig.var.spdf.df_16.core","pig.var.spdf.df_18.core",
                "pig.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar2)

# split elements into original dataframes
list2env(dfs, globalenv())


# save the SPDFs with the custom bins
write.csv(pig.var.spdf.df_10.core,
          file="Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_10.core.csv")
write.csv(pig.var.spdf.df_11.core,
          file="Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_11.core.csv")
write.csv(pig.var.spdf.df_13.core,
          file="Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_13.core.csv")
write.csv(pig.var.spdf.df_14.core,
          file="Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_14.core.csv")
write.csv(pig.var.spdf.df_16.core,
          file="Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_16.core.csv")
write.csv(pig.var.spdf.df_18.core,
          file="Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_18.core.csv")
write.csv(pig.var.spdf.df_20.core,
          file="Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_20.core.csv")


    # Discrete bins ####

# load spdfs (which has already had CV calculated and then put into bins)
pig.var.spdf.df_10.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_10.core.csv")
pig.var.spdf.df_11.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_11.core.csv")
pig.var.spdf.df_13.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_13.core.csv")
pig.var.spdf.df_14.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_14.core.csv")
pig.var.spdf.df_16.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_16.core.csv")
pig.var.spdf.df_18.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_18.core.csv")
pig.var.spdf.df_20.core <- read.csv("Results/PIG/Plots/variance/core_only/spdf/bins/custom/pig.var.spdf.df_20.core.csv")


# plot CV in bins
PIG_10_plot_bin_GS_var <- GSplotBin(pig.var.spdf.df_10.core,"CV",survey.area.core,"Variance","CV")
PIG_11_plot_bin_GS_var <- GSplotBin(pig.var.spdf.df_11.core,"CV",survey.area.core,"Variance","CV")
PIG_13_plot_bin_GS_var <- GSplotBin(pig.var.spdf.df_13.core,"CV",survey.area.core,"Variance","CV")
PIG_14_plot_bin_GS_var <- GSplotBin(pig.var.spdf.df_14.core,"CV",survey.area.core,"Variance","CV")
PIG_16_plot_bin_GS_var <- GSplotBin(pig.var.spdf.df_16.core,"CV",survey.area.core,"Variance","CV")
PIG_18_plot_bin_GS_var <- GSplotBin(pig.var.spdf.df_18.core,"CV",survey.area.core,"Variance","CV")
PIG_20_plot_bin_GS_var <- GSplotBin(pig.var.spdf.df_20.core,"CV",survey.area.core,"Variance","CV")

# save 
saveplot(PIG_10_plot_bin_GS_var,"Results/PIG/Plots/variance/core_only/greyscale/bins/PIG_10_plot_bin_GS_var.png")
saveplot(PIG_11_plot_bin_GS_var,"Results/PIG/Plots/variance/core_only/greyscale/bins/PIG_11_plot_bin_GS_var.png")
saveplot(PIG_13_plot_bin_GS_var,"Results/PIG/Plots/variance/core_only/greyscale/bins/PIG_13_plot_bin_GS_var.png")
saveplot(PIG_14_plot_bin_GS_var,"Results/PIG/Plots/variance/core_only/greyscale/bins/PIG_14_plot_bin_GS_var.png")
saveplot(PIG_16_plot_bin_GS_var,"Results/PIG/Plots/variance/core_only/greyscale/bins/PIG_16_plot_bin_GS_var.png")
saveplot(PIG_18_plot_bin_GS_var,"Results/PIG/Plots/variance/core_only/greyscale/bins/PIG_18_plot_bin_GS_var.png")
saveplot(PIG_20_plot_bin_GS_var,"Results/PIG/Plots/variance/core_only/greyscale/bins/PIG_20_plot_bin_GS_var.png")



    # Continous ####

# Load spatial dataframes
pig.var.spdf.df_10.core<-read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_10.core.csv")
#pig.var.spdf.df_11.core<-read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_11.core.csv")
#pig.var.spdf.df_13.core<-read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_13.core.csv")
#pig.var.spdf.df_14.core<-read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_14.core.csv")
#pig.var.spdf.df_16.core<-read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_16.core.csv")
#pig.var.spdf.df_18.core<-read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_18.core.csv")
pig.var.spdf.df_20.core<-read.csv("Results/PIG/Plots/variance/core_only/spdf/pig.var.spdf.df_20.core.csv")

# greyscale plots
PIG_varplot_final10.core.bw <- GSplotFun(pig.var.spdf.df_10.core, survey.area.core, "variance", "2010")
PIG_varplot_final11.core.bw <- GSplotFun(pig.var.spdf.df_11.core, survey.area.core, "variance")
PIG_varplot_final13.core.bw <- GSplotFun(pig.var.spdf.df_13.core, survey.area.core, "variance")
PIG_varplot_final14.core.bw <- GSplotFun(pig.var.spdf.df_14.core, survey.area.core, "variance")
PIG_varplot_final16.core.bw <- GSplotFun(pig.var.spdf.df_16.core, survey.area.core, "variance")
PIG_varplot_final18.core.bw <- GSplotFun(pig.var.spdf.df_18.core, survey.area.core, "variance")
PIG_varplot_final20.core.bw <- GSplotFun(pig.var.spdf.df_20.core, survey.area.core, "variance", "2020")

# save greyscale
saveplot(PIG_varplot_final10.core.bw, "Results/PIG/Plots/variance/core_only/greyscale/2010_PIG_var.core.bw.png")
saveplot(PIG_varplot_final11.core.bw, "Results/PIG/Plots/variance/core_only/greyscale/2011_PIG_var.core.bw.png")
saveplot(PIG_varplot_final13.core.bw, "Results/PIG/Plots/variance/core_only/greyscale/2013_PIG_var.core.bw.png")
saveplot(PIG_varplot_final14.core.bw, "Results/PIG/Plots/variance/core_only/greyscale/2014_PIG_var.core.bw.png")
saveplot(PIG_varplot_final16.core.bw, "Results/PIG/Plots/variance/core_only/greyscale/2016_PIG_var.core.bw.png")
saveplot(PIG_varplot_final18.core.bw, "Results/PIG/Plots/variance/core_only/greyscale/2018_PIG_var.core.bw.png")
saveplot(PIG_varplot_final20.core.bw, "Results/PIG/Plots/variance/core_only/greyscale/2020_PIG_var.core.bw.png")

# colour plots
PIG_varplot_final10.core.col <- CLplotFun(pig.var.spdf.df_10.core, survey.area.core, "variance")
PIG_varplot_final11.core.col <- CLplotFun(pig.var.spdf.df_11.core, survey.area.core, "variance")
PIG_varplot_final13.core.col <- CLplotFun(pig.var.spdf.df_13.core, survey.area.core, "variance")
PIG_varplot_final14.core.col <- CLplotFun(pig.var.spdf.df_14.core, survey.area.core, "variance")
PIG_varplot_final16.core.col <- CLplotFun(pig.var.spdf.df_16.core, survey.area.core, "variance")
PIG_varplot_final18.core.col <- CLplotFun(pig.var.spdf.df_18.core, survey.area.core, "variance")
PIG_varplot_final20.core.col <- CLplotFun(pig.var.spdf.df_20.core, survey.area.core, "variance")

# save colour
saveplot(PIG_varplot_final10.core.col, "Results/PIG/Plots/variance/core_only/colour/2010_PIG_var.core.col.png")
saveplot(PIG_varplot_final11.core.col, "Results/PIG/Plots/variance/core_only/colour/2011_PIG_var.core.col.png")
saveplot(PIG_varplot_final13.core.col, "Results/PIG/Plots/variance/core_only/colour/2013_PIG_var.core.col.png")
saveplot(PIG_varplot_final14.core.col, "Results/PIG/Plots/variance/core_only/colour/2014_PIG_var.core.col.png")
saveplot(PIG_varplot_final16.core.col, "Results/PIG/Plots/variance/core_only/colour/2016_PIG_var.core.col.png")
saveplot(PIG_varplot_final18.core.col, "Results/PIG/Plots/variance/core_only/colour/2018_PIG_var.core.col.png")
saveplot(PIG_varplot_final20.core.col, "Results/PIG/Plots/variance/core_only/colour/2020_PIG_var.core.col.png")


#### Gaur ####################################################################
## Load data ####

# Observation data. Unique to species 
gau_obsdata <- read.csv("Species_Data/GAU/R Data/obsdata.csv", header = TRUE)
gau_obsdata$object <- as.factor(gau_obsdata$object)
gau_obsdata$Sample.Label <- as.factor(gau_obsdata$Sample.Label)
str(gau_obsdata)
head(gau_obsdata)

# Transect data. Unique to species. DF model for gaur is Hr with no vars so stratum var is not needed
gau_distdata <- read.csv("Species_Data/GAU/R Data/distdata.csv", header = TRUE)
gau_distdata$object <- as.factor(gau_distdata$object)
gau_distdata$NameObserver <- as.factor(gau_distdata$NameObserver)
gau_distdata$trayonsect <- as.factor(gau_distdata$transect)
gau_distdata$year <- as.factor(gau_distdata$year)
gau_distdata$date <- as.Date(gau_distdata$date, format = "%d/%m/%Y")
str(gau_distdata)
head(gau_distdata)




## Plot the covariates across the grid with group sizes ####

# Warning - plots take a few minutes to run

# habitat
plot_GAUobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata200$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=gau_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/GAU/plot_GAUobs_habitat.png", plot_GAUobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_GAUobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys_200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=gau_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/GAU/plot_GAUobs_dstWater.png", plot_GAUobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_GAUobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=gau_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/GAU/plot_GAUobs_dstStlmnt.png", plot_GAUobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_GAUobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=gau_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/GAU/plot_GAUobs_dstRoad.png", plot_GAUobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_GAUobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=gau_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/GAU/plot_GAUobs_dstBorder.png", plot_GAUobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_GAUobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=gau_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/GAU/plot_GAUobs_dstStation.png", plot_GAUobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_GAUobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=gau_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/GAU/plot_GAUobs_dstELC.png", plot_GAUobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_GAUobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys_200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=gau_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/GAU/plot_GAUobs_elev.png", plot_GAUobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)


## Exploratory plots & linear models ####


## Checking that there are no large gaps in the range of variables for gaur observations. Only 34 observations so not that useful

# subset segdata to get only the segments with GAU observations
gau_varcheck <- segdata[match(gau_obsdata$Sample.Label,segdata$Sample.Label), ]

# habitat - fine
plot(segdata$habitat)
plot(gau_varcheck$habitat)

# dstWater - no observations beyond 1200 but variable range goes up to 2500
hist(segdata$dstWater)
hist(gau_varcheck$dstWater)

# dstStlmnt - actually ok - gap between 8km-10km but should be fine
hist(segdata$dstStlmnt)
hist(gau_varcheck$dstStlmnt)

# dstRoad - no observations beyond 1600 but variable range goes up to 4000
hist(segdata$dstRoad)
hist(gau_varcheck$dstRoad)

# dstBorder - no obs eyond 25000, variable goes up beyond 50000
hist(segdata$dstBorder)
hist(gau_varcheck$dstBorder)

# dstStation - no obs beyond 9000, variable goes up to 25000
hist(segdata$dstStation)
hist(gau_varcheck$dstStation)

# elevation - fine
hist(segdata$elevation)
hist(gau_varcheck$elevation)


## Histograms

# Distance 
gau_h1 <- ggplot(gau_distdata, aes(distance))+ geom_histogram(binwidth = 1)
gau_h2 <- ggplot(gau_distdata, aes(distance))+ geom_histogram(binwidth = 5)
gau_h3 <- ggplot(gau_distdata, aes(distance))+ geom_histogram(binwidth = 10)
gau_h4 <- ggplot(gau_distdata, aes(distance))+ geom_histogram(binwidth = 15)
gau_h5 <- ggplot(gau_distdata, aes(distance))+ geom_histogram(binwidth = 20)
gau_h6 <- ggplot(gau_distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(gau_h1,gau_h2,gau_h3,gau_h4,gau_h5,gau_h6)
# horrible distance histogram. frequency increases from the line up to a spike at around 25m before dropping off dramatically.  Too few data really.  I'm not sure we can attribute this to evasive modvement either because gaur are too large and noisy to "sneak" away from the line.


# cluster size, observer, habitat, year, month, transect
gau_h7 <- ggplot(gau_distdata, aes(size))+geom_histogram(binwidth = 0.5)
gau_h8 <- ggplot(gau_distdata, aes(NameObserver))+geom_histogram(stat="count")
gau_h9 <- ggplot(gau_distdata, aes(habitat))+geom_histogram(stat="count")
gau_h10 <- ggplot(gau_distdata, aes(year))+geom_histogram(stat="count")
gau_h11 <- ggplot(gau_distdata, aes(month))+geom_histogram(stat="count")
gau_h12 <- ggplot(gau_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(gau_h7,gau_h8,gau_h9,gau_h10,gau_h11,gau_h12)
# group size is tricky.  Most observations are of single individuals, although as a species they do aggregate into herds, although adult males are often solitary. Based on experience, I would be very doubtful that there are many proper herds left in Seima.  Probably safer to estimate abundance of groups.  Vadt majority of observations are in dense forest. Frequency of observations steadily decreasing over time, probably reflecting the species decline due to hunting. Transect 10 for some reason seems to get by far the most observations.  

## Plots of distance against variables

plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

gau_d1 <- ggplot(gau_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
gau_d2 <- ggplot(gau_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
gau_d3 <- ggplot(gau_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
gau_d4 <- ggplot(gau_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
gau_d5 <- ggplot(gau_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(gau_d1,gau_d2,gau_d3,gau_d4,gau_d5)
# No difference in distance for the different habitats. Not sure why distance seem to be larger in January but with so few data points I'd be wary of drawing any conclusions.

## Plots of cluster size against variables
gau_s1 <- ggplot(gau_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
gau_s2 <- ggplot(gau_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
gau_s3 <- ggplot(gau_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
gau_s4 <- ggplot(gau_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
gau_s5 <- ggplot(gau_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(gau_s1,gau_s2,gau_s3,gau_s4,gau_s5)
# Group size tends to be larger in open habitat. This could be a detectability issue. Largest group sizes all on transect 10.  Probably conincidence but I've asked Olly to double check the raw data to make sure there hasn't been a data entry issue.

## Linear models

# group size ~ distance
newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=gau_distdata)
plot(gau_distdata$size~gau_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# No significant relationship between group size and distance.  There is one outlier (distance = 100m, grp size = 1) which is pulling the trend line down, but which will be excluded anyway.

# remove the outlier
gau_distdata <- gau_distdata %>% filter(distance < 99)

newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=gau_distdata)
plot(gau_distdata$size~gau_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# no relationship between group size and distance

## Estimating the detection function ####

# DF model from CDS is hn with no adjustment

gauDF.hn <- ds(gau_distdata, truncation=50, key="hn")


    
## Fitting a spatial model ####

# The best detection function model for GAU is gauDF.hr 

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated.  

# I am not including dstELC in the models for GAU, as I don't think it is an appropriate variable, especially now that we are going to be predicting into the buffer zone

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
gau_obsdata <- gau_obsdata %>% filter(distance <= 50)



  ## Quasipoisson response ####

gauDSM.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+
                         s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                         habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.sat)
dev.off()
par(mfrow=c(2,3))
plot(gauDSM.sat, scale = 0)
gam.check(gauDSM.sat)
# DE = 39.2. only dstStation sig.  dstStation overfitted (but not sig)

# remove dstStlmnt
gauDSM.qp.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts")+
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.qp.2)
dev.off()
par(mfrow=c(2,2))
plot(gauDSM.qp.2, scale = 0)
gam.check(gauDSM.qp.2)
# DE = 35.4. no terms sig

# remove elevation from sat instead of dstStlmnt
gauDSM.qp.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.qp.3)
dev.off()
par(mfrow=c(2,3))
plot(gauDSM.qp.3, scale = 0)
gam.check(gauDSM.qp.3)
# DE = 31.1. NO terms sig

# compare qp.2 and qp.3
anova(gauDSM.qp.2,gauDSM.qp.3,test="Chisq")
# qp.2 is the better model

# as qp.2 but remove dstWater
gauDSM.qp.4 <- dsm(Nhat ~ s(elevation,bs="ts")+
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.qp.4)
dev.off()
par(mfrow=c(2,2))
plot(gauDSM.qp.4, scale = 0)
gam.check(gauDSM.qp.4)
# DE = 31.5. no terms sig

# remove elevation
gauDSM.qp.5 <- dsm(Nhat ~ s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.qp.5)
dev.off()
par(mfrow=c(2,1))
plot(gauDSM.qp.5, scale = 0)
gam.check(gauDSM.qp.5)
# DE = 27.5. dstBorder and dstStation sig. Smooth plots suggest they do very little! I will try and increase k for them both - gam.check suggests k is too low

# increase k
gauDSM.qp.6 <- dsm(Nhat ~ s(dstBorder,bs="ts",k=11)+ s(dstStation,bs="ts",k=11)+  
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.qp.6)
dev.off()
par(mfrow=c(2,1))
plot(gauDSM.qp.6, scale = 0)
gam.check(gauDSM.qp.6)
# DE = 27.5. DE the same. No full convergence.  


### gauDSM.qp.5 is the best QP model

  ## Tweedie response ####

gauDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                            habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.tw.sat)
dev.off()
par(mfrow=c(2,3))
plot(gauDSM.tw.sat, scale = 0)
gam.check(gauDSM.tw.sat)
# AIC = 12365, DE = 14.9. dstWater (and maybe dstBorder) sig. elevation and dstStation look shrunk

# remove elevation
gauDSM.tw.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+  
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                            habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.tw.2)
dev.off()
par(mfrow=c(2,2))
plot(gauDSM.tw.2, scale = 0)
gam.check(gauDSM.tw.2)
# AIC = 12363, DE = 10.8. no terms sig.

# remove dstStation
gauDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+  
                          s(dstBorder,bs="ts")+   
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.tw.3)
dev.off()
par(mfrow=c(2,2))
plot(gauDSM.tw.3, scale = 0)
gam.check(gauDSM.tw.3)
# AIC = 12363, DE = 10.8. Same as above. dstStlmnt has the lowest EDF and highest p value so will be removed

# remove dstStlmnt
gauDSM.tw.4 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstBorder,bs="ts")+   
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.tw.4)
dev.off()
par(mfrow=c(2,1))
plot(gauDSM.tw.4, scale = 0)
gam.check(gauDSM.tw.4)
# AIC = 12362, DE = 10.5. Better model. dstBorder now sig. Both remaining terms potentially linear

# first try and remove dstWater
gauDSM.tw.5 <- dsm(Nhat ~ s(dstBorder,bs="ts")+ habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.tw.5)
dev.off()
par(mfrow=c(2,1))
plot(gauDSM.tw.5, scale = 0)
gam.check(gauDSM.tw.5)
# AIC = 12361, DE = 9.7. Better model

# try and increase k for dstBorder
gauDSM.tw.6 <- dsm(Nhat ~ s(dstBorder,bs="ts",k=10)+ habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.tw.6)
dev.off()
par(mfrow=c(1,1))
plot(gauDSM.tw.6, scale = 0)
gam.check(gauDSM.tw.6)
# AIC = 12361, DE = 9.7. Doesn't make any difference.

# As tw.4 but with dstWater as a linear term
gauDSM.tw.7 <- dsm(Nhat ~ s(dstBorder,bs="ts")+   
                          habitat + dstWater,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.tw.7)
dev.off()
par(mfrow=c(1,1))
plot(gauDSM.tw.7, scale = 0)
gam.check(gauDSM.tw.7)
# AIC = 12363, DE = 10.5. No better


### best TW model is gauDSM.tw.5, although not great.  It wil be worth trying a GLM with habitat and dstBorder as linear terms

  ## Negative binomial response ####

# saturated model
gauDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+  
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                            habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.nb.sat)
dev.off()
par(mfrow=c(2,3))
plot(gauDSM.nb.sat, scale = 0)
gam.check(gauDSM.nb.sat)
# AIC = 472, DE = 41.8. Much better AIC and DE than QP and TW models. dstSlmnt & dstBorder sig. dstStlmnt looks linear

# remove dstStation
gauDSM.nb.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+  
                          s(dstBorder,bs="ts")+ s(elevation,bs="ts")+ 
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.nb.2)
dev.off()
par(mfrow=c(2,3))
plot(gauDSM.nb.2, scale = 0)
gam.check(gauDSM.nb.2)
# AIC = 482, DE = 37.5. Worse model than above.

# as sat but remove elevation instead
gauDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+   
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.nb.3)
dev.off()
par(mfrow=c(2,3))
plot(gauDSM.nb.3, scale = 0)
gam.check(gauDSM.nb.3)
# AIC = 492, DE = 29.9. even worse!

# as sat but remove dstWater
gauDSM.nb.4 <- dsm(Nhat ~ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+  
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.nb.4)
dev.off()
par(mfrow=c(2,2))
plot(gauDSM.nb.4, scale = 0)
gam.check(gauDSM.nb.4)
# AIC = 494, DE = 27.

# as sat but with dstWater, elevation, and dstStation removed
gauDSM.nb.5 <- dsm(Nhat ~  s(dstStlmnt,bs="ts")+ s(dstBorder,bs="ts")+   
                           habitat,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.nb.5)
dev.off()
par(mfrow=c(2,1))
plot(gauDSM.nb.5, scale = 0)
gam.check(gauDSM.nb.5)
# AIC = 472, DE = 41.8. Better! same AIC and DE as the saturated model. Both terms sig. dstStlmnt look linear

# dstStlmnt as linear term
gauDSM.nb.6 <- dsm(Nhat ~  s(dstBorder,bs="ts")+ habitat + dstStlmnt,
                  gauDF.hn, segdata, gau_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(gauDSM.nb.6)
dev.off()
par(mfrow=c(1,1))
plot(gauDSM.nb.6, scale = 0, xlab="Distance to border (m)")
gam.check(gauDSM.nb.6)
# AIC = 472, DE = 42.

# check if removing dstBorder makes a difference to predictions (due to gaps in variable range)
gauDSM.glm.1 <- dsm(Nhat ~ habitat + dstStlmnt,
                   gauDF.hn, segdata, gau_obsdata, method = "REML",
                   family = poisson(link = "log"), engine = "glm",
                   segment.area = segdata$Sample.Area, group=TRUE)
summary(gauDSM.glm.1)
# The dstBorder variable is important in restricting the presence of the species to areas close to the border, which we know is true. If we remove dstBorder the model predicts the species elsewhere in KSWS where we know it does not exist.  I am not sure exactly why this species is restricted to those areas, but dstBorder is clearly acting as a proxy for whatever it is. Therefore leaving dstBorder in is important, despite the gap in the variable range.

### gauDSM.nb.6 is the best NB model


## Model selection ####

# QP model is gauDSM.qp.5
# TW model is gauDSM.tw.5
# NB model is gauDSM.nb.6

summary(gauDSM.qp.5)
summary(gauDSM.tw.5)
summary(gauDSM.nb.6)

# NB model has a much lower AIC than TW, adnd the higest DE

# Compare QP model with NB model
anova(gauDSM.qp.5,gauDSM.nb.6,test="Chisq")
# NB model has much lower residual deviance

# Check Q-Q plot
dev.off()
par(mfrow=c(2,2))
gam.check(gauDSM.qp.5)
gam.check(gauDSM.nb.6)
# Q-Q plot much better for NB model

### Best model is gauDSM.nb.6

## Autocorrelation ####

par(mfrow=c(1,1))
dsm.cor(gauDSM.nb.6, max.lag=15, Segment.Label="Sample.Label")

# No significant autocorrelation

## Abundance estimation ####

# Predict over 2010 habitat 
gau.global.pred10.core <- predict(gauDSM.nb.6, preddata10_core, off.set = 40000)
write.csv(gau.global.pred10.core, file="Results/GAU/core_only/gau.pred10.core.csv")

# Predict over 2011 habitat 
gau.global.pred11.core <- predict(gauDSM.nb.6, preddata11_core, off.set = 40000)
write.csv(gau.global.pred11.core, file="Results/GAU/core_only/gau.pred11.core.csv")

# Predict over 2013 habitat 
gau.global.pred13.core <- predict(gauDSM.nb.6, preddata13_core, off.set = 40000)
write.csv(gau.global.pred13.core, file="Results/GAU/core_only/gau.pred13.core.csv")

# Predict over 2014 habitat 
gau.global.pred14.core <- predict(gauDSM.nb.6, preddata14_core, off.set = 40000)
write.csv(gau.global.pred14.core, file="Results/GAU/core_only/gau.pred14.core.csv")

# Predict over 2016 habitat 
gau.global.pred16.core <- predict(gauDSM.nb.6, preddata16_core, off.set = 40000)
write.csv(gau.global.pred16.core, file="Results/GAU/core_only/gau.pred16.core.csv")

# Predict over 2018 habitat 
gau.global.pred18.core <- predict(gauDSM.nb.6, preddata18_core, off.set = 40000)
write.csv(gau.global.pred18.core, file="Results/GAU/core_only/gau.pred18.core.csv")

# Predict over 2018 habitat 
gau.global.pred20.core <- predict(gauDSM.nb.6, preddata20_core, off.set = 40000)
write.csv(gau.global.pred20.core, file="Results/GAU/core_only/gau.pred20.core.csv")



# Create new dataframes for plotting
gau.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = gau.global.pred10.core)

gau.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = gau.global.pred11.core)

gau.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = gau.global.pred13.core)

gau.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = gau.global.pred14.core)

gau.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = gau.global.pred16.core)

gau.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = gau.global.pred18.core)

gau.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = gau.global.pred20.core)


## This creates a dataframe that can be plotted as a map
gau.spdf.df_10.core <- abunPlotDF(gau.df.Final10.core, pred.polys_200)
gau.spdf.df_11.core <- abunPlotDF(gau.df.Final11.core, pred.polys_200)
gau.spdf.df_13.core <- abunPlotDF(gau.df.Final13.core, pred.polys_200)
gau.spdf.df_14.core <- abunPlotDF(gau.df.Final14.core, pred.polys_200)
gau.spdf.df_16.core <- abunPlotDF(gau.df.Final16.core, pred.polys_200)
gau.spdf.df_18.core <- abunPlotDF(gau.df.Final18.core, pred.polys_200)
gau.spdf.df_20.core <- abunPlotDF(gau.df.Final20.core, pred.polys_200)

# save SPDFs
write.csv(gau.spdf.df_10.core,file="Results/GAU/Plots/core_only/spdf/gau.spdf.df_10.core.csv")
write.csv(gau.spdf.df_11.core,file="Results/GAU/Plots/core_only/spdf/gau.spdf.df_11.core.csv")
write.csv(gau.spdf.df_13.core,file="Results/GAU/Plots/core_only/spdf/gau.spdf.df_13.core.csv")
write.csv(gau.spdf.df_14.core,file="Results/GAU/Plots/core_only/spdf/gau.spdf.df_14.core.csv")
write.csv(gau.spdf.df_16.core,file="Results/GAU/Plots/core_only/spdf/gau.spdf.df_16.core.csv")
write.csv(gau.spdf.df_18.core,file="Results/GAU/Plots/core_only/spdf/gau.spdf.df_18.core.csv")
write.csv(gau.spdf.df_20.core,file="Results/GAU/Plots/core_only/spdf/gau.spdf.df_20.core.csv")


    ## Plotting continuous ####

# Load spatial dataframes
gau.spdf.df_10.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_10.core.csv")
#gau.spdf.df_11.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_11.core.csv") 
#gau.spdf.df_13.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_13.core.csv") 
#gau.spdf.df_14.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_14.core.csv") 
#gau.spdf.df_16.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_16.core.csv") 
#gau.spdf.df_18.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_18.core.csv")
gau.spdf.df_20.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_20.core.csv")

# greyscale plots
GAU_plot_10_core_gr <- GSplotFun(gau.spdf.df_10.core, survey.area.core, "abundance", "2010")
GAU_plot_11_core_gr <- GSplotFun(gau.spdf.df_11.core, survey.area.core, "abundance")
GAU_plot_13_core_gr <- GSplotFun(gau.spdf.df_13.core, survey.area.core, "abundance")
GAU_plot_14_core_gr <- GSplotFun(gau.spdf.df_14.core, survey.area.core, "abundance")
GAU_plot_16_core_gr <- GSplotFun(gau.spdf.df_16.core, survey.area.core, "abundance")
GAU_plot_18_core_gr <- GSplotFun(gau.spdf.df_18.core, survey.area.core, "abundance")
GAU_plot_20_core_gr <- GSplotFun(gau.spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(GAU_plot_10_core_gr,"Results/GAU/Plots/core_only/greyscale/GAU_plot_10_core_gr.png")
saveplot(GAU_plot_11_core_gr,"Results/GAU/Plots/core_only/greyscale/GAU_plot_11_core_gr.png")
saveplot(GAU_plot_13_core_gr,"Results/GAU/Plots/core_only/greyscale/GAU_plot_13_core_gr.png")
saveplot(GAU_plot_14_core_gr,"Results/GAU/Plots/core_only/greyscale/GAU_plot_14_core_gr.png")
saveplot(GAU_plot_16_core_gr,"Results/GAU/Plots/core_only/greyscale/GAU_plot_16_core_gr.png")
saveplot(GAU_plot_18_core_gr,"Results/GAU/Plots/core_only/greyscale/GAU_plot_18_core_gr.png")
saveplot(GAU_plot_20_core_gr,"Results/GAU/Plots/core_only/greyscale/GAU_plot_20_core_gr.png")

# colour plots
GAU_plot_10_core_col <- CLplotFun(gau.spdf.df_10.core, survey.area.core, "abundance")
GAU_plot_11_core_col <- CLplotFun(gau.spdf.df_11.core, survey.area.core, "abundance")
GAU_plot_13_core_col <- CLplotFun(gau.spdf.df_13.core, survey.area.core, "abundance")
GAU_plot_14_core_col <- CLplotFun(gau.spdf.df_14.core, survey.area.core, "abundance")
GAU_plot_16_core_col <- CLplotFun(gau.spdf.df_16.core, survey.area.core, "abundance")
GAU_plot_18_core_col <- CLplotFun(gau.spdf.df_18.core, survey.area.core, "abundance")
GAU_plot_20_core_col <- CLplotFun(gau.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(GAU_plot_10_core_col,"Results/GAU/Plots/core_only/colour/GAU_plot_10_core_col.png")
saveplot(GAU_plot_11_core_col,"Results/GAU/Plots/core_only/colour/GAU_plot_11_core_col.png")
saveplot(GAU_plot_13_core_col,"Results/GAU/Plots/core_only/colour/GAU_plot_13_core_col.png")
saveplot(GAU_plot_14_core_col,"Results/GAU/Plots/core_only/colour/GAU_plot_14_core_col.png")
saveplot(GAU_plot_16_core_col,"Results/GAU/Plots/core_only/colour/GAU_plot_16_core_col.png")
saveplot(GAU_plot_18_core_col,"Results/GAU/Plots/core_only/colour/GAU_plot_18_core_col.png")
saveplot(GAU_plot_20_core_col,"Results/GAU/Plots/core_only/colour/GAU_plot_20_core_col.png")


## plot grids (abundance and variance)

# greyscale 
gau_2yrs_gs <- 
  GAU_plot_10_core_gr + GAU_varplot_final10.core.bw  +
  GAU_plot_20_core_gr + GAU_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
gau_2yrs_gs[[1]] <- gau_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
gau_2yrs_gs[[2]] <- gau_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
gau_2yrs_gs[[2]] <- gau_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
gau_2yrs_gs[[4]] <- gau_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(gau_2yrs_gs, "Results/GAU/Plots/core_only/greyscale/plot_grids/gau_2yrs_gs.png")




    ## Plotting discrete ####
      # Add bins to SPDF - don't repeat ####

# this is the process of adding discrete bins to the abundance SPDFs. I will save them in a new folder so this only has to be done once

# Load original spatial dataframes
gau.spdf.df_10.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_10.core.csv")
gau.spdf.df_11.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_11.core.csv")
gau.spdf.df_13.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_13.core.csv")
gau.spdf.df_14.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_14.core.csv")
gau.spdf.df_16.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_16.core.csv")
gau.spdf.df_18.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_18.core.csv")
gau.spdf.df_20.core <- read.csv("Results/GAU/Plots/core_only/spdf/gau.spdf.df_20.core.csv")

# put spdf's into a list
dfs <- list(gau.spdf.df_10.core,gau.spdf.df_11.core,gau.spdf.df_13.core,gau.spdf.df_14.core,
            gau.spdf.df_16.core,gau.spdf.df_18.core,gau.spdf.df_20.core)

# name the elements
names(dfs) <- c("gau.spdf.df_10.core","gau.spdf.df_11.core","gau.spdf.df_13.core","gau.spdf.df_14.core",
                "gau.spdf.df_16.core","gau.spdf.df_18.core","gau.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunAbun)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(gau.spdf.df_10.core,file="Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_10.core.csv")
write.csv(gau.spdf.df_11.core,file="Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_11.core.csv")
write.csv(gau.spdf.df_13.core,file="Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_13.core.csv")
write.csv(gau.spdf.df_14.core,file="Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_14.core.csv")
write.csv(gau.spdf.df_16.core,file="Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_16.core.csv")
write.csv(gau.spdf.df_18.core,file="Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_18.core.csv")
write.csv(gau.spdf.df_20.core,file="Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_20.core.csv")


      # Plotting ####

# Load spatial dataframes (with bins)
gau.spdf.df_10.core <- read.csv("Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_10.core.csv")
gau.spdf.df_11.core <- read.csv("Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_11.core.csv")
gau.spdf.df_13.core <- read.csv("Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_13.core.csv")
gau.spdf.df_14.core <- read.csv("Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_14.core.csv")
gau.spdf.df_16.core <- read.csv("Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_16.core.csv")
gau.spdf.df_18.core <- read.csv("Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_18.core.csv")
gau.spdf.df_20.core <- read.csv("Results/GAU/Plots/core_only/spdf/bins/gau.spdf.df_20.core.csv")

## plot greyscale
GAU_10_plot_bin_GS <- GSplotBin(gau.spdf.df_10.core,"group2",survey.area.core,"Abundance","Relative abundance")
GAU_11_plot_bin_GS <- GSplotBin(gau.spdf.df_11.core,"group2",survey.area.core,"Abundance","Relative abundance")
GAU_13_plot_bin_GS <- GSplotBin(gau.spdf.df_13.core,"group2",survey.area.core,"Abundance","Relative abundance")
GAU_14_plot_bin_GS <- GSplotBin(gau.spdf.df_14.core,"group2",survey.area.core,"Abundance","Relative abundance")
GAU_16_plot_bin_GS <- GSplotBin(gau.spdf.df_16.core,"group2",survey.area.core,"Abundance","Relative abundance")
GAU_18_plot_bin_GS <- GSplotBin(gau.spdf.df_18.core,"group2",survey.area.core,"Abundance","Relative abundance")
GAU_20_plot_bin_GS <- GSplotBin(gau.spdf.df_20.core,"group2",survey.area.core,"Abundance","Relative abundance")


# save 
saveplot(GAU_10_plot_bin_GS,"Results/GAU/Plots/core_only/bins/GAU_10_plot_bin_GS.png")
saveplot(GAU_11_plot_bin_GS,"Results/GAU/Plots/core_only/bins/GAU_11_plot_bin_GS.png")
saveplot(GAU_13_plot_bin_GS,"Results/GAU/Plots/core_only/bins/GAU_13_plot_bin_GS.png")
saveplot(GAU_14_plot_bin_GS,"Results/GAU/Plots/core_only/bins/GAU_14_plot_bin_GS.png")
saveplot(GAU_16_plot_bin_GS,"Results/GAU/Plots/core_only/bins/GAU_16_plot_bin_GS.png")
saveplot(GAU_18_plot_bin_GS,"Results/GAU/Plots/core_only/bins/GAU_18_plot_bin_GS.png")
saveplot(GAU_20_plot_bin_GS,"Results/GAU/Plots/core_only/bins/GAU_20_plot_bin_GS.png")



## Variance estimation ####
  
# estimate variance
gau.var.Final10.core <- varEstfun(preddata10_core, gauDSM.nb.6)
gau.var.Final11.core <- varEstfun(preddata11_core, gauDSM.nb.6)
gau.var.Final13.core <- varEstfun(preddata13_core, gauDSM.nb.6)
gau.var.Final14.core <- varEstfun(preddata14_core, gauDSM.nb.6)
gau.var.Final16.core <- varEstfun(preddata16_core, gauDSM.nb.6)
gau.var.Final18.core <- varEstfun(preddata18_core, gauDSM.nb.6)
gau.var.Final20.core <- varEstfun(preddata20_core, gauDSM.nb.6)

# save variance estimates
write.csv(gau.var.Final10.core, file="Results/GAU/core_only/gau.var10.core.csv")
write.csv(gau.var.Final11.core, file="Results/GAU/core_only/gau.var11.core.csv")
write.csv(gau.var.Final13.core, file="Results/GAU/core_only/gau.var13.core.csv")
write.csv(gau.var.Final14.core, file="Results/GAU/core_only/gau.var14.core.csv")
write.csv(gau.var.Final16.core, file="Results/GAU/core_only/gau.var16.core.csv")
write.csv(gau.var.Final18.core, file="Results/GAU/core_only/gau.var18.core.csv")
write.csv(gau.var.Final20.core, file="Results/GAU/core_only/gau.var20.core.csv")

# create spdf's for plotting
gau.var.spdf.df_10.core <- varPlotDF(gau.var.Final10.core, pred.polys_200)
gau.var.spdf.df_11.core <- varPlotDF(gau.var.Final11.core, pred.polys_200)
gau.var.spdf.df_13.core <- varPlotDF(gau.var.Final13.core, pred.polys_200)
gau.var.spdf.df_14.core <- varPlotDF(gau.var.Final14.core, pred.polys_200)
gau.var.spdf.df_16.core <- varPlotDF(gau.var.Final16.core, pred.polys_200)
gau.var.spdf.df_18.core <- varPlotDF(gau.var.Final18.core, pred.polys_200)
gau.var.spdf.df_20.core <- varPlotDF(gau.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(gau.var.spdf.df_10.core,
          file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_10.core.csv")
write.csv(gau.var.spdf.df_11.core,
          file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_11.core.csv")
write.csv(gau.var.spdf.df_13.core,
          file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_13.core.csv")
write.csv(gau.var.spdf.df_14.core,
          file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_14.core.csv")
write.csv(gau.var.spdf.df_16.core,
          file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_16.core.csv")
write.csv(gau.var.spdf.df_18.core,
          file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_18.core.csv")
write.csv(gau.var.spdf.df_20.core,
          file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_20.core.csv") 



    # Calculate CV & add bins to SPDF - don't repeat ####

# Load spatial dataframes
gau.var.spdf.df_10.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_10.core.csv")
gau.var.spdf.df_11.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_11.core.csv")
gau.var.spdf.df_13.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_13.core.csv")
gau.var.spdf.df_14.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_14.core.csv")
gau.var.spdf.df_16.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_16.core.csv")
gau.var.spdf.df_18.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_18.core.csv")
gau.var.spdf.df_20.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_20.core.csv")


# first need to calculate CV from the variance (note: need the abundance spdf's loaded)
gau.var.spdf.df_10.core <- CVaddFun(gau.spdf.df_10.core,gau.var.spdf.df_10.core)
gau.var.spdf.df_11.core <- CVaddFun(gau.spdf.df_11.core,gau.var.spdf.df_11.core)
gau.var.spdf.df_13.core <- CVaddFun(gau.spdf.df_13.core,gau.var.spdf.df_13.core)
gau.var.spdf.df_14.core <- CVaddFun(gau.spdf.df_14.core,gau.var.spdf.df_14.core)
gau.var.spdf.df_16.core <- CVaddFun(gau.spdf.df_16.core,gau.var.spdf.df_16.core)
gau.var.spdf.df_18.core <- CVaddFun(gau.spdf.df_18.core,gau.var.spdf.df_18.core)
gau.var.spdf.df_20.core <- CVaddFun(gau.spdf.df_20.core,gau.var.spdf.df_20.core)


# Save the SPDF's with the CV value but no bins (for continuous plotting of the CV)
write.csv(gau.var.spdf.df_10.core,file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_10.core.csv")
write.csv(gau.var.spdf.df_11.core,file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_11.core.csv")
write.csv(gau.var.spdf.df_13.core,file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_13.core.csv")
write.csv(gau.var.spdf.df_14.core,file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_14.core.csv")
write.csv(gau.var.spdf.df_16.core,file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_16.core.csv")
write.csv(gau.var.spdf.df_18.core,file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_18.core.csv")
write.csv(gau.var.spdf.df_20.core,file="Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_20.core.csv")


### add bins

## Quartiles

# put spdf's into a list
dfs <- list(gau.var.spdf.df_10.core,gau.var.spdf.df_11.core,gau.var.spdf.df_13.core,gau.var.spdf.df_14.core,
            gau.var.spdf.df_16.core,gau.var.spdf.df_18.core,gau.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("gau.var.spdf.df_10.core","gau.var.spdf.df_11.core","gau.var.spdf.df_13.core",
                "gau.var.spdf.df_14.core","gau.var.spdf.df_16.core","gau.var.spdf.df_18.core",
                "gau.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(gau.var.spdf.df_10.core,file="Results/GAU/Plots/variance/core_only/spdf/bins/gau.var.spdf.df_10.core.csv")
write.csv(gau.var.spdf.df_11.core,file="Results/GAU/Plots/variance/core_only/spdf/bins/gau.var.spdf.df_11.core.csv")
write.csv(gau.var.spdf.df_13.core,file="Results/GAU/Plots/variance/core_only/spdf/bins/gau.var.spdf.df_13.core.csv")
write.csv(gau.var.spdf.df_14.core,file="Results/GAU/Plots/variance/core_only/spdf/bins/gau.var.spdf.df_14.core.csv")
write.csv(gau.var.spdf.df_16.core,file="Results/GAU/Plots/variance/core_only/spdf/bins/gau.var.spdf.df_16.core.csv")
write.csv(gau.var.spdf.df_18.core,file="Results/GAU/Plots/variance/core_only/spdf/bins/gau.var.spdf.df_18.core.csv")
write.csv(gau.var.spdf.df_20.core,file="Results/GAU/Plots/variance/core_only/spdf/bins/gau.var.spdf.df_20.core.csv")




## custom bins

# change column name from group2 to CV
gau.var.spdf.df_10.core <- gau.var.spdf.df_10.core %>% dplyr::rename(CV=group2)
gau.var.spdf.df_11.core <- gau.var.spdf.df_11.core %>% dplyr::rename(CV=group2)
gau.var.spdf.df_13.core <- gau.var.spdf.df_13.core %>% dplyr::rename(CV=group2)
gau.var.spdf.df_14.core <- gau.var.spdf.df_14.core %>% dplyr::rename(CV=group2)
gau.var.spdf.df_16.core <- gau.var.spdf.df_16.core %>% dplyr::rename(CV=group2)
gau.var.spdf.df_18.core <- gau.var.spdf.df_18.core %>% dplyr::rename(CV=group2)
gau.var.spdf.df_20.core <- gau.var.spdf.df_20.core %>% dplyr::rename(CV=group2)


# put spdf's into a list
dfs <- list(gau.var.spdf.df_10.core,gau.var.spdf.df_11.core,gau.var.spdf.df_13.core,gau.var.spdf.df_14.core,
            gau.var.spdf.df_16.core,gau.var.spdf.df_18.core,gau.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("gau.var.spdf.df_10.core","gau.var.spdf.df_11.core","gau.var.spdf.df_13.core",
                "gau.var.spdf.df_14.core","gau.var.spdf.df_16.core","gau.var.spdf.df_18.core",
                "gau.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar2)

# split elements into original dataframes
list2env(dfs, globalenv())


# save the SPDFs with the custom bins
write.csv(gau.var.spdf.df_10.core,
          file="Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_10.core.csv")
write.csv(gau.var.spdf.df_11.core,
          file="Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_11.core.csv")
write.csv(gau.var.spdf.df_13.core,
          file="Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_13.core.csv")
write.csv(gau.var.spdf.df_14.core,
          file="Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_14.core.csv")
write.csv(gau.var.spdf.df_16.core,
          file="Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_16.core.csv")
write.csv(gau.var.spdf.df_18.core,
          file="Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_18.core.csv")
write.csv(gau.var.spdf.df_20.core,
          file="Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_20.core.csv")


    # Discrete bins ####

# load spdfs (which has already had CV calculated and then put into bins)
gau.var.spdf.df_10.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_10.core.csv")
gau.var.spdf.df_11.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_11.core.csv")
gau.var.spdf.df_13.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_13.core.csv")
gau.var.spdf.df_14.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_14.core.csv")
gau.var.spdf.df_16.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_16.core.csv")
gau.var.spdf.df_18.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_18.core.csv")
gau.var.spdf.df_20.core <- read.csv("Results/GAU/Plots/variance/core_only/spdf/bins/custom/gau.var.spdf.df_20.core.csv")


# plot CV in bins
GAU_10_plot_bin_GS_var <- GSplotBin(gau.var.spdf.df_10.core,"CV",survey.area.core,"Variance","CV")
GAU_11_plot_bin_GS_var <- GSplotBin(gau.var.spdf.df_11.core,"CV",survey.area.core,"Variance","CV")
GAU_13_plot_bin_GS_var <- GSplotBin(gau.var.spdf.df_13.core,"CV",survey.area.core,"Variance","CV")
GAU_14_plot_bin_GS_var <- GSplotBin(gau.var.spdf.df_14.core,"CV",survey.area.core,"Variance","CV")
GAU_16_plot_bin_GS_var <- GSplotBin(gau.var.spdf.df_16.core,"CV",survey.area.core,"Variance","CV")
GAU_18_plot_bin_GS_var <- GSplotBin(gau.var.spdf.df_18.core,"CV",survey.area.core,"Variance","CV")
GAU_20_plot_bin_GS_var <- GSplotBin(gau.var.spdf.df_20.core,"CV",survey.area.core,"Variance","CV")

# save 
saveplot(GAU_10_plot_bin_GS_var,"Results/GAU/Plots/variance/core_only/greyscale/bins/GAU_10_plot_bin_GS_var.png")
saveplot(GAU_11_plot_bin_GS_var,"Results/GAU/Plots/variance/core_only/greyscale/bins/GAU_11_plot_bin_GS_var.png")
saveplot(GAU_13_plot_bin_GS_var,"Results/GAU/Plots/variance/core_only/greyscale/bins/GAU_13_plot_bin_GS_var.png")
saveplot(GAU_14_plot_bin_GS_var,"Results/GAU/Plots/variance/core_only/greyscale/bins/GAU_14_plot_bin_GS_var.png")
saveplot(GAU_16_plot_bin_GS_var,"Results/GAU/Plots/variance/core_only/greyscale/bins/GAU_16_plot_bin_GS_var.png")
saveplot(GAU_18_plot_bin_GS_var,"Results/GAU/Plots/variance/core_only/greyscale/bins/GAU_18_plot_bin_GS_var.png")
saveplot(GAU_20_plot_bin_GS_var,"Results/GAU/Plots/variance/core_only/greyscale/bins/GAU_20_plot_bin_GS_var.png")



    # Continuous ####

# Load spatial dataframes
gau.var.spdf.df_10.core<-read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_10.core.csv")
#gau.var.spdf.df_11.core<-read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_11.core.csv")
#gau.var.spdf.df_13.core<-read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_13.core.csv")
#gau.var.spdf.df_14.core<-read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_14.core.csv")
#gau.var.spdf.df_16.core<-read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_16.core.csv")
#gau.var.spdf.df_18.core<-read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_18.core.csv")
gau.var.spdf.df_20.core<-read.csv("Results/GAU/Plots/variance/core_only/spdf/gau.var.spdf.df_20.core.csv")

# greyscale plots
GAU_varplot_final10.core.bw <- GSplotFun(gau.var.spdf.df_10.core, survey.area.core, "variance", "2010")
GAU_varplot_final11.core.bw <- GSplotFun(gau.var.spdf.df_11.core, survey.area.core, "variance")
GAU_varplot_final13.core.bw <- GSplotFun(gau.var.spdf.df_13.core, survey.area.core, "variance")
GAU_varplot_final14.core.bw <- GSplotFun(gau.var.spdf.df_14.core, survey.area.core, "variance")
GAU_varplot_final16.core.bw <- GSplotFun(gau.var.spdf.df_16.core, survey.area.core, "variance")
GAU_varplot_final18.core.bw <- GSplotFun(gau.var.spdf.df_18.core, survey.area.core, "variance")
GAU_varplot_final20.core.bw <- GSplotFun(gau.var.spdf.df_20.core, survey.area.core, "variance", "2020")

# save greyscale
saveplot(GAU_varplot_final10.core.bw, "Results/GAU/Plots/variance/core_only/greyscale/2010_GAU_var.core.bw.png")
saveplot(GAU_varplot_final11.core.bw, "Results/GAU/Plots/variance/core_only/greyscale/2011_GAU_var.core.bw.png")
saveplot(GAU_varplot_final13.core.bw, "Results/GAU/Plots/variance/core_only/greyscale/2013_GAU_var.core.bw.png")
saveplot(GAU_varplot_final14.core.bw, "Results/GAU/Plots/variance/core_only/greyscale/2014_GAU_var.core.bw.png")
saveplot(GAU_varplot_final16.core.bw, "Results/GAU/Plots/variance/core_only/greyscale/2016_GAU_var.core.bw.png")
saveplot(GAU_varplot_final18.core.bw, "Results/GAU/Plots/variance/core_only/greyscale/2018_GAU_var.core.bw.png")
saveplot(GAU_varplot_final20.core.bw, "Results/GAU/Plots/variance/core_only/greyscale/2020_GAU_var.core.bw.png")

# colour plots
GAU_varplot_final10.core.col <- CLplotFun(gau.var.spdf.df_10.core, survey.area.core, "variance")
GAU_varplot_final11.core.col <- CLplotFun(gau.var.spdf.df_11.core, survey.area.core, "variance")
GAU_varplot_final13.core.col <- CLplotFun(gau.var.spdf.df_13.core, survey.area.core, "variance")
GAU_varplot_final14.core.col <- CLplotFun(gau.var.spdf.df_14.core, survey.area.core, "variance")
GAU_varplot_final16.core.col <- CLplotFun(gau.var.spdf.df_16.core, survey.area.core, "variance")
GAU_varplot_final18.core.col <- CLplotFun(gau.var.spdf.df_18.core, survey.area.core, "variance")
GAU_varplot_final20.core.col <- CLplotFun(gau.var.spdf.df_20.core, survey.area.core, "variance")

# save colour
saveplot(GAU_varplot_final10.core.col, "Results/GAU/Plots/variance/core_only/colour/2010_GAU_var.core.col.png")
saveplot(GAU_varplot_final11.core.col, "Results/GAU/Plots/variance/core_only/colour/2011_GAU_var.core.col.png")
saveplot(GAU_varplot_final13.core.col, "Results/GAU/Plots/variance/core_only/colour/2013_GAU_var.core.col.png")
saveplot(GAU_varplot_final14.core.col, "Results/GAU/Plots/variance/core_only/colour/2014_GAU_var.core.col.png")
saveplot(GAU_varplot_final16.core.col, "Results/GAU/Plots/variance/core_only/colour/2016_GAU_var.core.col.png")
saveplot(GAU_varplot_final18.core.col, "Results/GAU/Plots/variance/core_only/colour/2018_GAU_var.core.col.png")
saveplot(GAU_varplot_final20.core.col, "Results/GAU/Plots/variance/core_only/colour/2020_GAU_var.core.col.png")


#### Long-tailed macaque #####################################################
## Load data ####

# Observation data. Unique to species 
ltm_obsdata <- read.csv("Species_Data/LTM/R Data/obsdata.csv", header = TRUE)
ltm_obsdata$object <- as.factor(ltm_obsdata$object)
ltm_obsdata$Sample.Label <- as.factor(ltm_obsdata$Sample.Label)
str(ltm_obsdata)
head(ltm_obsdata)

# Transect data. Unique to species. Need stratum var for DF
ltm_distdata <- read.csv("Species_Data/LTM/R Data/distdata.csv", header = TRUE)
ltm_distdata$object <- as.factor(ltm_distdata$object)
ltm_distdata$NameObserver <- as.factor(ltm_distdata$NameObserver)
ltm_distdata$transect <- as.factor(ltm_distdata$transect)
ltm_distdata$stratum <- as.vector(scale(ltm_distdata$year, center=T, scale=T))
ltm_distdata$year <- as.factor(ltm_distdata$year)
ltm_distdata$date <- as.Date(ltm_distdata$date, format = "%d/%m/%Y")
str(ltm_distdata)
head(ltm_distdata)

## Plot the covariates across the grid with group sizes ####

## Warning - plots take a while to run

# habitat
plot_LTMobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata200$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=ltm_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/LTM/plot_LTMobs_habitat.png", plot_LTMobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_LTMobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys_200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=ltm_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/LTM/plot_LTMobs_dstWater.png", plot_LTMobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_LTMobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=ltm_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/LTM/plot_LTMobs_dstStlmnt.png", plot_LTMobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_LTMobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=ltm_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/LTM/plot_LTMobs_dstRoad.png", plot_LTMobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_LTMobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=ltm_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/LTM/plot_LTMobs_dstBorder.png", plot_LTMobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_LTMobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=ltm_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/LTM/plot_LTMobs_dstStation.png", plot_LTMobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_LTMobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=ltm_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/LTM/plot_LTMobs_dstELC.png", plot_LTMobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_LTMobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys_200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=ltm_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/LTM/plot_LTMobs_elev.png", plot_LTMobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)


## Exploratory plots and linear models ####


## Checking that there are no large gaps in the range of variables for LTM observations. 

# subset segdata to get only the segments with LTM observations
ltm_varcheck <- segdata[match(ltm_obsdata$Sample.Label,segdata$Sample.Label), ]

# habitat 
plot(segdata$habitat)
plot(ltm_varcheck$habitat)
# Good coverage

# dstWater 
hist(segdata$dstWater)
hist(ltm_varcheck$dstWater)
# gap between 1200 and 1700m. Not a big problem

# dstStlmnt 
hist(segdata$dstStlmnt)
hist(ltm_varcheck$dstStlmnt)
# Couple of small gaps between 7-8km and 9-10km.  

# dstRoad 
hist(segdata$dstRoad)
hist(ltm_varcheck$dstRoad)
# Not using this variable anymore 

# dstBorder 
hist(segdata$dstBorder)
hist(ltm_varcheck$dstBorder)
# Fine

# dstStation 
hist(segdata$dstStation)
hist(ltm_varcheck$dstStation)
# couple of small gaps - nothing major

# elevation 
hist(segdata$elevation)
hist(ltm_varcheck$elevation)
# one small gap between 650m and 700m. Not big problem

## Histograms

# Distance 
ltm_h1 <- ggplot(ltm_distdata, aes(distance))+ geom_histogram(binwidth = 1)
ltm_h2 <- ggplot(ltm_distdata, aes(distance))+ geom_histogram(binwidth = 5)
ltm_h3 <- ggplot(ltm_distdata, aes(distance))+ geom_histogram(binwidth = 10)
ltm_h4 <- ggplot(ltm_distdata, aes(distance))+ geom_histogram(binwidth = 15)
ltm_h5 <- ggplot(ltm_distdata, aes(distance))+ geom_histogram(binwidth = 20)
ltm_h6 <- ggplot(ltm_distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(ltm_h1,ltm_h2,ltm_h3,ltm_h4,ltm_h5,ltm_h6)
#  Evidence of evasive movement between perhaps 5 and 10, or 5 and 15m. Apart from that the histograms look pretty good


# cluster size, observer, habitat, year, month, transect
ltm_h7 <- ggplot(ltm_distdata, aes(size))+geom_histogram(binwidth = 0.5)
ltm_h8 <- ggplot(ltm_distdata, aes(NameObserver))+geom_histogram(stat="count")
ltm_h9 <- ggplot(ltm_distdata, aes(habitat))+geom_histogram(stat="count")
ltm_h10 <- ggplot(ltm_distdata, aes(year))+geom_histogram(stat="count")
ltm_h11 <- ggplot(ltm_distdata, aes(month))+geom_histogram(stat="count")
ltm_h12 <- ggplot(ltm_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(ltm_h7,ltm_h8,ltm_h9,ltm_h10,ltm_h11,ltm_h12)
# The majority of cluster sizes are between 1 and 10 individuals, but the maximum is 40. This species is known to exist in very large groups, so the recorded cluster sizes of <10 are likely underetimates of true group size. Nice coverage of observations in dense forest, open forest, and some in non-forest, highlighting that the species is pretty generalist in terms of habitat requirements. There was a major dip in observations in 2013 for some reason. I don't have an explanation for this. 

## Plots of distance against variables

plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

ltm_d1 <- ggplot(ltm_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
ltm_d2 <- ggplot(ltm_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
ltm_d3 <- ggplot(ltm_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
ltm_d4 <- ggplot(ltm_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
ltm_d5 <- ggplot(ltm_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(ltm_d1,ltm_d2,ltm_d3,ltm_d4,ltm_d5)
#  Slightly larger mean distance in open forest which we would expect. There may be a relationship between group size and distance, but we will see in the model below


## Plots of cluster size against variables
ltm_s1 <- ggplot(ltm_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
ltm_s2 <- ggplot(ltm_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
ltm_s3 <- ggplot(ltm_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
ltm_s4 <- ggplot(ltm_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
ltm_s5 <- ggplot(ltm_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(ltm_s1,ltm_s2,ltm_s3,ltm_s4,ltm_s5)
# mean group size in 2018 is much larger than the other years for some reason. Mean group size is also larger in NF than D or O forest.

## Linear models

# group size ~ distance
newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=ltm_distdata)
plot(ltm_distdata$size~ltm_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# There is no significant relationship between group size and distance

## Estimating the detection function ####


### the DF used in CDS is hn.strat.size.bin.

ltmDF.hn.strat.size.bin <- ds(ltm_distdata, truncation = 50, key="hn", formula=~stratum+size,
                              cutpoints = c(0,10,20,27,35,50))



 
## Fitting a spatial model ####

## The best DF is ltmDF.hn.cov8 

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated.

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
ltm_obsdata <- ltm_obsdata %>% filter(distance <= 50)


  ## Quasipoisson response ####

# saturated model
ltmDSM.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+
                         s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                         habitat,
                  ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.sat)
par(mfrow=c(2,3))
plot(ltmDSM.sat, scale = 0)
gam.check(ltmDSM.sat)
# DE = 23.9.  dstStlmnt not sig. everything looking overfitted

# reduce k
ltmDSM.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=5)+ 
                          s(elevation,bs="ts",k=5)+ s(dstBorder,bs="ts",k=5)+ 
                          s(dstStation,bs="ts",k=5)+  
                          habitat,
                  ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.sat2)
par(mfrow=c(2,3))
plot(ltmDSM.sat2, scale = 0)
gam.check(ltmDSM.sat2)
# DE = 19.5. dstStlmnt still only non-sig. Plots looking better

# remove dstStlmnt
ltmDSM.qp.3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                          s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+  
                          habitat,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.qp.3)
par(mfrow=c(2,2))
plot(ltmDSM.qp.3, scale = 0)
gam.check(ltmDSM.qp.3)
# DE = 19. All terms sig.


### The best QP model is ltmDSM.qp.3


  ## Tweedie Response ####

# saturated model
ltmDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                            habitat,
                  ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.tw.sat)
par(mfrow=c(2,3))
plot(ltmDSM.tw.sat, scale = 0)
gam.check(ltmDSM.tw.sat)
# AIC = 14077, DE = 18.9.  all terms sig excpet dstStlmnt. elevation and dstStation linear

# try increasing k for dstSltmnt, elevation, dstStation
ltmDSM.tw.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts",k=10)+ 
                             s(elevation,bs="ts",k=10)+ s(dstBorder,bs="ts")+ 
                             s(dstStation,bs="ts",k=10)+  
                             habitat,
                     ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                     family = tw(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.tw.sat2)
par(mfrow=c(2,3))
plot(ltmDSM.tw.sat2, scale = 0)
gam.check(ltmDSM.tw.sat2)
# AIC = 14077, DE = 18.9. no change

# remove dstStlmt
ltmDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+  s(elevation,bs="ts")+ 
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                      ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                      family = tw(), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.tw.3)
par(mfrow=c(2,2))
plot(ltmDSM.tw.3, scale = 0)
gam.check(ltmDSM.tw.3)
# AIC = 14075, DE = 18.4. Better model via AIC. All terms sig. Elevation & dstStation look linear

# dstStation as linear term
ltmDSM.tw.4 <- dsm(Nhat ~ s(dstWater,bs="ts")+  s(elevation,bs="ts")+ 
                          s(dstBorder,bs="ts")+  
                          habitat + dstStation,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.tw.4)
par(mfrow=c(2,2))
plot(ltmDSM.tw.4, scale = 0)
gam.check(ltmDSM.tw.4)
# AIC = 14078, DE = 16.9. Ran with warning, AIC gone up and DE gone down.

# put dstStation back as smooth term and elevation as linear term
ltmDSM.tw.5 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstBorder,bs="ts")+ 
                          s(dstStation,bs="ts")+  
                          habitat + elevation,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.tw.5)
par(mfrow=c(2,2))
plot(ltmDSM.tw.5, scale = 0)
gam.check(ltmDSM.tw.5)
# AIC = 14075, DE = 18.4. ran with warning but AIC and DE same as tw.3.

# try with both elevation and dstStation as linear terms
ltmDSM.tw.6 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstBorder,bs="ts")+ 
                          habitat + elevation + dstStation,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.tw.6)
par(mfrow=c(2,2))
plot(ltmDSM.tw.6, scale = 0)
gam.check(ltmDSM.tw.6)
# AIC = 14078, DE = 17.1. Worse model.

# back to tw.3 but with k increased for elevation & dstStation
ltmDSM.tw.7 <- dsm(Nhat ~ s(dstWater,bs="ts")+  s(elevation,bs="ts",k=15)+ 
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts",k=15)+  
                          habitat,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.tw.7)
par(mfrow=c(2,2))
plot(ltmDSM.tw.7, scale = 0)
gam.check(ltmDSM.tw.7)
# AIC = 14077, DE = 19.7. Worse model by AIC. dstStation still linear.

# reduce k for elevation
ltmDSM.tw.8 <- dsm(Nhat ~ s(dstWater,bs="ts")+  s(elevation,bs="ts",k=10)+ 
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts",k=15)+  
                          habitat,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.tw.8)
par(mfrow=c(2,2))
plot(ltmDSM.tw.8, scale = 0)
gam.check(ltmDSM.tw.8)
# AIC = 14075, DE = 18.4. Ran with warning but lower AIC. Making elevaiton and dstStation linear didn't improve the model, and I don't think increasing K did either. I'm going to stick with tw.3


### ltmDSM.tw.3 is the best TW model


  ## Negative binomial response ####

# saturated model
ltmDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+  
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                            habitat,
                  ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.nb.sat)
par(mfrow=c(2,3))
plot(ltmDSM.nb.sat, scale = 0)
gam.check(ltmDSM.nb.sat)
# AIC = 1633, DE = 25.2. dstStlmnt and elevation not sig. 

# remove dstStlmnt
ltmDSM.nb.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts")+  
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                     ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                     family = nb(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.nb.2)
par(mfrow=c(2,2))
plot(ltmDSM.nb.2, scale = 0)
gam.check(ltmDSM.nb.2)
# AIC = 1634, DE = 24.8. elevation not quite sig, and looking linear

# increase k for elevation
ltmDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts", k=15)+  
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.nb.3)
par(mfrow=c(2,2))
plot(ltmDSM.nb.3, scale = 0)
gam.check(ltmDSM.nb.3)
# AIC = 1634, DE = 24.8. increasing k makes no difference

# make elevation linear
ltmDSM.nb.4 <- dsm(Nhat ~ s(dstWater,bs="ts")+  s(dstBorder,bs="ts")+ 
                          s(dstStation,bs="ts")+  habitat + elevation,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.nb.4)
par(mfrow=c(2,2))
plot(ltmDSM.nb.4, scale = 0)
gam.check(ltmDSM.nb.4)
# AIC = 1635, DE = 24.8. Model slightly worse

# elevation as smooth, dstStation as linear
ltmDSM.nb.5 <- dsm(Nhat ~ s(dstWater,bs="ts")+  s(dstBorder,bs="ts")+ 
                          s(elevation,bs="ts")+  habitat + dstStation,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.nb.5)
par(mfrow=c(2,2))
plot(ltmDSM.nb.5, scale = 0)
gam.check(ltmDSM.nb.5)
# AIC = 1631, DE = 32.5. Ran with warning but AIC down and DE up.

# remove elevation
ltmDSM.nb.6 <- dsm(Nhat ~ s(dstWater,bs="ts")+  s(dstBorder,bs="ts")+ 
                          habitat + dstStation,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.nb.6)
par(mfrow=c(2,2))
plot(ltmDSM.nb.6, scale = 0)
gam.check(ltmDSM.nb.6)
# AIC = 1633, DE = 10.8. Woah, large drop in DE.



### ltmDSM.nb.5 is the best NB model

## Model selection ####

summary(ltmDSM.qp.3) # DE = 19
summary(ltmDSM.tw.3) # DE = 18.4
summary(ltmDSM.nb.5) # DE = 32.5

ltmDSM.tw.3$aic # 14075
ltmDSM.nb.5$aic # 1631
# NB model is better than TW

anova(ltmDSM.nb.5,ltmDSM.qp.3,test="Chisq")
# NB model has a lot less residual deviance

par(mfrow=c(2,2))
gam.check(ltmDSM.qp.3)
gam.check(ltmDSM.nb.5)
# Q-Q plot for NB model is much better

### ltmDSM.nb.5 is the best LTM modl

## Autocorrelation ####

par(mfrow=c(1,1))
dsm.cor(ltmDSM.nb.5, max.lag=15, Segment.Label="Sample.Label")
# there is some minor autocorrelation

# First recode sample labels and transect labels as numeric
segdata$sg.id <- as.numeric(sub("\\d+-","",segdata$Sample.Label))
segdata$tr.id <- as.numeric(segdata$Transect.Label)

# Try autoregressive structure
ltmDSM.nb.5.cor <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstBorder,bs="ts")+
                              s(elevation,bs="ts")+
                              habitat + dstStation,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   engine = "gamm", correlation=corAR1(form=~sg.id|tr.id), niterPQL=50,
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.nb.5.cor$lme)
summary(ltmDSM.nb.5.cor$gam)

# test autocorrelation for new model
dsm.cor(ltmDSM.nb.5.cor, max.lag=15, Segment.Label="Sample.Label")

par(mfrow=c(2,1))
dsm.cor(ltmDSM.nb.5, max.lag=15, Segment.Label="Sample.Label")
dsm.cor(ltmDSM.nb.5.cor, max.lag=15, Segment.Label="Sample.Label")

## not conviced that correlation strucutre has improved the autocor at all.

# try bivariate xy smooth
ltmDSM.nb.5.xy <- dsm(Nhat ~ s(dstWater,bs="ts")+  s(dstBorder,bs="ts")+ 
                             s(elevation,bs="ts")+ s(x,y,bs="ts")+  
                             habitat + dstStation,
                   ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.nb.5.xy)
par(mfrow=c(2,2))
plot(ltmDSM.nb.5.xy, scale = 0)
gam.check(ltmDSM.nb.5.xy)
# AIC = 1632, DE = 24.9. Worse model

# try univariate smooths
ltmDSM.nb.5.x.y <- dsm(Nhat ~ s(dstWater,bs="ts")+  s(dstBorder,bs="ts")+ 
                              s(elevation,bs="ts")+ s(x,bs="ts")+ s(y,bs="ts")+
                              habitat + dstStation,
                      ltmDF.hn.strat.size.bin, segdata, ltm_obsdata, method = "REML",
                      family = nb(), engine = "gam",
                      segment.area = segdata$Sample.Area, group = TRUE)
summary(ltmDSM.nb.5.x.y)
par(mfrow=c(2,2))
plot(ltmDSM.nb.5.x.y, scale = 0)
gam.check(ltmDSM.nb.5.x.y)
# AIC = 1632, DE = 24.9. 

## ok so the smooths aren't doing anything.  The autocorrelation is very very minor, and so I will jst move on the with ltmDSM.nb.5 model


## Abundance estimation ####


# Predict over 2010 habitat 
ltm.global.pred10.core <- predict(ltmDSM.nb.5, preddata10_core, off.set = 40000)
write.csv(ltm.global.pred10.core, file="Results/LTM/core_only/ltm.pred10.core.csv")

# Predict over 2011 habitat 
ltm.global.pred11.core <- predict(ltmDSM.nb.5, preddata11_core, off.set = 40000)
write.csv(ltm.global.pred11.core, file="Results/LTM/core_only/ltm.pred11.core.csv")

# Predict over 2013 habitat 
ltm.global.pred13.core <- predict(ltmDSM.nb.5, preddata13_core, off.set = 40000)
write.csv(ltm.global.pred13.core, file="Results/LTM/core_only/ltm.pred13.core.csv")

# Predict over 2014 habitat 
ltm.global.pred14.core <- predict(ltmDSM.nb.5, preddata14_core, off.set = 40000)
write.csv(ltm.global.pred14.core, file="Results/LTM/core_only/ltm.pred14.core.csv")

# Predict over 2016 habitat 
ltm.global.pred16.core <- predict(ltmDSM.nb.5, preddata16_core, off.set = 40000)
write.csv(ltm.global.pred16.core, file="Results/LTM/core_only/ltm.pred16.core.csv")

# Predict over 2018 habitat 
ltm.global.pred18.core <- predict(ltmDSM.nb.5, preddata18_core, off.set = 40000)
write.csv(ltm.global.pred18.core, file="Results/LTM/core_only/ltm.pred18.core.csv")

# Predict over 2020 habitat 
ltm.global.pred20.core <- predict(ltmDSM.nb.5, preddata20_core, off.set = 40000)
write.csv(ltm.global.pred20.core, file="Results/LTM/core_only/ltm.pred20.core.csv")


# Create new dataframes for plotting
ltm.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = ltm.global.pred10.core)

ltm.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = ltm.global.pred11.core)

ltm.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = ltm.global.pred13.core)

ltm.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = ltm.global.pred14.core)

ltm.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = ltm.global.pred16.core)

ltm.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = ltm.global.pred18.core)

ltm.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = ltm.global.pred20.core)

## This creates a dataframe that can be plotted as a map
ltm.spdf.df_10.core <- abunPlotDF(ltm.df.Final10.core, pred.polys_200)
ltm.spdf.df_11.core <- abunPlotDF(ltm.df.Final11.core, pred.polys_200)
ltm.spdf.df_13.core <- abunPlotDF(ltm.df.Final13.core, pred.polys_200)
ltm.spdf.df_14.core <- abunPlotDF(ltm.df.Final14.core, pred.polys_200)
ltm.spdf.df_16.core <- abunPlotDF(ltm.df.Final16.core, pred.polys_200)
ltm.spdf.df_18.core <- abunPlotDF(ltm.df.Final18.core, pred.polys_200)
ltm.spdf.df_20.core <- abunPlotDF(ltm.df.Final20.core, pred.polys_200)

# save SPDFs
write.csv(ltm.spdf.df_10.core,file="Results/LTM/Plots/core_only/spdf/ltm.spdf.df_10.core.csv")
write.csv(ltm.spdf.df_11.core,file="Results/LTM/Plots/core_only/spdf/ltm.spdf.df_11.core.csv")
write.csv(ltm.spdf.df_13.core,file="Results/LTM/Plots/core_only/spdf/ltm.spdf.df_13.core.csv")
write.csv(ltm.spdf.df_14.core,file="Results/LTM/Plots/core_only/spdf/ltm.spdf.df_14.core.csv")
write.csv(ltm.spdf.df_16.core,file="Results/LTM/Plots/core_only/spdf/ltm.spdf.df_16.core.csv")
write.csv(ltm.spdf.df_18.core,file="Results/LTM/Plots/core_only/spdf/ltm.spdf.df_18.core.csv")
write.csv(ltm.spdf.df_20.core,file="Results/LTM/Plots/core_only/spdf/ltm.spdf.df_20.core.csv")



    ## Plotting continuous ####

# Load spatial dataframes
ltm.spdf.df_10.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_10.core.csv")
#ltm.spdf.df_11.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_11.core.csv") 
#ltm.spdf.df_13.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_13.core.csv") 
#ltm.spdf.df_14.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_14.core.csv") 
#ltm.spdf.df_16.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_16.core.csv") 
#ltm.spdf.df_18.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_18.core.csv") 
ltm.spdf.df_20.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_20.core.csv")

# greyscale plots
LTM_plot_10_core_gr <- GSplotFun(ltm.spdf.df_10.core, survey.area.core, "abundance", "2010")
LTM_plot_11_core_gr <- GSplotFun(ltm.spdf.df_11.core, survey.area.core, "abundance")
LTM_plot_13_core_gr <- GSplotFun(ltm.spdf.df_13.core, survey.area.core, "abundance")
LTM_plot_14_core_gr <- GSplotFun(ltm.spdf.df_14.core, survey.area.core, "abundance")
LTM_plot_16_core_gr <- GSplotFun(ltm.spdf.df_16.core, survey.area.core, "abundance")
LTM_plot_18_core_gr <- GSplotFun(ltm.spdf.df_18.core, survey.area.core, "abundance")
LTM_plot_20_core_gr <- GSplotFun(ltm.spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(LTM_plot_10_core_gr,"Results/LTM/Plots/core_only/greyscale/LTM_plot_10_core_gr.png")
saveplot(LTM_plot_11_core_gr,"Results/LTM/Plots/core_only/greyscale/LTM_plot_11_core_gr.png")
saveplot(LTM_plot_13_core_gr,"Results/LTM/Plots/core_only/greyscale/LTM_plot_13_core_gr.png")
saveplot(LTM_plot_14_core_gr,"Results/LTM/Plots/core_only/greyscale/LTM_plot_14_core_gr.png")
saveplot(LTM_plot_16_core_gr,"Results/LTM/Plots/core_only/greyscale/LTM_plot_16_core_gr.png")
saveplot(LTM_plot_18_core_gr,"Results/LTM/Plots/core_only/greyscale/LTM_plot_18_core_gr.png")
saveplot(LTM_plot_20_core_gr,"Results/LTM/Plots/core_only/greyscale/LTM_plot_20_core_gr.png")

# colour plots
LTM_plot_10_core_col <- CLplotFun(ltm.spdf.df_10.core, survey.area.core, "abundance")
LTM_plot_11_core_col <- CLplotFun(ltm.spdf.df_11.core, survey.area.core, "abundance")
LTM_plot_13_core_col <- CLplotFun(ltm.spdf.df_13.core, survey.area.core, "abundance")
LTM_plot_14_core_col <- CLplotFun(ltm.spdf.df_14.core, survey.area.core, "abundance")
LTM_plot_16_core_col <- CLplotFun(ltm.spdf.df_16.core, survey.area.core, "abundance")
LTM_plot_18_core_col <- CLplotFun(ltm.spdf.df_18.core, survey.area.core, "abundance")
LTM_plot_20_core_col <- CLplotFun(ltm.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(LTM_plot_10_core_col,"Results/LTM/Plots/core_only/colour/LTM_plot_10_core_col.png")
saveplot(LTM_plot_11_core_col,"Results/LTM/Plots/core_only/colour/LTM_plot_11_core_col.png")
saveplot(LTM_plot_13_core_col,"Results/LTM/Plots/core_only/colour/LTM_plot_13_core_col.png")
saveplot(LTM_plot_14_core_col,"Results/LTM/Plots/core_only/colour/LTM_plot_14_core_col.png")
saveplot(LTM_plot_16_core_col,"Results/LTM/Plots/core_only/colour/LTM_plot_16_core_col.png")
saveplot(LTM_plot_18_core_col,"Results/LTM/Plots/core_only/colour/LTM_plot_18_core_col.png")
saveplot(LTM_plot_20_core_col,"Results/LTM/Plots/core_only/colour/LTM_plot_20_core_col.png")



## plot grids (abundance and variance)

# greyscale 
ltm_2yrs_gs <- 
  LTM_plot_10_core_gr + LTM_varplot_final10.core.bw  +
  LTM_plot_20_core_gr + LTM_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
ltm_2yrs_gs[[1]] <- ltm_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
ltm_2yrs_gs[[2]] <- ltm_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
ltm_2yrs_gs[[2]] <- ltm_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
ltm_2yrs_gs[[4]] <- ltm_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(ltm_2yrs_gs, "Results/LTM/Plots/core_only/greyscale/plot_grids/ltm_2yrs_gs.png")


    ## Plotting discrete ####
      # Add bins to SPDF - don't repeat ####

# this is the process of adding discrete bins to the abundance SPDFs. I will save them in a new folder so this only has to be done once

# Load original spatial dataframes
ltm.spdf.df_10.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_10.core.csv")
ltm.spdf.df_11.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_11.core.csv")
ltm.spdf.df_13.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_13.core.csv")
ltm.spdf.df_14.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_14.core.csv")
ltm.spdf.df_16.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_16.core.csv")
ltm.spdf.df_18.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_18.core.csv")
ltm.spdf.df_20.core <- read.csv("Results/LTM/Plots/core_only/spdf/ltm.spdf.df_20.core.csv")

# put spdf's into a list
dfs <- list(ltm.spdf.df_10.core,ltm.spdf.df_11.core,ltm.spdf.df_13.core,ltm.spdf.df_14.core,
            ltm.spdf.df_16.core,ltm.spdf.df_18.core,ltm.spdf.df_20.core)

# name the elements
names(dfs) <- c("ltm.spdf.df_10.core","ltm.spdf.df_11.core","ltm.spdf.df_13.core","ltm.spdf.df_14.core",
                "ltm.spdf.df_16.core","ltm.spdf.df_18.core","ltm.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunAbun)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(ltm.spdf.df_10.core,file="Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_10.core.csv")
write.csv(ltm.spdf.df_11.core,file="Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_11.core.csv")
write.csv(ltm.spdf.df_13.core,file="Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_13.core.csv")
write.csv(ltm.spdf.df_14.core,file="Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_14.core.csv")
write.csv(ltm.spdf.df_16.core,file="Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_16.core.csv")
write.csv(ltm.spdf.df_18.core,file="Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_18.core.csv")
write.csv(ltm.spdf.df_20.core,file="Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_20.core.csv")


      # Plotting ####

# Load spatial dataframes (with bins)
ltm.spdf.df_10.core <- read.csv("Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_10.core.csv")
ltm.spdf.df_11.core <- read.csv("Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_11.core.csv")
ltm.spdf.df_13.core <- read.csv("Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_13.core.csv")
ltm.spdf.df_14.core <- read.csv("Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_14.core.csv")
ltm.spdf.df_16.core <- read.csv("Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_16.core.csv")
ltm.spdf.df_18.core <- read.csv("Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_18.core.csv")
ltm.spdf.df_20.core <- read.csv("Results/LTM/Plots/core_only/spdf/bins/ltm.spdf.df_20.core.csv")

# change group2 (abundance) to factor and re-order. I've only done it for 2020 here but you can copy the code for the other years if you need to
ltm.spdf.df_20.core$group2 <- as.factor(ltm.spdf.df_20.core$group2)
ltm.spdf.df_20.core$group2 <- factor(ltm.spdf.df_20.core$group2, levels=c("High","Medium","Low","Very low"))

## plot greyscale
LTM_10_plot_bin_GS <- GSplotBin(ltm.spdf.df_10.core,"group2",survey.area.core,"Abundance","Relative abundance")
LTM_11_plot_bin_GS <- GSplotBin(ltm.spdf.df_11.core,"group2",survey.area.core,"Abundance","Relative abundance")
LTM_13_plot_bin_GS <- GSplotBin(ltm.spdf.df_13.core,"group2",survey.area.core,"Abundance","Relative abundance")
LTM_14_plot_bin_GS <- GSplotBin(ltm.spdf.df_14.core,"group2",survey.area.core,"Abundance","Relative abundance")
LTM_16_plot_bin_GS <- GSplotBin(ltm.spdf.df_16.core,"group2",survey.area.core,"Abundance","Relative abundance")
LTM_18_plot_bin_GS <- GSplotBin(ltm.spdf.df_18.core,"group2",survey.area.core,"Abundance","Relative abundance")
LTM_20_plot_bin_GS <- GSplotBin(ltm.spdf.df_20.core,"group2",survey.area.core,"Abundance","Relative abundance")


# save 
saveplot(LTM_10_plot_bin_GS,"Results/LTM/Plots/core_only/bins/LTM_10_plot_bin_GS.png")
saveplot(LTM_11_plot_bin_GS,"Results/LTM/Plots/core_only/bins/LTM_11_plot_bin_GS.png")
saveplot(LTM_13_plot_bin_GS,"Results/LTM/Plots/core_only/bins/LTM_13_plot_bin_GS.png")
saveplot(LTM_14_plot_bin_GS,"Results/LTM/Plots/core_only/bins/LTM_14_plot_bin_GS.png")
saveplot(LTM_16_plot_bin_GS,"Results/LTM/Plots/core_only/bins/LTM_16_plot_bin_GS.png")
saveplot(LTM_18_plot_bin_GS,"Results/LTM/Plots/core_only/bins/LTM_18_plot_bin_GS.png")
saveplot(LTM_20_plot_bin_GS,"Results/LTM/Plots/core_only/bins/LTM_20_plot_bin_GS.png")



## Variance estimation ####

# estimate variance
ltm.var.Final10.core <- varEstfun(preddata10_core, ltmDSM.nb.5)
ltm.var.Final11.core <- varEstfun(preddata11_core, ltmDSM.nb.5)
ltm.var.Final13.core <- varEstfun(preddata13_core, ltmDSM.nb.5)
ltm.var.Final14.core <- varEstfun(preddata14_core, ltmDSM.nb.5)
ltm.var.Final16.core <- varEstfun(preddata16_core, ltmDSM.nb.5)
ltm.var.Final18.core <- varEstfun(preddata18_core, ltmDSM.nb.5)
ltm.var.Final20.core <- varEstfun(preddata20_core, ltmDSM.nb.5)

# save variance estimates
write.csv(ltm.var.Final10.core, file="Results/LTM/core_only/ltm.var10.core.csv")
write.csv(ltm.var.Final11.core, file="Results/LTM/core_only/ltm.var11.core.csv")
write.csv(ltm.var.Final13.core, file="Results/LTM/core_only/ltm.var13.core.csv")
write.csv(ltm.var.Final14.core, file="Results/LTM/core_only/ltm.var14.core.csv")
write.csv(ltm.var.Final16.core, file="Results/LTM/core_only/ltm.var16.core.csv")
write.csv(ltm.var.Final18.core, file="Results/LTM/core_only/ltm.var18.core.csv")
write.csv(ltm.var.Final20.core, file="Results/LTM/core_only/ltm.var20.core.csv")

# create spdf's for plotting
ltm.var.spdf.df_10.core <- varPlotDF(ltm.var.Final10.core, pred.polys_200)
ltm.var.spdf.df_11.core <- varPlotDF(ltm.var.Final11.core, pred.polys_200)
ltm.var.spdf.df_13.core <- varPlotDF(ltm.var.Final13.core, pred.polys_200)
ltm.var.spdf.df_14.core <- varPlotDF(ltm.var.Final14.core, pred.polys_200)
ltm.var.spdf.df_16.core <- varPlotDF(ltm.var.Final16.core, pred.polys_200)
ltm.var.spdf.df_18.core <- varPlotDF(ltm.var.Final18.core, pred.polys_200)
ltm.var.spdf.df_20.core <- varPlotDF(ltm.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(ltm.var.spdf.df_10.core,
          file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_10.core.csv")
write.csv(ltm.var.spdf.df_11.core,
          file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_11.core.csv")
write.csv(ltm.var.spdf.df_13.core,
          file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_13.core.csv")
write.csv(ltm.var.spdf.df_14.core,
          file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_14.core.csv")
write.csv(ltm.var.spdf.df_16.core,
          file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_16.core.csv")
write.csv(ltm.var.spdf.df_18.core,
          file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_18.core.csv")
write.csv(ltm.var.spdf.df_20.core,
          file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_20.core.csv") 


    # Calculate CV & add bins to SPDF - don't repeat ####

# Load spatial dataframes
ltm.var.spdf.df_10.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_10.core.csv")
ltm.var.spdf.df_11.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_11.core.csv")
ltm.var.spdf.df_13.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_13.core.csv")
ltm.var.spdf.df_14.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_14.core.csv")
ltm.var.spdf.df_16.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_16.core.csv")
ltm.var.spdf.df_18.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_18.core.csv")
ltm.var.spdf.df_20.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_20.core.csv")


# first need to calculate CV from the variance (note: need the abundance spdf's loaded)
ltm.var.spdf.df_10.core <- CVaddFun(ltm.spdf.df_10.core,ltm.var.spdf.df_10.core)
ltm.var.spdf.df_11.core <- CVaddFun(ltm.spdf.df_11.core,ltm.var.spdf.df_11.core)
ltm.var.spdf.df_13.core <- CVaddFun(ltm.spdf.df_13.core,ltm.var.spdf.df_13.core)
ltm.var.spdf.df_14.core <- CVaddFun(ltm.spdf.df_14.core,ltm.var.spdf.df_14.core)
ltm.var.spdf.df_16.core <- CVaddFun(ltm.spdf.df_16.core,ltm.var.spdf.df_16.core)
ltm.var.spdf.df_18.core <- CVaddFun(ltm.spdf.df_18.core,ltm.var.spdf.df_18.core)
ltm.var.spdf.df_20.core <- CVaddFun(ltm.spdf.df_20.core,ltm.var.spdf.df_20.core)


# Save the SPDF's with the CV value but no bins (for continuous plotting of the CV)
write.csv(ltm.var.spdf.df_10.core,file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_10.core.csv")
write.csv(ltm.var.spdf.df_11.core,file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_11.core.csv")
write.csv(ltm.var.spdf.df_13.core,file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_13.core.csv")
write.csv(ltm.var.spdf.df_14.core,file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_14.core.csv")
write.csv(ltm.var.spdf.df_16.core,file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_16.core.csv")
write.csv(ltm.var.spdf.df_18.core,file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_18.core.csv")
write.csv(ltm.var.spdf.df_20.core,file="Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_20.core.csv")


### add bins

## Quartiles

# put spdf's into a list
dfs <- list(ltm.var.spdf.df_10.core,ltm.var.spdf.df_11.core,ltm.var.spdf.df_13.core,ltm.var.spdf.df_14.core,
            ltm.var.spdf.df_16.core,ltm.var.spdf.df_18.core,ltm.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("ltm.var.spdf.df_10.core","ltm.var.spdf.df_11.core","ltm.var.spdf.df_13.core",
                "ltm.var.spdf.df_14.core","ltm.var.spdf.df_16.core","ltm.var.spdf.df_18.core",
                "ltm.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(ltm.var.spdf.df_10.core,file="Results/LTM/Plots/variance/core_only/spdf/bins/ltm.var.spdf.df_10.core.csv")
write.csv(ltm.var.spdf.df_11.core,file="Results/LTM/Plots/variance/core_only/spdf/bins/ltm.var.spdf.df_11.core.csv")
write.csv(ltm.var.spdf.df_13.core,file="Results/LTM/Plots/variance/core_only/spdf/bins/ltm.var.spdf.df_13.core.csv")
write.csv(ltm.var.spdf.df_14.core,file="Results/LTM/Plots/variance/core_only/spdf/bins/ltm.var.spdf.df_14.core.csv")
write.csv(ltm.var.spdf.df_16.core,file="Results/LTM/Plots/variance/core_only/spdf/bins/ltm.var.spdf.df_16.core.csv")
write.csv(ltm.var.spdf.df_18.core,file="Results/LTM/Plots/variance/core_only/spdf/bins/ltm.var.spdf.df_18.core.csv")
write.csv(ltm.var.spdf.df_20.core,file="Results/LTM/Plots/variance/core_only/spdf/bins/ltm.var.spdf.df_20.core.csv")




## custom bins

# change column name from group2 to CV
ltm.var.spdf.df_10.core <- ltm.var.spdf.df_10.core %>% dplyr::rename(CV=group2)
ltm.var.spdf.df_11.core <- ltm.var.spdf.df_11.core %>% dplyr::rename(CV=group2)
ltm.var.spdf.df_13.core <- ltm.var.spdf.df_13.core %>% dplyr::rename(CV=group2)
ltm.var.spdf.df_14.core <- ltm.var.spdf.df_14.core %>% dplyr::rename(CV=group2)
ltm.var.spdf.df_16.core <- ltm.var.spdf.df_16.core %>% dplyr::rename(CV=group2)
ltm.var.spdf.df_18.core <- ltm.var.spdf.df_18.core %>% dplyr::rename(CV=group2)
ltm.var.spdf.df_20.core <- ltm.var.spdf.df_20.core %>% dplyr::rename(CV=group2)


# put spdf's into a list
dfs <- list(ltm.var.spdf.df_10.core,ltm.var.spdf.df_11.core,ltm.var.spdf.df_13.core,ltm.var.spdf.df_14.core,
            ltm.var.spdf.df_16.core,ltm.var.spdf.df_18.core,ltm.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("ltm.var.spdf.df_10.core","ltm.var.spdf.df_11.core","ltm.var.spdf.df_13.core",
                "ltm.var.spdf.df_14.core","ltm.var.spdf.df_16.core","ltm.var.spdf.df_18.core",
                "ltm.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar2)

# split elements into original dataframes
list2env(dfs, globalenv())


# save the SPDFs with the custom bins
write.csv(ltm.var.spdf.df_10.core,
          file="Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_10.core.csv")
write.csv(ltm.var.spdf.df_11.core,
          file="Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_11.core.csv")
write.csv(ltm.var.spdf.df_13.core,
          file="Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_13.core.csv")
write.csv(ltm.var.spdf.df_14.core,
          file="Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_14.core.csv")
write.csv(ltm.var.spdf.df_16.core,
          file="Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_16.core.csv")
write.csv(ltm.var.spdf.df_18.core,
          file="Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_18.core.csv")
write.csv(ltm.var.spdf.df_20.core,
          file="Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_20.core.csv")


    # Discrete bins ####

# load spdfs (which has already had CV calculated and then put into bins)
ltm.var.spdf.df_10.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_10.core.csv")
ltm.var.spdf.df_11.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_11.core.csv")
ltm.var.spdf.df_13.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_13.core.csv")
ltm.var.spdf.df_14.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_14.core.csv")
ltm.var.spdf.df_16.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_16.core.csv")
ltm.var.spdf.df_18.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_18.core.csv")
ltm.var.spdf.df_20.core <- read.csv("Results/LTM/Plots/variance/core_only/spdf/bins/custom/ltm.var.spdf.df_20.core.csv")

# make CV a factor and re-order. I have only done this for 2020 but you can copy the code for the other years if you need to
ltm.var.spdf.df_20.core$CV <- as.factor(ltm.var.spdf.df_20.core$CV)
ltm.var.spdf.df_20.core$CV <- factor(ltm.var.spdf.df_20.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))

# plot CV in bins
LTM_10_plot_bin_GS_var <- GSplotBin(ltm.var.spdf.df_10.core,"CV",survey.area.core,"Variance","CV")
LTM_11_plot_bin_GS_var <- GSplotBin(ltm.var.spdf.df_11.core,"CV",survey.area.core,"Variance","CV")
LTM_13_plot_bin_GS_var <- GSplotBin(ltm.var.spdf.df_13.core,"CV",survey.area.core,"Variance","CV")
LTM_14_plot_bin_GS_var <- GSplotBin(ltm.var.spdf.df_14.core,"CV",survey.area.core,"Variance","CV")
LTM_16_plot_bin_GS_var <- GSplotBin(ltm.var.spdf.df_16.core,"CV",survey.area.core,"Variance","CV")
LTM_18_plot_bin_GS_var <- GSplotBin(ltm.var.spdf.df_18.core,"CV",survey.area.core,"Variance","CV")
LTM_20_plot_bin_GS_var <- GSplotBin(ltm.var.spdf.df_20.core,"CV",survey.area.core,"Variance","CV")

# save 
saveplot(LTM_10_plot_bin_GS_var,"Results/LTM/Plots/variance/core_only/greyscale/bins/LTM_10_plot_bin_GS_var.png")
saveplot(LTM_11_plot_bin_GS_var,"Results/LTM/Plots/variance/core_only/greyscale/bins/LTM_11_plot_bin_GS_var.png")
saveplot(LTM_13_plot_bin_GS_var,"Results/LTM/Plots/variance/core_only/greyscale/bins/LTM_13_plot_bin_GS_var.png")
saveplot(LTM_14_plot_bin_GS_var,"Results/LTM/Plots/variance/core_only/greyscale/bins/LTM_14_plot_bin_GS_var.png")
saveplot(LTM_16_plot_bin_GS_var,"Results/LTM/Plots/variance/core_only/greyscale/bins/LTM_16_plot_bin_GS_var.png")
saveplot(LTM_18_plot_bin_GS_var,"Results/LTM/Plots/variance/core_only/greyscale/bins/LTM_18_plot_bin_GS_var.png")
saveplot(LTM_20_plot_bin_GS_var,"Results/LTM/Plots/variance/core_only/greyscale/bins/LTM_20_plot_bin_GS_var.png")



    # Continuous ####

# Load spatial dataframes
ltm.var.spdf.df_10.core<-read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_10.core.csv")
#ltm.var.spdf.df_11.core<-read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_11.core.csv")
#ltm.var.spdf.df_13.core<-read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_13.core.csv")
#ltm.var.spdf.df_14.core<-read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_14.core.csv")
#ltm.var.spdf.df_16.core<-read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_16.core.csv")
#ltm.var.spdf.df_18.core<-read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_18.core.csv")
ltm.var.spdf.df_20.core<-read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_20.core.csv")

# greyscale plots
LTM_varplot_final10.core.bw <- GSplotFun(ltm.var.spdf.df_10.core, survey.area.core, "variance", "2010")
LTM_varplot_final11.core.bw <- GSplotFun(ltm.var.spdf.df_11.core, survey.area.core, "variance")
LTM_varplot_final13.core.bw <- GSplotFun(ltm.var.spdf.df_13.core, survey.area.core, "variance")
LTM_varplot_final14.core.bw <- GSplotFun(ltm.var.spdf.df_14.core, survey.area.core, "variance")
LTM_varplot_final16.core.bw <- GSplotFun(ltm.var.spdf.df_16.core, survey.area.core, "variance")
LTM_varplot_final18.core.bw <- GSplotFun(ltm.var.spdf.df_18.core, survey.area.core, "variance")
LTM_varplot_final20.core.bw <- GSplotFun(ltm.var.spdf.df_20.core, survey.area.core, "variance", "2020")

# save greyscale
saveplot(LTM_varplot_final10.core.bw, "Results/LTM/Plots/variance/core_only/greyscale/2010_LTM_var.core.bw.png")
saveplot(LTM_varplot_final11.core.bw, "Results/LTM/Plots/variance/core_only/greyscale/2011_LTM_var.core.bw.png")
saveplot(LTM_varplot_final13.core.bw, "Results/LTM/Plots/variance/core_only/greyscale/2013_LTM_var.core.bw.png")
saveplot(LTM_varplot_final14.core.bw, "Results/LTM/Plots/variance/core_only/greyscale/2014_LTM_var.core.bw.png")
saveplot(LTM_varplot_final16.core.bw, "Results/LTM/Plots/variance/core_only/greyscale/2016_LTM_var.core.bw.png")
saveplot(LTM_varplot_final18.core.bw, "Results/LTM/Plots/variance/core_only/greyscale/2018_LTM_var.core.bw.png")
saveplot(LTM_varplot_final20.core.bw, "Results/LTM/Plots/variance/core_only/greyscale/2020_LTM_var.core.bw.png")

# colour plots
LTM_varplot_final10.core.col <- CLplotFun(ltm.var.spdf.df_10.core, survey.area.core, "variance")
LTM_varplot_final11.core.col <- CLplotFun(ltm.var.spdf.df_11.core, survey.area.core, "variance")
LTM_varplot_final13.core.col <- CLplotFun(ltm.var.spdf.df_13.core, survey.area.core, "variance")
LTM_varplot_final14.core.col <- CLplotFun(ltm.var.spdf.df_14.core, survey.area.core, "variance")
LTM_varplot_final16.core.col <- CLplotFun(ltm.var.spdf.df_16.core, survey.area.core, "variance")
LTM_varplot_final18.core.col <- CLplotFun(ltm.var.spdf.df_18.core, survey.area.core, "variance")
LTM_varplot_final20.core.col <- CLplotFun(ltm.var.spdf.df_20.core, survey.area.core, "variance")

# save colour
saveplot(LTM_varplot_final10.core.col, "Results/LTM/Plots/variance/core_only/colour/2010_LTM_var.core.col.png")
saveplot(LTM_varplot_final11.core.col, "Results/LTM/Plots/variance/core_only/colour/2011_LTM_var.core.col.png")
saveplot(LTM_varplot_final13.core.col, "Results/LTM/Plots/variance/core_only/colour/2013_LTM_var.core.col.png")
saveplot(LTM_varplot_final14.core.col, "Results/LTM/Plots/variance/core_only/colour/2014_LTM_var.core.col.png")
saveplot(LTM_varplot_final16.core.col, "Results/LTM/Plots/variance/core_only/colour/2016_LTM_var.core.col.png")
saveplot(LTM_varplot_final18.core.col, "Results/LTM/Plots/variance/core_only/colour/2018_LTM_var.core.col.png")
saveplot(LTM_varplot_final20.core.col, "Results/LTM/Plots/variance/core_only/colour/2020_LTM_var.core.col.png")



#### Pig-tailed macaque ######################################################
## Load data ####

# Observation data. Unique to species 
ptm_obsdata <- read.csv("Species_Data/PTM/R Data/obsdata.csv", header = TRUE)
ptm_obsdata$object <- as.factor(ptm_obsdata$object)
ptm_obsdata$Sample.Label <- as.factor(ptm_obsdata$Sample.Label)
str(ptm_obsdata)
head(ptm_obsdata)

# Transect data. Unique to species. stratum var needed for DF
ptm_distdata <- read.csv("Species_Data/PTM/R Data/distdata.csv", header = TRUE)
ptm_distdata$object <- as.factor(ptm_distdata$object)
ptm_distdata$NameObserver <- as.factor(ptm_distdata$NameObserver)
ptm_distdata$transect <- as.factor(ptm_distdata$transect)
ptm_distdata$stratum <- as.vector(scale(ptm_distdata$year, scale=T, center=T))
ptm_distdata$year <- as.factor(ptm_distdata$year)
ptm_distdata$date <- as.Date(ptm_distdata$date, format = "%d/%m/%Y")
str(ptm_distdata)
head(ptm_distdata)

# check there is no T20
ptm_distdata[ptm_distdata$transect=="20",]
# none


## Plot the covariates across the grid with group sizes ####

## Warning - plots take a while to run

# habitat
plot_PTMobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata200$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=ptm_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/PTM/plot_PTMobs_habitat.png", plot_PTMobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_PTMobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys_200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=ptm_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/PTM/plot_PTMobs_dstWater.png", plot_PTMobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_PTMobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=ptm_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/PTM/plot_PTMobs_dstStlmnt.png", plot_PTMobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_PTMobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=ptm_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/PTM/plot_PTMobs_dstRoad.png", plot_PTMobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_PTMobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=ptm_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/PTM/plot_PTMobs_dstBorder.png", plot_PTMobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_PTMobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=ptm_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/PTM/plot_PTMobs_dstStation.png", plot_PTMobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_PTMobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=ptm_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/PTM/plot_PTMobs_dstELC.png", plot_PTMobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_PTMobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys_200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=ptm_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/PTM/plot_PTMobs_elev.png", plot_PTMobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)


## Exploratory plots and linear models ####


## Checking that there are no large gaps in the range of variables for LTM observations. 

# subset segdata to get only the segments with LTM observations
ptm_varcheck <- segdata[match(ptm_obsdata$Sample.Label,segdata$Sample.Label), ]

# habitat 
plot(segdata$habitat)
plot(ptm_varcheck$habitat)
# Most observations in dense forest, but some in open 

# dstWater 
hist(segdata$dstWater)
hist(ptm_varcheck$dstWater)
# no gaps 

# dstStlmnt 
hist(segdata$dstStlmnt)
hist(ptm_varcheck$dstStlmnt)
# Gap between 8-10km, not a major problem   

# dstRoad 
hist(segdata$dstRoad)
hist(ptm_varcheck$dstRoad)
# This variable not being used anymore

# dstBorder 
hist(segdata$dstBorder)
hist(ptm_varcheck$dstBorder)
# No obs beyond 45km, but better than most other species  

# dstStation 
hist(segdata$dstStation)
hist(ptm_varcheck$dstStation)
# No gaps 

# elevation 
hist(segdata$elevation)
hist(ptm_varcheck$elevation)
# No obs beyond 650m - not a problem 

## Histograms

# Distance 
ptm_h1 <- ggplot(ptm_distdata, aes(distance))+ geom_histogram(binwidth = 1)
ptm_h2 <- ggplot(ptm_distdata, aes(distance))+ geom_histogram(binwidth = 5)
ptm_h3 <- ggplot(ptm_distdata, aes(distance))+ geom_histogram(binwidth = 10)
ptm_h4 <- ggplot(ptm_distdata, aes(distance))+ geom_histogram(binwidth = 15)
ptm_h5 <- ggplot(ptm_distdata, aes(distance))+ geom_histogram(binwidth = 20)
ptm_h6 <- ggplot(ptm_distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(ptm_h1,ptm_h2,ptm_h3,ptm_h4,ptm_h5,ptm_h6)
#  Not a great histogram.  Large 'cliff' in observations beyond 25m.  A hazard rate DF model is likely to work best here.

# cluster size, observer, habitat, year, month, transect
ptm_h7 <- ggplot(ptm_distdata, aes(size))+geom_histogram(binwidth = 0.5)
ptm_h8 <- ggplot(ptm_distdata, aes(NameObserver))+geom_histogram(stat="count")
ptm_h9 <- ggplot(ptm_distdata, aes(habitat))+geom_histogram(stat="count")
ptm_h10 <- ggplot(ptm_distdata, aes(year))+geom_histogram(stat="count")
ptm_h11 <- ggplot(ptm_distdata, aes(month))+geom_histogram(stat="count")
ptm_h12 <- ggplot(ptm_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(ptm_h7,ptm_h8,ptm_h9,ptm_h10,ptm_h11,ptm_h12)
#  vast majority of observations are of group sizes <10.  Likely to be underestimates. Most observations are in desnse forest but there are some obs in the other habitat types which is good. Number of observations have generally increased over the course of the study period.  

## Plots of distance against variables

plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

ptm_d1 <- ggplot(ptm_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
ptm_d2 <- ggplot(ptm_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
ptm_d3 <- ggplot(ptm_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
ptm_d4 <- ggplot(ptm_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
ptm_d5 <- ggplot(ptm_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(ptm_d1,ptm_d2,ptm_d3,ptm_d4,ptm_d5)
#  No real differences in distances between habitats. Not seeing much of a relationship between group size and distance


## Plots of cluster size against variables
ptm_s1 <- ggplot(ptm_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
ptm_s2 <- ggplot(ptm_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
ptm_s3 <- ggplot(ptm_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
ptm_s4 <- ggplot(ptm_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
ptm_s5 <- ggplot(ptm_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(ptm_s1,ptm_s2,ptm_s3,ptm_s4,ptm_s5)
# Nothing particulalry interesting in these plots

## Linear models

# group size ~ distance
newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=ptm_distdata)
plot(ptm_distdata$size~ptm_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# there is a positive relationship between size and distance but the effect size is small and the model is not signficant (p>0.05)

## Estimating the detection function ####

### DF taken from CDS. Binning was tested in the CDS analysis but rejected. The DF model to be used is Hr.strat.size

ptmDF.hr.strat.size <- ds(ptm_distdata, truncation=50, key="hr", formula=~stratum+size)


  
## Fitting a spatial model ####

## The best DF is ptmDF.hr.cov6 

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated.

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
ptm_obsdata <- ptm_obsdata %>% filter(distance <= 50)


  ## Quasipoisson response ####

# saturated model
ptmDSM.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+
                         s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                         habitat,
                  ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.sat)
par(mfrow=c(2,3))
plot(ptmDSM.sat, scale = 0)
gam.check(ptmDSM.sat)
# DE = 15.5.  dstStlmnt not sig. Plots look overfitted

# reduce k
ptmDSM.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(dstStlmnt,bs="ts",k=5)+ 
                          s(elevation,bs="ts",k=5)+ s(dstBorder,bs="ts",k=5)+ 
                          s(dstStation,bs="ts",k=5)+ habitat,
                  ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.sat2)
par(mfrow=c(2,3))
plot(ptmDSM.sat2, scale = 0)
gam.check(ptmDSM.sat2)
# DE = 11.8. Only dstStlmnt still not sig. Plots look better

# try dstsltmnt as linear term
ptmDSM.qp.3 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                          s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                          habitat + dstStlmnt,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.qp.3)
par(mfrow=c(2,2))
plot(ptmDSM.qp.3, scale = 0)
gam.check(ptmDSM.qp.3)
# DE = 11.8. Not sig as a linear term

# remove dstStlmnt
ptmDSM.qp.4 <- dsm(Nhat ~ s(dstWater,bs="ts",k=5)+ s(elevation,bs="ts",k=5)+ 
                          s(dstBorder,bs="ts",k=5)+ s(dstStation,bs="ts",k=5)+ 
                          habitat,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.qp.4)
par(mfrow=c(2,2))
plot(ptmDSM.qp.4, scale = 0)
gam.check(ptmDSM.qp.4)
# DE = 11.7. All terms sig. Plots look good


### The best QP model is ptmDSM.qp.4


  ## Tweedie Response ####

# saturated model
ptmDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                            habitat,
                  ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.tw.sat)
par(mfrow=c(2,3))
plot(ptmDSM.tw.sat, scale = 0)
gam.check(ptmDSM.tw.sat)
# AIC = 14427, DE = 13.1. dstStlmnt and dstStation not sig. dstBorder looks linear. Elevation possibly overfitted. 

# remove dstStation
ptmDSM.tw.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                          s(elevation,bs="ts")+ s(dstBorder,bs="ts")+   
                          habitat,
                     ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                     family = tw(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.tw.2)
par(mfrow=c(2,2))
plot(ptmDSM.tw.2, scale = 0)
gam.check(ptmDSM.tw.2)
# AIC = 14427, DE = 13.1. No change which suggest it wasn't doing anything. dstborder still looks linear, dstStlmnt still not sig. dstWater and elevation overfitted

# reduce k for dstWater and elevation
ptmDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts", k=5)+ s(dstStlmnt,bs="ts")+ 
                          s(elevation,bs="ts", k=5)+ s(dstBorder,bs="ts")+   
                          habitat,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.tw.3)
par(mfrow=c(2,2))
plot(ptmDSM.tw.3, scale = 0)
gam.check(ptmDSM.tw.3)
# AIC = 14439, DE = 12. elevation looks better but dstWater not - its gone linear which I don't think it is. dstStlmnt now sig

# increase k for dstwater
ptmDSM.tw.4 <- dsm(Nhat ~ s(dstWater,bs="ts", k=6)+ s(dstStlmnt,bs="ts")+ 
                          s(elevation,bs="ts", k=5)+ s(dstBorder,bs="ts")+   
                          habitat,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.tw.4)
par(mfrow=c(2,2))
plot(ptmDSM.tw.4, scale = 0)
gam.check(ptmDSM.tw.4)
# AIC = 14436, DE = 12.7. model improved. dstborder not linear, but looks overfitted

# remove dstSltmnt
ptmDSM.tw.5 <- dsm(Nhat ~ s(dstWater,bs="ts", k=6)+ s(elevation,bs="ts", k=5)+ 
                          s(dstBorder,bs="ts")+ habitat,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.tw.5)
par(mfrow=c(2,2))
plot(ptmDSM.tw.5, scale = 0)
gam.check(ptmDSM.tw.5)
# AIC = 14434, DE = 11.6. Improved model, although DE gone down. dstBorder linear

# dstBorder as linear term
ptmDSM.tw.6 <- dsm(Nhat ~ s(dstWater,bs="ts", k=6)+ s(elevation,bs="ts", k=5)+ 
                          habitat + dstBorder,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.tw.6)
par(mfrow=c(1,2))
plot(ptmDSM.tw.6, scale = 0)
gam.check(ptmDSM.tw.6)
# AIC = 14434, DE = 11.6. No real change in the model. dstBorder sig as linear term


### ptmDSM.tw.6 is the best TW model


  ## Negative binomial response ####

# saturated model
ptmDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+  
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                            habitat,
                  ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.nb.sat)
par(mfrow=c(2,3))
plot(ptmDSM.nb.sat, scale = 0)
gam.check(ptmDSM.nb.sat)
# AIC = 3605, DE = 19.5. dstStlmnt & dstStation not sig. dstBorder looks linear. elevation bit overfitted

# remove dstStlmnt
ptmDSM.nb.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts")+  
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                     ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                     family = nb(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.nb.2)
par(mfrow=c(2,2))
plot(ptmDSM.nb.2, scale = 0)
gam.check(ptmDSM.nb.2)
# AIC = 3605, DE = 19.5. no change in AIC or DE. dstStation not sig. dstBorder still linear

# remove dstStation
ptmDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts")+  
                          s(dstBorder,bs="ts")+ habitat,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.nb.3)
par(mfrow=c(2,2))
plot(ptmDSM.nb.3, scale = 0)
gam.check(ptmDSM.nb.3)
# AIC = 3605, DE = 19.5. All terms sig. 

# reduce k for elevation
ptmDSM.nb.4 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts", k=5)+  
                          s(dstBorder,bs="ts")+ habitat,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.nb.4)
par(mfrow=c(2,2))
plot(ptmDSM.nb.4, scale = 0)
gam.check(ptmDSM.nb.4)
# AIC = 3629, DE = 17.5. Elevation plot looks better. 

# try dstBorder as linear term
ptmDSM.nb.5 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts", k=5)+  
                          habitat + dstBorder,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.nb.5)
par(mfrow=c(2,2))
plot(ptmDSM.nb.5, scale = 0)
gam.check(ptmDSM.nb.5)
# AIC = 3629, DE = 17.5.


### ptmDSM.nb.5 is the best NB model

## Model selection ####

summary(ptmDSM.qp.4) # DE = 11.7
summary(ptmDSM.tw.6) # DE = 11.6
summary(ptmDSM.nb.5) # DE = 17.5

ptmDSM.tw.6$aic # 14434
ptmDSM.nb.5$aic # 3629
# NB model is better than TW

anova(ptmDSM.qp.4,ptmDSM.nb.5,test="Chisq")
# QP has significantly less residual deviance

par(mfrow=c(2,2))
gam.check(ptmDSM.qp.4)
gam.check(ptmDSM.nb.5)
# Q-Q plot and other diagnistics are better for NB. I have more faith in the NB model as both the NB and the TW model only have elevaiton & dtWater as sig terms, whereas QP reckons dstBorder and dstStation are sig. 

### ptmDSM.nb.5 is the best LTM modl

## Autocorrelation ####

par(mfrow=c(1,1))
dsm.cor(ptmDSM.nb.5, max.lag=15, Segment.Label="Sample.Label")
# there is some autocorrelation at lag 1 and 2

# Try autoregressive structure
ptmDSM.nb.5.cor <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation, bs="ts",k=5)+
                              habitat + dstBorder,
                   ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                   engine = "gamm", correlation=corAR1(form=~sg.id|tr.id), niterPQL=50,
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.nb.5.cor$lme)
summary(ptmDSM.nb.5.cor$gam)
# dstWater no longer sig.

# test autocorrelation for new model
dsm.cor(ptmDSM.nb.5.cor, max.lag=15, Segment.Label="Sample.Label")


## the corAR1 structure has removed the autocorrelation from lag 1, but not from lag 2. In fact it has slightly increased it at lag 2.

# I will re-model nb.5 but include univariate smooths of x and y
ptmDSM.nb.5.x.y <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(elevation,bs="ts",k=5)+  
                              s(x, bs="ts")+ s(y, bs="ts")+
                              habitat + dstBorder,
                  ptmDF.hr.strat.size, segdata, ptm_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(ptmDSM.nb.5.x.y)
par(mfrow=c(2,2))
plot(ptmDSM.nb.5.x.y, scale = 0)
gam.check(ptmDSM.nb.5.x.y)
# AIC = 3608, DE = 19.7. All terms sig. AIC down and DE up. Autocorr should be dealt with via the xy smooths


### ptmDSM.nb.5.x.y is the final PTM model

## Abundance estimation ####


# Predict over 2010 habitat 
ptm.global.pred10.core <- predict(ptmDSM.nb.5.x.y, preddata10_core, off.set = 40000)
write.csv(ptm.global.pred10.core, file="Results/PTM/core_only/ptm.pred10.core.csv")

# Predict over 2011 habitat 
ptm.global.pred11.core <- predict(ptmDSM.nb.5.x.y, preddata11_core, off.set = 40000)
write.csv(ptm.global.pred11.core, file="Results/PTM/core_only/ptm.pred11.core.csv")

# Predict over 2013 habitat 
ptm.global.pred13.core <- predict(ptmDSM.nb.5.x.y, preddata13_core, off.set = 40000)
write.csv(ptm.global.pred13.core, file="Results/PTM/core_only/ptm.pred13.core.csv")

# Predict over 2014 habitat 
ptm.global.pred14.core <- predict(ptmDSM.nb.5.x.y, preddata14_core, off.set = 40000)
write.csv(ptm.global.pred14.core, file="Results/PTM/core_only/ptm.pred14.core.csv")

# Predict over 2016 habitat 
ptm.global.pred16.core <- predict(ptmDSM.nb.5.x.y, preddata16_core, off.set = 40000)
write.csv(ptm.global.pred16.core, file="Results/PTM/core_only/ptm.pred16.core.csv")

# Predict over 2018 habitat 
ptm.global.pred18.core <- predict(ptmDSM.nb.5.x.y, preddata18_core, off.set = 40000)
write.csv(ptm.global.pred18.core, file="Results/PTM/core_only/ptm.pred18.core.csv")

# Predict over 2020 habitat 
ptm.global.pred20.core <- predict(ptmDSM.nb.5.x.y, preddata20_core, off.set = 40000)
write.csv(ptm.global.pred20.core, file="Results/PTM/core_only/ptm.pred20.core.csv")


# Create new dataframes for plotting
ptm.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = ptm.global.pred10.core)

ptm.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = ptm.global.pred11.core)

ptm.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = ptm.global.pred13.core)

ptm.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = ptm.global.pred14.core)

ptm.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = ptm.global.pred16.core)

ptm.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = ptm.global.pred18.core)

ptm.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = ptm.global.pred20.core)


## This creates a dataframe that can be plotted as a map
ptm.spdf.df_10.core <- abunPlotDF(ptm.df.Final10.core, pred.polys_200)
ptm.spdf.df_11.core <- abunPlotDF(ptm.df.Final11.core, pred.polys_200)
ptm.spdf.df_13.core <- abunPlotDF(ptm.df.Final13.core, pred.polys_200)
ptm.spdf.df_14.core <- abunPlotDF(ptm.df.Final14.core, pred.polys_200)
ptm.spdf.df_16.core <- abunPlotDF(ptm.df.Final16.core, pred.polys_200)
ptm.spdf.df_18.core <- abunPlotDF(ptm.df.Final18.core, pred.polys_200)
ptm.spdf.df_20.core <- abunPlotDF(ptm.df.Final20.core, pred.polys_200)

# save SPDFs
write.csv(ptm.spdf.df_10.core,file="Results/PTM/Plots/core_only/spdf/ptm.spdf.df_10.core.csv")
write.csv(ptm.spdf.df_11.core,file="Results/PTM/Plots/core_only/spdf/ptm.spdf.df_11.core.csv")
write.csv(ptm.spdf.df_13.core,file="Results/PTM/Plots/core_only/spdf/ptm.spdf.df_13.core.csv")
write.csv(ptm.spdf.df_14.core,file="Results/PTM/Plots/core_only/spdf/ptm.spdf.df_14.core.csv")
write.csv(ptm.spdf.df_16.core,file="Results/PTM/Plots/core_only/spdf/ptm.spdf.df_16.core.csv")
write.csv(ptm.spdf.df_18.core,file="Results/PTM/Plots/core_only/spdf/ptm.spdf.df_18.core.csv")
write.csv(ptm.spdf.df_20.core,file="Results/PTM/Plots/core_only/spdf/ptm.spdf.df_20.core.csv")



    ## Plotting continuous ####

# Load spatial dataframes
ptm.spdf.df_10.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_10.core.csv")
#ptm.spdf.df_11.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_11.core.csv") 
#ptm.spdf.df_13.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_13.core.csv") 
#ptm.spdf.df_14.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_14.core.csv") 
#ptm.spdf.df_16.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_16.core.csv") 
#ptm.spdf.df_18.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_18.core.csv") 
ptm.spdf.df_20.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_20.core.csv") 

# greyscale plots
PTM_plot_10_core_gr <- GSplotFun(ptm.spdf.df_10.core, survey.area.core, "abundance", "2010")
PTM_plot_11_core_gr <- GSplotFun(ptm.spdf.df_11.core, survey.area.core, "abundance")
PTM_plot_13_core_gr <- GSplotFun(ptm.spdf.df_13.core, survey.area.core, "abundance")
PTM_plot_14_core_gr <- GSplotFun(ptm.spdf.df_14.core, survey.area.core, "abundance")
PTM_plot_16_core_gr <- GSplotFun(ptm.spdf.df_16.core, survey.area.core, "abundance")
PTM_plot_18_core_gr <- GSplotFun(ptm.spdf.df_18.core, survey.area.core, "abundance")
PTM_plot_20_core_gr <- GSplotFun(ptm.spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(PTM_plot_10_core_gr,"Results/PTM/Plots/core_only/greyscale/PTM_plot_10_core_gr.png")
saveplot(PTM_plot_11_core_gr,"Results/PTM/Plots/core_only/greyscale/PTM_plot_11_core_gr.png")
saveplot(PTM_plot_13_core_gr,"Results/PTM/Plots/core_only/greyscale/PTM_plot_13_core_gr.png")
saveplot(PTM_plot_14_core_gr,"Results/PTM/Plots/core_only/greyscale/PTM_plot_14_core_gr.png")
saveplot(PTM_plot_16_core_gr,"Results/PTM/Plots/core_only/greyscale/PTM_plot_16_core_gr.png")
saveplot(PTM_plot_18_core_gr,"Results/PTM/Plots/core_only/greyscale/PTM_plot_18_core_gr.png")
saveplot(PTM_plot_20_core_gr,"Results/PTM/Plots/core_only/greyscale/PTM_plot_20_core_gr.png")

# colour plots
PTM_plot_10_core_col <- CLplotFun(ptm.spdf.df_10.core, survey.area.core, "abundance")
PTM_plot_11_core_col <- CLplotFun(ptm.spdf.df_11.core, survey.area.core, "abundance")
PTM_plot_13_core_col <- CLplotFun(ptm.spdf.df_13.core, survey.area.core, "abundance")
PTM_plot_14_core_col <- CLplotFun(ptm.spdf.df_14.core, survey.area.core, "abundance")
PTM_plot_16_core_col <- CLplotFun(ptm.spdf.df_16.core, survey.area.core, "abundance")
PTM_plot_18_core_col <- CLplotFun(ptm.spdf.df_18.core, survey.area.core, "abundance")
PTM_plot_20_core_col <- CLplotFun(ptm.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(PTM_plot_10_core_col,"Results/PTM/Plots/core_only/colour/PTM_plot_10_core_col.png")
saveplot(PTM_plot_11_core_col,"Results/PTM/Plots/core_only/colour/PTM_plot_11_core_col.png")
saveplot(PTM_plot_13_core_col,"Results/PTM/Plots/core_only/colour/PTM_plot_13_core_col.png")
saveplot(PTM_plot_14_core_col,"Results/PTM/Plots/core_only/colour/PTM_plot_14_core_col.png")
saveplot(PTM_plot_16_core_col,"Results/PTM/Plots/core_only/colour/PTM_plot_16_core_col.png")
saveplot(PTM_plot_18_core_col,"Results/PTM/Plots/core_only/colour/PTM_plot_18_core_col.png")
saveplot(PTM_plot_20_core_col,"Results/PTM/Plots/core_only/colour/PTM_plot_20_core_col.png")



## plot grids (abundance and variance)

# greyscale 
ptm_2yrs_gs <- 
  PTM_plot_10_core_gr + PTM_varplot_final10.core.bw  +
  PTM_plot_20_core_gr + PTM_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
ptm_2yrs_gs[[1]] <- ptm_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
ptm_2yrs_gs[[2]] <- ptm_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
ptm_2yrs_gs[[2]] <- ptm_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
ptm_2yrs_gs[[4]] <- ptm_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(ptm_2yrs_gs, "Results/PTM/Plots/core_only/greyscale/plot_grids/ptm_2yrs_gs.png")


    ## Plotting discrete ####
      # Add bins to SPDF - don't repeat ####

# this is the process of adding discrete bins to the abundance SPDFs. I will save them in a new folder so this only has to be done once

# Load original spatial dataframes
ptm.spdf.df_10.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_10.core.csv")
ptm.spdf.df_11.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_11.core.csv")
ptm.spdf.df_13.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_13.core.csv")
ptm.spdf.df_14.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_14.core.csv")
ptm.spdf.df_16.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_16.core.csv")
ptm.spdf.df_18.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_18.core.csv")
ptm.spdf.df_20.core <- read.csv("Results/PTM/Plots/core_only/spdf/ptm.spdf.df_20.core.csv")

# put spdf's into a list
dfs <- list(ptm.spdf.df_10.core,ptm.spdf.df_11.core,ptm.spdf.df_13.core,ptm.spdf.df_14.core,
            ptm.spdf.df_16.core,ptm.spdf.df_18.core,ptm.spdf.df_20.core)

# name the elements
names(dfs) <- c("ptm.spdf.df_10.core","ptm.spdf.df_11.core","ptm.spdf.df_13.core","ptm.spdf.df_14.core",
                "ptm.spdf.df_16.core","ptm.spdf.df_18.core","ptm.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunAbun)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(ptm.spdf.df_10.core,file="Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_10.core.csv")
write.csv(ptm.spdf.df_11.core,file="Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_11.core.csv")
write.csv(ptm.spdf.df_13.core,file="Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_13.core.csv")
write.csv(ptm.spdf.df_14.core,file="Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_14.core.csv")
write.csv(ptm.spdf.df_16.core,file="Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_16.core.csv")
write.csv(ptm.spdf.df_18.core,file="Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_18.core.csv")
write.csv(ptm.spdf.df_20.core,file="Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_20.core.csv")


      # Plotting ####

# Load spatial dataframes (with bins)
ptm.spdf.df_10.core <- read.csv("Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_10.core.csv")
ptm.spdf.df_11.core <- read.csv("Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_11.core.csv")
ptm.spdf.df_13.core <- read.csv("Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_13.core.csv")
ptm.spdf.df_14.core <- read.csv("Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_14.core.csv")
ptm.spdf.df_16.core <- read.csv("Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_16.core.csv")
ptm.spdf.df_18.core <- read.csv("Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_18.core.csv")
ptm.spdf.df_20.core <- read.csv("Results/PTM/Plots/core_only/spdf/bins/ptm.spdf.df_20.core.csv")

# change group2 (abundance) to factor and re-order. I've only done it for 2020 here but you can copy the code for the other years if you need to
ptm.spdf.df_20.core$group2 <- as.factor(ptm.spdf.df_20.core$group2)
ptm.spdf.df_20.core$group2 <- factor(ptm.spdf.df_20.core$group2, levels=c("High","Medium","Low","Very low"))


## plot greyscale
PTM_10_plot_bin_GS <- GSplotBin(ptm.spdf.df_10.core,"group2",survey.area.core,"Abundance","Relative abundance")
PTM_11_plot_bin_GS <- GSplotBin(ptm.spdf.df_11.core,"group2",survey.area.core,"Abundance","Relative abundance")
PTM_13_plot_bin_GS <- GSplotBin(ptm.spdf.df_13.core,"group2",survey.area.core,"Abundance","Relative abundance")
PTM_14_plot_bin_GS <- GSplotBin(ptm.spdf.df_14.core,"group2",survey.area.core,"Abundance","Relative abundance")
PTM_16_plot_bin_GS <- GSplotBin(ptm.spdf.df_16.core,"group2",survey.area.core,"Abundance","Relative abundance")
PTM_18_plot_bin_GS <- GSplotBin(ptm.spdf.df_18.core,"group2",survey.area.core,"Abundance","Relative abundance")
PTM_20_plot_bin_GS <- GSplotBin(ptm.spdf.df_20.core,"group2",survey.area.core,"Abundance","Relative abundance")


# save 
saveplot(PTM_10_plot_bin_GS,"Results/PTM/Plots/core_only/bins/PTM_10_plot_bin_GS.png")
saveplot(PTM_11_plot_bin_GS,"Results/PTM/Plots/core_only/bins/PTM_11_plot_bin_GS.png")
saveplot(PTM_13_plot_bin_GS,"Results/PTM/Plots/core_only/bins/PTM_13_plot_bin_GS.png")
saveplot(PTM_14_plot_bin_GS,"Results/PTM/Plots/core_only/bins/PTM_14_plot_bin_GS.png")
saveplot(PTM_16_plot_bin_GS,"Results/PTM/Plots/core_only/bins/PTM_16_plot_bin_GS.png")
saveplot(PTM_18_plot_bin_GS,"Results/PTM/Plots/core_only/bins/PTM_18_plot_bin_GS.png")
saveplot(PTM_20_plot_bin_GS,"Results/PTM/Plots/core_only/bins/PTM_20_plot_bin_GS.png")



## Variance estimation ####

# estimate variance
ptm.var.Final10.core <- varEstfun(preddata10_core, ptmDSM.nb.5.x.y)
ptm.var.Final11.core <- varEstfun(preddata11_core, ptmDSM.nb.5.x.y)
ptm.var.Final13.core <- varEstfun(preddata13_core, ptmDSM.nb.5.x.y)
ptm.var.Final14.core <- varEstfun(preddata14_core, ptmDSM.nb.5.x.y)
ptm.var.Final16.core <- varEstfun(preddata16_core, ptmDSM.nb.5.x.y)
ptm.var.Final18.core <- varEstfun(preddata18_core, ptmDSM.nb.5.x.y)
ptm.var.Final20.core <- varEstfun(preddata20_core, ptmDSM.nb.5.x.y)

# save variance estimates
write.csv(ptm.var.Final10.core, file="Results/PTM/core_only/ptm.var10.core.csv")
write.csv(ptm.var.Final11.core, file="Results/PTM/core_only/ptm.var11.core.csv")
write.csv(ptm.var.Final13.core, file="Results/PTM/core_only/ptm.var13.core.csv")
write.csv(ptm.var.Final14.core, file="Results/PTM/core_only/ptm.var14.core.csv")
write.csv(ptm.var.Final16.core, file="Results/PTM/core_only/ptm.var16.core.csv")
write.csv(ptm.var.Final18.core, file="Results/PTM/core_only/ptm.var18.core.csv")
write.csv(ptm.var.Final20.core, file="Results/PTM/core_only/ptm.var20.core.csv")

# create spdf's for plotting
ptm.var.spdf.df_10.core <- varPlotDF(ptm.var.Final10.core, pred.polys_200)
ptm.var.spdf.df_11.core <- varPlotDF(ptm.var.Final11.core, pred.polys_200)
ptm.var.spdf.df_13.core <- varPlotDF(ptm.var.Final13.core, pred.polys_200)
ptm.var.spdf.df_14.core <- varPlotDF(ptm.var.Final14.core, pred.polys_200)
ptm.var.spdf.df_16.core <- varPlotDF(ptm.var.Final16.core, pred.polys_200)
ptm.var.spdf.df_18.core <- varPlotDF(ptm.var.Final18.core, pred.polys_200)
ptm.var.spdf.df_20.core <- varPlotDF(ptm.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(ptm.var.spdf.df_10.core,
          file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_10.core.csv")
write.csv(ptm.var.spdf.df_11.core,
          file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_11.core.csv")
write.csv(ptm.var.spdf.df_13.core,
          file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_13.core.csv")
write.csv(ptm.var.spdf.df_14.core,
          file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_14.core.csv")
write.csv(ptm.var.spdf.df_16.core,
          file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_16.core.csv")
write.csv(ptm.var.spdf.df_18.core,
          file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_18.core.csv")
write.csv(ptm.var.spdf.df_20.core,
          file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_20.core.csv") 


    # Calculate CV & add bins to SPDF - don't repeat ####

# Load spatial dataframes
ptm.var.spdf.df_10.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_10.core.csv")
ptm.var.spdf.df_11.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_11.core.csv")
ptm.var.spdf.df_13.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_13.core.csv")
ptm.var.spdf.df_14.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_14.core.csv")
ptm.var.spdf.df_16.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_16.core.csv")
ptm.var.spdf.df_18.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_18.core.csv")
ptm.var.spdf.df_20.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_20.core.csv")


# first need to calculate CV from the variance (note: need the abundance spdf's loaded)
ptm.var.spdf.df_10.core <- CVaddFun(ptm.spdf.df_10.core,ptm.var.spdf.df_10.core)
ptm.var.spdf.df_11.core <- CVaddFun(ptm.spdf.df_11.core,ptm.var.spdf.df_11.core)
ptm.var.spdf.df_13.core <- CVaddFun(ptm.spdf.df_13.core,ptm.var.spdf.df_13.core)
ptm.var.spdf.df_14.core <- CVaddFun(ptm.spdf.df_14.core,ptm.var.spdf.df_14.core)
ptm.var.spdf.df_16.core <- CVaddFun(ptm.spdf.df_16.core,ptm.var.spdf.df_16.core)
ptm.var.spdf.df_18.core <- CVaddFun(ptm.spdf.df_18.core,ptm.var.spdf.df_18.core)
ptm.var.spdf.df_20.core <- CVaddFun(ptm.spdf.df_20.core,ptm.var.spdf.df_20.core)


# Save the SPDF's with the CV value but no bins (for continuous plotting of the CV)
write.csv(ptm.var.spdf.df_10.core,file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_10.core.csv")
write.csv(ptm.var.spdf.df_11.core,file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_11.core.csv")
write.csv(ptm.var.spdf.df_13.core,file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_13.core.csv")
write.csv(ptm.var.spdf.df_14.core,file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_14.core.csv")
write.csv(ptm.var.spdf.df_16.core,file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_16.core.csv")
write.csv(ptm.var.spdf.df_18.core,file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_18.core.csv")
write.csv(ptm.var.spdf.df_20.core,file="Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_20.core.csv")


### add bins

## Quartiles

# put spdf's into a list
dfs <- list(ptm.var.spdf.df_10.core,ptm.var.spdf.df_11.core,ptm.var.spdf.df_13.core,ptm.var.spdf.df_14.core,
            ptm.var.spdf.df_16.core,ptm.var.spdf.df_18.core,ptm.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("ptm.var.spdf.df_10.core","ptm.var.spdf.df_11.core","ptm.var.spdf.df_13.core",
                "ptm.var.spdf.df_14.core","ptm.var.spdf.df_16.core","ptm.var.spdf.df_18.core",
                "ptm.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar)

# split elements into original dataframes
list2env(dfs, globalenv())

# re-save the spdf's in new folder
write.csv(ptm.var.spdf.df_10.core,file="Results/PTM/Plots/variance/core_only/spdf/bins/ptm.var.spdf.df_10.core.csv")
write.csv(ptm.var.spdf.df_11.core,file="Results/PTM/Plots/variance/core_only/spdf/bins/ptm.var.spdf.df_11.core.csv")
write.csv(ptm.var.spdf.df_13.core,file="Results/PTM/Plots/variance/core_only/spdf/bins/ptm.var.spdf.df_13.core.csv")
write.csv(ptm.var.spdf.df_14.core,file="Results/PTM/Plots/variance/core_only/spdf/bins/ptm.var.spdf.df_14.core.csv")
write.csv(ptm.var.spdf.df_16.core,file="Results/PTM/Plots/variance/core_only/spdf/bins/ptm.var.spdf.df_16.core.csv")
write.csv(ptm.var.spdf.df_18.core,file="Results/PTM/Plots/variance/core_only/spdf/bins/ptm.var.spdf.df_18.core.csv")
write.csv(ptm.var.spdf.df_20.core,file="Results/PTM/Plots/variance/core_only/spdf/bins/ptm.var.spdf.df_20.core.csv")




## custom bins

# change column name from group2 to CV
ptm.var.spdf.df_10.core <- ptm.var.spdf.df_10.core %>% dplyr::rename(CV=group2)
ptm.var.spdf.df_11.core <- ptm.var.spdf.df_11.core %>% dplyr::rename(CV=group2)
ptm.var.spdf.df_13.core <- ptm.var.spdf.df_13.core %>% dplyr::rename(CV=group2)
ptm.var.spdf.df_14.core <- ptm.var.spdf.df_14.core %>% dplyr::rename(CV=group2)
ptm.var.spdf.df_16.core <- ptm.var.spdf.df_16.core %>% dplyr::rename(CV=group2)
ptm.var.spdf.df_18.core <- ptm.var.spdf.df_18.core %>% dplyr::rename(CV=group2)
ptm.var.spdf.df_20.core <- ptm.var.spdf.df_20.core %>% dplyr::rename(CV=group2)


# put spdf's into a list
dfs <- list(ptm.var.spdf.df_10.core,ptm.var.spdf.df_11.core,ptm.var.spdf.df_13.core,ptm.var.spdf.df_14.core,
            ptm.var.spdf.df_16.core,ptm.var.spdf.df_18.core,ptm.var.spdf.df_20.core)

# name the elements
names(dfs) <- c("ptm.var.spdf.df_10.core","ptm.var.spdf.df_11.core","ptm.var.spdf.df_13.core",
                "ptm.var.spdf.df_14.core","ptm.var.spdf.df_16.core","ptm.var.spdf.df_18.core",
                "ptm.var.spdf.df_20.core")

# run the binning function for all spdf's
dfs <- lapply(dfs, binFunVar2)

# split elements into original dataframes
list2env(dfs, globalenv())


# save the SPDFs with the custom bins
write.csv(ptm.var.spdf.df_10.core,
          file="Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_10.core.csv")
write.csv(ptm.var.spdf.df_11.core,
          file="Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_11.core.csv")
write.csv(ptm.var.spdf.df_13.core,
          file="Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_13.core.csv")
write.csv(ptm.var.spdf.df_14.core,
          file="Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_14.core.csv")
write.csv(ptm.var.spdf.df_16.core,
          file="Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_16.core.csv")
write.csv(ptm.var.spdf.df_18.core,
          file="Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_18.core.csv")
write.csv(ptm.var.spdf.df_20.core,
          file="Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_20.core.csv")


    # Discrete bins ####

# load spdfs (which has already had CV calculated and then put into bins)
ptm.var.spdf.df_10.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_10.core.csv")
ptm.var.spdf.df_11.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_11.core.csv")
ptm.var.spdf.df_13.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_13.core.csv")
ptm.var.spdf.df_14.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_14.core.csv")
ptm.var.spdf.df_16.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_16.core.csv")
ptm.var.spdf.df_18.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_18.core.csv")
ptm.var.spdf.df_20.core <- read.csv("Results/PTM/Plots/variance/core_only/spdf/bins/custom/ptm.var.spdf.df_20.core.csv")

# make CV into factor and re-order
ptm.var.spdf.df_10.core$CV <- as.factor(ptm.var.spdf.df_10.core$CV)
ptm.var.spdf.df_10.core$CV <- factor(ptm.var.spdf.df_10.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
ptm.var.spdf.df_11.core$CV <- as.factor(ptm.var.spdf.df_11.core$CV)
ptm.var.spdf.df_11.core$CV <- factor(ptm.var.spdf.df_11.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
ptm.var.spdf.df_13.core$CV <- as.factor(ptm.var.spdf.df_13.core$CV)
ptm.var.spdf.df_13.core$CV <- factor(ptm.var.spdf.df_13.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
ptm.var.spdf.df_14.core$CV <- as.factor(ptm.var.spdf.df_14.core$CV)
ptm.var.spdf.df_14.core$CV <- factor(ptm.var.spdf.df_14.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
ptm.var.spdf.df_16.core$CV <- as.factor(ptm.var.spdf.df_16.core$CV)
ptm.var.spdf.df_16.core$CV <- factor(ptm.var.spdf.df_16.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
ptm.var.spdf.df_18.core$CV <- as.factor(ptm.var.spdf.df_18.core$CV)
ptm.var.spdf.df_18.core$CV <- factor(ptm.var.spdf.df_18.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))
ptm.var.spdf.df_20.core$CV <- as.factor(ptm.var.spdf.df_20.core$CV)
ptm.var.spdf.df_20.core$CV <- factor(ptm.var.spdf.df_20.core$CV, 
                                     levels=c("< 10%","11-20%","21-30%","31-40%","41-50%","51-60%","> 60%"))

# plot CV in bins
PTM_10_plot_bin_GS_var <- GSplotBin(ptm.var.spdf.df_10.core,"CV",survey.area.core,"Variance","CV")
PTM_11_plot_bin_GS_var <- GSplotBin(ptm.var.spdf.df_11.core,"CV",survey.area.core,"Variance","CV")
PTM_13_plot_bin_GS_var <- GSplotBin(ptm.var.spdf.df_13.core,"CV",survey.area.core,"Variance","CV")
PTM_14_plot_bin_GS_var <- GSplotBin(ptm.var.spdf.df_14.core,"CV",survey.area.core,"Variance","CV")
PTM_16_plot_bin_GS_var <- GSplotBin(ptm.var.spdf.df_16.core,"CV",survey.area.core,"Variance","CV")
PTM_18_plot_bin_GS_var <- GSplotBin(ptm.var.spdf.df_18.core,"CV",survey.area.core,"Variance","CV")
PTM_20_plot_bin_GS_var <- GSplotBin(ptm.var.spdf.df_20.core,"CV",survey.area.core,"Variance","CV")

# save 
saveplot(PTM_10_plot_bin_GS_var,"Results/PTM/Plots/variance/core_only/greyscale/bins/PTM_10_plot_bin_GS_var.png")
saveplot(PTM_11_plot_bin_GS_var,"Results/PTM/Plots/variance/core_only/greyscale/bins/PTM_11_plot_bin_GS_var.png")
saveplot(PTM_13_plot_bin_GS_var,"Results/PTM/Plots/variance/core_only/greyscale/bins/PTM_13_plot_bin_GS_var.png")
saveplot(PTM_14_plot_bin_GS_var,"Results/PTM/Plots/variance/core_only/greyscale/bins/PTM_14_plot_bin_GS_var.png")
saveplot(PTM_16_plot_bin_GS_var,"Results/PTM/Plots/variance/core_only/greyscale/bins/PTM_16_plot_bin_GS_var.png")
saveplot(PTM_18_plot_bin_GS_var,"Results/PTM/Plots/variance/core_only/greyscale/bins/PTM_18_plot_bin_GS_var.png")
saveplot(PTM_20_plot_bin_GS_var,"Results/PTM/Plots/variance/core_only/greyscale/bins/PTM_20_plot_bin_GS_var.png")



    # Continuous ####

# Load spatial dataframes
ptm.var.spdf.df_10.core<-read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_10.core.csv")
#ptm.var.spdf.df_11.core<-read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_11.core.csv")
#ptm.var.spdf.df_13.core<-read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_13.core.csv")
#ptm.var.spdf.df_14.core<-read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_14.core.csv")
#ptm.var.spdf.df_16.core<-read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_16.core.csv")
#ptm.var.spdf.df_18.core<-read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_18.core.csv")
ptm.var.spdf.df_20.core<-read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_20.core.csv")

# greyscale plots
PTM_varplot_final10.core.bw <- GSplotFun(ptm.var.spdf.df_10.core, survey.area.core, "variance", "2010")
PTM_varplot_final11.core.bw <- GSplotFun(ptm.var.spdf.df_11.core, survey.area.core, "variance")
PTM_varplot_final13.core.bw <- GSplotFun(ptm.var.spdf.df_13.core, survey.area.core, "variance")
PTM_varplot_final14.core.bw <- GSplotFun(ptm.var.spdf.df_14.core, survey.area.core, "variance")
PTM_varplot_final16.core.bw <- GSplotFun(ptm.var.spdf.df_16.core, survey.area.core, "variance")
PTM_varplot_final18.core.bw <- GSplotFun(ptm.var.spdf.df_18.core, survey.area.core, "variance")
PTM_varplot_final20.core.bw <- GSplotFun(ptm.var.spdf.df_20.core, survey.area.core, "variance", "2020")

# save greyscale
saveplot(PTM_varplot_final10.core.bw, "Results/PTM/Plots/variance/core_only/greyscale/2010_PTM_var.core.bw.png")
saveplot(PTM_varplot_final11.core.bw, "Results/PTM/Plots/variance/core_only/greyscale/2011_PTM_var.core.bw.png")
saveplot(PTM_varplot_final13.core.bw, "Results/PTM/Plots/variance/core_only/greyscale/2013_PTM_var.core.bw.png")
saveplot(PTM_varplot_final14.core.bw, "Results/PTM/Plots/variance/core_only/greyscale/2014_PTM_var.core.bw.png")
saveplot(PTM_varplot_final16.core.bw, "Results/PTM/Plots/variance/core_only/greyscale/2016_PTM_var.core.bw.png")
saveplot(PTM_varplot_final18.core.bw, "Results/PTM/Plots/variance/core_only/greyscale/2018_PTM_var.core.bw.png")
saveplot(PTM_varplot_final20.core.bw, "Results/PTM/Plots/variance/core_only/greyscale/2020_PTM_var.core.bw.png")

# colour plots
PTM_varplot_final10.core.col <- CLplotFun(ptm.var.spdf.df_10.core, survey.area.core, "variance")
PTM_varplot_final11.core.col <- CLplotFun(ptm.var.spdf.df_11.core, survey.area.core, "variance")
PTM_varplot_final13.core.col <- CLplotFun(ptm.var.spdf.df_13.core, survey.area.core, "variance")
PTM_varplot_final14.core.col <- CLplotFun(ptm.var.spdf.df_14.core, survey.area.core, "variance")
PTM_varplot_final16.core.col <- CLplotFun(ptm.var.spdf.df_16.core, survey.area.core, "variance")
PTM_varplot_final18.core.col <- CLplotFun(ptm.var.spdf.df_18.core, survey.area.core, "variance")
PTM_varplot_final20.core.col <- CLplotFun(ptm.var.spdf.df_20.core, survey.area.core, "variance")

# save colour
saveplot(PTM_varplot_final10.core.col, "Results/PTM/Plots/variance/core_only/colour/2010_PTM_var.core.col.png")
saveplot(PTM_varplot_final11.core.col, "Results/PTM/Plots/variance/core_only/colour/2011_PTM_var.core.col.png")
saveplot(PTM_varplot_final13.core.col, "Results/PTM/Plots/variance/core_only/colour/2013_PTM_var.core.col.png")
saveplot(PTM_varplot_final14.core.col, "Results/PTM/Plots/variance/core_only/colour/2014_PTM_var.core.col.png")
saveplot(PTM_varplot_final16.core.col, "Results/PTM/Plots/variance/core_only/colour/2016_PTM_var.core.col.png")
saveplot(PTM_varplot_final18.core.col, "Results/PTM/Plots/variance/core_only/colour/2018_PTM_var.core.col.png")
saveplot(PTM_varplot_final20.core.col, "Results/PTM/Plots/variance/core_only/colour/2020_PTM_var.core.col.png")


#### Stump-tailed macaque #####################################################
## Load data ####

# Observation data. Unique to species 
stm_obsdata <- read.csv("Species_Data/STM/R Data/obsdata.csv", header = TRUE)
stm_obsdata$object <- as.factor(stm_obsdata$object)
stm_obsdata$Sample.Label <- as.factor(stm_obsdata$Sample.Label)
str(stm_obsdata)
head(stm_obsdata)

# Transect data. Unique to species. Stratum not needed in DF
stm_distdata <- read.csv("Species_Data/STM/R Data/distdata.csv", header = TRUE)
stm_distdata$object <- as.factor(stm_distdata$object)
stm_distdata$NameObserver <- as.factor(stm_distdata$NameObserver)
stm_distdata$transect <- as.factor(stm_distdata$transect)
stm_distdata$year <- as.factor(stm_distdata$year)
stm_distdata$date <- as.Date(stm_distdata$date, format = "%d/%m/%Y")
str(stm_distdata)
head(stm_distdata)

# check for T20 obs
stm_distdata[stm_distdata$transect=="20",]
# none

## Plot the covariates across the grid with group sizes ####

## Warning - plots take a while to run

# habitat
plot_STMobs_habitat <- ggplot() + 
                    grid_plot_obj(preddata200$habitat, "habitat", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Habitat",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), 
                    data=stm_distdata, colour="red", alpha=I(0.7))+
                    gg.opts
ggsave("Plots/STM/plot_STMobs_habitat.png", plot_STMobs_habitat, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstWater
plot_STMobs_dstWater <- ggplot() + 
                     grid_plot_obj(preddata200$dstWater, "dstWater", pred.polys_200) + 
                     coord_equal()+
                     labs(fill="Distance to water",x="x",y="y",size="Group size")+
                     geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                     geom_point(aes(x, y, size=size), 
                     data=stm_distdata, colour="red", alpha=I(0.7))+
                     gg.opts

ggsave("Plots/STM/plot_STMobs_dstWater.png", plot_STMobs_dstWater, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStlmnt
plot_STMobs_dstStlmnt <- ggplot() + 
                    grid_plot_obj(preddata200$dstStlmnt, "dstStlmnt", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to settlement",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=stm_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts
ggsave("Plots/STM/plot_STMobs_dstStlmnt.png", plot_STMobs_dstStlmnt, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstRoad
plot_STMobs_dstRoad <- ggplot() + 
                    grid_plot_obj(preddata200$dstRoad, "dstRoad", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to road",x="x",y="y",size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=stm_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/STM/plot_STMobs_dstRoad.png", plot_STMobs_dstRoad, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstBorder
plot_STMobs_dstBorder <- ggplot() + 
                      grid_plot_obj(preddata200$dstBorder, "dstBorder", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to VN border",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=stm_distdata, 
                      colour="red", alpha=I(0.7))+
                      gg.opts

ggsave("Plots/STM/plot_STMobs_dstBorder.png", plot_STMobs_dstBorder, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstStation
plot_STMobs_dstStation <- ggplot() + 
                    grid_plot_obj(preddata200$dstStation, "dstStation", pred.polys_200) + 
                    coord_equal()+
                    labs(fill="Distance to ranger station",x="x",y="y",
                    size="Group size")+
                    geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                    geom_point(aes(x, y, size=size), data=stm_distdata, colour="red", 
                    alpha=I(0.7))+
                    gg.opts

ggsave("Plots/STM/plot_STMobs_dstStation.png", plot_STMobs_dstStation, width = 20, 
       height = 20, units = "cm", dpi = 300)

# dstELC
plot_STMobs_dstELC <- ggplot() + 
                      grid_plot_obj(preddata200$dstELC, "dstELC", pred.polys_200) + 
                      coord_equal()+
                      labs(fill="Distance to ELC",x="x",y="y",size="Group size")+
                      geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                      geom_point(aes(x, y, size=size), data=stm_distdata, colour="red", 
                      alpha=I(0.7))+
                      gg.opts

ggsave("Plots/STM/plot_STMobs_dstELC.png", plot_STMobs_dstELC, width = 20, 
       height = 20, units = "cm", dpi = 300)

# elevation
plot_STMobs_elev <- ggplot() + 
                  grid_plot_obj(preddata200$elevation, "elevation", pred.polys_200) + 
                  coord_equal()+
                  labs(fill="Elevation (m)",x="x",y="y",size="Group size")+
                  geom_line(aes(x, y, group=Transect.Label), data=segdata)+
                  geom_point(aes(x, y, size=size), data=stm_distdata, colour="red", 
                  alpha=I(0.7))+
                  gg.opts

ggsave("Plots/STM/plot_STMobs_elev.png", plot_STMobs_elev, width = 20, height = 20, 
       units = "cm", dpi = 300)


## Exploratory plots and linear models ####


## Checking that there are no large gaps in the range of variables for LTM observations. 

# subset segdata to get only the segments with LTM observations
stm_varcheck <- segdata[match(stm_obsdata$Sample.Label,segdata$Sample.Label), ]
par(mfrow=c(1,2))

# habitat 
plot(segdata$habitat)
plot(stm_varcheck$habitat)
# vast majority of observations are in dense forest (matches with known ecology (IUCN redlist) 

# dstWater 
hist(segdata$dstWater)
hist(stm_varcheck$dstWater)
# No gaps 

# dstStlmnt 
hist(segdata$dstStlmnt)
hist(stm_varcheck$dstStlmnt)
# Small gap between 8-10km   

# dstRoad 
hist(segdata$dstRoad)
hist(stm_varcheck$dstRoad)
# No obs beyond 2km - this variable will not be used 

# dstBorder 
hist(segdata$dstBorder)
hist(stm_varcheck$dstBorder)
# Gap between 10-15km, and no observations beyond 30km 

# dstStation 
hist(segdata$dstStation)
hist(stm_varcheck$dstStation)
# no observations beyond 11km 

# elevation 
hist(segdata$elevation)
hist(stm_varcheck$elevation)
# gap between 350-400m. Not a major problem 

## Histograms

# Distance 
stm_h1 <- ggplot(stm_distdata, aes(distance))+ geom_histogram(binwidth = 1)
stm_h2 <- ggplot(stm_distdata, aes(distance))+ geom_histogram(binwidth = 5)
stm_h3 <- ggplot(stm_distdata, aes(distance))+ geom_histogram(binwidth = 10)
stm_h4 <- ggplot(stm_distdata, aes(distance))+ geom_histogram(binwidth = 15)
stm_h5 <- ggplot(stm_distdata, aes(distance))+ geom_histogram(binwidth = 20)
stm_h6 <- ggplot(stm_distdata, aes(distance))+ geom_histogram(binwidth = 40)
plot_grid(stm_h1,stm_h2,stm_h3,stm_h4,stm_h5,stm_h6)
#  Not a horrible distance histogram.  There is potentially some minor evidence of clumping - spikes in observations at 0, 6, 11, 21, 23m.  Slightly odd distances for there to be clumps though, so I am not convinced this is a problem with distance estimation by the teams. If they were rounding up or down they would tend to do it to the nearest 5 or 10m.  There are only 32 observations so this pattern is more likely just an artefact of few observations i.e. it would probably dissapear with more observations

# cluster size, observer, habitat, year, month, transect
stm_h7 <- ggplot(stm_distdata, aes(size))+geom_histogram(binwidth = 0.5)
stm_h8 <- ggplot(stm_distdata, aes(NameObserver))+geom_histogram(stat="count")
stm_h9 <- ggplot(stm_distdata, aes(habitat))+geom_histogram(stat="count")
stm_h10 <- ggplot(stm_distdata, aes(year))+geom_histogram(stat="count")
stm_h11 <- ggplot(stm_distdata, aes(month))+geom_histogram(stat="count")
stm_h12 <- ggplot(stm_distdata, aes(transect))+geom_histogram(stat="count")
plot_grid(stm_h7,stm_h8,stm_h9,stm_h10,stm_h11,stm_h12)
#  vast majority of observations have cluster sizes of <5.  These will be underetimates.  Major drop in observations in 2013, which also occured for LTM (but not PTM). I am not sure why this is.

## Plots of distance against variables

plotlabs <- function(title,x,y) {
  
  title = title
  xlab = x
  ylab = y
  
  list(labs(x = x, y=y, title=title))
}

stm_d1 <- ggplot(stm_distdata, aes(x=habitat, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by habitat","Habitat","Distance (m)")
stm_d2 <- ggplot(stm_distdata, aes(x=NameObserver, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by observer","Observer","Distance (m)")
stm_d3 <- ggplot(stm_distdata, aes(x=month, y=distance))+geom_boxplot()+ 
      plotlabs("Distance by month","Month","Distance (m)")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
stm_d4 <- ggplot(stm_distdata, aes(x=size, y=distance))+geom_point()+
      plotlabs("Distance by size","Group size","Distance (m)")
stm_d5 <- ggplot(stm_distdata, aes(x=transect, y=distance))+geom_point()+
      plotlabs("Distance by transect","Transect","Distance (m)")
plot_grid(stm_d1,stm_d2,stm_d3,stm_d4,stm_d5)
#  distance is larger in open forest but too few observations in open forest to take anything away. Too few observations to get much from these plots


## Plots of cluster size against variables
stm_s1 <- ggplot(stm_distdata, aes(x=habitat, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by habitat","Habitat","Group size")
stm_s2 <- ggplot(stm_distdata, aes(x=NameObserver, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by observer","observer","Group size")
stm_s3 <- ggplot(stm_distdata, aes(x=month, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by month","month","Group size")+
      scale_x_discrete(limits=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul"))
stm_s4 <- ggplot(stm_distdata, aes(x=year, y=size))+geom_boxplot()+ 
      plotlabs("Grp size by year","year","Group size")
stm_s5 <- ggplot(stm_distdata, aes(x=as.factor(transect), y=size))+geom_boxplot()+ 
      plotlabs("Grp size by transect","transect","Group size")
plot_grid(stm_s1,stm_s2,stm_s3,stm_s4,stm_s5)
# There seems to be a general decrease in group size over the course of the stud period, but hard to say anythiing really with only 32 observations. 

## Linear models

par(mfrow=c(1,1))
# group size ~ distance
newdist <- data.frame(distance=seq(0,100,len=10))
lm1 <- lm(size~distance, data=stm_distdata)
plot(stm_distdata$size~stm_distdata$distance)
lines(newdist$distance, as.vector(predict(lm1,newdist)))
summary(lm1)
# There appears to be a weak positive relationship, but with small effect size and not significant to 0.05.  Also this relationship is being driven by one outlier at around 47m.   

## Estimating the detection function ####

### the DF model from CDS is HR with no adjustment.
stmDF.hr <- ds(stm_distdata, truncation=50, key="hr")


 
## Fitting a spatial model ####

## The best DF is stmDF.hn

# I am setting group = TRUE which means abundance of groups rather than individuals will be estimated.

# We need to define segment.area = "Sample.Area"

# Use method=REML

# Need to test quasipoisson, tweedie, negative binomial distributions

# Need to test for autocorrelation. If present add a covariance structure to the model

# Need to remove observations from 'obsdata' that have a distance greater than the truncation distance used in the detection function (50m)
stm_obsdata <- stm_obsdata %>% filter(distance <= 50)


  ## Quasipoisson response ####

# saturated model
stmDSM.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+
                         s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                         habitat,
                  stmDF.hr, segdata, stm_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.sat)
par(mfrow=c(2,3))
plot(stmDSM.sat, scale = 0)
gam.check(stmDSM.sat)
# DE = 35.3. Ran with warning. Nothing significant, not even habitat!

# increase k for all
stmDSM.sat2 <- dsm(Nhat ~ s(dstWater,bs="ts",k=15)+ s(dstStlmnt,bs="ts",k=15)+ 
                          s(elevation,bs="ts",k=15)+ s(dstBorder,bs="ts",k=15)+ 
                          s(dstStation,bs="ts",k=15)+  
                          habitat,
                  stmDF.hr, segdata, stm_obsdata, method = "REML",
                  family = quasipoisson(link = "log"), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.sat2)
par(mfrow=c(2,3))
plot(stmDSM.sat2, scale = 0)
gam.check(stmDSM.sat2)
# DE = 46.9. Ran with multiple warnings. Nowt of signifcance

# Remove dstWater and reduce k
stmDSM.qp.3 <- dsm(Nhat ~ s(dstStlmnt,bs="ts",k=10)+ s(elevation,bs="ts",k=10)+ 
                          s(dstBorder,bs="ts",k=10)+ s(dstStation,bs="ts",k=10)+  
                          habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.qp.3)
par(mfrow=c(2,2))
plot(stmDSM.qp.3, scale = 0)
gam.check(stmDSM.qp.3)
# DE = 33. ran with warning. Nothing sig

# remove dstSltmnt and remove k setting
stmDSM.qp.4 <- dsm(Nhat ~  s(elevation,bs="ts")+ s(dstBorder,bs="ts")+ 
                           s(dstStation,bs="ts")+  
                           habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.qp.4)
par(mfrow=c(2,2))
plot(stmDSM.qp.4, scale = 0)
gam.check(stmDSM.qp.4)
# DE = 32. Ran with warnings. dstBorder close to sig. 

# remove dstStation
stmDSM.qp.5 <- dsm(Nhat ~ s(elevation,bs="ts")+ s(dstBorder,bs="ts")+ 
                          habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.qp.5)
par(mfrow=c(2,2))
plot(stmDSM.qp.5, scale = 0)
gam.check(stmDSM.qp.5)
# DE = 27.9. Ran with warnings. Really surprised habitat isn't significant

# remove dstBorder
stmDSM.qp.6 <- dsm(Nhat ~ s(elevation,bs="ts")+
                          habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = quasipoisson(link = "log"), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.qp.6)
par(mfrow=c(1,1))
plot(stmDSM.qp.6, scale = 0)
gam.check(stmDSM.qp.6)
# DE = 15.5. Ran without warnings. elevation and habiatat (D,O) sig.


### The best QP model is stmDSM.qp.6


  ## Tweedie Response ####

# saturated model
stmDSM.tw.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+ 
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                            habitat,
                     stmDF.hr, segdata, stm_obsdata, method = "REML",
                  family = tw(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.tw.sat)
par(mfrow=c(2,3))
plot(stmDSM.tw.sat, scale = 0)
gam.check(stmDSM.tw.sat)
# AIC = 13965, DE = 12.2. Some habitat and dstBorder sig.

# remove dstStation
stmDSM.tw.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                          s(elevation,bs="ts")+ s(dstBorder,bs="ts")+   
                          habitat,
                     stmDF.hr, segdata, stm_obsdata, method = "REML",
                     family = tw(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.tw.2)
par(mfrow=c(2,2))
plot(stmDSM.tw.2, scale = 0)
gam.check(stmDSM.tw.2)
# AIC = 13965, DE = 12.2. still only dstBorder. elevation looks the worst despite it being sig in the QP models

# remove elevation
stmDSM.tw.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ 
                          s(dstBorder,bs="ts")+   
                          habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.tw.3)
par(mfrow=c(2,2))
plot(stmDSM.tw.3, scale = 0)
gam.check(stmDSM.tw.3)
# AIC = 13965, DE = 12.2.

# remove dstWater
stmDSM.tw.4 <- dsm(Nhat ~ s(dstStlmnt,bs="ts")+ s(dstBorder,bs="ts")+   
                          habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.tw.4)
par(mfrow=c(2,2))
plot(stmDSM.tw.4, scale = 0)
gam.check(stmDSM.tw.4)
# AIC = 13964, DE = 12.

# remove dstStlmnt
stmDSM.tw.5 <- dsm(Nhat ~ s(dstBorder,bs="ts")+ habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = tw(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.tw.5)
par(mfrow=c(1,1))
plot(stmDSM.tw.5, scale = 0)
gam.check(stmDSM.tw.5)
# AIC = 13963, DE = 11.8. Perhaps a GLM would be worth exploring as dstBorder is only term and has EDF of 0.92


### stmDSM.tw.5 is the best TW model


  ## Negative binomial response ####

# saturated model
stmDSM.nb.sat <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+ s(elevation,bs="ts")+  
                            s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                            habitat,
                  stmDF.hr, segdata, stm_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.nb.sat)
par(mfrow=c(2,3))
plot(stmDSM.nb.sat, scale = 0)
gam.check(stmDSM.nb.sat)
# AIC = 533, DE = 30.8. dstborder (and some hab) only sig term.

# remove elevation
stmDSM.nb.2 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+  
                          s(dstBorder,bs="ts")+ s(dstStation,bs="ts")+  
                          habitat,
                     stmDF.hr, segdata, stm_obsdata, method = "REML",
                     family = nb(), engine = "gam",
                     segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.nb.2)
par(mfrow=c(2,2))
plot(stmDSM.nb.2, scale = 0)
gam.check(stmDSM.nb.2)
# AIC = 533, DE = 30.8. No change

# remove dstStation
stmDSM.nb.3 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstStlmnt,bs="ts")+  
                          s(dstBorder,bs="ts")+ +  
                          habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.nb.3)
par(mfrow=c(2,2))
plot(stmDSM.nb.3, scale = 0)
gam.check(stmDSM.nb.3)
# AIC = 533, DE = 30.8.

# remove dstStlmnt
stmDSM.nb.4 <- dsm(Nhat ~ s(dstWater,bs="ts")+ s(dstBorder,bs="ts")+
                          habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.nb.4)
par(mfrow=c(1,2))
plot(stmDSM.nb.4, scale = 0)
gam.check(stmDSM.nb.4)
# AIC = 533, DE = 30.8.

# remove dstWater
stmDSM.nb.5 <- dsm(Nhat ~ s(dstBorder,bs="ts")+ 
                          habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   family = nb(), engine = "gam",
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.nb.5)
par(mfrow=c(1,1))
plot(stmDSM.nb.5, scale = 0)
gam.check(stmDSM.nb.5)
# AIC = 533, DE = 30.2. AIC and DE have barely changed from the saturated model suggesting the other covars contribute very little. I think testing a GLM is a good idea


### stmDSM.nb.5 is the best NB model



  ## GLM ####

# try all terms as linear terms in a GLM
stmDSM.glm.1 <- dsm(Nhat ~ habitat+dstBorder+dstStlmnt+dstStation+dstWater+elevation,
                    stmDF.hr, segdata, stm_obsdata, method = "REML",
                  family = poisson(link='log'), engine = "glm",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.glm.1)
# dstBorder sig. dstStlmnt & dstWater close.

# remove dstStation
stmDSM.glm.2 <- dsm(Nhat ~ habitat+dstBorder+dstStlmnt+dstWater+elevation,
                    stmDF.hr, segdata, stm_obsdata, method = "REML",
                  family = poisson(link='log'), engine = "glm",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.glm.2)
# As above

# remove elevation
stmDSM.glm.3 <- dsm(Nhat ~ habitat+dstBorder+dstStlmnt+dstWater,
                    stmDF.hr, segdata, stm_obsdata, method = "REML",
                    family = poisson(link='log'), engine = "glm",
                    segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.glm.3)
# dstStlmnt & dtWater still quite weak

# remove dstStlmnt
stmDSM.glm.4 <- dsm(Nhat ~ habitat+dstBorder+dstWater,
                    stmDF.hr, segdata, stm_obsdata, method = "REML",
                    family = poisson(link='log'), engine = "glm",
                    segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.glm.4)

# I am tempted to predict using GLM.3 as it has more terms than any other model. I would be interested to see what the predictions look like, becuase the predictions just using habitat and dstBorder don't seem particulalry realistic to me.


## Model selection ####

summary(stmDSM.qp.6) # DE = 5.5
summary(stmDSM.tw.5) # DE = 11.8
summary(stmDSM.nb.5) # DE = 30.2
summary(stmDSM.glm.3) # includes the most terms

stmDSM.tw.5$aic # 13963
stmDSM.nb.5$aic # 533
# NB model is better than TW

anova(stmDSM.nb.5,stmDSM.qp.6,test="Chisq")
# NB model has a lot less residual deviance

par(mfrow=c(2,2))
gam.check(stmDSM.qp.6)
gam.check(stmDSM.nb.5)
# Neither Q-Q plot is particlarly nice but NB model is much better

# compate nb model to glm
anova(stmDSM.nb.5,stmDSM.glm.3,test="Chisq")

### stmDSM.nb.5 is the best STM model

## Autocorrelation ####

par(mfrow=c(1,1))
dsm.cor(stmDSM.nb.5, max.lag=15, Segment.Label="Sample.Label")
# there is some autocorrelation at lag 1 

# Try autoregressive structure
stmDSM.nb.5.cor <- dsm(Nhat ~ s(dstBorder,bs="ts")+ habitat,
                   stmDF.hr, segdata, stm_obsdata, method = "REML",
                   engine = "gamm", correlation=corAR1(form=~sg.id|tr.id), niterPQL=50,
                   segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.nb.5.cor$lme)
summary(stmDSM.nb.5.cor$gam)
# No convergence (too few observations?)

# I will re-model nb.5 but include univariate smooths of x and y
stmDSM.nb.5.xy <- dsm(Nhat ~ s(dstBorder,bs="ts")+ s(x,bs="ts")+ s(y,bs="ts")+
                             habitat,
                  stmDF.hr, segdata, stm_obsdata, method = "REML",
                  family = nb(), engine = "gam",
                  segment.area = segdata$Sample.Area, group = TRUE)
summary(stmDSM.nb.5.xy)
dev.off()
par(mfrow=c(2,2))
plot(stmDSM.nb.5.xy, scale = 0)
gam.check(stmDSM.nb.5.xy)
# AIC = 533, DE = 30.8. 

par(mfrow=c(1,1))
dsm.cor(stmDSM.nb.5.xy, max.lag=15, Segment.Label="Sample.Label")
# The autocor is still there at lag 1, but it has been reduced. Enither of the smooth term are that significna tthough, so I am unsure whether the autocor is being dealt with properly. This is just not a very good set of models as there are too few data points.  I'm not sure these models will make it into the paper.

### stmDSM.nb.5.xy is the final STM model

## Abundance estimation ####
 

# Predict over 2010 habitat 
stm.global.pred10.core <- predict(stmDSM.nb.5.xy, preddata10_core, off.set = 40000)
write.csv(stm.global.pred10.core, file="Results/STM/core_only/stm.pred10.core.csv")

# Predict over 2011 habitat 
stm.global.pred11.core <- predict(stmDSM.nb.5.xy, preddata11_core, off.set = 40000)
write.csv(stm.global.pred11.core, file="Results/STM/core_only/stm.pred11.core.csv")

# Predict over 2013 habitat 
stm.global.pred13.core <- predict(stmDSM.nb.5.xy, preddata13_core, off.set = 40000)
write.csv(stm.global.pred13.core, file="Results/STM/core_only/stm.pred13.core.csv")

# Predict over 2014 habitat 
stm.global.pred14.core <- predict(stmDSM.nb.5.xy, preddata14_core, off.set = 40000)
write.csv(stm.global.pred14.core, file="Results/STM/core_only/stm.pred14.core.csv")

# Predict over 2016 habitat 
stm.global.pred16.core <- predict(stmDSM.nb.5.xy, preddata16_core, off.set = 40000)
write.csv(stm.global.pred16.core, file="Results/STM/core_only/stm.pred16.core.csv")

# Predict over 2018 habitat 
stm.global.pred18.core <- predict(stmDSM.nb.5.xy, preddata18_core, off.set = 40000)
write.csv(stm.global.pred18.core, file="Results/STM/core_only/stm.pred18.core.csv")

# Predict over 2020 habitat 
stm.global.pred20.core <- predict(stmDSM.nb.5.xy, preddata20_core, off.set = 40000)
write.csv(stm.global.pred20.core, file="Results/STM/core_only/stm.pred20.core.csv")



# Create new dataframes for plotting
stm.df.Final10.core <- data.frame(id = 1:47801,
                         abundance = stm.global.pred10.core)

stm.df.Final11.core <- data.frame(id = 1:47801,
                         abundance = stm.global.pred11.core)

stm.df.Final13.core <- data.frame(id = 1:47801,
                         abundance = stm.global.pred13.core)

stm.df.Final14.core <- data.frame(id = 1:47801,
                         abundance = stm.global.pred14.core)

stm.df.Final16.core <- data.frame(id = 1:47801,
                         abundance = stm.global.pred16.core)

stm.df.Final18.core <- data.frame(id = 1:47801,
                         abundance = stm.global.pred18.core)

stm.df.Final20.core <- data.frame(id = 1:47801,
                         abundance = stm.global.pred20.core)


## This creates a dataframe that can be plotted as a map
stm.spdf.df_10.core <- abunPlotDF(stm.df.Final10.core, pred.polys_200)
stm.spdf.df_11.core <- abunPlotDF(stm.df.Final11.core, pred.polys_200)
stm.spdf.df_13.core <- abunPlotDF(stm.df.Final13.core, pred.polys_200)
stm.spdf.df_14.core <- abunPlotDF(stm.df.Final14.core, pred.polys_200)
stm.spdf.df_16.core <- abunPlotDF(stm.df.Final16.core, pred.polys_200)
stm.spdf.df_18.core <- abunPlotDF(stm.df.Final18.core, pred.polys_200)
stm.spdf.df_20.core <- abunPlotDF(stm.df.Final20.core, pred.polys_200)

# save SPDFs
write.csv(stm.spdf.df_10.core,file="Results/STM/Plots/core_only/spdf/stm.spdf.df_10.core.csv")
write.csv(stm.spdf.df_11.core,file="Results/STM/Plots/core_only/spdf/stm.spdf.df_11.core.csv")
write.csv(stm.spdf.df_13.core,file="Results/STM/Plots/core_only/spdf/stm.spdf.df_13.core.csv")
write.csv(stm.spdf.df_14.core,file="Results/STM/Plots/core_only/spdf/stm.spdf.df_14.core.csv")
write.csv(stm.spdf.df_16.core,file="Results/STM/Plots/core_only/spdf/stm.spdf.df_16.core.csv")
write.csv(stm.spdf.df_18.core,file="Results/STM/Plots/core_only/spdf/stm.spdf.df_18.core.csv")
write.csv(stm.spdf.df_20.core,file="Results/STM/Plots/core_only/spdf/stm.spdf.df_20.core.csv")



    ## Plotting ####

# Load spatial dataframes
stm.spdf.df_10.core <- read.csv("Results/STM/Plots/core_only/spdf/stm.spdf.df_10.core.csv")
#stm.spdf.df_11.core <- read.csv("Results/STM/Plots/core_only/spdf/stm.spdf.df_11.core.csv") 
#stm.spdf.df_13.core <- read.csv("Results/STM/Plots/core_only/spdf/stm.spdf.df_13.core.csv") 
#stm.spdf.df_14.core <- read.csv("Results/STM/Plots/core_only/spdf/stm.spdf.df_14.core.csv") 
#stm.spdf.df_16.core <- read.csv("Results/STM/Plots/core_only/spdf/stm.spdf.df_16.core.csv") 
#stm.spdf.df_18.core <- read.csv("Results/STM/Plots/core_only/spdf/stm.spdf.df_18.core.csv")
stm.spdf.df_20.core <- read.csv("Results/STM/Plots/core_only/spdf/stm.spdf.df_20.core.csv")

# greyscale plots
STM_plot_10_core_gr <- GSplotFun(stm.spdf.df_10.core, survey.area.core, "abundance", "2010")
STM_plot_11_core_gr <- GSplotFun(stm.spdf.df_11.core, survey.area.core, "abundance")
STM_plot_13_core_gr <- GSplotFun(stm.spdf.df_13.core, survey.area.core, "abundance")
STM_plot_14_core_gr <- GSplotFun(stm.spdf.df_14.core, survey.area.core, "abundance")
STM_plot_16_core_gr <- GSplotFun(stm.spdf.df_16.core, survey.area.core, "abundance")
STM_plot_18_core_gr <- GSplotFun(stm.spdf.df_18.core, survey.area.core, "abundance")
STM_plot_20_core_gr <- GSplotFun(stm.spdf.df_20.core, survey.area.core, "abundance", "2020")

# save greyscale
saveplot(STM_plot_10_core_gr,"Results/STM/Plots/core_only/greyscale/STM_plot_10_core_gr.png")
saveplot(STM_plot_11_core_gr,"Results/STM/Plots/core_only/greyscale/STM_plot_11_core_gr.png")
saveplot(STM_plot_13_core_gr,"Results/STM/Plots/core_only/greyscale/STM_plot_13_core_gr.png")
saveplot(STM_plot_14_core_gr,"Results/STM/Plots/core_only/greyscale/STM_plot_14_core_gr.png")
saveplot(STM_plot_16_core_gr,"Results/STM/Plots/core_only/greyscale/STM_plot_16_core_gr.png")
saveplot(STM_plot_18_core_gr,"Results/STM/Plots/core_only/greyscale/STM_plot_18_core_gr.png")
saveplot(STM_plot_20_core_gr,"Results/STM/Plots/core_only/greyscale/STM_plot_20_core_gr.png")

# colour plots
STM_plot_10_core_col <- CLplotFun(stm.spdf.df_10.core, survey.area.core, "abundance")
STM_plot_11_core_col <- CLplotFun(stm.spdf.df_11.core, survey.area.core, "abundance")
STM_plot_13_core_col <- CLplotFun(stm.spdf.df_13.core, survey.area.core, "abundance")
STM_plot_14_core_col <- CLplotFun(stm.spdf.df_14.core, survey.area.core, "abundance")
STM_plot_16_core_col <- CLplotFun(stm.spdf.df_16.core, survey.area.core, "abundance")
STM_plot_18_core_col <- CLplotFun(stm.spdf.df_18.core, survey.area.core, "abundance")
STM_plot_20_core_col <- CLplotFun(stm.spdf.df_20.core, survey.area.core, "abundance")

# save colour
saveplot(STM_plot_10_core_col,"Results/STM/Plots/core_only/colour/STM_plot_10_core_col.png")
saveplot(STM_plot_11_core_col,"Results/STM/Plots/core_only/colour/STM_plot_11_core_col.png")
saveplot(STM_plot_13_core_col,"Results/STM/Plots/core_only/colour/STM_plot_13_core_col.png")
saveplot(STM_plot_14_core_col,"Results/STM/Plots/core_only/colour/STM_plot_14_core_col.png")
saveplot(STM_plot_16_core_col,"Results/STM/Plots/core_only/colour/STM_plot_16_core_col.png")
saveplot(STM_plot_18_core_col,"Results/STM/Plots/core_only/colour/STM_plot_18_core_col.png")
saveplot(STM_plot_20_core_col,"Results/STM/Plots/core_only/colour/STM_plot_20_core_col.png")


## plot grids (abundance and variance)

# greyscale 
stm_2yrs_gs <- 
  STM_plot_10_core_gr + STM_varplot_final10.core.bw  +
  STM_plot_20_core_gr + STM_varplot_final20.core.bw

# remove x axis labels and text for plots 1 and 2
stm_2yrs_gs[[1]] <- stm_2yrs_gs[[1]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())
stm_2yrs_gs[[2]] <- stm_2yrs_gs[[2]] + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank())

# remove y axis labels and text for plots 2 and 4
stm_2yrs_gs[[2]] <- stm_2yrs_gs[[2]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())
stm_2yrs_gs[[4]] <- stm_2yrs_gs[[4]] + theme(axis.title.y = element_blank()) + 
  theme(axis.text.y = element_blank())

# save
saveplot(stm_2yrs_gs, "Results/STM/Plots/core_only/greyscale/plot_grids/stm_2yrs_gs.png")


## Variance estimation ####

# estimate variance
stm.var.Final10.core <- varEstfun(preddata10_core, stmDSM.nb.5.xy)
stm.var.Final11.core <- varEstfun(preddata11_core, stmDSM.nb.5.xy)
stm.var.Final13.core <- varEstfun(preddata13_core, stmDSM.nb.5.xy)
stm.var.Final14.core <- varEstfun(preddata14_core, stmDSM.nb.5.xy)
stm.var.Final16.core <- varEstfun(preddata16_core, stmDSM.nb.5.xy)
stm.var.Final18.core <- varEstfun(preddata18_core, stmDSM.nb.5.xy)
stm.var.Final20.core <- varEstfun(preddata20_core, stmDSM.nb.5.xy)

# save variance estimates
write.csv(stm.var.Final10.core, file="Results/STM/core_only/stm.var10.core.csv")
write.csv(stm.var.Final11.core, file="Results/STM/core_only/stm.var11.core.csv")
write.csv(stm.var.Final13.core, file="Results/STM/core_only/stm.var13.core.csv")
write.csv(stm.var.Final14.core, file="Results/STM/core_only/stm.var14.core.csv")
write.csv(stm.var.Final16.core, file="Results/STM/core_only/stm.var16.core.csv")
write.csv(stm.var.Final18.core, file="Results/STM/core_only/stm.var18.core.csv")
write.csv(stm.var.Final20.core, file="Results/STM/core_only/stm.var20.core.csv")

# create spdf's for plotting
stm.var.spdf.df_10.core <- varPlotDF(stm.var.Final10.core, pred.polys_200)
stm.var.spdf.df_11.core <- varPlotDF(stm.var.Final11.core, pred.polys_200)
stm.var.spdf.df_13.core <- varPlotDF(stm.var.Final13.core, pred.polys_200)
stm.var.spdf.df_14.core <- varPlotDF(stm.var.Final14.core, pred.polys_200)
stm.var.spdf.df_16.core <- varPlotDF(stm.var.Final16.core, pred.polys_200)
stm.var.spdf.df_18.core <- varPlotDF(stm.var.Final18.core, pred.polys_200)
stm.var.spdf.df_20.core <- varPlotDF(stm.var.Final20.core, pred.polys_200)

# save spdf's
write.csv(stm.var.spdf.df_10.core,
          file="Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_10.core.csv")
write.csv(stm.var.spdf.df_11.core,
          file="Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_11.core.csv")
write.csv(stm.var.spdf.df_13.core,
          file="Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_13.core.csv")
write.csv(stm.var.spdf.df_14.core,
          file="Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_14.core.csv")
write.csv(stm.var.spdf.df_16.core,
          file="Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_16.core.csv")
write.csv(stm.var.spdf.df_18.core,
          file="Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_18.core.csv")
write.csv(stm.var.spdf.df_20.core,
          file="Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_20.core.csv") 


    ## Plotting variance ####

# Load spatial dataframes
stm.var.spdf.df_10.core<-read.csv("Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_10.core.csv")
#stm.var.spdf.df_11.core<-read.csv("Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_11.core.csv")
#stm.var.spdf.df_13.core<-read.csv("Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_13.core.csv")
#stm.var.spdf.df_14.core<-read.csv("Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_14.core.csv")
#stm.var.spdf.df_16.core<-read.csv("Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_16.core.csv")
#stm.var.spdf.df_18.core<-read.csv("Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_18.core.csv")
stm.var.spdf.df_20.core<-read.csv("Results/STM/Plots/variance/core_only/spdf/stm.var.spdf.df_20.core.csv")

# greyscale plots
STM_varplot_final10.core.bw <- GSplotFun(stm.var.spdf.df_10.core, survey.area.core, "variance","2010")
STM_varplot_final11.core.bw <- GSplotFun(stm.var.spdf.df_11.core, survey.area.core, "variance")
STM_varplot_final13.core.bw <- GSplotFun(stm.var.spdf.df_13.core, survey.area.core, "variance")
STM_varplot_final14.core.bw <- GSplotFun(stm.var.spdf.df_14.core, survey.area.core, "variance")
STM_varplot_final16.core.bw <- GSplotFun(stm.var.spdf.df_16.core, survey.area.core, "variance")
STM_varplot_final18.core.bw <- GSplotFun(stm.var.spdf.df_18.core, survey.area.core, "variance")
STM_varplot_final20.core.bw <- GSplotFun(stm.var.spdf.df_20.core, survey.area.core, "variance", "2020")

# save greyscale
saveplot(STM_varplot_final10.core.bw, "Results/STM/Plots/variance/core_only/greyscale/2010_STM_var.core.bw.png")
saveplot(STM_varplot_final11.core.bw, "Results/STM/Plots/variance/core_only/greyscale/2011_STM_var.core.bw.png")
saveplot(STM_varplot_final13.core.bw, "Results/STM/Plots/variance/core_only/greyscale/2013_STM_var.core.bw.png")
saveplot(STM_varplot_final14.core.bw, "Results/STM/Plots/variance/core_only/greyscale/2014_STM_var.core.bw.png")
saveplot(STM_varplot_final16.core.bw, "Results/STM/Plots/variance/core_only/greyscale/2016_STM_var.core.bw.png")
saveplot(STM_varplot_final18.core.bw, "Results/STM/Plots/variance/core_only/greyscale/2018_STM_var.core.bw.png")
saveplot(STM_varplot_final20.core.bw, "Results/STM/Plots/variance/core_only/greyscale/2020_STM_var.core.bw.png")

# colour plots
STM_varplot_final10.core.col <- CLplotFun(stm.var.spdf.df_10.core, survey.area.core, "variance")
STM_varplot_final11.core.col <- CLplotFun(stm.var.spdf.df_11.core, survey.area.core, "variance")
STM_varplot_final13.core.col <- CLplotFun(stm.var.spdf.df_13.core, survey.area.core, "variance")
STM_varplot_final14.core.col <- CLplotFun(stm.var.spdf.df_14.core, survey.area.core, "variance")
STM_varplot_final16.core.col <- CLplotFun(stm.var.spdf.df_16.core, survey.area.core, "variance")
STM_varplot_final18.core.col <- CLplotFun(stm.var.spdf.df_18.core, survey.area.core, "variance")
STM_varplot_final20.core.col <- CLplotFun(stm.var.spdf.df_20.core, survey.area.core, "variance")

# save colour
saveplot(STM_varplot_final10.core.col, "Results/STM/Plots/variance/core_only/colour/2010_STM_var.core.col.png")
saveplot(STM_varplot_final11.core.col, "Results/STM/Plots/variance/core_only/colour/2011_STM_var.core.col.png")
saveplot(STM_varplot_final13.core.col, "Results/STM/Plots/variance/core_only/colour/2013_STM_var.core.col.png")
saveplot(STM_varplot_final14.core.col, "Results/STM/Plots/variance/core_only/colour/2014_STM_var.core.col.png")
saveplot(STM_varplot_final16.core.col, "Results/STM/Plots/variance/core_only/colour/2016_STM_var.core.col.png")
saveplot(STM_varplot_final18.core.col, "Results/STM/Plots/variance/core_only/colour/2018_STM_var.core.col.png")
saveplot(STM_varplot_final20.core.col, "Results/STM/Plots/variance/core_only/colour/2020_STM_var.core.col.png")


#### Plotting all species #####################################################
  ## separate plots for arboreal and ground species, with CV plots ####

# YCG, BSD, GSL, 

# RMJ, GPF, LTM, PTM 

arb2 <- YCG_20_plot_bin_GS + YCG_20_plot_bin_GS_var + 
        BSD_20_plot_bin_GS + BSD_20_plot_bin_GS_var +
        GSL_20_plot_bin_GS + GSL_20_plot_bin_GS_var +
        plot_layout(ncol=2)

# remove y axis labels from all plots
arb2[[1]] <- arb2[[1]] + theme(axis.title.y = element_blank())
arb2[[2]] <- arb2[[2]] + theme(axis.title.y = element_blank())
arb2[[3]] <- arb2[[3]] + theme(axis.title.y = element_blank())
arb2[[4]] <- arb2[[4]] + theme(axis.title.y = element_blank())
arb2[[5]] <- arb2[[5]] + theme(axis.title.y = element_blank())
arb2[[6]] <- arb2[[6]] + theme(axis.title.y = element_blank())

# remove y axis numbers and all ticks from all plots
arb2[[1]] <- arb2[[1]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) 
arb2[[2]] <- arb2[[2]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
arb2[[3]] <- arb2[[3]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
arb2[[4]] <- arb2[[4]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
arb2[[5]] <- arb2[[5]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
arb2[[6]] <- arb2[[6]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())

# remove x axis labels from all plots
arb2[[1]] <- arb2[[1]] + theme(axis.title.x = element_blank())
arb2[[2]] <- arb2[[2]] + theme(axis.title.x = element_blank())
arb2[[3]] <- arb2[[3]] + theme(axis.title.x = element_blank())
arb2[[4]] <- arb2[[4]] + theme(axis.title.x = element_blank())
arb2[[5]] <- arb2[[5]] + theme(axis.title.x = element_blank())
arb2[[6]] <- arb2[[6]] + theme(axis.title.x = element_blank())

# remove x axis numbers and ticks from all plots
arb2[[1]] <- arb2[[1]] + theme(axis.text.x = element_blank())
arb2[[2]] <- arb2[[2]] + theme(axis.text.x = element_blank())
arb2[[3]] <- arb2[[3]] + theme(axis.text.x = element_blank())
arb2[[4]] <- arb2[[4]] + theme(axis.text.x = element_blank())
arb2[[5]] <- arb2[[5]] + theme(axis.text.x = element_blank())
arb2[[6]] <- arb2[[6]] + theme(axis.text.x = element_blank())

# remove legend title for plots 3,4,5,6
#arb2[[3]] <- arb2[[3]] + theme(legend.title = element_blank())
#arb2[[4]] <- arb2[[4]] + theme(legend.title = element_blank())
#arb2[[5]] <- arb2[[5]] + theme(legend.title = element_blank())
#arb2[[6]] <- arb2[[6]] + theme(legend.title = element_blank())

# remove legends for plots 1,2,5,6
arb2[[1]] <- arb2[[1]] + theme(legend.position = "none")
arb2[[2]] <- arb2[[2]] + theme(legend.position = "none")
arb2[[5]] <- arb2[[5]] + theme(legend.position = "none")
arb2[[6]] <- arb2[[6]] + theme(legend.position = "none")

# change legend position for plot 3
arb2[[3]] <- arb2[[3]] + theme(legend.position = "left")

# increase the size of the legend
arb2[[3]] <- arb2[[3]] + theme(legend.key.size = unit(1.0,"cm"), legend.key.width = unit(1.0,"cm"), 
                               legend.title = element_text(size = 15))
arb2[[4]] <- arb2[[4]] + theme(legend.key.size = unit(1.0,"cm"), legend.key.width = unit(1.0,"cm"),
                               legend.title = element_text(size = 15))

# add titles for plots 1,3,5 ad align right
#arb2[[1]] <- arb2[[1]] + ggtitle("Yellow-cheeked crested gibbon") + theme(plot.title = element_text(hjust=1))
#arb2[[3]] <- arb2[[3]] + ggtitle("Black-shanked douc") + theme(plot.title = element_text(hjust=1))
#arb2[[5]] <- arb2[[5]] + ggtitle("Germain's silver langur") + theme(plot.title = element_text(hjust=1))

ggsave("Results/ALL_PLOTS/arb_spp_1.png",arb2,dpi = 300, height = 30, width = 30, units = "cm")




#### ground based species

gb <- GPF_20_plot_bin_GS + GPF_20_plot_bin_GS_var + 
      RMJ_20_plot_bin_GS + RMJ_20_plot_bin_GS_var +
      LTM_20_plot_bin_GS + LTM_20_plot_bin_GS_var +
      PTM_20_plot_bin_GS + PTM_20_plot_bin_GS_var +
        plot_layout(ncol=2)

# remove y axis labels from all plots
gb[[1]] <- gb[[1]] + theme(axis.title.y = element_blank())
gb[[2]] <- gb[[2]] + theme(axis.title.y = element_blank())
gb[[3]] <- gb[[3]] + theme(axis.title.y = element_blank())
gb[[4]] <- gb[[4]] + theme(axis.title.y = element_blank())
gb[[5]] <- gb[[5]] + theme(axis.title.y = element_blank())
gb[[6]] <- gb[[6]] + theme(axis.title.y = element_blank())
gb[[7]] <- gb[[7]] + theme(axis.title.y = element_blank())
gb[[8]] <- gb[[8]] + theme(axis.title.y = element_blank())

# remove y axis numbers and all ticks from all plots
gb[[1]] <- gb[[1]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) 
gb[[2]] <- gb[[2]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
gb[[3]] <- gb[[3]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
gb[[4]] <- gb[[4]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
gb[[5]] <- gb[[5]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
gb[[6]] <- gb[[6]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
gb[[7]] <- gb[[7]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
gb[[8]] <- gb[[8]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())

# remove x axis labels from all plots
gb[[1]] <- gb[[1]] + theme(axis.title.x = element_blank())
gb[[2]] <- gb[[2]] + theme(axis.title.x = element_blank())
gb[[3]] <- gb[[3]] + theme(axis.title.x = element_blank())
gb[[4]] <- gb[[4]] + theme(axis.title.x = element_blank())
gb[[5]] <- gb[[5]] + theme(axis.title.x = element_blank())
gb[[6]] <- gb[[6]] + theme(axis.title.x = element_blank())
gb[[7]] <- gb[[7]] + theme(axis.title.x = element_blank())
gb[[8]] <- gb[[8]] + theme(axis.title.x = element_blank())

# remove x axis numbers and ticks from all plots
gb[[1]] <- gb[[1]] + theme(axis.text.x = element_blank())
gb[[2]] <- gb[[2]] + theme(axis.text.x = element_blank())
gb[[3]] <- gb[[3]] + theme(axis.text.x = element_blank())
gb[[4]] <- gb[[4]] + theme(axis.text.x = element_blank())
gb[[5]] <- gb[[5]] + theme(axis.text.x = element_blank())
gb[[6]] <- gb[[6]] + theme(axis.text.x = element_blank())
gb[[7]] <- gb[[7]] + theme(axis.text.x = element_blank())
gb[[8]] <- gb[[8]] + theme(axis.text.x = element_blank())

# remove legends for all plots 
gb[[1]] <- gb[[1]] + theme(legend.position = "none")
gb[[2]] <- gb[[2]] + theme(legend.position = "none")
gb[[3]] <- gb[[3]] + theme(legend.position = "none")
gb[[4]] <- gb[[4]] + theme(legend.position = "none")
gb[[5]] <- gb[[5]] + theme(legend.position = "none")
gb[[6]] <- gb[[6]] + theme(legend.position = "none")
gb[[7]] <- gb[[7]] + theme(legend.position = "none")
gb[[8]] <- gb[[8]] + theme(legend.position = "none")

# change legend position for plot 5
#gb[[5]] <- gb[[5]] + theme(legend.position = "left")

# increase the size of the legend
#gb[[5]] <- gb[[5]] + theme(legend.key.size = unit(0.8,"cm"), legend.key.width = unit(0.8,"cm"), 
                               #legend.title = element_text(size = 13))
#gb[[6]] <- gb[[6]] + theme(legend.key.size = unit(0.8,"cm"), legend.key.width = unit(0.8,"cm"),
                               #legend.title = element_text(size = 13))

ggsave("Results/ALL_PLOTS/gb_spp_1.png",gb,dpi = 300, height = 30, width = 30, units = "cm")

  ## Single plot for all species - no CV plots ####


# YCG, BSD, GSL, RMJ, GPF, LTM, PTM 

all.p <- PTM_20_plot_bin_GS + GPF_20_plot_bin_GS +
         YCG_20_plot_bin_GS + BSD_20_plot_bin_GS +
         LTM_20_plot_bin_GS + GSL_20_plot_bin_GS +
         RMJ_20_plot_bin_GS + guide_area() +
         plot_layout(ncol=2)

# remove y axis labels from all plots
all.p[[1]] <- all.p[[1]] + theme(axis.title.y = element_blank())
all.p[[2]] <- all.p[[2]] + theme(axis.title.y = element_blank())
all.p[[3]] <- all.p[[3]] + theme(axis.title.y = element_blank())
all.p[[4]] <- all.p[[4]] + theme(axis.title.y = element_blank())
all.p[[5]] <- all.p[[5]] + theme(axis.title.y = element_blank())
all.p[[6]] <- all.p[[6]] + theme(axis.title.y = element_blank())
all.p[[7]] <- all.p[[7]] + theme(axis.title.y = element_blank())

# remove y axis numbers and all ticks from all plots
all.p[[1]] <- all.p[[1]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) 
all.p[[2]] <- all.p[[2]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[3]] <- all.p[[3]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[4]] <- all.p[[4]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[5]] <- all.p[[5]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[6]] <- all.p[[6]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[7]] <- all.p[[7]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())

# remove x axis labels from all plots
all.p[[1]] <- all.p[[1]] + theme(axis.title.x = element_blank())
all.p[[2]] <- all.p[[2]] + theme(axis.title.x = element_blank())
all.p[[3]] <- all.p[[3]] + theme(axis.title.x = element_blank())
all.p[[4]] <- all.p[[4]] + theme(axis.title.x = element_blank())
all.p[[5]] <- all.p[[5]] + theme(axis.title.x = element_blank())
all.p[[6]] <- all.p[[6]] + theme(axis.title.x = element_blank())
all.p[[7]] <- all.p[[7]] + theme(axis.title.x = element_blank())

# remove x axis numbers and ticks from all plots
all.p[[1]] <- all.p[[1]] + theme(axis.text.x = element_blank())
all.p[[2]] <- all.p[[2]] + theme(axis.text.x = element_blank())
all.p[[3]] <- all.p[[3]] + theme(axis.text.x = element_blank())
all.p[[4]] <- all.p[[4]] + theme(axis.text.x = element_blank())
all.p[[5]] <- all.p[[5]] + theme(axis.text.x = element_blank())
all.p[[6]] <- all.p[[6]] + theme(axis.text.x = element_blank())
all.p[[7]] <- all.p[[7]] + theme(axis.text.x = element_blank())

# remove legend title for plots 3,4,5,6
#arb2[[3]] <- arb2[[3]] + theme(legend.title = element_blank())
#arb2[[4]] <- arb2[[4]] + theme(legend.title = element_blank())
#arb2[[5]] <- arb2[[5]] + theme(legend.title = element_blank())
#arb2[[6]] <- arb2[[6]] + theme(legend.title = element_blank())

# remove legends for plots 2,3,4,5,6,7
all.p[[2]] <- all.p[[2]] + theme(legend.position = "none")
all.p[[3]] <- all.p[[3]] + theme(legend.position = "none")
all.p[[4]] <- all.p[[4]] + theme(legend.position = "none")
all.p[[5]] <- all.p[[5]] + theme(legend.position = "none")
all.p[[6]] <- all.p[[6]] + theme(legend.position = "none")
all.p[[7]] <- all.p[[7]] + theme(legend.position = "none")

# put legend in the empty space
all.p <- all.p + plot_layout(guides = 'collect')

# increase the size of the legend
all.p[[1]] <- all.p[[1]] + theme(legend.key.size = unit(1.5,"cm"), legend.key.width = unit(2.0,"cm"), 
                                 legend.title = element_text(size = 19),
                                 legend.text = element_text(size = 15))

# add plot titles
all.p[[1]] <- all.p[[1]] + ggtitle("Pig-tailed macaque")
all.p[[2]] <- all.p[[2]] + ggtitle("Green peafowl")
all.p[[3]] <- all.p[[3]] + ggtitle("Yellow-cheeked crested gibbon")
all.p[[4]] <- all.p[[4]] + ggtitle("Black-shanked douc")
all.p[[5]] <- all.p[[5]] + ggtitle("Long-tailed macaque")
all.p[[6]] <- all.p[[6]] + ggtitle("Germain's silver langur")
all.p[[7]] <- all.p[[7]] + ggtitle("Northern red muntjac")


ggsave("Results/ALL_PLOTS/all_spp_N.png",all.p,dpi = 300, height = 30, width = 30, units = "cm")




### plot with 3 columns

# YCG, BSD, GSL, RMJ, GPF, LTM, PTM 

all.p <- PTM_20_plot_bin_GS + GPF_20_plot_bin_GS +
         YCG_20_plot_bin_GS + BSD_20_plot_bin_GS +
         LTM_20_plot_bin_GS + GSL_20_plot_bin_GS +
         RMJ_20_plot_bin_GS + guide_area() +
         plot_layout(ncol=3)

# remove y axis labels from all plots
all.p[[1]] <- all.p[[1]] + theme(axis.title.y = element_blank())
all.p[[2]] <- all.p[[2]] + theme(axis.title.y = element_blank())
all.p[[3]] <- all.p[[3]] + theme(axis.title.y = element_blank())
all.p[[4]] <- all.p[[4]] + theme(axis.title.y = element_blank())
all.p[[5]] <- all.p[[5]] + theme(axis.title.y = element_blank())
all.p[[6]] <- all.p[[6]] + theme(axis.title.y = element_blank())
all.p[[7]] <- all.p[[7]] + theme(axis.title.y = element_blank())

# remove y axis numbers and all ticks from all plots
all.p[[1]] <- all.p[[1]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) 
all.p[[2]] <- all.p[[2]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[3]] <- all.p[[3]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[4]] <- all.p[[4]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[5]] <- all.p[[5]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[6]] <- all.p[[6]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.p[[7]] <- all.p[[7]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())

# remove x axis labels from all plots
all.p[[1]] <- all.p[[1]] + theme(axis.title.x = element_blank())
all.p[[2]] <- all.p[[2]] + theme(axis.title.x = element_blank())
all.p[[3]] <- all.p[[3]] + theme(axis.title.x = element_blank())
all.p[[4]] <- all.p[[4]] + theme(axis.title.x = element_blank())
all.p[[5]] <- all.p[[5]] + theme(axis.title.x = element_blank())
all.p[[6]] <- all.p[[6]] + theme(axis.title.x = element_blank())
all.p[[7]] <- all.p[[7]] + theme(axis.title.x = element_blank())

# remove x axis numbers and ticks from all plots
all.p[[1]] <- all.p[[1]] + theme(axis.text.x = element_blank())
all.p[[2]] <- all.p[[2]] + theme(axis.text.x = element_blank())
all.p[[3]] <- all.p[[3]] + theme(axis.text.x = element_blank())
all.p[[4]] <- all.p[[4]] + theme(axis.text.x = element_blank())
all.p[[5]] <- all.p[[5]] + theme(axis.text.x = element_blank())
all.p[[6]] <- all.p[[6]] + theme(axis.text.x = element_blank())
all.p[[7]] <- all.p[[7]] + theme(axis.text.x = element_blank())

# remove legends for plots 2,3,4,5,6,7
all.p[[2]] <- all.p[[2]] + theme(legend.position = "none")
all.p[[3]] <- all.p[[3]] + theme(legend.position = "none")
all.p[[4]] <- all.p[[4]] + theme(legend.position = "none")
all.p[[5]] <- all.p[[5]] + theme(legend.position = "none")
all.p[[6]] <- all.p[[6]] + theme(legend.position = "none")
all.p[[7]] <- all.p[[7]] + theme(legend.position = "none")

# put legend in the empty space
all.p <- all.p + plot_layout(guides = 'collect')

# increase the size of the legend
all.p[[1]] <- all.p[[1]] + theme(legend.key.size = unit(1.5,"cm"), legend.key.width = unit(2.0,"cm"), 
                                 legend.title = element_text(size = 19),
                                 legend.text = element_text(size = 15))

# add plot titles
all.p[[1]] <- all.p[[1]] + ggtitle("Northern pig-tailed macaque")
all.p[[2]] <- all.p[[2]] + ggtitle("Green Peafowl")
all.p[[3]] <- all.p[[3]] + ggtitle("Southern yellow-cheeked crested gibbon")
all.p[[4]] <- all.p[[4]] + ggtitle("Black-shanked douc")
all.p[[5]] <- all.p[[5]] + ggtitle("Long-tailed macaque")
all.p[[6]] <- all.p[[6]] + ggtitle("Germain's silver langur")
all.p[[7]] <- all.p[[7]] + ggtitle("Northern red muntjac")

ggsave("Results/ALL_PLOTS/all_spp_N_wide.png",all.p,dpi = 300, height = 30, width = 30, units = "cm")

  ## Separate CV plots for all species ####


### plot with 3 columns

# YCG, BSD, GSL, RMJ, GPF, LTM, PTM 

all.cv <- PTM_20_plot_bin_GS_var + GPF_20_plot_bin_GS_var +
         YCG_20_plot_bin_GS_var + BSD_20_plot_bin_GS_var +
         LTM_20_plot_bin_GS_var + GSL_20_plot_bin_GS_var +
         RMJ_20_plot_bin_GS_var + guide_area() +
         plot_layout(ncol=3)

# remove y axis labels from all plots
all.cv[[1]] <- all.cv[[1]] + theme(axis.title.y = element_blank())
all.cv[[2]] <- all.cv[[2]] + theme(axis.title.y = element_blank())
all.cv[[3]] <- all.cv[[3]] + theme(axis.title.y = element_blank())
all.cv[[4]] <- all.cv[[4]] + theme(axis.title.y = element_blank())
all.cv[[5]] <- all.cv[[5]] + theme(axis.title.y = element_blank())
all.cv[[6]] <- all.cv[[6]] + theme(axis.title.y = element_blank())
all.cv[[7]] <- all.cv[[7]] + theme(axis.title.y = element_blank())

# remove y axis numbers and all ticks from all plots
all.cv[[1]] <- all.cv[[1]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank()) 
all.cv[[2]] <- all.cv[[2]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.cv[[3]] <- all.cv[[3]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.cv[[4]] <- all.cv[[4]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.cv[[5]] <- all.cv[[5]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.cv[[6]] <- all.cv[[6]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())
all.cv[[7]] <- all.cv[[7]] + theme(axis.text.y = element_blank(), axis.ticks = element_blank())

# remove x axis labels from all plots
all.cv[[1]] <- all.cv[[1]] + theme(axis.title.x = element_blank())
all.cv[[2]] <- all.cv[[2]] + theme(axis.title.x = element_blank())
all.cv[[3]] <- all.cv[[3]] + theme(axis.title.x = element_blank())
all.cv[[4]] <- all.cv[[4]] + theme(axis.title.x = element_blank())
all.cv[[5]] <- all.cv[[5]] + theme(axis.title.x = element_blank())
all.cv[[6]] <- all.cv[[6]] + theme(axis.title.x = element_blank())
all.cv[[7]] <- all.cv[[7]] + theme(axis.title.x = element_blank())

# remove x axis numbers and ticks from all plots
all.cv[[1]] <- all.cv[[1]] + theme(axis.text.x = element_blank())
all.cv[[2]] <- all.cv[[2]] + theme(axis.text.x = element_blank())
all.cv[[3]] <- all.cv[[3]] + theme(axis.text.x = element_blank())
all.cv[[4]] <- all.cv[[4]] + theme(axis.text.x = element_blank())
all.cv[[5]] <- all.cv[[5]] + theme(axis.text.x = element_blank())
all.cv[[6]] <- all.cv[[6]] + theme(axis.text.x = element_blank())
all.cv[[7]] <- all.cv[[7]] + theme(axis.text.x = element_blank())

# remove legends for plots 2,3,4,5,6,7
all.cv[[2]] <- all.cv[[2]] + theme(legend.position = "none")
all.cv[[3]] <- all.cv[[3]] + theme(legend.position = "none")
all.cv[[4]] <- all.cv[[4]] + theme(legend.position = "none")
all.cv[[5]] <- all.cv[[5]] + theme(legend.position = "none")
all.cv[[6]] <- all.cv[[6]] + theme(legend.position = "none")
all.cv[[7]] <- all.cv[[7]] + theme(legend.position = "none")

# put legend in the empty space
all.cv <- all.cv + plot_layout(guides = 'collect')

# increase the size of the legend
all.cv[[1]] <- all.cv[[1]] + theme(legend.key.size = unit(1.0,"cm"), legend.key.width = unit(1.5,"cm"), 
                                 legend.title = element_text(size = 19),
                                 legend.text = element_text(size = 12))

# add plot titles
all.cv[[1]] <- all.cv[[1]] + ggtitle("Pig-tailed macaque")
all.cv[[2]] <- all.cv[[2]] + ggtitle("Green peafowl")
all.cv[[3]] <- all.cv[[3]] + ggtitle("Yellow-cheeked crested gibbon")
all.cv[[4]] <- all.cv[[4]] + ggtitle("Black-shanked douc")
all.cv[[5]] <- all.cv[[5]] + ggtitle("Long-tailed macaque")
all.cv[[6]] <- all.cv[[6]] + ggtitle("Germain's silver langur")
all.cv[[7]] <- all.cv[[7]] + ggtitle("Northern red muntjac")

ggsave("Results/ALL_PLOTS/all_spp_CV_wide.png",all.cv,dpi = 300, height = 30, width = 30, units = "cm")

#### Post-submission calculations --------------------------------------------------

# Reviewer 3 wants the precision of the spatial models to be discussed in the manuscript. Currently the CV plots are in the Supporting Information. I will add a sentence or two in the Results about the range of CV values.

# load spdf's with the CV values
gpf_cv <- read.csv("Results/GPF/Plots/variance/core_only/spdf/gpf.var.spdf.df_20.core.csv")
bsd_cv <- read.csv("Results/BSD/Plots/variance/core_only/spdf/bsd.var.spdf.df_20.core.csv")
rmj_cv <- read.csv("Results/RMJ/Plots/variance/core_only/spdf/rmj.var.spdf.df_20.core.csv")
ptm_cv <- read.csv("Results/PTM/Plots/variance/core_only/spdf/ptm.var.spdf.df_20.core.csv")
ltm_cv <- read.csv("Results/LTM/Plots/variance/core_only/spdf/ltm.var.spdf.df_20.core.csv")
gsl_cv <- read.csv("Results/GSL/Plots/variance/core_only/spdf/gsl.var.spdf.df_20.core.csv")
ycg_cv <- read.csv("Results/YCG/Plots/variance/core_only/spdf/ycg.var.spdf.df_20.core.csv")

gpf_sum <- summary(gpf_cv$group2)
bsd_sum <- summary(bsd_cv$group2)
rmj_sum <- summary(rmj_cv$group2)
ptm_sum <- summary(ptm_cv$group2)
ltm_sum <- summary(ltm_cv$group2)
gsl_sum <- summary(gsl_cv$group2)
ycg_sum <- summary(ycg_cv$group2)

lst <- list(gpf_sum,bsd_sum,rmj_sum,ptm_sum,ltm_sum,gsl_sum,ycg_sum)

medians <- sapply(lst, function(x) x[3])
