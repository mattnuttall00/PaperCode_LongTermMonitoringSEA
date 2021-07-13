#### This script is for the fitting of population trends for the key species surveyed on line transects in Keo Seima Wildlife Sanctuary between 2010 and 2020.  These trends use the population point estimates that were produced using the script "Conventional_DistSamp_abundance_estimates.R". These results in combination with the density surface models (See "Density_surface_models.R") using the same data will form one chapter of my PhD and have been submitted to Conservation Science and Practice.  The bootstrapping code ("Rfunc.R" source file) was written by Prof Rachel Fewster, who will be a co-author on this chapter/paper


#### Source script (libraries, data, functions) ####


## sources R code and libraries, and loads the allData data-frame
source("Rfunc.R")   

## Create effort data "eff_dat.2":
eff_dat.2 <- createEffortData.func(filename="Input/sample.table.csv")

#### Run bootstrap ####


## Run bootstrap: this generates n replicates using method 2.
## Method 2 = stop when overall effort in each habitat across years reaches or exceeds that in the real data.
set.seed(123)
boot.res <- bootstrap.func(Nrep=2000, method=2) # 


# explore the data
head(boot.res[[1]]$repData)
head(boot.res[[1]]$sampleInfo)

# which "new" transect label corresponds to original transect 12 in this particular BS replicate
boot.res[[5]]$sampleInfo[12,]

unique(boot.res[[5]]$repData$obs.transectOrig[boot.res[[5]]$repData$Transect==12])

  ## Edit strata from repData[[i]]$stratum ####


### THIS IS NO LONGER NECESSARY AS THE EFFORT STRATA HAVE BEEN REMOVED FROM ALLDATA PRIOR TO RUNNING BOOTSTRAPPING


## The "stratum" column in each dataframe of boot.res still has the priority strata, which we need to remove and then make numeric for incusion as a linear predictor in some detection function models

# change the column from factor to character
boot.res <- lapply(boot.res, function(x) {x$repData$stratum <- as.character(x$repData$stratum);x})

# remove strata
boot.res <- lapply(boot.res, function(x) {
  x$repData$stratum[x$repData$stratum=="2011_H"] <- "2011"
  x$repData$stratum[x$repData$stratum=="2011_L"] <- "2011"
  x$repData$stratum[x$repData$stratum=="2014_H"] <- "2014"
  x$repData$stratum[x$repData$stratum=="2014_L"] <- "2014"
  x$repData$stratum[x$repData$stratum=="2016_H"] <- "2016"
  x$repData$stratum[x$repData$stratum=="2016_L"] <- "2016"
;x})

boot.res <- lapply(boot.res, function(x) {x$repData$stratum <- as.numeric(x$repData$stratum);x})

str(boot.res[[20]]$repData)



  ## Plot bootstrap results ####

## Plot the number of sightings among the bootstrap replicates and real data and print a mean summary:
bootstrap.plot(boot.res, habitat="All", plotwhat="Sightings")

## Plot (say) the number of sightings in Dense habitat:
bootstrap.plot(boot.res, habitat="Dense", plotwhat="Sightings")

## Plot the effort per year instead of the number of sightings:
bootstrap.plot(boot.res, habitat="All", plotwhat="Effort")

## Loop to plot number of sightings overall, and in each of the habitats:
for(hab in c("All", "Dense", "Open", "Nonf")){
  cat("\n", hab, "habitat:\n")
  print(bootstrap.plot(boot.res, habitat=hab, plotwhat="Sightings"))
  readline("Enter for next plot...")
}

## Loop to plot amount of effort overall, and in each of the habitats:
for(hab in c("All", "Dense", "Open", "Nonf")){
  cat("\n", hab, "habitat:\n")
  print(bootstrap.plot(boot.res, habitat=hab, plotwhat="Effort"))
  readline("Enter for next plot...")
}

#### Trend estimation ####
  ## Trends in abundance - GAMs ####

# load in abundance estimates for all species
ycg.dat <- read.csv("Output/Results/YCG_results_final.csv",        header = TRUE)
bsd.dat <- read.csv("Output/Results/BSD_results_final.csv",        header = TRUE)
btg.dat <- read.csv("Output/Results/BTG_results_final.csv",        header = TRUE)
gau.dat <- read.csv("Output/Results/GAU_results_final.csv",        header = TRUE)
gpf.dat <- read.csv("Output/Results/GPF_results_final.csv",        header = TRUE)
gsl.dat <- read.csv("Output/Results/GSL_results_final_binned.csv", header = TRUE)
ltm.dat <- read.csv("Output/Results/LTM_results_final_binned.csv", header = TRUE)
pig.dat <- read.csv("Output/Results/PIG_results_final_binned.csv", header = TRUE)
ptm.dat <- read.csv("Output/Results/PTM_results_final.csv",        header = TRUE)
rmj.dat <- read.csv("Output/Results/RMJ_results_final_binned.csv", header = TRUE)
stm.dat <- read.csv("Output/Results/STM_results_final.csv",        header = TRUE)

# load in quantiles for all spcies
ycg.quants <- read.csv("Output/Results/Plots/BS_quants/ycg.quants.csv")
bsd.quants <- read.csv("Output/Results/Plots/BS_quants/bsd.quants.csv")
#btg.quants <- read.csv("Output/Results/Plots/BS_quants/btg.quants.csv")
gau.quants <- read.csv("Output/Results/Plots/BS_quants/gau.quants.csv")
gpf.quants <- read.csv("Output/Results/Plots/BS_quants/gpf.quants.csv")
gsl.quants <- read.csv("Output/Results/Plots/BS_quants/gsl.quants.csv")
ltm.quants <- read.csv("Output/Results/Plots/BS_quants/ltm.quants.csv")
pig.quants <- read.csv("Output/Results/Plots/BS_quants/pig.quants.csv")
ptm.quants <- read.csv("Output/Results/Plots/BS_quants/ptm.quants.csv")
rmj.quants <- read.csv("Output/Results/Plots/BS_quants/rmj.quants.csv")
stm.quants <- read.csv("Output/Results/Plots/BS_quants/stm.quants.csv")


  ## Functions (plots & trends) ####

# 95% CIs only, colour
plot95fun <- function(dat,quants,species,label,ylab,ymax){
  
  ggplot()+
    ylim(0,ymax)+
    geom_point(data=dat[dat$Label==label,], aes(x=Year, y=Estimate), size=3.5)+
    geom_errorbar(data=dat[dat$Label==label,],aes(x=Year, ymin=lcl, ymax=ucl),width=0.2)+
    geom_line(data=quants, aes(x=year, y=Q50), colour="red", size=1)+
    geom_line(data=quants, aes(x=year, y=Q2.5), colour="red", linetype="dashed", size=0.8)+
    geom_line(data=quants, aes(x=year, y=Q97.5), colour="red", linetype="dashed", size=0.8)+
    scale_x_continuous(breaks = c(2010,2011,2013,2014,2016,2018,2020))+
    theme(panel.background = element_blank())+
    theme(axis.line = element_line(colour = "black"))+
    ggtitle(species)+
    ylab(ylab)+
    xlab("Year")
}

# 95 & 85% CIs, colour
plot85fun <- function(dat,quants,species,label,ylab,ymax){
  
  ggplot()+
    ylim(0,ymax)+
    geom_point(data=dat[dat$Label==label,], aes(x=Year, y=Estimate), size=3.5)+
    geom_errorbar(data=dat[dat$Label==label,],aes(x=Year, ymin=lcl, ymax=ucl),width=0.2)+
    geom_line(data=quants, aes(x=year, y=Q50), colour="red", size=1)+
    geom_line(data=quants, aes(x=year, y=Q2.5), colour="red", linetype="dashed", size=0.8)+
    geom_line(data=quants, aes(x=year, y=Q97.5), colour="red", linetype="dashed", size=0.8)+
    geom_line(data=quants, aes(x=year, y=Q7.5), colour="blue", linetype="dashed", size=0.8)+
    geom_line(data=quants, aes(x=year, y=Q92.5), colour="blue", linetype="dashed", size=0.8)+
    scale_x_continuous(breaks = c(2010,2011,2013,2014,2016,2018,2020))+
    theme(panel.background = element_blank())+
    theme(axis.line = element_line(colour = "black"))+
    ggtitle(species)+
    ylab(ylab)+
    xlab("Year")
}

# 95 & 85% CIs, black and white
plot85Grfun <- function(dat,quants,species,label,ylab,ymax){
  
  ggplot()+
    ylim(0,ymax)+
    geom_point(data=dat[dat$Label==label,], aes(x=Year, y=Estimate), size=3.5)+
    geom_errorbar(data=dat[dat$Label==label,],aes(x=Year, ymin=lcl, ymax=ucl),width=0.2)+
    geom_line(data=quants, aes(x=year, y=Q50), size=1)+
    geom_line(data=quants, aes(x=year, y=Q2.5), linetype="dotted", size=0.8)+
    geom_line(data=quants, aes(x=year, y=Q97.5), linetype="dotted", size=0.8)+
    geom_line(data=quants, aes(x=year, y=Q7.5), linetype="dashed", size=0.8)+
    geom_line(data=quants, aes(x=year, y=Q92.5), linetype="dashed", size=0.8)+
    scale_x_continuous(breaks = c(2010,2011,2013,2014,2016,2018,2020))+
    theme(panel.background = element_blank())+
    theme(axis.line = element_line(colour = "black"))+
    ggtitle(species)+
    ylab(ylab)+
    xlab("Year")
}

# 95% CIs, black and white
plot95Grfun <- function(dat,quants,species,label,ylab,ymax){
  
  ggplot()+
    ylim(0,ymax)+
    geom_point(data=dat[dat$Label==label,], aes(x=Year, y=Estimate), size=3.5)+
    geom_errorbar(data=dat[dat$Label==label,],aes(x=Year, ymin=lcl, ymax=ucl),width=0.2)+
    geom_line(data=quants, aes(x=year, y=Q50), size=1)+
    geom_line(data=quants, aes(x=year, y=Q2.5), linetype="dotted", size=0.8)+
    geom_line(data=quants, aes(x=year, y=Q97.5), linetype="dotted", size=0.8)+
    scale_x_continuous(breaks = c(2010,2011,2013,2014,2016,2018,2020))+
    theme(panel.background = element_blank())+
    theme(axis.line = element_line(colour = "black"))+
    ggtitle(species)+
    ylab(ylab)+
    xlab("Year")
  
}

# 85 and 95% CIs, black and white, faded points, ribbon CIs
plot95Grfun2 <- function(dat,quants,species,label,ylab,ymax){
  
  ggplot() +
    ylim(0,ymax)+
    geom_point(data=dat[dat$Label==label,], aes(x=Year, y=Estimate), size=2, colour = "grey70") +
    geom_errorbar(data=dat[dat$Label==label,],aes(x=Year, ymin=lcl, ymax=ucl),width=0.1, colour = "grey80") +
    geom_ribbon(data=quants, aes(x=year, ymin=Q2.5, ymax=Q97.5), fill="grey60", alpha=0.3) + # 95% CIs
    geom_ribbon(data=quants, aes(x=year, ymin=Q7.5, ymax=Q92.5), fill="grey70", alpha=0.3) + # 85% CIs
    geom_line(data=quants, aes(x=year, y=Q50), size=1) +
    scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2020)) +
    ggtitle(species) +
    ylab(ylab) +
    xlab("Year") +
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "grey20"),
          axis.title = element_text(colour = "grey20"),
          axis.text = element_text(colour = "grey20"), 
          axis.ticks = element_line(colour = "grey20"),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          plot.title = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))
  
  
}

# 85 and 95% CIS, black and white, no points, ribbon CIs
plot95Grfun3 <- function(quants,species,label,ylab,ymax){
  
  ggplot() +
    ylim(0,ymax)+
    geom_ribbon(data=quants, aes(x=year, ymin=Q2.5, ymax=Q97.5), fill="grey60", alpha=0.3) + # 95% CIs
    geom_ribbon(data=quants, aes(x=year, ymin=Q7.5, ymax=Q92.5), fill="grey70", alpha=0.3) + # 85% CIs
    geom_line(data=quants, aes(x=year, y=Q50), size=1) +
    scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2020)) +
    ggtitle(species) +
    ylab(ylab) +
    xlab("Year") +
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "grey20"),
          axis.title = element_text(colour = "grey20"),
          axis.text = element_text(colour = "grey20"), 
          axis.ticks = element_line(colour = "grey20"),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          plot.title = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))
  
  
}

# 95% CIs, black and white, faded points, ribbon CIs, solid line 
plot95Grfun4 <- function(dat,quants,species,label,ylab,ymax){
  
  ggplot() +
    ylim(0,ymax)+
    geom_point(data=dat[dat$Label==label,], aes(x=Year, y=N), size=2, colour = "grey50") +
    geom_errorbar(data=dat[dat$Label==label,],aes(x=Year, ymin=n_lcl, ymax=n_ucl),width=0.1, colour = "grey50") +
    geom_ribbon(data=quants, aes(x=year, ymin=Q2.5, ymax=Q97.5), fill="grey60", alpha=0.3) + # 95% CIs
    geom_line(data=quants, aes(x=year, y=Q50), size=1) +
    scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2020)) +
    ggtitle(species) +
    ylab(ylab) +
    xlab("Year") +
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "grey20"),
          axis.title = element_text(colour = "grey20"),
          axis.text = element_text(colour = "grey20"), 
          axis.ticks = element_line(colour = "grey20"),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          plot.title = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))
  
  
}

# 95% CIs, black and white, faded points, ribbon CIs, dashed line 
plot95Grfun5 <- function(dat,quants,species,label,ylab,ymax){
  
  ggplot() +
    ylim(0,ymax)+
    geom_point(data=dat[dat$Label==label,], aes(x=Year, y=Estimate), size=2, colour = "grey70") +
    geom_errorbar(data=dat[dat$Label==label,],aes(x=Year, ymin=lcl, ymax=ucl),width=0.1, colour = "grey80") +
    geom_ribbon(data=quants, aes(x=year, ymin=Q2.5, ymax=Q97.5), fill="grey60", alpha=0.3) + # 95% CIs
    geom_line(data=quants, aes(x=year, y=Q50), size=1, linetype="dashed") +
    scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2020)) +
    ggtitle(species) +
    ylab(ylab) +
    xlab("Year") +
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "grey20"),
          axis.title = element_text(colour = "grey20"),
          axis.text = element_text(colour = "grey20"), 
          axis.ticks = element_line(colour = "grey20"),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          plot.title = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))
  
  
}



### trend functions

# calculates whether each replicate is positive or negative trend (bootstrapping method)
trendFunc <- function(x){
  ifelse(x[1] < x[100], "positive", "negative")
}

# function to apply trendFunc to a single dataframe within a list (see "Trends" section)
calcfun <- function(c){
  df <- data.frame(apply(c,2,trendFunc))
  tbl <- table(df)
  return(tbl)
}


# 85% CI overlapping method
trendFunc2 <- function(x){
  ifelse(x$Q7.5[1] > x$Q92.5[100], "negative",
         ifelse(x$Q92.5[1] < x$Q7.5[100], "positive", "stable"))
}


    ## YCG ####

ycg.dat <- select(ycg.dat, -X)
ycg.dat$year <- as.numeric(ycg.dat$year)
str(ycg.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
ycg.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  gib.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), 
                    data=ycg.dat[ycg.dat$Label=="Grp",])
  ycg.aic.res[paste0("df", dfval)] <- summary(gib.gam.df)$aic }


# Choose the fit with lowest AIC:
which.min(ycg.aic.res) # df=1


## Get confidence intervals from the bootstrapped replicates

# Function for each species. These functions extract the species-specific transect data from a single bootstrapping replicate, create the required files for ds(), fits a detection function (using the model formulation from the CDS analysis for that species), fits a GAM to the annual point estimates, then predicts estimates across a range of "year" values (length=100, year as numeric) which is used to fit the smooth lines

# the detection funcion model formulation for each species is the same as the final DF model in the CDS analysis (see "CDS_allSpp_2010-2018.R"), except no covariates will be included as this causes covariance issues due to the nature of the bootstrapping

# The GAM model formulation in each function for each species is the same as the model formulation used in getting the main trend line. 

# function to fit a GAM to each replicate of boot.res
fitspecies.func.YCG <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="YCG",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model
  try(detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 60, key = "hn"))
  
  
  # extract estimates
  estimates <- detfunc$dht$individuals$N[1:7, 1:2]
  estimates <- estimates %>% dplyr::rename(Year = Label) 
  estimates$Year <- as.numeric(estimates$Year)
  
  # fit a GAM & predict
  gamfit <- gam(Estimate ~ s(Year, df=1), family=gaussian(link="identity"), data = estimates)
  newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100))
  pred <- predict.Gam(gamfit, newdata = newdata, type = "response")
  gampred <- data.frame(pred)
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(ycg.bs.gams <- lapply(1:length(boot.res), fitspecies.func.YCG)) #70 mins

# put output list into a dataframe
ycg.bs.gams.df <- data.frame(matrix(unlist(ycg.bs.gams), nrow=100, byrow = FALSE))

# save
write.csv(ycg.bs.gams.df, file="Output/Results/Trends/Bootstraps/ycg.bs.gams.df.csv")

# read
#ycg.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/ycg.bs.gams.df.csv")
#ycg.bs.gams.df <- ycg.bs.gams.df %>% select(-X)

# test what percentage of the last BS gam estimates are above the first. This is to test for a significant upward trend 
trend.df.ycg <- data.frame(apply(ycg.bs.gams.df,2,trendFunc))
table(trend.df.ycg)
# 89% of BS estimates suggest a positive trend, so the trend is not signficant. 

# extract 2.5, 7.5, 50, 92.5 and 97.5 quantiles (85% and 95% CIs)
ycg.bs_ints <- data.frame(apply(ycg.bs.gams.df, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
ycg.bs_ints <- ycg.bs_ints %>% rownames_to_column("quant")


# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
ycg.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(ycg.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(ycg.bs_ints[2,2:101]),
                         Q50 = as.numeric(ycg.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(ycg.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(ycg.bs_ints[5,2:101]))
# save quantiles
write.csv(ycg.quants, "Output/Results/Plots/BS_quants/ycg.quants.csv")

# load quantiles
ycg.quants <- read.csv("Output/Results/Plots/BS_quants/ycg.quants.csv")
str(ycg.quants)
ycg.quants <- ycg.quants[ ,-1]



### plots

# 95% CIs, colour
ycg_plot_95 <- plot95fun(ycg.dat,ycg.quants,"Yellow-cheeked crested Gibbon","Grp","Group abundance",1400)

# 95 & 85% CIs, black and white
ycg_plot_85_gr <- plot85Grfun(ycg.dat,ycg.quants,"Yellow-cheeked crested Gibbon","Grp",
                              "Group abundance",1400)

# 95% black and white
ycg_plot_95_gr <- plot95Grfun(ycg.dat,ycg.quants,"Yellow-cheeked crested gibbon","Ind",
                              "Individual abundance",3100)


    ## GSL ####

gsl.dat <- select(gsl.dat, -X)
gsl.dat$year <- as.numeric(gsl.dat$year)
str(gsl.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
gsl.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  gsl.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), 
                    data=gsl.dat[gsl.dat$Label=="Grp",])
  gsl.aic.res[paste0("df", dfval)] <- summary(gsl.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(gsl.aic.res) # df=3


## Get confidence intervals from the bootstrapped replicates

# Function for each species. For more details see this section in the YCG section above. 

# function to fit a GAM to each replicate of boot.res. The ke ymodel used in the CDS is HR, but the function keeps failing after a varying number of model fits with an error about the variance-covariance matrix.  I wonder whether this is to do with the HR model. I have checked in the CDS analysis and HN/Unif were preferred models by AIC, and I had concluded that HR was the most ecologically accurate. BUt I had also concluded that the resutls were pretty much the same. So I will test the below function with HN to see if it still fails.
fitspecies.func.GSL <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="GSL",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model
  try(detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 60, key = "hn",cutpoints = c(0,10,19,30,50,60)))
  
  
  # extract estimates
  estimates <- detfunc$dht$individuals$N[1:7, 1:2]
  estimates <- estimates %>% dplyr::rename(Year = Label) 
  estimates$Year <- as.numeric(estimates$Year)
  
  # fit a GAM & predict
  gamfit <- gam(Estimate ~ s(Year, df=1), family=gaussian(link="identity"), data = estimates)
  newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100))
  pred <- predict.Gam(gamfit, newdata = newdata, type = "response")
  gampred <- data.frame(pred)
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(gsl.bs.gams <- lapply(1:length(boot.res), fitspecies.func.GSL)) # 2 hrs 21 mins

# put output list into a dataframe
gsl.bs.gams.df <- data.frame(matrix(unlist(gsl.bs.gams), nrow=100, byrow = FALSE))

# save (to avoid having to re-run)
write.csv(gsl.bs.gams.df, file="Output/Results/Trends/Bootstraps/gsl.bs.gams.df.csv")

# load
#gsl.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/gsl.bs.gams.df.csv")
#gsl.bs.gams.df <- gsl.bs.gams.df[ ,-1]

# test what percentage of the last BS gam estimates are above the first. This is to test for a significant upward trend 
trend.df <- data.frame(apply(gsl.bs.gams.df,2,trendFunc))
table(trend.df)
# 93% of replicates suggest negative trend

# extract 2.5, 7.5, 50, 92.5 and 97.5 quantiles (85 & 95% CIs)
gsl.bs_ints <- data.frame(apply(gsl.bs.gams.df, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
gsl.bs_ints <- gsl.bs_ints %>% rownames_to_column("quant")

# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
gsl.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(gsl.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(gsl.bs_ints[2,2:101]),
                         Q50 = as.numeric(gsl.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(gsl.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(gsl.bs_ints[5,2:101]))
# save quantiles
write.csv(gsl.quants, "Output/Results/Plots/BS_quants/gsl.quants.csv")

# load quantiles
gsl.quants <- read.csv("Output/Results/Plots/BS_quants/gsl.quants.csv")
gsl.quants <- gsl.quants[ ,-1]


### plots

# 95% colour
gsl_plot_95 <- plot95fun(gsl.dat,gsl.quants,"Silver langur","Grp","Group abundance",2500)

# 95 & 85% black and white
gsl_plot_85_gr <- plot85Grfun(gsl.dat,gsl.quants,"Silver langur","Grp","Group abundance",2500)

# 95% black and white
gsl_plot_95_gr <- plot95Grfun(gsl.dat,gsl.quants,"Silver langur","Ind","Group abundance",10500)


    ## LTM ####

ltm.dat <- select(ltm.dat, -X)
ltm.dat$year <- as.numeric(ltm.dat$year)
str(ltm.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
ltm.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  ltm.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), data=ltm.dat[ltm.dat$Label=="Grp",])
  ltm.aic.res[paste0("df", dfval)] <- summary(ltm.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(ltm.aic.res) #df=1



## Get confidence intervals from the bootstrapped replicates

# Function for each species. For more details see this section in the YCG section above. 

# function to fit a GAM to each replicate of boot.res
fitspecies.func.LTM <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="LTM",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model
  try(detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 50, key = "hn",cutpoints = c(0,10,20,27,35,50)))
  
  
  # extract estimates
  estimates <- detfunc$dht$individuals$N[1:7, 1:2]
  estimates <- estimates %>% dplyr::rename(Year = Label) 
  estimates$Year <- as.numeric(estimates$Year)
  
  # fit a GAM & predict
  gamfit <- gam(Estimate ~ s(Year, df=1), family=gaussian(link="identity"), data = estimates)
  newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100))
  pred <- predict.Gam(gamfit, newdata = newdata, type = "response")
  gampred <- data.frame(pred)
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(ltm.bs.gams <- lapply(1:length(boot.res), fitspecies.func.LTM))

# put output list into a dataframe
ltm.bs.gams.df <- data.frame(matrix(unlist(ltm.bs.gams), nrow=100, byrow = FALSE))

# save (to avoid having to re-run)
write.csv(ltm.bs.gams.df, file="Output/Results/Trends/Bootstraps/ltm.bs.gams.df.csv")

# load (if required)
#ltm.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/ltm.bs.gams.df.csv")
#ltm.bs.gams.df <- ltm.bs.gams.df[ ,-1]


# test what percentage of the last BS gam estimates are above the first. This is to test for a significant upward trend 
trend.df <- data.frame(apply(ltm.bs.gams.df,2,trendFunc))
table(trend.df)
# 75% of replicates suggest negative trend


# extract 2.5,50 and 97.5 quantiles
ltm.bs_ints <- data.frame(apply(ltm.bs.gams.df, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
ltm.bs_ints <- ltm.bs_ints %>% rownames_to_column("quant")

# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
ltm.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(ltm.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(ltm.bs_ints[2,2:101]),
                         Q50 = as.numeric(ltm.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(ltm.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(ltm.bs_ints[5,2:101]))

# save quantiles
write.csv(ltm.quants, "Output/Results/Plots/BS_quants/ltm.quants.csv")

# load quantiles
ltm.quants <- read.csv("Output/Results/Plots/BS_quants/ltm.quants.csv")
ltm.quants <- ltm.quants[ ,-1]


### plots 

# 95% colour
ltm_plot_95 <- plot95fun(ltm.dat,ltm.quants,"Long-tailed macaque","Grp","Group abundance",2000)

# 95 & 85% black and white
ltm_plot_85_gr <- plot85Grfun(ltm.dat,ltm.quants,"Long-tailed macaque","Grp","Group abundance",2000)

# 95% black and white
ltm_plot_95_gr <- plot95Grfun(ltm.dat,ltm.quants,"Long-tailed macaque","Grp","Group abundance",2000)


    ## PTM ####

ptm.dat <- select(ptm.dat, -X)
ptm.dat$year <- as.numeric(ptm.dat$year)
str(ptm.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
ptm.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  ptm.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), 
                    data=ptm.dat[ptm.dat$Label=="Grp",])
  ptm.aic.res[paste0("df", dfval)] <- summary(ptm.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(ptm.aic.res) # df=1


## Get confidence intervals from the bootstrapped replicates

# Function for each species. For more details see this section in the YCG section above. 

# HR models keep failing and causing the function to fail. Switching to HN

# function to fit a GAM to each replicate of boot.res
fitspecies.func.PTM <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="PTM",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model
  try(
    detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 50, key = "hr"))
  
  # extract estimates
  try(estimates <- detfunc$dht$individuals$N[1:7, 1:2])
  try(estimates <- estimates %>% dplyr::rename(Year = Label)) 
  try(estimates$Year <- as.numeric(estimates$Year))
  
  # fit a GAM & predict
  try(gamfit <- gam(Estimate ~ s(Year, df=1), family=gaussian(link="identity"), data = estimates))
  try(newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100)))
  try(pred <- predict.Gam(gamfit, newdata = newdata, type = "response"))
  try(gampred <- data.frame(pred))
  
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(ptm.bs.gams <- lapply(1:length(boot.res), fitspecies.func.PTM)) # 3 hr 50 mins

# Some of the DF models failed and so I think there are some empty list elements which are causing issues
sapply(ptm.bs.gams, min)
which(err=="Error in data.frame(pred) : object 'pred' not found\n")
# 182, 611, 738, 781, 1531, 1631, 1917

# make copy
ptm.bs.gams2 <- ptm.bs.gams

# remove elements
ptm.bs.gams2 <- ptm.bs.gams2[- c(182, 611, 738, 781, 1531, 1631, 1917)]


# put output list into a dataframe
ptm.bs.gams.df <- data.frame(matrix(unlist(ptm.bs.gams2), nrow=100, byrow = FALSE))

# save (to avoid having to re-run)
write.csv(ptm.bs.gams.df, file="Output/Results/Trends/Bootstraps/ptm.bs.gams.df.csv")

# load
#ptm.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/ptm.bs.gams.df.csv")
#ptm.bs.gams.df <- ptm.bs.gams.df[ ,-1]


# test what percentage of the last BS gam estimates are above the first. This is to test for a significant upward trend 
trend.df <- data.frame(apply(ptm.bs.gams.df,2,trendFunc))
table(trend.df)
# 96.5% of replicates suggest positive trend

# extract 2.5, 7.5, 50, 92.5, and 97.5 quantiles
ptm.bs_ints <- data.frame(apply(ptm.bs.gams.df, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
ptm.bs_ints <- ptm.bs_ints %>% rownames_to_column("quant")

# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
ptm.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(ptm.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(ptm.bs_ints[2,2:101]),
                         Q50 = as.numeric(ptm.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(ptm.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(ptm.bs_ints[5,2:101]))

# save quantiles
write.csv(ptm.quants, "Output/Results/Plots/BS_quants/ptm.quants.csv")

# load quantiles
ptm.quants <- read.csv("Output/Results/Plots/BS_quants/ptm.quants.csv")
ptm.quants <- ptm.quants[ ,-1]


### plots

# 95% colour
ptm_plot_95 <- plot95fun(ptm.dat,ptm.quants,"Pig-tailed macaque","Grp","Group abundance",2700)

# 95 & 85% black and white
ptm_plot_85_gr <- plot85Grfun(ptm.dat,ptm.quants,"Pig-tailed macaque","Grp","Group abundance",2700)

# 95% black and white
ptm_plot_95_gr <- plot95Grfun(ptm.dat,ptm.quants,"Pig-tailed macaque","Ind","Group abundance",8750)



    ## STM ####

stm.dat <- select(stm.dat, -X)
stm.dat$year <- as.numeric(stm.dat$year)
str(stm.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
stm.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  stm.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), 
                    data=stm.dat[stm.dat$Label=="Grp",])
  stm.aic.res[paste0("df", dfval)] <- summary(stm.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(stm.aic.res) # 1


## Get confidence intervals from the bootstrapped replicates

# I have changed the key function for the DF in the function below from HR to HN. This is because the models were struggling to produce estimates from the replicates, which was causing the bootstrapped CIs to be ridiculous. I re-ran the bootstrapping with 10 replicates, and then tried with HN key, and the results were better. After consultation with Rachel, we have decided to re-run the 2000 reps with a HN key function to see what improvements there are. 

# function to fit a GAM to each replicate of boot.res
fitspecies.func.STM <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="STM",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model
  try(detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 35, key = "hn"))
  
  
  # extract estimates
  try(estimates <- detfunc$dht$individuals$N[1:7, 1:2])
  try(estimates <- estimates %>% dplyr::rename(Year = Label)) 
  try(estimates$Year <- as.numeric(estimates$Year))
  
  # fit a GAM & predict
  try(gamfit <- gam(Estimate ~ s(Year, df=1), family=gaussian(link="identity"), data = estimates))
  try(newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100)))
  try(pred <- predict.Gam(gamfit, newdata = newdata, type = "response"))
  try(gampred <- data.frame(pred))
  
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(stm.bs.gams <- lapply(1:length(boot.res), fitspecies.func.STM)) 


# When running the bootstrapping for individuals, some of the replicates failed because the detection functions couldn't be fitted. That means that some of the elements in the output list were not dataframes of predictions, but were elements of class "try-error" (from the try() code in the main funcion above).  
# Loop through list and identify elements that are not data.frames (i.e. they are try-errors)
res <- data.frame()
for(i in 1:length(stm.bs.gams)){
  result <- class(stm.bs.gams[[i]])
  res <- rbind(res,result)
}

# add index and filter
res$num <- 1:2000
rem <- res %>% filter(X.data.frame.=="try-error") %>% select(num)
# elements to remove are: 21,180,271,305,796,1049,1582,1665,1669,1775,1841,1886


# make a copy of original list
stm.bs.gams2 <- stm.bs.gams

# remove elements with errors
stm.bs.gams2 <- stm.bs.gams2[-c(21, 180, 271, 305, 796, 1049, 1582, 1665, 1669, 1775, 1841, 1886)]

# put into df
stm.bs.gams.df <- data.frame(matrix(unlist(stm.bs.gams2), nrow=100, byrow = FALSE))

# save GAM predictions (to avoid having to re-run)
write.csv(stm.bs.gams.df, file="Output/Results/Trends/Bootstraps/stm.bs.gams.df.csv")

# put Estimate output list into a dataframe
stm.bs.est.df <- data.frame(matrix(unlist(stm.est.res), nrow=7, byrow = FALSE))
stm.bs.est.df <- stm.bs.est.df[ ,seq(2,length(stm.bs.est.df), by=2)]
stm.bs.est.df$year <- c("2010","2011","2013","2014","2016","2018","2020")

# save estimates
write.csv(stm.bs.est.df, file="Output/Results/Trends/Bootstraps/stm.bs.ests.df.csv")


# load GAMS
#stm.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/stm.bs.gams.df.csv")
#stm.bs.gams.df <- stm.bs.gams.df[ ,-1]

# load estimates
#stm.est.res <- read.csv("Output/Results/Trends/Bootstraps/stm.bs.ests.df.csv")




### The bootstrapped quantile CIs were coming out pretty wild, so below are some diagnostics

# split the GAM predictions from the point estimates from each replicate
stm.pred.res <- lapply(stm.bs.gams, function(x)x[[1]])
stm.est.res <- lapply(stm.bs.gams, function(x)x[[2]])

# extract highest annual point estimate from each replicate
maxEst <- unlist(lapply(stm.est.res, function(x) max(x$Estimate)))
maxEstdf <- data.frame(rep = 1:2000,
                       estimate = maxEst,
                       label="Maximum (any year)")

plot(maxEstdf$rep, maxEstdf$estimate)

# other diagnostic plots for maximum estimates per replicate
par(mfrow=c(2,2))
density(maxEst)
hist(maxEst)
boxplot(maxEst)
stripchart(maxEst)
plot(sort(maxEst), type="l")


# extract only 2011 estiamtes
est11 <- unlist(lapply(stm.est.res, function(x) x$Estimate[2]))
est11df <- data.frame(rep=1:2000,
                      estimate = est11,
                      label = "2011 only")

plot(est11df$rep,est11df$est11)

# merge the two dataframes (highest estimate and 2011 estimate)
combdf <- rbind(maxEstdf,est11df)

ggplot(combdf, aes(x=rep, y=estimate, group=label, colour=label))+
  geom_point()+
  theme_bw()
 

# extract all other years estimates
est10 <- unlist(lapply(stm.est.res, function(x) x$Estimate[1]))
est13 <- unlist(lapply(stm.est.res, function(x) x$Estimate[3]))
est14 <- unlist(lapply(stm.est.res, function(x) x$Estimate[4]))
est16 <- unlist(lapply(stm.est.res, function(x) x$Estimate[5]))
est18 <- unlist(lapply(stm.est.res, function(x) x$Estimate[6]))
est20 <- unlist(lapply(stm.est.res, function(x) x$Estimate[7]))

# create dataframe for them all
est10df <- data.frame(rep=1:2000, estimate=est10, label="2010")
est11df <- data.frame(rep=1:2000, estimate=est11, label="2011")
est13df <- data.frame(rep=1:2000, estimate=est13, label="2013")
est14df <- data.frame(rep=1:2000, estimate=est14, label="2014")
est16df <- data.frame(rep=1:2000, estimate=est16, label="2016")
est18df <- data.frame(rep=1:2000, estimate=est18, label="2018")
est20df <- data.frame(rep=1:2000, estimate=est20, label="2020")

alldf <- rbind(est10df,est11df,est13df,est14df,est16df,est18df,est20df)

ggplot(alldf, aes(x=rep, y=estimate, group=label, colour=label))+
  geom_point()+
  theme_bw()


### based on all of the above (primarly the diagnostic plots for the highest estimate from each replicate), we have chosen 700 as the point above which the estimates are outliers and ridiculous. Therefore we will remove all estimates from the bootstrap replicates above 700 before doing the rest of the process below

# This just changes all values above 700 to 0, so that they can't feature in the 95% and 85% quantiles
stm.bs.gams.df.sub <- data.frame(apply(stm.bs.gams.df, 2, function(x) ifelse(x >700, 0, x)))
str(stm.bs.gams.df.sub)





# test what percentage of the last BS gam estimates are above the first. This is to test for a significant upward trend 
trend.df.stm <- data.frame(apply(stm.bs.gams.df,2,trendFunc))
table(trend.df.stm)
# 100% of replicates suggest negative trend


# extract 2.5, 7.5, 50, 92.5, and 97.5 quantiles
stm.bs_ints <- data.frame(apply(stm.bs.gams.df, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
stm.bs_ints <- stm.bs_ints %>% rownames_to_column("quant")

# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
stm.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(stm.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(stm.bs_ints[2,2:101]),
                         Q50 = as.numeric(stm.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(stm.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(stm.bs_ints[5,2:101]))

# save quantiles
write.csv(stm.quants, "Output/Results/Plots/BS_quants/stm.quants.csv")

# load quantiles
stm.quants <- read.csv("Output/Results/Plots/BS_quants/stm.quants.csv")
stm.quants <- stm.quants[ ,-1]

# change negative 2.5% quantiles to 0
stm.quants <- stm.quants %>% mutate(Q2.5 = ifelse(Q2.5 <0, 0, Q2.5))


### plots

# 95% colour
stm_plot_95 <- plot95fun(stm.dat,stm.quants,"Stump-tailed macaque","Grp","Group abundance",800)

# 95 & 85% black and white
stm_plot_85_gr <- plot85Grfun(stm.dat,stm.quants,"Stump-tailed macaque","Grp","Group abundance",800)



# as there are no observations in 2018, I want to plot the point differently (as with BTG). Therefore I will make this plot manually

# add grouping level to annual estimates so that I can change the pooint type for 2013 and 2020
stm.dat$Label2 <- ifelse(stm.dat$Estimate==0,1,0)
stm.dat$Label2 <- as.factor(stm.dat$Label2)

# plot 95 black and white
stm_plot_95_gr <- ggplot()+
                  ylim(0,4000)+
                  geom_point(data=stm.dat[stm.dat$Label=="Ind",], 
                    aes(x=Year, y=Estimate, group=Label2, shape=Label2), size=3.5)+
                  scale_shape_manual(values = c(16,1))+
                  geom_errorbar(data=stm.dat[stm.dat$Label=="Ind",],aes(x=Year, ymin=lcl, ymax=ucl),width=0.2)+
                  geom_line(data=stm.quants, aes(x=year, y=Q50), size=1)+
                  geom_line(data=stm.quants, aes(x=year, y=Q2.5), linetype="dotted", size=0.8)+
                  geom_line(data=stm.quants, aes(x=year, y=Q97.5), linetype="dotted", size=0.8)+
                  scale_x_continuous(breaks = c(2010,2011,2013,2014,2016,2018,2020))+
                  theme(panel.background = element_blank())+
                  theme(axis.line = element_line(colour = "black"))+
                  theme(legend.position="none")+
                  ggtitle("Stump-tailed macaque")+
                  ylab("Group abundance")+
                  xlab("Year")


# add grouping level to annual estimates so that I can change the pooint type for 2013 and 2020
stm.dat$Label2 <- ifelse(stm.dat$N==0,1,0)
stm.dat$Label2 <- as.factor(stm.dat$Label2)

# plot 95% with ribbons, faded points, and hollow points for no obs
stm.95 <- ggplot() +
    ylim(0,3800)+
    geom_point(data=stm.dat[stm.dat$Label=="Ind",], aes(x=Year, y=N, group=Label2, shape=Label2), 
               size=2, colour = "grey50") +
    scale_shape_manual(values=c(16,1))+
    geom_errorbar(data=stm.dat[stm.dat$Label=="Ind",],aes(x=Year, ymin=n_lcl, ymax=n_ucl),width=0.1, colour = "grey50") +
    geom_ribbon(data=stm.quants, aes(x=year, ymin=Q2.5, ymax=Q97.5), fill="grey60", alpha=0.3) + # 95% CIs
    geom_line(data=stm.quants, aes(x=year, y=Q50), size=1) +
    scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2020)) +
    ggtitle("Stump-tailed macaque") +
    ylab("Abundance") +
    xlab("Year") +
    theme(legend.position="none",
          panel.background = element_blank(), 
          axis.line = element_line(colour = "grey20"),
          axis.title = element_text(colour = "grey20"),
          axis.text = element_text(colour = "grey20"), 
          axis.ticks = element_line(colour = "grey20"),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          plot.title = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))
  
  




    ## BSD ####

bsd.dat <- select(bsd.dat, -X)
bsd.dat$year <- as.numeric(bsd.dat$year)
str(bsd.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
bsd.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  bsd.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), 
                    data=bsd.dat[bsd.dat$Label=="Grp",])
  bsd.aic.res[paste0("df", dfval)] <- summary(bsd.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(bsd.aic.res) 


## Get confidence intervals from the bootstrapped replicates

# Function for each species. For more details see this section in the YCG section above. 

# for BSD, all of the annual CDS det func models were HN, and so HN will be used. All truncation distances were 50m apart from 2018 which was 60. Therefoe 50 will be used.

# function to fit a GAM to each replicate of boot.res
fitspecies.func.BSD <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="BSD",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model. For BSD, I looked at the most common key function from all the annual analyses, and Hn is the most common. Therefore that is what will be used here.I have also not used any bins, because in the CDS anlaysis there are more years that are not binned than are binned. The most common truncation distance is 50m and so that will be used.
  try(detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 50, key = "hn"))
  
  
  # extract estimates
  estimates <- detfunc$dht$individuals$N[1:7, 1:2]
  estimates <- estimates %>% dplyr::rename(Year = Label) 
  estimates$Year <- as.numeric(estimates$Year)
  
  # fit a GAM & predict
  gamfit <- gam(Estimate ~ s(Year, df=2), family=gaussian(link="identity"), data = estimates)
  newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100))
  pred <- predict.Gam(gamfit, newdata = newdata, type = "response")
  gampred <- data.frame(pred)
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(bsd.bs.gams <- lapply(1:length(boot.res), fitspecies.func.BSD))

# put output list into a dataframe
bsd.bs.gams.df <- data.frame(matrix(unlist(bsd.bs.gams), nrow=100, byrow = FALSE))

# save (to avoid having to re-run)
write.csv(bsd.bs.gams.df, file="Output/Results/Trends/Bootstraps/bsd.bs.gams.df.csv")

# load
bsd.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/bsd.bs.gams.df.csv")
bsd.bs.gams.df <- bsd.bs.gams.df[ ,-1]

# test what percentage of the last BS gam estimates are above the first. This is to test for a significant upward trend
trend.df <- data.frame(apply(bsd.bs.gams.df,2,trendFunc))
table(trend.df)
# 87% of replicates suggest positive trend


# extract 2.5, 7.5, 50, 92.5, and 97.5 quantiles
bsd.bs_ints <- data.frame(apply(bsd.bs.gams.df, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
bsd.bs_ints <- bsd.bs_ints %>% rownames_to_column("quant")

# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
bsd.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(bsd.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(bsd.bs_ints[2,2:101]),
                         Q50 = as.numeric(bsd.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(bsd.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(bsd.bs_ints[5,2:101]))

# save quantiles
write.csv(bsd.quants, "Output/Results/Plots/BS_quants/bsd.quants.csv")

# load quantiles
bsd.quants <- read.csv("Output/Results/Plots/BS_quants/bsd.quants.csv")


### plots

# 95% colour
bsd_plot_95 <- plot95fun(bsd.dat,bsd.quants,"Black-shanked douc","Grp","Group abundance",17000)

# 95 & 85% black and white
bsd_plot_85_gr <- plot85Grfun(bsd.dat,bsd.quants,"Black-shanked douc","Grp","Group abundance",17000)

# 95% black and white
bsd_plot_95_gr <- plot95Grfun(bsd.dat,bsd.quants,"Black-shanked douc","Ind","Group abundance",50500)


    ## BTG ####

btg.dat <- select(btg.dat, -X)
btg.dat$year <- as.numeric(btg.dat$year)
btg.dat <- filter(btg.dat, data=="wcs")
str(btg.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
btg.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  btg.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), 
                    data=btg.dat[btg.dat$Label=="Ind",])
  btg.aic.res[paste0("df", dfval)] <- summary(btg.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(btg.aic.res) # df= 3

# Bootstrapping does not work for BTG, and therefore I will just fit a single GAM to the point estimates to get a trend line, with no CIs

# final model
btg.gam <- gam(N ~ s(Year, df=3), family=gaussian(link="identity"), data=btg.dat[btg.dat$Label=="Ind",])
plot.Gam(btg.gam)

# predict to get trend line
btg.newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100))
btg.gam.pred <- predict.Gam(btg.gam, newdata = btg.newdata, type = "response")
btg.pred.results <- cbind(btg.newdata,btg.gam.pred) 

# change negative trend values to 0
btg.pred.results <- btg.pred.results %>% mutate(btg.gam.pred = ifelse(btg.gam.pred<0,0,btg.gam.pred))

# add grouping level to annual estimates so that I can change the pooint type for 2013 and 2020
btg.dat$Label2 <- ifelse(btg.dat$N==0,1,0)
btg.dat$Label2 <- as.factor(btg.dat$Label2)

# plot
btg_plot <- ggplot()+
  ylim(0,2000)+
  geom_point(data=btg.dat[btg.dat$Label=="Ind",], 
             aes(x=Year, y=Estimate, group=Label2, shape=Label2), size=3.5)+
  scale_shape_manual(values = c(16,1))+
  geom_errorbar(data=btg.dat[btg.dat$Label=="Ind",],aes(x=Year, ymin=lcl, ymax=ucl),width=0.2)+
  geom_line(data=btg.pred.results, aes(x=Year, y=btg.gam.pred), colour="black", size=1)+
  scale_x_continuous(breaks = c(2010,2011,2013,2014,2016,2018,2020))+
  theme(legend.position="none")+
  ggtitle("Banteng")+
  ylab("Individual abundance")+
  xlab("Year")
  


# plot with no points
btg_plot <- ggplot()+
  ylim(0,350)+
  geom_line(data=btg.pred.results, aes(x=Year, y=btg.gam.pred), colour="black", size=1)+
  scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2020))+
  theme(legend.position="none")+
  ggtitle("Banteng")+
  ylab("Individual abundance")+
  xlab("Year")+
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "grey20"),
        axis.title = element_text(colour = "grey20"),
        axis.text = element_text(colour = "grey20"), 
        axis.ticks = element_line(colour = "grey20"),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))

# plot with faded points, 95 & 85% CIs
btg_plot2 <- ggplot()+
  ylim(0,2000)+
  geom_point(data=btg.dat[btg.dat$Label=="Ind",], aes(x=Year, y=Estimate, group=Label2, shape=Label2), 
             size=2, colour = "grey70") +
  scale_shape_manual(values = c(16,1))+
  geom_errorbar(data=btg.dat[btg.dat$Label=="Ind",],aes(x=Year, ymin=lcl, ymax=ucl),width=0.1, colour = "grey80") +
  geom_line(data=btg.pred.results, aes(x=Year, y=btg.gam.pred), colour="black", size=1)+
  scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2020))+
  theme(legend.position="none")+
  ggtitle("Banteng")+
  ylab("Abundance")+
  xlab("Year")+
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "grey20"),
        axis.title = element_text(colour = "grey20"),
        axis.text = element_text(colour = "grey20"), 
        axis.ticks = element_line(colour = "grey20"),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))


# plot with faded points, 95% CIs
btg.95 <- ggplot()+
  ylim(0,2000)+
  geom_point(data=btg.dat[btg.dat$Label=="Ind",], aes(x=Year, y=N, group=Label2, shape=Label2), 
             size=2, colour = "grey50") +
  scale_shape_manual(values=c(16,1))+
  geom_errorbar(data=btg.dat[btg.dat$Label=="Ind",],aes(x=Year, ymin=n_lcl, ymax=n_ucl),width=0.1, colour = "grey50") +
  geom_line(data=btg.pred.results, aes(x=Year, y=btg.gam.pred), colour="black", size=1)+
  scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2020))+
  theme(legend.position="none")+
  ggtitle("Banteng")+
  ylab("Individual abundance")+
  xlab("Year")+
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "grey20"),
        axis.title = element_text(colour = "grey20"),
        axis.text = element_text(colour = "grey20"), 
        axis.ticks = element_line(colour = "grey20"),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.title = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))


    ## GAU ####


gau.dat <- select(gau.dat, -X)
gau.dat$year <- as.numeric(gau.dat$year)
str(gau.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
gau.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  gau.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), data=gau.dat[gau.dat$Label=="Ind",])
  gau.aic.res[paste0("df", dfval)] <- summary(gau.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(gau.aic.res) # df = 1


## Get confidence intervals from the bootstrapped replicates

# Function for each species. For more details see this section in the YCG section above. 

# function to fit a GAM to each replicate of boot.res
fitspecies.func.GAU <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="GAU",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model
  try(detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 40, key = "hn"))
  
  
  # extract estimates
  estimates <- detfunc$dht$individuals$N[1:7, 1:2]
  estimates <- estimates %>% dplyr::rename(Year = Label) 
  estimates$Year <- as.numeric(estimates$Year)
  
  # fit a GAM & predict
  gamfit <- gam(Estimate ~ s(Year, df=1), family=gaussian(link="identity"), data = estimates)
  newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100))
  pred <- predict.Gam(gamfit, newdata = newdata, type = "response")
  gampred <- data.frame(pred)
  
  return(list(gampred, estimates))
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(gau.bs.gams <- lapply(1:length(boot.res), fitspecies.func.GAU)) # 1hr 8 mins


# split the GAM predictions from the point estimates from each replicate
gau.pred.res <- lapply(gau.bs.gams, function(x)x[[1]])
gau.est.res <- lapply(gau.bs.gams, function(x)x[[2]])


# put GAM output list into a dataframe
gau.bs.gams.df <- data.frame(matrix(unlist(gau.pred.res), nrow=100, byrow = FALSE))

# save (to avoid having to re-run)
write.csv(gau.bs.gams.df, file="Output/Results/Trends/Bootstraps/gau.bs.gams.df.csv")


# put Estimate output list into a dataframe
gau.bs.est.df <- data.frame(matrix(unlist(gau.est.res), nrow=7, byrow = FALSE))
gau.bs.est.df <- gau.bs.est.df[ ,seq(2,length(gau.bs.est.df), by=2)]
gau.bs.est.df$year <- c("2010","2011","2013","2014","2016","2018","2020")

# save estimates
write.csv(gau.bs.est.df, file="Output/Results/Trends/Bootstraps/gau.bs.est.df.csv")

# load GAMs
#gau.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/gau.bs.gams.df.csv")
#gau.bs.gams.df <- gau.bs.gams.df[ ,-1]


### The bootstrap CIs are pretty wild, so below are some diagnostics

# extract highest annual point estimate from each replicate
maxEst <- unlist(lapply(gau.est.res, function(x) max(x$Estimate)))
maxEstdf <- data.frame(rep = 1:2000,
                       estimate = maxEst,
                       label="Maximum (any year)")

plot(maxEstdf$rep, maxEstdf$estimate)

# other diagnostic plots for maximum estimates per replicate
par(mfrow=c(2,2))
density(maxEst)
hist(maxEst) # 1400
boxplot(maxEst) # 2100
stripchart(maxEst) # 2000
plot(sort(maxEst), type="l") # 1700


# extract all estimates from each replicate
est10 <- unlist(lapply(gau.est.res, function(x) x$Estimate[1]))
est11 <- unlist(lapply(gau.est.res, function(x) x$Estimate[2]))
est13 <- unlist(lapply(gau.est.res, function(x) x$Estimate[3]))
est14 <- unlist(lapply(gau.est.res, function(x) x$Estimate[4]))
est16 <- unlist(lapply(gau.est.res, function(x) x$Estimate[5]))
est18 <- unlist(lapply(gau.est.res, function(x) x$Estimate[6]))
est20 <- unlist(lapply(gau.est.res, function(x) x$Estimate[7]))

# create dataframe for them all
est10df <- data.frame(rep=1:2000, estimate=est10, label="2010")
est11df <- data.frame(rep=1:2000, estimate=est11, label="2011")
est13df <- data.frame(rep=1:2000, estimate=est13, label="2013")
est14df <- data.frame(rep=1:2000, estimate=est14, label="2014")
est16df <- data.frame(rep=1:2000, estimate=est16, label="2016")
est18df <- data.frame(rep=1:2000, estimate=est18, label="2018")
est20df <- data.frame(rep=1:2000, estimate=est20, label="2020")

alldf <- rbind(est10df,est11df,est13df,est14df,est16df,est18df,est20df)

ggplot(alldf, aes(x=rep, y=estimate, group=label, colour=label))+
  geom_point()+
  theme_bw()

### Based on all of the above, it is slightly tricky to decide.  Some of the diagnostic plots (boxplot, stripchart) the cutoff point would be around 2000, whereas the histogram and the ordered plot suggest lower than that.  The ordered plot suggests ~1700, and this figure sits in between the others. Therefore I think this is an appropriate compromise.

# This just changes all values above 1700 to 0, so that they can't feature in the 95% and 85% quantiles
gau.bs.gams.df.sub <- data.frame(apply(gau.bs.gams.df, 2, function(x) ifelse(x >1700, 0, x)))
str(stm.bs.gams.df.sub)



# test what percentage of the last BS gam estimates are above the first. This is to test for a significant trend 
trend.df <- data.frame(apply(gau.bs.gams.df.sub,2,trendFunc))
table(trend.df)
# 96% of replicates suggest negative trend

# extract 2.5, 7.5, 50, 92.5, and 97.5 quantiles
gau.bs_ints <- data.frame(apply(gau.bs.gams.df.sub, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
gau.bs_ints <- gau.bs_ints %>% rownames_to_column("quant")

# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
gau.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(gau.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(gau.bs_ints[2,2:101]),
                         Q50 = as.numeric(gau.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(gau.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(gau.bs_ints[5,2:101]))

# change negative 2.5, 50, and 7.5 quantile values to 0
gau.quants$Q2.5[gau.quants$Q2.5 < 0] <- 0
gau.quants$Q7.5[gau.quants$Q7.5 < 0] <- 0
gau.quants$Q50[gau.quants$Q50 < 0] <- 0

# save quantiles
write.csv(gau.quants, "Output/Results/Plots/BS_quants/gau.quants.csv")

# load quantiles
gau.quants <- read.csv("Output/Results/Plots/BS_quants/gau.quants.csv")
gau.quants <- gau.quants[ ,-1]


### plots

# 95% colour
gau_plot_95 <- plot95fun(gau.dat,gau.quants,"Gaur","Ind","Individual abundance",3000)

# 95 & 85% black and white
gau_plot_85_gr <- plot85Grfun(gau.dat,gau.quants,"Gaur","Ind","Individual abundance",3000)

# because there are years with no observation I will make those points different and make the plot manually

# add grouping level to annual estimates so that I can change the point type for 2013 and 2020
gau.dat$Label2 <- ifelse(gau.dat$N==0,1,0)
gau.dat$Label2 <- as.factor(gau.dat$Label2)

# change negative 50% quantile values to 0
gau.quants <- gau.quants %>% mutate(Q50 = ifelse(Q50 < 0, 0, Q50))


# plot 95% with ribbons, faded points, and hollow points for no obs
gau.95 <- ggplot() +
    ylim(0,3250)+
    geom_point(data=gau.dat[gau.dat$Label=="Ind",], aes(x=Year, y=N, group=Label2, shape=Label2), 
               size=2, colour = "grey50") +
    scale_shape_manual(values=c(16,1))+
    geom_errorbar(data=gau.dat[gau.dat$Label=="Ind",],aes(x=Year, ymin=n_lcl, ymax=n_ucl),width=0.1, colour = "grey50") +
    geom_ribbon(data=gau.quants, aes(x=year, ymin=Q2.5, ymax=Q97.5), fill="grey60", alpha=0.3) + # 95% CIs
    geom_line(data=gau.quants, aes(x=year, y=Q50), size=1) +
    scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2020)) +
    ggtitle("Gaur") +
    ylab("Abundance") +
    xlab("Year") +
    theme(legend.position="none",
          panel.background = element_blank(), 
          axis.line = element_line(colour = "grey20"),
          axis.title = element_text(colour = "grey20"),
          axis.text = element_text(colour = "grey20"), 
          axis.ticks = element_line(colour = "grey20"),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          plot.title = element_text(margin = margin(t = 0, r = 0, b = 20, l = 0)))


    ## PIG ####

pig.dat <- select(pig.dat, -X)
pig.dat$year <- as.numeric(pig.dat$year)
str(pig.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
pig.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  pig.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), data=pig.dat[pig.dat$Label=="Ind",])
  pig.aic.res[paste0("df", dfval)] <- summary(pig.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(pig.aic.res) # df = 3


## Get confidence intervals from the bootstrapped replicates

# Function for each species. For more details see this section in the YCG section above. 

# function to fit a GAM to each replicate of boot.res
fitspecies.func.PIG <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="PIG",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model
  try(detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 60, key = "hr",cutpoints = c(0,5,12,16,25,32,45,60)))
  
  
  # extract estimates
  estimates <- detfunc$dht$individuals$N[1:7, 1:2]
  estimates <- estimates %>% dplyr::rename(Year = Label) 
  estimates$Year <- as.numeric(estimates$Year)
  
  # fit a GAM & predict
  gamfit <- gam(Estimate ~ s(Year, df=3), family=gaussian(link="identity"), data = estimates)
  newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100))
  pred <- predict.Gam(gamfit, newdata = newdata, type = "response")
  gampred <- data.frame(pred)
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(pig.bs.gams <- lapply(1:length(boot.res), fitspecies.func.PIG)) #48mins

# put output list into a dataframe
pig.bs.gams.df <- data.frame(matrix(unlist(pig.bs.gams), nrow=100, byrow = FALSE))

# save (to avoid having to re-run)
write.csv(pig.bs.gams.df, file="Output/Results/Trends/Bootstraps/pig.bs.gams.df.csv")

# load
#pig.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/pig.bs.gams.df.csv")
#pig.bs.gams.df <- pig.bs.gams.df[ ,-1]


# test what percentage of the last BS gam estimates are above the first. This is to test for a significant trend 
trend.df <- data.frame(apply(pig.bs.gams.df,2,trendFunc))
table(trend.df)
# 97% of replicates suggest negative trend

# extract 2.5, 7.5, 50, 92.5, and 97.5 quantiles
pig.bs_ints <- data.frame(apply(pig.bs.gams.df, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
pig.bs_ints <- pig.bs_ints %>% rownames_to_column("quant")

# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
pig.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(pig.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(pig.bs_ints[2,2:101]),
                         Q50 = as.numeric(pig.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(pig.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(pig.bs_ints[5,2:101]))

# save quantiles
write.csv(pig.quants, "Output/Results/Plots/BS_quants/pig.quants.csv")

# load quantiles
pig.quants <- read.csv("Output/Results/Plots/BS_quants/pig.quants.csv")
pig.quants <- pig.quants[ ,-1]


### plots

# 95% colour
pig_plot_95 <- plot95fun(pig.dat,pig.quants,"Wild pig","Ind","Individual abundance",7000)

# 95 & 85% black and white
pig_plot_85_gr <- plot85Grfun(pig.dat,pig.quants,"Wild pig","Ind","Individual abundance",7000)

# 95% black and white
pig_plot_95_gr <- plot95Grfun(pig.dat,pig.quants,"Wild pig","Ind","Individual abundance",7000)



    ## GPF ####

gpf.dat <- select(gpf.dat, -X)
gpf.dat$year <- as.numeric(gpf.dat$year)
str(gpf.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
gpf.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  gpf.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), data=gpf.dat[gpf.dat$Label=="Ind",])
  gpf.aic.res[paste0("df", dfval)] <- summary(gpf.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(gpf.aic.res) # df = 3



## Get confidence intervals from the bootstrapped replicates

# Function for each species. For more details see this section in the YCG section above. 

# function to fit a GAM to each replicate of boot.res
fitspecies.func.GPF <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="GPF",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model
  try(detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 80, key = "hn"))
  
  
  # extract estimates
  estimates <- detfunc$dht$individuals$N[1:7, 1:2]
  estimates <- estimates %>% dplyr::rename(Year = Label) 
  estimates$Year <- as.numeric(estimates$Year)
  
  # fit a GAM & predict
  gamfit <- gam(Estimate ~ s(Year, df=3), family=gaussian(link="identity"), data = estimates)
  newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100))
  pred <- predict.Gam(gamfit, newdata = newdata, type = "response")
  gampred <- data.frame(pred)
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(gpf.bs.gams <- lapply(1:length(boot.res), fitspecies.func.GPF)) # 31 mins

# put output list into a dataframe
gpf.bs.gams.df <- data.frame(matrix(unlist(gpf.bs.gams), nrow=100, byrow = FALSE))

# save (to avoid having to re-run)
write.csv(gpf.bs.gams.df, file="Output/Results/Trends/Bootstraps/gpf.bs.gams.df.csv")

# load
#gpf.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/gpf.bs.gams.df.csv")
#gpf.bs.gams.df <- gpf.bs.gams.df[ ,-1]


# test what percentage of the last BS gam estimates are above the first. This is to test for a significant trend 
trend.df <- data.frame(apply(gpf.bs.gams.df,2,trendFunc))
table(trend.df)
# 99% of replicates suggest positive trend


# extract 2.5, 7.5, 50, 92.5, and 97.5 quantiles
gpf.bs_ints <- data.frame(apply(gpf.bs.gams.df, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
gpf.bs_ints <- gpf.bs_ints %>% rownames_to_column("quant")

# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
gpf.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(gpf.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(gpf.bs_ints[2,2:101]),
                         Q50 = as.numeric(gpf.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(gpf.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(gpf.bs_ints[5,2:101]))

# save quantiles
write.csv(gpf.quants, "Output/Results/Plots/BS_quants/gpf.quants.csv")

# load quantiles
gpf.quants <- read.csv("Output/Results/Plots/BS_quants/gpf.quants.csv")
gpf.quants <- gpf.quants[ , -1]

### plots

# 95% colour
gpf_plot_95 <- plot95fun(gpf.dat,gpf.quants,"Green peafowl","Ind","Individual abundance",3000)

# 95 & 85% black and white
gpf_plot_85_gr <- plot85Grfun(gpf.dat,gpf.quants,"Green peafowl","Ind","Individual abundance",3000)

# 95% black and white
gpf_plot_95_gr <- plot95Grfun(gpf.dat,gpf.quants,"Green peafowl","Ind","Individual abundance",3000)



    ## RMJ ####

rmj.dat <- select(rmj.dat, -X)
rmj.dat$year <- as.numeric(rmj.dat$year)
str(rmj.dat)

# fit gams to real estiamtes (with varying degrees of freedom)
rmj.aic.res <- c(df1=NA, df2=NA, df3=NA)

for(dfval in 1:3){
  rmj.gam.df <- gam(Estimate ~ s(Year, df=dfval), family=gaussian(link="identity"), 
                    data=rmj.dat[rmj.dat$Label=="Ind",])
  rmj.aic.res[paste0("df", dfval)] <- summary(rmj.gam.df)$aic }

# Choose the fit with lowest AIC:
which.min(rmj.aic.res) # df = 3


## Get confidence intervals from the bootstrapped replicates

# Function for each species. For more details see this section in the YCG section above. 

# For the annual DF models for RMJ, the majority are uniform, so that is what I will use in the function. The majority of trunc distances are 60 so that is what will be used.

# function to fit a GAM to each replicate of boot.res
fitspecies.func.RMJ <- function(bootrep){
  
  
  ## Extract the bootstrap replicate data for the desired species:
  repData <- boot.res[[bootrep]]$repData
  repDataSpecies <- repData[repData$species=="RED",]
  repDataSpecies <- as.data.frame(repDataSpecies)
  sampleInfo <- boot.res[[bootrep]]$sampleInfo
  
  
  # create sample.table using SampleInfo
  sample.table <- data.frame(Sample.Label = rep(sampleInfo$Transect,times=7),
                             Region.Label = rep(c("2010","2011","2013","2014","2016","2018","2020"),
                                                each=nrow(sampleInfo)),
                             Effort = c(sampleInfo[,3],sampleInfo[,4],sampleInfo[,5],
                                          sampleInfo[,6],sampleInfo[,7],sampleInfo[,8],sampleInfo[,9]))
  
  
  # create obs.table
  obs.table <- data.frame(object =  repDataSpecies$object,
                          Region.Label =  repDataSpecies$stratum,
                          Sample.Label =  repDataSpecies$Transect)
  
  # create region.table
    # no need to create this now as the full.region.table is in the correct format
  
  
  ## fit the detection function model
  try(detfunc <- ds(repDataSpecies, region.table = full.region.table, 
                sample.table = sample.table, obs.table = obs.table,
                truncation = 60, key = "hn", adjustment= "cos", 
                cutpoints = c(0,7,14,21,28,35,42,49,56,60)))
  
  
  # extract estimates
  estimates <- detfunc$dht$individuals$N[1:7, 1:2]
  estimates <- estimates %>% dplyr::rename(Year = Label) 
  estimates$Year <- as.numeric(estimates$Year)
  
  # fit a GAM & predict
  gamfit <- gam(Estimate ~ s(Year, df=3), family=gaussian(link="identity"), data = estimates)
  newdata <- data.frame(Year = seq(from=2010, to=2020, length.out = 100))
  pred <- predict.Gam(gamfit, newdata = newdata, type = "response")
  gampred <- data.frame(pred)
  
}

## Call fitspecies.func to fit GAM to all replicates in boot.res
system.time(rmj.bs.gams <- lapply(1:length(boot.res), fitspecies.func.RMJ))

# put output list into a dataframe
rmj.bs.gams.df <- data.frame(matrix(unlist(rmj.bs.gams), nrow=100, byrow = FALSE))

# save (to avoid having to re-run)
write.csv(rmj.bs.gams.df, file="Output/Results/Trends/Bootstraps/rmj.bs.gams.df.csv")

# load
rmj.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/rmj.bs.gams.df.csv")

# test what percentage of the last BS gam estimates are above the first. This is to test for a significant trend 
trend.df <- data.frame(apply(rmj.bs.gams.df,2,trendFunc))
table(trend.df)
# 100% of replicates suggest negative trend


# extract 2.5, 7.5, 50, 92.5, and 97.5 quantiles
rmj.bs_ints <- data.frame(apply(rmj.bs.gams.df, 1, quantile, probs=c(0.025, 0.075, 0.5, 0.925, 0.975)))
rmj.bs_ints <- rmj.bs_ints %>% rownames_to_column("quant")

# Put into tidy format
#BS_ints_tidy <- gather(BS_ints, key = "year", value = "confits", -quant)

# quantiles into vectors
rmj.quants <- data.frame(year = seq(from=2010, to=2020, length.out = 100),
                         Q2.5 = as.numeric(rmj.bs_ints[1,2:101]),
                         Q7.5 = as.numeric(rmj.bs_ints[2,2:101]),
                         Q50 = as.numeric(rmj.bs_ints[3,2:101]),
                         Q92.5 = as.numeric(rmj.bs_ints[4,2:101]),
                         Q97.5 = as.numeric(rmj.bs_ints[5,2:101]))

# save quantiles
write.csv(rmj.quants, "Output/Results/Plots/BS_quants/rmj.quants.csv")

# load quantiles
rmj.quants <- read.csv("Output/Results/Plots/BS_quants/rmj.quants.csv")
rmj.quants <- rmj.quants[ ,-1]

### plots

# 95% colour
rmj_plot_95 <- plot95fun(rmj.dat,rmj.quants,"Red muntjac","Ind","Individual abundance",7000)

# 95 & 85% black and white
rmj_plot_85_gr <- plot85Grfun(rmj.dat,rmj.quants,"Red muntjac","Ind","Individual abundance",7000)

# 95% black and white
rmj_plot_95_gr <- plot95Grfun(rmj.dat,rmj.quants,"Red muntjac","Ind","Individual abundance",7000)



  ## Trends & estimates ####

    
# load all BS estimates
ycg.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/ycg.bs.gams.df.csv")
ycg.bs.gams.df <- ycg.bs.gams.df %>% select(-X)
gsl.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/gsl.bs.gams.df.csv")
gsl.bs.gams.df <- gsl.bs.gams.df %>% select(-X)
ltm.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/ltm.bs.gams.df.csv")
ltm.bs.gams.df <- ltm.bs.gams.df %>% select(-X)
ptm.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/ptm.bs.gams.df.csv")
ptm.bs.gams.df <- ptm.bs.gams.df %>% select(-X)
stm.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/stm.bs.gams.df.csv")
stm.bs.gams.df <- stm.bs.gams.df %>% select(-X)
bsd.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/bsd.bs.gams.df.csv")
bsd.bs.gams.df <- bsd.bs.gams.df %>% select(-X)
gau.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/gau.bs.gams.df.csv")
gau.bs.gams.df <- gau.bs.gams.df %>% select(-X)
pig.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/pig.bs.gams.df.csv")
pig.bs.gams.df <- pig.bs.gams.df %>% select(-X)
gpf.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/gpf.bs.gams.df.csv")
gpf.bs.gams.df <- gpf.bs.gams.df %>% select(-X)
rmj.bs.gams.df <- read.csv("Output/Results/Trends/Bootstraps/rmj.bs.gams.df.csv")
rmj.bs.gams.df <- rmj.bs.gams.df %>% select(-X)

# load in quantiles for all spcies
ycg.quants <- read.csv("Output/Results/Plots/BS_quants/ycg.quants.csv")
bsd.quants <- read.csv("Output/Results/Plots/BS_quants/bsd.quants.csv")
#btg.quants <- read.csv("Output/Results/Plots/BS_quants/btg.quants.csv")
gau.quants <- read.csv("Output/Results/Plots/BS_quants/gau.quants.csv")
gpf.quants <- read.csv("Output/Results/Plots/BS_quants/gpf.quants.csv")
gsl.quants <- read.csv("Output/Results/Plots/BS_quants/gsl.quants.csv")
ltm.quants <- read.csv("Output/Results/Plots/BS_quants/ltm.quants.csv")
pig.quants <- read.csv("Output/Results/Plots/BS_quants/pig.quants.csv")
ptm.quants <- read.csv("Output/Results/Plots/BS_quants/ptm.quants.csv")
rmj.quants <- read.csv("Output/Results/Plots/BS_quants/rmj.quants.csv")
stm.quants <- read.csv("Output/Results/Plots/BS_quants/stm.quants.csv")


    # Trends assessments ####
      # Bootstrap method ####


### This method looks at each bootstrap replicate and checks whether the last estimate is higher (positive) or lower (negative) than the first estimate. If 95% of them agree, then you can say that there is a positive/negative trend

# put all into a list
trend.list <- list(ycg.bs.gams.df,gsl.bs.gams.df,ltm.bs.gams.df,ptm.bs.gams.df,stm.bs.gams.df,
                     bsd.bs.gams.df,gau.bs.gams.df,pig.bs.gams.df,gpf.bs.gams.df,rmj.bs.gams.df)

# name elements
names(trend.list) <- c("ycg","gsl","ltm","ptm","stm","bsd","gau","pig","gpf","rmj")

# apply the trend function to each element
trend.df <- as.data.frame(lapply(trend.list,calcfun))
table(trend.df)

# manually create dataframe as can't get the bastard pivot_longer() to work
trend.df.all <- data.frame(species=c("ycg","gsl","ltm","ptm","stm","bsd","gau","pig","gpf","rmj"),
                           positive_btsrp = c(1783, 1084, 580,  1924, 0, 1749,78,70,1971,0),
                           negative_btsrp = c(217,  916,  1420, 69,   1988, 251, 1922,1930,29,2000))

trend.df.all$total <- apply(trend.df.all[,2:3],1,sum)
trend.df.all$trend <- c("positive","positive","negative","positive","negative","positive","negative",
                        "negative","positive","negative")

trend.df.all$perc <- ifelse(trend.df.all$trend=="positive",
                            trend.df.all$positive_btsrp/trend.df.all$total*100,
                            trend.df.all$negative_btsrp/trend.df.all$total*100)  

trend.df.all$sig <- ifelse(trend.df.all$perc > 95, "Yes", "no")

# save
write.csv(trend.df.all, file="OUtput/Results/Trends/trend_significance.csv")

# this not working, so did it manually above
trend.df.l <- trend.df %>% pivot_longer(ends_with(".df"), names_to = "species", values_to = "trend") %>% 
                            pivot_longer(ends_with(".Freq"),names_to = "sp.freq", values_to = "freq")


      # Overlapping CI method ####

### This method uses the 85% CIs and sees if the 2010 lcl is higher than the 2020 ucl then there is a negative trend, and if the 2010 ucl is lower than the 2020 lcl then there is a positive trend. 

# put all into a list
trend.list <- list(ycg.quants,bsd.quants,gau.quants,gpf.quants,gsl.quants,ltm.quants,pig.quants,ptm.quants,
                   rmj.quants,stm.quants)

# name elements
names(trend.list) <- c("ycg","bsd","gau","gpf","gsl","ltm","pig","ptm","rmj","stm")

# run function over list
trends.85.list <- lapply(trend.list, trendFunc2)

#
trends.85 <- as.data.frame(do.call(rbind,trends.85.list))


    # Extract predicted estimates ####

### create results table and export to excel using the GAM predictions

# put all species quantile dataframes into a list
trend.list <- list(ycg.quants,bsd.quants,gau.quants,gpf.quants,gsl.quants,ltm.quants,pig.quants,ptm.quants,
              rmj.quants,stm.quants)

# name the list elements
names(trend.list) <- c("ycg","bsd","gau","gpf","gsl","ltm","pig","ptm","rmj","stm")

# function to create dataframe using the required data from the quantiles
extrFun <- function(x){
  
  df <- data.frame(Year = c("2010","2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"),
                   Estimate = c(x$Q50[x$year==Closest(x$year,2010)],x$Q50[x$year==Closest(x$year,2011)],
                                x$Q50[x$year==Closest(x$year,2012)],x$Q50[x$year==Closest(x$year,2013)],
                                x$Q50[x$year==Closest(x$year,2014)],x[51,5],
                                x$Q50[x$year==Closest(x$year,2016)],x$Q50[x$year==Closest(x$year,2017)],
                                x$Q50[x$year==Closest(x$year,2018)],x$Q50[x$year==Closest(x$year,2019)],
                                x$Q50[x$year==Closest(x$year,2020)]),
                   lcl_95 = c(x$Q2.5[x$year==Closest(x$year,2010)],x$Q2.5[x$year==Closest(x$year,2011)],
                              x$Q2.5[x$year==Closest(x$year,2012)],x$Q2.5[x$year==Closest(x$year,2013)],
                              x$Q2.5[x$year==Closest(x$year,2014)],x[51,3],
                              x$Q2.5[x$year==Closest(x$year,2016)],x$Q2.5[x$year==Closest(x$year,2017)],
                              x$Q2.5[x$year==Closest(x$year,2018)],x$Q2.5[x$year==Closest(x$year,2019)],
                              x$Q2.5[x$year==Closest(x$year,2020)]),
                   ucl_95 = c(x$Q97.5[x$year==Closest(x$year,2010)],x$Q97.5[x$year==Closest(x$year,2011)],
                              x$Q97.5[x$year==Closest(x$year,2012)],x$Q97.5[x$year==Closest(x$year,2013)],
                              x$Q97.5[x$year==Closest(x$year,2014)],x[51,7],
                              x$Q97.5[x$year==Closest(x$year,2016)],x$Q97.5[x$year==Closest(x$year,2017)],
                              x$Q97.5[x$year==Closest(x$year,2018)],x$Q97.5[x$year==Closest(x$year,2019)],
                              x$Q97.5[x$year==Closest(x$year,2020)]),
                   lcl_85 = c(x$Q7.5[x$year==Closest(x$year,2010)],x$Q7.5[x$year==Closest(x$year,2011)],
                              x$Q7.5[x$year==Closest(x$year,2012)],x$Q7.5[x$year==Closest(x$year,2013)],
                              x$Q7.5[x$year==Closest(x$year,2014)],x[51,4],
                              x$Q7.5[x$year==Closest(x$year,2016)],x$Q7.5[x$year==Closest(x$year,2017)],
                              x$Q7.5[x$year==Closest(x$year,2018)],x$Q7.5[x$year==Closest(x$year,2019)],
                              x$Q7.5[x$year==Closest(x$year,2020)]),
                   ucl_85 = c(x$Q92.5[x$year==Closest(x$year,2010)],x$Q92.5[x$year==Closest(x$year,2011)],
                              x$Q92.5[x$year==Closest(x$year,2012)],x$Q92.5[x$year==Closest(x$year,2013)],
                              x$Q92.5[x$year==Closest(x$year,2014)],x[51,6],
                              x$Q92.5[x$year==Closest(x$year,2016)],x$Q92.5[x$year==Closest(x$year,2017)],
                              x$Q92.5[x$year==Closest(x$year,2018)],x$Q92.5[x$year==Closest(x$year,2019)],
                              x$Q92.5[x$year==Closest(x$year,2020)]))
  

}

# Apply function to the list
trend.estimates.list <- lapply(trend.list, extrFun)

# rbind the resulting dataframes together
trend.estimates <- do.call(rbind,trend.estimates.list)

# add species column
trend.estimates$Species <- rep(c("ycg","bsd","gau","gpf","gsl","ltm","pig","ptm","rmj","stm"),each=11)

# remove annoying rownames
rownames(trend.estimates) <- c()

# reorder
trend.estimates <- trend.estimates %>% select(Species,Year,Estimate,lcl_95,ucl_95,lcl_85,ucl_85)

# save 
write.csv(trend.estimates, file="Output/Results/Trends/trend_estimates.csv")




  ## Plot all ####
    ## 95% CIs only (black and white), bold points, dotted line CIs ####

# final trend plots for all species

# all species  
ycg_plot_95_gr + bsd_plot_95_gr + gsl_plot_95_gr + ltm_plot_95_gr + ptm_plot_95_gr + stm_plot_95_gr + 
  btg_plot + gau_plot_95_gr + pig_plot_95_gr + rmj_plot_95_gr + gpf_plot_95_gr


# final trend plots for primates only
final_trend_plot_prims <- ycg_plot_95_gr + bsd_plot_95_gr + gsl_plot_95_gr + ltm_plot_95_gr + 
                          ptm_plot_95_gr + stm_plot_95_gr

# remove y axis labels from plots 2,3,5,6
final_trend_plot_prims[[2]] <- final_trend_plot_prims[[2]] + theme(axis.title.y = element_blank())
final_trend_plot_prims[[3]] <- final_trend_plot_prims[[3]] + theme(axis.title.y = element_blank())
final_trend_plot_prims[[5]] <- final_trend_plot_prims[[5]] + theme(axis.title.y = element_blank())
final_trend_plot_prims[[6]] <- final_trend_plot_prims[[6]] + theme(axis.title.y = element_blank())

# remove x axis labels from plots 1,2,3
final_trend_plot_prims[[1]] <- final_trend_plot_prims[[1]] + theme(axis.title.x = element_blank())
final_trend_plot_prims[[2]] <- final_trend_plot_prims[[2]] + theme(axis.title.x = element_blank())
final_trend_plot_prims[[3]] <- final_trend_plot_prims[[3]] + theme(axis.title.x = element_blank())

# adjust angle of x axis tick mark text and remove from plots 1,2,3
final_trend_plot_prims[[1]] <- final_trend_plot_prims[[1]] + theme(axis.text.x = element_blank())
final_trend_plot_prims[[2]] <- final_trend_plot_prims[[2]] + theme(axis.text.x = element_blank())
final_trend_plot_prims[[3]] <- final_trend_plot_prims[[3]] + theme(axis.text.x = element_blank())
final_trend_plot_prims[[4]] <- final_trend_plot_prims[[4]] + theme(axis.text.x = element_text(angle=45,hjust = 1))
final_trend_plot_prims[[5]] <- final_trend_plot_prims[[5]] + theme(axis.text.x = element_text(angle=45,hjust = 1))
final_trend_plot_prims[[6]] <- final_trend_plot_prims[[6]] + theme(axis.text.x = element_text(angle=45,hjust = 1))

# adjust size of tick mark numbers
final_trend_plot_prims[[1]] <- final_trend_plot_prims[[1]] + theme(axis.text = element_text(size=12))
final_trend_plot_prims[[2]] <- final_trend_plot_prims[[2]] + theme(axis.text = element_text(size=12))
final_trend_plot_prims[[3]] <- final_trend_plot_prims[[3]] + theme(axis.text = element_text(size=12))
final_trend_plot_prims[[4]] <- final_trend_plot_prims[[4]] + theme(axis.text = element_text(size=12))
final_trend_plot_prims[[5]] <- final_trend_plot_prims[[5]] + theme(axis.text = element_text(size=12))
final_trend_plot_prims[[6]] <- final_trend_plot_prims[[6]] + theme(axis.text = element_text(size=12))

# adjust size of axis titles
final_trend_plot_prims[[1]] <- final_trend_plot_prims[[1]] + theme(axis.title = element_text(size=15))
final_trend_plot_prims[[4]] <- final_trend_plot_prims[[4]] + theme(axis.title = element_text(size=15))
final_trend_plot_prims[[5]] <- final_trend_plot_prims[[5]] + theme(axis.title = element_text(size=15))
final_trend_plot_prims[[6]] <- final_trend_plot_prims[[6]] + theme(axis.title = element_text(size=15))


# save
ggsave("Output/Results/Plots/final_trend_plot_prims.png", final_trend_plot_prims,
       width = 30, height = 20, units = "cm", dpi = 300)





# final trend plots for non-primates
final_trend_plot_ungs <- btg_plot + gau_plot_95_gr + pig_plot_95_gr + rmj_plot_95_gr + gpf_plot_95_gr

# remove y axis labels from plots 2,3,5
final_trend_plot_ungs[[2]] <- final_trend_plot_ungs[[2]] + theme(axis.title.y = element_blank())
final_trend_plot_ungs[[3]] <- final_trend_plot_ungs[[3]] + theme(axis.title.y = element_blank())
final_trend_plot_ungs[[5]] <- final_trend_plot_ungs[[5]] + theme(axis.title.y = element_blank())

# remove x axis labels from plots 1,2
final_trend_plot_ungs[[1]] <- final_trend_plot_ungs[[1]] + theme(axis.title.x = element_blank())
final_trend_plot_ungs[[2]] <- final_trend_plot_ungs[[2]] + theme(axis.title.x = element_blank())

# adjust angle of x axis tick mark text and remove from plots 1,2
final_trend_plot_ungs[[1]] <- final_trend_plot_ungs[[1]] + theme(axis.text.x = element_blank())
final_trend_plot_ungs[[2]] <- final_trend_plot_ungs[[2]] + theme(axis.text.x = element_blank())
final_trend_plot_ungs[[3]] <- final_trend_plot_ungs[[3]] + theme(axis.text.x = element_text(angle=45,hjust = 1))
final_trend_plot_ungs[[4]] <- final_trend_plot_ungs[[4]] + theme(axis.text.x = element_text(angle=45,hjust = 1))
final_trend_plot_ungs[[5]] <- final_trend_plot_ungs[[5]] + theme(axis.text.x = element_text(angle=45,hjust = 1))

# adjust size of tick mark numbers
final_trend_plot_ungs[[1]] <- final_trend_plot_ungs[[1]] + theme(axis.text = element_text(size=12))
final_trend_plot_ungs[[2]] <- final_trend_plot_ungs[[2]] + theme(axis.text = element_text(size=12))
final_trend_plot_ungs[[3]] <- final_trend_plot_ungs[[3]] + theme(axis.text = element_text(size=12))
final_trend_plot_ungs[[4]] <- final_trend_plot_ungs[[4]] + theme(axis.text = element_text(size=12))
final_trend_plot_ungs[[5]] <- final_trend_plot_ungs[[5]] + theme(axis.text = element_text(size=12))

# adjust size of axis titles
final_trend_plot_ungs[[1]] <- final_trend_plot_ungs[[1]] + theme(axis.title = element_text(size=15))
final_trend_plot_ungs[[3]] <- final_trend_plot_ungs[[3]] + theme(axis.title = element_text(size=15))
final_trend_plot_ungs[[4]] <- final_trend_plot_ungs[[4]] + theme(axis.title = element_text(size=15))
final_trend_plot_ungs[[5]] <- final_trend_plot_ungs[[5]] + theme(axis.title = element_text(size=15))


# save
ggsave("Output/Results/Plots/final_trend_plot_ungs.png", final_trend_plot_ungs,
       width = 30, height = 20, units = "cm", dpi = 300)







    ## 85% and 95% CIs (black and white), no points, ribbon CIs ####

# first make all the individual plots
ycg.p <- plot95Grfun3(ycg.quants,"Yellow-cheeked crested gibbon","Grp","Group abundance", 950)
gsl.p <- plot95Grfun3(gsl.quants,"Germain's silver langur","Grp","Group abundance", 1600)
ltm.p <- plot95Grfun3(ltm.quants,"Long-tailed macaque","Grp","Group abundance", 800)
ptm.p <- plot95Grfun3(ptm.quants,"Pig-tailed macaque","Grp","Group abundance", 2100)
stm.p <- plot95Grfun3(stm.quants,"Stump-tailed macaque","Grp","Group abundance", 550)
bsd.p <- plot95Grfun3(bsd.quants,"Black-shanked douc langur","Grp","Group abundance", 11550)
btg.p <- btg_plot
gau.p <- plot95Grfun3(gau.quants,"Gaur","Ind","Individual abundance", 1250)
pig.p <- plot95Grfun3(pig.quants,"Wild pig","Ind","Individual abundance", 3600)
gpf.p <- plot95Grfun3(gpf.quants,"Green peafowl","Ind","Individual abundance", 1950)
rmj.p <- plot95Grfun3(rmj.quants,"Red muntjac","Ind","Individual abundance", 4800)


## Primates

# make plot grid
prims.p <- ycg.p+ gsl.p + bsd.p + ltm.p + ptm.p + stm.p 

# remove y axis labels from plots 2,3,5,6
prims.p[[2]] <- prims.p[[2]] + theme(axis.title.y = element_blank())
prims.p[[3]] <- prims.p[[3]] + theme(axis.title.y = element_blank())
prims.p[[5]] <- prims.p[[5]] + theme(axis.title.y = element_blank())
prims.p[[6]] <- prims.p[[6]] + theme(axis.title.y = element_blank())

# remove x axis labels and text from plots 1,2,3 and increase size for 4,5,6
prims.p[[1]] <- prims.p[[1]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
prims.p[[2]] <- prims.p[[2]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
prims.p[[3]] <- prims.p[[3]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
prims.p[[4]] <- prims.p[[4]] + theme(axis.title.x = element_text(size=14))
prims.p[[5]] <- prims.p[[5]] + theme(axis.title.x = element_text(size=14))
prims.p[[6]] <- prims.p[[6]] + theme(axis.title.x = element_text(size=14))

# adjust y axis label size
prims.p[[1]] <- prims.p[[1]] + theme(axis.title.y = element_text(size=14))
prims.p[[4]] <- prims.p[[4]] + theme(axis.title.y = element_text(size=14))

# adjust size of tick mark numbers
prims.p[[1]] <- prims.p[[1]] + theme(axis.text = element_text(size=11))
prims.p[[2]] <- prims.p[[2]] + theme(axis.text = element_text(size=11))
prims.p[[3]] <- prims.p[[3]] + theme(axis.text = element_text(size=11))
prims.p[[4]] <- prims.p[[4]] + theme(axis.text = element_text(size=11))
prims.p[[5]] <- prims.p[[5]] + theme(axis.text = element_text(size=11))
prims.p[[6]] <- prims.p[[6]] + theme(axis.text = element_text(size=11))


ggsave("Output/Results/Plots/Trends/trend_plot_prims.png", prims.p,
       width = 30, height = 20, units = "cm", dpi = 300)


## non-primates

# make plot grid
ungs.p <- btg.p+ gau.p + pig.p + rmj.p + gpf.p 

# remove y axis labels from plots 2,3,5,6
ungs.p[[2]] <- ungs.p[[2]] + theme(axis.title.y = element_blank())
ungs.p[[3]] <- ungs.p[[3]] + theme(axis.title.y = element_blank())
ungs.p[[5]] <- ungs.p[[5]] + theme(axis.title.y = element_blank())

# remove x axis labels and text from plots 1,2,3 and increase size for 4,5,6
ungs.p[[1]] <- ungs.p[[1]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p[[2]] <- ungs.p[[2]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p[[3]] <- ungs.p[[3]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p[[4]] <- ungs.p[[4]] + theme(axis.title.x = element_text(size=14))
ungs.p[[5]] <- ungs.p[[5]] + theme(axis.title.x = element_text(size=14))

# adjust y axis label size
ungs.p[[1]] <- ungs.p[[1]] + theme(axis.title.y = element_text(size=14))
ungs.p[[4]] <- ungs.p[[4]] + theme(axis.title.y = element_text(size=14))

# adjust size of tick mark numbers
ungs.p[[1]] <- ungs.p[[1]] + theme(axis.text = element_text(size=11))
ungs.p[[2]] <- ungs.p[[2]] + theme(axis.text = element_text(size=11))
ungs.p[[3]] <- ungs.p[[3]] + theme(axis.text = element_text(size=11))
ungs.p[[4]] <- ungs.p[[4]] + theme(axis.text = element_text(size=11))
ungs.p[[5]] <- ungs.p[[5]] + theme(axis.text = element_text(size=11))


ggsave("Output/Results/Plots/Trends/trend_plot_ungs.png", ungs.p,
       width = 30, height = 20, units = "cm", dpi = 300)

    ## 85% and 95% CIs (black and white), faded points, ribbon CIs ####

# first make all the individual plots
ycg.p2 <- plot95Grfun2(ycg.dat, ycg.quants,"Yellow-cheeked crested gibbon","Grp","Group abundance", 1550)
gsl.p2 <- plot95Grfun2(gsl.dat, gsl.quants,"Germain's silver langur","Grp","Group abundance", 2500)
ltm.p2 <- plot95Grfun2(ltm.dat, ltm.quants,"Long-tailed macaque","Grp","Group abundance", 1950)
ptm.p2 <- plot95Grfun2(ptm.dat, ptm.quants,"Pig-tailed macaque","Grp","Group abundance", 2500)
stm.p2 <- plot95Grfun2(stm.dat, stm.quants,"Stump-tailed macaque","Grp","Group abundance", 650)
bsd.p2 <- plot95Grfun2(bsd.dat, bsd.quants,"Black-shanked douc langur","Grp","Group abundance", 13000)
btg.p2 <- btg_plot2
gau.p2 <- plot95Grfun2(gau.dat, gau.quants,"Gaur","Ind","Individual abundance", 3250)
pig.p2 <- plot95Grfun2(pig.dat, pig.quants,"Wild pig","Ind","Individual abundance", 6200)
gpf.p2 <- plot95Grfun2(gpf.dat, gpf.quants,"Green peafowl","Ind","Individual abundance", 2800)
rmj.p2 <- plot95Grfun2(rmj.dat, rmj.quants,"Red muntjac","Ind","Individual abundance", 6150)


## Primates

# make plot grid
prims.p2 <- ycg.p2+ gsl.p2 + bsd.p2 + ltm.p2 + ptm.p2 + stm.p2 

# remove y axis labels from plots 2,3,5,6
prims.p2[[2]] <- prims.p2[[2]] + theme(axis.title.y = element_blank())
prims.p2[[3]] <- prims.p2[[3]] + theme(axis.title.y = element_blank())
prims.p2[[5]] <- prims.p2[[5]] + theme(axis.title.y = element_blank())
prims.p2[[6]] <- prims.p2[[6]] + theme(axis.title.y = element_blank())

# remove x axis labels and text from plots 1,2,3 and increase size for 4,5,6
prims.p2[[1]] <- prims.p2[[1]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
prims.p2[[2]] <- prims.p2[[2]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
prims.p2[[3]] <- prims.p2[[3]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
prims.p2[[4]] <- prims.p2[[4]] + theme(axis.title.x = element_text(size=14))
prims.p2[[5]] <- prims.p2[[5]] + theme(axis.title.x = element_text(size=14))
prims.p2[[6]] <- prims.p2[[6]] + theme(axis.title.x = element_text(size=14))

# adjust y axis label size
prims.p2[[1]] <- prims.p2[[1]] + theme(axis.title.y = element_text(size=14))
prims.p2[[4]] <- prims.p2[[4]] + theme(axis.title.y = element_text(size=14))

# adjust size of tick mark numbers
prims.p2[[1]] <- prims.p2[[1]] + theme(axis.text = element_text(size=11))
prims.p2[[2]] <- prims.p2[[2]] + theme(axis.text = element_text(size=11))
prims.p2[[3]] <- prims.p2[[3]] + theme(axis.text = element_text(size=11))
prims.p2[[4]] <- prims.p2[[4]] + theme(axis.text = element_text(size=11))
prims.p2[[5]] <- prims.p2[[5]] + theme(axis.text = element_text(size=11))
prims.p2[[6]] <- prims.p2[[6]] + theme(axis.text = element_text(size=11))


ggsave("Output/Results/Plots/Trends/trend_plot_prims_points.png", prims.p2,
       width = 30, height = 20, units = "cm", dpi = 300)


## non-primates

# make plot grid
ungs.p2 <- btg.p2 + gau.p2 + pig.p2 + rmj.p2 + gpf.p2 

# remove y axis labels from plots 2,3,5,6
ungs.p2[[2]] <- ungs.p2[[2]] + theme(axis.title.y = element_blank())
ungs.p2[[3]] <- ungs.p2[[3]] + theme(axis.title.y = element_blank())
ungs.p2[[5]] <- ungs.p2[[5]] + theme(axis.title.y = element_blank())

# remove x axis labels and text from plots 1,2,3 and increase size for 4,5,6
ungs.p2[[1]] <- ungs.p2[[1]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p2[[2]] <- ungs.p2[[2]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p2[[3]] <- ungs.p2[[3]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p2[[4]] <- ungs.p2[[4]] + theme(axis.title.x = element_text(size=14))
ungs.p2[[5]] <- ungs.p2[[5]] + theme(axis.title.x = element_text(size=14))

# adjust y axis label size
ungs.p2[[1]] <- ungs.p2[[1]] + theme(axis.title.y = element_text(size=14))
ungs.p2[[4]] <- ungs.p2[[4]] + theme(axis.title.y = element_text(size=14))

# adjust size of tick mark numbers
ungs.p2[[1]] <- ungs.p2[[1]] + theme(axis.text = element_text(size=11))
ungs.p2[[2]] <- ungs.p2[[2]] + theme(axis.text = element_text(size=11))
ungs.p2[[3]] <- ungs.p2[[3]] + theme(axis.text = element_text(size=11))
ungs.p2[[4]] <- ungs.p2[[4]] + theme(axis.text = element_text(size=11))
ungs.p2[[5]] <- ungs.p2[[5]] + theme(axis.text = element_text(size=11))


ggsave("Output/Results/Plots/Trends/trend_plot_ungs_points.png", ungs.p2,
       width = 30, height = 20, units = "cm", dpi = 300)

    ## 95% CIs (black and white), faded points, ribbon CIs ####


# first make all the individual plots
ycg.95 <- plot95Grfun4(ycg.dat, ycg.quants,"Southern yellow-cheeked crested gibbon","Ind","Abundance", 3100)
gsl.95 <- plot95Grfun4(gsl.dat, gsl.quants,"Germain's silver langur","Ind","Abundance", 10450)
ltm.95 <- plot95Grfun4(ltm.dat, ltm.quants,"Long-tailed macaque","Ind","Abundance", 8450)
ptm.95 <- plot95Grfun4(ptm.dat, ptm.quants,"Northern pig-tailed macaque","Ind","Abundance", 8750)
stm.95 <- stm.95
bsd.95 <- plot95Grfun4(bsd.dat, bsd.quants,"Black-shanked douc","Ind","Abundance", 50500)
btg.95 <- btg.95
gau.95 <- gau.95
pig.95 <- plot95Grfun4(pig.dat, pig.quants,"Wild pig","Ind","Abundance", 6200)
gpf.95 <- plot95Grfun4(gpf.dat, gpf.quants,"Green Peafowl","Ind","Abundance", 2800)
rmj.95 <- plot95Grfun4(rmj.dat, rmj.quants,"Northern red muntjac","Ind","Abundance", 6150)


      # Arranged by primates/non-primates ####

## Primates

## Primates

# make plot grid
prims.p2 <- ycg.95+ gsl.95 + bsd.95 + ltm.95 + ptm.95 + stm.95 

# remove y axis labels from plots 2,3,5,6
prims.p2[[2]] <- prims.p2[[2]] + theme(axis.title.y = element_blank())
prims.p2[[3]] <- prims.p2[[3]] + theme(axis.title.y = element_blank())
prims.p2[[5]] <- prims.p2[[5]] + theme(axis.title.y = element_blank())
prims.p2[[6]] <- prims.p2[[6]] + theme(axis.title.y = element_blank())

# remove x axis labels and text from plots 1,2,3 and increase size for 4,5,6
prims.p2[[1]] <- prims.p2[[1]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
prims.p2[[2]] <- prims.p2[[2]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
prims.p2[[3]] <- prims.p2[[3]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
prims.p2[[4]] <- prims.p2[[4]] + theme(axis.title.x = element_text(size=14))
prims.p2[[5]] <- prims.p2[[5]] + theme(axis.title.x = element_text(size=14))
prims.p2[[6]] <- prims.p2[[6]] + theme(axis.title.x = element_text(size=14))

# adjust y axis label size
prims.p2[[1]] <- prims.p2[[1]] + theme(axis.title.y = element_text(size=14))
prims.p2[[4]] <- prims.p2[[4]] + theme(axis.title.y = element_text(size=14))

# adjust size of tick mark numbers
prims.p2[[1]] <- prims.p2[[1]] + theme(axis.text = element_text(size=11))
prims.p2[[2]] <- prims.p2[[2]] + theme(axis.text = element_text(size=11))
prims.p2[[3]] <- prims.p2[[3]] + theme(axis.text = element_text(size=11))
prims.p2[[4]] <- prims.p2[[4]] + theme(axis.text = element_text(size=11))
prims.p2[[5]] <- prims.p2[[5]] + theme(axis.text = element_text(size=11))
prims.p2[[6]] <- prims.p2[[6]] + theme(axis.text = element_text(size=11))

ggsave("Output/Results/Plots/Trends/trend_plot95_prims_points.png", prims.p2,
       width = 30, height = 20, units = "cm", dpi = 300)


## non-primates


## 2 rows, 3 columns

# make plot grid
ungs.p2 <- btg.95 + gau.95 + pig.95 + rmj.95 + gpf.95 

# remove y axis labels from plots 2,3,5,6
ungs.p2[[2]] <- ungs.p2[[2]] + theme(axis.title.y = element_blank())
ungs.p2[[3]] <- ungs.p2[[3]] + theme(axis.title.y = element_blank())
ungs.p2[[5]] <- ungs.p2[[5]] + theme(axis.title.y = element_blank())

# remove x axis labels and text from plots 1,2 and increase size for 3,4,5
ungs.p2[[1]] <- ungs.p2[[1]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p2[[2]] <- ungs.p2[[2]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p2[[3]] <- ungs.p2[[3]] + theme(axis.title.x = element_text(size=14))
ungs.p2[[4]] <- ungs.p2[[4]] + theme(axis.title.x = element_text(size=14))
ungs.p2[[5]] <- ungs.p2[[5]] + theme(axis.title.x = element_text(size=14))

# adjust y axis label size
ungs.p2[[1]] <- ungs.p2[[1]] + theme(axis.title.y = element_text(size=14))
ungs.p2[[4]] <- ungs.p2[[4]] + theme(axis.title.y = element_text(size=14))

# adjust size of tick mark numbers
ungs.p2[[1]] <- ungs.p2[[1]] + theme(axis.text = element_text(size=11))
ungs.p2[[2]] <- ungs.p2[[2]] + theme(axis.text = element_text(size=11))
ungs.p2[[3]] <- ungs.p2[[3]] + theme(axis.text = element_text(size=11))
ungs.p2[[4]] <- ungs.p2[[4]] + theme(axis.text = element_text(size=11))
ungs.p2[[5]] <- ungs.p2[[5]] + theme(axis.text = element_text(size=11))

ggsave("Output/Results/Plots/Trends/trend_plot95_ungsWIDE_points.png", ungs.p2,
       width = 30, height = 20, units = "cm", dpi = 300)



## 3 rows, 2 columns

# make plot grid
ungs.p3 <- btg.95 + gau.95 + pig.95 + rmj.95 + gpf.95 + plot_layout(ncol=2)

# remove y axis labels from plots 2,4
ungs.p3[[2]] <- ungs.p3[[2]] + theme(axis.title.y = element_blank())
ungs.p3[[4]] <- ungs.p3[[4]] + theme(axis.title.y = element_blank())

# remove x axis labels and text from plots 1,2,3 and increase size for 4,5
ungs.p3[[1]] <- ungs.p3[[1]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p3[[2]] <- ungs.p3[[2]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p3[[3]] <- ungs.p3[[3]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
ungs.p3[[4]] <- ungs.p3[[4]] + theme(axis.title.x = element_text(size=14))
ungs.p3[[5]] <- ungs.p3[[5]] + theme(axis.title.x = element_text(size=14))

# adjust y axis label size
ungs.p3[[1]] <- ungs.p3[[1]] + theme(axis.title.y = element_text(size=14))
ungs.p3[[3]] <- ungs.p3[[3]] + theme(axis.title.y = element_text(size=14))
ungs.p3[[5]] <- ungs.p3[[5]] + theme(axis.title.y = element_text(size=14))

# adjust size of tick mark numbers
ungs.p3[[1]] <- ungs.p3[[1]] + theme(axis.text = element_text(size=11))
ungs.p3[[2]] <- ungs.p3[[2]] + theme(axis.text = element_text(size=11))
ungs.p3[[3]] <- ungs.p3[[3]] + theme(axis.text = element_text(size=11))
ungs.p3[[4]] <- ungs.p3[[4]] + theme(axis.text = element_text(size=11))
ungs.p3[[5]] <- ungs.p3[[5]] + theme(axis.text = element_text(size=11))

ggsave("Output/Results/Plots/Trends/trend_plot95_ungsLONG_points.png", ungs.p3,
       width = 30, height = 20, units = "cm", dpi = 300)

      # Arranged by trend ####

### Stable or increasing species 

# make plot grid
stab.p2 <- ptm.95 + gpf.95 + ycg.95 + bsd.95 + ltm.95 + gsl.95   

# remove y axis labels from plots 2,3,5,6
stab.p2[[2]] <- stab.p2[[2]] + theme(axis.title.y = element_blank())
stab.p2[[3]] <- stab.p2[[3]] + theme(axis.title.y = element_blank())
stab.p2[[5]] <- stab.p2[[5]] + theme(axis.title.y = element_blank())
stab.p2[[6]] <- stab.p2[[6]] + theme(axis.title.y = element_blank())

# remove x axis labels and text from plots 1,2,3 and increase size for 4,5,6
stab.p2[[1]] <- stab.p2[[1]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
stab.p2[[2]] <- stab.p2[[2]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
stab.p2[[3]] <- stab.p2[[3]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
stab.p2[[4]] <- stab.p2[[4]] + theme(axis.title.x = element_text(size=14))
stab.p2[[5]] <- stab.p2[[5]] + theme(axis.title.x = element_text(size=14))
stab.p2[[6]] <- stab.p2[[6]] + theme(axis.title.x = element_text(size=14))

# adjust y axis label size
stab.p2[[1]] <- stab.p2[[1]] + theme(axis.title.y = element_text(size=14))
stab.p2[[4]] <- stab.p2[[4]] + theme(axis.title.y = element_text(size=14))

# adjust size of tick mark numbers
stab.p2[[1]] <- stab.p2[[1]] + theme(axis.text = element_text(size=11))
stab.p2[[2]] <- stab.p2[[2]] + theme(axis.text = element_text(size=11))
stab.p2[[3]] <- stab.p2[[3]] + theme(axis.text = element_text(size=11))
stab.p2[[4]] <- stab.p2[[4]] + theme(axis.text = element_text(size=11))
stab.p2[[5]] <- stab.p2[[5]] + theme(axis.text = element_text(size=11))
stab.p2[[6]] <- stab.p2[[6]] + theme(axis.text = element_text(size=11))


ggsave("Output/Results/Plots/Trends/trend_plot95_stable.png", stab.p2,
       width = 30, height = 20, units = "cm", dpi = 300)

### decreasing species

# make plot grid
decr.p2 <- stm.95 + btg.95 + gau.95 + rmj.95 + pig.95

# remove y axis labels from plots 2,3,5
decr.p2[[2]] <- decr.p2[[2]] + theme(axis.title.y = element_blank())
decr.p2[[3]] <- decr.p2[[3]] + theme(axis.title.y = element_blank())
decr.p2[[5]] <- decr.p2[[5]] + theme(axis.title.y = element_blank())

# remove x axis labels and text from plots 1 and 2 and increase size for 3, 4 and 5
decr.p2[[1]] <- decr.p2[[1]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
decr.p2[[2]] <- decr.p2[[2]] + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())
decr.p2[[3]] <- decr.p2[[3]] + theme(axis.title.x = element_text(size=14))
decr.p2[[4]] <- decr.p2[[4]] + theme(axis.title.x = element_text(size=14))
decr.p2[[5]] <- decr.p2[[5]] + theme(axis.title.x = element_text(size=14))

# adjust y axis label size for 1 and 4
decr.p2[[1]] <- decr.p2[[1]] + theme(axis.title.y = element_text(size=14))
decr.p2[[4]] <- decr.p2[[4]] + theme(axis.title.y = element_text(size=14))

# adjust size of tick mark numbers
decr.p2[[1]] <- decr.p2[[1]] + theme(axis.text = element_text(size=11))
decr.p2[[2]] <- decr.p2[[2]] + theme(axis.text = element_text(size=11))
decr.p2[[3]] <- decr.p2[[3]] + theme(axis.text = element_text(size=11))
decr.p2[[4]] <- decr.p2[[4]] + theme(axis.text = element_text(size=11))
decr.p2[[5]] <- decr.p2[[5]] + theme(axis.text = element_text(size=11))

ggsave("Output/Results/Plots/Trends/trend_plot95_decrease.png", decr.p2,
       width = 30, height = 20, units = "cm", dpi = 300)



## plot the two plots together
stab.p2 / decr.p2


# arboreal primates: YCG, BSD, GSL, 
# Semi-arb species:  GPF, PTM, LTM, STM, 
# ground based spec: BTG, GAU, RMJ, PIG

all.p <- ycg.95+bsd.95+gsl.95+gpf.95+ptm.95+ltm.95+stm.95+btg.95+gau.95+rmj.95+pig.95 + plot_layout(ncol=4)




ggsave("Output/Results/Plots/Trends/trend_plot95_prims_points.png", prims.p2,
       width = 30, height = 20, units = "cm", dpi = 300)