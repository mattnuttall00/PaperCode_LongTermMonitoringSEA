## Rachel's code file for bootstrap resampling.

## Load libraries and data ####

library('tidyverse')
library('Distance')
library('data.table')
library('gam')
library('patchwork')
library('DescTools')

## load data (allData). This is the master data from the CDS analysis
load("./Output/Data/KSWS_MASTER.Rdata")
## In this dataframe, stratum == year.

## Remove T20 from allData and obs.table
allData <- allData[allData$obs.transect != 20, ]
obs.table <- obs.table[obs.table$Sample.Label != 20, ]

## ===================================================
## Code in the if(F){} brackets below is for copy-and-paste: it isn't run when the R code is sourced.
## ===================================================
if(F){
        ## Running:
        source("Rfunc.R")   ## sources R code and libraries, and loads the allData data-frame

        ## Create effort data "eff_dat.2":
        eff_dat.2 <- createEffortData.func(filename="Input/sample.table.csv")  ## *** Rachel - changed .txt to .csv

        ## ---------------------------------------------------------------------------------------------------------------
        ## Run bootstrap: this generates 500 replicates using method 2.
        ## Method 2 = stop when overall effort in each habitat across years reaches or exceeds that in the real data.
        boot.res <- bootstrap.func(Nrep=500, method=2)

        ## ---------------------------------------------------------------------------------------------------------------
        ## Plot the number of sightings among the bootstrap replicates and real data and print a mean summary:
        bootstrap.plot(boot.res, habitat="All", plotwhat="Sightings")

        ## Plot (say) the number of sightings in Dense habitat:
        bootstrap.plot(boot.res, habitat="Dense", plotwhat="Sightings")

        ## Plot the effort per year instead of the number of sightings:
        bootstrap.plot(boot.res, habitat="All", plotwhat="Effort")

        ## Loop to plot number of sightings overall, and in each of the habitats:
        for(hab in c("All", "Dense", "Open")){## , "Nonf")){  *** Rachel - removed Nonf category 20/7/2020
                cat("\n", hab, "habitat:\n")
                print(bootstrap.plot(boot.res, habitat=hab, plotwhat="Sightings"))
                readline("Enter for next plot...")
        }

        ## Loop to plot amount of effort overall, and in each of the habitats:
        for(hab in c("All", "Dense", "Open")){## , "Nonf")){  *** Rachel - removed Nonf category 20/7/2020
                cat("\n", hab, "habitat:\n")
                print(bootstrap.plot(boot.res, habitat=hab, plotwhat="Effort"))
                readline("Enter for next plot...")
        }

        ## ---------------------------------------------------------------------------------------------------------------
        ## The required bootstrap replicates are in the list boot.res.
        ## For bootstrap replicate i, the sightings data is in boot.res[[i]]$repData, and the sample-table
        ## is in boot.res[[i]]$sampleInfo.  Check out some replicates as follows:
        ##
        head(boot.res[[5]]$repData)
        boot.res[[5]]$sampleInfo
        ## The column "Transect" in both objects corresponds to the bootstrap transect label.  To see which
        ## transect each sighting was originally on in the real data, use TransectOrig in the sampleInfo data,
        ## and obs.transectOrig in the repData.  These should be consistent across the two data-frames:
        ## e.g. all instances of Transect=12 in either data-frame should have the same result for TransectOrig
        ## and for obs.transectOrig.  Check this out as so:
        boot.res[[5]]$sampleInfo[12,]
        ## The line above tells you which TransectOrig corresponds to Transect 12 in the bootstrap replicate.
        ## The line below tells you all values of obs.transectOrig in the resampled data corresponding to
        ## Transect label 12: it should be a single answer the same as TransectOrig above.
        unique(boot.res[[5]]$repData$obs.transectOrig[boot.res[[5]]$repData$Transect==12])
        ## Note that the habitat type for this transect is most easily found from boot.res[[5]]$sampleInfo[12,]:
        ## it will be easiest to stick with the bootstrap datasets for all downstream analyses; don't go back to the
        ## original data for some lookups as this could get confusing.

        ## ---------------------------------------------------------------------------------------------------------------
        ## Each bootstrap replicate contains data for all species.  To isolate data for a single species, say BSD,
        ##  for bootstrap replicate 5, use:
        boot.res[[5]]$repData[boot.res[[5]]$repData$species=="BSD", ]
        ## Or more likely, you'll be writing a function something like this:
        fitspecies.func <- function(bootrep){
                ## Extract the bootstrap replicate data for the desired species:
                repData <- boot.res[[bootrep]]$repData
                repDataSpecies <- repData[repData$species=="BSD",]
                sampleInfo <- boot.res[[bootrep]]$sampleInfo
                ## fit the model and return the results here
        }
        ## Call fitspecies.func using something like:
        lapply(1:length(boot.res), fitspecies.func)

}
## ===================================================
## End of the if(F) copy-and-paste example code.
## ===================================================


####################################################
## FUNCTIONS BELOW
####################################################

createEffortData.func <- function(filename){
        ## createEffortData.func: Matt's code to create the effort data-frame "eff_dat.2" put into a function
        ## for easy use by Rachel.
        ## filename set to Matt's directory structure by default.

        ## load in data that has effort
        eff_dat <- read.csv(filename, header = TRUE)
        eff_dat <- eff_dat[,-1]
        eff_dat <- eff_dat[eff_dat$Sample.Label !=20,] # remove T20
        eff_dat$Region.Label <- as.character(eff_dat$Region.Label)

        ## rename the high and low strata transects with just the year
        eff_dat$Region.Label[eff_dat$Region.Label=="2011_H"] <- "2011"
        eff_dat$Region.Label[eff_dat$Region.Label=="2011_L"] <- "2011"
        eff_dat$Region.Label[eff_dat$Region.Label=="2014_H"] <- "2014"
        eff_dat$Region.Label[eff_dat$Region.Label=="2014_L"] <- "2014"
        eff_dat$Region.Label[eff_dat$Region.Label=="2016_H"] <- "2016"
        eff_dat$Region.Label[eff_dat$Region.Label=="2016_L"] <- "2016"

        ## "widen" the data so that each year has its own column and "effort" is the value
        eff_dat.2 <- pivot_wider(eff_dat, id_cols = c("Sample.Label","Region.Label"), names_from = Region.Label,
                                 values_from = Effort)

        ## I need to allocate a broad habitat class to each transect.  Obviously this will not capture the fine-scale habitat hetergeneity within transects, but the only reason we are doing this is to account for the spatial trends in species for the bootstrapping (see Rachel's e-mail).  By splitting the transects by broad habitat category, we will be capturing all of the spatial trends amongst species. i.e. all of the species have either low density in dense forest (D) areas and high density in open forest (O) areas, or vice versa.  Therefore by splitting the transects by those categories, we are are capturing the transects that have either high or low density for each species.  I will allocate habitat category based on the land cover layer in QGIS. The results are as follows:

        ## Dense forest - 1:8, 10:12, 15, 18:19, 28, 32:33, 36,
        ## Open forest - 9, 13:14, 16:17, 21:27, 30:31, 34:35, 37:40
        ## Non-forest - 20, 29

        ## add habitat category ## *** Rachel - changed the last remaining "nonf" transect (29) to "dense" below:
        dense <- c(1:8, 10:12, 15, 18:19, 28, 29, 32:33, 36) ## *** Rachel - added transect 29 here 20/7/2020
        open <- c(9, 13:14, 16:17, 21:27, 30:31, 34:35, 37:40)
        ## nonf <- 29 *** Rachel - deleted this 20/7/2020 as above

        eff_dat.2 <- eff_dat.2 %>% mutate(Habitat = case_when(Sample.Label %in% dense ~ "D",
                                          Sample.Label %in% open ~ "O"))## *** Rachel: deleted next line 20/7/2020
        ## Sample.Label %in% nonf ~ "NF"))

        ## re-name columns and re-order
        eff_dat.2 <- eff_dat.2 %>% dplyr::select(Habitat, Transect=Sample.Label,
                                                 eff_10="2010",eff_11="2011",eff_13="2013",eff_14="2014",
                                                 eff_16="2016",eff_18="2018",eff_20="2020")

        ## Return finalised effort data:
        return(eff_dat.2)

}

##########################################################

bootstrap.func <- function(obsdat=allData, effdat=eff_dat.2, Nrep=1000, method=2){
        ## bootstrap.func
        ## Rachel's bootstrapping function - 5/5/2020
        ## obsdat is the distance-sampling observations by transect: allData by default.
        ## effdat is the data on effort by transect: eff_dat.2 by default.
        ## Nrep is the number of bootstrap replicates required.
        ##
        ## If method=1, we keep resampling transects until the first time that the effort for each individual year
        ## first equals or exceeds that in the real data.
        ##
        ## If method=2, we keep resampling transects until the total effort in all years put together
        ## first equals or exceeds that in the real data. This might mean that the effort in some years is less than
        ## that in the real data.
        ##
        ## Set method=2 by default.

        ## -------------------------------------------------------------------------------------------------
        ## Housekeeping:
        if(!(method %in% c(1, 2))) stop("Method should be 1 or 2: see the blurb of bootstrap.func for explanation.")

        ## -------------------------------------------------------------------------------------------------
        ## Rename Transect in effdat as TransectOrig, so it is not confused with the bootstrap replicate transect labels.
        ## Similarly for columns object and obs.transect in obsdat.
        ## *** Rachel new for 20/7/2020: when transect 20 was removed from the data, obsdat$object was no
        ## longer the same as the row-numbers in obsdat. Instead of indexing obsdat by object number,
        ## the code is now changed to create indexOrig (the original row numbers for each object)
        ## for indexing purposes. Before transect 20 was removed, objectOrig and indexOrig were the same thing.
        names(effdat)[names(effdat)=="Transect"] <- "TransectOrig"
        names(obsdat)[names(obsdat)=="object"] <- "objectOrig"
        names(obsdat)[names(obsdat)=="obs.transect"] <- "obs.transectOrig"
        obsdat$indexOrig <- 1:nrow(obsdat) ## *** Rachel added 20/7/2020 : indexOrig gives original row numbers

        ## -------------------------------------------------------------------------------------------------
        ## Preliminary calculations of total effort per year in each habitat category:
        effnames <- c("eff_10", "eff_11", "eff_13", "eff_14", "eff_16", "eff_18", "eff_20")
        ## Create a list which contains the rows of effdat organised by habitat: these are the equivalents of Matt's
        ## df_dense, df_open, and df_nonf.  Putting them in a list means we can access them in a loop over habitats
        ## instead of having to go through each habitat separately.
        effdatList <- list(Dense=effdat[effdat$Habitat=="D",],
                           Open=effdat[effdat$Habitat=="O",])##, *** Rachel removed Nonf category 20/7/2020
                           ## Nonf=effdat[effdat$Habitat=="NF",])
        ## Total effort in Dense, Open, and NonForest habitats by year: call these "lims" to match Matt's naming,
        ## and make into a list with names Dense, Open, Nonf.
        limsList <- lapply(effdatList, function(x) colSums(x[, effnames]))
        ## To save computation later, find the limsSums: the total effort across years in each habitat:
        limsSum <- lapply(limsList, sum)

        ## -------------------------------------------------------------------------------------------------
        ## Create a reference list of which objects in obsdat belong to each transect: objectTransectRefList is a
        ## list of length 40 where the i'th element gives all object IDs in obsdat that were observed on transect i.
        ## We will use objectTransectRefList inside the bootstrapping function to quickly find the required
        ## rows of obsdat for each bootstrap replicate.
        ## *** Matt edited line below to convert object and obs.transect into objectOrig and obs.transectOrig: good.
        ## Presumably the original worked because of partial matching of column-names that has been deprecated.
        ## *** Rachel further editing it 20/7/2020 to change objectOrig to indexOrig.
        ## *** Additionally, transect 20 is no longer in effdat, but the transects are still numbered 1 to 40 with
        ## 20 missing.  We want to be able to index the TransectRefList by the transect number, so that
        ## list element 40 corresponds to transect label 40.  So replace lapply(1:nrow(effdat) in the line below with
        ## lapply(1:max(effdat$TransectOrig). Also change the name from objectTransectRefList to
        ## indexTransectRefList. The result is that indexTransectRefList[[i]] is the rows of obsdat
        ## that correspond to the transect with label i, for i=1, ..., 19, 21, ..., 40.
        indexTransectRefList <- lapply(1:max(effdat$TransectOrig),
                                        function(i) obsdat$indexOrig[obsdat$obs.transectOrig==i])

        ## -------------------------------------------------------------------------------------------------
        ## Create an inner function for a single bootstrap replicate.
        oneRep.func <- function(bootrep){
                ## oneRep.func is one bootstrap replicate. "bootrep" specifies which replicate it is: rep=1, 2, ..., Nrep.
                if(bootrep%%10==0) cat("Bootstrap replicate", bootrep, "\n")
                ## Resample transects until the desired effort is obtained, according to method 1 or method 2.
                ## First create an empty list "transectRep" with components Dense, Open, Nonf:
                ## this will contain the list of transects and effort for each habitat (like Matt's objects "first", "first_open",
                ## "first_nf").
                transectRep <- vector("list", 2)  ## *** Rachel changed 3 to 2 to remove Nonf 20/7/2020
                names(transectRep) <- c("Dense", "Open")## , "Nonf") ## *** Rachel removed Nonf 20/7/2020
                ## Now create the sample for each habitat and fill in the three data-frames inside transectRep:
                for(habitat in c("Dense", "Open")){## , "Nonf")){  ## *** Rachel removed Nonf 20/7/2020
                        ## effortReached is equivalent of Matt's "exc": changed name as we don't want to "exceed" any more
                        ## but just be greater than or equal to.
                        effortReached <- FALSE
                        while(!effortReached){
                                transectRep[[habitat]] <- rbind(transectRep[[habitat]],
                                                                sample_n(tbl = effdatList[[habitat]], size = 1))
                                ## Find the total effort in this habitat so far:
                                colsum <- colSums(transectRep[[habitat]][,effnames])
                                ## Decide if effort is reached or not according to method 1 or method 2.
                                ## Method 1: all elements of colsum should be >= the corresponding members of lims;
                                ## Method 2: the sum of colsum should be >= the corresponding limsum:
                                effortReached <- ifelse(method==1,
                                                        all(colsum >= limsList[[habitat]]),  ## method 1
                                                        (sum(colsum) >= limsSum[[habitat]])) ## method 2
                        }  ## End of effortReached while loop
                } ## End of looping over habitat dense, open, nonf.  The list transectRep is now filled for each habitat.

                ## Combine the three habitats into a single data-frame, sampleInfoRep: ## *** Rachel removed Nonf
                sampleInfoRep <- as.data.frame(rbind(transectRep$Dense, transectRep$Open))## , transectRep$Nonf))

                ## Now label the transects in sampleInfoRep as 1:nrow(sampleInfoRep): this will be the transect ID
                ## for analysing the bootstrap data, and we want it to be distinct for each replicate transect, whether
                ## or not that replicate transect is a duplicate of a transect from the original data.
                ## sampleInfoRep$Transect gives the bootstrap transect label (1:nrow(sampleInfoRep)).
                ## The original (real-data) transect that this is a copy of is in sampleInfoRep$TransectOrig.
                sampleInfoRep$Transect <- 1:nrow(sampleInfoRep)
                ## End of code for resampling transects until the desired effort is reached.

                ## -----------------------------------------------------------------------------------------------------
                ## Now assemble the bootstrap sightings data: obsdatRep.
                ## *** Rachel 20/7/2020 : modified the lines below to talk about indexes (row-numbers) instead of
                ## objects:
                ## Find all the rows corresponding to sampleInfoRep$TransectOrig using the reference list
                ## defined in indexTransectRefList at the top of the function:
                indexTransectRep <- indexTransectRefList[sampleInfoRep$TransectOrig]
                indexesRep <- unlist(indexTransectRep)

                ## Now we just need the corresponding rows of the original obsdat:
                obsdatRep <- obsdat[indexesRep,]

                ## *** Rachel 20/7/2020 - change objectTransectRep to indexTransectRep below:
                ## The lengths of the elements of indexTransectRep give us the new transect labels in the replicate data:
                ## e.g. if indexTransectRep[[5]] has 100 objects, this means that 100 rows of obsdatRep should have
                ## the transect label 5:
                obsdatRep$Transect <- rep(sampleInfoRep$Transect, times=unlist(lapply(indexTransectRep, length)))
                ## Lastly create a new column for "object" in obsdatRep which goes from 1 to nrow(obsdatRep):
                ## we need to treat each row as if it is a new object, even if it is a duplicate of another row.
                obsdatRep$object <- 1:nrow(obsdatRep)  ## *** This is OK, we can label the replicate objects 1:nrow
                ## Return the observed data and the effort data for this replicate:
                return(list(repData=obsdatRep, sampleInfo=sampleInfoRep))
        }  ## End of oneRep.func

        ## -------------------------------------------------------------------------------------------------------------
        ## Run onerep.func Nrep times and return the results.
        ## bootres contains all Nrep data-frames (observations of all species) and transect information.
        ## The corresponding items from (say) replicate 51 are in bootres[[51]]$repData and
        ## bootres[[51]]$sampleInfo.
        bootres <- lapply(1:Nrep, oneRep.func)
        ## Add the bootstrapping method to the attributes of bootres:
        attributes(bootres)$method <- method
        return(bootres)
}

##########################################################

bootstrap.plot <- function(bootres, plotwhat="Sightings", habitat="All", obsdat=allData, effdat=eff_dat.2){
        ## bootstrap.plot 9/5/2020
        ## Takes bootres, the output from bootstrap.func, and plots the bootstrap distribution compared with the
        ## real-data values by year of the requested output.
        ## plotwhat can be either "Sightings" or "Effort", respectively plotting the number of sightings per year
        ## or the effort per year.
        ## habitat can be "All", "Dense", "Open", or "Nonf". If "All" it gives the total across all habitat types.
        ##
        ## The real-data values are generated from obsdat and effdat, so these should match the choices that
        ## were used to generate bootres in bootstrap.func. This step could be made more bullet-proof, but I
        ## haven't done so because I assume that allData and eff_dat.2 are always going to be the right choices.
        ##
        ## Uses Matt's nifty ggplot approach for the plot.

        ## ------------------------------------------------------------------------------------------------------------------
        ## Housekeeping:
        if(!(plotwhat %in% c("Sightings", "Effort"))) stop("plotwhat should be either 'Sightings' or 'Effort'.")
        if(!(habitat %in% c("All","Dense","Open")))## ,"Nonf"))) ## *** Rachel removed Nonf 20/7/2020
                stop("habitat should be one of: 'All', 'Dense', or 'Open'.")## , or 'Nonf'.")  ## *** Rachel removed Nonf
        effnames <- c("eff_10", "eff_11", "eff_13", "eff_14", "eff_16", "eff_18", "eff_20") ## *** Rachel added eff_20
        ## Create a named vector, habitatCode, so we can easily reference "Dense" to "D", "Open" to "O, and
        ## "Nonf" to "NF".  For example if we want to convert "Dense" to "D", ask for habitatCode["Dense"].
        ## In general, habitatCode[habitat] gives the right thing for the habitat asked for in the function arguments.
        habitatCode <- c(All="All", Dense="D", Open="O", Nonf="NF")  ## *** Rachel - this is OK, no change needed
        ## We also need a column called Transect in the real data to match the format in the bootstrap replicates:
        obsdat$Transect <- obsdat$obs.transect

        ## ------------------------------------------------------------------------------------------------------------------
        ## Define separate data-processing functions according to whether plotwhat="Sightings" or "Effort".
        ## If bootrep=0 then perform the same operations on the real-data, obsdat and effdat.
        processSightings.func <- function(bootrep){
                ## For a single bootstrap replicate, use bootres[[bootrep]]$repData and extract a short data-frame of
                ## year and number of sightings.
                if(bootrep>0) repData <- bootres[[bootrep]]$repData
                else repData <- obsdat
                ## Focus on the habitat required: if habitat is not "All", we need to identify it using the Transect numbers
                ## in sampleInfo for that replicate.  We identify transects by their bootstrap labels (i.e. Transect, not
                ## TransectOrig) in both data-frames.
                if(habitat!="All"){
                        if(bootrep>0) sampleInfo <- bootres[[bootrep]]$sampleInfo
                        else sampleInfo <- effdat
                        transectsHab <- sampleInfo$Transect[sampleInfo$Habitat==habitatCode[habitat]]
                        repData <- repData[repData$Transect %in% transectsHab,]
                }
                ## Generate "year" by using the stratum column, but strip out any instances of _H or _L so that
                ## (for example) 2011_H becomes just 2011.
                repData$year <- gsub(pattern="_H", replacement="", x=repData$stratum)
                repData$year <- gsub(pattern="_L", replacement="", x=repData$year)
                ## The info required is just the table of repData$year within the habitat-specific repData:
                as.data.frame(table(repData$year), stringsAsFactors=FALSE)
        }
        ## ------------------------------------------------------------------------------------------------------------------
        ## Function to process the effort in each bootstrap replicate.  Supply bootrep=0 for the real-data values.
        processEffort.func <- function(bootrep){
                ## For a single bootstrap replicate, use bootres[[bootrep]]$repData and extract a short data-frame of
                ## year and amount of effort.
                if(bootrep>0) sampleInfo <- bootres[[bootrep]]$sampleInfo
                else sampleInfo <- effdat
                ## Focus on the habitat required: if habitat is not "All", subset the sampleInfo data by habitat code:
                if(habitat!="All")
                        sampleInfo <- sampleInfo[sampleInfo$Habitat == habitatCode[habitat],]
                ## The info required is just the colSums of sampleInfo by year within the habitat-specific sample data:
                effortRep <- colSums(sampleInfo[effnames])
                ## Rename years from "eff_10" to "2010" etc:
                yearnames <- gsub(pattern="eff_", replacement="20", x=names(effortRep))
                ## Return the result:
                data.frame(year=yearnames, effort=effortRep, row.names=NULL, stringsAsFactors=FALSE)
        }

        ## ------------------------------------------------------------------------------------------------------------------
        ## Create the plots
        ## ------------------------------------------------------------------------------------------------------------------
        ## Sightings plot:
        if(plotwhat=="Sightings"){
                ## Count number of sightings by year within each bootstrap replicate:
                ## the code below is clunky but gets the job done for creating a single long-form dataframe with
                ## columns "year" and "Nsightings".
                sightings.df <- NULL
                for(i in 1:length(bootres)) sightings.df <- rbind(sightings.df, processSightings.func(i))
                names(sightings.df) <- c("year", "Nsightings")

                ## Find the equivalent results for the real data by putting bootrep=0:
                original.df <- processSightings.func(0)
                names(original.df) <- c("year", "Nsightings")

                ## Plot sightings:
                print(
                      ggplot(sightings.df, aes(x=year, y=Nsightings)) +
                      geom_point() +
                      geom_point(data=original.df, aes(x=year, y=Nsightings), shape=4, size=6, colour="red") +
                      ggtitle(paste0("Method ", attributes(bootres)$method, " Sightings: ", habitat, " habitat",
                                     ifelse(habitat=="All", "s", "")))
                      )

                ## Return the summary of mean sightings by year:
                bootMeans <- unlist(lapply(split(sightings.df, sightings.df$year), function(x)mean(x$Nsightings)))
                summary.df <- original.df
                summary.df$BootstrapMean <- bootMeans[summary.df$year]
                summary.df$PercentDiff <-
                        (summary.df$BootstrapMean - summary.df$Nsightings)/summary.df$Nsightings*100
                return(summary.df)

        } ## End of plotwhat="Sightings"

        ## ------------------------------------------------------------------------------------------------------------------
        ## Effort plot:
        ## Sightings plot:
        if(plotwhat=="Effort"){
                ## Find effort by year within each bootstrap replicate:
                ## the code below creates a single long-form dataframe with columns "year" and "effort".
                effort.df <- NULL
                for(i in 1:length(bootres)) effort.df <- rbind(effort.df, processEffort.func(i))
                names(effort.df) <- c("year", "effort")

                ## Find the equivalent results for the real data by putting bootrep=0:
                original.df <- processEffort.func(0)
                names(original.df) <- c("year", "effort")

                ## Plot effort:
                print(
                      ggplot(effort.df, aes(x=year, y=effort)) +
                      geom_point() +
                      geom_point(data=original.df, aes(x=year, y=effort), shape=4, size=6, colour="blue") +
                      ggtitle(paste0("Method ", attributes(bootres)$method, " Effort: ", habitat, " habitat",
                                     ifelse(habitat=="All", "s", "")))
                      )

                ## Return the summary of mean sightings by year:
                bootMeans <- unlist(lapply(split(effort.df, effort.df$year), function(x)mean(x$effort)))
                summary.df <- original.df
                summary.df$BootstrapMean <- bootMeans[summary.df$year]
                summary.df$PercentDiff <- (summary.df$BootstrapMean - summary.df$effort)/summary.df$effort*100
                return(summary.df)
        }

}

###############################################################
