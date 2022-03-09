setwd("/Users/ariglenn/Desktop/Ratledge Work/IMRproject")

library(raster) #extract rasters
library(MetBrewer) #color palet
library(sf) #helps switch from different type of shape files
library(tidyverse)
library(dplyr)
library(fixest)
library(glmnet) #machine learning weighted regression
library(ggplot2) #plotter
library(ggthemes) #makes figures look nice
library(reshape2) #reshapes data
library(utils)
library(foreign) #interpret different types of data frames
library(MCPanel) #matrix completion

#might be able to remove
#library(rgdal)
#library(lfe)

#communties from all over Uganda - anonimize communities by jittering them by 5km in radius
#means village could appear up to 5km radius away from electric grid line
#2km buffer, anything within 2km is included
#3km buffer - #X control villages, #Y treatment villages
#X and Y values will differ between 2km, 3km, and 4km
#
#4km buffer - 

# loading treatment and control data
TC_3km <- read.csv("treat_control_3km.csv") #csv is in the drive

#takes data from all rows w/2014 as a data value under year by View(TC_2km)
#964 total = 888 control and 76 treated
TC_3km_2014sub <- subset(TC_3km, year == "2014")

#creates a matrix with 964 rows, 1 col, filled w/NA
#number of rows = length of data = number of villages
numVillages <- 825 #964 #825
count <- matrix(NA,numVillages,1)
#labels each row from 1 to 964
count[,1] <- c(1:numVillages)

#assigns each row in the data set of the 2014 subset with a unqiue identifies (1-964)
TC_3km_2014sub_count <- cbind(count,TC_3km_2014sub)

#deletes everything except the longitude and latitude values from the data set (from columns 4 and 5)
#needs to be long, lat
TC_3km_2014sub_count <- TC_3km_2014sub_count[,c(5,4)]

# here is how to load and prep infant mortality data i just sent

#creates a box of latitude and longitude that roughly focus on Uganda
#this is just setting the extent of the raster we want to extract
# i.e. just lat and long box around uganda, not the whole world!
e <- as(extent(29, 34, -2, 5), 'SpatialPolygons')

#ensures that all data is in the same coordinate system
crs(e) <- "EPSG:4326"

#imports TIF file that has many layers with stack
#infant data from ghdx data set
ghdx_infant <- stack("IHME_LMICS_U5M_2000_2017_D_INFANT_MEAN_Y2019M10D16.TIF")

#gets infant data from the defined Uganda box (e)
ghdx_trial_ug <- crop(ghdx_infant, e)

#ensures that all data is in the same coordinate system
pr_ghdx_trial_ug <- projectRaster(ghdx_trial_ug,crs="+proj=longlat +datum=NAD27")

#from the TIF file, extracts the points related to the TC_2km_2014sub_count
#gets infant data related to each village (coordinates for each TC_2km_2014sub_count -- 964 villages)
rast_pr_ghdx_trial_ug <- raster::extract(pr_ghdx_trial_ug, TC_3km_2014sub_count)

#ensures that all data is in the same data frame
rast_pr_ghdx_trial_ug <- as.data.frame(rast_pr_ghdx_trial_ug)

#renames each column for years 2000-2017
colnames(rast_pr_ghdx_trial_ug) <- c("imr_2000", "imr_2001", "imr_2002", "imr_2003", "imr_2004", 
                                     "imr_2005", "imr_2006", "imr_2007", "imr_2008", "imr_2009", 
                                     "imr_2010", "imr_2011", "imr_2012", "imr_2013", "imr_2014", 
                                     "imr_2015", "imr_2016", "imr_2017")

#matrix that contains imr for years 2000-2017 for each 964 villages in Uganda


# this is the website for the nighttime light data. 
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/YGIVCD

count <- matrix(NA,numVillages,1) # just creating a unique identifier
count[,1] <- c(1:numVillages)

# now starting the raster extract for 2003 to 2016

NTL_2003 <- raster("LongNTL_2003.tif")
ug_NTL_2003 <- crop(NTL_2003, e)
pr_NTL_2003 <- projectRaster(ug_NTL_2003,crs="+proj=longlat +datum=NAD27")
# plot just to spot check that there are not import problems
plot(pr_NTL_2003)

rast_val_NTL_2003_TC <- raster::extract(pr_NTL_2003, TC_3km_2014sub_count)
rast_val_NTL_2003_TC <- as.data.frame(rast_val_NTL_2003_TC)
colnames(rast_val_NTL_2003_TC) <- "ntl_2003"

# this save the data so we don't have to do the extraction again
write.csv(rast_val_NTL_2003_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2003_TC.csv")

# this combines the extracted data to the original data frame
TC_3km_2014sub_count <- cbind(TC_3km_2014sub, rast_val_NTL_2003_TC)

# 2004
NTL_2004 <- raster("LongNTL_2004.tif")
ug_NTL_2004 <- crop(NTL_2004, e)
pr_NTL_2004 <- projectRaster(ug_NTL_2004,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2004)

rast_val_NTL_2004_TC <- raster::extract(pr_NTL_2004, TC_3km_2014sub_count)
rast_val_NTL_2004_TC <- as.data.frame(rast_val_NTL_2004_TC)
colnames(rast_val_NTL_2004_TC) <- "ntl_2004"

write.csv(rast_val_NTL_2004_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2004_TC.csv")

TC_3km_2014sub_count <- cbind(TC_3km_2014sub, rast_val_NTL_2004_TC)


# 2005
NTL_2005 <- raster("LongNTL_2005.tif")
ug_NTL_2005 <- crop(NTL_2005, e)
pr_NTL_2005 <- projectRaster(ug_NTL_2005,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2005)

rast_val_NTL_2005_TC <- raster::extract(pr_NTL_2005, TC_2km_2014sub_count)
rast_val_NTL_2005_TC <- as.data.frame(rast_val_NTL_2005_TC)
colnames(rast_val_NTL_2005_TC) <- "ntl_2005"

write.csv(rast_val_NTL_2005_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2005_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2005_TC)


# 2006
NTL_2006 <- raster("LongNTL_2006.tif")
ug_NTL_2006 <- crop(NTL_2006, e)
pr_NTL_2006 <- projectRaster(ug_NTL_2006,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2006)

rast_val_NTL_2006_TC <- raster::extract(pr_NTL_2006, TC_2km_2014sub_count)
rast_val_NTL_2006_TC <- as.data.frame(rast_val_NTL_2006_TC)
colnames(rast_val_NTL_2006_TC) <- "ntl_2006"

write.csv(rast_val_NTL_2006_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2006_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2006_TC)


# 2007
NTL_2007 <- raster("LongNTL_2007.tif")
ug_NTL_2007 <- crop(NTL_2007, e)
pr_NTL_2007 <- projectRaster(ug_NTL_2007,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2007)

rast_val_NTL_2007_TC <- raster::extract(pr_NTL_2007, TC_2km_2014sub_count)
rast_val_NTL_2007_TC <- as.data.frame(rast_val_NTL_2007_TC)
colnames(rast_val_NTL_2007_TC) <- "ntl_2007"

write.csv(rast_val_NTL_2007_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2007_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2007_TC)


# 2008
NTL_2008 <- raster("LongNTL_2008.tif")
ug_NTL_2008 <- crop(NTL_2008, e)
pr_NTL_2008 <- projectRaster(ug_NTL_2008,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2008)

rast_val_NTL_2008_TC <- raster::extract(pr_NTL_2008, TC_2km_2014sub_count)
rast_val_NTL_2008_TC <- as.data.frame(rast_val_NTL_2008_TC)
colnames(rast_val_NTL_2008_TC) <- "ntl_2008"

write.csv(rast_val_NTL_2008_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2008_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2008_TC)


# 2009
NTL_2009 <- raster("LongNTL_2009.tif")
ug_NTL_2009 <- crop(NTL_2009, e)
pr_NTL_2009 <- projectRaster(ug_NTL_2009,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2009)

rast_val_NTL_2009_TC <- raster::extract(pr_NTL_2009, TC_2km_2014sub_count)
rast_val_NTL_2009_TC <- as.data.frame(rast_val_NTL_2009_TC)
colnames(rast_val_NTL_2009_TC) <- "ntl_2009"

write.csv(rast_val_NTL_2009_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2009_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2009_TC)


# 2010
NTL_2010 <- raster("LongNTL_2010.tif")
ug_NTL_2010 <- crop(NTL_2010, e)
pr_NTL_2010 <- projectRaster(ug_NTL_2010,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2010)

rast_val_NTL_2010_TC <- raster::extract(pr_NTL_2010, TC_2km_2014sub_count)
rast_val_NTL_2010_TC <- as.data.frame(rast_val_NTL_2010_TC)
colnames(rast_val_NTL_2010_TC) <- "ntl_2010"

write.csv(rast_val_NTL_2010_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2010_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2010_TC)


# 2011
NTL_2011 <- raster("LongNTL_2011.tif")
ug_NTL_2011 <- crop(NTL_2011, e)
pr_NTL_2011 <- projectRaster(ug_NTL_2011,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2011)

rast_val_NTL_2011_TC <- raster::extract(pr_NTL_2011, TC_2km_2014sub_count)
rast_val_NTL_2011_TC <- as.data.frame(rast_val_NTL_2011_TC)
colnames(rast_val_NTL_2011_TC) <- "ntl_2011"

write.csv(rast_val_NTL_2011_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2011_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2011_TC)


# 2012
NTL_2012 <- raster("LongNTL_2012.tif")
ug_NTL_2012 <- crop(NTL_2012, e)
pr_NTL_2012 <- projectRaster(ug_NTL_2012,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2012)

rast_val_NTL_2012_TC <- raster::extract(pr_NTL_2012, TC_2km_2014sub_count)
rast_val_NTL_2012_TC <- as.data.frame(rast_val_NTL_2012_TC)
colnames(rast_val_NTL_2012_TC) <- "ntl_2012"

write.csv(rast_val_NTL_2012_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2012_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2012_TC)


# 2013
NTL_2013 <- raster("LongNTL_2013.tif")
ug_NTL_2013 <- crop(NTL_2013, e)
pr_NTL_2013 <- projectRaster(ug_NTL_2013,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2013)

rast_val_NTL_2013_TC <- raster::extract(pr_NTL_2013, TC_2km_2014sub_count)
rast_val_NTL_2013_TC <- as.data.frame(rast_val_NTL_2013_TC)
colnames(rast_val_NTL_2013_TC) <- "ntl_2013"

write.csv(rast_val_NTL_2013_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2013_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2013_TC)


# 2014
NTL_2014 <- raster("LongNTL_2014.tif")
ug_NTL_2014 <- crop(NTL_2014, e)
pr_NTL_2014 <- projectRaster(ug_NTL_2014,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2014)

rast_val_NTL_2014_TC <- raster::extract(pr_NTL_2014, TC_2km_2014sub_count)
rast_val_NTL_2014_TC <- as.data.frame(rast_val_NTL_2014_TC)
colnames(rast_val_NTL_2014_TC) <- "ntl_2014"

write.csv(rast_val_NTL_2014_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2014_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2014_TC)


# 2015
NTL_2015 <- raster("LongNTL_2015.tif")
ug_NTL_2015 <- crop(NTL_2015, e)
pr_NTL_2015 <- projectRaster(ug_NTL_2015,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2015)

rast_val_NTL_2015_TC <- raster::extract(pr_NTL_2015, TC_2km_2014sub_count)
rast_val_NTL_2015_TC <- as.data.frame(rast_val_NTL_2015_TC)
colnames(rast_val_NTL_2015_TC) <- "ntl_2015"

write.csv(rast_val_NTL_2015_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2015_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2015_TC)


# 2016
NTL_2016 <- raster("LongNTL_2016.tif")
ug_NTL_2016 <- crop(NTL_2016, e)
pr_NTL_2016 <- projectRaster(ug_NTL_2016,crs="+proj=longlat +datum=NAD27")
plot(pr_NTL_2016)

rast_val_NTL_2016_TC <- raster::extract(pr_NTL_2016, TC_2km_2014sub_count)
rast_val_NTL_2016_TC <- as.data.frame(rast_val_NTL_2016_TC)
colnames(rast_val_NTL_2016_TC) <- "ntl_2016"

write.csv(rast_val_NTL_2016_TC, "rast_val_NTL_csvFiles/rast_val_NTL_2016_TC.csv")

raster_TC_2km_2014sub <- cbind(raster_TC_2km_2014sub, rast_val_NTL_2016_TC)





raster_TC_2km_2014sub_zeros <- raster_TC_2km_2014sub

raster_TC_2km_2014sub_zeros[is.na(raster_TC_2km_2014sub_zeros)] <- 0