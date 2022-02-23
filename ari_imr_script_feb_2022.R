
setwd("/Users/ariglenn/Desktop/Ratledge Work/IMRproject")

#not sure we still need all these, but...
library(raster)
library(MetBrewer)
library(sf)
library(tidyverse)
library(rgdal)
library(fixest)
library(ggplot2)
library(utils)
library(dplyr)
library(fixest)
library(reshape2)
library(lfe)
library(foreign) 
#library(MCPanel) 
library(glmnet)
library(ggthemes)


# loading treatment and control data
TC_2km <- read.csv("TC_2km.csv") #csv is in the drive

TC_2km_2014sub <- subset(TC_2km, year == "2014") #964 total = 888 control and 76 treated

count <- matrix(NA,964,1)
count[,1] <- c(1:964)

TC_2km_2014sub_count <- cbind(count,TC_2km_2014sub)

TC_2km_2014sub_count <- TC_2km_2014sub_count[,c(5,4)] #needs to be long, lat

# here is how to load and prep infant mortality data i just sent

e <- as(extent(29, 34, -2, 5), 'SpatialPolygons') #this is just setting the extent of the raster we want to extract
# i.e. just lat and long box around uganda, not the whole world!
crs(e) <- "EPSG:4326"

ghdx_infant <- stack("IHME_LMICS_U5M_2000_2017_D_INFANT_MEAN_Y2019M10D16.TIF")
ghdx_trial_ug <- crop(ghdx_infant, e)
pr_ghdx_trial_ug <- projectRaster(ghdx_trial_ug,crs="+proj=longlat +datum=NAD27")

rast_pr_ghdx_trial_ug <- raster::extract(pr_ghdx_trial_ug, TC_2km_2014sub_count)
rast_pr_ghdx_trial_ug <- as.data.frame(rast_pr_ghdx_trial_ug)

colnames(rast_pr_ghdx_trial_ug) <- c("imr_2000", "imr_2001", "imr_2002", "imr_2003", "imr_2004", 
                                     "imr_2005", "imr_2006", "imr_2007", "imr_2008", "imr_2009", 
                                     "imr_2010", "imr_2011", "imr_2012", "imr_2013", "imr_2014", 
                                     "imr_2015", "imr_2016", "imr_2017")

# i will add more later

