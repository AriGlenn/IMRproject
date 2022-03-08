# The is the 3 and 4 km script for the 5 penalty

# NEED TO CHANGE NAMES IN LINE 225

# This script should produce the causal estimates for DiD, MC, SC-EN, and SC-ENt. 
# It also runs pre-trend tests for DiD, and MC and SC-ENt


# Not sure all these are needed anymore, some may be legacy, but...
library(utils)
library(dplyr)
library(sf)
library(fixest)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(wesanderson)
library(fixest)
library(rgdal)
library(lfe)
library(foreign) #for importing certain file types
library(MCPanel) #for MC
library(glmnet)
library(ggthemes) 

# Add your WD here
setwd("/Users/ariglenn/Desktop/Ratledge Work/IMRproject")

ug_dist_2010 <- st_read(dsn = "Ug_Aug_2010_shapefile.shp", layer = "Ug_Aug_2010_shapefile", stringsAsFactors = F)

ug_grid_2010_shape_3km <- ug_dist_2010  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 3000) %>%
  st_transform("+init=epsg:4326")

ug_grid_2010_shape_4km <- ug_dist_2010  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 4000) %>%
  st_transform("+init=epsg:4326")

ug_dist_2013  <- st_read(dsn = "UG_2013_Shapefile_Aug30.shp", layer = "UG_2013_Shapefile_Aug30", stringsAsFactors = F)

ug_grid_2013_shape_3km <- ug_dist_2013  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 3000) %>%
  st_transform("+init=epsg:4326")

ug_grid_2013_shape_4km <- ug_dist_2013  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 4000) %>%
  st_transform("+init=epsg:4326")

ug_dist_2016  <- st_read(dsn = "Distribution_Lines_2016_Operational.shp", layer = "Distribution_Lines_2016_Operational", stringsAsFactors = F)

ug_grid_2016_shape_3km <- ug_dist_2016  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 3000) %>%
  st_transform("+init=epsg:4326")

ug_grid_2016_shape_4km <- ug_dist_2016  %>%
  st_transform("+init=epsg:3358") %>%
  st_buffer(dist = 4000) %>%
  st_transform("+init=epsg:4326")

atlas_cluster_data_w_WI_scenarios_causal <- read.csv("atlas_cluster_data_nov_12_2020.csv")


atlas_cluster_data_w_WI_scenarios_causal <- atlas_cluster_data_w_WI_scenarios_causal[,-1]

dhs_1_9_84_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_108_dhs.csv") 
colnames(dhs_1_9_84_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_216_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_240_dhs.csv")
colnames(dhs_1_9_216_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_348_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_372_dhs.csv")
colnames(dhs_1_9_348_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_480_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_504_dhs.csv")
colnames(dhs_1_9_480_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_612_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_636_dhs.csv")
colnames(dhs_1_9_612_DHSpanel)[1] <- "DHSID_year"





###########
# Here is the 3km single / non-Ensemble version # 
###########
##### CI_V0_1_7_config_id_42 #####

# here i am making the ENS df
config_3q_42_ENS <- as.data.frame(matrix(NA,30583,2))
config_3q_42_ENS[,1] <- dhs_1_9_84_DHSpanel[,1]
config_3q_42_ENS[,2] <- (dhs_1_9_84_DHSpanel[,2] + dhs_1_9_84_DHSpanel[,2] + dhs_1_9_84_DHSpanel[,2] + dhs_1_9_84_DHSpanel[,2] + dhs_1_9_84_DHSpanel[,2]) / 5

sapply(config_3q_42_ENS, class)
as.numeric(config_3q_42_ENS$X..DHSID)
colnames(config_3q_42_ENS)[1] <- "dhsID_year"
colnames(config_3q_42_ENS)[2] <- "3q_1CL_ENS_pred"

config_3q_42_ENS$DHSID <- substr(config_3q_42_ENS$dhsID_year, 1, 14)
config_3q_42_ENS$year <- substr(config_3q_42_ENS$dhsID_year, 16, 19)
config_3q_42_ENS$survey_year <- substr(config_3q_42_ENS$dhsID_year, 3, 6)

merged_3q_42_ENS <- merge(atlas_cluster_data_w_WI_scenarios_causal, config_3q_42_ENS, by = "DHSID")
merged_3q_42_ENS_final <- merged_3q_42_ENS[,c(1,20:22,26:30)]  #need to understand why this changes some times. 
colnames(merged_3q_42_ENS_final)[8] <- "year"

sf_3q_42_ENS_dhs <- st_as_sf(merged_3q_42_ENS_final, coords = c("long","lat"), crs = ("+init=epsg:4326"))#make sure this crs approach works, works when i plot

#using 2011 as the unique year.  
merged_3q_42_ENS_final_2011 <- subset(merged_3q_42_ENS_final, year == 2011)#1,798 units #weird that it is one diff from before, ie 1799

sf_3q_42_ENS_2011 <- st_as_sf(merged_3q_42_ENS_final_2011, coords = c("long","lat"), crs = ("+init=epsg:4326"))

#do this one first
village_intersects_2013_line_3 <- st_intersects(sf_3q_42_ENS_2011, ug_grid_2013_shape_3km, sparse = FALSE) 
true_index_village_13_3 <- which(apply(village_intersects_2013_line_3,1,any))
true_index_village_13_3_df <- as.data.frame(true_index_village_13_3)
length(true_index_village_13_3_df[,1])
true_matrix_village_13_3 <- matrix(NA,858,2)
true_matrix_village_13_3[,1] <- true_index_village_13_3
#then 2010, then next line

#i find the 2010 intersection
village_intersects_2010_line_3 <- st_intersects(sf_3q_42_ENS_2011, ug_grid_2010_shape_3km, sparse = FALSE) 
true_index_village_10_3 <- which(apply(village_intersects_2010_line_3,1,any))
true_index_village_10_3_df <- as.data.frame(true_index_village_10_3)
length(true_index_village_10_3_df[,1]) #
true_matrix_village_10_3 <- matrix(NA,858,3) #same as 2013
true_matrix_village_10_3[1:758,1] <- true_index_village_10_3
true_matrix_village_10_3[759:853,2] <- 0  

true_matrix_village_13_3[,2] <- true_matrix_village_10_3[,1]
true_matrix_village_13_3 <- as.data.frame(true_matrix_village_13_3)

village_intersects_2016_line_3 <- st_intersects(sf_3q_42_ENS_2011, ug_grid_2016_shape_3km, sparse = FALSE) #this is the line that gives you treated
true_index_village_16_3 <- which(apply(village_intersects_2016_line_3,1,any))
true_index_village_16_3_df <- as.data.frame(true_index_village_16_3)
length(true_index_village_16_3_df[,1]) #
true_matrix_village_16_3 <- matrix(NA,1053,2)
true_matrix_village_16_3[1:1053,1] <- true_index_village_16_3

# Here i am getting a list of non-electrified communities in 2016 (control)
# before applying the rural masks
#sapply(merged_v8_37_ENS_final_2011, class)

control_for_village_WI_3 <- matrix(NA,1798,4)
colnames(control_for_village_WI_3) <- c("row_num", "name", "lat", "long")
control_for_village_WI_3[,1] <- rep(1:1798) 
control_for_village_WI_3[,2] <- as.character(merged_3q_42_ENS_final_2011[,1])
control_for_village_WI_3[,3] <- merged_3q_42_ENS_final_2011[,2]
control_for_village_WI_3[,4] <- merged_3q_42_ENS_final_2011[,3]

# Connected communities in 2016, including those that are treated
connected_village_WI_3 <- control_for_village_WI_3[control_for_village_WI_3[,1] %in% true_matrix_village_16_3[,1], ]
connected_village_WI_3 <- as.data.frame(connected_village_WI_3)

# These are the never connected villages. 
control_for_village_WI_3_list <- control_for_village_WI_3[!control_for_village_WI_3[,1] %in% true_matrix_village_16_3[,1], ]
control_for_village_WI_3_list <- as.data.frame(control_for_village_WI_3_list)

# Now dropping duplicates, which should leave just the newly electrified
treatment_villages_3 <- matrix(NA,858+758,1)
treatment_villages_3[1:858,] <- true_matrix_village_13_3[1:858,1]
treatment_villages_3[859:1616,] <- true_matrix_village_13_3[1:758,2]
treatment_villages_3_list <- as.data.frame(as.numeric(names(which(table(treatment_villages_3)==1))))
colnames(treatment_villages_3_list) <- c("treated_units")
#100 units

treatment_for_village_WI_3 <- matrix(NA,1798,4)
colnames(treatment_for_village_WI_3) <- c("row_num", "name", "lat", "long")
treatment_for_village_WI_3[,1] <- rep(1:1798) 
treatment_for_village_WI_3[,2] <- as.character(merged_3q_42_ENS_final_2011[,1])
treatment_for_village_WI_3[,3] <- merged_3q_42_ENS_final_2011[,2]
treatment_for_village_WI_3[,4] <- merged_3q_42_ENS_final_2011[,3]

#treatment_for_village_WI[,2] <- village_names[,1] add this later
treatment_for_village_WI_3 <- treatment_for_village_WI_3[treatment_for_village_WI_3[,1] %in% treatment_villages_3_list[,1], ]
treatment_for_village_WI_3 <- as.data.frame(treatment_for_village_WI_3)

treat_control_3q_42_ENS_3km <- rbind(control_for_village_WI_3_list, treatment_for_village_WI_3)
#checking for any duplicates - have 745+100 = 845
treat_control_3q_42_ENS_3km <- treat_control_3q_42_ENS_3km[!(duplicated(treat_control_3q_42_ENS_3km) | duplicated(treat_control_3q_42_ENS_3km, fromLast = TRUE)), ]
#now we are at 825 total, dropping 10 from the treat and control so 735 and 90

treat_control_3q_42_ENS_3km$treated <- 0
treat_control_3q_42_ENS_3km[736:825,5] <- 1
colnames(treat_control_3q_42_ENS_3km)[2] <- "DHSID"

c_3q_42_ENS_treat_control_3km <- merge(merged_3q_42_ENS_final, treat_control_3q_42_ENS_3km, by = "DHSID")

c_3q_42_ENS_treat_control_3km <- c_3q_42_ENS_treat_control_3km[,c(1:10,13)]
colnames(c_3q_42_ENS_treat_control_3km)[2] <- "lat"
colnames(c_3q_42_ENS_treat_control_3km)[3] <- "long"
colnames(c_3q_42_ENS_treat_control_3km)[11] <- "treated_unit"

c_3q_42_ENS_treat_control_3km <- c_3q_42_ENS_treat_control_3km[!(c_3q_42_ENS_treat_control_3km$year=="2003"),]
c_3q_42_ENS_treat_control_3km <- c_3q_42_ENS_treat_control_3km[!(c_3q_42_ENS_treat_control_3km$year=="2004"),]
c_3q_42_ENS_treat_control_3km <- c_3q_42_ENS_treat_control_3km[!(c_3q_42_ENS_treat_control_3km$year=="2005"),]
c_3q_42_ENS_treat_control_3km <- c_3q_42_ENS_treat_control_3km[!(c_3q_42_ENS_treat_control_3km$year=="2017"),]
c_3q_42_ENS_treat_control_3km <- c_3q_42_ENS_treat_control_3km[!(c_3q_42_ENS_treat_control_3km$year=="2018"),]
c_3q_42_ENS_treat_control_3km <- c_3q_42_ENS_treat_control_3km[!(c_3q_42_ENS_treat_control_3km$year=="2019"),]

treated <- function(x) { 
  if(x == 2003 | x == 2004 | x == 2005 | x == 2006 | x == 2007 | x == 2008 | x == 2009 | x == 2010) y <- 0
  if(x == 2011 | x == 2012 | x == 2013 | x == 2014 | x == 2015 | x == 2016) y <- 1
  return(y)
}

c_3q_42_ENS_treat_control_3km$treat_year <- sapply(c_3q_42_ENS_treat_control_3km$year,treated)

c_3q_42_ENS_treat_control_3km$treated_id <- c_3q_42_ENS_treat_control_3km$treated_unit*c_3q_42_ENS_treat_control_3km$treat_year

colnames(c_3q_42_ENS_treat_control_3km)[7] <- "cnn_3q_42_ENS_pred"

#11 year base DiD
base_WI_did_3km <- feols(cnn_3q_42_ENS_pred ~ treated_id | DHSID + year, data=c_3q_42_ENS_treat_control_3km)
summary(base_WI_did_3km)

#.053, p = .13  -- aug 31, getting .09, ad .02

c_3q_42_ENS_treat_control_C_3km <- subset(c_3q_42_ENS_treat_control_3km, treated_unit == 0)
c_3q_42_ENS_treat_control_T_3km <- subset(c_3q_42_ENS_treat_control_3km, treated_unit == 1)


# Here are the five splits for our preferred CNN run
# You can download these from the dropbox
# Let me know if you need me to send them
dhs_1_9_84_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_108_dhs.csv")
colnames(dhs_1_9_84_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_216_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_240_dhs.csv")
colnames(dhs_1_9_216_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_348_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_372_dhs.csv")
colnames(dhs_1_9_348_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_480_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_504_dhs.csv")
colnames(dhs_1_9_480_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_612_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_636_dhs.csv")
colnames(dhs_1_9_612_DHSpanel)[1] <- "DHSID_year"

# Here is just some code that needs to be updated with above code. 
dhs_1_9_84_DHSpanel_col1 <- dhs_1_9_84_DHSpanel
colnames(dhs_1_9_84_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_216_DHSpanel_col1 <- dhs_1_9_216_DHSpanel
colnames(dhs_1_9_216_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_348_DHSpanel_col1 <- dhs_1_9_348_DHSpanel
colnames(dhs_1_9_348_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_480_DHSpanel_col1 <- dhs_1_9_480_DHSpanel
colnames(dhs_1_9_480_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_612_DHSpanel_col1 <- dhs_1_9_612_DHSpanel
colnames(dhs_1_9_612_DHSpanel_col1)[1] <- "dhsID_year"

# Here I add the splits to the above DF
c_3q_42_ENS_randomized_treat_control_C_3km <- c_3q_42_ENS_treat_control_C_3km
c_3q_42_ENS_randomized_treat_control_C_3km <- merge(c_3q_42_ENS_randomized_treat_control_C_3km, dhs_1_9_84_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C_3km <- merge(c_3q_42_ENS_randomized_treat_control_C_3km, dhs_1_9_216_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C_3km <- merge(c_3q_42_ENS_randomized_treat_control_C_3km, dhs_1_9_348_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C_3km <- merge(c_3q_42_ENS_randomized_treat_control_C_3km, dhs_1_9_480_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C_3km <- merge(c_3q_42_ENS_randomized_treat_control_C_3km, dhs_1_9_612_DHSpanel_col1[,1:2], by = "dhsID_year")

#c_3q_42_ENS_randomized_treat_control_C_3km <- c_3q_42_ENS_randomized_treat_control_C_3km[,-2]
c_3q_42_ENS_randomized_treat_control_C_3km$meanWI <- rowMeans(c_3q_42_ENS_randomized_treat_control_C_3km[,14:18])

c_3q_42_ENS_randomized_treat_control_T_3km <- c_3q_42_ENS_treat_control_T_3km
c_3q_42_ENS_randomized_treat_control_T_3km <- merge(c_3q_42_ENS_randomized_treat_control_T_3km, dhs_1_9_84_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T_3km <- merge(c_3q_42_ENS_randomized_treat_control_T_3km, dhs_1_9_216_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T_3km <- merge(c_3q_42_ENS_randomized_treat_control_T_3km, dhs_1_9_348_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T_3km <- merge(c_3q_42_ENS_randomized_treat_control_T_3km, dhs_1_9_480_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T_3km <- merge(c_3q_42_ENS_randomized_treat_control_T_3km, dhs_1_9_612_DHSpanel_col1[,1:2], by = "dhsID_year")

#c_3q_42_ENS_randomized_treat_control_T_3km <- c_3q_42_ENS_randomized_treat_control_T_3km[,-2]
c_3q_42_ENS_randomized_treat_control_T_3km$meanWI <- rowMeans(c_3q_42_ENS_randomized_treat_control_T_3km[,14:18])


# Here is the massive for loop. 
# Warning: this started as a basic for loop, and has been added to many times, so it's not optimized for order.

ptm <- proc.time()

set.seed(789) #456 is the second one I've been using

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
p = 90 # number of treated units for the SC-ENt
r = 90 # number of treated units for the SC-EN
full_loop_output_5q_CL1_3km <- matrix(NA, 100, 20) # Column headers are at the bottom.
inside_boot_run <- matrix(NA,1,p) # this is the inside df to capture the bootwise SC-ENt estimates. 
inside_boot_run_2 <- matrix(NA,1,r) # this is the inside df to capture the bootwise SC-EN estimates. 
#coeff_save <- matrix(NA,899,(n*r)+1) # this is a relic from old coeff saving, don't worry about now
#coeff_save[,1] <- c(1:899)
inside_boot_run_pretrend <- matrix(NA,1,p) # this is the inside df to capture the bootwise SC-ENt estimates for the pretrend analysis

for(m in 1:n) {
  
  sample_C <- c_3q_42_ENS_randomized_treat_control_C_3km # import the control df
  randomized_C <- t(apply(sample_C[,14:18], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[20] <- "randomized"
  
  #same as above for the treatment
  sample_T <- c_3q_42_ENS_randomized_treat_control_T_3km 
  randomized_T <- t(apply(sample_T[,14:18], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[20] <- "randomized"
  
  unique_C_units <- unique(sample_C[,2]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 735, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,2]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 90, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here i am running the Did and saving the coeff and p-value to the loop df
  sample_did <- summary(feols(randomized ~ treated_id | DHSID + year, data=sample_TC)) 
  sample_coeff <- as.data.frame(coeftable(sample_did)[1])
  full_loop_output_5q_CL1_3km[m,1] <- sample_coeff[1,1]
  sample_p <- as.data.frame(pvalue(sample_did))
  full_loop_output_5q_CL1_3km[m,2] <- sample_p[1,1]
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated_unit * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  full_loop_output_5q_CL1_3km[m,10] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  full_loop_output_5q_CL1_3km[m,11] <- placebo_p[1,1]
  
  # here i am making the long control df wide, so that i can do MC and SC
  c_df_to_reshape_boot <- sample_control[,c(1,8,20)]
  c_df_to_reshape_boot$id <- rep(1:735, each = 11)
  c_df_to_reshape_boot$re_id <- paste(c_df_to_reshape_boot[,1], c_df_to_reshape_boot[,4], sep="_")
  c_df_to_reshape_boot <- c_df_to_reshape_boot[,c(2,3,5)]
  wide_c_df_boot <- reshape(c_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_c_df_boot  <- wide_c_df_boot  %>% arrange(year)
  wide_c_df_boot  <- wide_c_df_boot [,-1]
  
  # here i am making the long treatment df wide, so that i can do MC and SC
  t_df_to_reshape_boot <- sample_treatment[,c(1,8,20)]
  t_df_to_reshape_boot$id <- rep(1:90, each = 11)
  t_df_to_reshape_boot$re_id <- paste(t_df_to_reshape_boot[,1], t_df_to_reshape_boot[,4], sep="_")
  t_df_to_reshape_boot <- t_df_to_reshape_boot[,c(2,3,5)]
  wide_t_df_boot <- reshape(t_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_t_df_boot <- wide_t_df_boot %>% arrange(year)
  wide_t_df_boot <- wide_t_df_boot[,-1]
  
  # here i cbind them together
  wide_TC_df_boot <- cbind(wide_c_df_boot,wide_t_df_boot)
  
  # here i save the sample average treated values, so that i can compare it later with the MC and SC estimates
  full_loop_output_5q_CL1_3km[m,3] <- rowMeans(wide_TC_df_boot[11,736:825])
  
  # here i am doing the MC code, note MC imputes the whole matrix it predicts values in literally every cell. 
  sample_matrix <- wide_TC_df_boot # input wide df created above 
  
  N_basic_control_test <- 825 # = N
  T_basic_control_test <- 11 # = T
  M_basic_control_test <- sample_matrix # input wide df
  M_basic_control_test <- as.matrix(M_basic_control_test) # has to be a matrix for this code
  
  mask_basic_control_test <- matrix(1,11,825) # this is the zero / one mask
  mask_basic_control_test[6:11,736:825] <- 0 # we put in zeros for every unit we want to predict, i.e every treated unit. 
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1) # 1 = time and unit fixed effects, 0s would not include fixed effects
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  # this is the predicted matrix
  
  # here i am putting the average imputed values in the year we are looking at into the output df
  # this is the counterfactual
  # if i wanted to look at 2015, i would change 11 to 10. 
  full_loop_output_5q_CL1_3km[m,4] <- mean(model_with_both_basic_control_test$est[11,736:825])
  
  # here i subtract the predicted mean from the observed mean to tell how far off this loops prediction was. 
  # this is ATE
  full_loop_output_5q_CL1_3km[m,5] <- full_loop_output_5q_CL1_3km[m,3] - full_loop_output_5q_CL1_3km[m,4]
  
  # now i'm doing the SC-ENt part
  for(o in 1:p) {  
    
    X_prod <- t(as.matrix(sample_matrix[1:5,1:735])) # because we are doing the transpose version here, i select all control values in the pre-years
    Y_prod <- t(as.matrix(sample_matrix[11,1:735])) # because we are doing the transpose version here, i select the post year we are interested in, 11 = 2016
    
    en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100) # this is the panel regression stage
    coeff_prod <- coef(en_prod, s = (en_prod$lambda.min)) # this save the year coeffs in each loop
    
    post_prod <- as.matrix(t(sample_matrix[1:5,735+o])) # these are the units to apply to the coeffs to, note this is where the loop comes in. 
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = (en_prod$lambda.min)) # pred_prod is each units counter factual in 2016
    
    inside_boot_run[1,o] <- pred_prod # saving each counterfactual to the inside boot run 
    
  }
  
  full_loop_output_5q_CL1_3km[m,6] <- mean(inside_boot_run[1,])  # here i store the average counteractual for each boot run
  full_loop_output_5q_CL1_3km[m,7] <- full_loop_output_5q_CL1_3km[m,3]  - full_loop_output_5q_CL1_3km[m,6] # here i find the ATE for each boot run, as with MC 
  
  # here i am creating a df to store the SC-EN (unit) coeffs. 
  # i'm not using it now, but left it here
  inside_coeff <- matrix(NA,735,r+1)
  inside_coeff[,1] <- c(1:735)
  inside_coeff <- as.data.frame(inside_coeff)
  colnames(inside_coeff)[1] <- "row_num"
  
  
  # NOW I am starting the pre-trend test for mc and sc-ent
  
  # as above, i am putting in the average treated unit value for 2010
  full_loop_output_5q_CL1_3km[m,12] <- rowMeans(wide_TC_df_boot[5,736:825])
  
  # i then follow the same MC procedure as above except that it stops at 2010. 
  pre_trend_sample_matrix <- sample_matrix[1:5,]
  
  N_basic_control_test <- 825
  T_basic_control_test <- 5
  M_basic_control_test <- pre_trend_sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,5,825)
  mask_basic_control_test[5,736:825] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  model_with_both_basic_control_test$err <- model_with_both_basic_control_test$est - M_basic_control_test
  
  model_with_both_basic_control_test$msk_err <- model_with_both_basic_control_test$err*(1-mask_basic_control_test)
  
  model_with_both_basic_control_test$test_RMSE <- 
    sqrt((1/sum(1-mask_basic_control_test)) * sum(model_with_both_basic_control_test$msk_err^2))
  
  full_loop_output_5q_CL1_3km[m,13] <- mean(model_with_both_basic_control_test$est[5,736:825]) # this is the average 2010 counterfactual
  full_loop_output_5q_CL1_3km[m,14] <- full_loop_output_5q_CL1_3km[m,12] - full_loop_output_5q_CL1_3km[m,13] # this is the pretrend estimated effect
  
  # here i take the imputed values 2010 treated from the proceeding step and slice them into the observed values df 
  # i then recreate a long DF and run a faux Did using the imputed units
  # the goal is to see if the MC imputed units have a better pretrend test than the DiD pretrend test above. 
  new_pretrend_df_mc <- pre_trend_sample_matrix
  new_pretrend_df_mc[5,736:825] <- model_with_both_basic_control_test$est[5,736:825]
  new_pretrend_df_mc_long <- reshape(new_pretrend_df_mc, 
                                     direction = "long", list(names(new_pretrend_df_mc)[1:825]))
  colnames(new_pretrend_df_mc_long)[1] <- "unit_id"
  colnames(new_pretrend_df_mc_long)[2] <- "randomized"
  colnames(new_pretrend_df_mc_long)[3] <- "year_id"
  new_pretrend_df_mc_long$year <- c(2006, 2007, 2008, 2009, 2010)
  new_pretrend_df_mc_long$treated_unit <- 0
  new_pretrend_df_mc_long[3676:4125,5] <- 1
  new_pretrend_df_mc_long$treated_year <- c(0,0,0,0,1)
  new_pretrend_df_mc_long$treated_id <- new_pretrend_df_mc_long$treated_unit * new_pretrend_df_mc_long$treated_year
  
  # i'm storing the MC / DiD pre-trend output here. 
  MC_did <- summary(feols(randomized ~ treated_id | unit_id + year, data=new_pretrend_df_mc_long))
  sample_coeff_mc <- as.data.frame(coeftable(MC_did)[1])
  full_loop_output_5q_CL1_3km[m,15] <- sample_coeff_mc[1,1]
  MC_p <- as.data.frame(pvalue(MC_did))
  full_loop_output_5q_CL1_3km[m,16] <- MC_p[1,1]
  
  # here i follow the same steps as above to predict treated values in 2010 via SC-ENt
  for(o in 1:p) {  
    
    X_prod <- t(as.matrix(pre_trend_sample_matrix[1:4,1:735]))
    Y_prod <- t(as.matrix(pre_trend_sample_matrix[5,1:735]))
    
    en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
    coeff_prod <- coef(en_prod, s = (en_prod$lambda.min))
    
    post_prod <- as.matrix(t(sample_matrix[1:4,735+o]))
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = (en_prod$lambda.min))
    
    inside_boot_run_pretrend[1,o] <- pred_prod
    
  }
  
  full_loop_output_5q_CL1_3km[m,17] <- mean(inside_boot_run_pretrend[1,])  # this is the average 2010 counterfactual
  full_loop_output_5q_CL1_3km[m,18] <- full_loop_output_5q_CL1_3km[m,12]  - full_loop_output_5q_CL1_3km[m,17] # this is the pretrend estimated effect
  
  # as with mc above, i am splicing in the imputed treated units from the pretrend SC analysis
  # and then running a faux DiD
  new_pretrend_df_SC <- pre_trend_sample_matrix
  new_pretrend_df_SC[5,736:825] <- inside_boot_run_pretrend[1,1:90]
  new_pretrend_df_SC_long <- reshape(new_pretrend_df_SC, 
                                     direction = "long", list(names(new_pretrend_df_SC)[1:825]))
  colnames(new_pretrend_df_SC_long)[1] <- "unit_id"
  colnames(new_pretrend_df_SC_long)[2] <- "randomized"
  colnames(new_pretrend_df_SC_long)[3] <- "year_id"
  new_pretrend_df_SC_long$year <- c(2006, 2007, 2008, 2009, 2010)
  new_pretrend_df_SC_long$treated_unit <- 0
  new_pretrend_df_SC_long[3676:4125,5] <- 1
  new_pretrend_df_SC_long$treated_year <- c(0,0,0,0,1)
  new_pretrend_df_SC_long$treated_id <- new_pretrend_df_SC_long$treated_unit * new_pretrend_df_SC_long$treated_year
  
  # i'm storing the SC / DiD pre-trend output here. 
  SC_did <- summary(feols(randomized ~ treated_id | unit_id + year, data=new_pretrend_df_SC_long))
  sample_coeff_SC <- as.data.frame(coeftable(SC_did)[1])
  full_loop_output_5q_CL1_3km[m,19] <- sample_coeff_SC[1,1]
  SC_p <- as.data.frame(pvalue(SC_did))
  full_loop_output_5q_CL1_3km[m,20] <- SC_p[1,1]
  
}


# i've include an explanation for the colnames below
colnames(full_loop_output_5q_CL1_3km) <- c("full_did_coeff", "full_did_Pvalue", "rando_2016_ave", "mc_2016_ave", "mc_2016_ATE", 
                                           "scT_2016_ave", "scT_2016_ATE", "sc_2016_Ave", "sc_2016_ATE", "placebo_did_coeff", 
                                           "placebo_did_Pvalue", "rando_2010_ave", "mc_2010_ave", "mc_2010_ATE", "mc_did_coeff", 
                                           "mc_did_Pvalue", "sc_2010_ave", "sc_2010_ATE", "sc_did_coeff", "sc_did_Pvalue")

proc.time() - ptm   

#write.csv(full_loop_output_5q_CL1_3km, "full_loop_with_pretrends_and_DiDs_3km_5penalty_cvglment_aug31.csv")

# Column name explanations
#1 "full_did_coeff", this is the 11 year, base DiD estimate
#2 "full_did_Pvalue", this is the corresponding p value
#3 "rando_2016_ave", this is when you realize how bad your column names are :), and is also the randomized mean treatment value in 2016
#4 "mc_2016_ave", this is the mean MC counterfactual for 2016
#5 "mc_2016_ATE", this is the ATE as predicted by MC
#6 "scT_2016_ave", this is the mean SCt counterfactual for 2016
#7 "scT_2016_ATE", this is the ATE as predicted by SCt
#8 "sc_2016_Ave", not used any longer
#9 "sc_2016_ATE", not used any longer
#10 "placebo_did_coeff", this is the pretrend DiD test coeff for each loop
#11 "placebo_did_Pvalue", this is the corresponding p value
#12 "rando_2010_ave", this is the randomized mean treatment value in 2010
#13 "mc_2010_ave", this is the mean MC counterfactual for 2010
#14 "mc_2010_ATE", this is the ATE as predicted by MC for 2010 and used for pretrend comparison
#15 "mc_did_coeff", this is the 2010 pretend DiD test when i splice in the imputed MC estimates
#16 "mc_did_Pvalue", this is the corresponding p value
#17 "sc_2010_ave", this is the mean SC counterfactual for 2010
#18 "sc_2010_ATE", this is the ATE as predicted by SC for 2010 and used for pretrend comparison
#19 "sc_did_coeff", this is the 2010 pretend DiD test when i splice in the imputed SC estimates
#20 "sc_did_Pvalue", this is the corresponding p value

# Here are the analysis steps
## means, sds, etc etc, 

# Here is a basic output chart.  Not all the analysis is formatted / mechanized yet, but this gives the basic output 
Full_Boot_Analysis_3km <- matrix(NA,20,5)
colnames(Full_Boot_Analysis_3km) <- c("run", "mean", "sd", "- 95%", "+ 95%" )
Full_Boot_Analysis_3km[1,1] <- "full_did_coeff"
Full_Boot_Analysis_3km[1,2] <-  round(mean(full_loop_output_5q_CL1_3km[,1]), 3)
Full_Boot_Analysis_3km[1,3] <-  round(sd(full_loop_output_5q_CL1_3km[,1]), 3)
rank_did <- as.data.frame(full_loop_output_5q_CL1_3km[,1])
rank_did <- as.data.frame(rank_did[order(rank_did[,1],decreasing = FALSE) , ])
Full_Boot_Analysis_3km[1,4] <- round(mean(rank_did[2:3,]), 3)
Full_Boot_Analysis_3km[1,5] <- round(mean(rank_did[98:99,]), 3)

Full_Boot_Analysis_3km[2,1] <- "mc_2016_ATE"
Full_Boot_Analysis_3km[2,2] <-  round(mean(full_loop_output_5q_CL1_3km[,5]), 3)
Full_Boot_Analysis_3km[2,3] <-  round(sd(full_loop_output_5q_CL1_3km[,5]), 3)
rank_mc <- as.data.frame(full_loop_output_5q_CL1_3km[,5])
rank_mc <- as.data.frame(rank_mc[order(rank_mc[,1],decreasing = FALSE) , ])
Full_Boot_Analysis_3km[2,4] <- round(mean(rank_mc[2:3,]), 3)
Full_Boot_Analysis_3km[2,5] <- round(mean(rank_mc[98:99,])  , 3)

Full_Boot_Analysis_3km[3,1] <- "scT_2016_ATE"
Full_Boot_Analysis_3km[3,2] <-  round(mean(full_loop_output_5q_CL1_3km[,7]), 3)
Full_Boot_Analysis_3km[3,3] <-  round(sd(full_loop_output_5q_CL1_3km[,7]), 3)
rank_sc <- as.data.frame(full_loop_output_5q_CL1_3km[,7])
rank_sc <- as.data.frame(rank_sc[order(rank_sc[,1],decreasing = FALSE) , ])
Full_Boot_Analysis_3km[3,4] <- round(mean(rank_sc[2:3,]), 3)
Full_Boot_Analysis_3km[3,5] <- round(mean(rank_sc[98:99,]) , 3)

Full_Boot_Analysis_3km[4,1] <- "Pretrend Analysis Below"
Full_Boot_Analysis_3km[4,4] <- "Pvalue < .1"

Full_Boot_Analysis_3km[5,1] <- "placebo_did_coeff"
Full_Boot_Analysis_3km[5,2] <-  round(mean(full_loop_output_5q_CL1_3km[,10]), 3)
Full_Boot_Analysis_3km[5,3] <-  round(sd(full_loop_output_5q_CL1_3km[,10]), 3)
Full_Boot_Analysis_3km[5,4] <-  length(which(full_loop_output_5q_CL1_3km[,11] < .1))
Full_Boot_Analysis_3km[5,5] <-  NA    

Full_Boot_Analysis_3km[6,1] <- "mc_did_coeff"
Full_Boot_Analysis_3km[6,2] <-  round(mean(full_loop_output_5q_CL1_3km[,15]), 3)
Full_Boot_Analysis_3km[6,3] <-  round(sd(full_loop_output_5q_CL1_3km[,15]), 3)
Full_Boot_Analysis_3km[6,4] <-  length(which(full_loop_output_5q_CL1_3km[,16] < .1))
Full_Boot_Analysis_3km[6,5] <-  NA

Full_Boot_Analysis_3km[7,1] <- "sc_did_coeff"
Full_Boot_Analysis_3km[7,2] <-  round(mean(full_loop_output_5q_CL1_3km[,19]), 3)
Full_Boot_Analysis_3km[7,3] <-  round(sd(full_loop_output_5q_CL1_3km[,19]), 3)
Full_Boot_Analysis_3km[7,4] <-  length(which(full_loop_output_5q_CL1_3km[,20] < .1))
Full_Boot_Analysis_3km[7,5] <-  NA

# here is a plot depicting dd, mc and sc-en precision


full_loop_output_5q_CL1_3km <- as.data.frame(full_loop_output_5q_CL1_3km)

pre_prediction_comparison <- ggplot(data = full_loop_output_5q_CL1_3km) +
  geom_point(aes(x = 1, y = placebo_did_coeff)) +
  geom_point(aes(x = 2, y = mc_2010_ATE), color = "springgreen") +
  geom_point(aes(x = 3, y = sc_2010_ATE), color = "dodgerblue") +
  theme_clean()
pre_prediction_comparison


# annaul means
control_year_means <- c_3q_42_ENS_randomized_treat_control_C_3km %>%
  dplyr::group_by(year) %>% 
  dplyr::summarise(year_mean = mean(meanWI, na.rm = TRUE))

treat_year_means <- c_3q_42_ENS_randomized_treat_control_T_3km %>%
  dplyr::group_by(year) %>% 
  dplyr::summarise(year_mean = mean(meanWI, na.rm = TRUE))


###########
# here is the 4km analysis
##########

###########
# Here is the 4km single / non-Ensemble version # 
###########
##### CI_V0_1_7_config_id_42 #####

#do this one first
village_intersects_2013_line_4 <- st_intersects(sf_3q_42_ENS_2011, ug_grid_2013_shape_4km, sparse = FALSE) 
true_index_village_13_4 <- which(apply(village_intersects_2013_line_4,1,any))
true_index_village_13_4_df <- as.data.frame(true_index_village_13_4)
length(true_index_village_13_4_df[,1])
true_matrix_village_13_4 <- matrix(NA,955,2)
true_matrix_village_13_4[,1] <- true_index_village_13_4
#then 2010, then next line

#i find the 2010 intersection
village_intersects_2010_line_4 <- st_intersects(sf_3q_42_ENS_2011, ug_grid_2010_shape_4km, sparse = FALSE) 
true_index_village_10_4 <- which(apply(village_intersects_2010_line_4,1,any))
true_index_village_10_4_df <- as.data.frame(true_index_village_10_4)
length(true_index_village_10_4_df[,1]) #
true_matrix_village_10_4 <- matrix(NA,955,3) #same as 2013
true_matrix_village_10_4[1:850,1] <- true_index_village_10_4
true_matrix_village_10_4[851:955,2] <- 0  

true_matrix_village_13_4[,2] <- true_matrix_village_10_4[,1]
true_matrix_village_13_4 <- as.data.frame(true_matrix_village_13_4)

village_intersects_2016_line_4 <- st_intersects(sf_3q_42_ENS_2011, ug_grid_2016_shape_4km, sparse = FALSE) #this is the line that gives you treated
true_index_village_16_4 <- which(apply(village_intersects_2016_line_4,1,any))
true_index_village_16_4_df <- as.data.frame(true_index_village_16_4)
length(true_index_village_16_4_df[,1]) #
true_matrix_village_16_4 <- matrix(NA,1173,2)
true_matrix_village_16_4[1:1173,1] <- true_index_village_16_4

# Here i am getting a list of non-electrified communities in 2016 (control)
# before applying the rural masks
#sapply(merged_v8_37_ENS_final_2011, class)

control_for_village_WI_4 <- matrix(NA,1798,4)
colnames(control_for_village_WI_4) <- c("row_num", "name", "lat", "long")
control_for_village_WI_4[,1] <- rep(1:1798) 
control_for_village_WI_4[,2] <- as.character(merged_3q_42_ENS_final_2011[,1])
control_for_village_WI_4[,3] <- merged_3q_42_ENS_final_2011[,2]
control_for_village_WI_4[,4] <- merged_3q_42_ENS_final_2011[,3]

# Connected communities in 2016, including those that are treated
connected_village_WI_4 <- control_for_village_WI_4[control_for_village_WI_4[,1] %in% true_matrix_village_16_4[,1], ]
connected_village_WI_4 <- as.data.frame(connected_village_WI_4)

# These are the never connected villages. 
control_for_village_WI_4_list <- control_for_village_WI_4[!control_for_village_WI_4[,1] %in% true_matrix_village_16_4[,1], ]
control_for_village_WI_4_list <- as.data.frame(control_for_village_WI_4_list)

# Now dropping duplicates, which should leave just the newly electrified
treatment_villages_4 <- matrix(NA,955+850,1)
treatment_villages_4[1:955,] <- true_matrix_village_13_4[1:955,1]
treatment_villages_4[956:1805,] <- true_matrix_village_13_4[1:850,2]
treatment_villages_4_list <- as.data.frame(as.numeric(names(which(table(treatment_villages_4)==1))))
colnames(treatment_villages_4_list) <- c("treated_units")
#100 units

treatment_for_village_WI_4 <- matrix(NA,1798,4)
colnames(treatment_for_village_WI_4) <- c("row_num", "name", "lat", "long")
treatment_for_village_WI_4[,1] <- rep(1:1798) 
treatment_for_village_WI_4[,2] <- as.character(merged_3q_42_ENS_final_2011[,1])
treatment_for_village_WI_4[,3] <- merged_3q_42_ENS_final_2011[,2]
treatment_for_village_WI_4[,4] <- merged_3q_42_ENS_final_2011[,3]

#treatment_for_village_WI[,2] <- village_names[,1] add this later
treatment_for_village_WI_4 <- treatment_for_village_WI_4[treatment_for_village_WI_4[,1] %in% treatment_villages_4_list[,1], ]
treatment_for_village_WI_4 <- as.data.frame(treatment_for_village_WI_4)

treat_control_3q_42_ENS_4km <- rbind(control_for_village_WI_4_list, treatment_for_village_WI_4)
#checking for any duplicates - have 625+105 = 730
treat_control_3q_42_ENS_4km <- treat_control_3q_42_ENS_4km[!(duplicated(treat_control_3q_42_ENS_4km) | duplicated(treat_control_3q_42_ENS_4km, fromLast = TRUE)), ]
#now we are at 706 total, dropping 12 from the treat and control so 613 and 93

treat_control_3q_42_ENS_4km$treated <- 0
treat_control_3q_42_ENS_4km[614:706,5] <- 1
colnames(treat_control_3q_42_ENS_4km)[2] <- "DHSID"

c_3q_42_ENS_treat_control_4km <- merge(merged_3q_42_ENS_final, treat_control_3q_42_ENS_4km, by = "DHSID")

c_3q_42_ENS_treat_control_4km <- c_3q_42_ENS_treat_control_4km[,c(1:10,13)]
colnames(c_3q_42_ENS_treat_control_4km)[2] <- "lat"
colnames(c_3q_42_ENS_treat_control_4km)[3] <- "long"
colnames(c_3q_42_ENS_treat_control_4km)[11] <- "treated_unit"

c_3q_42_ENS_treat_control_4km <- c_3q_42_ENS_treat_control_4km[!(c_3q_42_ENS_treat_control_4km$year=="2003"),]
c_3q_42_ENS_treat_control_4km <- c_3q_42_ENS_treat_control_4km[!(c_3q_42_ENS_treat_control_4km$year=="2004"),]
c_3q_42_ENS_treat_control_4km <- c_3q_42_ENS_treat_control_4km[!(c_3q_42_ENS_treat_control_4km$year=="2005"),]
c_3q_42_ENS_treat_control_4km <- c_3q_42_ENS_treat_control_4km[!(c_3q_42_ENS_treat_control_4km$year=="2017"),]
c_3q_42_ENS_treat_control_4km <- c_3q_42_ENS_treat_control_4km[!(c_3q_42_ENS_treat_control_4km$year=="2018"),]
c_3q_42_ENS_treat_control_4km <- c_3q_42_ENS_treat_control_4km[!(c_3q_42_ENS_treat_control_4km$year=="2019"),]

treated <- function(x) { 
  if(x == 2003 | x == 2004 | x == 2005 | x == 2006 | x == 2007 | x == 2008 | x == 2009 | x == 2010) y <- 0
  if(x == 2011 | x == 2012 | x == 2013 | x == 2014 | x == 2015 | x == 2016) y <- 1
  return(y)
}

c_3q_42_ENS_treat_control_4km$treat_year <- sapply(c_3q_42_ENS_treat_control_4km$year,treated)

c_3q_42_ENS_treat_control_4km$treated_id <- c_3q_42_ENS_treat_control_4km$treated_unit*c_3q_42_ENS_treat_control_4km$treat_year

colnames(c_3q_42_ENS_treat_control_4km)[7] <- "cnn_3q_42_ENS_pred"

#11 year base DiD
base_WI_did_4km <- feols(cnn_3q_42_ENS_pred ~ treated_id | DHSID + year, data=c_3q_42_ENS_treat_control_4km)
summary(base_WI_did_4km)

#.046, p = .17

c_3q_42_ENS_treat_control_C_4km <- subset(c_3q_42_ENS_treat_control_4km, treated_unit == 0)
c_3q_42_ENS_treat_control_T_4km <- subset(c_3q_42_ENS_treat_control_4km, treated_unit == 1)


# Here are the five splits for our preferred CNN run
# You can download these from the dropbox
# Let me know if you need me to send them
dhs_1_9_84_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_108_dhs.csv")
colnames(dhs_1_9_84_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_216_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_240_dhs.csv")
colnames(dhs_1_9_216_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_348_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_372_dhs.csv")
colnames(dhs_1_9_348_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_480_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_504_dhs.csv")
colnames(dhs_1_9_480_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_612_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_636_dhs.csv")
colnames(dhs_1_9_612_DHSpanel)[1] <- "DHSID_year"

# Here is just some code that needs to be updated with above code. 
dhs_1_9_84_DHSpanel_col1 <- dhs_1_9_84_DHSpanel
colnames(dhs_1_9_84_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_216_DHSpanel_col1 <- dhs_1_9_216_DHSpanel
colnames(dhs_1_9_216_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_348_DHSpanel_col1 <- dhs_1_9_348_DHSpanel
colnames(dhs_1_9_348_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_480_DHSpanel_col1 <- dhs_1_9_480_DHSpanel
colnames(dhs_1_9_480_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_612_DHSpanel_col1 <- dhs_1_9_612_DHSpanel
colnames(dhs_1_9_612_DHSpanel_col1)[1] <- "dhsID_year"

# Here I add the splits to the above DF
c_3q_42_ENS_randomized_treat_control_C_4km <- c_3q_42_ENS_treat_control_C_4km
c_3q_42_ENS_randomized_treat_control_C_4km <- merge(c_3q_42_ENS_randomized_treat_control_C_4km, dhs_1_9_84_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C_4km <- merge(c_3q_42_ENS_randomized_treat_control_C_4km, dhs_1_9_216_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C_4km <- merge(c_3q_42_ENS_randomized_treat_control_C_4km, dhs_1_9_348_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C_4km <- merge(c_3q_42_ENS_randomized_treat_control_C_4km, dhs_1_9_480_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C_4km <- merge(c_3q_42_ENS_randomized_treat_control_C_4km, dhs_1_9_612_DHSpanel_col1[,1:2], by = "dhsID_year")

#c_3q_42_ENS_randomized_treat_control_C_4km <- c_3q_42_ENS_randomized_treat_control_C_4km[,-2]
c_3q_42_ENS_randomized_treat_control_C_4km$meanWI <- rowMeans(c_3q_42_ENS_randomized_treat_control_C_4km[,14:18])

c_3q_42_ENS_randomized_treat_control_T_4km <- c_3q_42_ENS_treat_control_T_4km
c_3q_42_ENS_randomized_treat_control_T_4km <- merge(c_3q_42_ENS_randomized_treat_control_T_4km, dhs_1_9_84_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T_4km <- merge(c_3q_42_ENS_randomized_treat_control_T_4km, dhs_1_9_216_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T_4km <- merge(c_3q_42_ENS_randomized_treat_control_T_4km, dhs_1_9_348_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T_4km <- merge(c_3q_42_ENS_randomized_treat_control_T_4km, dhs_1_9_480_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T_4km <- merge(c_3q_42_ENS_randomized_treat_control_T_4km, dhs_1_9_612_DHSpanel_col1[,1:2], by = "dhsID_year")

#c_3q_42_ENS_randomized_treat_control_T_4km <- c_3q_42_ENS_randomized_treat_control_T_4km[,-2]
c_3q_42_ENS_randomized_treat_control_T_4km$meanWI <- rowMeans(c_3q_42_ENS_randomized_treat_control_T_4km[,14:18])


# Here is the massive for loop. 
# Warning: this started as a basic for loop, and has been added to many times, so it's not optimized for order.

ptm <- proc.time()

set.seed(789) #456 is the second one I've been using

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
p = 93 # number of treated units for the SC-ENt
r = 93 # number of treated units for the SC-EN
full_loop_output_5q_CL1_4km <- matrix(NA, 100, 20) # Column headers are at the bottom.
inside_boot_run <- matrix(NA,1,p) # this is the inside df to capture the bootwise SC-ENt estimates. 
inside_boot_run_2 <- matrix(NA,1,r) # this is the inside df to capture the bootwise SC-EN estimates. 
#coeff_save <- matrix(NA,899,(n*r)+1) # this is a relic from old coeff saving, don't worry about now
#coeff_save[,1] <- c(1:899)
inside_boot_run_pretrend <- matrix(NA,1,p) # this is the inside df to capture the bootwise SC-ENt estimates for the pretrend analysis

for(m in 1:n) {
  
  sample_C <- c_3q_42_ENS_randomized_treat_control_C_4km # import the control df
  randomized_C <- t(apply(sample_C[,14:18], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[20] <- "randomized"
  
  #same as above for the treatment
  sample_T <- c_3q_42_ENS_randomized_treat_control_T_4km 
  randomized_T <- t(apply(sample_T[,14:18], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[20] <- "randomized"
  
  unique_C_units <- unique(sample_C[,2]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 613, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,2]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 93, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here i am running the Did and saving the coeff and p-value to the loop df
  sample_did <- summary(feols(randomized ~ treated_id | DHSID + year, data=sample_TC)) 
  sample_coeff <- as.data.frame(coeftable(sample_did)[1])
  full_loop_output_5q_CL1_4km[m,1] <- sample_coeff[1,1]
  sample_p <- as.data.frame(pvalue(sample_did))
  full_loop_output_5q_CL1_4km[m,2] <- sample_p[1,1]
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated_unit * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  full_loop_output_5q_CL1_4km[m,10] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  full_loop_output_5q_CL1_4km[m,11] <- placebo_p[1,1]
  
  # here i am making the long control df wide, so that i can do MC and SC
  c_df_to_reshape_boot <- sample_control[,c(1,8,20)]
  c_df_to_reshape_boot$id <- rep(1:613, each = 11)
  c_df_to_reshape_boot$re_id <- paste(c_df_to_reshape_boot[,1], c_df_to_reshape_boot[,4], sep="_")
  c_df_to_reshape_boot <- c_df_to_reshape_boot[,c(2,3,5)]
  wide_c_df_boot <- reshape(c_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_c_df_boot  <- wide_c_df_boot  %>% arrange(year)
  wide_c_df_boot  <- wide_c_df_boot [,-1]
  
  # here i am making the long treatment df wide, so that i can do MC and SC
  t_df_to_reshape_boot <- sample_treatment[,c(1,8,20)]
  t_df_to_reshape_boot$id <- rep(1:93, each = 11)
  t_df_to_reshape_boot$re_id <- paste(t_df_to_reshape_boot[,1], t_df_to_reshape_boot[,4], sep="_")
  t_df_to_reshape_boot <- t_df_to_reshape_boot[,c(2,3,5)]
  wide_t_df_boot <- reshape(t_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_t_df_boot <- wide_t_df_boot %>% arrange(year)
  wide_t_df_boot <- wide_t_df_boot[,-1]
  
  # here i cbind them together
  wide_TC_df_boot <- cbind(wide_c_df_boot,wide_t_df_boot)
  
  # here i save the sample average treated values, so that i can compare it later with the MC and SC estimates
  full_loop_output_5q_CL1_4km[m,3] <- rowMeans(wide_TC_df_boot[11,614:706])
  
  # here i am doing the MC code, note MC imputes the whole matrix it predicts values in literally every cell. 
  sample_matrix <- wide_TC_df_boot # input wide df created above 
  
  N_basic_control_test <- 706 # = N
  T_basic_control_test <- 11 # = T
  M_basic_control_test <- sample_matrix # input wide df
  M_basic_control_test <- as.matrix(M_basic_control_test) # has to be a matrix for this code
  
  mask_basic_control_test <- matrix(1,11,706) # this is the zero / one mask
  mask_basic_control_test[6:11,614:706] <- 0 # we put in zeros for every unit we want to predict, i.e every treated unit. 
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1) # 1 = time and unit fixed effects, 0s would not include fixed effects
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  # this is the predicted matrix
  
  # here i am putting the average imputed values in the year we are looking at into the output df
  # this is the counterfactual
  # if i wanted to look at 2015, i would change 11 to 10. 
  full_loop_output_5q_CL1_4km[m,4] <- mean(model_with_both_basic_control_test$est[11,614:706])
  
  # here i subtract the predicted mean from the observed mean to tell how far off this loops prediction was. 
  # this is ATE
  full_loop_output_5q_CL1_4km[m,5] <- full_loop_output_5q_CL1_4km[m,3] - full_loop_output_5q_CL1_4km[m,4]
  
  # now i'm doing the SC-ENt part
  for(o in 1:p) {  
    
    X_prod <- t(as.matrix(sample_matrix[1:5,1:613])) # because we are doing the transpose version here, i select all control values in the pre-years
    Y_prod <- t(as.matrix(sample_matrix[11,1:613])) # because we are doing the transpose version here, i select the post year we are interested in, 11 = 2016
    
    en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100) # this is the panel regression stage
    coeff_prod <- coef(en_prod, s = (en_prod$lambda.min)) # this save the year coeffs in each loop
    
    post_prod <- as.matrix(t(sample_matrix[1:5,613+o])) # these are the units to apply to the coeffs to, note this is where the loop comes in. 
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = (en_prod$lambda.min)) # pred_prod is each units counter factual in 2016
    
    inside_boot_run[1,o] <- pred_prod # saving each counterfactual to the inside boot run 
    
  }
  
  full_loop_output_5q_CL1_4km[m,6] <- mean(inside_boot_run[1,])  # here i store the average counteractual for each boot run
  full_loop_output_5q_CL1_4km[m,7] <- full_loop_output_5q_CL1_4km[m,3]  - full_loop_output_5q_CL1_4km[m,6] # here i find the ATE for each boot run, as with MC 
  
  # here i am creating a df to store the SC-EN (unit) coeffs. 
  # i'm not using it now, but left it here
  inside_coeff <- matrix(NA,613,r+1)
  inside_coeff[,1] <- c(1:613)
  inside_coeff <- as.data.frame(inside_coeff)
  colnames(inside_coeff)[1] <- "row_num"
  
  
  # NOW I am starting the pre-trend test for mc and sc-ent
  
  # as above, i am putting in the average treated unit value for 2010
  full_loop_output_5q_CL1_4km[m,12] <- rowMeans(wide_TC_df_boot[5,614:706])
  
  # i then follow the same MC procedure as above except that it stops at 2010. 
  pre_trend_sample_matrix <- sample_matrix[1:5,]
  
  N_basic_control_test <- 706
  T_basic_control_test <- 5
  M_basic_control_test <- pre_trend_sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,5,706)
  mask_basic_control_test[5,614:706] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  model_with_both_basic_control_test$err <- model_with_both_basic_control_test$est - M_basic_control_test
  
  model_with_both_basic_control_test$msk_err <- model_with_both_basic_control_test$err*(1-mask_basic_control_test)
  
  model_with_both_basic_control_test$test_RMSE <- 
    sqrt((1/sum(1-mask_basic_control_test)) * sum(model_with_both_basic_control_test$msk_err^2))
  
  full_loop_output_5q_CL1_4km[m,13] <- mean(model_with_both_basic_control_test$est[5,614:706]) # this is the average 2010 counterfactual
  full_loop_output_5q_CL1_4km[m,14] <- full_loop_output_5q_CL1_4km[m,12] - full_loop_output_5q_CL1_4km[m,13] # this is the pretrend estimated effect
  
  # here i take the imputed values 2010 treated from the proceeding step and slice them into the observed values df 
  # i then recreate a long DF and run a faux Did using the imputed units
  # the goal is to see if the MC imputed units have a better pretrend test than the DiD pretrend test above. 
  new_pretrend_df_mc <- pre_trend_sample_matrix
  new_pretrend_df_mc[5,614:706] <- model_with_both_basic_control_test$est[5,614:706]
  new_pretrend_df_mc_long <- reshape(new_pretrend_df_mc, 
                                     direction = "long", list(names(new_pretrend_df_mc)[1:706]))
  colnames(new_pretrend_df_mc_long)[1] <- "unit_id"
  colnames(new_pretrend_df_mc_long)[2] <- "randomized"
  colnames(new_pretrend_df_mc_long)[3] <- "year_id"
  new_pretrend_df_mc_long$year <- c(2006, 2007, 2008, 2009, 2010)
  new_pretrend_df_mc_long$treated_unit <- 0
  new_pretrend_df_mc_long[3066:3530,5] <- 1
  new_pretrend_df_mc_long$treated_year <- c(0,0,0,0,1)
  new_pretrend_df_mc_long$treated_id <- new_pretrend_df_mc_long$treated_unit * new_pretrend_df_mc_long$treated_year
  
  # i'm storing the MC / DiD pre-trend output here. 
  MC_did <- summary(feols(randomized ~ treated_id | unit_id + year, data=new_pretrend_df_mc_long))
  sample_coeff_mc <- as.data.frame(coeftable(MC_did)[1])
  full_loop_output_5q_CL1_4km[m,15] <- sample_coeff_mc[1,1]
  MC_p <- as.data.frame(pvalue(MC_did))
  full_loop_output_5q_CL1_4km[m,16] <- MC_p[1,1]
  
  # here i follow the same steps as above to predict treated values in 2010 via SC-ENt
  for(o in 1:p) {  
    
    X_prod <- t(as.matrix(pre_trend_sample_matrix[1:4,1:613]))
    Y_prod <- t(as.matrix(pre_trend_sample_matrix[5,1:613]))
    
    en_prod <- cv.glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
    coeff_prod <- coef(en_prod, s = (en_prod$lambda.min))
    
    post_prod <- as.matrix(t(sample_matrix[1:4,613+o]))
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = (en_prod$lambda.min))
    
    inside_boot_run_pretrend[1,o] <- pred_prod
    
  }
  
  full_loop_output_5q_CL1_4km[m,17] <- mean(inside_boot_run_pretrend[1,])  # this is the average 2010 counterfactual
  full_loop_output_5q_CL1_4km[m,18] <- full_loop_output_5q_CL1_4km[m,12]  - full_loop_output_5q_CL1_4km[m,17] # this is the pretrend estimated effect
  
  # as with mc above, i am splicing in the imputed treated units from the pretrend SC analysis
  # and then running a faux DiD
  new_pretrend_df_SC <- pre_trend_sample_matrix
  new_pretrend_df_SC[5,614:706] <- inside_boot_run_pretrend[1,1:93]
  new_pretrend_df_SC_long <- reshape(new_pretrend_df_SC, 
                                     direction = "long", list(names(new_pretrend_df_SC)[1:706]))
  colnames(new_pretrend_df_SC_long)[1] <- "unit_id"
  colnames(new_pretrend_df_SC_long)[2] <- "randomized"
  colnames(new_pretrend_df_SC_long)[3] <- "year_id"
  new_pretrend_df_SC_long$year <- c(2006, 2007, 2008, 2009, 2010)
  new_pretrend_df_SC_long$treated_unit <- 0
  new_pretrend_df_SC_long[3066:3530,5] <- 1
  new_pretrend_df_SC_long$treated_year <- c(0,0,0,0,1)
  new_pretrend_df_SC_long$treated_id <- new_pretrend_df_SC_long$treated_unit * new_pretrend_df_SC_long$treated_year
  
  # i'm storing the SC / DiD pre-trend output here. 
  SC_did <- summary(feols(randomized ~ treated_id | unit_id + year, data=new_pretrend_df_SC_long))
  sample_coeff_SC <- as.data.frame(coeftable(SC_did)[1])
  full_loop_output_5q_CL1_4km[m,19] <- sample_coeff_SC[1,1]
  SC_p <- as.data.frame(pvalue(SC_did))
  full_loop_output_5q_CL1_4km[m,20] <- SC_p[1,1]
  
}


# i've include an explanation for the colnames below
colnames(full_loop_output_5q_CL1_4km) <- c("full_did_coeff", "full_did_Pvalue", "rando_2016_ave", "mc_2016_ave", "mc_2016_ATE", 
                                           "scT_2016_ave", "scT_2016_ATE", "sc_2016_Ave", "sc_2016_ATE", "placebo_did_coeff", 
                                           "placebo_did_Pvalue", "rando_2010_ave", "mc_2010_ave", "mc_2010_ATE", "mc_did_coeff", 
                                           "mc_did_Pvalue", "sc_2010_ave", "sc_2010_ATE", "sc_did_coeff", "sc_did_Pvalue")

proc.time() - ptm   

#write.csv(full_loop_output_5q_CL1_4km, "full_loop_with_pretrends_and_DiDs_4km_5penalty_cvglment_aug31.csv")

# Column name explanations
#1 "full_did_coeff", this is the 11 year, base DiD estimate
#2 "full_did_Pvalue", this is the corresponding p value
#3 "rando_2016_ave", this is when you realize how bad your column names are :), and is also the randomized mean treatment value in 2016
#4 "mc_2016_ave", this is the mean MC counterfactual for 2016
#5 "mc_2016_ATE", this is the ATE as predicted by MC
#6 "scT_2016_ave", this is the mean SCt counterfactual for 2016
#7 "scT_2016_ATE", this is the ATE as predicted by SCt
#8 "sc_2016_Ave", not used any longer
#9 "sc_2016_ATE", not used any longer
#10 "placebo_did_coeff", this is the pretrend DiD test coeff for each loop
#11 "placebo_did_Pvalue", this is the corresponding p value
#12 "rando_2010_ave", this is the randomized mean treatment value in 2010
#13 "mc_2010_ave", this is the mean MC counterfactual for 2010
#14 "mc_2010_ATE", this is the ATE as predicted by MC for 2010 and used for pretrend comparison
#15 "mc_did_coeff", this is the 2010 pretend DiD test when i splice in the imputed MC estimates
#16 "mc_did_Pvalue", this is the corresponding p value
#17 "sc_2010_ave", this is the mean SC counterfactual for 2010
#18 "sc_2010_ATE", this is the ATE as predicted by SC for 2010 and used for pretrend comparison
#19 "sc_did_coeff", this is the 2010 pretend DiD test when i splice in the imputed SC estimates
#20 "sc_did_Pvalue", this is the corresponding p value

# Here are the analysis steps
## means, sds, etc etc, 

# Here is a basic output chart.  Not all the analysis is formatted / mechanized yet, but this gives the basic output 
Full_Boot_Analysis_4km <- matrix(NA,20,5)
colnames(Full_Boot_Analysis_4km) <- c("run", "mean", "sd", "- 95%", "+ 95%" )
Full_Boot_Analysis_4km[1,1] <- "full_did_coeff"
Full_Boot_Analysis_4km[1,2] <-  round(mean(full_loop_output_5q_CL1_4km[,1]), 3)
Full_Boot_Analysis_4km[1,3] <-  round(sd(full_loop_output_5q_CL1_4km[,1]), 3)
rank_did <- as.data.frame(full_loop_output_5q_CL1_4km[,1])
rank_did <- as.data.frame(rank_did[order(rank_did[,1],decreasing = FALSE) , ])
Full_Boot_Analysis_4km[1,4] <- round(mean(rank_did[2:3,]), 3)
Full_Boot_Analysis_4km[1,5] <- round(mean(rank_did[98:99,]), 3)

Full_Boot_Analysis_4km[2,1] <- "mc_2016_ATE"
Full_Boot_Analysis_4km[2,2] <-  round(mean(full_loop_output_5q_CL1_4km[,5]), 3)
Full_Boot_Analysis_4km[2,3] <-  round(sd(full_loop_output_5q_CL1_4km[,5]), 3)
rank_mc <- as.data.frame(full_loop_output_5q_CL1_4km[,5])
rank_mc <- as.data.frame(rank_mc[order(rank_mc[,1],decreasing = FALSE) , ])
Full_Boot_Analysis_4km[2,4] <- round(mean(rank_mc[2:3,]), 3)
Full_Boot_Analysis_4km[2,5] <- round(mean(rank_mc[98:99,])  , 3)

Full_Boot_Analysis_4km[3,1] <- "scT_2016_ATE"
Full_Boot_Analysis_4km[3,2] <-  round(mean(full_loop_output_5q_CL1_4km[,7]), 3)
Full_Boot_Analysis_4km[3,3] <-  round(sd(full_loop_output_5q_CL1_4km[,7]), 3)
rank_sc <- as.data.frame(full_loop_output_5q_CL1_4km[,7])
rank_sc <- as.data.frame(rank_sc[order(rank_sc[,1],decreasing = FALSE) , ])
Full_Boot_Analysis_4km[3,4] <- round(mean(rank_sc[2:3,]), 3)
Full_Boot_Analysis_4km[3,5] <- round(mean(rank_sc[98:99,]) , 3)

Full_Boot_Analysis_4km[4,1] <- "Pretrend Analysis Below"
Full_Boot_Analysis_4km[4,4] <- "Pvalue < .1"

Full_Boot_Analysis_4km[5,1] <- "placebo_did_coeff"
Full_Boot_Analysis_4km[5,2] <-  round(mean(full_loop_output_5q_CL1_4km[,10]), 3)
Full_Boot_Analysis_4km[5,3] <-  round(sd(full_loop_output_5q_CL1_4km[,10]), 3)
Full_Boot_Analysis_4km[5,4] <-  length(which(full_loop_output_5q_CL1_4km[,11] < .1))
Full_Boot_Analysis_4km[5,5] <-  NA    

Full_Boot_Analysis_4km[6,1] <- "mc_did_coeff"
Full_Boot_Analysis_4km[6,2] <-  round(mean(full_loop_output_5q_CL1_4km[,15]), 3)
Full_Boot_Analysis_4km[6,3] <-  round(sd(full_loop_output_5q_CL1_4km[,15]), 3)
Full_Boot_Analysis_4km[6,4] <-  length(which(full_loop_output_5q_CL1_4km[,16] < .1))
Full_Boot_Analysis_4km[6,5] <-  NA

Full_Boot_Analysis_4km[7,1] <- "sc_did_coeff"
Full_Boot_Analysis_4km[7,2] <-  round(mean(full_loop_output_5q_CL1_4km[,19]), 3)
Full_Boot_Analysis_4km[7,3] <-  round(sd(full_loop_output_5q_CL1_4km[,19]), 3)
Full_Boot_Analysis_4km[7,4] <-  length(which(full_loop_output_5q_CL1_4km[,20] < .1))
Full_Boot_Analysis_4km[7,5] <-  NA

# here is a plot depicting dd, mc and sc-en precision


full_loop_output_5q_CL1_4km <- as.data.frame(full_loop_output_5q_CL1_4km)

pre_prediction_comparison <- ggplot(data = full_loop_output_5q_CL1_4km) +
  geom_point(aes(x = 1, y = placebo_did_coeff)) +
  geom_point(aes(x = 2, y = mc_2010_ATE), color = "springgreen") +
  geom_point(aes(x = 3, y = sc_2010_ATE), color = "dodgerblue") +
  theme_clean()
pre_prediction_comparison


# annaul means
control_year_means <- c_3q_42_ENS_randomized_treat_control_C_4km %>%
  dplyr::group_by(year) %>% 
  dplyr::summarise(year_mean = mean(meanWI, na.rm = TRUE))

treat_year_means <- c_3q_42_ENS_randomized_treat_control_T_4km %>%
  dplyr::group_by(year) %>% 
  dplyr::summarise(year_mean = mean(meanWI, na.rm = TRUE))


#supp fig 3km and 4km
pdf("supp_results_3km_4km_5penalty_august2021.pdf", width=6, height=6) 
supp_results_3km_4km_5penalty_august2021 <- ggplot() +
  geom_point(aes(x = .1, y = .222), color = "springgreen", size = 2) +
  geom_segment(aes(x = .1, y = 0.082, xend = .1, yend = 0.371), color = "springgreen", size = .2) +
  geom_point(aes(x = .15, y = 0.279), color = "dodgerblue", size = 2) +
  geom_segment(aes(x = .15, y = 0.136, xend = .15, yend = 0.421), color = "dodgerblue", size = .2) +
  
  geom_point(aes(x = .25, y = 0.25), color = "springgreen", size = 2) +
  geom_segment(aes(x = .25, y = .091, xend = .25, yend = 0.405), color = "springgreen", size = .2) +
  geom_point(aes(x = .30, y = 0.302), color = "dodgerblue", size = 2) +
  geom_segment(aes(x = .30, y = 0.157, xend = .30, yend = 0.421), color = "dodgerblue", size = .2) +
  
  geom_point(aes(x = .4, y = 0.2), color = "springgreen", size = 2) +
  geom_segment(aes(x = .4, y = .072, xend = .4, yend = .36), color = "springgreen", size = .2) +
  geom_point(aes(x = .45, y = 0.241), color = "dodgerblue", size = 2) +
  geom_segment(aes(x = .45, y = 0.105, xend = .45, yend = 0.421), color = "dodgerblue", size = .2) +
  
  ylim(-.25,.75) +
  xlim(.05,.5) +
  ylab("Estimated Causal Effect") +
  xlab("") +
  theme_clean() +
  theme(legend.position = "none", plot.background = element_rect(color = "white"), axis.title = element_text(size = 16), axis.text = element_text(size = 16), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
supp_results_3km_4km_5penalty_august2021
dev.off()

save.image(file = "july_15_main.Rdata")
save.image(file = "july_31a_main.Rdata")


load("july_15_main.Rdata")

load("july_31a_main.Rdata")






































######
these are the other splits





# Here are the five splits for our preferred CNN run
# You can download these from the dropbox
# Let me know if you need me to send them
dhs_1_9_96_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_96_dhs.csv")
colnames(dhs_1_9_96_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_228_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_228_dhs.csv")
colnames(dhs_1_9_228_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_360_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_360_dhs.csv")
colnames(dhs_1_9_360_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_492_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_492_dhs.csv")
colnames(dhs_1_9_492_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_624_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_624_dhs.csv")
colnames(dhs_1_9_624_DHSpanel)[1] <- "DHSID_year"

# Here is just some code that needs to be updated with above code. 
dhs_1_9_96_DHSpanel_col1 <- dhs_1_9_96_DHSpanel
colnames(dhs_1_9_96_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_228_DHSpanel_col1 <- dhs_1_9_228_DHSpanel
colnames(dhs_1_9_228_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_360_DHSpanel_col1 <- dhs_1_9_360_DHSpanel
colnames(dhs_1_9_360_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_492_DHSpanel_col1 <- dhs_1_9_492_DHSpanel
colnames(dhs_1_9_492_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_624_DHSpanel_col1 <- dhs_1_9_624_DHSpanel
colnames(dhs_1_9_624_DHSpanel_col1)[1] <- "dhsID_year"

#loading the 0 CL
dhs_1_9_0_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_0_dhs.csv")
colnames(dhs_1_9_0_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_228_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_228_dhs.csv")
colnames(dhs_1_9_228_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_360_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_360_dhs.csv")
colnames(dhs_1_9_360_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_492_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_492_dhs.csv")
colnames(dhs_1_9_492_DHSpanel)[1] <- "DHSID_year"
dhs_1_9_624_DHSpanel <- read.csv("preds_CI-V0-1-9_config_id_624_dhs.csv")
colnames(dhs_1_9_624_DHSpanel)[1] <- "DHSID_year"

# Here is just some code that needs to be updated with above code. 
dhs_1_9_96_DHSpanel_col1 <- dhs_1_9_0_DHSpanel
colnames(dhs_1_9_96_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_228_DHSpanel_col1 <- dhs_1_9_132_DHSpanel
colnames(dhs_1_9_228_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_360_DHSpanel_col1 <- dhs_1_9_264_DHSpanel
colnames(dhs_1_9_360_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_492_DHSpanel_col1 <- dhs_1_9_396_DHSpanel
colnames(dhs_1_9_492_DHSpanel_col1)[1] <- "dhsID_year"
dhs_1_9_624_DHSpanel_col1 <- dhs_1_9_528_DHSpanel
colnames(dhs_1_9_624_DHSpanel_col1)[1] <- "dhsID_year"


# Here I add the splits to the above DF
# w 3 cl
c_3q_42_ENS_randomized_treat_control_C <- c_3q_42_ENS_treat_control_C
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_96_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_228_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_360_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_492_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_624_DHSpanel_col1[,1:2], by = "dhsID_year")

c_3q_42_ENS_randomized_treat_control_C$meanWI <- rowMeans(c_3q_42_ENS_randomized_treat_control_C[,14:18])
c_3q_42_ENS_randomized_treat_control_C <- c_3q_42_ENS_randomized_treat_control_C[,-2]

c_3q_42_ENS_randomized_treat_control_T <- c_3q_42_ENS_treat_control_T
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_96_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_228_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_360_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_492_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_624_DHSpanel_col1[,1:2], by = "dhsID_year")

c_3q_42_ENS_randomized_treat_control_T$meanWI <- rowMeans(c_3q_42_ENS_randomized_treat_control_T[,14:18])
c_3q_42_ENS_randomized_treat_control_T <- c_3q_42_ENS_randomized_treat_control_T[,-2]


# with 0 cl
c_3q_42_ENS_randomized_treat_control_C <- c_3q_42_ENS_treat_control_C
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_0_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_132_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_264_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_396_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_C <- merge(c_3q_42_ENS_randomized_treat_control_C, dhs_1_9_528_DHSpanel_col1[,1:2], by = "dhsID_year")

c_3q_42_ENS_randomized_treat_control_C$meanWI <- rowMeans(c_3q_42_ENS_randomized_treat_control_C[,14:18])
c_3q_42_ENS_randomized_treat_control_C <- c_3q_42_ENS_randomized_treat_control_C[,-2]

c_3q_42_ENS_randomized_treat_control_T <- c_3q_42_ENS_treat_control_T
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_0_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_132_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_264_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_396_DHSpanel_col1[,1:2], by = "dhsID_year")
c_3q_42_ENS_randomized_treat_control_T <- merge(c_3q_42_ENS_randomized_treat_control_T, dhs_1_9_528_DHSpanel_col1[,1:2], by = "dhsID_year")

c_3q_42_ENS_randomized_treat_control_T$meanWI <- rowMeans(c_3q_42_ENS_randomized_treat_control_T[,14:18])
c_3q_42_ENS_randomized_treat_control_T <- c_3q_42_ENS_randomized_treat_control_T[,-2]

# Here is the massive for loop. 
# Warning: this started as a basic for loop, and has been added to many times, so it's not optimized for order.

ptm <- proc.time()

set.seed(789) #456 is the second one I've been using

n = 100 # Number of total loops. I do 100 because 500 or 1,000 take too long on my computer.  Feel free to up it!
p = 76 # number of treated units for the SC-ENt
r = 76 # number of treated units for the SC-EN
full_loop_output_5q_CL1 <- matrix(NA, 100, 20) # Column headers are at the bottom.
inside_boot_run <- matrix(NA,1,p) # this is the inside df to capture the bootwise SC-ENt estimates. 
inside_boot_run_2 <- matrix(NA,1,r) # this is the inside df to capture the bootwise SC-EN estimates. 
coeff_save <- matrix(NA,899,(n*r)+1) # this is a relic from old coeff saving, don't worry about now
coeff_save[,1] <- c(1:899)
inside_boot_run_pretrend <- matrix(NA,1,p) # this is the inside df to capture the bootwise SC-ENt estimates for the pretrend analysis

for(m in 1:n) {
  
  sample_C <- c_3q_42_ENS_randomized_treat_control_C # import the control df
  randomized_C <- t(apply(sample_C[,14:18], 1, function(d) sample(d, 1))) # randomly select which of the five splits estimates to use per unit. 
  # is the t() here eating tons of memory again?
  randomized_C <- as.data.frame(randomized_C)
  sample_C<- cbind(sample_C, t(randomized_C)) # adding the radomized column to the df
  colnames(sample_C)[20] <- "randomized"
  
  #same as above for the treatment
  sample_T <- c_3q_42_ENS_randomized_treat_control_T 
  randomized_T <- t(apply(sample_T[,14:18], 1, function(d) sample(d, 1)))
  randomized_T <- as.data.frame(randomized_T)
  sample_T<- cbind(sample_T, t(randomized_T))
  colnames(sample_T)[20] <- "randomized"
  
  unique_C_units <- unique(sample_C[,2]) #here i am just simplifying the df for randomization
  unique_C_units <- as.data.frame(unique_C_units)
  sample_unique_C_units <- unique_C_units[sample(nrow(unique_C_units), size = 888, replace = TRUE),] #here i randomized units for the current loop. 
  sample_unique_C_units <- as.data.frame(sample_unique_C_units)
  colnames(sample_unique_C_units) <- c("DHSID") 
  sample_control <- merge(sample_unique_C_units, sample_C, by = "DHSID") # here i merge the randomized selection with the full df
  
  #same as above for the treatment
  unique_T_units <- unique(sample_T[,2]) 
  unique_T_units <- as.data.frame(unique_T_units)
  sample_unique_T_units <- unique_T_units[sample(nrow(unique_T_units), size = 76, replace = TRUE),]
  sample_unique_T_units <- as.data.frame(sample_unique_T_units)
  colnames(sample_unique_T_units) <- c("DHSID") 
  sample_treatment <- merge(sample_unique_T_units, sample_T, by = "DHSID")
  
  sample_TC <- rbind(sample_control, sample_treatment) # here i bind the randomized control and treatment together for DiD analysis
  
  # here i am running the Did and saving the coeff and p-value to the loop df
  sample_did <- summary(feols(randomized ~ treated_id | DHSID + year, data=sample_TC)) 
  sample_coeff <- as.data.frame(coeftable(sample_did)[1])
  full_loop_output_5q_CL1[m,1] <- sample_coeff[1,1]
  sample_p <- as.data.frame(pvalue(sample_did))
  full_loop_output_5q_CL1[m,2] <- sample_p[1,1]
  
  # here is the DiD pretrend test
  placebo_test <- subset(sample_TC, year < 2011)
  
  fake_treated <- function(x) { 
    if(x == 2006 | x == 2007 | x == 2008 | x == 2009 ) y <- 0
    if(x == 2010 ) y <- 1
    return(y)
  }
  
  placebo_test$fake_treat_year <- sapply(placebo_test$year, fake_treated)
  placebo_test$fake_treat_id <- placebo_test$treated_unit * placebo_test$fake_treat_year
  
  # as above, i am running the pretrend DiD and saving the coeff and p-value for each loop.   
  placebo_did <- summary(feols(randomized ~ fake_treat_id | DHSID + year, data=placebo_test))
  sample_coeff <- as.data.frame(coeftable(placebo_did)[1])
  full_loop_output_5q_CL1[m,10] <-sample_coeff[1,1]
  placebo_p <- as.data.frame(pvalue(placebo_did))
  full_loop_output_5q_CL1[m,11] <- placebo_p[1,1]
  
  # here i am making the long control df wide, so that i can do MC and SC
  c_df_to_reshape_boot <- sample_control[,c(1,8,20)]
  c_df_to_reshape_boot$id <- rep(1:888, each = 11)
  c_df_to_reshape_boot$re_id <- paste(c_df_to_reshape_boot[,1], c_df_to_reshape_boot[,4], sep="_")
  c_df_to_reshape_boot <- c_df_to_reshape_boot[,c(2,3,5)]
  wide_c_df_boot <- reshape(c_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_c_df_boot  <- wide_c_df_boot  %>% arrange(year)
  wide_c_df_boot  <- wide_c_df_boot [,-1]
  
  # here i am making the long treatment df wide, so that i can do MC and SC
  t_df_to_reshape_boot <- sample_treatment[,c(1,8,20)]
  t_df_to_reshape_boot$id <- rep(1:76, each = 11)
  t_df_to_reshape_boot$re_id <- paste(t_df_to_reshape_boot[,1], t_df_to_reshape_boot[,4], sep="_")
  t_df_to_reshape_boot <- t_df_to_reshape_boot[,c(2,3,5)]
  wide_t_df_boot <- reshape(t_df_to_reshape_boot, idvar = "year", timevar = "re_id", direction = "wide")
  wide_t_df_boot <- wide_t_df_boot %>% arrange(year)
  wide_t_df_boot <- wide_t_df_boot[,-1]
  
  # here i cbind them together
  wide_TC_df_boot <- cbind(wide_c_df_boot,wide_t_df_boot)
  
  # here i save the sample average treated values, so that i can compare it later with the MC and SC estimates
  full_loop_output_5q_CL1[m,3] <- rowMeans(wide_TC_df_boot[11,889:964])
  
  # here i am doing the MC code, note MC imputes the whole matrix it predicts values in literally every cell. 
  sample_matrix <- wide_TC_df_boot # input wide df created above 
  
  N_basic_control_test <- 964 # = N
  T_basic_control_test <- 11 # = T
  M_basic_control_test <- sample_matrix # input wide df
  M_basic_control_test <- as.matrix(M_basic_control_test) # has to be a matrix for this code
  
  mask_basic_control_test <- matrix(1,11,964) # this is the zero / one mask
  mask_basic_control_test[6:11,889:964] <- 0 # we put in zeros for every unit we want to predict, i.e every treated unit. 
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1) # 1 = time and unit fixed effects, 0s would not include fixed effects
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  # this is the predicted matrix
  
  # we don't need any of these but i leave them in just in case
  model_with_both_basic_control_test$err <- model_with_both_basic_control_test$est - M_basic_control_test
  model_with_both_basic_control_test$msk_err <- model_with_both_basic_control_test$err*(1-mask_basic_control_test)
  model_with_both_basic_control_test$test_RMSE <- 
    sqrt((1/sum(1-mask_basic_control_test)) * sum(model_with_both_basic_control_test$msk_err^2))
  
  # here i am putting the average imputed values in the year we are looking at into the output df
  # this is the counterfactual
  # if i wanted to look at 2015, i would change 11 to 10. 
  full_loop_output_5q_CL1[m,4] <- mean(model_with_both_basic_control_test$est[11,889:964])
  
  # here i subtract the predicted mean from the observed mean to tell how far off this loops prediction was. 
  # this is ATE
  full_loop_output_5q_CL1[m,5] <- full_loop_output_5q_CL1[m,3] - full_loop_output_5q_CL1[m,4]
  
  # now i'm doing the SC-ENt part
  for(o in 1:p) {  
    
    X_prod <- t(as.matrix(sample_matrix[1:5,1:888])) # because we are doing the transpose version here, i select all control values in the pre-years
    Y_prod <- t(as.matrix(sample_matrix[11,1:888])) # because we are doing the transpose version here, i select the post year we are interested in, 11 = 2016
    
    en_prod <- glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100) # this is the panel regression stage
    coeff_prod <- coef(en_prod, s = min(en_prod$lambda)) # this save the year coeffs in each loop
    
    post_prod <- as.matrix(t(sample_matrix[1:5,888+o])) # these are the units to apply to the coeffs to, note this is where the loop comes in. 
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = min(en_prod$lambda)) # pred_prod is each units counter factual in 2016
    
    inside_boot_run[1,o] <- pred_prod # saving each counterfactual to the inside boot run 
    
  }
  
  full_loop_output_5q_CL1[m,6] <- mean(inside_boot_run[1,])  # here i store the average counteractual for each boot run
  full_loop_output_5q_CL1[m,7] <- full_loop_output_5q_CL1[m,3]  - full_loop_output_5q_CL1[m,6] # here i find the ATE for each boot run, as with MC 
  
  # here i am creating a df to store the SC-EN (unit) coeffs. 
  # i'm not using it now, but left it here
  inside_coeff <- matrix(NA,888,r+1)
  inside_coeff[,1] <- c(1:888)
  inside_coeff <- as.data.frame(inside_coeff)
  colnames(inside_coeff)[1] <- "row_num"
  
  # here is the normal SC-EN run. i'm not using it in the analysis, but left it hear just to see how it goes 
  # the steps are mirror to above
  for(q in 1:r) {
    X_prod <- (as.matrix(sample_matrix[1:5,1:888]))
    Y_prod <- (as.matrix(sample_matrix[1:5,888+q]))
    
    en_prod <- glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
    coeff_prod <- coef(en_prod, s = min(en_prod$lambda))
    
    post_prod <- as.matrix((sample_matrix[11,1:888]))
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = min(en_prod$lambda))
    
    inside_boot_run_2[1,q] <- pred_prod
    
    coeff_prod_df <- as.data.frame(summary(coeff_prod))
    coeff_prod_df  <- coeff_prod_df[-1,] #to get rid of the intercept
    colnames(coeff_prod_df)[1] <- "row_num"
    coeff_prod_df <- coeff_prod_df[,c(1,3)]
    
    inside_coeff <- merge(coeff_prod_df, inside_coeff, by = "row_num", all = TRUE)
    
  }
  assign( paste("inside_coeff", m, sep = "_") , inside_coeff , envir = globalenv() ) 
  
  full_loop_output_5q_CL1[m,8] <- mean(inside_boot_run_2[1,]) 
  full_loop_output_5q_CL1[m,9] <- full_loop_output_5q_CL1[m,3]  - full_loop_output_5q_CL1[m,8] 
  
  # NOW I am starting the pre-trend test for mc and sc-ent
  
  # as above, i am putting in the average treated unit value for 2010
  full_loop_output_5q_CL1[m,12] <- rowMeans(wide_TC_df_boot[5,889:964])
  
  # i then follow the same MC procedure as above except that it stops at 2010. 
  pre_trend_sample_matrix <- sample_matrix[1:5,]
  
  N_basic_control_test <- 964
  T_basic_control_test <- 5
  M_basic_control_test <- pre_trend_sample_matrix
  M_basic_control_test <- as.matrix(M_basic_control_test)
  
  mask_basic_control_test <- matrix(1,5,964)
  mask_basic_control_test[5,889:964] <- 0
  mask_basic_control_test <- as.matrix(mask_basic_control_test)
  
  model_with_both_basic_control_test <- mcnnm_cv(M_basic_control_test, mask_basic_control_test, to_estimate_u = 1, 
                                                 to_estimate_v = 1)
  
  model_with_both_basic_control_test$est <- model_with_both_basic_control_test$L + 
    replicate(N_basic_control_test,model_with_both_basic_control_test$u) + t(replicate(T_basic_control_test,model_with_both_basic_control_test$v))
  
  model_with_both_basic_control_test$err <- model_with_both_basic_control_test$est - M_basic_control_test
  
  model_with_both_basic_control_test$msk_err <- model_with_both_basic_control_test$err*(1-mask_basic_control_test)
  
  model_with_both_basic_control_test$test_RMSE <- 
    sqrt((1/sum(1-mask_basic_control_test)) * sum(model_with_both_basic_control_test$msk_err^2))
  
  full_loop_output_5q_CL1[m,13] <- mean(model_with_both_basic_control_test$est[5,889:964]) # this is the average 2010 counterfactual
  full_loop_output_5q_CL1[m,14] <- full_loop_output_5q_CL1[m,12] - full_loop_output_5q_CL1[m,13] # this is the pretrend estimated effect
  
  # here i take the imputed values 2010 treated from the proceeding step and slice them into the observed values df 
  # i then recreate a long DF and run a faux Did using the imputed units
  # the goal is to see if the MC imputed units have a better pretrend test than the DiD pretrend test above. 
  new_pretrend_df_mc <- pre_trend_sample_matrix
  new_pretrend_df_mc[5,889:964] <- model_with_both_basic_control_test$est[5,889:964]
  new_pretrend_df_mc_long <- reshape(new_pretrend_df_mc, 
                                     direction = "long", list(names(new_pretrend_df_mc)[1:964]))
  colnames(new_pretrend_df_mc_long)[1] <- "unit_id"
  colnames(new_pretrend_df_mc_long)[2] <- "randomized"
  colnames(new_pretrend_df_mc_long)[3] <- "year_id"
  new_pretrend_df_mc_long$year <- c(2006, 2007, 2008, 2009, 2010)
  new_pretrend_df_mc_long$treated_unit <- 0
  new_pretrend_df_mc_long[4441:4820,5] <- 1
  new_pretrend_df_mc_long$treated_year <- c(0,0,0,0,1)
  new_pretrend_df_mc_long$treated_id <- new_pretrend_df_mc_long$treated_unit * new_pretrend_df_mc_long$treated_year
  
  # i'm storing the MC / DiD pre-trend output here. 
  MC_did <- summary(feols(randomized ~ treated_id | unit_id + year, data=new_pretrend_df_mc_long))
  sample_coeff_mc <- as.data.frame(coeftable(MC_did)[1])
  full_loop_output_5q_CL1[m,15] <- sample_coeff_mc[1,1]
  MC_p <- as.data.frame(pvalue(MC_did))
  full_loop_output_5q_CL1[m,16] <- MC_p[1,1]
  
  # here i follow the same steps as above to predict treated values in 2010 via SC-ENt
  for(o in 1:p) {  
    
    X_prod <- t(as.matrix(pre_trend_sample_matrix[1:4,1:888]))
    Y_prod <- t(as.matrix(pre_trend_sample_matrix[5,1:888]))
    
    en_prod <- glmnet(X_prod, Y_prod, alpha = .5, nlambda = 100)
    coeff_prod <- coef(en_prod, s = min(en_prod$lambda))
    
    post_prod <- as.matrix(t(sample_matrix[1:4,888+o]))
    pred_prod <- predict(en_prod, newx = post_prod, 
                         s = min(en_prod$lambda))
    
    inside_boot_run_pretrend[1,o] <- pred_prod
    
  }
  
  full_loop_output_5q_CL1[m,17] <- mean(inside_boot_run_pretrend[1,])  # this is the average 2010 counterfactual
  full_loop_output_5q_CL1[m,18] <- full_loop_output_5q_CL1[m,12]  - full_loop_output_5q_CL1[m,17] # this is the pretrend estimated effect
  
  # as with mc above, i am splicing in the imputed treated units from the pretrend SC analysis
  # and then running a faux DiD
  new_pretrend_df_SC <- pre_trend_sample_matrix
  new_pretrend_df_SC[5,889:964] <- inside_boot_run_pretrend[1,1:76]
  new_pretrend_df_SC_long <- reshape(new_pretrend_df_SC, 
                                     direction = "long", list(names(new_pretrend_df_SC)[1:964]))
  colnames(new_pretrend_df_SC_long)[1] <- "unit_id"
  colnames(new_pretrend_df_SC_long)[2] <- "randomized"
  colnames(new_pretrend_df_SC_long)[3] <- "year_id"
  new_pretrend_df_SC_long$year <- c(2006, 2007, 2008, 2009, 2010)
  new_pretrend_df_SC_long$treated_unit <- 0
  new_pretrend_df_SC_long[4441:4820,5] <- 1
  new_pretrend_df_SC_long$treated_year <- c(0,0,0,0,1)
  new_pretrend_df_SC_long$treated_id <- new_pretrend_df_SC_long$treated_unit * new_pretrend_df_SC_long$treated_year
  
  # i'm storing the SC / DiD pre-trend output here. 
  SC_did <- summary(feols(randomized ~ treated_id | unit_id + year, data=new_pretrend_df_SC_long))
  sample_coeff_SC <- as.data.frame(coeftable(SC_did)[1])
  full_loop_output_5q_CL1[m,19] <- sample_coeff_SC[1,1]
  SC_p <- as.data.frame(pvalue(SC_did))
  full_loop_output_5q_CL1[m,20] <- SC_p[1,1]
  
}


# i've include an explanation for the colnames below
colnames(full_loop_output_5q_CL1) <- c("full_did_coeff", "full_did_Pvalue", "rando_2016_ave", "mc_2016_ave", "mc_2016_ATE", 
                                       "scT_2016_ave", "scT_2016_ATE", "sc_2016_Ave", "sc_2016_ATE", "placebo_did_coeff", 
                                       "placebo_did_Pvalue", "rando_2010_ave", "mc_2010_ave", "mc_2010_ATE", "mc_did_coeff", 
                                       "mc_did_Pvalue", "sc_2010_ave", "sc_2010_ATE", "sc_did_coeff", "sc_did_Pvalue")

proc.time() - ptm   

#write.csv(full_loop_output_5q_CL1, "full_loop_with_pretrends_and_DiDs_3CL.csv")

# Column name explanations
#1 "full_did_coeff", this is the 11 year, base DiD estimate
#2 "full_did_Pvalue", this is the corresponding p value
#3 "rando_2016_ave", this is when you realize how bad your column names are :), and is also the randomized mean treatment value in 2016
#4 "mc_2016_ave", this is the mean MC counterfactual for 2016
#5 "mc_2016_ATE", this is the ATE as predicted by MC
#6 "scT_2016_ave", this is the mean SCt counterfactual for 2016
#7 "scT_2016_ATE", this is the ATE as predicted by SCt
#8 "sc_2016_Ave", not used any longer
#9 "sc_2016_ATE", not used any longer
#10 "placebo_did_coeff", this is the pretrend DiD test coeff for each loop
#11 "placebo_did_Pvalue", this is the corresponding p value
#12 "rando_2010_ave", this is the randomized mean treatment value in 2010
#13 "mc_2010_ave", this is the mean MC counterfactual for 2010
#14 "mc_2010_ATE", this is the ATE as predicted by MC for 2010 and used for pretrend comparison
#15 "mc_did_coeff", this is the 2010 pretend DiD test when i splice in the imputed MC estimates
#16 "mc_did_Pvalue", this is the corresponding p value
#17 "sc_2010_ave", this is the mean SC counterfactual for 2010
#18 "sc_2010_ATE", this is the ATE as predicted by SC for 2010 and used for pretrend comparison
#19 "sc_did_coeff", this is the 2010 pretend DiD test when i splice in the imputed SC estimates
#20 "sc_did_Pvalue", this is the corresponding p value

# Here are the analysis steps
## means, sds, etc etc, 

# Here is a basic output chart.  Not all the analysis is formatted / mechanized yet, but this gives the basic output 
Full_Boot_Analysis <- matrix(NA,20,5)
colnames(Full_Boot_Analysis) <- c("run", "mean", "sd", "- 95%", "+ 95%" )
Full_Boot_Analysis[1,1] <- "full_did_coeff"
Full_Boot_Analysis[1,2] <-  round(mean(full_loop_output_5q_CL1[,1]), 3)
Full_Boot_Analysis[1,3] <-  round(sd(full_loop_output_5q_CL1[,1]), 3)
rank_did <- as.data.frame(full_loop_output_5q_CL1[,1])
rank_did <- as.data.frame(rank_did[order(rank_did[,1],decreasing = FALSE) , ])
Full_Boot_Analysis[1,4] <- round(mean(rank_did[2:3,]), 3)
Full_Boot_Analysis[1,5] <- round(mean(rank_did[98:99,]), 3)

Full_Boot_Analysis[2,1] <- "mc_2016_ATE"
Full_Boot_Analysis[2,2] <-  round(mean(full_loop_output_5q_CL1[,5]), 3)
Full_Boot_Analysis[2,3] <-  round(sd(full_loop_output_5q_CL1[,5]), 3)
rank_mc <- as.data.frame(full_loop_output_5q_CL1[,5])
rank_mc <- as.data.frame(rank_mc[order(rank_mc[,1],decreasing = FALSE) , ])
Full_Boot_Analysis[2,4] <- round(mean(rank_mc[2:3,]), 3)
Full_Boot_Analysis[2,5] <- round(mean(rank_mc[98:99,])  , 3)

Full_Boot_Analysis[3,1] <- "scT_2016_ATE"
Full_Boot_Analysis[3,2] <-  round(mean(full_loop_output_5q_CL1[,7]), 3)
Full_Boot_Analysis[3,3] <-  round(sd(full_loop_output_5q_CL1[,7]), 3)
rank_sc <- as.data.frame(full_loop_output_5q_CL1[,7])
rank_sc <- as.data.frame(rank_sc[order(rank_sc[,1],decreasing = FALSE) , ])
Full_Boot_Analysis[3,4] <- round(mean(rank_sc[2:3,]), 3)
Full_Boot_Analysis[3,5] <- round(mean(rank_sc[98:99,]) , 3)

Full_Boot_Analysis[4,1] <- "Pretrend Analysis Below"
Full_Boot_Analysis[4,4] <- "Pvalue < .1"

Full_Boot_Analysis[5,1] <- "placebo_did_coeff"
Full_Boot_Analysis[5,2] <-  round(mean(full_loop_output_5q_CL1[,10]), 3)
Full_Boot_Analysis[5,3] <-  round(sd(full_loop_output_5q_CL1[,10]), 3)
Full_Boot_Analysis[5,4] <-  length(which(full_loop_output_5q_CL1[,11] < .1))
Full_Boot_Analysis[5,5] <-  NA    

Full_Boot_Analysis[6,1] <- "mc_did_coeff"
Full_Boot_Analysis[6,2] <-  round(mean(full_loop_output_5q_CL1[,15]), 3)
Full_Boot_Analysis[6,3] <-  round(sd(full_loop_output_5q_CL1[,15]), 3)
Full_Boot_Analysis[6,4] <-  length(which(full_loop_output_5q_CL1[,16] < .1))
Full_Boot_Analysis[6,5] <-  NA

Full_Boot_Analysis[7,1] <- "sc_did_coeff"
Full_Boot_Analysis[7,2] <-  round(mean(full_loop_output_5q_CL1[,19]), 3)
Full_Boot_Analysis[7,3] <-  round(sd(full_loop_output_5q_CL1[,19]), 3)
Full_Boot_Analysis[7,4] <-  length(which(full_loop_output_5q_CL1[,20] < .1))
Full_Boot_Analysis[7,5] <-  NA