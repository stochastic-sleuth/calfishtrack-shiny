##%######################################################%##
#                                                          #
####           Detections downloads for Shiny           ####
#                                                          #
##%######################################################%##

## LAST UPDATED: 12/12/2023 - Jessica Frey

library(tidyverse)
library(rerddap) # read tagged fish from Calfishtrack 
library(lubridate)
library(dbplyr)


# Set working directory **IMPORTANT TO SET TO THE LATEST FOLDER**
setwd("C:/Users/Cyril.Michel/Desktop/AT_shiny/calfishtrack-shiny/data/detections")

# Clear cached data in order to retrieve all latest data
cache_delete_all()

# Establish ERDDAP url and database name
my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo <- info('FED_JSATS_taggedfish', url = my_url)
TaggedFish <- tabledap(JSATSinfo, url = my_url) %>% 
  as_tibble()

my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo <- info('FED_JSATS_receivers', url = my_url)
ReceiverDeployments <- tabledap(JSATSinfo, url = my_url)

my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo <- info('FED_JSATS_detects', url = my_url)

#sanity check
JSATSinfo[["variables"]][["variable_name"]]

# Function to download the studyID in the proper format for Shiny
download_detections <- function(studyID) {
  
  # Retrieve detection data from ERDDAP
  #
  # Arguments:
  #  studyID: StudyID name to retrieve data for
  #     
  # Return:
  #  df of detection data formatted correctly, add in RKM, Region, Lat, Lon, 
  # Release RKM, Release Lat, Release Lon, format types, rename cols
  # Get detection data from ERDDAP
  df <- tabledap(JSATSinfo,
                 fields = c('study_id', 'receiver_serial_number', 'fish_id', 
                            'receiver_general_location', 'receiver_location','first_time'),
                 paste0('study_id=', '"',studyID, '"'),
                 url = my_url,
                 distinct = T) 
  
  # Merge detections with receiver deployments regardless of 'receiver_last_valid'
  df <- df %>% 
        left_join(., ReceiverDeployments %>% 
                     select(receiver_serial_number, receiver_location, receiver_general_location, receiver_general_river_km, receiver_region,
                            receiver_general_latitude, receiver_general_longitude, latitude, longitude, receiver_start, receiver_end, receiver_last_valid),
                  by = c("receiver_serial_number", "receiver_general_location", "receiver_location"), relationship = "many-to-many")
  
  # Filter data to make sure first_time was within the start and last_valid datefirst_time
  df <- as.data.frame(df) # for some reason the following command only worked when converted to data.frame without tabledap class
  df <- df %>%
        mutate(first_time = as.POSIXct(first_time, format = '%Y-%m-%d %H:%M:%S', tz = 'Etc/GMT+8'),
               receiver_start = as.POSIXct(receiver_start, format = '%m/%d/%Y %H:%M:%S', tz = 'Etc/GMT+8'),
               receiver_last_valid = as.POSIXct(receiver_last_valid, format = '%m/%d/%Y %H:%M:%S', tz = 'Etc/GMT+8')) %>%
        mutate(fish_id = ifelse((first_time >= receiver_start & first_time <= receiver_last_valid), fish_id, NA)) %>% # if detection first_time was outside of valid datefirst_time for a receiver then assign NA
        filter(!is.na(fish_id))
  
  # Merge with tagged fish dataframe
  df <- df %>%
           left_join(., TaggedFish %>%
                        select(fish_id, release_river_km, release_latitude, release_longitude))
  
  # Rename columns and change column types as ERDDAP returns data all in 
  # character format
  df <- as.data.frame(df) %>% 
    rename(
      studyID = study_id,
      SN = receiver_serial_number,
      FishID = fish_id,
      GEN = receiver_general_location,
      LAT = latitude,
      LON = longitude,
      GenRKM = receiver_general_river_km,
      Region = receiver_region,
      GenLat = receiver_general_latitude,
      GenLon = receiver_general_longitude,
      RelRKM = release_river_km,
      local_time = first_time
    ) %>% 
    mutate(
      GenLat = ifelse(is.na(GenLat), release_latitude, GenLat),
      GenLon = ifelse(is.na(GenLon), release_longitude, GenLon),
      GenLat = as.numeric(GenLat),
      GenLon = as.numeric(GenLon),
      GenRKM = as.numeric(GenRKM),
      RelRKM = as.numeric(RelRKM),
      GenRKM = ifelse(is.na(GenRKM), RelRKM, GenRKM),
      time_pst = ymd_hms(local_time)
    ) %>% 
    select(c(studyID, SN, FishID, GEN, local_time, LAT, LON, GenRKM, Region, GenLat, GenLon, time_pst))
  
  write_csv(df, paste0(studyID, ".csv"))
  
}

# Look for new studies
unique(TaggedFish$study_id)[grep("2023", unique(TaggedFish$study_id))]

new_studies <- c('BC_Jumpstart_2023', 'Butte_Sink_2023', 'DeerCk_Wild_STH_2023',
                 'FR_Fall_Chinook_2023', 'FR_Spring_2023', 'Lower_Yuba_FRH_Chinook_2023',
                 'MillCk_Wild_STH_2023', 'Mok_Fall_2023', 'Nimbus_Fall_2023',
                 'Putah_Creek_CHN_2023', 'SacRiverSpringJPE_2023', 'SB_Spring_2023',
                 'Seasonal_Survival_2023', 'Spring_Pulse_2023', 'Stan_Steelhead_2023',
                 'Sutter_Bypass_Rice_Subyearling_2023', 'Sutter_Bypass_Rice_Yearling_2023',
                 'Winter_H_2023')

# Provide the studyID name you want to download here
new_studies <- "Wild_stock_Chinook_RBDD_2022"

# Call the download function and provide the studyID you want and it will
# retrieve the data from ERDDAP and save it to csv in the detections folder
for(i in new_studies){
  print(i)
  download_detections(i)
}


# Wild_stock_Chinook_RBDD_2022<- read.csv("/Users/jessefrey/Desktop/calfishtrack/calfishtrack-shiny/data/detections/Wild_stock_Chinook_RBDD_2022.csv")
# head(Wild_stock_Chinook_RBDD_2022)
# Wild_stock_Chinook_RBDD_2021<- read.csv("/Users/jessefrey/Desktop/calfishtrack/calfishtrack-shiny/data/detections/Wild_stock_Chinook_RBDD_2021.csv")
# head(Wild_stock_Chinook_RBDD_2021)
# Juv_Green_Sturgeon_2017<- read.csv("/Users/jessefrey/Desktop/calfishtrack/calfishtrack-shiny/data/detections/Juv_Green_Sturgeon_2017.csv")
# head(Juv_Green_Sturgeon_2017)
# 
# 
# Winter_H_2021<- read.csv("/Users/jessefrey/Desktop/calfishtrack/calfishtrack-shiny/data/detections/Winter_H_2021.csv")
# head(Winter_H_2021)

# We need to provide a fresh vemco deployment csv re: line 125 in app.R 

library(DBI)
library(odbc) #Read in Access database off of SQL server

con <- dbConnect(odbc(),
                 Driver = "ODBC Driver 17 for SQL Server",
                 Server = "calfishtrack-server.database.windows.net",
                 Database = "JSATS_Database",
                 UID = "jsats_user",
                 PWD = "Pass@123",
                 Port = 1433)
odbcListObjects(con)

Vemco_Dep<-dbReadTable(con, "VemcoReceiverDeployments")
#remove 999 named location
Vemco_Dep<- Vemco_Dep[Vemco_Dep$GPSname != "999", ]

Rec_Loc<-dbReadTable(con, "ReceiverLocations")
#remove 999 named location
Rec_Loc<- Rec_Loc[Rec_Loc$GPSname != "999", ]

vloc<-subset(Rec_Loc, select = c("GPSname","GeneralLocation", "GenLat", "GenLon"))
vloc<- vloc %>% rename(GEN = GeneralLocation)

Vemco_Dep<-merge(vloc,Vemco_Dep, all.y = TRUE)

Vemco_Dep<- subset(Vemco_Dep, select = c("GPSname", "GEN", "VemcoSN", "StartTime", "EndTime", "GenLat", "GenLon"))


write.csv(Vemco_Dep, "~/Desktop/calfishtrack/calfishtrack-shiny/data/vemco_rec_deployments.csv")
