##%######################################################%##
#                                                          #
####           Detections downloads for Shiny           ####
#                                                          #
##%######################################################%##

## LAST UPDATED: 03/07/2022 - Paul Carvalho

library(tidyverse)
library(rerddap)
library(lubridate)

# Set working directory **IMPORTANT TO SET TO THE LATEST FOLDER**
setwd("C:/Users/pgcar/Google Drive/1 Work/1 NOAA UCSC AT/1 Projects/Shiny App/shiny_012023/data/detections")

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
                            'receiver_general_location', 'receiver_location','local_time'),
                 paste0('study_id=', '"',studyID, '"'),
                 url = my_url,
                 distinct = T) 
  
  # Merge detections with receiver deployments regardless of 'receiver_last_valid'
  df <- df %>% 
        left_join(., ReceiverDeployments %>% 
                     select(receiver_serial_number, receiver_location, receiver_general_location, receiver_general_river_km, receiver_region,
                            receiver_general_latitude, receiver_general_longitude, latitude, longitude, receiver_start, receiver_end, receiver_last_valid),
                  by = c("receiver_serial_number", "receiver_general_location", "receiver_location"))
  
  # Filter data to make sure local_time was within the start and last_valid datetime
  df <- as.data.frame(df) # for some reason the following command only worked when converted to data.frame without tabledap class
  df <- df %>%
        mutate(local_time = as.POSIXct(local_time, format = '%Y-%m-%d %H:%M:%S', tz = 'Etc/GMT+8'),
               receiver_start = as.POSIXct(receiver_start, format = '%m/%d/%Y %H:%M:%S', tz = 'Etc/GMT+8'),
               receiver_last_valid = as.POSIXct(receiver_last_valid, format = '%m/%d/%Y %H:%M:%S', tz = 'Etc/GMT+8')) %>%
        mutate(fish_id = ifelse((local_time >= receiver_start & local_time <= receiver_last_valid), fish_id, NA)) %>% # if detection time was outside of valid datetime for a receiver then assign NA
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
      RelRKM = release_river_km
    ) %>% 
    mutate(
      GenLat = ifelse(is.na(GenLat), release_latitude, GenLat),
      GenLon = ifelse(is.na(GenLon), release_longitude, GenLon),
      GenLat = as.numeric(GenLat),
      GenLon = as.numeric(GenLon),
      GenRKM = as.numeric(GenRKM),
      RelRKM = as.numeric(RelRKM),
      time_pst = ymd_hms(local_time),
      GenRKM = ifelse(is.na(GenRKM), RelRKM, GenRKM)
    ) %>% 
    select(c(studyID, SN, FishID, GEN, local_time, LAT, LON, GenRKM, Region, GenLat, GenLon, time_pst))
  
  write_csv(df, paste0(studyID, ".csv"))
  
}

# Look for new studies
unique(TaggedFish$study_id)[grep("2022", unique(TaggedFish$study_id))]

# Provide the studyID name you want to download here
studyID <- "SJ_steelhead_2022"

# Call the download function and provide the studyID you want and it will
# retrieve the data from ERDDAP and save it to csv in the detections folder
download_detections(studyID)


