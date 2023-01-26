##%######################################################%##
#                                                          #
####      Survival Analysis for Shiny                   ####
#                                                          #
##%######################################################%##
# Author: Tom Pham
# Updated: Paul Carvalho (6/10/2022)
# Script and data info: This script performs survival analysis using CJS model
# and outputs both reach per 10km and cumulative estimates to csv for use in
# the Shiny app
# Data source: Data extracted from NOAA ERDDAP data server (oceanview)

library(RMark)
library(tidyverse)
library(rerddap)
library(lubridate)
library(clusterPower) # DO NOT UPDATE TO v0.7, removes function expit()
library(cowplot)
library(leaflet)
library(vroom)

memory.limit(size=56000)

# Set working directory
setwd("~/github/calfishtrack-shiny/data") # PC 
studyid_list <- read_csv("studyid_names.csv")

# Clear cached data in order to retrieve all latest data
cache_delete_all()

# Need to create logit() and expit() functions because they're not updated in clusterPower package
expit <- function(x){
   out <- exp(x) / (1 + exp(x))
   return(out)
}
logit <- function(p){
   out <- log(p/(1-p))
   return(out)
}

### Load TaggedFish and ReceiverDeployments tables through ERDDAP ---------------------------------------------------------------------

# Establish ERDDAP url and database name
my_url     <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo  <- info('FED_JSATS_taggedfish', url = my_url)
TaggedFish <- tabledap(JSATSinfo,
                       fields = c('fish_id','study_id','tag_id_hex','release_location','release_river_km','fish_release_date','release_latitude',
                                  'release_longitude'),
                       url = my_url) %>% as_tibble()

JSATSinfo           <- info('FED_JSATS_receivers', url = my_url)
ReceiverDeployments <- tabledap(JSATSinfo, 
                                fields = c('receiver_serial_number','receiver_general_location','receiver_location','receiver_river_km',
                                           'receiver_start','receiver_end','receiver_last_valid','receiver_general_river_km','receiver_region',
                                           'receiver_general_latitude','receiver_general_longitude'),
                                url = my_url) %>% as_tibble()

my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo <- info('FED_JSATS_detects', url = my_url)

# Retrieve list of all studyIDs on FED_JSATS
studyid_list <- tabledap(JSATSinfo,
                         fields = c('study_id'),
                         url = my_url,
                         distinct = TRUE) %>% 
                  filter(study_id != "2017_BeaconTag")


# Load functions ----------------------------------------------------------
get_detections <- function(studyID) {
  # Retrieve detection data from ERDDAP
  #
  # Arguments:
  #  studyID: StudyID name to retrieve data for
  #     
  # Return:
  #  df of detection data formatted correctly, add in RKM, Region, Lat, Lon, 
  # Release RKM, Release Lat, Release Lon, format types, rename cols
  
  df <- tabledap(JSATSinfo,
                 fields = c('study_id', 'fish_id', 'receiver_general_location',
                            'time'),
                 paste0('study_id=', '"',studyID, '"'),
                 url = my_url,
                 distinct = T) %>% 
    left_join(ReceiverDeployments %>% select(receiver_general_location, receiver_general_river_km, receiver_region, 
                                             receiver_general_latitude, receiver_general_longitude)) %>% 
    left_join(TaggedFish %>% select(fish_id, release_river_km, release_latitude, release_longitude, release_location) %>% distinct()) %>% 
    distinct()
  
  # Rename columns and change column types as ERDDAP returns data all in 
  # character format
  df <- df %>% 
    rename(
      StudyID = study_id,
      FishID = fish_id,
      GEN = receiver_general_location,
      GenRKM = receiver_general_river_km,
      Region = receiver_region,
      GenLat = receiver_general_latitude,
      GenLon =receiver_general_longitude,
      RelRKM = release_river_km,
      Rel_loc = release_location
    ) %>% 
    mutate(
      GenLat = ifelse(is.na(GenLat), release_latitude, GenLat),
      GenLon = ifelse(is.na(GenLon), release_longitude, GenLon),
      GenLat = as.numeric(GenLat),
      GenLon = as.numeric(GenLon),
      GenRKM = as.numeric(GenRKM),
      RelRKM = as.numeric(RelRKM),
      time = ymd_hms(time),
      GenRKM = ifelse(is.na(GenRKM), RelRKM, GenRKM)
    ) %>% 
    as_tibble() # ERDDAP by default returns a table.dap object which does not play nice with
  # maggittr (pipes) so convert to tibble
  
  # Check for duplicate GEN with different GenRKM, if found replace with mean 
  # GenRKM, GenLat, GenLon
  dup_GEN <- df %>% 
    distinct(GEN, GenRKM, GenLat, GenLon) %>% 
    group_by(GEN) %>% 
    summarise_at(
      c("GenRKM", "GenLat", "GenLon"), mean
    )
  
  # Replace any duplicated GEN with mean values
  df <- df %>% 
    rowwise() %>% 
    mutate(
      GenRKM = ifelse(GEN %in% dup_GEN$GEN, dup_GEN$GenRKM[dup_GEN$GEN == GEN], 
                      GenRKM),
      GenLat = ifelse(GEN %in% dup_GEN$GEN, dup_GEN$GenLat[dup_GEN$GEN == GEN], 
                      GenLat),
      GenLon = ifelse(GEN %in% dup_GEN$GEN, dup_GEN$GenLon[dup_GEN$GEN == GEN], 
                      GenLon)
    )
  
}

get_receiver_GEN <- function(all_detections) {
  # Get a list of all receiver sites and metadata for a given detections df
  #
  # Arguments:
  #  all_detections: detections df 
  #
  # Return:
  #  df of receiver sites along with RKM, Lat, Lon, Region
  
  reach.meta <- all_detections %>% 
    bind_rows() %>% 
    distinct(GEN, GenRKM, GenLat, GenLon, Region) %>% 
    # Necessary because detections files shows differing RKM, Lat, Lon for some 
    # GEN sometimes
    group_by(GEN) %>% 
    summarise(
      GenRKM = mean(GenRKM),
      GenLat = mean(GenLat),
      GenLon = mean(GenLon),
      Region = first(Region)
    ) %>% 
    arrange(desc(GenRKM))
}

aggregate_GEN <- function(detections, 
                          replace_dict = list(replace_with = list(c("Chipps"),
                                                                  c("Benicia"),
                                                                  c("SacTrawl")),
                                              replace_list = list(c("ChippsE", "ChippsW"),
                                                                  c("BeniciaE", "BeniciaW"),
                                                                  c("SacTrawl1", "SacTrawl2")))) {
  # Replace GEN in detections df according to replace_list, basically aggregate
  # sites into one. By default this is done for ChippsE/W, BeniciaE/W, and
  # SacTrawl1/2
  #
  # Arguments:
  #  detections: a detections df
  #  replace_dict: list of receiver locations to aggregate and, the aggregated
  #  name
  #     
  # Return:
  #  a detections df that has replaced list of GEN with the aggregated GEN,
  #  mean RKM, mean Lat, mean Lon. Creates reach.meta.aggregate which is the list
  #  of receiver sites with new aggregated GEN's, along with RKM, Lat, Lon
  
  # Make a copy of reach.meta (receiver metadata)
  reach.meta.aggregate <<- reach.meta
  
  # Walk through each key/pair value
  for (i in 1:length(replace_dict$replace_with)) {
    # Unlist for easier to use format
    replace_list <- unlist(replace_dict[[2]][i])
    replace_with <- unlist(replace_dict[[1]][i])
    
    # Gather receiver data for the replace_list and replace, get the mean genrkm. 
    # This will be used to replace in the detections
    replace <- reach.meta %>% 
      select(GEN, GenRKM, GenLat, GenLon, Region) %>% 
      filter(GEN %in% c(replace_list, replace_with)) %>%
      distinct() %>% 
      select(-GEN) %>% 
      group_by(Region) %>% 
      summarise_all(mean)
    
    # Replace replace_list GENs name with replace_with GEN, and replace all of 
    # their genrkm with the averaged val
    detections <- detections %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
      )
    
    # This new df shows receiver metadata and reflects the aggregation done
    reach.meta.aggregate <<- reach.meta.aggregate %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
        Region = ifelse(GEN == "End", "End", ifelse(GEN %in% c(replace_list, replace_with), replace$Region, Region))
      ) %>% 
      distinct()
  }
  detections
}


make_EH <- function(detections_tmp) {
  # Make an encounter history df
  #
  # Arguments:
  #  detections: a detections df
  #     
  # Return:
  #  Encounter history df. A matrix of every fish tagged for a given studyID
  #  at every given receiver site (that is in reach.meta.aggregate) and whether
  #  it was present 1 or absent 0 in the detection df
  
  # Get earliest detection for each fish at each GEN
  min_detects <- detections_tmp %>% 
    filter(GEN %in% reach.meta.aggregate$GEN) %>% 
    group_by(FishID, GEN, GenRKM) %>% 
    summarise(min_time = min(time)) %>% 
    arrange(FishID, min_time) %>%
    ungroup()
  
  # Get list of all tagged fish for the studyID
  fish <- TaggedFish %>% 
    filter(study_id == detections_tmp$StudyID[1]) %>% 
    arrange(fish_id) %>% 
    pull(fish_id)
  
  # Create matrix of all combinations of fish and GEN
  EH <- expand.grid(
    fish,
    reach.meta.aggregate$GEN, stringsAsFactors = FALSE 
  )
  
  names(EH) <- c('FishID', 'GEN')  
  
  # Add col detect to min_detects, these fish get a 1
  min_detects$detect <- 1
  
  # Join in detections to the matrix, fish detected a GEN will be given a 1
  # otherwise it will be given a 0
  EH <- EH %>% 
    left_join(
      min_detects %>% 
        select(
          FishID, GEN, detect
        ), by = c("FishID", "GEN")
    ) %>% 
    # Replace NA with 0 https://stackoverflow.com/questions/28992362/dplyr-join-define-na-values
    mutate_if(
      is.numeric, coalesce, 0
    )
  
  # Reshape the df wide, so that columns are GEN locations, rows are fish, 
  # values are 1 or 0 for presence/absence
  EH <- reshape(EH, idvar = 'FishID', timevar = 'GEN', direction = 'wide')
  colnames(EH) <- gsub('detect.', '', colnames(EH))
  # Manually make the release column a 1 because all fish were released there
  # sometimes detections df does not reflect that accurately
  EH[2] <- 1
  EH
}

create_inp <- function(detections_tmp, EH) { 
  # Create an inp df
  #
  # Arguments:
  #  detections: a detections df (make sure to use the aggregated version)
  #  EH: an encounter history df
  #     
  # Return:
  #  inp df i.e. Fish01 | 11101, a record of a fish and it's presence/absence
  #  at each given receiver location. 
  
  EH.inp <- EH %>% 
    # Collapse the encounter columns into a single column of 1's and 0's
    unite("ch", 2:(length(EH)), sep ="") %>% 
    # Use the detections df to get the StudyID assignment
    mutate(StudyID = unique(detections_tmp$StudyID))
  EH.inp
}


get_mark_model <- function(all.inp, standardized, multiple, multi_model = times) {
  # Run a CJS Mark model
  #
  # Arguments:
  #  all.inp: inp df, can be more than one studyID
  #  standardized: TRUE or FALSE, if you want outputs to be standardized to 
  #     per10km or not
  #  multiple: TRUE or FALSE, if you have multiple studyIDs or not   
  #  multi_model: model type to use for cases of multiple StudyIDs. By default
  #     multi_model = times runs reach*StudyID, alternatively 
  #     multi_model = plus runs reach+StudyID
  #
  # Return:
  #  the outputs of running a CJS Mark model, df with phi and p estimates, LCI
  #  UCI, SE
  
  # For single studyID
  if (multiple == F) {
    # If standardized, set time.intervals to reach_length to get per 10km
    if (standardized) {
      all.process <- process.data(all.inp, model="CJS", begin.time=1, 
                                  time.intervals = reach_length)
    } else {
      all.process <- process.data(all.inp, model="CJS", begin.time=1)
    }
    # For multiple studyID
  } else {
    # If multiple studyIDs, set groups to "StudyID
    if (standardized) {
      all.process <- process.data(all.inp, model="CJS", begin.time=1, 
                                  time.intervals = reach_length, groups = "StudyID")
    } else {
      all.process <- process.data(all.inp, model="CJS", begin.time=1,
                                  groups = "StudyID")
    }
  }
  
  
  all.ddl <- make.design.data(all.process)
  rm(list=ls(pattern="p.t.x.y"))
  rm(list=ls(pattern="Phi.t.x.y"))
  
  # Set the model up differently depending on single vs multiple StudyIDs
  # and multi_model argument (times vs plus)
  if (multiple) {
    if (multi_model == times) {
      p.t.x.y <- list(formula= ~time*StudyID)
      Phi.t.x.y <- list(formula= ~time*StudyID)
    }else {
      p.t.x.y <- list(formula= ~time+StudyID)
      Phi.t.x.y <- list(formula= ~time+StudyID)
    }
    
  }else {
    p.t.x.y <- list(formula= ~time)
    Phi.t.x.y <- list(formula= ~time)
  }
  
  cml = create.model.list("CJS")
  model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl) 
  outputs <- model.outputs$Phi.t.x.y.p.t.x.y$results$real
  
}

get_cum_survival <- function(all.inp, add_release) {
  # Run a CJS Mark model for cumulative survival
  #
  # Arguments:
  #  all.inp: inp df, can be more than one studyID
  #  add_release: TRUE or FALSE, if you wish to add an extra dummy row at the top
  #  to show 100% survival at the release location
  #
  # Return:
  #  Cumulative survival outputs of CJS Mark model
  
  all.process <- process.data(all.inp, model = "CJS", begin.time = 1)
  all.ddl <- make.design.data(all.process)
  
  rm(list=ls(pattern="Phi.t"))
  rm(list=ls(pattern="p.t"))
  
  p.t <- list(formula= ~time) 
  Phi.t <- list(formula= ~time)
  
  cml = create.model.list("CJS")
  
  model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl, realvcv = TRUE)
  
  reaches <- nchar(all.inp$ch[1]) - 1
  
  phi.t <- model.outputs$Phi.t.p.t$results$real$estimate[1:reaches] 
  phi.t.vcv <- model.outputs$Phi.t.p.t$results$real.vcv
  
  cum.phi <- cumprod(phi.t)
  
  # calculate standard errors for the cumulative product. 
  cum.phi.se <- deltamethod.special("cumprod", phi.t[1:reaches], 
                                    phi.t.vcv[1:(reaches),1:(reaches)])
  
  
  ### Output estimate, SE, LCI, UCI to a dataframe
  cumulative <- data.frame(cum.phi = cum.phi, 
                           cum.phi.se = cum.phi.se,
                           # LCI = cum.phi - 1.96 * cum.phi.se,
                           # UCI = cum.phi + 1.96 * cum.phi.se)
                           LCI = expit(logit(cum.phi) - 1.96 * sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))),
                           UCI = expit(logit(cum.phi) + 1.96 * sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))))
 
  
  # Round to 3 digits
  cumulative <- round(cumulative,3)
  
  # If add_release TRUE, add in the dummy row to the top which just represents
  # survival of 100% at release
  if (add_release == T) {
    cumulative <- cumulative %>% 
      add_row(
        .before = 1,
        cum.phi = 1,
        cum.phi.se = NA,
        LCI = NA, 
        UCI = NA
      ) %>% 
      mutate(
        StudyID = all.inp$StudyID[1]
      )
  }else {
    cumulative <- cumulative %>% 
      mutate(
        StudyID = all.inp$StudyID[1]
      )
  }
  
}

format_phi <- function(outputs, multiple) {
  # Format phi outputs for plotting and table outputs
  #
  # Arguments:
  #  output: output df from Mark model
  #  mutliple: TRUE/FALSE if there were multiple StudyIDs in the outputs
  #
  # Return:
  #  properly formatted df for phi outputs, now ready to plot
  
  
  # Identify number of fish detected at each GEN
  unique_detects <- get_unique_detects(aggregated) 
  # %>% 
  #   filter(StudyID %in% unique(outputs$StudyID))
  
  # Format for single
  if (multiple == F) {
    outputs %>%
      # Grab first half of Mark outputs which represent Phi values
      slice(1:(nrow(outputs) / 2)) %>%
      select(
        -c("fixed", "note")
      ) %>% 
      add_column(
        StudyID = studyID,
        reach_start = reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)],
        reach_end = reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))],
        rkm_start = reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)],
        rkm_end = reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))]
      ) %>% 
      left_join(
        reach.meta.aggregate %>%
          select(GEN, Region) %>%
          distinct(),
        by = c("reach_start" = "GEN")
      ) %>% 
      mutate(
        Reach = paste0(reach_start, " to \n", reach_end),
        RKM = paste0(rkm_start, " to ", rkm_end)
      ) %>% 
      left_join(unique_detects %>% 
                  select(GEN, count_at_start = count),
                by = c("reach_start" = "GEN")) %>% 
      mutate(count_at_start = ifelse(is.na(count_at_start), 0, count_at_start)) %>% 
      left_join(unique_detects %>% 
                  select(GEN, count_at_end = count),
                by = c("reach_end" = "GEN")) %>% 
      mutate(count_at_end = ifelse(is.na(count_at_end), 0, count_at_end)) %>% 
      add_column(
        GenLat_start = reach.meta.aggregate$GenLat[1:(length(reach.meta.aggregate$GenLat)-1)],
        GenLon_start = reach.meta.aggregate$GenLon[1:(length(reach.meta.aggregate$GenLon)-1)],
        GenLat_end = reach.meta.aggregate$GenLat[2:(length(reach.meta.aggregate$GenLat))],
        GenLon_end = reach.meta.aggregate$GenLon[2:(length(reach.meta.aggregate$GenLon))]
      )
  } else {
    outputs %>%
      slice(1:(nrow(outputs) / 2)) %>%
      select(
        -c("fixed", "note")
      ) %>% 
      rownames_to_column(var = "StudyID") %>%
      add_column(
        reach_start = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN)-1)], 2),
        reach_end = rep(reach.meta.aggregate$GEN[2:(length(reach.meta.aggregate$GEN))], 2),
        rkm_start = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM)-1)], 2),
        rkm_end = rep(reach.meta.aggregate$GenRKM[2:(length(reach.meta.aggregate$GenRKM))], 2)
      ) %>% 
      left_join(
        reach.meta.aggregate %>%
          select(GEN, Region) %>%
          distinct(),
        by = c("reach_start" = "GEN")
      ) %>% 
      mutate(
        Reach = paste0(reach_start, " to \n", reach_end),
        RKM = paste0(rkm_start, " to ", rkm_end)
      ) %>% 
      rowwise() %>%
      mutate(
        StudyID = strsplit(strsplit(StudyID, "Phi g")[[1]][2], " ")[[1]][1]
      ) %>% 
      left_join(unique_detects %>% 
                  select(GEN, StudyID, count_at_start = count),
                by = c("reach_start" = "GEN", "StudyID")) %>% 
      mutate(count_at_start = ifelse(is.na(count_at_start), 0, count_at_start)) %>% 
      left_join(unique_detects %>% 
                  select(GEN, StudyID, count_at_end = count),
                by = c("reach_end" = "GEN", "StudyID")) %>% 
      mutate(count_at_end = ifelse(is.na(count_at_end), 0, count_at_end)) %>% 
      add_column(
        StudyID = studyID,
        GenLat_start = reach.meta.aggregate$GenLat[1:(length(reach.meta.aggregate$GenLat)-1)],
        GenLon_start = reach.meta.aggregate$GenLon[1:(length(reach.meta.aggregate$GenLon)-1)],
        GenLat_end = reach.meta.aggregate$GenLat[2:(length(reach.meta.aggregate$GenLat))],
        GenLon_end = reach.meta.aggregate$GenLon[2:(length(reach.meta.aggregate$GenLon))]
      )
  }
  
}

get_unique_detects <- function(all_aggregated){
  # Get raw number of unique fish detected at each GEN 
  #
  # Arguments:
  #  all_aggregated: df of detections that have been replaced with aggregating
  #  receiver locations
  #
  # Return:
  #  df of each GEN in a detections df and the raw number of unique fish detected
  
  all_aggregated %>% 
    bind_rows() %>% 
    select(StudyID, FishID, GEN, GenRKM) %>% 
    distinct() %>% 
    group_by(StudyID, GEN, GenRKM) %>% 
    summarise(
      count = n()
    ) %>% 
    arrange(StudyID, desc(GenRKM)) %>% 
    ungroup()
}

make_phi_table <- function(phi, standardized = T) {
  # Format phi outputs further to be ready to save as a csv
  #
  # Arguments:
  #  phi: phi outputs from Mark model, must be formatted with (format_phi) first
  #
  # Return:
  #  phi df formatted the way I want to be saved as csv
  ifelse(standardized, label <-  'Survival rate per 10km (SE)',
         label <- 'Survival rate (SE)')
  
  phi %>% 
    select(studyID, reach_num, Reach, RKM, Region, Estimate = estimate, SE = se, 
           LCI = lcl, UCI = ucl) %>% 
    mutate(
      Reach = str_remove_all(Reach, "\n"),
      Estimate = round(Estimate, 2),
      SE = round(SE, 2),
      LCI = round(LCI, 2),
      UCI = round(UCI, 2),
      Estimate2 = paste0(Estimate, " (", as.character(SE), ")")
    ) %>% 
    rename(!!label := Estimate2,
           'Reach #' = reach_num)
}


format_cum_surv <- function(cum_survival_all) {
  # Format cumulative survival outputs for plotting and table outputs
  #
  # Arguments:
  #  cum_survival_all: output df from Mark model cumulative survival
  #
  # Return:
  #  properly formatted df for phi outputs, now ready to plot
  StudyID = unique(cum_survival_all$StudyID)
   
  cum_survival_all <- cum_survival_all %>% 
    add_column(
      GEN = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN))], 
                length(StudyID)),
      RKM = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM))], 
                length(StudyID)),
      reach_num = rep(seq(0, (nrow(reach.meta.aggregate))-1, 1), length(StudyID))
    ) %>% 
    left_join(
      reach.meta.aggregate %>%
        select(GEN, Region) %>% 
        distinct(),
      by = c("GEN")
    ) %>% 
    mutate(
      cum.phi = round(cum.phi, digits = 2),
      cum.phi.se = round(cum.phi.se, digits = 2),
      LCI = round(LCI, digits = 2),
      UCI = round(UCI, digits = 2),
      'Survival estimate (SE)' = paste0(cum.phi, " (", as.character(cum.phi.se), ")"),
      'Reach #' = reach_num,
    ) %>% 
    filter(
      GEN != "GoldenGateW"
    ) 
  
}

plot_cum_surv <- function(cum_survival_all, add_breaks, multiple, padding = 5.5,
                          text_size = 20) {
  # Plot cumulative survival outputs from Mark model
  #
  # Arguments:
  #  cum_survival_all: cumulative survival outputs from Mark model, 
  #     must be formatted with (format_cum_surv) first
  #  add_breaks: TRUE/FALSE whether to add vertical line breaks to represent
  #     regions
  #  multiple: TRUE/FALSE, whether there are multiple studyIDs or not
  #  padding: leftside plot margin, default set to 5.5 good for most, but can
  #     be adjusted of the xaxis label too long and gets cut off
  #
  # Return:
  #  plot of cumulative survival with estimate and error bars representing LCI, UCI
  
  # Create the levels order 
  lvls <- cum_survival_all %>% 
    select(GEN) %>% 
    distinct() %>% 
    pull()
  
  if (multiple) {
    cum_survival_all %>% 
      mutate(
        GEN = factor(GEN, levels = lvls)
      ) %>%
      ggplot(mapping = aes(x = GEN, y = cum.phi, group = studyID)) +
      geom_point(size = 2, aes(color = studyID)) +
      geom_errorbar(mapping = aes(x= GEN, ymin = LCI, ymax = UCI, 
                                  color = studyID),  width = .1) +
      geom_line(size = 0.7, aes(color = studyID)) +
      {if(add_breaks)geom_vline(xintercept = region_breaks, linetype = "dotted")} +
      ylab("Cumulative survival") +
      xlab("Receiver Location") +
      scale_y_continuous(breaks = seq(0, 1, 0.1)) +
      # scale_x_continuous(breaks = 0:max(cum_survival_all$reach_num)) +
      # NEEDS TO BE FIXED FOR IF THERE ARE MORE THEN 2 STUDYIDS
      scale_color_manual(values=c("#007EFF", "#FF8100")) +
      theme_classic() +
      theme(
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.margin = margin(5.5, 5.5, 5.5, padding, "pt"),
        legend.position = "top",
        text = element_text(size=text_size)
      ) 
  } else {
    cum_survival_all %>% 
      mutate(
        # GEN = factor(GEN, levels = reach.meta.aggregate$GEN,
        #              labels = paste0(reach.meta.aggregate$GEN, " (",
        #                              reach.meta.aggregate$GenRKM, ")"))
        GEN = factor(GEN, levels = lvls)
      ) %>%
      ggplot(mapping = aes(x = GEN, y = cum.phi)) +
      geom_point(size = 2) +
      geom_errorbar(mapping = aes(x= GEN, ymin = LCI, ymax = UCI),  width = .1) +
      geom_line(size = 0.7, group = 1) +
      {if(add_breaks)geom_vline(xintercept = region_breaks, linetype = "dotted")} +
      ylab("Cumulative survival") +
      xlab("Receiver Location") +
      scale_y_continuous(breaks = seq(0, 1, 0.1)) +
      # scale_x_continuous(breaks = 0:max(cum_survival_all$GEN)) +
      theme_classic() +
      theme(
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.margin=unit(c(.5,.5,.5,.5), "cm"),
        text = element_text(size=text_size)
      ) 
  }
}

pivot_cum_surv <- function (p){
  df <- p %>% 
    pivot_wider(names_from = studyID, values_from = c("Survival estimate (SE)", 
                                                      "LCI", "UCI"),
                names_glue = "{studyID} {.value}")
  
  prefixes <- unique(cum_survival_all$studyID)
  
  names_to_order <- map(prefixes, ~ names(df)[grep(paste0(.x, " "), names(df))]) %>% unlist
  names_id <- setdiff(names(df), names_to_order)
  
  df <- df %>%
    select(names_id, names_to_order)
}

# Identify StudyIDs with multiple release locations -----------------------

find_multi_release_loc <- function(studyid_list) {
  rel_loc <- TaggedFish %>% 
    filter(study_id %in% studyid_list) %>% 
    select(study_id, release_location) %>% 
    distinct() %>% 
    group_by(study_id) %>% 
    filter(n() > 1)

}

# Identify the studyIDs with multiple releases and their locations
multi_rel_loc <- find_multi_release_loc(studyid_list$study_id)

# find_multi_release_date <- function(studyid_list) {
#   rel_loc <- TaggedFish %>% 
#     filter(study_id %in% studyid_list) %>% 
#     select(study_id, fish_release_date) %>% 
#     mutate(fish_release_date = as.Date(mdy_hms(fish_release_date))) %>% 
#     distinct() %>% 
#     group_by(study_id) %>% 
#     filter(n() > 1)
# }
# 
# multi_rel_date <- find_multi_release_date(studyid_list)


### Output survival estimates for multi-release location groups -------------
# "FR_Spring_2013", "FR_Spring_2014", "FR_Spring_2015", "FR_Spring_2019" ,"FR_Spring_2020" ,
# "Nimbus_Fall_2018", "ColemanAltRel_2021", "ColemanAltRel_2022", "SJ_steelhead_2021", "SJ_steelhead_2022"

# Select StudyID
studyID <- "SJ_Steelhead_2021"

# Get detections
detections <- get_detections(studyID)

# Get list of GEN
reach.meta <- get_receiver_GEN(detections)

# Remove some GEN used for survival analysis, should be linear path, no upstream
# In general, ignore Delta sites, more specific removal according to each studyID

# Feather River studies - remove any Sac GEN above Blw_FRConf
# if (str_detect(studyID, "FR_Spring")) {
#   reach.meta <- reach.meta %>% 
#     filter(
#       !Region %in% c("North Delta", "East Delta", "West Delta", "South Delta", 
#                      "Yolo Bypass", "Franks Tract", "Upper Sac R") |
#         GEN %in% c("ChippsE", "ChippsW"),
#       !(Region == "Lower Sac R" & GenRKM > 203.46)
#     )
# } else if (studyID == "CNFH_FMR_2019") {
#   reach.meta <- reach.meta %>% 
#     filter(
#       !Region %in% c("North Delta", "East Delta", "West Delta", "South Delta", 
#                      "Yolo Bypass", "Franks Tract") |
#         GEN %in% c("ChippsE", "ChippsW"),
#       GEN != "Butte1"
#     )
# } else if (studyID == "Nimbus_Fall_2018") {
#   reach.meta <- reach.meta %>% 
#     filter(
#       !Region %in% c("North Delta", "East Delta", "West Delta", "South Delta", 
#                      "Yolo Bypass", "Franks Tract", "Upper Sac R") |
#         GEN %in% c("ChippsE", "ChippsW"),
#       !(Region == "Lower Sac R" & GenRKM > 172.00)
#     )
# } else {
#   reach.meta <- reach.meta %>%
#     filter(
#       !Region %in% c("North Delta", "East Delta", "West Delta", "South Delta",
#                      "Yolo Bypass", "Franks Tract") |
#         GEN %in% c("ChippsE", "ChippsW"),
#       GenRKM <= TaggedFish %>%
#         filter(study_id == studyID) %>%
#         mutate(release_river_km = as.numeric(release_river_km)) %>%
#         group_by(study_id) %>%
#         summarise(max_rkm = max(release_river_km)) %>%
#         pull(max_rkm)
#     )
# }

# Filter out receiver locations for SJ Steelhead studies
reach.meta <- reach.meta %>%
              filter(Region %in% c("Carquinez Strait", "San Joaquin River ", "SF Bay") |
                     GEN %in% c("OR_HOR_US", "OR_HOR_DS", "ChippsE", "ChippsW", "Stockton_Rel", "Head_of_Old_River_Rel")) %>%
              filter(GEN != "Las Palmas" & GEN != "DurhamFerryUS_1")
reach.meta <- rbind(reach.meta, data.frame(GEN = "Head_of_Old_River_Rel", GenRKM = 156, GenLat = 37.80789, GenLon = -121.3295, Region = NA))

# # Filter out receiver locations in the Delta, and any that are above the 
# # release RKM
# reach.meta <- reach.meta %>%
#   filter(
#     !Region %in% c("North Delta", "East Delta", "West Delta", "South Delta",
#                    "Yolo Bypass", "Franks Tract") |
#       GEN %in% c("ChippsE", "ChippsW"),
#     GenRKM <= TaggedFish %>%
#       filter(study_id == studyID) %>%
#       mutate(release_river_km = as.numeric(release_river_km)) %>%
#       group_by(study_id) %>%
#       summarise(max_rkm = max(release_river_km)) %>%
#       pull(max_rkm)
#   )

# Visually inspect receiver locations, determine if sites need to be removed
leaflet(data = reach.meta) %>% 
  addTiles() %>% 
  addMarkers(lng = ~GenLon, lat = ~GenLat, label = ~GEN, 
             labelOptions = labelOptions(noHide = T))

# Identify the unique release locations
rel_loc <- TaggedFish %>% 
  filter(study_id == detections$StudyID[1]) %>% 
  select(release_location, release_river_km) %>% 
  distinct() %>% 
  arrange(desc(as.numeric(release_river_km))) %>% 
  pull(release_location)

run_multi_survival <- function(release_loc, type) {
  # If the current release location is not the first begin at its starting
  # position in reach.meta
  if(reach.meta$GenRKM[reach.meta$GEN == release_loc] < max(reach.meta$GenRKM)){
    reach.meta <- reach.meta %>% filter(GenRKM <= reach.meta$GenRKM[reach.meta$GEN == release_loc])
  }

  detects <- detections %>% 
             filter(Rel_loc == release_loc,
                    GEN %in% reach.meta$GEN)

  aggregated <- aggregate_GEN(detects)
  
  EH <- make_EH(aggregated)
  
  inp <- create_inp(aggregated, EH)
  
  # Get river KM
  KM <- reach.meta.aggregate$GenRKM
  
  # Get the reach lengths per 10km
  reach_length <- round((abs(diff(KM))/10), digits = 3)
  
  if (type == "reach") {
    reach_surv <- get_mark_model(inp, standardized = T, multiple = F)
    reach_surv <- format_phi(reach_surv, multiple = F) %>% 
      filter(reach_end != "GoldenGateW") %>% 
      mutate(release = release_loc)
    
    fish_count <- get_unique_detects(aggregated)
    reach_surv <- reach_surv %>% 
      left_join(
        fish_count %>% 
          select(reach_start = GEN, count_at_start = count)
      ) %>% 
      left_join(
        fish_count %>% 
          select(reach_end = GEN, count_at_end = count)
      ) %>% 
      mutate_if(is.numeric, coalesce, 0)
  } else if (type == "cumulative") {
    cum_surv <- get_cum_survival(inp, add_release = T)
    cum_surv <- format_cum_surv(cum_surv)
    cum_surv <- cum_surv %>% 
      mutate(release = release_loc)
    
    fish_count <- get_unique_detects(aggregated)
    cum_surv <- cum_surv %>%
      left_join(
        fish_count %>%
          select(GEN, count)
      ) %>%
      left_join(
        reach.meta.aggregate %>% 
          select(GEN, GenLat, GenLon)
      ) %>% 
      mutate_if(is.numeric, coalesce, 0)
  } else {
    print("Type must be 'reach' or 'cumulative")
  }
  
}

reach_surv <- lapply(rel_loc, run_multi_survival, type = "reach") %>% 
  bind_rows()

### CODE FOR RUNNING SAN JOAQUIN STEELHEAD STUDIES ####
reach.meta$GenRKM[which(reach.meta$GEN == "Head_of_Old_River_Rel")] <- reach.meta$GenRKM[which(reach.meta$GEN == "Mossdale")] - 4
reach.meta$GenRKM[which(reach.meta$GEN == "OR_HOR_US")] <- reach.meta$GenRKM[which(reach.meta$GEN == "Head_of_Old_River_Rel")] - 1.16
reach.meta$GenRKM[which(reach.meta$GEN == "OR_HOR_DS")] <- reach.meta$GenRKM[which(reach.meta$GEN == "OR_HOR_US")] - 0.2
reach.meta <- reach.meta %>% arrange(desc(GenRKM))
reach.meta.main <- reach.meta 
reach.meta.1 <- reach.meta %>%
                filter(!Region %in% c("South Delta", NA))
reach.meta.2 <- reach.meta %>%
                filter(!Region %in% c("San Joaquin River ") &
                       !GEN %in% c("Stockton_Rel"))
reach.meta.3 <- reach.meta %>%
                filter(!Region %in% c("San Joaquin River ") &
                       !GEN %in% c("OR_HOR_US", "OR_HOR_DS", "Head_of_Old_River_Rel"))

reach.meta <- reach.meta.main
# reach_surv1 <- reach_surv
# reach_surv2 <- reach_surv
# reach_surv3 <- reach_surv
reach_surv <- rbind(reach_surv1, reach_surv2, reach_surv3)
# cum_surv1 <- cum_surv
# cum_surv2 <- cum_surv
# cum_surv3 <- cum_surv
cum_surv <- rbind(cum_surv1, cum_surv2, cum_surv3)

# Rerun first release group just to restore reach.meta.aggregate to full version
run_multi_survival(rel_loc[1], type = "reach")

reach.meta.aggregate <- reach.meta.aggregate %>%
   mutate(reach_num = 1:length(GEN))

# Fix reach_num for releases after the first
reach_surv <- reach_surv %>% 
  # rowwise() %>% 
  # mutate(reach_num = ifelse(release != rel_loc[1], 
  #                           reach_surv$reach_num[reach_surv$reach_start == reach_start & reach_surv$release == rel_loc[1]],
  #                           1:length(reach_surv$Reach)# reach_num
  #   )
  # ) %>% 
  left_join(
    reach.meta.aggregate %>% 
      select(GEN, GenLat_start = GenLat, GenLon_start = GenLon, reach_num),
    by = c("reach_start" = "GEN")
  ) %>% 
  left_join(
    reach.meta.aggregate %>% 
      select(GEN, GenLat_end = GenLat, GenLon_end = GenLon),
    by = c("reach_end" = "GEN")
  )

# Format SJ reach survival
reach_surv <- reach_surv %>%
   select(estimate, se, lcl, ucl, StudyID, reach_start, reach_end, rkm_start, rkm_end, Region, Reach, RKM, reach_num, release, 
          count_at_start, count_at_end, GenLat_start = GenLat_start.x, GenLon_start = GenLon_start.x, 
          GenLat_end = GenLat_end.x, GenLon_end = GenLon_end.x)

write_csv(reach_surv, paste0("Survival/Reach Survival Per 10km/", studyID, "_reach_survival.csv"))

cum_surv <- lapply(rel_loc, run_multi_survival, type = "cumulative") %>% 
  bind_rows()

run_multi_survival(rel_loc[1], type = "cumulative")

cum_surv <- cum_surv %>% 
  left_join(reach.meta.aggregate %>% 
              select(GEN, GenLat, GenLon, reach_n = reach_num))
cum_surv <- cum_surv %>% 
   select(cum.phi, cum.phi.se, LCI, UCI, StudyID, GEN, RKM, reach_num = reach_n, Region, release, count, GenLat, GenLon)

# Fix reach_num for releases after the first
cum_surv <- cum_surv %>% 
  rowwise() %>%
  mutate(
    reach_num = ifelse(
      release != rel_loc[1],
      cum_surv$reach_num[cum_surv$GEN == GEN
                              & cum_surv$release == rel_loc[1]],
      reach_num
    )
  )

write_csv(cum_surv, paste0("Survival/Cumulative Survival/", studyID, "_cum_survival.csv"))

cleanup(ask = FALSE)


#### Get Reach per 10km Survival estimates ------------------------------------------

# 1. Provide StudyID

c("DeerCk_Wild_STH_2019", "DeerCk_Wild_STH_2020",
  "MillCk_Wild_STH_2015", "MillCk_Wild_STH_2019",
  "MillCk_Wild_STH_2020", "Upper_Butte_2019", "Upper_Butte_2020")

studyID <- "Winter_H_2022" 

# 2. Read the detections and get the receiver GEN data; if detections file
# has not been saved for the studyID yet, run save_new_detections.R first

# detect_file <- paste0('../detections/', studyID, ".csv")
# detections <- vroom(detect_file)
detections <- get_detections(studyID)

reach.meta <- get_receiver_GEN(detections)

# QA the detections for 
# -duplicate GenRKM with different names
# -GEN that is above the release RKM

# Check for any GenRKM above release RKM for the studyID
reach.meta %>% 
  filter(
    GenRKM > TaggedFish %>% 
      filter(study_id == studyID) %>% 
      distinct(release_river_km) %>% 
      mutate(release_river_km = as.numeric(release_river_km)) %>% 
      pull()
  )

# Visually inspect receiver locations, determine if sites need to be removed
leaflet(data = reach.meta) %>% 
  addTiles() %>% 
  addMarkers(lng = ~GenLon, lat = ~GenLat, label = ~GEN, 
             labelOptions = labelOptions(noHide = T))

# Check for any duplicate GenRKM
reach.meta %>% 
  group_by(GenRKM) %>% 
  filter(n() > 1)

# Remove GEN as necessary
reach.meta <- reach.meta %>%
  filter(!GEN %in% c("BattleCk3","BattleCk4", "SutterWestSacRiver"))

reach.meta <- reach.meta %>% 
  filter(
    !GEN %in% c("Mill_Ck_Conf", "GCID", "GCID_Main", "GCID_Conf"),
    !is.na(Region),
    !Region %in% c("North Delta", "South Delta", "East Delta", "West Delta", 
                   "Yolo Bypass", "Lower Mok R", "Delta") |
      GEN %in% c("ChippsE", "ChippsW")
  )

# 3. Using the replacement dictionary, aggregate the detections
aggregated <- aggregate_GEN(detections)

# 5. Create an encounter history from the aggregated detections
EH <- make_EH(aggregated)

# 6. Create an inp from the given encounter history
inp <- create_inp(aggregated, EH)

# 7. Run reach per 10km survival in Mark using the inp

# Get river KM
KM <- reach.meta.aggregate$GenRKM

# Get the reach lengths per 10km
reach_length <- round((abs(diff(KM))/10), digits = 3)

# Get the estimates
reach_surv <- get_mark_model(inp, standardized = T, multiple = F)
cleanup(ask = F)

# Format the estimates table
reach_surv_formatted <- format_phi(reach_surv, multiple = F)

write_csv(reach_surv_formatted, paste0("Reach Survival per 10km/", studyID, "_reach_survival.csv"))

# Get Cumulative Survival Estimates ---------------------------------------

# 1. Get the cumulative survival estimates
cum_surv <- get_cum_survival(inp, add_release = T)

# 2. Format the estimates table
cum_surv_formatted <- format_cum_surv(cum_surv)

# Add in unique fish counts at each GEN and lat/lon
fish_count <- get_unique_detects(aggregated)

cum_surv_formatted <- cum_surv_formatted %>% 
  left_join(
    fish_count %>% 
      select(GEN, count) %>% 
      distinct()
  ) %>% 
  mutate_if(is.numeric, coalesce, 0) %>% 
  left_join(
    reach.meta.aggregate %>% 
      select(GEN, GenLat, GenLon)
  )

write_csv(cum_surv_formatted, paste0("Cumulative Survival/", studyID, "_cum_survival.csv"))

cleanup(ask = F)




# Check num fish per studyID ----------------------------------------------

studyIDs <- c("DeerCk_Wild_STH_2019", "DeerCk_Wild_STH_2020",
              "MillCk_Wild_STH_2015", "MillCk_Wild_STH_2019",
              "MillCk_Wild_STH_2020", "Upper_Butte_2019", "Upper_Butte_2020")


TaggedFish %>% 
  filter(study_id %in% studyIDs) %>% 
  group_by(study_id) %>% 
  count()

studyID <- "DeerCk_Wild_STH_2019" 

# 2. Read the detections and get the receiver GEN data; if detections file
# has not been saved for the studyID yet, run save_new_detections.R first

detect_file <- paste0('./data/detections/', studyID, ".csv")
detections <- vroom(detect_file)

reach.meta <- get_receiver_GEN(detections)

# QA the detections for 
# -duplicate GenRKM with different names
# -GEN that is above the release RKM

# Check for any GenRKM above release RKM for the studyID
reach.meta %>% 
  filter(
    GenRKM > TaggedFish %>% 
      filter(study_id == studyID) %>% 
      distinct(release_river_km) %>% 
      mutate(release_river_km = as.numeric(release_river_km)) %>% 
      pull()
  )

# Visually inspect receiver locations, determine if sites need to be removed
leaflet(data = reach.meta) %>% 
  addTiles() %>% 
  addMarkers(lng = ~GenLon, lat = ~GenLat, label = ~GEN, 
             labelOptions = labelOptions(noHide = T))


