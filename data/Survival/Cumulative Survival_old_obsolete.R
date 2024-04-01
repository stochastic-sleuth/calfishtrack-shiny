##### Cumulative Survival
# Create cumulative survival for 'study groups', i.e. Winter_H, ColemanFall
# using bomber sites. Format and output results to csv for use in Shiny app

# 1. Load in proper bomber sites
# 2. Get detections data for study group through ERDDAP
# 3. Create EH from detections, using the bomber sites
# 4. Run MARK Phi.t.x.y p.t.x.y model 

library(RMark)
library(tidyverse)
library(RODBC)
library(rerddap)
library(lubridate)
library(clusterPower)


### 1. Load in proper bomber sites ------------------------------------------

# Choose the 'group' to process, i.e. 'Winter', 'ColemanFall'
group <- "Wild_stock_Chinook_Rbdd_2021"

## Load river km and receiver names.  These should be all the receiver locations in common for the specific inp and study group you are running
reach.meta <- read_csv(paste0("./outputs/Bomber Sites/", group, "_bomber_sites.csv"))

### 2.  Get detections data for study group through ERDDAP ---------------------------------------------------------------------
# Get TaggedFish and ReceiverDeployments tables

my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo <- info('FED_JSATS_taggedfish', url = my_url)
TaggedFish <- tabledap(JSATSinfo,
                       fields = c('fish_id','study_id','tag_id_hex','release_location','release_river_km','fish_release_date'),
                       url = my_url) %>% as_tibble()

JSATSinfo <- info('FED_JSATS_receivers', url = my_url)
ReceiverDeployments <- tabledap(JSATSinfo, 
                                fields = c('receiver_serial_number','receiver_general_location','receiver_location','receiver_river_km',
                                           'receiver_start','receiver_end','receiver_last_valid','receiver_general_river_km'),
                                url = my_url) %>% as_tibble()

# Establish ERDDAP url and database name
JSATSinfo <- info('FED_JSATS_detects', url = my_url)

# Retrieve list of all studyIDs on FED_JSATS
studyid_list <- tabledap(JSATSinfo,
                         fields = c('study_id'),
                         url = my_url,
                         distinct = TRUE
) %>% 
  filter(study_id != "2017_BeaconTag") %>% 
  pull(study_id)

# Function that accepts a studyid and calls tabledap, allows repeated calls since
# You can not simply do an %in% query
get_erddap <- function(studyid) {
  df <- tabledap(JSATSinfo,
                 fields = c('fish_id','local_time','time','receiver_serial_number','receiver_river_km','study_id','receiver_location','receiver_general_location'),
                 paste0('study_id=', '"',studyid, '"'),
                 url = my_url) %>% 
    left_join(ReceiverDeployments %>% select(receiver_general_location, receiver_general_river_km)) %>% 
    left_join(TaggedFish %>% select(fish_id, release_river_km,release_location))
  
  # Rename columns and change column types as ERDDAP returns data all in character format
  df <- df %>% 
    rename(
      StudyID = study_id,
      FishID = fish_id,
      GEN = receiver_general_location,
      GenRKM = receiver_general_river_km,
      RelRKM = release_river_km,
      Rel_loc = release_location
    ) %>% 
    mutate(
      GenRKM = as.numeric(GenRKM),
      time = ymd_hms(time), # ERDDAP returns everything in character so parse time back into POSIXct
      RelRKM = as.numeric(RelRKM)
    )
  
  
  # For Winter aggregate: Caldwell Park_Rel, Bonnyview_Rel
  if (group == "Winter") {
    df <- df %>% 
    mutate(GEN = ifelse(GEN == "Caldwell Park_Rel", "Bonnyview_Rel", GEN),
           GenRKM = ifelse(GenRKM == 551.288, 540.2, GenRKM))
  }
  
  df
}

studyGroup <- studyid_list[str_detect(studyid_list, group)]

# Manually remove "ColemanLateFall_2020" which is online but not finished yet
studyGroup <- studyGroup[!studyGroup %in% c("DeerCk_SH_Wild_2018",
                                            "RBDD_WR_2018", "ColemanLateFall_2020",
                                            "Winter_H_2013_TE")]  

# Retreive ERDDAP detection data 
all_detections <- lapply(studyGroup, get_erddap)

## Figure out list of sites to use. 
## 1. For max # of studyID find all common GEN 
## 2. If one ends before the others, find common GEN among the remaining studyIDs
## 3. Continue until only 1 left
## 4. If last studyID does not make it to GG, use the remaining bomber sites

studyID_num <- length(studyGroup)
furthest_rkm <- 600

#all_detections_bind <- as_tibble(bind_rows(all_detections)) 
get_GEN <- function(detections) {
  detections %>% 
    mutate(GenRKM = ifelse(is.na(GenRKM), RelRKM, GenRKM)) %>% 
    select(StudyID, GEN, GenRKM) %>% 
    distinct()
}

all_detections_bind <- as_tibble(lapply(all_detections, get_GEN) %>% 
  bind_rows())

# Continue until there are no studyIDs remaining
while (studyID_num != 0) {
  # Get GEN, GenRKM for all detections for given number of studyID
  temp <- all_detections_bind %>% 
    select(StudyID, GEN, GenRKM) %>% 
    distinct() %>%
    filter(
      #GEN %in% reach.meta$GEN,
      GenRKM < furthest_rkm
    ) %>% 
    group_by(GEN) %>% 
    filter(n() == studyID_num) %>% 
    arrange(StudyID, desc(GenRKM))
  
  # Add these GEN to master list  
  # If first iteration create dataframe, else add on to it
  if (studyID_num == length(studyGroup)) {
    GEN_list <- temp %>% select(GEN, GenRKM) %>% distinct()
  }else {
    GEN_list <- rbind(GEN_list, temp %>% select(GEN, GenRKM) %>% distinct())
  }
  
  # Now identify the furthest GenRKM that was common among them 
  furthest_rkm <- min(GEN_list$GenRKM)
  studyID_num <- studyID_num - 1 
}

# If none of the StudyID's ever make it to GG use the remaining bombersites 
if (min(GEN_list$GenRKM) > 2) {
  GEN_list <- rbind(GEN_list, reach.meta %>%
                      filter(GenRKM < furthest_rkm) %>%
                      select(GEN, GenRKM))
}

GEN_list <- GEN_list %>% 
  left_join(
    ReceiverDeployments %>% 
      select(GEN = receiver_general_location, Region = receiver_region) %>% 
      distinct()
  ) %>% 
  filter(
    !(Region %in% c("East Delta", "North Delta", "South Delta", "West Delta", 
                    "Yolo Bypass")) | GEN %in% c("ChippsE", "ChippsW")
  )

release <- all_detections[[1]]$Rel_loc[1]

# Get raw number of fish detected at each bomber site
get_raw_detects <- function(detections) {
  detections %>% 
    # bind_rows() %>% 
    filter(
      GEN %in% c(release, GEN_list$GEN)
    ) %>% 
    select(
      StudyID, FishID, GEN, GenRKM
    ) %>% 
    distinct() %>% 
    group_by(StudyID, GEN, GenRKM) %>% 
    count() %>% 
    arrange(desc(GenRKM))
    
}

raw_survival <- lapply(all_detections, get_raw_detects) %>% 
  bind_rows()


### 3. Create EH from detections, using the bomber sites ----------------------------------------------------------------------

create_EH <- function(detections) {
  print(unique(detections$StudyID))
  
  
  # ### Aggregate Caldwell Park Rel with Bonnyview Rel for Winter group only
  # if (str_detect(i, "Winter")) {
  #   detections <- detections %>% 
  #     mutate(GEN = ifelse(GEN == "Caldwell Park_Rel", "Bonnyview_Rel", GEN),
  #            GenRKM = ifelse(GenRKM == 551.288, 540.2, GenRKM))
  # }
  
  
  # Get earliest detection for each fish at each GEN
  min_detects <- detections %>% 
    filter(GEN %in% c(GEN_list$GEN)) %>% 
    group_by(FishID, GEN, GenRKM) %>% 
    summarise(
      min_time = min(time)
    ) %>% 
    arrange(
      FishID, min_time
    )
  
  # Quick way to check if all the studyIDs have the same release site
  # lapply(all_detections, function(x) x %>% 
  #          filter(GenRKM == max(GenRKM)) %>% 
  #          pull(GEN) %>% 
  #          unique())
  
  # Find the release site = max GenRKM, so that we can add it into the bomber sites
  release <- unique(detections$Rel_loc)
  
  # Save the release GenRKM so that reach lengths can be calculated starting from release
  release_genrkm <<- unique(detections$RelRKM)
  
  EH <- expand.grid(
    sort(unique(detections$FishID)),
    GEN_list$GEN, stringsAsFactors = FALSE 
  )
  
  names(EH) <- c('FishID', 'GEN')  
  
  min_detects$detect <- 1
  
  # Add in detections 
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
  
  EH <- reshape(EH, idvar = 'FishID', timevar = 'GEN', direction = 'wide')
  colnames(EH) <- gsub('detect.', '', colnames(EH))
  
  EH[2] <- 1
  
  # Collapse reach encounters into a single column string, necessary for MARK
  EH.inp <- EH %>% 
    # nrow(GEN_info) + 1 is number of GEN locations plus 1 because column 1 is the FishID
    unite("ch", 2:(length(EH)), sep ="") %>% 
    mutate(studyID = unique(detections$StudyID))
  
  
}

# Create EH for all studyIDs in the study group
all_inp <- lapply(all_detections, create_EH)


### 4. Run MARK Phi.t p.t model  --------------------------------------------

get_cum_survival <- function(EH) {
  # Take the input data frame and the user-defined arguments and 
  # creates a list (processed data) containing the data and numerous 
  # defined attributes that the remaining functions use in defining the analysis models  
  all.process <- process.data(EH, model = "CJS", begin.time = 1)
  
  # Create the design data and PIM structure which depends on the selected type of
  # analysis model (e.g., CJS or Multistrata), number of occasions, grouping variables 
  # and other attributes of the data that were defined in the processed data
  all.ddl <- make.design.data(all.process)
  
  # remove previously saved lists
  rm(list=ls(pattern="Phi"))
  rm(list=ls(pattern="p.t"))
  rm(list=ls(pattern="p.dot"))
  
  ## Now, set up basic model structures
  # This is setting detection probability (p) to be a function of time x year, which is typically the best for simple survival analysis
  p.t <- list(formula= ~time) 
  
  # This is setting the survival model to be a function of time (which is actually reach) x year
  Phi.t <- list(formula= ~time)
  
  # Create a model list of the above models 
  cml = create.model.list("CJS")
  
  # Run mark.wrapper for all model structures (the combination of all phi and p models possible). Warning, this could take awhile if you have many models
  model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl, realvcv = TRUE)
  
  
  reaches <- nchar(EH$ch[1]) - 1
  
  # Single
  phi.t <- model.outputs$Phi.t.p.t$results$real$estimate[1:(reaches - 1)] # reaches - 1 because the last reach isn't meaningful
  phi.t.vcv <- model.outputs$Phi.t.p.t$results$real.vcv
  
  cum.phi <- cumprod(phi.t)
  
  # calculate standard errors for the cumulative product. Again, change the "17"s to the # reaches-1 that you have
  cum.phi.se <- deltamethod.special("cumprod", phi.t[1:reaches - 1], phi.t.vcv[1:(reaches - 1),1:(reaches - 1)])
  
  ### Output estimate, SE, LCI, UCI to a dataframe
  cumulative <- data.frame(cum.phi = cum.phi, cum.phi.se = cum.phi.se, #RKM = KM[1:reaches]
                           LCI = expit(logit(cum.phi)-1.96*sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))),
                           UCI = expit(logit(cum.phi)+1.96*sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))))
  
  # Round to 3 digits
  cumulative <- round(cumulative,3)
  
  # Add in reach start GEN, GenRKM, and StudyID
  cumulative <- cumulative %>% 
    mutate(
      GEN = GEN_list$GEN[2:(length(GEN_list$GEN)-1)], # Rel, and bomber sites minus the last one
      GenRKM = GEN_list$GenRKM[2:(length(GEN_list$GenRKM)-1)],
      StudyID = EH$studyID[1],
    )

  return(cumulative)
  
}

cumulative_survivals <- lapply(all_inp, get_cum_survival) %>% 
  bind_rows() 

cumulative_survivals <- cumulative_survivals %>% 
  left_join(
    raw_survival %>% 
      ungroup %>% 
      select(StudyID, GEN, n_fish = n), 
    by = c("StudyID", "GEN")
  )

# Build the file name and path, then write to csv
file_name <- paste0("./outputs/Cumulative Survival/", group, "_cum_survival.csv")
write_csv(cumulative_survivals, file_name)



# Clean up mark files
cleanup(ask=FALSE)

