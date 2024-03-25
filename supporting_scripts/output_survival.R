##%######################################################%##
#                                                          #
####      Survival Analysis for Shiny                   ####
#                                                          #
##%######################################################%##
# Author: Tom Pham
# Updated: Paul Carvalho (6/10/2022)
# Major overhaul by Cyril Michel 3/6/2023, major changes made:
# to automatically drop bad receiver sites 
# to automatically drop lower Sac River sites when weirs are overtopping
# to automatically rerun the mark model with profile intervals when parameter estimates have oversized CIs
# to automatically reduce the # of reaches to a max of 25 per study to make shiny plots more legible. drops lease useful recv sites preferentially (ones that create super short reaches)
# automatically create release groups and estimate per reach and cumulative survival separately per release group
# and if a release group has <5 fish released, drop the release group completely


# Script and data info: This script performs survival analysis using CJS model
# and outputs both reach per 10km and cumulative estimates to csv for use in
# the Shiny app
# Data source: Data extracted from NOAA ERDDAP data server (oceanview)


# ### ALWAYS RUN ----------------------------------------------------------
library(RMark)
library(tidyverse)
library(rerddap)
library(lubridate)
library(cowplot)
library(leaflet)
library(vroom)
library(data.table)
library(cder)

## this removes annoying messages when running summarise function
options(dplyr.summarise.inform = FALSE)


# memory.limit(size=56000)

# Set working directory to 'data' folder
setwd("C:/Users/Cyril.Michel/Desktop/AT_shiny/calfishtrack-shiny/data")
studyid_list <- fread("studyid_names.csv")

# Clear cached data in order to retrieve all latest data
cache_delete_all()

# Load TaggedFish and ReceiverDeployments tables through ERDDAP ---------------------------------------------------------------------

# Establish ERDDAP url and database name
my_url     <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo  <- info('FED_JSATS_taggedfish', url = my_url)
TaggedFish <- tabledap(JSATSinfo,
                       fields = c('fish_id','study_id','tag_id_hex','release_location','release_river_km','fish_release_date','release_latitude',
                                  'release_longitude'),
                       url = my_url) %>% as_tibble()

TaggedFish$fish_release_date_PST <- as.POSIXct(TaggedFish$fish_release_date, tz = "Etc/GMT+8", format = "%m/%d/%Y %H:%M:%S")
  
JSATSinfo           <- info('FED_JSATS_receivers', url = my_url)
ReceiverDeployments <- tabledap(JSATSinfo, 
                                # fields = c('receiver_serial_number','receiver_general_location','receiver_location','receiver_river_km',
                                #            'receiver_start','receiver_end','receiver_last_valid','receiver_general_river_km','receiver_region',
                                #            'receiver_general_latitude','receiver_general_longitude'),
                                url = my_url) %>% as_tibble()

## read in weir overtopping data
TIS <- cdec_query(stations = c("TIS"), sensors = 20, durations = "H", start.date = as.Date(min(TaggedFish$fish_release_date_PST)), end.date = as.Date(max(TaggedFish$fish_release_date_PST) + (7*60*60*24)))
#FRE_by_day <- cdec_query(stations = c("FRE"), sensors = 41, durations = "D", start.date = min(TaggedFish$fish_release_date_PST), end.date = max(TaggedFish$fish_release_date_PST) + (7*60*60*24))
#FRE_by_day$day <- as.Date(FRE_by_day$DateTime)

TIS$day <- as.Date(TIS$DateTime)
TIS_by_day <- aggregate(list(Value = TIS$Value), by = list(day= TIS$day), mean,na.rm = T)

# Load functions ----------------------------------------------------------
# Need to create logit() and expit() functions because they're not updated in clusterPower package
# naming expit "expiit" because for some unknown reason R (or Rstudio) drops the expit function when named "expit" after some time has passed in a session?!?
expiit <- function(x){
  out <- exp(x) / (1 + exp(x))
  return(out)
}
logit <- function(p){
  out <- log(p/(1-p))
  return(out)
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
                                                                  c("SacTrawl"),
                                                                  c("SacBlwGeo"),
                                                                  c("SacBlwSteamboat"),
                                                                  c("Sac_Rio_Vista"),
                                                                  c("Freeport"),
                                                                  c("GCID"),
                                                                  c("Altube")),
                                              replace_list = list(c("ChippsE", "ChippsW"),
                                                                  c("BeniciaE", "BeniciaW"),
                                                                  c("SacTrawl1", "SacTrawl2"),
                                                                  c("SacBlwGeo_1", "SacBlwGeo_2"),
                                                                  c("SacBlwSteamboat_1", "SacBlwSteamboat_2"),
                                                                  c("Sac_Rio_Vista_1", "Sac_Rio_Vista_2"),
                                                                  c("Freeport_1", "Freeport_2"),
                                                                  c("GCID_abv", "GCID_blw"),
                                                                  c("Abv_Altube1", "Abv_Altube2")))) {
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
  
  reach.meta <- get_receiver_GEN(detections)
  # Make a copy of reach.meta (receiver metadata)
  reach.meta.aggregate <- reach.meta
  
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
    reach.meta.aggregate <- reach.meta.aggregate %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
        Region = ifelse(GEN == "End", "End", ifelse(GEN %in% c(replace_list, replace_with), replace$Region, Region))
      ) %>% 
      distinct()
  }
  ## also, use mean GenRKM for all GENs, since some GEN recvs apear to have different GenRKMs
  gens <- aggregate(data = detections, GenRKM ~ GEN, mean, na.rm=T)
  detections$GenRKM <- NULL
  detections <- merge(detections, gens, by = "GEN")
  return(detections)
}

assign_releases <- function(detections_tmp, day_gap = 6, rkm_gap = 20, minimum_rel_size = 5) {
  # Make an encounter history df
  #
  # Arguments:
  #  detections: a detections df
  #     
  # Return:
  #  Encounter history df. A matrix of every fish tagged for a given studyID
  #  at every given receiver site (that is in reach.meta.aggregate) and whether
  #  it was present 1 or absent 0 in the detection df
  
  rels <- detections_tmp[substr(detections_tmp$GEN,nchar(detections_tmp$GEN)-3,nchar(detections_tmp$GEN))%in% c("_Rel","_rel"),]
  detections_tmp <- merge(detections_tmp,data.frame(FishID = rels$FishID, RelRKM = rels$GenRKM, fish_release_date = rels$local_time))
  releases <- aggregate(list(count = detections_tmp$FishID), by = list(fish_release_date = detections_tmp$fish_release_date, RelRKM = as.numeric(detections_tmp$RelRKM)), function(x){length(unique(x))})
  releases <- releases[order(releases$fish_release_date),]
  releases$days_since_last <- as.numeric(difftime(releases$fish_release_date, shift(releases$fish_release_date, type = "lag"), units = "days"))
  #releases$distance_since_last <- releases$RelRKM- shift(releases$RelRKM, type = "lag")
  
  releases$gap_time <- 0
  releases[which(releases$days_since_last >=day_gap), "gap_time"] <- 1
  #releases[which(releases$days_since_last >=day_gap | abs(releases$distance_since_last) > rkm_gap), "gap"] <- 1
  releases$group <- cumsum(releases$gap_time)+1
  
  releases <- releases[order(releases$group,releases$RelRKM),]
  releases$distance_since_last <- releases$RelRKM- shift(releases$RelRKM, type = "lag")
  
  releases$gap_space <- 0
  releases[which(abs(releases$distance_since_last) > rkm_gap), "gap_space"] <- 1
  releases$gap_final <- rowSums(releases[c("gap_space","gap_time")])
  releases[releases$gap_final >1,"gap_final"] <- 1
  
  #releases <- releases[order(releases$group,releases$RelRKM),]
  releases$group <- cumsum(releases$gap_final)+1
  #releases$group <- LETTERS[releases$group]
  
  detections_tmp <- merge(detections_tmp, releases)
  
  rel_size <- aggregate(data = detections_tmp, FishID ~ group, FUN = function(x){length(unique(x))})
  
  ## drop any fish release groups that had less fish than the minimum rel size
  detections_tmp <- detections_tmp[detections_tmp$group %in% rel_size[rel_size$FishID>minimum_rel_size,"group"],]
  ## give warning
  if(nrow(rel_size[rel_size$FishID<=minimum_rel_size,])>0){
    print(paste("dropping",nrow(rel_size[rel_size$FishID<=minimum_rel_size,]),"release groups out of",nrow(rel_size),"due to insufficient numbers"))
  }
  
  ## now reassign group so that there is a group 1, rmark doesn't seem to like groups starting at 2
  names <- unique(detections_tmp$group)
  if(length(names)>1){
    group_names <- data.frame(group = names[order(names)])
    group_names$group_new <- 1:nrow(group_names)
    detections_tmp <- merge(detections_tmp,group_names, by = "group")
    detections_tmp$group <- detections_tmp$group_new
    detections_tmp$group_new <- NULL
  }
  return(detections_tmp)
}


make_EH <- function(detections_tmp, rkm_gap = 20) {
  
  rels <- detections_tmp[substr(detections_tmp$GEN,nchar(detections_tmp$GEN)-3,nchar(detections_tmp$GEN)) %in% c("_Rel","_rel"),]
  rels_unique <- unique(rels[,c("group","GEN","GenRKM")])
  reach.meta.group <- get_receiver_GEN(detections_tmp)
  
  ## in the rare situation that there are release group(s) downstream of other release group(s), we have to remove the release groups
  ## as inp sites, and instead assign this release as 1s at the nearest upstream actual receiver site
  if(length(unique(detections_tmp$group)) > 1 & any(abs(diff(rels_unique$GenRKM))>rkm_gap)){
    downsteam_rels <- rels_unique[rels_unique$GenRKM < max(rels_unique$GenRKM),]
    ## find closest upstream receiver site
    for (j in 1:nrow(downsteam_rels)){
      reach.meta.upstream <- reach.meta.group[reach.meta.group$GenRKM >downsteam_rels[j,"GenRKM"],]
      detections_tmp[detections_tmp$GEN == downsteam_rels[j,"GEN"],"GenRKM"] <- reach.meta.upstream[which.min(reach.meta.upstream$GenRKM - downsteam_rels[j,"GenRKM"]),"GenRKM"]
      detections_tmp[detections_tmp$GEN == downsteam_rels[j,"GEN"],"GEN"] <- reach.meta.upstream[which.min(reach.meta.upstream$GenRKM - downsteam_rels[j,"GenRKM"]),"GEN"]
    }
  }
  
  # Get earliest detection for each fish at each GEN
  min_detects <- detections_tmp %>% 
    #filter(GEN %in% reach.meta.aggregate$GEN) %>% 
    group_by(FishID, GEN, GenRKM, group) %>% 
    summarise(min_time = min(local_time)) %>% 
    arrange(FishID, min_time) %>%
    ungroup()
  
  # Create matrix of all combinations of fish and GEN
  EH <- table(min_detects$FishID,min_detects$GEN) # A will be rows, B will be columns
  
  EH <- as.data.frame(EH)
  names(EH) <- c('FishID', 'GEN',"detect")  
  EH <- merge(EH, unique(min_detects[,c("FishID", "group")]))
  EH <- merge(EH, unique(min_detects[,c("GEN", "GenRKM")]))
  
  EH <- EH[order(EH$GenRKM,decreasing = T),]
  EH$GenRKM <- NULL
  
  # Reshape the df wide, so that columns are GEN locations, rows are fish, 
  # values are 1 or 0 for presence/absence
  EH <- reshape(EH, idvar = c('FishID','group'), timevar = 'GEN', direction = 'wide')
  colnames(EH) <- gsub('detect.', '', colnames(EH))
  # Manually make the release column a 1 because all fish were released there
  # sometimes detections df does not reflect that accurately
  #EH[2] <- 1
  return(EH)
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
    unite("ch", 3:(length(EH)), sep ="") %>% 
    # Use the detections df to get the StudyID assignment
    mutate(group = EH$group,studyID = unique(detections_tmp$studyID))
  EH.inp
}


get_mark_model <- function(all.inp, standardized = T, multi_model = "times", detections) {
  # Run a CJS Mark model
  #
  # Arguments:
  #  all.inp: inp df, can be more than one studyID
  #  standardized: TRUE or FALSE, if you want outputs to be standardized to 
  #     per10km or not
  #  multi_model: model type to use for cases of multiple StudyIDs. By default
  #     multi_model = times runs reach*StudyID, alternatively 
  #     multi_model = plus runs reach+StudyID
  #
  # Return:
  #  the outputs of running a CJS Mark model, df with phi and p estimates, LCI
  #  UCI, SE
  # Run a CJS Mark model for cumulative survival
  
  reach.meta.all <- get_receiver_GEN(detections)
  ## remove any downstream release locations before calculating reach_length
  reach.meta.all <- rbind(reach.meta.all[1,],
                            reach.meta.all[!(substr(reach.meta.all$GEN, nchar(reach.meta.all$GEN)-3,nchar(reach.meta.all$GEN))%in% c("_Rel","_rel")),])
  reach.meta.all$reach_num <- 0:(nrow(reach.meta.all)-1)
  
  reach_length <- round((abs(diff(reach.meta.all$GenRKM))/10), digits = 3) 

  if(length(unique(all.inp$group))==1){

    all.process <- process.data(all.inp, model="CJS", begin.time=1)

    all.ddl <- make.design.data(all.process)
    
    ## fix survival between benicia gates. too complicated to automate for now, just reducing Benicia down to 1 line
    # if(nrow(reach.meta.all[reach.meta.all$GEN %in% c("BeniciaE","BeniciaW"),])>1){
    #   all.ddl$Phi$fix <- NA
    #   beniciaw_reach <- reach.meta.all[reach.meta.all$GEN == "BeniciaW","reach_num"]
    #   all.ddl$Phi$fix[all.ddl$Phi$time == beniciaw_reach$reach_num] <- 1
    # }
    
    rm(list=ls(pattern="Phi.t"))
    rm(list=ls(pattern="p.t"))
    
    p.t <- list(formula= ~time) 
    Phi.t <- list(formula= ~time)
    
    cml = create.model.list("CJS")
    
    model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl,brief = T, output = F, silent = T)
    ## identify parameters estimates that have large CIs, rerun with profile.int
    model.outputs$Phi.t.p.t$results$real$CI <- model.outputs$Phi.t.p.t$results$real$ucl - model.outputs$Phi.t.p.t$results$real$lcl
    profile_int <- which(model.outputs$Phi.t.p.t$results$real$CI>0.75)
    if(length(profile_int)>0){
      print(paste("Using profile intervals on",length(profile_int),"parameters"))
      model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl,brief = T, output = F, silent = T, profile.int = profile_int)
    }
    
    # surv_index <- grep("Phi g",rownames(model.outputs$Phi.t.p.t$results$real))
    # p_index <- grep("p g",rownames(model.outputs$Phi.t.p.t$results$real))
    # ## for releases made further downstream than others, subset based in inp length
    # surv_index_length <- nrow(reach.meta.group)-1
    # surv_index <- tail(surv_index, n = surv_index_length)
    
    #phi.t <- model.outputs$Phi.t.p.t$results$real 
    surv_index <- grep("Phi g",rownames(model.outputs$Phi.t.p.t$results$real))
    unique_detects <- get_unique_detects(detections)
    # ## for releases made further downstream than others, subset based in inp length
    surv_index_length <- nrow(reach.meta.all)-1
    surv_index <- tail(surv_index, n = surv_index_length)
    phi.t <- model.outputs$Phi.t.p.t$results$real[surv_index,]
    if (standardized) {
      phi.t$se <- sqrt(((phi.t$estimate^((1-reach_length)/reach_length))/reach_length)^2 * phi.t$se^2)
      phi.t$estimate <- phi.t$estimate^(1/reach_length)
      phi.t$lcl <- expiit(logit(phi.t$estimate)-1.96*sqrt(phi.t$se^2/((exp(logit(phi.t$estimate))/(1+exp(logit(phi.t$estimate)))^2)^2)))
      phi.t$ucl <- expiit(logit(phi.t$estimate)+1.96*sqrt(phi.t$se^2/((exp(logit(phi.t$estimate))/(1+exp(logit(phi.t$estimate)))^2)^2)))
      # if se is 0, above code produces NaNs for ucl and lcl, but Rmark produces the estimate for these, so mimic
      phi.t[is.na(phi.t$lcl),c("lcl","ucl")] <- phi.t[is.na(phi.t$lcl),c("estimate")]
    }
    phi.t$reach_start <- reach.meta.all$GEN[1:nrow(reach.meta.all)-1]
    phi.t$reach_end <- reach.meta.all$GEN[2:nrow(reach.meta.all)]
    phi.t$rkm_start <- reach.meta.all$GenRKM[1:nrow(reach.meta.all)-1]
    phi.t$rkm_end <- reach.meta.all$GenRKM[2:nrow(reach.meta.all)]
    phi.t$Region <- reach.meta.all$Region[2:nrow(reach.meta.all)]
    phi.t$Reach <- paste(phi.t$reach_start,"to",phi.t$reach_end)
    phi.t$RKM <- paste(phi.t$rkm_start,"to",phi.t$rkm_end)
    phi.t$count_at_start <- unique_detects$count[1:nrow(unique_detects)-1]
    phi.t$count_at_end <- unique_detects$count[2:nrow(unique_detects)]
    phi.t$GenLat_start <- reach.meta.all$GenLat[1:nrow(reach.meta.all)-1]
    phi.t$GenLon_start <- reach.meta.all$GenLon[1:nrow(reach.meta.all)-1]
    phi.t$GenLat_end <- reach.meta.all$GenLat[2:nrow(reach.meta.all)]
    phi.t$GenLon_end <- reach.meta.all$GenLon[2:nrow(reach.meta.all)]
    phi.t$release <- reach.meta.all$GEN[1]
    phi.t$group <- 1
    phi.t$reach_num <- reach.meta.all$reach_num[2:nrow(reach.meta.all)]
    phi.t$mean_rel_date <- unique_detects$mean_rel_date[1]

    ### Output estimate, SE, LCI, UCI to a dataframe
    cumulative <- phi.t
    
  }else{
    
    ## we standardize after the fact when there are multiple release groups since reach lengths for downstream releases are wrong (since those detections get assigned to a rel loc nearby)
    all.inp$group <- as.factor(all.inp$group)
    all.process <- process.data(all.inp, model="CJS", begin.time=1,
                                  groups = "group")

    all.ddl <- make.design.data(all.process)
    
    # ## fix survival between benicia gates. too complicated to automate for now, just reducing Benicia down to 1 line
    # if(nrow(reach.meta.all[reach.meta.all$GEN %in% c("BeniciaE","BeniciaW"),])>1){
    #   all.ddl$Phi$fix <- NA
    #   beniciaw_reach <- reach.meta.all[reach.meta.all$GEN == "BeniciaW","reach_num"]
    #   all.ddl$Phi$fix[all.ddl$Phi$time == beniciaw_reach$reach_num] <- 1
    # }
    
    rm(list=ls(pattern="Phi.t"))
    rm(list=ls(pattern="p.t"))
    
    # Set the model up differently depending on single vs multiple StudyIDs
    # and multi_model argument (times vs plus)

    if (multi_model == "times") {
      p.t <- list(formula= ~time)
      Phi.txr <- list(formula= ~time*group)
    }else {
      p.t <- list(formula= ~time)
      Phi.txr <- list(formula= ~time+group)
    } 
      
    cml = create.model.list("CJS")
    
    model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl,brief = T, output = F, silent = T)
    #model.outputs$Phi.txr.p.t$results$real$reach <- as.numeric(sub('.*t', '', row.names(model.outputs$Phi.txr.p.t$results$real)))
    ## identify parameters estimates that have large CIs, rerun with profile.int
    model.outputs$Phi.txr.p.t$results$real$CI <- model.outputs$Phi.txr.p.t$results$real$ucl - model.outputs$Phi.txr.p.t$results$real$lcl
    profile_int <- which(model.outputs$Phi.txr.p.t$results$real$CI>0.9)
    if(length(profile_int)>0){
      print(paste("Using profile intervals on",length(profile_int),"parameters"))
      model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl,brief = T, output = F, silent = T, profile.int = profile_int)
    }
    
    
    ## create empty df to fill
    cumulative <- data.frame(model.outputs$Phi.txr.p.t$results$real[0,],
                             reach_start = as.character(),
                             reach_end = as.character(),
                             rkm_start = as.numeric(),
                             rkm_end = as.numeric(),
                             Region = as.character(),
                             Reach = as.character(),
                             RKM = as.character(),
                             count_at_start = as.numeric(),
                             count_at_end = as.numeric(),
                             GenLat_start = as.numeric(),
                             GenLon_start = as.numeric(),
                             GenLat_end = as.numeric(),
                             GenLon_end = as.numeric(),
                             release = as.character(), 
                             group = as.numeric())
    
    
    #reach.meta.all.w.rels <- get_receiver_GEN(detections)
    
    for (j in unique(all.inp$group)){
      reach.meta.group <- data.frame(get_receiver_GEN(detections[detections$group == j,]),group = j)
      reach.meta.group.all <- merge(reach.meta.group,reach.meta.all, all = T, by = "GEN")
      #reach.meta.group.all <- reach.meta.group.all[order(reach.meta.group.all$reach_num),]
      ## these next steps may seem clunky but come from a need of being able to get all the GEN metadata, but sometimes the data exists in reach.meta.group and not in reach.meta.all and vice versa
      reach.meta.group.all$GenRKM <- apply(reach.meta.group.all[,c("GenRKM.x","GenRKM.y")], MARGIN = 1, FUN = function(x){first(na.omit(x))})
      reach.meta.group.all$GenLat <- apply(reach.meta.group.all[,c("GenLat.x","GenLat.y")], MARGIN = 1, FUN = function(x){first(na.omit(x))})
      reach.meta.group.all$GenLon <- apply(reach.meta.group.all[,c("GenLon.x","GenLon.y")], MARGIN = 1, FUN = function(x){first(na.omit(x))})
      reach.meta.group.all$Region <- apply(reach.meta.group.all[,c("Region.x","Region.y")], MARGIN = 1, FUN = function(x){first(na.omit(x))})
      
      surv_index <- grep(paste("Phi g",j,sep =""),rownames(model.outputs$Phi.txr.p.t$results$real))
      unique_detects <- get_unique_detects(detections[detections$group == j,])
      unique_detects$count_at_start <- unique_detects$count
      unique_detects$count_at_end <- unique_detects$count
      # ## for releases made further downstream than others, subset based in inp length
      # surv_index_length <- nrow(reach.meta.group)-1
      # surv_index <- tail(surv_index, n = surv_index_length)
      reach.meta.group.all <- reach.meta.group.all[order(reach.meta.group.all$GenRKM,decreasing = T),]
      ## here, if release was downstream of other releases, we need to subset surv_index differently (due to having assigned release to existing recv site)
      if(min(which(is.na(reach.meta.group.all[,"group"])==F)) ==1){
        surv_index_index <- min(which(is.na(reach.meta.group.all[,"group"])==F)):(max(which(is.na(reach.meta.group.all[,"group"])==F))-1)
      }else{
        surv_index_index <- (min(which(is.na(reach.meta.group.all[,"group"])==F))-1):(max(which(is.na(reach.meta.group.all[,"group"])==F))-2)
      }
      #surv_index_index <- surv_index_index[surv_index_index >0]
      surv_index <- surv_index[surv_index_index]
      reach.meta.group.final <- reach.meta.group.all[min(which(is.na(reach.meta.group.all[,"group"])==F)):max(which(is.na(reach.meta.group.all[,"group"])==F)),]
      ## make sure release loc has a reach_num, this sometimes doesn't happen automatically when a release is downstream of another
      #reach.meta.group.final[1,"reach_num"] <- min(reach.meta.group.final$reach_num, na.rm = T)-1
      reach_length_group <- abs(diff(reach.meta.group.final$GenRKM))
      
      phi.t <- model.outputs$Phi.txr.p.t$results$real[surv_index,]
      if (standardized) {
          phi.t$se <- sqrt(((phi.t$estimate^((1-reach_length_group)/reach_length_group))/reach_length_group)^2 * phi.t$se^2)
          phi.t$estimate <- phi.t$estimate^(1/reach_length_group)
          phi.t$lcl <- expiit(logit(phi.t$estimate)-1.96*sqrt(phi.t$se^2/((exp(logit(phi.t$estimate))/(1+exp(logit(phi.t$estimate)))^2)^2)))
          phi.t$ucl <- expiit(logit(phi.t$estimate)+1.96*sqrt(phi.t$se^2/((exp(logit(phi.t$estimate))/(1+exp(logit(phi.t$estimate)))^2)^2)))
          # if se is 0, above code produces NaNs for ucl and lcl, but Rmark produces the estimate for these, so mimic
          phi.t[is.na(phi.t$lcl),c("lcl","ucl")] <- phi.t[is.na(phi.t$lcl),c("estimate")]
      }
      phi.t$reach_start <- reach.meta.group.final$GEN[1:nrow(reach.meta.group.final)-1]
      phi.t$reach_end <- reach.meta.group.final$GEN[2:nrow(reach.meta.group.final)]
      phi.t$rkm_start <- round(reach.meta.group.final$GenRKM[1:nrow(reach.meta.group.final)-1],3)
      phi.t$rkm_end <- round(reach.meta.group.final$GenRKM[2:nrow(reach.meta.group.final)],3)
      phi.t$Region <- reach.meta.group.final$Region[2:nrow(reach.meta.group.final)]
      phi.t$Reach <- paste(phi.t$reach_start,"to",phi.t$reach_end)
      phi.t$RKM <- paste(phi.t$rkm_start,"to",phi.t$rkm_end)
      phi.t$GenLat_start <- reach.meta.group.final$GenLat[1:nrow(reach.meta.group.final)-1]
      phi.t$GenLon_start <- reach.meta.group.final$GenLon[1:nrow(reach.meta.group.final)-1]
      phi.t$GenLat_end <- reach.meta.group.final$GenLat[2:nrow(reach.meta.group.final)]
      phi.t$GenLon_end <- reach.meta.group.final$GenLon[2:nrow(reach.meta.group.final)]
      phi.t$release <- reach.meta.group.final$GEN[1]
      phi.t$mean_rel_date <- unique_detects$mean_rel_date[1]
      
      phi.t <- merge(phi.t, unique_detects[,c("GEN","count_at_start")], by.x = "reach_start", by.y = "GEN", all.x = T)
      phi.t <- merge(phi.t, unique_detects[,c("GEN","count_at_end")], by.x = "reach_end", by.y = "GEN", all.x = T)
      phi.t[is.na(phi.t$count_at_start),"count_at_start"] <- 0
      phi.t[is.na(phi.t$count_at_end),"count_at_end"] <- 0
      
      # phi.t.vcv <- model.outputs$Phi.txr.p.t$results$real.vcv[surv_index,surv_index]
      # cum.phi <- cumprod(phi.t$estimate)
      # # calculate standard errors for the cumulative product. 
      # cum.phi.se <- deltamethod.special("cumprod", phi.t$estimate, phi.t.vcv)
      cumulative <- rbind(cumulative,data.frame(phi.t, group = j))
    }
    # p_index <- grep("p g1",rownames(model.outputs$Phi.txr.p.t$results$real))
    # surv_index_length <- nrow(reach.meta.all)-1
    # p_index <- tail(p_index, n = surv_index_length)
    # p.t <- model.outputs$Phi.txr.p.t$results$real[p_index,] 
    # cumulative <- rbind(cumulative,data.frame(p.t,GEN = reach.meta.all[2:nrow(reach.meta.all),"GEN"],group = "All"))
    cumulative <- merge(cumulative, reach.meta.all[,c("GEN","reach_num")], by.x = "reach_end", by.y = "GEN")
  }
  # Round to 3 digits
  cumulative[,c("estimate","se","lcl","ucl")] <- round(cumulative[,c("estimate","se","lcl","ucl")],3)
  cumulative$fixed <- NULL
  cumulative$note <- NULL
  
  ## in rare instances, profile likelihood still doesn't fix incorrectly large CIs, but if surv is 1 and SE is 0, we can assume ucl and lcl of 1
  cumulative[cumulative$estimate == 1 & cumulative$se < 0.03 & cumulative$lcl == 0,c("lcl","ucl")] <- 1
  
  cumulative <- cumulative[order(cumulative$group,cumulative$reach_num),]
  
  ## remove last reach due to unidentifiability
  cumulative <- cumulative[cumulative$reach_num < max(cumulative$reach_num),]
  return(cumulative)
  
}


remove_bad_sites <- function(all.inp, EH, p_cutoff = 0.7, detections) {
  
  reach.meta.all <- get_receiver_GEN(detections)
  ## remove any downstream release locations before calculating reach_length
  reach.meta.all <- rbind(reach.meta.all[1,],
                          reach.meta.all[!(substr(reach.meta.all$GEN, nchar(reach.meta.all$GEN)-3,nchar(reach.meta.all$GEN))%in% c("_Rel","_rel")),])
  reach.meta.all$reach_num <- 0:(nrow(reach.meta.all)-1)
  
  
  if(length(unique(all.inp$group))==1){
    
    all.process <- process.data(all.inp, model="CJS", begin.time=1)
    
    all.ddl <- make.design.data(all.process)
    
    rm(list=ls(pattern="Phi.t"))
    rm(list=ls(pattern="p.t"))
    
    p.t <- list(formula= ~time) 
    Phi.t <- list(formula= ~time)
    
    cml = create.model.list("CJS")
    
    model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl,brief = T, output = F, silent = T)
    
    #phi.t <- model.outputs$Phi.t.p.t$results$real 
    surv_index <- grep("p g",rownames(model.outputs$Phi.t.p.t$results$real))
    # ## for releases made further downstream than others, subset based in inp length
    surv_index_length <- nrow(reach.meta.all)-1
    surv_index <- tail(surv_index, n = surv_index_length)
    phi.t <- model.outputs$Phi.t.p.t$results$real[surv_index,]

    phi.t$GEN <- reach.meta.all$GEN[2:nrow(reach.meta.all)]
    
    
    ### Output estimate, SE, LCI, UCI to a dataframe
    cumulative <- phi.t
    
  }else{
    
    all.inp$group <- as.factor(all.inp$group)
    
    ## we standardize after the fact when there are multiple release groups since reach lengths for downstream releases are wrong (since those detections get assigned to a rel loc nearby)
    all.process <- process.data(all.inp, model="CJS", begin.time=1,
                                groups = "group")
    
    all.ddl <- make.design.data(all.process)
    
    rm(list=ls(pattern="Phi.t"))
    rm(list=ls(pattern="p.t"))
    
    # Set the model up differently depending on single vs multiple StudyIDs
    # and multi_model argument (times vs plus)
    # if lots of parameters, may need to keep survival just time variable to speed up
    if((ncol(EH)-3)*length(unique(EH$group))>100){
      p.t <- list(formula= ~time*group)
      Phi.txr <- list(formula= ~time)
    }else{
      p.t <- list(formula= ~time*group)
      Phi.txr <- list(formula= ~time*group) 
    }

    cml = create.model.list("CJS")
    
    model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl,brief = T, output = F, silent = T)
    #model.outputs$Phi.txr.p.t$results$real$reach <- as.numeric(sub('.*t', '', row.names(model.outputs$Phi.txr.p.t$results$real)))
    
    ## create empty df to fill
    cumulative <- model.outputs$Phi.txr.p.t$results$real[0,]
    
    
    #reach.meta.all.w.rels <- get_receiver_GEN(detections)
    
    for (j in unique(all.inp$group)){
      reach.meta.group <- data.frame(get_receiver_GEN(detections[detections$group == j,]),group = j)
      reach.meta.group.all <- merge(reach.meta.group,reach.meta.all, all = T, by = "GEN")
      #reach.meta.group.all <- reach.meta.group.all[order(reach.meta.group.all$reach_num),]
      ## these next steps may seem clunky but come from a need of being able to get all the GEN metadata, but sometimes the data exists in reach.meta.group and not in reach.meta.all and vice versa
      reach.meta.group.all$GenRKM <- apply(reach.meta.group.all[,c("GenRKM.x","GenRKM.y")], MARGIN = 1, FUN = function(x){first(na.omit(x))})
      reach.meta.group.all$GenLat <- apply(reach.meta.group.all[,c("GenLat.x","GenLat.y")], MARGIN = 1, FUN = function(x){first(na.omit(x))})
      reach.meta.group.all$GenLon <- apply(reach.meta.group.all[,c("GenLon.x","GenLon.y")], MARGIN = 1, FUN = function(x){first(na.omit(x))})
      reach.meta.group.all$Region <- apply(reach.meta.group.all[,c("Region.x","Region.y")], MARGIN = 1, FUN = function(x){first(na.omit(x))})
      
      surv_index <- grep(paste("p g",j,sep =""),rownames(model.outputs$Phi.txr.p.t$results$real))
      # ## for releases made further downstream than others, subset based in inp length
      # surv_index_length <- nrow(reach.meta.group)-1
      # surv_index <- tail(surv_index, n = surv_index_length)
      reach.meta.group.all <- reach.meta.group.all[order(reach.meta.group.all$GenRKM,decreasing = T),]
      ## here, if release was downstream of other releases, we need to subset surv_index differently (due to having assigned release to existing recv site)
      if(min(which(is.na(reach.meta.group.all[,"group"])==F)) ==1){
        surv_index_index <- min(which(is.na(reach.meta.group.all[,"group"])==F)):(max(which(is.na(reach.meta.group.all[,"group"])==F))-1)
      }else{
        surv_index_index <- (min(which(is.na(reach.meta.group.all[,"group"])==F))-1):(max(which(is.na(reach.meta.group.all[,"group"])==F))-2)
      }
      #surv_index_index <- surv_index_index[surv_index_index >0]
      surv_index <- surv_index[surv_index_index]
      reach.meta.group.final <- reach.meta.group.all[min(which(is.na(reach.meta.group.all[,"group"])==F)):max(which(is.na(reach.meta.group.all[,"group"])==F)),]
      
      phi.t <- model.outputs$Phi.txr.p.t$results$real[surv_index,]

      phi.t$GEN <- reach.meta.group.final$GEN[2:nrow(reach.meta.group.final)]
      phi.t$GEN <- reach.meta.group.final$GEN[2:nrow(reach.meta.group.final)]
      
      # phi.t.vcv <- model.outputs$Phi.txr.p.t$results$real.vcv[surv_index,surv_index]
      # cum.phi <- cumprod(phi.t$estimate)
      # # calculate standard errors for the cumulative product. 
      # cum.phi.se <- deltamethod.special("cumprod", phi.t$estimate, phi.t.vcv)
      cumulative <- rbind(cumulative,data.frame(phi.t, group = j))
    }

  }
  
  ## here, select sites that have a p < p_cutoff (default 0.7)
  bad_sites <- unique(cumulative[cumulative$estimate < p_cutoff,"GEN"])
  ## don't remove rel sites or GG sites
  bad_sites <- bad_sites[!(substr(bad_sites,nchar(bad_sites)-3,nchar(bad_sites))%in% c("_Rel","_rel"))]
  bad_sites <- bad_sites[substr(bad_sites,1,6)!= "Golden"]
  if(length(bad_sites)>0){
    print(paste("dropped", bad_sites))
    detections <- detections[!detections$GEN %in% bad_sites,]
  }
  return(detections)

}  
 
get_cum_survival <- function(all.inp, add_release, detections) {
  # Run a CJS Mark model for cumulative survival
  #
  # Arguments:
  #  all.inp: inp df, can be more than one studyID
  #  add_release: TRUE or FALSE, if you wish to add an extra dummy row at the top
  #  to show 100% survival at the release location
  #
  # Return:
  #  Cumulative survival outputs of CJS Mark model
  
  reach.meta.all <- get_receiver_GEN(detections)
  ## remove any downstream release locations before calculating reach_length
  reach.meta.all <- rbind(reach.meta.all[1,],
                          reach.meta.all[!(substr(reach.meta.all$GEN, nchar(reach.meta.all$GEN)-3,nchar(reach.meta.all$GEN))%in% c("_Rel","_rel")),])
  reach.meta.all$reach_num <- 0:(nrow(reach.meta.all)-1)
  
  if(length(unique(all.inp$group))==1){
    
    reach.meta.group <- get_receiver_GEN(detections)
    
    all.process <- process.data(all.inp, model = "CJS", begin.time = 1)
    all.ddl <- make.design.data(all.process)
    
    rm(list=ls(pattern="Phi.t"))
    rm(list=ls(pattern="p.t"))
    
    p.t <- list(formula= ~time) 
    Phi.t <- list(formula= ~time)
    
    cml = create.model.list("CJS")
    
    model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl, realvcv = TRUE,brief = T, output = F, silent = T)
    ## identify parameters estimates that have large CIs, rerun with profile.int
    model.outputs$Phi.t.p.t$results$real$CI <- model.outputs$Phi.t.p.t$results$real$ucl - model.outputs$Phi.t.p.t$results$real$lcl
    profile_int <- which(model.outputs$Phi.t.p.t$results$real$CI>0.9)
    if(length(profile_int)>0){
      print(paste("Using profile intervals on",length(profile_int),"parameters"))
      model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl,brief = T, realvcv = TRUE, output = F, silent = T, profile.int = profile_int)
    }
    
    model.outputs$Phi.t.p.t$results$real$reach <- as.numeric(sub('.*t', '', row.names(model.outputs$Phi.t.p.t$results$real)))
    
    unique_detects <- get_unique_detects(detections)
    surv_index <- grep("Phi g",rownames(model.outputs$Phi.t.p.t$results$real))
    ## for releases made further downstream than others, subset based in inp length
    surv_index_length <- nrow(reach.meta.group)-1
    surv_index <- tail(surv_index, n = surv_index_length)
    
    phi.t <- model.outputs$Phi.t.p.t$results$real[surv_index,] 
    phi.t.vcv <- model.outputs$Phi.t.p.t$results$real.vcv[surv_index,surv_index]
    
    cum.phi <- cumprod(phi.t$estimate)
    
    # calculate standard errors for the cumulative product. 
    cum.phi.se <- deltamethod.special("cumprod", phi.t$estimate, phi.t.vcv)
    
    
    ### Output estimate, SE, LCI, UCI to a dataframe

    cumulative <- rbind(data.frame(reach = min(phi.t$reach)-0.5,
                                   cum.phi = 1, 
                                   cum.phi.se = NA,
                                   LCI = NA,
                                   UCI = NA,
                                   GEN = reach.meta.group[1,"GEN"],
                                   group = 1,      
                                   release = reach.meta.group$GEN[1],
                                   mean_rel_date = unique_detects$mean_rel_date[1]),
                        data.frame(reach = phi.t$reach,
                                   cum.phi = cum.phi, 
                                   cum.phi.se = cum.phi.se,
                                   LCI = expiit(logit(cum.phi) - 1.96 * sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))),
                                   UCI = expiit(logit(cum.phi) + 1.96 * sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))),
                                   GEN = reach.meta.group[2:nrow(reach.meta.group),"GEN"],
                                   group = 1,      
                                   release = reach.meta.group$GEN[1],
                                   mean_rel_date = unique_detects$mean_rel_date[1]))
    
    
    # If add_release TRUE, add in the dummy row to the top which just represents
    # survival of 100% at release
    if (add_release == T) {
      cumulative <- cumulative %>% 
        add_row(
          .before = 1,
          group = 1,
          reach = 0.5,
          cum.phi = 1,
          cum.phi.se = NA,
          LCI = NA, 
          UCI = NA
        ) %>% 
        mutate(
          studyID = all.inp$studyID[1]
        )
    }else {
      cumulative <- cumulative %>% 
        mutate(
          studyID = all.inp$studyID[1]
        )
    }
  }else{
    
    all.inp$group <- as.factor(all.inp$group)
    all.process <- process.data(all.inp, model = "CJS", begin.time = 1, groups = "group")
    all.ddl <- make.design.data(all.process)
    
    rm(list=ls(pattern="Phi.t"))
    rm(list=ls(pattern="p.t"))
    
    p.t <- list(formula= ~time) 
    Phi.txr <- list(formula= ~time*group)
    
    cml = create.model.list("CJS")
    
    model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl, realvcv = TRUE,brief = T, output = F, silent = T)
    ## identify parameters estimates that have large CIs, rerun with profile.int
    model.outputs$Phi.txr.p.t$results$real$CI <- model.outputs$Phi.txr.p.t$results$real$ucl - model.outputs$Phi.txr.p.t$results$real$lcl
    profile_int <- which(model.outputs$Phi.txr.p.t$results$real$CI>0.9)
    if(length(profile_int)>0){
      print(paste("Using profile intervals on",length(profile_int),"parameters"))
      model.outputs <- mark.wrapper(cml, data=all.process, ddl=all.ddl,brief = T, realvcv = TRUE, output = F, silent = T, profile.int = profile_int)
    }
    
    model.outputs$Phi.txr.p.t$results$real$reach <- as.numeric(sub('.*t', '', row.names(model.outputs$Phi.txr.p.t$results$real)))
    
    
    ### Output estimate, SE, LCI, UCI to a dataframe
    cumulative <- data.frame(reach = as.numeric(),
                             cum.phi = as.numeric(), 
                             cum.phi.se = as.numeric(), 
                             # LCI = cum.phi - 1.96 * cum.phi.se,
                             # UCI = cum.phi + 1.96 * cum.phi.se)
                             LCI = as.numeric(), 
                             UCI = as.numeric(),
                             GEN = as.character(),
                             group = as.character(),
                             release = as.character(),
                             mean_rel_date = as.character())
    
    for (j in unique(all.inp$group)){

      reach.meta.group <- data.frame(get_receiver_GEN(detections[detections$group == j,]),group = j)
      reach.meta.group.all <- merge(reach.meta.group,reach.meta.all, by = "GEN", all = T)
      reach.meta.group.all <- reach.meta.group.all[order(reach.meta.group.all$reach_num),]
      reach.meta.group.all$GenRKM <- apply(reach.meta.group.all[,c("GenRKM.x","GenRKM.y")], MARGIN = 1, FUN = function(x){first(na.omit(x))})
      
      unique_detects <- get_unique_detects(detections[detections$group == j,])
      
      surv_index <- grep(paste("Phi g",j,sep =""),rownames(model.outputs$Phi.txr.p.t$results$real))
      #reach.meta.group.all <- reach.meta.group.all[order(reach.meta.group.all$GenRKM,decreasing = T),]
      reach.meta.group.all <- reach.meta.group.all[order(reach.meta.group.all$GenRKM,decreasing = T),]
      
      ## here, if release was downstream of other releases, we need to subset surv_index differently (due to having assigned release to existing recv site)
      if(min(which(is.na(reach.meta.group.all[,"group"])==F)) ==1){
        surv_index_index <- min(which(is.na(reach.meta.group.all[,"group"])==F)):(max(which(is.na(reach.meta.group.all[,"group"])==F))-1)
      }else{
        surv_index_index <- (min(which(is.na(reach.meta.group.all[,"group"])==F))-1):(max(which(is.na(reach.meta.group.all[,"group"])==F))-2)
      }
      #surv_index_index <- surv_index_index[surv_index_index >0]
      surv_index <- surv_index[surv_index_index]
      reach.meta.group.final <- reach.meta.group.all[min(which(is.na(reach.meta.group.all[,"group"])==F)):max(which(is.na(reach.meta.group.all[,"group"])==F)),]
      phi.t <- model.outputs$Phi.txr.p.t$results$real[surv_index,] 
      phi.t.vcv <- model.outputs$Phi.txr.p.t$results$real.vcv[surv_index,surv_index]
      cum.phi <- cumprod(phi.t$estimate)
      # calculate standard errors for the cumulative product. 
      if(length(surv_index)==1){
        cum.phi.se <- phi.t$se
      }else{
        cum.phi.se <- deltamethod.special("cumprod", phi.t$estimate, phi.t.vcv)
      }
      cumulative <- rbind(cumulative,rbind(data.frame(reach = min(phi.t$reach)-0.5,
                                                cum.phi = 1, 
                                                cum.phi.se = NA,
                                                LCI = NA,
                                                UCI = NA,
                                                GEN = reach.meta.group[1,"GEN"],
                                                group = j,
                                                release = reach.meta.group$GEN[1],
                                                mean_rel_date = unique_detects$mean_rel_date[1]),
                                          data.frame(reach = phi.t$reach,
                                                cum.phi = cum.phi, 
                                                cum.phi.se = cum.phi.se,
                                                LCI = expiit(logit(cum.phi) - 1.96 * sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))),
                                                UCI = expiit(logit(cum.phi) + 1.96 * sqrt(cum.phi.se^2/((exp(logit(cum.phi))/(1+exp(logit(cum.phi)))^2)^2))),
                                                GEN = reach.meta.group.final[2:nrow(reach.meta.group.final),"GEN"],
                                                group = j,      
                                                release = reach.meta.group$GEN[1],
                                                mean_rel_date = unique_detects$mean_rel_date[1])))
    }
  }
  # Round to 3 digits
  cumulative[,c("cum.phi", "cum.phi.se", "LCI", "UCI")] <- round(cumulative[,c("cum.phi", "cum.phi.se", "LCI", "UCI")],3)
  return(cumulative)
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
    select(FishID, GEN, GenRKM,group,fish_release_date) %>% 
    distinct() %>% 
    group_by(group, GEN, GenRKM) %>% 
    summarise(
      count = n(),
      mean_rel_date = mean(as.Date(fish_release_date))
    ) %>% 
    arrange(desc(GenRKM)) %>% 
    ungroup()
}


format_cum_surv <- function(cum_survival_all, detections) {
  # Format cumulative survival outputs for plotting and table outputs
  #
  # Arguments:
  #  cum_survival_all: output df from Mark model cumulative survival
  #
  # Return:
  #  properly formatted df for phi outputs, now ready to plot
  
  reach.meta.aggregate <- get_receiver_GEN(detections)
  
  
  cum_survival_all <- merge(cum_survival_all,reach.meta.aggregate, all.x = T)
  
  
  cum_survival_all <- cum_survival_all %>% 
    # add_column(
    #   GEN = rep(reach.meta.aggregate$GEN[1:(length(reach.meta.aggregate$GEN))], 
    #             length(groups)),
    #   RKM = rep(reach.meta.aggregate$GenRKM[1:(length(reach.meta.aggregate$GenRKM))], 
    #             length(groups)),
    #   reach_num = rep(seq(0, (nrow(reach.meta.aggregate))-1, 1), length(groups))
    # ) %>% 
    # left_join(
    #   reach.meta.aggregate %>%
    #     select(GEN, Region) %>% 
    #     distinct(),
    #   by = c("GEN")
    # ) %>% 
    mutate(
      cum.phi = round(cum.phi, digits = 3),
      cum.phi.se = round(cum.phi.se, digits = 3),
      LCI = round(LCI, digits = 2),
      UCI = round(UCI, digits = 2),
      'Survival estimate (SE)' = paste0(cum.phi, " (", as.character(cum.phi.se), ")"),
      reach_num = reach,
    )
  
  ## remove last reach due to unidentifiability
  cum_survival_all <- cum_survival_all[cum_survival_all$reach_num < max(cum_survival_all$reach_num),]
  
  return(cum_survival_all)
}


#### Get Reach per 10km Survival estimates ####------------------------------------------

# 1. Run all studyids

studyid_to_run <- as.data.frame(studyid_list[studyid_list$incl_survival == "YES",])

## here, uncomment code to instead run just a few select studies
#studyid_to_run <- studyid_to_run[studyid_to_run$study_id %in% c("Wild_stock_Chinook_Rbdd_2021","Wild_stock_Chinook_RBDD_2022","Lower_Yuba_FRH_Chinook_2023" ),]

for(i in 1:nrow(studyid_to_run)){
  print(studyid_to_run[i,"study_id"])
  day_gap_custom <- 6
  rkm_gap_custom <- 20
  p_cutoff_custom <- 0.7
  if(is.na(studyid_to_run[i,"custom_rkm_gap"])==F){
    rkm_gap_custom <- studyid_to_run[i,"custom_rkm_gap"]
    print(paste("Custom rkm gap",rkm_gap_custom))
  }
  if(is.na(studyid_to_run[i,"custom_day_gap"])==F){
    day_gap_custom <- studyid_to_run[i,"custom_day_gap"]
    print(paste("Custom day gap",day_gap_custom))
  }
  if(is.na(studyid_to_run[i,"p_cutoff_custom"])==F){
    p_cutoff_custom <- studyid_to_run[i,"p_cutoff_custom"]
    print(paste("Custom p cutoff for bad sites",p_cutoff_custom))
  }
  # 2. Read the detections and get the receiver GEN data; if detections file
  # has not been saved for the studyID yet, run save_new_detections.R first
  try({  
    detect_file <- paste0('detections/', studyid_to_run[i,"study_id"], ".csv")
    detections <- vroom(detect_file, show_col_types=FALSE)
    #detections <- get_detections(studyID)
    
    ## Head_of_Old_River is a release site but doesn't contain "_Rel"
    detections[detections$GEN == "Head_of_Old_River","GEN"] <- "Head_of_Old_River_Rel"
    
    ## Altube Island is a release site but doesn't contain "_Rel"
    detections[detections$GEN == "Altube Island","GEN"] <- "Altube Island_Rel"
    
    # 3. Remove Delta receivers, since we never use these due to braiding
    #detections <- detections[detections$Region %in% c("Carquinez Strait","Battle Ck","Lower Sac R","SF Bay","Upper Sac R") | detections$GEN %in% c("ChippsE",	"ChippsW"),]
    aggregated <- detections[!(detections$Region %in% c("South Delta","East Delta","North Delta","West Delta","Delta","Yolo Bypass","Sutter Bypass","Suisun Bay")),]
    aggregated <- aggregated[!(aggregated$Region %in% c("Lower Mok R") & aggregated$GenRKM< 139),]
    
    ## add back in chipps since this is downstream of braiding
    aggregated <- rbind(aggregated, detections[detections$GEN %in% c("ChippsE",	"ChippsW"),])
    ## add back in any release locs in the delta
    aggregated <- rbind(aggregated,detections[substr(detections$GEN,nchar(detections$GEN)-3,nchar(detections$GEN)) %in% c("_Rel","_rel") &
                                                detections$Region %in% c("South Delta","East Delta","North Delta","West Delta","Delta","Yolo Bypass","Sutter Bypass"),])
    
    ## remove lower river receiver sites if bypasses were flooding
    min_detects <- aggregated %>%
      #filter(GEN %in% reach.meta.aggregate$GEN) %>%
      group_by(Region, FishID, GEN, GenRKM) %>%
      summarise(min_time = min(local_time)) %>%
      arrange(GenRKM, min_time) %>%
      ungroup() %>%
      group_by(Region, GEN, GenRKM) %>%
      summarise(min_arrival = min(min_time),max_arrival = max(min_time)) %>%
      arrange(GenRKM, min_arrival) %>%
      ungroup()
    river_min_detects <- min_detects[min_detects$Region == "Lower Sac R",]
    if(nrow(river_min_detects)>0){
    ## using Tisdale weir because it is first to spill so most conservative
      if(sum(TIS_by_day[TIS_by_day$day >= min(river_min_detects$min_arrival) & TIS_by_day$day <= max(river_min_detects$max_arrival),"Value"],na.rm = T)>100){
        print("dropping lower Sac locations due to weir spilling")
        ## RKM 331 is just upstream of Moulton Weir, so most conservative since Moulton is last to spill
        aggregated <- aggregated[!(aggregated$Region == "Lower Sac R" & aggregated$GenRKM < 331),]
      }
    }
    # 3. Using the replacement dictionary, aggregate the detections
    aggregated <- aggregate_GEN(aggregated)
    
    # 4. assign release groups automatically
    aggregated <- assign_releases(aggregated, day_gap = day_gap_custom, rkm_gap = rkm_gap_custom)
    
    # 5. Remove any fish detections upstream of release location
    aggregated <- aggregated[aggregated$RelRKM >= aggregated$GenRKM,]
    
    # Similarly, remove any fish that were released in Feather and run up Sac (or Moke), or vice versa
    rel_river <- unique(aggregated[substr(aggregated$GEN,nchar(aggregated$GEN)-3,nchar(aggregated$GEN)) %in% c("_Rel","_rel"),"Region"])
    if(any(rel_river == "Feather_R")){
      aggregated <- aggregated[!(aggregated$GenRKM > 205 & aggregated$Region == "Lower Sac R"),]
      aggregated <- aggregated[!(aggregated$Region == "Lower Mok R"),]
    }
    if(any(rel_river %in% c("Upper Sac R","Battle Ck"))){
      aggregated <- aggregated[!(aggregated$GenRKM > 205 & aggregated$Region == "Feather_R"),]
      aggregated <- aggregated[!(aggregated$Region == "Lower Mok R"),]
    }
    if(any(rel_river %in% c("Lower Mok R"))){
      aggregated <- aggregated[!(aggregated$Region == "Lower Sac R"),]
    }
    if(any(rel_river %in% c("South Delta","San Joaquin River"))){
      aggregated <- aggregated[!(aggregated$Region == "Lower Sac R"),]
    }
    
    ## Also, remove any "_Rel_Rec" receivers. these are meant to double check that all tags are on at release, but gum up the works of this code
    aggregated <- aggregated[substr(aggregated$GEN, nchar(aggregated$GEN)-7,nchar(aggregated$GEN)) != "_Rel_Rec",]
    ## these next 2 are also a release receiver but not labeled as such
    aggregated <- aggregated[aggregated$GEN != "UpperButte_RST",]
    aggregated <- aggregated[aggregated$GEN != "FR_Gridley_Rel_DS",]
    aggregated <- aggregated[aggregated$GEN != "BoydsPump",]
    
    ## finally, a back door to remove pervasively erroneous sites
    if (studyid_to_run[i,"study_id"] == "SacRiverSpringJPE_2023"){
      aggregated <- aggregated[aggregated$GEN != "Blw_Salt",]
      aggregated <- aggregated[aggregated$GEN != "BlwOrd",]
    }
    if (studyid_to_run[i,"study_id"] == "ColemanLateFall_2019"){
      aggregated <- aggregated[aggregated$GEN != "Abv_Otter_Island",]
    }
    
    # 5. Create an encounter history from the aggregated detections
    EH <- make_EH(aggregated, rkm_gap = rkm_gap_custom) #summarise by FishID, GEN
    
    # 6. Create an inp from the given encounter history
    inp <- create_inp(detections_tmp = aggregated, EH = EH)
    
    # 7. Remove bad sites
    aggregated <- remove_bad_sites(all.inp = inp, EH = EH, detections = aggregated, p_cutoff = p_cutoff_custom)
    
    reach.meta <- get_receiver_GEN(aggregated)
    
    # QA the detections for 
    # -duplicate GenRKM with different names
    # -GEN that is above the release RKM
    
    
    # Check for any duplicate GenRKM
    dupe<- reach.meta %>% 
      group_by(GenRKM) %>% 
      filter(n() > 1)
    # Three Mile Slough is a dupe of Three Mile North
    # Freeport doesn't show up as a dupe because RKM dist is different, but it looks like one point on the map - will aggregate
    
    # Visually inspect receiver locations, determine if sites need to be removed
    leaflet(data = reach.meta) %>% 
      addTiles() %>% 
      addMarkers(lng = ~GenLon, lat = ~GenLat, label = ~GEN, 
                 labelOptions = labelOptions(noHide = T))
    
    ## remove detections from sites that were removed manually above
    aggregated <- aggregated[aggregated$GEN %in% reach.meta$GEN,]
    
    ## if study has a lot of reaches (>25), drop some short reaches
    reach.meta.all <- as.data.frame(get_receiver_GEN(aggregated))
    if(length(unique(reach.meta.all$GEN))>25){
      ## make sure not to remove the receiver site just upstream or downstream of rel locs
      reach.meta.all$keep <- 0
      reach.meta.all$reach_num <- 0:(nrow(reach.meta.all)-1)
      rels <- reach.meta.all[(substr(reach.meta.all$GEN, nchar(reach.meta.all$GEN)-3,nchar(reach.meta.all$GEN))%in% c("_Rel","_rel")),"reach_num"]
      reach.meta.all[reach.meta.all$reach_num %in% (rels[rels>0]-1),"keep"] <- 1
      reach.meta.all[reach.meta.all$reach_num %in% (rels+1),"keep"] <- 1
      ## make sure to also not drop GG_W or Benicia w
      reach.meta.all[reach.meta.all$GEN %in% c("GoldenGateW","BeniciaW"),"keep"] <- 1
      ## remove any release locations before calculating reach_length
      reach.meta.all <- reach.meta.all[!(substr(reach.meta.all$GEN, nchar(reach.meta.all$GEN)-3,nchar(reach.meta.all$GEN))%in% c("_Rel","_rel")),]
      reach.meta.all$reach_length <- c(NA,round((abs(diff(reach.meta.all$GenRKM))/10), digits = 3))
      reach.meta.all <- reach.meta.all[order(reach.meta.all$reach_length, decreasing = F),]
      num_to_drop <- length(unique(reach.meta.all$GEN))-25
      print(paste("dropping",num_to_drop,"sites to shorten to 25 reaches"))
      reach.meta.all <- reach.meta.all[which(reach.meta.all$keep == 0),]
      
      aggregated <- aggregated[!(aggregated$GEN %in% reach.meta.all[1:num_to_drop,"GEN"]),]
      reach.meta.all <- as.data.frame(get_receiver_GEN(aggregated))
      leaflet(data = reach.meta.all) %>% 
        addTiles() %>% 
        addMarkers(lng = ~GenLon, lat = ~GenLat, label = ~GEN, 
                   labelOptions = labelOptions(noHide = T))
    }
    

    # 8. Remake EH and inp
    EH <- make_EH(aggregated, rkm_gap = rkm_gap_custom)
    inp <- create_inp(detections_tmp = aggregated, EH = EH)
    
    # 7. Run reach per 10km survival in Mark using the inp
    
    # Get the estimates
    reach_surv <- get_mark_model(inp, standardized = T, detections = aggregated)
    cleanup(ask = F)
    
    # Format the estimates table
    #reach_surv_formatted <- format_phi(outputs = reach_surv, multiple = F, detections = aggregated)
    reach_surv$release_site <- reach_surv$release
    reach_surv$release <- paste(reach_surv$mean_rel_date,reach_surv$release_site, sep = "_")
    
    write_csv(reach_surv, paste0("Survival/Reach Survival per 10km/", studyid_to_run[i,"study_id"], "_reach_survival.csv"))
    
    # Get Cumulative Survival Estimates ---------------------------------------
    
    # 1. Get the cumulative survival estimates
    cum_surv <- get_cum_survival(inp, add_release = F, detections = aggregated)
    #cum_surv <- get_cum_survival(detections)
    
    # 2. Format the estimates table
    cum_surv_formatted <- format_cum_surv(cum_surv, detections = aggregated)
    
    # Add in unique fish counts at each GEN and lat/lon
    fish_count <- get_unique_detects(aggregated) #Summarise by StudyID, GEN
    
    cum_surv_formatted$group <- as.integer(cum_surv_formatted$group)
    
    cum_surv_formatted <- cum_surv_formatted %>% 
      left_join(
        fish_count %>% 
          select(GEN,group, count) %>% 
          distinct(),by = join_by(GEN, group)
      ) %>% 
      mutate_if(is.numeric, coalesce, 0) %>% 
      left_join(
        reach.meta %>% 
          select(GEN, GenLat, GenLon),
        by = join_by(GEN, GenLat, GenLon)
      )
    
    cum_surv_formatted <- cum_surv_formatted[order(cum_surv_formatted$group,cum_surv_formatted$reach),]
    
    cum_surv_formatted$release_site <- cum_surv_formatted$release
    cum_surv_formatted$release <- paste(cum_surv_formatted$mean_rel_date,cum_surv_formatted$release_site, sep = "_")
    
    write_csv(cum_surv_formatted, paste0("Survival/Cumulative Survival/", studyid_to_run[i,"study_id"], "_cum_survival.csv"))
    
    cleanup(ask = F)
  })
}

