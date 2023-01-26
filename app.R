# Author: Tom Pham
#
# Updated: 12/08/2022 by Paul Carvalho (paul.carvalho@noaa.gov)
#
# Delta route-specific survival developed by: Jake Kelley (jrkelley@usgs.gov) and Russ Perry (rperry@usgs.gov)
#
# Purpose: Shiny app for UCSC/NOAA Acoustic Telemetry Project
# This app visualizes data for the acoustic project in several ways
# -Background project data
# -Receiver Deployment map/table
# -Sacramento River hydrology graph at several CDEC stations
# -Outmigration animation of downstream movement of fish
# -Time of day in which fish are detected
# -Survival (reach per 10km, cumulative, and delta routing/survival (provided by USGS)) plots/tables
# -Movement: travel time from rel to GEN

# Load data
# --- ReceiverDeployments
# --- VemcoReceiverDeployments
# --- StudyIDs
# --- TaggedFish
# --- CDEC Flows

library(shinydashboard) # shiny tabs
library(shinythemes)
library(leaflet) # web mapping
library(leaflet.extras)
library(leaflet.minicharts) # web mapping animation
library(dygraphs) # interactive time series graphing
library(plotly) # interactive plotting
library(cder)
# May need to update CDEC package to library(cder)
library(xts) # formatting for dygraphs
library(lubridate) # date time madness
library(tidyverse)
library(suncalc) # getting sunset sunrise times
library(arsenal) # nice summary stats tables
library(gt) # fancy tables
library(rerddap) # ERDDAP data retrievals
library(httr) # Check HTTP status for CDEC/ERDDAP
library(vroom) # Fastest way to read csv
library(sf) # To display gis layers
library(DT)
library(DiagrammeR) # for Delta route survival diagram
library(kableExtra) # for Delta route survival table


# Global ------------------------------------------------------------------------------------------------------------------------------

## Load Receiver Deployments -----------------------------------------------
# Clear cached data in order to retrieve all latest data
cache_delete_all()

# If last dl was > 90 days update files, else read in csv

# Download updates every 90 days
last_checked_date <- read_rds("./data/last_checked_date.RDS")

# If last update check was < 90 days read in CSVs (much faster load times)
if (as.numeric(Sys.Date() - last_checked_date) < 90) {
   ReceiverDeployments <- vroom("./data/rec_deployments.csv", col_types = cols())
   
} else { 
  # Else check if ERDDAP is online, x returns TRUE if database is down or "Timeout"
  # if the http check timeouts out 
  x <- tryCatch(http_error("oceanview.pfeg.noaa.gov/erddap/tabledap/FED_JSATS_receivers.html", 
                           timeout(3)), error=function(e) print("Timeout"))
  
  # If the database isn't working then read csv
  if (x == TRUE | x == "Timeout") {
     ReceiverDeployments <- vroom("./data/rec_deployments.csv")
     
  } else {
    # If database is working then check for updates
    
    ## Download ReceiverDeployments
    my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
    JSATSinfo <- info('FED_JSATS_receivers', url = my_url)
    ReceiverDeployments <- tabledap(JSATSinfo, url = my_url)  
    
    # Fix column names and correct column types
    ReceiverDeployments <- ReceiverDeployments %>% 
      rename(
        SN = receiver_serial_number,
        GEN = receiver_general_location,
        Region = receiver_region,
        GPSname = receiver_location,
        LAT = latitude,
        LON = longitude,
        RKM = receiver_river_km,
        GenLat = receiver_general_latitude,
        GenLon = receiver_general_longitude,
        GenRKM = receiver_general_river_km,
        RecMake = receiver_make,
        StartTime = receiver_start,
        EndTime = receiver_end
      ) %>% 
      mutate_at(vars(SN, LAT, LON, RKM, GenLat, GenLon, GenRKM), as.numeric) %>% 
      mutate(
        StartTime = mdy_hms(StartTime),
        EndTime = mdy_hms(EndTime),
        water_year = ifelse(month(StartTime) <= 9, year(StartTime),
                            year(StartTime) + 1)
      ) %>% 
      filter(
        SN != 1
      )
    
    # Save latest update to file
    write_csv(ReceiverDeployments, "./data/rec_deployments.csv")
    
    # Change the last saved date to today
    last_checked_date <- Sys.Date()
    saveRDS(last_checked_date, "./data/last_checked_date.RDS")
  }
}

## Load Vemco Receiver Deployments -----------------------------------------

# This must be manually updated
VemcoReceiverDeployments <- vroom("./data/vemco_rec_deployments.csv",
                                  col_types = cols()) %>%
  left_join(ReceiverDeployments %>% 
              select(GPSname, GEN, GenLat, GenLon) %>% 
              distinct(),
            by = c("GPSname", "GEN", "GenLat", "GenLon") ) %>%
  mutate(
    water_year = ifelse(month(StartTime) <= 9, year(StartTime),
                        year(StartTime) + 1)
  )

# Assigns ReceiverDeployments an ID to be able to be identifiable when clicked
ReceiverDeployments$uid <- 1:nrow(ReceiverDeployments)
VemcoReceiverDeployments$uid <- 1:nrow(VemcoReceiverDeployments)


## Load StudyIDs -----------------------------------------------------------

# Table of studyid, descriptive name, whether it should be included in 
# survival tab, whether it should be included in non-survival tabs
studyid_dictionary <- vroom("./data/studyid_names.csv", col_types = cols())

# Descriptive names for studyids to be used in survival
surv_studyid_descript_name <- studyid_dictionary %>% 
  filter(incl_survival == "YES") %>% 
  pull(descript_name)

# Descriptive names for studyids to be used in survival
deltasurv_studyid_descript_name <- studyid_dictionary %>% 
   filter(incl_deltasurv == "YES") %>% 
   pull(descript_name)

# Descriptive names for studyids to be used in outmigration animation, 
# data explorer, time of day, movement
other_studyid_descript_name <- studyid_dictionary %>% 
  filter(incl_other == "YES") %>% 
  pull(descript_name)


## Load TaggedFish ---------------------------------------------------------

TaggedFish <- vroom("./data/tagged_fish.csv", col_types = cols())

setwd("./data")
# Update TaggedFish table from ERDDAP every 30 days
if (as.numeric(Sys.Date() - last_checked_date) < 30) { 
  my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
  JSATSinfo <- as.info('FED_JSATS_taggedfish', url = my_url)
  
  
  TaggedFish <- tabledap(JSATSinfo, url = my_url,
                         fields = c('fish_id', 'study_id', 'fish_type', 'fish_origin', 'fish_release_date', 'tag_id_hex', 'tag_id_decimal', 'tag_weight',
                                    'tag_model', 'tag_pulse_rate_interval_nominal', 'tag_warranty_life', 'fish_length_type', 'fish_length', 'fish_weight',
                                    'release_location', 'release_latitude', 'release_longitude', 'release_river_km'))

  # Fix column names and correct column types
  TaggedFish <- TaggedFish %>% 
    mutate_at(vars(tag_weight, tag_pulse_rate_interval_nominal, tag_warranty_life,
                   fish_length, fish_weight, release_latitude, release_longitude,
                   release_river_km), as.numeric) 
  
  closeAllConnections()
  # Save latest update to file
  write.csv(TaggedFish, "tagged_fish.csv")
}

TaggedFish <- vroom("tagged_fish.csv", col_types = cols())
setwd("..")

## Load CDEC hydrology data ------------------------------------------------
# Read saved CDEC flow data
comb_flow <- vroom("./data/comb_flow.csv", col_types = c(Index = "D", 
                                                         KES = "d", BND = "d", BTC = "d", WLK = "d",
                                                         MEN = "d", NEW = "d", VNS = "d", MSD = "d",
                                                         BND.Natural = "d"))

## ADD A NEW STATION OR SENSOR ##
## NOTE: Uncomment this chunk of code to add a new station or sensor
# new_station <- "BND"
# new_sensor  <- 8
# new_data <- cdec_query(stations = new_station, sensors = new_sensor, "D", min(comb_flow$Index), max(comb_flow$Index)) %>%
#             select(StationID, Index = DateTime, Value) %>%
#             mutate(Index = as.Date(format(Index, "%Y-%m-%d"))) %>%
#             distinct() %>%
#             pivot_wider(names_from = StationID, values_from = Value)
# names(new_data) <- c("Index", "BND.Natural") # Change names because we are adding full natural flow sensor for BND
# comb_flow <- comb_flow %>%
#              left_join(., new_data, by = "Index")

# # Update file with new flow data  from CDEC if it has been over 7 days since 
# # last download
if (as.numeric(Sys.Date() - max(comb_flow$Index)) > 7) {
  # Choose CDEC gauges to display
  gauges <- c("KES", "BND", "BTC", "WLK", "MEN", "NEW", "VNS", "MSD", "BND")

  # The last date of downloaded data in comb_flow + 1
  last_date <- max(comb_flow$Index) + 1

  hydro_err <- tryCatch({
     # apply the list of gauges to function that queries CDEC to get daily flow
     # then turns it into an xts (time series object) object which is needed for dygraphs
     flows <- cdec_query(stations = gauges, c(23, rep(41, 7), 8), "D", last_date, Sys.Date()) %>%
              mutate(StationID = ifelse(SensorNumber == 8, "BND.Natural", StationID)) %>%    
              select(StationID, Index = DateTime, Value) %>%
              mutate(Index = format(Index, "%Y-%m-%d")) %>%
              distinct() %>%
              pivot_wider(names_from = StationID, values_from = Value)
     
     write_csv(flows, "./data/comb_flow.csv", append = T)
     rm(flows)
     
  }, error = function(err) {
     return(1)
  })
} else {
   hydro_err <- 0
}

comb_flow <- as.xts(read.csv.zoo("./data/comb_flow.csv"))

cdec_stations <- vroom("./data/cdec_stations.csv", col_types = cols())

# Read saved CDEC temperature data
temps <- vroom("./data/water_temps.csv", col_types = c(Index = "D", BND = "d", KWK = "d", MSD = "d", WLK = "d")) %>%
         select(Index, BND, KWK, MSD, WLK)

# Update file with new water temp data from CDEC if it has been over 7 days since last download
if (as.numeric(Sys.Date() - max(temps$Index)) > 7) {
   # Choose CDEC gauges to display
   wtemp_gauges <- c("BND", "KWK", "MSD", "WLK")
   
   # The last date of downloaded data in temps + 1
   last_date <- max(temps$Index) + 1
   
   wtemp_err <- tryCatch({
      # apply the list of gauges to function that queries CDEC to get daily temperature
      # then turns it into an xts (time series object) object which is needed for dygraphs
      new_temps <- cdec_query(stations = wtemp_gauges, sensors = rep(25, length(wtemp_gauges)), start.date = last_date, end.date = Sys.Date()) %>%
                   mutate(Index = format(DateTime, "%Y-%m-%d")) %>%
                   group_by(StationID, Index) %>%
                   summarise(temp = mean(Value, na.rm = TRUE)) %>%
                   mutate(temp = (temp - 32)/1.8) %>%
                   pivot_wider(names_from = StationID, values_from = temp)
      
      write_csv(rbind(temps, new_temps), "./data/water_temps.csv")
      rm(new_temps)
      
   }, error = function(err) {
      return(1)
   })
} else {
   wtemp_err <- 0
}

temps <- vroom("./data/water_temps.csv", col_types = c(Index = "D", BND = "d", KWK = "d", MSD = "d", WLK = "d")) %>%
         select(Index, BND, KWK, MSD, WLK)


##################### SECTION FOR RUSS ET AL #####################
# Use this section for loading data and data manipulation

##################################################################


## Helper functions --------------------------------------------------------

# Convert a studyid descriptive name to short studyid name
convert_descriptive_to_studyid <- function(descript_studyid) {
  studyid_dictionary$study_id[studyid_dictionary$descript_name == descript_studyid]
}


# UI ----------------------------------------------------------------------

ui <- fluidPage(theme = shinytheme("flatly"),
                navbarPage("Central Valley Enhanced Acoustic Tagging Project",
                           tabPanel("Background", 
                                    mainPanel(
                                      uiOutput("background")
                                    )
                           ),
                           tabPanel("Receiver Deployments",
                                    sidebarLayout(
                                      sidebarPanel(
                                        radioButtons("receiverType", 
                                                     "Receiver type",
                                                     c("JSATS", 
                                                       "Real-time", "Vemco")
                                        ),
                                        selectInput("water_year",
                                                    "Water Year:",
                                                    ""),
                                        helpText("Note: water year is defined as the 12 month period starting 
                                   October 1 to September 30 of the following calendar year. 
                                   The water year is designated by the calendar year in which 
                                   it ends and which includes 9 of the 12 months. Thus, the year
                                   ending September 30, 1999 is called the 1999 water year."),
                                        htmlOutput("map_marker_click")
                                      ),
                                      mainPanel(
                                        leafletOutput("map", width = "100%", 
                                                      height = "650"),
                                        dataTableOutput("receiver_table")
                                      )
                                    )
                           ),
                           tabPanel("Hydrology",
                                    fluidRow(style='margin: 0px',
                                             column(6, style = "padding: 0px", align = "center",
                                                    box(dygraphOutput("plot1", width = "100%"), width = 12
                                                    ),
                                                    fluidRow(style = "margin: 0px", align = "center",
                                                             box(dygraphOutput("plot3", width = "100%"), width = 12
                                                             )
                                                    )
                                             ),
                                             column(6, style = "padding: 0px", align = "center",
                                                    box(dygraphOutput("plot2", width = "100%"), width = 12
                                                        ),
                                                    fluidRow(style = "margin: 0px", align = "center",
                                                             box(dygraphOutput("plot4", width = "100%"), width = 12
                                                             )
                                                    )
                                             )
                                    ),
                                    fluidRow(style="margin: 5px",
                                             column(12, style = "padding: 0px", align = "center",
                                                    box(leafletOutput("hydromap", width = "75%", height = "650px"), width = 12
                                                        )
                                                    )
                                    ),
                                    textOutput('hydrology_error'),
                                    tags$head(tags$style("#hydrology_error{color: red;
                                                         font-size: 20px;
                                                         }")
                                    ),
                                    textOutput('temperature_error'),
                                    tags$head(tags$style("#temperature_error{color: red;
                                                         font-size: 20px;
                                                         }")
                                    )
                                    # textOutput("text1"),
                                    # textOutput("text2"),
                                    # textOutput("text3"),
                                    # downloadButton("downloadData", "Download")

                           ),
                           tabPanel("Outmigration Animation",
                                    sidebarLayout(
                                      sidebarPanel(
                                        selectInput("anim_dataset", 
                                                    "Choose a study population",
                                                    choices = other_studyid_descript_name),
                                        helpText("This tool visualizes study group detection data into an animated time series",
                                                 "of fish outmigration. Unique fish detections are identified at each general",
                                                 "receiver location and summed by day. Note that is it possible for a fish to",
                                                 "be detected at more than one general receiver location in a single day.")
                                      ),
                                      mainPanel(
                                        leafletOutput("timemap", width = "100%", 
                                                      height = "650")
                                      )
                                    )
                           ),
                           tabPanel("Data Explorer",
                                    headerPanel("Select Options"),
                                    sidebarPanel(
                                      selectizeInput("data_explorer_datasets", 
                                                     "Study Populations",
                                                     # Adding empty string first 
                                                     # allows no default selection
                                                     choices = c("", other_studyid_descript_name),
                                                     options = list(maxItems = 4)),
                                      selectInput("variable", "Variable",
                                                  choices = c("Weight", "Length")),
                                      selectInput("plot_type", "Plot type", 
                                                  choices = c("boxplot", 
                                                              "histogram", 
                                                              "density")),
                                      # Show this bin width slider only if 
                                      # histogram is chosen
                                      conditionalPanel(
                                        condition = "input.plot_type == 'histogram'",
                                        # Place holder for bin slider input, 
                                        # created in server side
                                        uiOutput("bin_slider")
                                      )
                                    ),
                                    mainPanel(
                                      plotlyOutput("plotly_plot"),
                                      tableOutput("summarytable")
                                    )
                           ),
                           tabPanel("Time of Day",
                                    headerPanel("Select Options"),
                                    mainPanel(
                                      selectInput("time_of_day_input", 
                                                  "Choose a study population",
                                                  choices = c("", other_studyid_descript_name),
                                      ),
                                      radioButtons("time_of_day_radio", 
                                                   "Choose all detections or by receiver GEN",
                                                   c("All Detections", "By General Location")
                                      ),
                                      conditionalPanel(
                                        condition = "input.time_of_day_radio == 'By General Location'",
                                        selectInput("time_of_day_GEN", "Select a General Location",
                                                    "")
                                      ),
                                      plotlyOutput("time_of_day_plot"),
                                      uiOutput("time_of_day_caption")
                                    )
                           ),
                           tabPanel("Survival",
                                    tabsetPanel(
                                      tabPanel("Cumulative Survival",
                                               headerPanel("Select Options"),
                                               sidebarPanel(
                                                 uiOutput("cumSurvSelect"),
                                                 helpText("Note: These survival results are
                                          preliminary and for discussion purposes
                                          only. Detection data has not been
                                          filtered for predator detections, and
                                          survival estimates have not been
                                          adjusted for any potential premature
                                          tag failures."),
                                                 uiOutput("cumsurvival_text"),
                                                 radioButtons("cumsurvival_radio",
                                                              "View",
                                                              c("Plot", "Table")
                                                 ),
                                                 leafletOutput("survival_map2")
                                               ),
                                               mainPanel(
                                                 conditionalPanel(
                                                   condition = "input.cumsurvival_radio == 'Plot'",
                                                   plotlyOutput("plotly_survival_output")
                                                 ),
                                                 conditionalPanel(
                                                   condition = "input.cumsurvival_radio == 'Table'",
                                                   dataTableOutput('cumSurvDT'),
                                                   downloadButton("cumSurvDownload")
                                                 )
                                               )
                                      ),
                                      tabPanel("Reach Survival",
                                               headerPanel("Select Options"),
                                               sidebarPanel(
                                                 uiOutput("reachSurvSelect"),
                                                 helpText("Note: These survival results are
                                          preliminary and for discussion purposes
                                          only. Detection data has not been
                                          filtered for predator detections, and
                                          survival estimates have not been
                                          adjusted for any potential premature
                                          tag failures."),
                                                 uiOutput("reachSurv_text"),
                                                 radioButtons("reachSurvRadio",
                                                              "View",
                                                              c("Plot", "Table")
                                                 ),
                                                 leafletOutput("reachSurvMap")
                                               ),
                                               mainPanel(
                                                 conditionalPanel(
                                                   condition = "input.reachSurvRadio == 'Plot'",
                                                   plotlyOutput("reachSurvPlotly")
                                                 ),
                                                 conditionalPanel(
                                                   condition = "input.reachSurvRadio == 'Table'",
                                                   dataTableOutput('reachSurvDT'),
                                                   downloadButton("reachSurvDownload")
                                                 )
                                               )
                                      ),
                                      
                                      
                                      ##################### USER-INTERFACE FOR DELTA ROUTE-SPECIFIC SURVIVAL #####################
                                      tabPanel("Delta route-specific survival",
                                               tags$head(tags$style("svg { pointer-events: none;}")),
                                               headerPanel("Select Options"),
                                               sidebarPanel(
                                                  uiOutput("deltaSurvSelect"),
                                                  helpText("Note: These survival results are
                                          preliminary and for discussion purposes
                                          only. Detection data has not been
                                          filtered for predator detections, and
                                          survival estimates have not been
                                          adjusted for any potential premature
                                          tag failures."),
                                          radioButtons("deltaSurvRadio",
                                                       "View",
                                                       c("Diagram", "Barplot", "Table")
                                          ),
                                          leafletOutput("deltasurv_map")
                                               ),
                                          mainPanel(
                                             conditionalPanel(
                                                condition = "input.deltaSurvRadio == 'Diagram'",
                                                textOutput("deltasurv_diagramtext"),
                                                tags$head(tags$style("#deltasurv_diagramtext{font-size: 20px;}")),
                                                grVizOutput("deltasurv_diagram", width = "900", height = "650px")),
                                             conditionalPanel(
                                                condition = "input.deltaSurvRadio == 'Barplot'",
                                                textOutput("deltasurv_barplottext"),
                                                tags$head(tags$style("#deltasurv_barplottext{font-size: 16px;}")),
                                                plotOutput("deltasurv_barplot", width = "400px", height = "400px")),
                                             conditionalPanel(
                                                condition = "input.deltaSurvRadio == 'Table'",
                                                tableOutput("deltasurv_table"))
                                             
                                          )
                                      )
                                      
                                      ############################################################################################
                                    )
                                    
                           ),
                           tabPanel("Movement",
                                    tabsetPanel(
                                      tabPanel("Travel Time Table",
                                               headerPanel("Select Options"),
                                               sidebarPanel(
                                                 selectInput(
                                                   "movement_dataset", "Choose a study population",
                                                   choices = other_studyid_descript_name
                                                 ),
                                                 uiOutput(
                                                   "movement_second_select"
                                                 )
                                               ),
                                               mainPanel(
                                                 gt_output(
                                                   "movement_gt"
                                                 )
                                               )
                                      )
                                      # ,
                                      # tabPanel("Movement Plot",
                                      #          headerPanel("Select Options"),
                                      #          sidebarPanel(
                                      #            selectInput(
                                      #              "movement_plot_dataset", "Choose a study population",
                                      #              choices = studyid_list
                                      #            )
                                      #          ),
                                      #          mainPanel(
                                      #            plotlyOutput(
                                      #              "movement_plot"
                                      #            )
                                      #          )
                                      # )
                                    )
                                    
                           )
                )
)


# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  
  ##%######################################################%##
  #                                                          #
  ####                    Background                      ####
  #                                                          #
  ##%######################################################%##
  
  output$background <- renderUI({
    withTags({
      div(class="header", checked=NA,
          h1("Background"),
          p("Since 2012, juvenile salmon have been tagged and tracked throughout 
            California's Central Valley (CCV) using Juvenile Salmon Acoustic 
            Telemetry System (JSATS) technology. This technology allows 
            researchers to monitor the movement and survival rates of various 
            populations over 500 river kilometers, from Redding to the Pacific 
            Ocean. The data compiled here is a result of the hard work and 
            coordination between various federal and state agencies, universities, 
            and private consultants, with the goal of conserving and restoring 
            California's once abundant, but now imperiled, Chinook salmon 
            populations. This data is open access and is managed and hosted by 
            the National Marine Fisheries Service here:", 
            a(href="https://oceanview.pfeg.noaa.gov/erddap/tabledap/FED_JSATS_detects.html", 
              "https://oceanview.pfeg.noaa.gov/erddap/tabledap/FED_JSATS_detects.html"),
            "."),
          br(),
          h3("Receiver Deployments"),
          p("An interactive map of acoustic receivers (JSATS, Vemco, Real-time) 
            deployed throughout the CCV by water year, which is defined as the 
            12 month period beginning October 1 and ending September 30 of the 
            following calendar year (i.e. water year 2017 spans 10/1/17-9/30/18). 
            The deployment period of individual receivers is displayed in a 
            table below the map, once a general location is clicked. This map is 
            useful to identify a particular site of interest, and to see when 
            coverage existed in that area and which receiver(s) potentially 
            recorded migrating fish. "),
          br(),
          h3("Hydrology"),
          p("The hydrology of the Sacramento River is highly variable, and is 
            largely driven by storm events and associated runoff in the winter, 
            followed by dam controlled releases for agricultural purposes in the 
            summer. These interactive graphs of Sacramento (SR) and San Joaquin River (SJR) flows (cubic feet 
            per second) at four different ", 
            a(href="https://cdec.water.ca.gov/index.html", "CDEC"), 
            " (California Data Exchange Center) stations on each river: ", 
            a(href="https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=KES", "KES"), 
            " (Keswick Reservoir, SR), ", 
            a(href="https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=BND", "BND"), 
            " (Bend Bridge, SR), ",
            a(href="https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=BTC", "BTC"), 
            " (Butte City, SR), ",
            a(href="https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=WLK", "WLK"), 
            " (Wilkins Slough, SR), ",
            a(href="https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=MEN", "MEN"), 
            " (Mendota, SJR), ",
            a(href="https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=NEW", "NEW"), 
            " (Newman, SJR), ",
            a(href="https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=VNS", "VNS"), 
            " (Vernalis, SJR), and",
            a(href="https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=MSD", "MSD"), 
            " (Mossdale Bridge, SJR) displays the annual and inter-annual variation 
            in flows from upstream to downstream gauging stations. In addition, the temperatures at KES, BND,
            WLK, and MSD are shown. These plots were created using the R package ",
            a(href="https://rstudio.github.io/dygraphs/", "dygraphs"),"."),
          br(),
          h3("Outmigration Animation"),
          p("The Outmigration Animation tab visualizes fish outmigration using 
            detection data over time. Unique fish detections at each receiver 
            location are plotted by day, and numbers of unique detections per 
            receiver are overlaid in each point. This animation helps to 
            visualize the movement of fish across a broad geographic landscape, 
            and demonstrates the differences in outmigration timing across space 
            and time among specific populations. This animation was created 
            using the R packages  ", 
            a(href="https://rstudio.github.io/leaflet/", "leaflet"),
            "and ",
            a(href="https://github.com/bhaskarvk/leaflet.extras", "leaflet.extras"),
            "."),
          br(),
          h3("Data Explorer"),
          p("Data explorer allows users to look at the tagged fish data 
            individually or across different study populations in a number of 
            ways. Plots can be created interactively, choosing variables to 
            explore (length, weight), and with different plot types (boxplot, 
            histogram, density) as well as summary tables."),
          br(),
          h3("Time of Day"),
          p("Time of day allows users to visually explore behavioral differences 
            in fish movement between night and day throughout the migratory 
            corridor. Using detection times we can look at the distribution in 
            movement times for an entire study group, or the movement times for 
            a specific receiver location. "),
          br(),
          h3("Survival"),
          p("Survival estimates are calculated by summarizing the detection 
            history of study populations, and assuming a fish has died when it 
            is not detected by subsequent downstream receivers. We used a CJS 
            survival model in RMark to estimate reach specific (survival per 10 
            river kilometers) and cumulative (from release location to last 
            downstream detection) survival rates. We used a simple survival model 
            (survival and detection efficiency is a function of time) to derive 
            these estimates. For each survival tab, a map is displayed with all 
            receiver locations used to generate survival for each population, 
            and a figure and table displays the estimates and associated error. 
            The unique number of fish detected at each receiver location is also 
            provided."),
          br(),
          p("For the Delta route-specific survival tab, we ran a Bayesian 
            multi-state mark-recapture model in JAGS to obtain survival and 
            route-entrainment probability estimates for the Delta between 
            Sacramento and Benicia. We modeled routing and survival in four 
            possible routes through the Delta: Sacramento River mainstem, 
            Sutter Slough, Steamboat Slough, and Georgiana Slough. In this tab, 
            there is a diagram showing survival and routing estimates for each 
            reach and junction in the Delta, and a barplot and table summarizing 
            survival and migration route probabilities for each of the four 
            migratory pathways. This section was developed by the Quantitative 
            Fisheries Ecology Section of USGS Western Fisheries Research Center. 
            For questions or comments, send an email to our group at 
            gs-b-crrl_qfes@usgs.gov ."),
          br(),
          i("*These survival results are preliminary and for discussion purposes 
            only. Detection data has not been filtered for predator detections, 
            and survival estimates have not been adjusted for any potential 
            premature tag failures."),
          br(),
          h3("Movement"),
          p("Fish movement is summarized by study population in a table format. 
            For each receiver location, the minimum, median, and maximum travel 
            time is calculated in days and kilometers per day. The number of 
            unique fish detected at each receiver location is displayed as well. 
            We will update this page to display travel times in an interactive 
            plot soon, so stay tuned."),
          h3("Questions or comments?"),
          p("This Shiny App was developed by UCSC/NOAA, Southwest Science 
            Fisheries Center, Fisheries Ecology Division. Source code used to create 
            this site can be viewed", 
            a(href="https://github.com/pham-thomas/shiny-telemetry", "here"),
            ".", " If you have any
            questions or comments please feel free to send us an ",
            a(href="mailto:jeremy.notch@noaa.gov", "email"),
            ".")
      )
    })
  })
  
  ##%######################################################%##
  #                                                          #
  ####               Receiver Deployments                 ####
  #                                                          #
  ##%######################################################%##
  
  # Updates allowable water_year selection input based on radio button selection
  observe({
    updateSelectInput(session, "water_year",
                      choices = outVar()
    )
  })
  
  # Reactive list of years depending on radio button selection
  outVar <-  reactive({
    receiver_type_selection <- as.character(input$receiverType)
    
    # If the radio button is "Real-time" get the years they were deployed
    if (receiver_type_selection == "Real-time") {
      vars <- ReceiverDeployments %>%
        filter(RecMake %in% c("ATS SR3017", "Tekno RT") & 
                 SN != "999" & GPSname != "999") %>% 
        select(water_year) %>% 
        distinct() %>% 
        drop_na() %>% 
        pull(water_year)
      # If the radio button is "Autonomous" get the years they were deployed
    }else if(receiver_type_selection == "JSATS") {
      vars <- ReceiverDeployments %>%
        filter(!(RecMake %in% c("ATS SR3017", "Tekno RT")) &
                 SN != "999" & GPSname != "999") %>% 
        select(water_year) %>% 
        distinct() %>% 
        drop_na() %>% 
        pull(water_year)
    }else if(receiver_type_selection == "Vemco") {
      vars <- VemcoReceiverDeployments %>%
        select(water_year) %>% 
        distinct() %>% 
        drop_na() %>% 
        pull(water_year)
    }
    sort(vars)
  })
  
  # Create receiver deployments map
  output$map <- renderLeaflet({
    
    # Use user input to query receiver type and assign colors
    if (input$receiverType == "JSATS") {
      receivers <- subset(ReceiverDeployments, water_year == input$water_year & 
                            !(RecMake %in% c("ATS SR3017", "Tekno RT")))
      color <- "purple"
    }else if (input$receiverType == "Real-time"){
      receivers <- subset(ReceiverDeployments, water_year == input$water_year & 
                            RecMake %in% c("ATS SR3017", "Tekno RT"))
      color <- "red"
    }else if(input$receiverType == "Vemco"){ 
      receivers <- subset(VemcoReceiverDeployments, water_year == input$water_year)
      color <- "orange"
    }
    
    receivers <- subset(receivers, GPSname != "999")
    
    # Create the leaflet map
    leaflet(data = receivers) %>% 
      addTiles(group = "Open Street Map") %>% 
      # Provide optional basemap tiles
      addProviderTiles(providers$Stamen.Terrain, group = "Stamen Terrain",
                       options = providerTileOptions(noWrap = TRUE)) %>% 
      addProviderTiles(providers$Esri.NatGeoWorldMap, group = "Esri Nat Geo",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      # Set the default view
      setView(-122.159729,39.119407, 7) %>% 
      # Display receiver deploments as circle markers
      addCircleMarkers(~GenLon, ~GenLat, 
                       popup = paste(receivers$GEN),
                       layerId = ~uid,
                       radius = 3,
                       stroke = FALSE,
                       color = color, 
                       fillOpacity = 1) %>% 
      # Give control for selecting basemaps
      addLayersControl(
        baseGroups = c("Open Street Map", "Stamen Terrain", "Esri Nat Geo", 
                       "Esri World Imagery"),
        options = layersControlOptions(collapsed = FALSE)
      ) %>% 
      # Button to reset to original view
      addResetMapButton()
  })
  
  prevRecType <- reactiveVal("")
  prevYear <- reactiveVal("")
  
  # Produces table output in response to clicks on receiver deployment map markers
  observeEvent(input$map_marker_click, {
    # Essentially the attributes of the clicked receiver marker is assigned to p
    p <- input$map_marker_click
    
    if (input$receiverType == "Vemco") {
      
      # Give me the GEN for the clicked marker
      GEN_click = VemcoReceiverDeployments$GEN[VemcoReceiverDeployments$uid == p$id]
      
      prevRecType(input$receiverType)
      prevYear(input$water_year)
      
      df <- VemcoReceiverDeployments %>%
        filter(
          GEN == GEN_click & water_year == input$water_year
        ) %>%
        select(
          GPSname, VemcoSN, StartTime, EndTime, GenLat, GenLon
        ) %>%
        arrange(StartTime) %>%
        mutate(
          VemcoSN = as.character(VemcoSN),
          StartTime = as.character(StartTime),
          EndTime = as.character(EndTime)
        )
      
      output$receiver_table <- renderDataTable(df)
      
      
    } else {
      GEN_click = ReceiverDeployments$GEN[ReceiverDeployments$uid == p$id]
      
      prevRecType(input$receiverType)
      prevYear(input$water_year)
      
      # Display all the GPSnames that are associated with that GEN for the given
      # water year along with SN, StartTime, and EndTime
      df <-  ReceiverDeployments %>%
        filter(
          GEN == GEN_click & water_year == input$water_year
        ) %>%
        arrange(StartTime) %>%
        mutate(
          SN = as.character(SN),
          StartTime = as.character(StartTime),
          # If the receiver is RT and its EndTime is NA make it say 'Active' instead
          EndTime = ifelse(RecMake %in% c("ATS SR3017", "Tekno RT") & is.na(EndTime),
                           "Active", as.character(EndTime))
        ) %>%
        select(
          GPSname, SN, StartTime, EndTime, GenLat, GenLon
        )
      
      output$receiver_table <- renderDataTable(df)
    }
    
  })
  
  # If user switches receiver type display an empty table
  observeEvent(input$receiverType, {
    if (!is_empty(prevRecType()) & input$receiverType != prevRecType()) {
      df <- datatable(
        data.frame(GPSname = NA, 
                   SN = NA, 
                   StartTime = NA, 
                   EndTime = NA,
                   GenLat = NA,
                   GenLon = NA)
      )
      
      output$receiver_table <- renderDataTable(df)
      
    }
  })
  
  # If user switches year display an empty table
  observeEvent(input$water_year, {
    if (!is_empty(prevYear()) & input$water_year != prevYear()) {
      df <- datatable(
        data.frame(GPSname = NA, 
                   SN = NA, 
                   StartTime = NA, 
                   EndTime = NA,
                   GenLat = NA,
                   GenLon = NA)
      )
      
      output$receiver_table <- renderDataTable(df)
      
    }
  })
  
  
  ##%######################################################%##
  #                                                          #
  ####                     Hydrology                      ####
  #                                                          #
  ##%######################################################%##
  
  # Create a dygraph - interactive flow graph of different CDEC gauges: KES, BND, BTC, WLK
  # output$dygraph <- renderDygraph({
  #   dygraph(comb_flow, main = "Sacramento River Flows (CFS)") %>%
  #     # Adds date range selector tool, use a year from today through todays date as the default range
  #     dyRangeSelector(dateWindow = c(Sys.Date() - 365, Sys.Date())) %>%
  #     # Individual series customizations
  #     dySeries("KES", strokeWidth = 1.5, strokePattern = "solid", color = "#008033") %>%
  #     dySeries("BND", strokeWidth = 1.5, strokePattern = "solid", color = "#668000") %>%
  #     dySeries("BTC", strokeWidth = 1.5, strokePattern = "solid", color = "#003380") %>%
  #     dySeries("WLK", strokeWidth = 1.5, strokePattern = "solid", color = "#660080") %>%
  #     dyOptions(fillGraph = TRUE, fillAlpha = 0.2) %>%
  #     # Make legend wide enough to fit single line
  #     dyLegend(width = 400) %>%
  #     # Makes it so that mouse hovered series is highlighted
  #     dyHighlight(highlightSeriesOpts = list(strokeWidth = 3))
  # })
  
  if(hydro_err == 1){
     output$hydrology_error <- renderText({'ERROR RETRIEVING DATA FROM CDEC'})
  } else {
     output$hydrology_error <- renderText({''})
  }
  
  if(wtemp_err == 1){
     output$temperature_error <- renderText({'ERROR RETRIEVING DATA FROM CDEC'})
  } else {
     output$temperature_error <- renderText({''})
  }
  
  
  comb_flow1 <- comb_flow[,-c(5:8)]
  output$plot1 <- renderDygraph({
     dygraph(comb_flow1, main = "Sacramento River Flows (CFS)", group = "sac") %>%
        # Adds date range selector tool, use a year from today through todays date as the default range
        dyRangeSelector(dateWindow = c(Sys.Date() - 365, Sys.Date())) %>%
        # Individual series customizations
        dySeries("KES", strokeWidth = 1.5, strokePattern = "solid", color = "#008033") %>%
        dySeries("BND", strokeWidth = 1.5, strokePattern = "solid", color = "#668000") %>%
        dySeries("BND.Natural", strokeWidth = 1.5, strokePattern = "dashed", color = "#668000") %>%
        dySeries("BTC", strokeWidth = 1.5, strokePattern = "solid", color = "#003380") %>%
        dySeries("WLK", strokeWidth = 1.5, strokePattern = "solid", color = "#660080") %>%
        dyOptions(fillGraph = TRUE, fillAlpha = 0.2) %>%
        # Make legend wide enough to fit single line
        dyLegend(width = 500) %>%
        # Makes it so that mouse hovered series is highlighted
        dyHighlight(highlightSeriesOpts = list(strokeWidth = 3))
  })

  temps <- as.data.frame(temps)
  row.names(temps) <- temps[,1]
  temps1 <- temps[,c(2,3,5)]
  names(temps1) <- c("BND", "KES", "WLK")
  output$plot3 <- renderDygraph({
     dygraph(temps1, main = "Sacramento River Temperature (C)", group = "sac") %>%
        # Adds date range selector tool, use a year from today through todays date as the default range
        dyRangeSelector(dateWindow = c(Sys.Date() - 365, Sys.Date())) %>%
        # Individual series customizations
        dySeries("KES", strokeWidth = 1.5, strokePattern = "solid", color = "#008033") %>%
        dySeries("BND", strokeWidth = 1.5, strokePattern = "solid", color = "#668000") %>%
        dySeries("WLK", strokeWidth = 1.5, strokePattern = "solid", color = "#660080") %>%
        dyOptions(fillGraph = TRUE, fillAlpha = 0.2) %>%
        # Make legend wide enough to fit single line
        dyLegend(width = 500) %>%
        # Makes it so that mouse hovered series is highlighted
        dyHighlight(highlightSeriesOpts = list(strokeWidth = 3))
  })

  comb_flow2 <- comb_flow[,-c(1:4,9)]
  output$plot2 <- renderDygraph({
     dygraph(comb_flow2, main = "San Joaquin River Flows (CFS)", group = "sj") %>%
        # Adds date range selector tool, use a year from today through todays date as the default range
        dyRangeSelector(dateWindow = c(Sys.Date() - 365, Sys.Date())) %>%
        # Individual series customizations
        dySeries("MEN", strokeWidth = 1.5, strokePattern = "solid", color = "#801500") %>%
        dySeries("NEW", strokeWidth = 1.5, strokePattern = "solid", color = "#007180") %>%
        dySeries("VNS", strokeWidth = 1.5, strokePattern = "solid", color = "#800071") %>%
        dySeries("MSD", strokeWidth = 1.5, strokePattern = "solid", color = "#919F27") %>%
        dyOptions(fillGraph = TRUE, fillAlpha = 0.2) %>%
        # Make legend wide enough to fit single line
        dyLegend(width = 500) %>%
        # Makes it so that mouse hovered series is highlighted
        dyHighlight(highlightSeriesOpts = list(strokeWidth = 3))
  })
  
  temps2 <- temps[,c(1,4)]
  temps2 <- xts(temps2[,-1], order.by = temps2[,1])
  names(temps2) <- "MSD"
  output$plot4 <- renderDygraph({
     dygraph(temps2, main = "San Joaquin River Temperature (C)", group = "sj") %>%
        # Adds date range selector tool, use a year from today through todays date as the default range
        dyRangeSelector(dateWindow = c(Sys.Date() - 365, Sys.Date())) %>%
        # Individual series customizations
        dySeries("MSD", strokeWidth = 1.5, strokePattern = "solid", color = "#919F27") %>%
        dyOptions(fillGraph = TRUE, fillAlpha = 0.2) %>%
        # Make legend wide enough to fit single line
        dyLegend(width = 500) %>%
        # Makes it so that mouse hovered series is highlighted
        dyHighlight(highlightSeriesOpts = list(strokeWidth = 3))
  })
  

  # Create leaflet map that corresponds to the dygraph, it displays the gauge locations
  output$hydromap <- renderLeaflet({
    # Set colors for each gauge
    pal <- colorFactor(
      palette = c("#668000", "#003380", "#008033", "#801500","#919F27","#007180", "#800071", "#660080"),
      domain = c( "BND","BTC","KES","MEN","MSD","NEW","VNS","WLK")
    )

    # Bring in a Sacramento River shapefile
    rivers   <- st_read("./data/GIS/sac_river_dissolve.shp")
    # Transform datum
    rivers   <- st_transform(rivers, crs = '+proj=longlat +datum=WGS84')
    # Need to drop z dimension: https://gis.stackexchange.com/questions/253898/adding-a-linestring-by-st-read-in-shiny-leaflet
    rivers   <- st_zm(rivers, drop = T, what = "ZM")

    leaflet(data = cdec_stations) %>%
      addPolygons(
        data = rivers, color = "#668db3",
        weight = 3, opacity = 1
      ) %>%
      addProviderTiles(providers$Stamen.Terrain, group = "Stamen Terrain",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      addTiles(group = "Open Street Map") %>%
      addProviderTiles(providers$Esri.NatGeoWorldMap, group = "Esri Nat Geo",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      setView(mean(cdec_stations$longitude), mean(cdec_stations$latitude), 7) %>%
      addCircleMarkers(
        ~longitude, ~latitude,
        popup = ~name,
        color = ~pal(station_id),
        fillOpacity = 1.5,
        label = ~station_id,
        # Make labels always appear and offset them to the right of the markr
        labelOptions = labelOptions(
          noHide = T,
          direction = "right",
          offset = c(10,0))) %>%
      addLayersControl(
        baseGroups = c("Stamen Terrain", "Open Street Map", "Esri Nat Geo",
                       "Esri World Imagery"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })

  # Explanatory text for the download button
  output$text1 <- renderText("Click the button below to download the hydrology data,
                             in CSV, currently in view based on the date range selector.")
  output$text2 <- renderText(paste0("Start Date: ", as.Date(ymd_hms(input$dygraph_date_window[[1]]))))
  output$text3 <- renderText(paste0("End Date: ", as.Date(ymd_hms(input$dygraph_date_window[[2]]))))


  # Reactive value of hydrology data based on date range selection
  hydroDataset <- reactive({
    start_date <- as.Date(ymd_hms(input$dygraph_date_window[[1]]))
    end_date <- as.Date(ymd_hms(input$dygraph_date_window[[2]]))

    # First convert comb_flow xts to a dataframe and extract index from rowname into a
    # column
    comb_flow_subset <- data.frame(date=index(comb_flow), coredata(comb_flow))

    # Filter the data based on date range selection values
    comb_flow_subset %>%
      filter(
        date >= start_date,
        date <= end_date
      )
  })

  # Adds download button and allows user to download csv of the hydrology data
  # Currently in, view filtered by user's date range selected
  output$downloadData <- downloadHandler(
    filename = function() {
      start_date <- as.Date(ymd_hms(input$dygraph_date_window[[1]]))
      end_date <- as.Date(ymd_hms(input$dygraph_date_window[[2]]))
      paste(start_date, "_", end_date, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(hydroDataset(), file, row.names = FALSE)
    }
  )
  
  ##%######################################################%##
  #                                                          #
  ####               Outmigration Animation               ####
  #                                                          #
  ##%######################################################%##
  
  # Reactive expression that retrieves detections from ERDDAP and formats it for 
  # the time step animation
  timestepVar<-  reactive({
    req(input$anim_dataset)
    
    files <- list.files("./detections/", full.names = T)
    
    # Choose the file that matches the studyID and read it in
    original_studyid <- convert_descriptive_to_studyid(input$anim_dataset)
    file <- files[str_detect(files, original_studyid)]
    detections <- vroom(file, col_types = cols())
    
    # Get release info from TaggedFish
    release <- TaggedFish %>% 
      filter(study_id == original_studyid) %>% 
      select(release_latitude, release_longitude, release_river_km) %>% 
      distinct()
    
    detections <- detections %>% 
      mutate(
        date = as.Date(time_pst),
        GenLat = ifelse(SN == 1, LAT, GenLat),
        GenLon = ifelse(SN == 1, LON, GenLon),
        GenRKM = ifelse(SN == 1, release$release_river_km, GenRKM)
      )
    
    # Get the list of unique receiver general locations for the study group
    receiver_loc <- detections %>% 
      select(GEN, GenRKM, GenLat, GenLon) %>% 
      arrange(desc(GenRKM)) %>% 
      distinct()
    
    # Create a grid of all receiver general locations by all dates from the 
    # earliest detection date to the latest detection date
    timestep <- expand.grid(receiver_loc$GEN, seq(min(detections$date), 
                                                  max(detections$date), by = 1),
                            stringsAsFactors = F)
    colnames(timestep) <- c("GEN", "date")
    
    # Summarise detections counts by GEN and date and join into the grid
    timestep <- timestep %>% 
      left_join(detections %>% 
                  select(FishID, GEN, date) %>% 
                  distinct() %>% 
                  group_by(GEN, date) %>% 
                  summarise(num_fish = n()), by = c('GEN', 'date')
      ) %>% 
      mutate(
        num_fish = ifelse(is.na(num_fish), 0, num_fish)
      ) %>% 
      left_join(receiver_loc %>% 
                  select(-GenRKM), 
                by = 'GEN')
  })
  
  # Output leaflet Outmigration Animation of fish outmigration
  output$timemap <- renderLeaflet({
    data <- timestepVar()
    
    leaflet(data = data, width = "100%", height = "800px") %>%
      addProviderTiles(providers$Stamen.Terrain, group = "Stamen Terrain",
                       options = providerTileOptions(noWrap = TRUE)) %>% 
      addTiles(group = "Open Street Map") %>% 
      addProviderTiles(providers$Esri.NatGeoWorldMap, group = "Esri Nat Geo",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      setView(mean(data$GenLon), mean(data$GenLat), 7) %>%
      addMinicharts(
        data$GenLon, data$GenLat,
        chartdata = data$num_fish,
        time = data$date,
        fillColor = "blue",
        width = 60, height = 60,
        popup = popupArgs(
          showValues = FALSE,
          supValues = data %>% select(GEN, num_fish),
          supLabels = c("GEN", "N = ")
        ),
        showLabels = TRUE,
        opacity = .7
      ) %>% 
      # Give control for selecting basemaps
      addLayersControl(
        baseGroups = c("Stamen Terrain", "Open Street Map", "Esri Nat Geo", 
                       "Esri World Imagery"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  ##%######################################################%##
  #                                                          #
  ####                   Data Explorer                    ####
  #                                                          #
  ##%######################################################%##
  
  # Reactive function that returns the filtered TaggedFish table by selected studyID
  taggedfishVar <-  reactive({
    req(input$data_explorer_datasets)
    
    # List of studyIDs to query for
    original_studyids <- studyid_dictionary$study_id[studyid_dictionary$descript_name 
                                                     %in% input$data_explorer_datasets]
    
    # Filter TaggedFish by studyID's
    TaggedFish %>% 
      filter(
        study_id %in% original_studyids,
        fish_weight > 0 # filter out placeholder values
      ) %>% 
      rename(
        StudyID = study_id,
        FishID = fish_id,
        Length = fish_length,
        Weight = fish_weight
      )
    
  })
  
  # Finds the optimal bin size using Freedman-Diaconis rule
  # https://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram
  # Sets that optimal bin width to the default setting in the sliderInput
  output$bin_slider <- renderUI({
    tagged_fish <- taggedfishVar()
    
    bins <- nclass.FD(unlist(tagged_fish[,input$variable]))
    sliderInput("bin_width", "Bin width", min = 1, max = 50, value = bins)
  })
  
  # Output a summary table - using package arsenal
  output$summarytable <- renderTable({
    
    mylabels <- list(Length = "Length (mm)", Weight ="Weight (g)")
    as.data.frame(summary(tableby(StudyID ~ Length + Weight, 
                                  data = taggedfishVar(),
                                  control = tableby.control(digits = 2)), 
                          labelTranslations = mylabels, text = "html"))
    
  }, sanitize.text.function = function(x) x)
  
  # Create plotly plot
  output$plotly_plot <- renderPlotly({
    tagged_fish <- taggedfishVar()
    
    # Based on the user selected plot type, use different ggplot options 
    plot_type <- switch(input$plot_type,
                        "boxplot" 	= 	geom_boxplot(),
                        "histogram" =	geom_histogram(bins = input$bin_width, 
                                                     alpha=0.5,
                                                     position="identity"),
                        "density" 	=	geom_density(alpha=.75)
    )
    
    # Based on the user selected var type, use different labels 
    var_lab <- switch(input$variable,
                      "Length" 	= " (mm)",
                      "Weight" =	" (g)"
    )
    
    # Create the plotly plot differently depending on plot type
    if (input$plot_type == "boxplot") {
      ggplotly(ggplot(tagged_fish, 
                      aes(
                        x = StudyID,
                        y = get(input$variable),
                        fill = StudyID
                      )
      ) + 
        plot_type +
        ylab(paste0(input$variable, var_lab)) +
        scale_fill_brewer(palette="Accent") +
        theme_classic() +
        theme(legend.position="top", axis.text=element_text(size=12))) 
    } else {
      ggplotly(ggplot(tagged_fish, 
                      aes(
                        x = get(input$variable),
                        fill = StudyID
                      )
      ) + 
        plot_type +
        scale_y_continuous(expand = c(0, 0)) +
        xlab(paste0(input$variable, var_lab)) +
        ylab(ifelse(input$plot_type == 'histogram', 'Count', 'Density')) +
        scale_fill_brewer(palette="Accent") +
        theme_classic() +
        theme(legend.position="top", axis.text=element_text(size=12))) 
    }
    
  })
  
  ##%######################################################%##
  #                                                          #
  ####                    Time of Day                     ####
  #                                                          #
  ##%######################################################%##
  
  # Reactive list of GEN depending on StudyID selected
  GEN_list <-  reactive({
    detections <- timeofdayVar()
    
    unique(detections$GEN)
  })
  
  # Updates allowable GEN selection input using GEN_list()
  observe({
    updateSelectInput(session, "time_of_day_GEN",
                      choices = GEN_list()
    )
  })
  
  # Reactive expression that retrieves detections from ERDDAP and formats it for 
  # the time step animation
  timeofdayVar<-  reactive({
    req(input$time_of_day_input)
    
    # Directory of detection files
    files <- list.files("./detections/", full.names = T)
    
    original_studyid <- convert_descriptive_to_studyid(input$time_of_day_input)
    file <- files[str_detect(files, original_studyid)]
    detections <- vroom(file, col_types = cols()) %>% 
      mutate(hour = hour(time_pst))
    
  })
  
  output$time_of_day_plot <- renderPlotly({
    # Retrieve the file name in the detections folder that corresponds to the 
    # StudyID selected by the user, Look at all files in the detections folder, 
    # choose the one that contains the StudyID string in its name
    
    detections <- timeofdayVar()
    
    # Get first detection time of each fish at each GEN
    filtered <- detections  %>% 
      group_by(FishID, GEN) %>% 
      summarise(first_detection = min(hour))
    
    # If user selected to examine time of arrivals by GEN filter the detection 
    # file by the GEN they chose
    if (input$time_of_day_radio == "By General Location") {
      filtered <- filtered %>% 
        filter(GEN == input$time_of_day_GEN)
    }
    
    # Bins data by hour explicitly from 0:23 hours
    # Necessary to do it this way instead of dplyr group_by b/c not all hours 
    # may be represented
    
    hour_freq <- data.frame("hour" = 0:23) %>% 
      left_join(filtered %>% 
                  group_by(first_detection) %>% 
                  count(), by = c("hour" = "first_detection")) %>% 
      # Replace NA with 0
      rename(freq = n) %>% 
      mutate(freq = ifelse(is.na(freq), 0, freq),
             percent_pass = ((freq / sum(freq)) * 100)
      )
    
    # Find Sunset, Sunrise for earliest and latest date using mean lat/lon
    receivers <- detections %>% 
      select(GEN, LAT, LON) %>% 
      distinct()
    
    avg_lat = mean(receivers$LAT)
    avg_lon = mean(receivers$LON)
    
    # Get earliest date and latest dates to compare differences in these times
    earliest_date <- as.Date(min(detections$time_pst))
    latest_date <- as.Date(max(detections$time_pst))
    
    earliest_sr_ss <- getSunlightTimes(earliest_date, 
                                       avg_lat, 
                                       avg_lon, 
                                       keep = c("sunrise", "sunset"), 
                                       tz = "America/Los_Angeles")
    latest_sr_ss <- getSunlightTimes(latest_date, 
                                     avg_lat, 
                                     avg_lon, 
                                     keep = c("sunrise", "sunset"), 
                                     tz = "America/Los_Angeles")
    
    earliest_sunrise <- hour(earliest_sr_ss$sunrise) + 
      (minute(earliest_sr_ss$sunrise)/60) + 
      (second(earliest_sr_ss$sunrise)/3600)
    
    earliest_sunset <- hour(earliest_sr_ss$sunset) + 
      (minute(earliest_sr_ss$sunset)/60) + (second(earliest_sr_ss$sunset)/3600)
    
    latest_sunrise <- hour(latest_sr_ss$sunrise) + 
      (minute(latest_sr_ss$sunrise)/60) + (second(latest_sr_ss$sunrise)/3600)
    
    latest_sunset <- hour(latest_sr_ss$sunset) + 
      (minute(latest_sr_ss$sunset)/60) + (second(latest_sr_ss$sunset)/3600)
    
    # # Convert hour to factor and give it specific levels, so that it will appear 
    # # as 12:23, 0:11 on the x-axis
    # hour_freq$hour <- factor(hour_freq$hour, levels = c(12:23, 0:11))
    
    hour_freq$hour <- factor(hour_freq$hour, levels = c(0:23))
    
    # Create bar plot to show proportion of detections by hour
    p <- ggplot(data = hour_freq, mapping = aes(x = hour, y = percent_pass)) +
      # # Add rectangle representing "nighttime", using the earliest sunrise and earliest sunset values, had to use 
      # # ymax = 999999 because plotly won't take Inf
      # geom_rect(data=hour_freq, aes(NULL,NULL,xmin=earliest_sunrise,xmax=earliest_sunset),
      #           ymin=0,ymax=999999, size=0.5, alpha=0.2) +
      
      # Have to add +1 to geom_rect and geom_vline because x-axis set as factor
      # and is offset by 1 for some reason
      geom_rect(aes(NULL,NULL,xmin=0,xmax=min(c(earliest_sunrise, latest_sunrise)) + 1),
                ymin=0,ymax=999999, size=0.5, alpha=0.2) +
      geom_rect(aes(NULL,NULL,xmin=max(c(earliest_sunset, latest_sunset)) + 1,xmax=25),
                ymin=0,ymax=999999, size=0.5, alpha=0.2) +
      geom_col() +
      # Add lines to represent the latest sunrise and latest sunset values
      geom_vline(xintercept = max(c(earliest_sunrise, latest_sunrise)) + 1, linetype = "dashed", size = 1.25) +
      geom_vline(xintercept = min(c(earliest_sunset, latest_sunset)) + 1, linetype = "dashed", size = 1.25) +
      ylab("% Fish Passed") +
      xlab("Time of day (h)") + 
      scale_y_continuous(expand = c(0, 0)) +
      theme_classic() +
      theme(
        text = element_text(size = 15),
        # Make caption text small and left justify
        plot.caption = element_text(size = 10, hjust = 0),
        plot.margin = margin(l = 20, b = 20)
      )
    
    if (input$time_of_day_radio == "By General Location") {
      # Get the GenRKM value for this GEN being examined, for use in labeling in plot
      genrkm <- detections %>% 
        filter(GEN == input$time_of_day_GEN) %>% 
        select(GenRKM) %>% 
        distinct() %>% 
        unlist()
      
      p <- (p + labs(
        title = paste0("% Fish arrivals at ", input$time_of_day_GEN),
        subtitle = paste0("GenRKM = (", genrkm, ")")))
    } else {
      p <- (p + labs(
        title = "% Fish arrivals"))
    }
    
    ggplotly(p)
  })
  
  # Output the plot caption using renderUI instead of in ggplot because it 
  # will autowidth the same length as the plot
  # Couldn't figure out how to autowidth in ggplot
  output$time_of_day_caption <- renderUI({
    
    detections <- timeofdayVar() %>% 
      mutate(time_pst = ymd_hms(time_pst))
    # Get earliest date and latest dates to compare differences in these times
    earliest_date <- as.Date(min(detections$time_pst))
    latest_date <- as.Date(max(detections$time_pst))
    
    paste0("The histogram displays the frequency of fish detections as a function of hour of day. The gray box represents \nhours of nighttime between the earliest sunrise and latest sunset. The dotted lines represent \nhours of nighttime between the latest sunrise and earliest sunset. The first date of detection was at: ", earliest_date, " and the last date of detection was at: ", latest_date, ".")
  })
  
  
  ##%######################################################%##
  #                                                          #
  ####                      Survival                      ####
  #                                                          #
  ##%######################################################%##
  
  # Synchronize studyID selection between tabs
  # https://stackoverflow.com/questions/44516768/r-shiny-synchronize-filters-on-multiple-tabs
  studyIDSelect <- reactiveVal("")  
  
  output$cumSurvSelect <- renderUI({
    selectInput(inputId = "id1", 
                label = "Select", 
                choices = surv_studyid_descript_name, 
                selected = studyIDSelect())
  })
  
  output$reachSurvSelect <- renderUI({
    selectInput(inputId = "id2", 
                label = "Select", 
                choices = surv_studyid_descript_name, 
                selected = studyIDSelect())
  })
  
  output$deltaSurvSelect <- renderUI({
     selectInput(inputId = "id3", 
                 label = "Select", 
                 choices = deltasurv_studyid_descript_name, 
                 selected = studyIDSelect())
  })
  
  observeEvent(input$id2,{
    studyIDSelect(input$id2)
  })
  
  observeEvent(input$id1,{
    studyIDSelect(input$id1)
  })
  
  observeEvent(input$id3,{
     studyIDSelect(input$id3)
  })
  
  # Reactive expression for cumulative survival input
  cumsurvivalVar<-  reactive({
    req(studyIDSelect())
    
    name <- studyIDSelect() 
    original_studyid <- convert_descriptive_to_studyid(name)
    
    file <- list.files("./data/Survival/Cumulative Survival/", original_studyid, full.names = T)
    df <- vroom(file, col_types = cols())
    
    df <- df %>% 
      mutate_at(c("cum.phi", "cum.phi.se", "LCI", "UCI"), round, digits = 3) %>% 
      mutate(
        RKM = round(RKM, digits = 2),
        id = seq.int(nrow(df))
      ) %>% 
      dplyr::rename(
        Survival = cum.phi,
        SE = cum.phi.se,
        Count = count
        
      )
  })
  
  
  output$cumsurvival_text <- renderUI({
    sample_size <- cumsurvivalVar() %>% 
      select("Count") %>% 
      slice(1)
    
    tags$h4(paste0("Number of fish tagged: ", sample_size))
  })
  
  output$plotly_survival_output <- renderPlotly({
    
    df <- cumsurvivalVar()
    
    # If survival estimates contain the column 'release'
    if (any(str_detect(colnames(df), "release"))) {
      # Set the order of release groups
      df$release <- factor(df$release, unique(df$release))
       
      release_colors <- c("#E69F00","#56B4E9","#009E73", "#0072B2", "#D55E00", "#CC79A7")
      release_colors <- release_colors[1:length(unique(df$release))]
       
      p <- df %>% 
        ggplot(mapping = aes(x = RKM, y = Survival, group = release)) +
        geom_point(position = position_dodge(width = 1), 
                   aes(color = release,
                       text = paste('</br> Survival: ', Survival,
                                    '</br> LCI: ', LCI,
                                    '</br> UCI: ', UCI,
                                    '</br> GEN: ', GEN,
                                    '</br> RKM: ', RKM,
                                    '</br> Count: ', Count,
                                    '</br> Release: ', release
                       ))) +
        geom_line(aes(color = release)) +
        geom_errorbar(mapping = aes(x = RKM, ymin = LCI, ymax = UCI, 
                                    color = release),  width = .1,
                      position = position_dodge(1)) +
        scale_color_manual(values=release_colors) +
        labs(
          x = "RKM",
          y = "Survival",
          color = "Release"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5)
        ) 
      
      ggplotly(p, tooltip = "text") %>% 
        layout(xaxis = list(autorange = "reversed"))
    } else {
      df %>% 
        plot_ly(x = ~RKM, y = ~Survival, type = 'scatter', mode = 'lines+markers',
                hoverinfo = 'text',
                text = ~paste('</br> Survival: ', Survival,
                              '</br> LCI: ', LCI,
                              '</br> UCI: ', UCI,
                              '</br> GEN: ', GEN,
                              '</br> GenRKM: ', RKM,
                              '</br> Raw number of fish to site: ', Count
                ),
                # https://rpubs.com/chelsea/68601
                error_y = ~list(
                  type = "data", 
                  symmetric = FALSE, 
                  arrayminus = Survival - LCI,
                  array = UCI - Survival)) %>% 
        layout(title = paste0("Cumulative Survival for ", studyIDSelect())) %>%
        layout(xaxis = list(autorange = "reversed",
                            showline = FALSE,
                            zeroline = FALSE)) %>% 
        layout(yaxis = list(showline = FALSE,
                            zeroline = FALSE))
    }
  })
  
  output$survival_map2 <- renderLeaflet({
    
    df <- cumsurvivalVar()
    
    df %>% 
      slice(1:(max(df$reach_num) + 1)) %>% 
      # filter(reach_num != 0) %>% 
      leaflet() %>% 
      addProviderTiles(providers$Stamen.Terrain, group = "Stamen Terrain",
                       options = providerTileOptions(noWrap = TRUE)) %>% 
      addTiles(group = "Open Street Map") %>% 
      addProviderTiles(providers$Esri.NatGeoWorldMap, group = "Esri Nat Geo",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      addMarkers(
        ~GenLon, ~GenLat, 
        layerId = as.character(df$id),
        popup = paste0( "<b>Receiver Location: </b>", 
                        df$GEN,
                        "</br>",
                        "<b>RKM: </b>", 
                        df$RKM,
                        "</br>"
        ),
        label = ~GEN
      ) %>% 
      # Give control for selecting basemaps
      addLayersControl(
        baseGroups = c("Stamen Terrain", "Open Street Map", "Esri Nat Geo", 
                       "Esri World Imagery"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  output$cumSurvDT <- renderDataTable({
    df <- cumsurvivalVar()
    
    # If survival estimates contain the column 'release'
    if (any(str_detect(colnames(df), "release"))) {
      
      df <- cumsurvivalVar() %>% 
        select('Receiver Location' =  GEN, RKM, Region, Survival, SE, LCI, UCI, 
               Count, release, id)
      # manual method here: https://stackoverflow.com/questions/60399441/tidyrpivot-wider-reorder-column-names-grouping-by-name-from
      prefixes <- unique(df$release)
      
      # Widen data so that numbers for each release represented in a single row
      df <- df %>%
        select(-id) %>% 
        pivot_wider(
          names_from = release,
          values_from = c("Survival", "SE", "LCI", "UCI", 
                          "Count"),
          names_glue = "{release} {.value}"
        ) %>%
        mutate(id = row_number())
      
      names_to_order <- map(prefixes, ~ names(df)[grep(paste0(.x, " "), 
                                                       names(df))]) %>% unlist
      names_id <- setdiff(names(df), names_to_order)
      
      df <- df %>%
        select(names_id, names_to_order) %>% 
        select(-id, id) 
      
      
      datatable(df, selection = "single", 
                options=list(stateSave = TRUE,
                             dom = 'Bfrtip',
                             scrollX = TRUE, # allows horizontal scrolling 
                             #fixedHeader = TRUE, # freezes header,
                             scrollY = 500,
                             rownames = FALSE),
                fillContainer = T)
    } else{
      dat <- cumsurvivalVar() %>% 
        select('Receiver Location' =  GEN, RKM, Region, Survival, LCI, UCI, Count, id)
      
      datatable(dat, selection = "single", 
                options=list(stateSave = TRUE,
                             dom = 'Bfrtip',
                             rownames = FALSE))
    }
  })
  
  # Create download button from cumulative file
  output$cumSurvDownload <- downloadHandler(
    
    # Use the file name of the cumulative survival csv 
    filename = function() {
      name <- studyIDSelect() 
      original_studyid <- convert_descriptive_to_studyid(name)
      
      file <- list.files("./data/Survival/Cumulative Survival/", original_studyid)
    },
    content = function(con) {
      name <- studyIDSelect() 
      original_studyid <- convert_descriptive_to_studyid(name)
      
      file <- list.files("./data/Survival/Cumulative Survival/", original_studyid, 
                         full.names = T)
      
      file.copy(file, con)
    }
  )
  
  # to keep track of previously selected row
  prev_row <- reactiveVal()
  
  # new icon style for selected rows
  my_icon = makeAwesomeIcon(icon = 'flag', markerColor = 'red', iconColor = 'white')
  
  # new icon style for release site
  my_icon2 = makeAwesomeIcon(icon = 'star', markerColor = 'purple', iconColor = 'white')
  
  # https://stackoverflow.com/questions/48781380/shiny-how-to-highlight-an-object-on-a-leaflet-map-when-selecting-a-record-in-a
  observeEvent(input$cumSurvDT_rows_selected, {
    row_selected <-  cumsurvivalVar()[input$cumSurvDT_rows_selected,]
    proxy <- leafletProxy('survival_map2')
    proxy %>%
      addAwesomeMarkers(
        popup = paste0( "<b>Receiver Location: </b>", 
                        row_selected$GEN,
                        "</br>",
                        "<b>RKM: </b>", 
                        row_selected$RKM,
                        "</br>"
        ),
        label = row_selected$GEN,
        layerId = as.character(row_selected$id),
        lng=row_selected$GenLon,
        lat=row_selected$GenLat,
        icon = my_icon)
    
    # Reset previously selected marker
    if(!is.null(prev_row()))
    {
      proxy %>%
        addMarkers(popup=as.character(prev_row()$GEN), 
                   layerId = as.character(prev_row()$id),
                   lng=prev_row()$GenLon, 
                   lat=prev_row()$GenLat)
    }
    # set new value to reactiveVal 
    prev_row(row_selected)
  })
  
  observeEvent(input$survival_map2_marker_click, {
    clickId <- input$survival_map2_marker_click$id
    dataTableProxy("cumSurvDT") %>%
      selectRows(which(cumsurvivalVar()$id == clickId)) %>%
      selectPage(which(input$cumSurvDT_rows_all == clickId) %/% input$cumSurvDT_state$length + 1)
  })
  
  # Reset markers when switching from table view to plot view
  observeEvent(input$cumsurvival_radio, {
    # Reset previously selected marker
    if(!is.null(prev_row()) & input$cumsurvival_radio == "Plot")
    {
      proxy <- leafletProxy('survival_map2')
      proxy %>%
        addMarkers(popup=as.character(prev_row()$GEN), 
                   layerId = as.character(prev_row()$id),
                   lng=prev_row()$GenLon, 
                   lat=prev_row()$GenLat)
    }
  })
  
  # Reactive expression for reach survival input
  reachSurvVar<-  reactive({
    req(studyIDSelect())
    
    name <- studyIDSelect()
    original_studyid <- convert_descriptive_to_studyid(name)
    
    file <- list.files("./data/Survival/Reach Survival Per 10km", original_studyid, full.names = T)
    df <- vroom(file, col_types = cols())
    
    # If studyID contains multi rel loc, rearrange data to wide format
    if(length(unique(df$release)) > 1) {
      
      df <- df %>% 
        mutate_at(c("estimate", "se", "lcl", "ucl"), round, digits = 3) %>% 
        mutate_at(c("rkm_start", "rkm_end"), round, digits = 2) %>% 
        mutate(id = seq.int(nrow(df))) %>%
        rename(
          Survival = estimate,
          SE = se,
          LCI = lcl,
          UCI = ucl
        ) %>% 
        filter(reach_end != "GoldenGateW")
      
    } else {
      df <- df %>% 
        mutate_at(c("estimate", "se", "lcl", "ucl"), round, digits = 3) %>% 
        mutate_at(c("rkm_start", "rkm_end"), round, digits = 2) %>% 
        mutate(id = seq.int(nrow(df))) %>%
        rename(
          Survival = estimate,
          SE = se,
          LCI = lcl,
          UCI = ucl
        ) %>% 
        filter(reach_end != "GoldenGateW")
    }
  })
  
  output$reachSurv_text <- renderUI({
    sample_size <- cumsurvivalVar() %>% 
      select("Count") %>% 
      slice(1)
    
    tags$h4(paste0("Number of fish tagged: ", sample_size))
  })
  
  output$reachSurvMap <- renderLeaflet({
    
    df <- reachSurvVar()
    
    df %>% 
      # Do this b/c if multiple releases there will be duplictate rows only 
      # need one set
      dplyr::slice(1:max(df$reach_num)) %>% 
      leaflet() %>% 
      addProviderTiles(providers$Stamen.Terrain, group = "Stamen Terrain",
                       options = providerTileOptions(noWrap = TRUE)) %>% 
      addTiles(group = "Open Street Map") %>% 
      addProviderTiles(providers$Esri.NatGeoWorldMap, group = "Esri Nat Geo",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery",
                       options = providerTileOptions(noWrap = TRUE)) %>%
      addMarkers(
        ~GenLon_end, ~GenLat_end, 
        layerId = df$id,
        popup = paste0( "<b>Reach: </b>", 
                        df$reach_start, " to ", df$reach_end,
                        "</br>",
                        "<b>RKM: </b>", 
                        df$rkm_start, " to ", df$rkm_end,
                        "</br>"
        ),
        label = ~Reach
      ) %>% 
      # Add the release site as a unique marker
      addAwesomeMarkers(
        data = df %>% filter(id == 1),
        lng = ~GenLon_start,
        lat = ~GenLat_start,
        popup = paste0( "<b>Release Location: </b>"
                        , df$reach_start[1]
                        , "<br>"
                        , "<b># Fish Tagged: </b>"
                        , TaggedFish %>% 
                          filter(study_id == df$StudyID[1]) %>% 
                          count() %>% 
                          unlist()
        ),
        label = ~reach_start,
        icon = my_icon2
      ) %>%
      addLayersControl(
        baseGroups = c("Stamen Terrain", "Open Street Map", "Esri Nat Geo", 
                       "Esri World Imagery"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  output$reachSurvPlotly <- renderPlotly({
    
    name <- studyIDSelect()
    
    df <- reachSurvVar()
    
    # Plot differently if there is a single release vs multi
    # Check if df has the colname 'release
    if(any(sapply(colnames(df), str_detect, "release"))) {
      # Set the order of release groups
      df$release <- factor(df$release, unique(df$release))
      
      release_colors <- c("#E69F00","#56B4E9","#009E73", "#0072B2", "#D55E00", "#CC79A7")
      release_colors <- release_colors[1:length(unique(df$release))]
      
      p <- df %>% 
        mutate(reach_end = factor(reach_end, levels = unique(reach_end))) %>%
        mutate(reach_end = reorder(reach_end, -rkm_end)) %>%
        ggplot(mapping = aes(x = reach_end, y = Survival, group = release)) +
        geom_point(position = position_dodge(width = 0.5), 
                   aes(color = release,
                       text = paste('</br> Survival: ', Survival,
                                    '</br> LCI: ', LCI,
                                    '</br> UCI: ', UCI,
                                    '</br> Reach: ', paste0(reach_start, " to ", reach_end),
                                    '</br> RKM: ', paste0(rkm_start, " to " , rkm_end),
                                    '</br> Count Start: ', count_at_start,
                                    '</br> Count End: ', count_at_end
                       ))) +
        geom_errorbar(mapping = aes(x = reach_end, ymin = LCI, ymax = UCI, 
                                    color = release),  width = .1,
                      position = position_dodge(.5)) +
        scale_color_manual(values=release_colors) +
        labs(
          title = paste0("Reach Survival for ", name),
          x = "Receiver Location",
          y = "Survival per 10km",
          color = "Release"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)
        )
      
      ggplotly(p, tooltip = "text")
    } else {
      p <- df %>% 
        mutate(reach_end = factor(reach_end, levels = unique(reach_end))) %>% 
        ggplot(mapping = aes(x = reach_end, y = Survival)) +
        geom_point(position = position_dodge(width = 0.5), 
                   aes(text = paste('</br> Survival: ', Survival,
                                    '</br> LCI: ', LCI,
                                    '</br> UCI: ', UCI,
                                    '</br> Reach: ', paste0(reach_start, " to ", reach_end),
                                    '</br> RKM: ', paste0(rkm_start, " to " , rkm_end),
                                    '</br> Count Start: ', count_at_start,
                                    '</br> Count End: ', count_at_end
                   )),
                   color = "#007EFF") +
        geom_errorbar(mapping = aes(x = reach_end, ymin = LCI, ymax = UCI),  
                      width = .1, position = position_dodge(.5),
                      color = "#007EFF") +
        labs(
          title = paste0("Reach Survival for ", name),
          x = "Receiver Location",
          y = "Survival per 10km"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)
        )
      
      ggplotly(p, tooltip = "text")
    }
    
  })
  
  output$reachSurvDT <- renderDataTable({
    df <- reachSurvVar()
    
    # If multiple release
    if(length(unique(df$release)) > 1) {
      # Pivot_wider does not natively allow ordering columns by names
      # manual method here: https://stackoverflow.com/questions/60399441/tidyrpivot-wider-reorder-column-names-grouping-by-name-from
      prefixes <- unique(df$release)
      
      # Widen data so that numbers for each release represented in a single row
      df <- df %>%
        select(-id) %>% 
        dplyr::rename(
          'Count Start' = 'count_at_start',
          'Count End' = 'count_at_end'
        ) %>% 
        pivot_wider(
          names_from = release,
          values_from = c("Survival", "SE", "LCI", "UCI", 
                          "Count Start", "Count End"),
          names_glue = "{release} {.value}"
        ) %>%
        mutate(id = row_number())
      
      names_to_order <- map(prefixes, ~ names(df)[grep(paste0(.x, " "), names(df))]) %>% unlist
      names_id <- setdiff(names(df), names_to_order)
      
      df <- df %>%
        select(names_id, names_to_order) %>% 
        select(-id, id)
      
      
      df %>% 
        select(-c(reach_start, reach_end, rkm_start, rkm_end, Region, reach_num,
                  GenLon_start, GenLon_end, GenLat_start, GenLat_end)) %>% 
        datatable(selection = "single", 
                  options=list(stateSave = TRUE,
                               dom = 'Bfrtip',
                               scrollX = TRUE, # allows horizontal scrolling 
                               #fixedHeader = TRUE, # freezes header,
                               scrollY = 500,
                               rownames = FALSE),
                  # This is allows x scrolling and freezer headers to work simultaneously
                  # https://github.com/rstudio/DT/issues/643
                  fillContainer = T
        )
      
    } else {
      df <- reachSurvVar() %>% 
        select('Reach Start' = reach_start, 'Reach End' = reach_end, 
               'RKM Start' = rkm_start, 'RKM End' = rkm_end, Region, 
               Survival, LCI, UCI, 'Count Start' = count_at_start,
               'Count End' = count_at_end, id)
      
      datatable(df, selection = "single", #extensions = 'Buttons',
                options=list(stateSave = TRUE,
                             dom = 'Bfrtip',
                             rownames = FALSE))
    }
    
    
  })
  
  # Create download button from cumulative file
  output$reachSurvDownload <- downloadHandler(
    
    # Use the file name of the cumulative survival csv 
    filename = function() {
      name <- studyIDSelect() 
      original_studyid <- convert_descriptive_to_studyid(name)
      
      
      file <- list.files("./data/Survival/Reach Survival Per 10km/", original_studyid)
    },
    content = function(con) {
      name <- studyIDSelect() 
      original_studyid <- convert_descriptive_to_studyid(name)
      
      file <- list.files("./data/Survival/Reach Survival Per 10km/", original_studyid, 
                         full.names = T)
      
      file.copy(file, con)
    }
  )
  
  # to keep track of previously selected row
  prev_row2 <- reactiveVal()
  
  observeEvent(input$reachSurvDT_rows_selected, {
    row_selected <-  reachSurvVar()[input$reachSurvDT_rows_selected,]
    proxy <- leafletProxy('reachSurvMap')
    proxy %>%
      addAwesomeMarkers(
        popup = paste0( "<b>Reach: </b>", 
                        row_selected$reach_start, " to ", row_selected$reach_end,
                        "</br>",
                        "<b>RKM: </b>", 
                        row_selected$rkm_start, " to ", row_selected$rkm_end,
                        "</br>"
        ),
        label = row_selected$Reach,
        layerId = as.character(row_selected$id),
        lng=row_selected$GenLon_end,
        lat=row_selected$GenLat_end,
        icon = my_icon)
    
    # Reset previously selected marker
    if(!is.null(prev_row2()))
    {
      proxy %>%
        addMarkers(popup = paste0( "<b>Reach: </b>", 
                                   prev_row2()$reach_start, " to ", prev_row2()$reach_end,
                                   "</br>",
                                   "<b>RKM: </b>", 
                                   prev_row2()$rkm_start, " to ", prev_row2()$rkm_end,
                                   "</br>"
        ), 
        layerId = as.character(prev_row2()$id),
        lng=prev_row2()$GenLon_end, 
        lat=prev_row2()$GenLat_end)
    }
    # set new value to reactiveVal 
    prev_row2(row_selected)
  })
  
  observeEvent(input$reachSurvMap_marker_click, {
    clickId <- input$reachSurvMap_marker_click$id
    dataTableProxy("reachSurvDT") %>%
      selectRows(which(reachSurvVar()$id == clickId)) %>%
      selectPage(which(input$reachSurvDT_rows_all == clickId) %/% input$reachSurvDT_state$length + 1)
  })
  
  # Reset markers when switching from table view to plot view
  observeEvent(input$reachSurvRadio, {
    # Reset previously selected marker
    if(!is.null(prev_row2()) & input$reachSurvRadio == "Plot")
    {
      proxy <- leafletProxy('reachSurvMap')
      proxy %>%
        addMarkers(popup = paste0( 
          "<b>Reach: </b>", 
          paste0(prev_row2()$reach_start, ' to ', prev_row2()$reach_end),
          "</br>",
          "<b>RKM: </b>", 
          paste0(prev_row2()$rkm_start, ' to ', prev_row2()$rkm_end),
          "</br>",
          prev_row2()$count_at_end), 
          layerId = as.character(prev_row2()$id),
          lng=prev_row2()$GenLon_end, 
          lat=prev_row2()$GenLat_end)
    }
  })
  
  ##################### SERVER-SIDE FOR DELTA ROUTE-SPECIFIC SURVIVAL #####################
  # Here is a reactive expression that will grab the Rdata file based on user selection from the drop down menu
  
  # Reactive expression for cumulative survival input
  deltasurvivalVar<-  reactive({
     req(studyIDSelect())
     
     name <- studyIDSelect() 
     original_studyid <- convert_descriptive_to_studyid(name)
     
     file <- list.files("./data/Survival/Delta route specific survival/", original_studyid, full.names = T)
     
     summarydata <- vroom(file, col_types = cols())
     
     
  })
  
  ## leaflet map
  
  output$deltasurv_map <- renderLeaflet({
     
     geo_col <- "#8F009A"
     sac_col <- "#0071B7"
     stm_col <- "#180F9B"
     sut_col  <- "#F99B00"
     
     
     sac <- st_read("./data/Survival/Delta route specific survival/KMLs/Sacramento.kml") %>%
        st_coordinates() %>%
        subset(, select = -c(Z, L1, L2))
     
     sut <- st_read("./data/Survival/Delta route specific survival/KMLs/Sutter full route.kml") %>%
        st_coordinates() %>%
        subset(, select = -c(Z, L1, L2))
     
     sutb <- st_read("./data/Survival/Delta route specific survival/KMLs/Minor via Sutter.kml") %>%
        st_coordinates() %>%
        subset(, select = -c(Z, L1, L2))
     
     stm <- st_read("./data/Survival/Delta route specific survival/KMLs/Steamboat full route.kml") %>%
        st_coordinates() %>%
        subset(, select = -c(Z, L1, L2))
     
     geo <- st_read("./data/Survival/Delta route specific survival/KMLs/Georgiana full route.kml") %>%
        st_coordinates() %>%
        subset(, select = -c(Z, L1, L2))
     
     int <- st_read("./data/Survival/Delta route specific survival/KMLs/Interior Delta.kml") %>%
        st_coordinates() %>%
        subset(, select = -c(Z, L1, L2))
     
     ## Setup USGS Leaflet Map
     grp <- c("USGS Topo", "USGS Imagery Only", "USGS Imagery Topo",
              "USGS Shaded Relief", "Hydrography")
     
     att <- paste0("<a href='https://www.usgs.gov/'>",
                   "U.S. Geological Survey</a> | ",
                   "<a href='https://www.usgs.gov/laws/policies_notices.html'>",
                   "Policies</a>")
     
     GetURL <- function(service, host = "basemap.nationalmap.gov") {
        sprintf("https://%s/arcgis/services/%s/MapServer/WmsServer", host, service)
     }
     
     # Code below creates marker locations and popup labels
     lon.sac <- c(-121.50862,-121.51629,-121.56578,-121.56178,
                  -121.52522,-121.53536,
                  -122.12338,-122.12695)
     
     lat.sac <- c(38.58044,38.57117,38.296288,38.29145,
                  38.23956,38.23922,
                  38.04337,38.04494)
     
     label.sac <- c("Sacramento River at Tower Bridge",
                    "Sacramento River at I80-50 Bridge",
                    "Sacramento River below Steamboat Slough US",
                    "Sacramento River below Steamboat Slough DS",
                    "Sacramento River below Georgiana Slough US",
                    "Sacramento River below Georgiana Slough DS",
                    "Sacramento River at Benicia E", 
                    "Sacramento River at Benicia W")
     
     
     lon.sut <- c(-121.58522,-121.58583)
     
     lat.sut <- c(38.330406,38.327988)
     
     label.sut <- c("Sutter Slough Entrance US", 
                    "Sutter Slough Entrance DS")
     
     lon.stm <- c(-121.58695,-121.58916)
     
     lat.stm <- c(38.285114,38.28016)
     
     label.stm <- c("Steamboat Slough Entrance US", 
                    "Steamboat Slough Entrance DS")
     
     lon.geo <- c(-121.5196, -121.52287)
     
     lat.geo <- c(38.23417,38.23112)
     
     label.geo <- c("Georgiana Slough Entrance US", 
                    "Georgiana Slough Entrance DS")
     
     
     #parameterize the icons for each group
     icon.sac <- makeAwesomeIcon(icon= 'dot', markerColor = 'lightblue', iconColor = 'lightblue')
     icon.sut <- makeAwesomeIcon(icon= 'dot', markerColor = 'orange', iconColor = 'orange')
     icon.stm <- makeAwesomeIcon(icon= 'dot', markerColor = 'darkblue', iconColor = 'darkblue')
     icon.geo <- makeAwesomeIcon(icon= 'dot', markerColor = 'purple', iconColor = 'purple')
     
     ############################################################################################
     
     #Create map
     leaflet() %>%
        addWMSTiles(GetURL("USGSShadedReliefOnly"), attribution = att, layers = "0") %>%
        addWMSTiles(GetURL("USGSHydroCached"), options = WMSTileOptions(format = "image/png", transparent = TRUE), layers = "0") %>%
        #overlay Groups
        addAwesomeMarkers(lng=lon.sac, lat=lat.sac, icon=icon.sac, label = label.sac, group = "Telemetry stations") %>%
        addAwesomeMarkers(lng=lon.sut, lat=lat.sut, icon=icon.sut, label = label.sut, group = "Telemetry stations") %>%
        addAwesomeMarkers(lng=lon.stm, lat=lat.stm, icon=icon.stm, label = label.stm, group = "Telemetry stations") %>%
        addAwesomeMarkers(lng=lon.geo, lat=lat.geo, icon=icon.geo, label = label.geo, group = "Telemetry stations") %>%
        addPolygons(data = sut, stroke = TRUE, weight = 2, color = sut_col, group = "Sutter Slough", label = "Sutter Slough") %>%
        addPolygons(data = sutb, stroke = TRUE, weight = 2, color = sut_col, group = "Sutter Slough", label = "Miner Slough") %>%
        addPolygons(data = stm, stroke = TRUE, weight = 2, color = stm_col, group = "Steamboat Slough", label = "Steamboat Slough") %>%
        addPolygons(data = geo, stroke = TRUE, weight = 2, color = geo_col, group = "Georgiana Slough", label = "Interior Delta via Georgiana Slough") %>%
        addPolygons(data = int, stroke = FALSE, weight = 2, color = geo_col, group = "Georgiana Slough", label = "Interior Delta") %>%
        addPolygons(data = sac, stroke = TRUE, weight = 2, color = sac_col, group = "Sacramento River", label = "Sacramento River") %>%
        
        #Layers control
        addLayersControl(
           overlayGroups = c("Telemetry stations", ##For markers
                             "Sacramento River",
                             "Sutter Slough","Steamboat Slough",
                             "Georgiana Slough"
           ),
           options = layersControlOptions(collapsed = TRUE)
        )
     
  })
  
  
  
  
  ## diagram
  
  output$deltasurv_diagram <- renderGrViz({
     
     summarydata <- deltasurvivalVar()
     
     
     # prep data for diagram
     
     release_nodelabel <- c(paste0("Release"))
     
     p_nodelabel <- c(paste0("Sacramento River\nat Sacramento\np=",summarydata$est[summarydata$parameter=="P_3"], " (",summarydata$l95[summarydata$parameter=="P_3"], "-",summarydata$u95[summarydata$parameter=="P_3"], ")"),
                      paste0("Sutter Slough\np=",summarydata$est[summarydata$parameter=="P_14"], " (",summarydata$l95[summarydata$parameter=="P_14"], "-",summarydata$u95[summarydata$parameter=="P_14"], ")"),
                      paste0("Steamboat Slough\np=",summarydata$est[summarydata$parameter=="P_15"], " (",summarydata$l95[summarydata$parameter=="P_15"], "-",summarydata$u95[summarydata$parameter=="P_15"], ")"),
                      paste0("Sacramento River\nbelow Steamboat\np=",summarydata$est[summarydata$parameter=="P_8"], " (",summarydata$l95[summarydata$parameter=="P_8"], "-",summarydata$u95[summarydata$parameter=="P_8"], ")"),
                      paste0("Georgiana Slough\np=",summarydata$est[summarydata$parameter=="P_21"], " (",summarydata$l95[summarydata$parameter=="P_21"], "-",summarydata$u95[summarydata$parameter=="P_21"], ")"),
                      paste0("Sacramento River\nbelow Georgiana\np=",summarydata$est[summarydata$parameter=="P_22"], " (",summarydata$l95[summarydata$parameter=="P_22"], "-",summarydata$u95[summarydata$parameter=="P_22"], ")"),
                      paste0("Sacramento River\nat Benicia\np=",summarydata$est[summarydata$parameter=="P_30"], " (",summarydata$l95[summarydata$parameter=="P_30"], "-",summarydata$u95[summarydata$parameter=="P_30"], ")"))
     
     r_edgelabel <- c(paste0("r=",summarydata$est[summarydata$parameter=="r1_sut"], " (",summarydata$l95[summarydata$parameter=="r1_sut"], "-",summarydata$u95[summarydata$parameter=="r1_sut"], ")"),
                      paste0("r=",summarydata$est[summarydata$parameter=="r1_ste"], " (",summarydata$l95[summarydata$parameter=="r1_ste"], "-",summarydata$u95[summarydata$parameter=="r1_ste"], ")"),
                      paste0("r=",summarydata$est[summarydata$parameter=="r1_sac"], " (",summarydata$l95[summarydata$parameter=="r1_sac"], "-",summarydata$u95[summarydata$parameter=="r1_sac"], ")"),
                      paste0("r=",summarydata$est[summarydata$parameter=="r2_sac"], " (",summarydata$l95[summarydata$parameter=="r2_sac"], "-",summarydata$u95[summarydata$parameter=="r2_sac"], ")"),
                      paste0("r=",summarydata$est[summarydata$parameter=="r2_geo"], " (",summarydata$l95[summarydata$parameter=="r2_geo"], "-",summarydata$u95[summarydata$parameter=="r2_geo"], ")"))
     
     #release_edgelabel <- c(paste0("S=",summarydata$est[summarydata$parameter=="S_2"]))
     
     S_edgelabel <- c("", #paste0("S=",summarydata$est[summarydata$parameter=="S_2"], " (",summarydata$l95[summarydata$parameter=="S_2"], "-",summarydata$u95[summarydata$parameter=="S_2"], ")"),
                      paste0("S=",summarydata$est[summarydata$parameter=="S_4"], " (",summarydata$l95[summarydata$parameter=="S_4"], "-",summarydata$u95[summarydata$parameter=="S_4"], ")"),
                      paste0("S=",summarydata$est[summarydata$parameter=="S_17"], " (",summarydata$l95[summarydata$parameter=="S_17"], "-",summarydata$u95[summarydata$parameter=="S_17"], ")"),
                      paste0("S=",summarydata$est[summarydata$parameter=="S_16"], " (",summarydata$l95[summarydata$parameter=="S_16"], "-",summarydata$u95[summarydata$parameter=="S_16"], ")"),
                      paste0("S=",summarydata$est[summarydata$parameter=="S_12"], " (",summarydata$l95[summarydata$parameter=="S_12"], "-",summarydata$u95[summarydata$parameter=="S_12"], ")"),
                      paste0("S=",summarydata$est[summarydata$parameter=="S_23"], " (",summarydata$l95[summarydata$parameter=="S_23"], "-",summarydata$u95[summarydata$parameter=="S_23"], ")"),
                      paste0("S=",summarydata$est[summarydata$parameter=="S_24"], " (",summarydata$l95[summarydata$parameter=="S_24"], "-",summarydata$u95[summarydata$parameter=="S_24"], ")"))
     
     
     ###### Generate Diagram
     # create nodes
     # release site
     nodes_1 <-
        create_node_df(
           n = 1,
           label = release_nodelabel,
           color = "black",
           shape = "oval",
           fixedsize=FALSE,
           fillcolor="white",
           fontcolor = "black",
           penwidth=2,
           style="bold",
           fontsize=18
        )
     # detection sites
     nodes_2 <-
        create_node_df(
           n = 7,
           label = p_nodelabel,
           type = "lower",
           color = "black",
           shape = "rectangle",
           fixedsize=FALSE,
           fillcolor="white",
           fontcolor = "black",
           penwidth=2,
           style="bold",
           fontsize=18
        )
     # junctions and invisible nodes for spacing
     nodes_3 <-
        create_node_df(
           n = 6,
           label = c("","","","","",""),
           shape = "circle",
           fixedsize = TRUE,
           width = 0.07,
           color = c("lightblue3","lightblue3","orange","blue4","lightblue3","purple"),
           fillcolor=c("lightblue3","lightblue3","orange","blue4","lightblue3","purple"),
           fontsize = 0.1
        )
     # legend
     nodes_4 <-
        create_node_df(
           n = 1,
           label = c("S = survival probability (95% credible intervals)\nr = route probability (95% credible intervals)\np = detection probability (95% credible intervals)"),
           color = "white",
           shape = "rectangle",
           fixedsize=FALSE,
           fillcolor="white",
           fontcolor = "black",
           fontsize=18
        )
     # S labels
     nodes_5 <-
        create_node_df(
           n = 7,
           label = S_edgelabel,
           color = "white",
           shape = "rectangle",
           fixedsize=TRUE,
           width = 0.1,
           fillcolor="white",
           fontcolor = "black",
           fontsize=18
        )
     # r labels
     nodes_6 <-
        create_node_df(
           n = 5,
           label = r_edgelabel,
           color = "white",
           shape = "rectangle",
           fixedsize=TRUE,
           width = 0.1,
           fillcolor="white",
           fontcolor = "black",
           fontsize=18
        )
     # blank node to set width of diagram so labels are not cut off
     nodes_7 <-
        create_node_df(
           n = 1,
           label = "",
           color = "white",
           shape = "rectangle",
           fixedsize=TRUE,
           width = 0.1,
           fillcolor="white",
           fontcolor = "black",
           fontsize=18
        )
     # blank nodes for route lines in legend
     nodes_8 <-
        create_node_df(
           n = 8,
           label = "",
           color = "white",
           shape = "rectangle",
           fixedsize=TRUE,
           width = 0.1,
           fillcolor="white",
           fontcolor = "black",
           fontsize=18
        )
     # legend - route text
     nodes_9 <-
        create_node_df(
           n = 1,
           label = c("\nSacramento River\nSutter Slough\nSteamboat Slough\n Georgiana Slough"),
           color = "white",
           shape = "rectangle",
           fixedsize=FALSE,
           fillcolor="white",
           fontcolor = "black",
           fontsize=18
        )
     # combine nodes
     all_nodes <- combine_ndfs(nodes_1, nodes_2, nodes_3, nodes_4, nodes_5, nodes_6, nodes_7, nodes_8, nodes_9)
     # create edges
     # release survival edge
     edges_0 <-
        create_edge_df(
           from = c(1),
           to = c(2),
           #label = release_edgelabel,
           fontsize=20,
           arrowsize=1,
           penwidth=5,
           color="lightblue3",
           arrowhead="normal",
           tooltip=" "
        )
     # survival edges
     edges_1 <-
        create_edge_df(
           from = c(2,3,4,5,6,7),
           to = c(9,11,12,10,14,13),
           #label = S_edgelabel,
           fontsize=20,
           arrowsize=0.2,
           penwidth=5,
           color=c("lightblue3","orange","blue4","lightblue3","purple","lightblue3"),
           arrowhead="none"
        )
     # routing edges
     edges_2 <-
        create_edge_df(
           from = c(9,9,9,10,10),
           to =    c(3,4,5,6,7),
           #label = r_edgelabel,
           fontsize=20,
           arrowsize=1,
           penwidth=5,
           color=c("orange","blue4","lightblue3","purple","lightblue3"),
           arrowhead="normal"
        )
     # spacing edges
     edges_3 <-
        create_edge_df(
           from = c(11,12,13,14),
           to =    c(8,8,8,8),
           label=c("","","",""),
           fontsize=20,
           arrowsize=1,
           penwidth=5,
           color=c("orange","blue4","lightblue3","purple"),
           arrowhead="normal"
        )
     # colored lines in legend
     edges_4 <-
        create_edge_df(
           from = c(29,30,31,32),
           to = c(33,34,35,36),
           fontsize=20,
           arrowsize=0.2,
           penwidth=5,
           color=c("lightblue3","orange","blue4","purple"),
           arrowhead="none"
        )
     # combine edges
     all_edges <- combine_edfs(edges_0, edges_1, edges_2, edges_3, edges_4)
     
     # create graph
     graph <- create_graph(nodes_df=all_nodes, edges_df=all_edges)
     # set node positions
     graph <-
        graph %>%
        # Release
        set_node_position(
           node = 1,
           x = 5, y = 11.5) %>% 
        # Sacramento
        set_node_position(
           node = 2,
           x = 5, y = 9) %>%
        # Below Steamboat
        set_node_position(
           node = 5,
           x = 10, y = 6) %>% 
        # Sutter
        set_node_position(
           node = 3,
           x = -2, y = 6) %>% 
        # Steamboat
        set_node_position(
           node = 4,
           x = 3, y = 6) %>% 
        # Below Georgiana
        set_node_position(
           node = 7,
           x = 8, y = 3) %>% 
        # Georgiana
        set_node_position(
           node = 6,
           x = 12, y = 3) %>% 
        # Benicia
        set_node_position(
           node = 8,
           x = 5, y = 0) %>% 
        # junction 1 (Sutter/Steamboat)
        set_node_position(
           node = 9,
           x = 5, y = 7) %>%
        # junction 2 (Georgiana)
        set_node_position(
           node = 10,
           x = 10, y = 4) %>%  
        # Sutter spacing
        set_node_position(
           node = 11,
           x = -2, y = 1) %>%   
        # Steamboat spacing
        set_node_position(
           node = 12,
           x = 3, y = 1) %>%  
        # Below Georgiana spacing
        set_node_position(
           node = 13,
           x = 8, y = 1) %>%  
        # Georgiana spacing
        set_node_position(
           node = 14,
           x = 12, y = 1) %>% 
        # legend
        set_node_position(
           node = 15,
           x = 0, y = 10.3) %>% 
        # release survival label
        set_node_position(
           node = 16,
           x = 6.4, y = 10.25) %>%   
        # Sacramento survival label
        set_node_position(
           node = 17,
           x = 6.4, y = 8) %>%  
        # Sutter survival label
        set_node_position(
           node = 19,
           x = -0.6, y = 4) %>%  
        # Steamboat survival label
        set_node_position(
           node = 18,
           x = 4.4, y = 4) %>% 
        # Below Steamboat survival label
        set_node_position(
           node = 20,
           x = 11.4, y = 5) %>%
        # Georgiana survival label
        set_node_position(
           node = 21,
           x = 13.4, y = 2) %>% 
        # Below Georgiana survival label
        set_node_position(
           node = 22,
           x = 9.4, y = 2) %>% 
        # Sutter route label
        set_node_position(
           node = 23,
           x = 2, y = 6.9) %>%  
        # Steamboat route label
        set_node_position(
           node = 24,
           x = 5.5, y = 6.5) %>% 
        # Below Steamboat route label
        set_node_position(
           node = 25,
           x = 7.5, y = 6.9) %>%
        # Below Georgiana route label
        set_node_position(
           node = 26,
           x = 8.2, y = 3.85) %>% 
        # Georgiana route label
        set_node_position(
           node = 27,
           x = 11.8, y = 3.85) %>% 
        # node to make visible width wider
        set_node_position(
           node = 28,
           x = 15, y = 1) %>% 
        # legend lines
        set_node_position(
           node = 29,
           x = 0, y = 9.3) %>% 
        set_node_position(
           node = 30,
           x = 0, y = 9) %>% 
        set_node_position(
           node = 31,
           x = 0, y = 8.7) %>% 
        set_node_position(
           node = 32,
           x = 0, y = 8.4) %>% 
        set_node_position(
           node = 33,
           x = 2, y = 9.3) %>% 
        set_node_position(
           node = 34,
           x = 2, y = 9) %>% 
        set_node_position(
           node = 35,
           x = 2, y = 8.7) %>% 
        set_node_position(
           node = 36,
           x = 2, y = 8.4) %>% 
        # legend - route text
        set_node_position(
           node = 37,
           x = -1, y = 9)
     
     render_graph(graph)
  })
  
  # diagram text
  output$deltasurv_diagramtext <- renderText({ paste("Diagram of routing, detection, and survival from Sacramento to Benicia with 95% credible intervals") })
  
  
  ## barplot
  
  output$deltasurv_barplot <- renderPlot({
     
     summarydata <- deltasurvivalVar()
     
     # data prep
     route_survival_summary <- data.frame(Route = c("Sacramento River", "Sutter Slough", "Steamboat Slough", "Georgiana Slough"),
                                          Entrainment = c(summarydata$est[summarydata$parameter=="sac_route"],
                                                          summarydata$est[summarydata$parameter=="sutter_route"],
                                                          summarydata$est[summarydata$parameter=="steam_route"],
                                                          summarydata$est[summarydata$parameter=="georg_route"]),
                                          Survival = c(summarydata$est[summarydata$parameter=="sac_surv"],
                                                       summarydata$est[summarydata$parameter=="sutter_surv"],
                                                       summarydata$est[summarydata$parameter=="steam_surv"],
                                                       summarydata$est[summarydata$parameter=="georg_surv"]))
     route_arrow_x <- c(0,route_survival_summary$Entrainment[1],route_survival_summary$Entrainment[1]+route_survival_summary$Entrainment[2],
                        route_survival_summary$Entrainment[1]+route_survival_summary$Entrainment[2]+route_survival_summary$Entrainment[3],
                        route_survival_summary$Entrainment[1]+route_survival_summary$Entrainment[2]+route_survival_summary$Entrainment[3]+route_survival_summary$Entrainment[4])
     route_label_x <- c(route_survival_summary$Entrainment[1]/2,route_survival_summary$Entrainment[1]+(route_survival_summary$Entrainment[2])/2,
                        route_survival_summary$Entrainment[1]+route_survival_summary$Entrainment[2]+(route_survival_summary$Entrainment[3])/2,
                        route_survival_summary$Entrainment[1]+route_survival_summary$Entrainment[2]+route_survival_summary$Entrainment[3]+(route_survival_summary$Entrainment[4])/2)
     
     # make plot
     
     par(xpd=NA, mar = c(9.6, 4.1, 0.6, 0.6))
     routesurv_barplot <- barplot(height = route_survival_summary$Survival, 
                                  width = route_survival_summary$Entrainment, 
                                  names = route_survival_summary$Entrainment,
                                  col = c("lightblue","orange","blue4","purple"),
                                  ylab = "Survival Probability",
                                  xlab = "",
                                  space = 0, ylim = c(0,1), xaxs="i")
     arrows(x0 = routesurv_barplot,
            y0 = c(summarydata$l95[summarydata$parameter=="sac_surv"],summarydata$l95[summarydata$parameter=="sutter_surv"],summarydata$l95[summarydata$parameter=="steam_surv"],summarydata$l95[summarydata$parameter=="georg_surv"]),
            y1 = c(summarydata$u95[summarydata$parameter=="sac_surv"],summarydata$u95[summarydata$parameter=="sutter_surv"],summarydata$u95[summarydata$parameter=="steam_surv"],summarydata$u95[summarydata$parameter=="georg_surv"]),
            angle = 90,
            code = 3,
            length = 0.1)
     arrows(y0 = -0.03,
            x0 = route_arrow_x[1:4],
            x1 = route_arrow_x[2:5],
            angle = 90,
            code = 3,
            length = 0.08,
            lwd = 2)
     box()
     text(x = route_label_x,
          ## Move labels to just below bottom of chart.
          y = -0.15,
          ## Use names from the data list.
          labels = c("Sacramento River", "Sutter Slough", "Steamboat Slough", "Georgiana Slough"),
          ## Change the clipping region.
          xpd = NA,
          ## Rotate the labels by 35 degrees.
          srt = 35,
          ## Adjust the labels to almost 100% right-justified.
          adj = 0.965)
     title(xlab = "Route Probability", line = 8)
  })
  
  # barplot text
  output$deltasurv_barplottext <- renderText({ paste("Barplot of route-specific survival 
                                                     and entrainment probabilities from 
                                                     Sacramento to Benicia. Width of bars 
                                                     represent entrainment and height of 
                                                     bars represent survival with 95% 
                                                     credible intervals. Entrainment 
                                                     comprises the cumulative routing 
                                                     probability at each junction along 
                                                     each migration pathway.") })
  
  
  ## table
  
  output$deltasurv_table <- function() {
     
     summarydata <- deltasurvivalVar()
     
     # route specific entrainment and survival starting from Sacramento
     # with 95% credible intervals
     route_survival_summary_table <- data.frame(
        Route = c("Sacramento River", "Sutter Slough", "Steamboat Slough", "Georgiana Slough"), 
        Entrainment = c(paste0(summarydata$est[summarydata$parameter=="sac_route"]," (", 
                               summarydata$l95[summarydata$parameter=="sac_route"], "-", 
                               summarydata$u95[summarydata$parameter=="sac_route"], ")"), 
                        paste0(summarydata$est[summarydata$parameter=="sutter_route"]," (", 
                               summarydata$l95[summarydata$parameter=="sutter_route"], "-", 
                               summarydata$u95[summarydata$parameter=="sutter_route"], ")"), 
                        paste0(summarydata$est[summarydata$parameter=="steam_route"]," (", 
                               summarydata$l95[summarydata$parameter=="steam_route"], "-", 
                               summarydata$u95[summarydata$parameter=="steam_route"], ")"), 
                        paste0(summarydata$est[summarydata$parameter=="georg_route"]," (", 
                               summarydata$l95[summarydata$parameter=="georg_route"], "-", 
                               summarydata$u95[summarydata$parameter=="georg_route"], ")")), 
        Survival = c(paste0(summarydata$est[summarydata$parameter=="sac_surv"]," (", 
                            summarydata$l95[summarydata$parameter=="sac_surv"], "-", 
                            summarydata$u95[summarydata$parameter=="sac_surv"], ")"), 
                     paste0(summarydata$est[summarydata$parameter=="sutter_surv"]," (", 
                            summarydata$l95[summarydata$parameter=="sutter_surv"], "-", 
                            summarydata$u95[summarydata$parameter=="sutter_surv"], ")"), 
                     paste0(summarydata$est[summarydata$parameter=="steam_surv"]," (", 
                            summarydata$l95[summarydata$parameter=="steam_surv"], "-", 
                            summarydata$u95[summarydata$parameter=="steam_surv"], ")"), 
                     paste0(summarydata$est[summarydata$parameter=="georg_surv"]," (", 
                            summarydata$l95[summarydata$parameter=="georg_surv"], "-", 
                            summarydata$u95[summarydata$parameter=="georg_surv"], ")")))
     
     kable(route_survival_summary_table, row.names = F, "html", caption = "Sacramento to Benicia route-specific entrainment and survival probabilities with 95% credible intervals") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive", "bordered"), full_width = F, position = "left")}
  
  #########################################################################################
  
  
  ##%######################################################%##
  #                                                          #
  ####                      Movement                      ####
  #                                                          #
  ##%######################################################%##
  
  # This is used to get the release locations for the selected studyid
  movement_rel_locs <- reactive({
    req(input$movement_dataset)
    
    name <- input$movement_dataset
    original_studyid <- convert_descriptive_to_studyid(name)
    
    TaggedFish %>% 
      filter(study_id == original_studyid) %>% 
      select(release_location) %>% 
      distinct() %>% 
      pull(release_location)
    
  })
  
  
  movementVar <-  reactive({
    req(input$movement_dataset)
    
    files <- list.files("./detections/", full.names = T)
    
    # Choose the file that matches the studyID and read it in
    original_studyid <- convert_descriptive_to_studyid(input$movement_dataset)
    file <- files[str_detect(files, original_studyid)]
    df <- vroom(file, col_types = cols())
    
    df <- df %>% 
      left_join(
        TaggedFish %>% 
          select(FishID = fish_id, fish_release_date, release_location,
                 release_river_km),
        by = "FishID"
      ) %>% 
      mutate(fish_release_date = substr(fish_release_date, 1, nchar(fish_release_date)-3)) %>%
      mutate(fish_release_date = mdy_hm(fish_release_date)) %>% 
      # Remove any detections above the release RKM
      filter(GenRKM <= release_river_km) %>% 
      # If dataset has multi release filter for the selected one otherwise
      # leave it untouched https://stackoverflow.com/a/49411556
      filter(if(length(movement_rel_locs()) > 1 )  
        (release_location == input$movement_rel_select) else TRUE)
    
    # Calculate time from release to receiver and dist to receiver
    min_detects <- df %>% 
      group_by(release_location, FishID, GEN, GenRKM, fish_release_date, 
               release_river_km) %>% 
      summarize(
        min_time = min(time_pst)
      ) %>% 
      mutate(
        km_from_rel = release_river_km - GenRKM,
        days_since_rel = min_time - fish_release_date #difftime(min_time, fish_release_date, units = "days")
      ) %>% 
      arrange(FishID, min_time) %>% 
      mutate(km_day = km_from_rel / as.numeric(days_since_rel, units='days')) %>% 
      filter(km_from_rel > 0,
             abs(km_day) < 200) # DOUBLE CHECK IF THIS IS OK TO USE
    
    speeds <- min_detects %>% 
      group_by(release_location, GEN, GenRKM) %>% 
      summarize(
        N = n(),
        min_travel_days = min(as.numeric(days_since_rel, units='days'), na.rm=T),
        median_travel_days = median(days_since_rel),
        max_travel_days = max(days_since_rel),
        mean_km_day = mean(km_day),
        median_km_day = median(km_day)
      ) %>% 
      ungroup() %>% 
      mutate_if(is.difftime, as.numeric) %>% 
      arrange(desc(GenRKM))
    
    
  })
  
  
  output$movement_second_select <- renderUI({
    
    # Create the input selection for release location if there are multi release
    if (length(movement_rel_locs()) > 1 ) {
      selectInput("movement_rel_select", "Select Release Location", 
                  choices = movement_rel_locs())
    }
  })
  
  output$movement_gt <- render_gt({
    
    df <- movementVar()
    
    # If > 1 rel loc, filter data by selected rel loc
    if (length(movement_rel_locs()) > 1 ) {
      
      df <- df %>% 
        filter(release_location == input$movement_rel_select)
    }
    
    
    gt_tbl <- gt(data = df %>% 
                   select(-release_location), rowname_col = "GEN")
    
    gt_tbl %>% 
      tab_header(
        title = md("**Travel time summary statistics**"),
        subtitle = "Days and travel rates calculated from point of release to general
        receiver location"
      ) %>%
      tab_stubhead(label = "General Receiver Location") %>% 
      fmt_number(
        decimals = 1,
        columns = c(min_travel_days, median_travel_days, max_travel_days,
                    mean_km_day, median_km_day)
      ) %>% 
      cols_label(
        min_travel_days = "Min",
        median_travel_days = "Median", 
        max_travel_days = "Max", 
        mean_km_day = "Mean", 
        median_km_day = "Median"
      ) %>% 
      tab_spanner(
        label = "Time (days)",
        columns = c(min_travel_days, median_travel_days, max_travel_days)
      ) %>%
      tab_spanner(
        label = "Travel rate (km/day)",
        columns = c(mean_km_day, median_km_day)
      )
  })  
  
  movement_plot_min_detects <- reactive({
    req(input$movement_plot_dataset)
    
    # Directory of detection files
    files <- list.files("./detections/", full.names = T)
    original_studyid <- convert_descriptive_to_studyid(input$movement_plot_dataset)
    file <- files[str_detect(files, original_studyid)]
    detections <- vroom(file, col_types = cols())
    
    detections <- detections %>% 
      left_join(
        TaggedFish %>% 
          select(FishID = fish_id, fish_release_date, release_location,
                 release_river_km),
        by = "FishID"
      ) %>% 
      mutate(fish_release_date = mdy_hm(fish_release_date)) %>% 
      # Remove any detections above the release RKM
      filter(GenRKM <= release_river_km) 
    
    unique_detects <- detections %>% 
      group_by(
        FishID, GEN, GenRKM
      ) %>% 
      summarise(
        first_detect = min(time_pst),
        last_detect = max(time_pst)
      ) %>% 
      arrange(
        FishID, first_detect
      )
    
    # Vectorized way of getting swim speed
    reach_movement <- unique_detects %>% 
      group_by(FishID) %>% 
      mutate(
        reach_start = lag(GEN),
        reach_end = GEN,
        reach = paste0(reach_start, "_to_", reach_end),
        start_detect = dplyr::lag(last_detect),
        end_detect = first_detect,
        swim_time = difftime(end_detect, start_detect, units = "days"),
        rkm_start = lag(GenRKM),
        rkm_end = GenRKM,
        genrkm_diff = rkm_start- rkm_end,
        km_day = (genrkm_diff / as.numeric(swim_time)),
      ) %>% 
      filter(!is.na(rkm_start)) %>% 
      select(-c(start_detect, end_detect)) %>% 
      select(FishID, reach_start, reach_end, reach, first_detect, last_detect, 
             swim_time, genrkm_diff, rkm_start, rkm_end, km_day)
    
    # Build sequential reach names
    reach_names <- detections %>% 
      select(reach_start = GEN, GenRKM) %>% 
      distinct() %>% 
      arrange(desc(GenRKM)) %>% 
      mutate(
        reach_end = lead(reach_start),
        reach = ifelse(!is.na(reach_end), paste0(reach_start, "_to_", reach_end), NA)
      ) %>% 
      filter(!is.na(reach)) %>% 
      pull(reach)
    
    # Use only sequential reaches
    reach_movement <- reach_movement %>% 
      filter(reach %in% reach_names) %>% 
      mutate(midreach = (rkm_start + rkm_end)/2)
    
    reach_movement <- reach_movement %>% 
      mutate(
        midreach = factor(midreach, levels = unique(sort(reach_movement$midreach, decreasing = T)))
      )
    
  })
  
  output$movement_plot <- renderPlotly({
    df <- movement_plot_min_detects()
    
    p <- ggplot(data = df, aes(x = factor(midreach), y = km_day)) +
      geom_boxplot() +
      geom_jitter(shape=16, position=position_jitter(0.2)) +
      xlab("Midreach (RKM)") +
      ylab("Swim speed (KM/day)")
    
    ggplotly(p)
  })
}

# source('global.R')
# source('ui.R', local = TRUE)
# source('server.R', local = TRUE)

shinyApp(ui = ui, server = server)
