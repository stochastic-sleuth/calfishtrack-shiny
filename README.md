# Central Valley Enhanced Acoustic Tagging Project

## Background
Since 2012, juvenile salmon have been tagged and tracked throughout California's Central Valley (CCV) using Juvenile Salmon Acoustic Telemetry System (JSATS) technology. This technology allows researchers to monitor the movement and survival rates of various populations over 500 river kilometers, from Redding to the Pacific Ocean. The data compiled here is a result of the hard work and coordination between various federal and state agencies, universities, and private consultants, with the goal of conserving and restoring California's once abundant, but now imperiled, Chinook salmon populations. This data is open access and is managed and hosted by the National Marine Fisheries Service here: https://oceanview.pfeg.noaa.gov/erddap/tabledap/FED_JSATS_detects.html.

### Receiver Deployments
An interactive map of acoustic receivers (JSATS, Vemco, Real-time) deployed throughout the CCV by water year, which is defined as the 12 month period beginning October 1 and ending September 30 of the following calendar year (i.e. water year 2017 spans 10/1/17-9/30/18). The deployment period of individual receivers is displayed in a table below the map, once a general location is clicked. This map is useful to identify a particular site of interest, and to see when coverage existed in that area and which receiver(s) potentially recorded migrating fish.

### Hydrology
The hydrology of the Sacramento River is highly variable, and is largely driven by storm events and associated runoff in the winter, followed by dam controlled releases for agricultural purposes in the summer. These interactive graphs of Sacramento (SR) and San Joaquin River (SJR) flows (cubic feet per second) at four different CDEC (California Data Exchange Center) stations on each river: KES (Keswick Reservoir, SR), BND (Bend Bridge, SR), BTC (Butte City, SR), WLK (Wilkins Slough, SR), MEN (Mendota, SJR), NEW (Newman, SJR), VNS (Vernalis, SJR), and MSD (Mossdale Bridge, SJR) displays the annual and inter-annual variation in flows from upstream to downstream gauging stations. In addition, the temperatures at KES, BND, WLK, and MSD are shown. These plots were created using the R package dygraphs.

### Outmigration animation
The Outmigration Animation tab visualizes fish outmigration using detection data over time. Unique fish detections at each receiver location are plotted by day, and numbers of unique detections per receiver are overlaid in each point. This animation helps to visualize the movement of fish across a broad geographic landscape, and demonstrates the differences in outmigration timing across space and time among specific populations. This animation was created using the R packages leaflet and leaflet.extras. 

### Data Explorer
Data explorer allows users to look at the tagged fish data individually or across different study populations in a number of ways. Plots can be created interactively, choosing variables to explore (length, weight), and with different plot types (boxplot, histogram, density) as well as summary tables.

### Time of Day
Time of day allows users to visually explore behavioral differences in fish movement between night and day throughout the migratory corridor. Using detection times we can look at the distribution in movement times for an entire study group, or the movement times for a specific receiver location. 

### Survival
Survival estimates are calculated by summarizing the detection history of study populations, and assuming a fish has died when it is not detected by subsequent downstream receivers. We used a CJS survival model in RMark to estimate reach specific (survival per 10 river kilometers) and cumulative (from release location to last downstream detection) survival rates. We used a simple survival model (survival and detection efficiency is a function of time) to derive these estimates. For each survival tab, a map is displayed with all receiver locations used to generate survival for each population, and a figure and table displays the estimates and associated error. The unique number of fish detected at each receiver location is also provided.

For the Delta route-specific survival tab, we ran a Bayesian multi-state mark-recapture model in JAGS to obtain survival and route-entrainment probability estimates for the Delta between Sacramento and Benicia. We modeled routing and survival in four possible routes through the Delta: Sacramento River mainstem, Sutter Slough, Steamboat Slough, and Georgiana Slough. In this tab, there is a diagram showing survival and routing estimates for each reach and junction in the Delta, and a barplot and table summarizing survival and migration route probabilities for each of the four migratory pathways. This section was developed by the Quantitative Fisheries Ecology Section of USGS Western Fisheries Research Center. For questions or comments, send an email to our group at gs-b-crrl_qfes@usgs.gov .

*These survival results are preliminary and for discussion purposes only. Detection data has not been filtered for predator detections, and survival estimates have not been adjusted for any potential premature tag failures.*

### Movement
Fish movement is summarized by study population in a table format. For each receiver location, the minimum, median, and maximum travel time is calculated in days and kilometers per day. The number of unique fish detected at each receiver location is displayed as well. We will update this page to display travel times in an interactive plot soon, so stay tuned.

### Questions or comments?
This Shiny App was developed by UCSC/NOAA, Southwest Science Fisheries Center, Fisheries Ecology Division. Source code used to create this site can be viewed here . If you have any questions or comments please feel free to send us an email . 
