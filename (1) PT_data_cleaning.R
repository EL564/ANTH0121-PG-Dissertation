################################################################################
#
#
#(1) PT_data_cleaning
#script includes: merging and cleaning data, plus interactive map
#author: Ella Lipscombe 
#
#
################################################################################

#start session 

#load data

#set working directory
setwd("D:/MSc HEB 2023-2025/Modules 2024-2025/ANTH0121 Dissertation/Methodology/R script") #set work directory 

#shortcut (load csvs created in this script)
SU_combined <- read.csv ("SU_combined.csv", header = TRUE) #code to import master data set

#load raw ENDOW data and distances from QGIS 
PT_HH <- read.csv("HH_data_pantanal.csv", header=TRUE) #household egos
PT_indiv <- read.csv("PT_Indiv.csv", header=TRUE) #individual-level data
PT_distances <- read.csv("D:/MSc HEB 2023-2025/Modules 2024-2025/ANTH0121 Dissertation/Methodology/R script/Distances_SU.csv") #distances
PT_indiv_subset <- read.csv("PT_indiv_subset.csv", header=TRUE) #subsetted household egos from the full individual data frame
  
#packages 
install.packages(c("dplyr","sf","INLA","fmesher","mapview"))
library(dplyr)
library(sf)
library(INLA)
library(fmesher)
library(mapview)
library(ggplot2)

################################################################################

#clean data
distance_subset <- PT_distances [,c(1,5)] #subset distances data to include relevant columns only 

names(distance_subset)[names(distance_subset) == "SharingUni"] <- "SharingUnitID" #standardising names
names(distance_subset)[names(distance_subset) == "HubDist"] <- "DistancefromPA"

#calculate no. members per SU 
SU_members <- PT_indiv %>%
  group_by(SharingUnitID) %>%
  summarise(SU_members = n()) #count the number of individuals per sharing unit

PT_indiv_subset <- PT_indiv_subset %>%
  left_join(SU_members, by = "SharingUnitID") #merge with subset data of isolated egos/households of interest 

View(PT_indiv_subset)

#combine data at the SU level
indiv_combined <- left_join(PT_indiv_subset, PT_HH, by="SharingUnitID") #adding the coordinates and to PT_ind_data_MD dataset -- values align by matching SharingUnitID
indiv_combined <- left_join (indiv_combined, distance_subset, by="SharingUnitID") #do the same with distance from PA
indiv_combined <- indiv_combined[which(indiv_combined$Longitude!="NA"),] #remove missing coordinates

SU_combined <- indiv_combined #subset data into heads of households and focals only -- these are the observations we are interested in
SU_combined$WealthMeasure <- NULL #remove empty variable WealthMeasure to leave us with 62 observations and 24 variables

#standardising the capitalisation of values
SU_combined$GroupID[SU_combined$GroupID == "Catho"] <- "Catholic" #rename values in column GroupID
SU_combined$GroupID[SU_combined$GroupID == "catho"] <- "Catholic" 
SU_combined$GroupID[SU_combined$GroupID == "spiritualist"] <- "Spiritualist"

SU_combined$MaritalStatus[SU_combined$MaritalStatus == "widower"] <- "Widower" #rename values in column MaritalStatus
SU_combined$MaritalStatus[SU_combined$MaritalStatus == "single"] <- "Single"

SU_combined$Occupation[SU_combined$Occupation == "healthassistant"] <- "Health Assistant" #rename values in column Occupation 
SU_combined$Occupation[SU_combined$Occupation == "ParkRanger"] <- "Park Ranger"


#re-leveling categorical variable for future analysis
SU_combined = SU_combined %>%  
  mutate(GroupID_numerical = case_when(
    GroupID == "Protestant" ~ 1, 
    GroupID == "Catholic" ~ 2,
    GroupID == "Adventist" ~ 3,
    GroupID == "Spiritualist" ~4))


#renaming variable names 
names(SU_combined)[names(SU_combined) == "displaced"] <- "DisplacementStatus" #standardising the capitalisation of variable names 
names(SU_combined)[names(SU_combined) == "OtherNoetic"] <- "NoeticWealth"
names(SU_combined)[names(SU_combined) == "GroupID"] <- "Religion" 
names(SU_combined)[names(SU_combined) == "Status"] <- "PositionofPower" 
names(SU_combined)[names(SU_combined) == "WEALTH"] <- "MaterialWealth"
names(SU_combined)[names(SU_combined) == "GroupID_numerical"] <- "ReligionLevel" 
names(SU_combined)[names(SU_combined)== "IndivID"] <- "Ego"

#save clean master dataset as CSV 
write.csv(SU_combined, "SU_combined.csv", row.names = FALSE) #save cleaned master dataset at the SU level
SU_combined <- read.csv ("SU_combined.csv", header = TRUE) #code to import master data set

rm(food_SU_subset,indiv_combined,PT_food_sec,PT_HH, PT_indiv, SU_members) #remove all the other data sets and subsets so one remains
View(SU_combined) #view clean master data set at the SU level


################################################################################


#interactive map of study site (not used in thesis)

PT_sites <- st_as_sf(
  SU_combined, coords = c("Longitude","Latitude"), crs = 4326) #create an interactive map of the displaced and non-displaced communities 

map <- mapview(PT_sites, 
             col.regions = "darkred", lwd = 1,
        alpha.regions = 0.2, #set the transparency level of the regions to 20%. The lower the value, the more transparent the regions will be, allowing for better visibility of overlapping areas.
        alpha = 0.3) #this sets the transparency for the points or features themselves to 30%. Similar to alpha.regions, this helps in viHHalizing the data without overwhelming the viewer.
map


#saving map
leaflet_map <- map@map  #extract the leaflet object

#save as interactive HTML
htmlwidgets::saveWidget(leaflet_map, file = "map.html", selfcontained = TRUE)

#code citation: Pirani, M., 2019. R Studio 4.4.3.

################################################################################

#end session

################################################################################