# Load libraries
library(tidyverse) # For data wrangling
library(data.table) # To load multiple files at once
library(Hmisc) # For the %nin% operator 
library(CoordinateCleaner) # To flag records with erroneous coordinates
library(countrycode) # Associated with coordinate cleaner
library(rgdal) # To handle vector data 
library(sp) # To handle vector data 

setwd("F:/Connectivity/data/01_cleanOccurrences")

# 1. Load occurrences of all species ----------------------------------------
data <- 
  list.files(pattern = "*.txt") %>% 
  map_df(~fread(., encoding = 'UTF-8'))

# 2. Remove some species that had wrong information about home ranges ----------------------------------
data <- data %>% 
  filter(`SCIENTIFIC NAME` %nin% c("Campephilus pollens", "Campylorhynchus nuchalis"))


# 4. Delete records with certain protocols -------------------------------------

# See all the protocols 
table(data$`PROTOCOL TYPE`)

# Remove Historical and Incidental protocols 
try <- data %>% filter(`PROTOCOL TYPE` %in% c("Traveling", "Stationary", "Area")) %>% droplevels()

# 3. Delete records with routes longer than 2 km and periods longer than 2 hours -------------------------------------
data <- data %>% 
  filter(`DURATION MINUTES` <= 120 | is.na(`DURATION MINUTES`)) %>% 
  filter(`EFFORT DISTANCE KM` <= 2 | is.na(`DURATION MINUTES`)) %>% 
  droplevels()


#plot data to get an overview
wm <- borders("world", colour="gray50", fill="gray50")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = data, aes(x = LONGITUDE, y = LATITUDE),
             colour = "darkred", size = 0.5)+
  theme_bw()

# 4. Remove duplicates -----------------------------------------------------------------

data <- data %>%
  
  # Remove records that have the same combination of the following variables 
  distinct(LONGITUDE,LATITUDE,`SCIENTIFIC NAME`,`OBSERVATION DATE`,`TIME OBSERVATIONS STARTED`, `DURATION MINUTES`, `OBSERVATION COUNT`, .keep_all = TRUE) 

# 5. Clean coordinates --------------------------------------------------------------

#convert country code from ISO2c to ISO3c
data$ISO3 <-  countrycode(data$`COUNTRY CODE`, origin =  'iso2c', destination = 'iso3c')


#flag problems

flags <- clean_coordinates(x = data, 
                           lon = "LONGITUDE", 
                           lat = "LATITUDE",
                           countries = "ISO3",
                           species = "SCIENTIFIC NAME",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros", "countries", "seas")) # most test are on by default

summary(flags)
plot(flags, lon = "LONGITUDE", lat = "LATITUDE")

# Check the flagged records
flagged_coordinates <- flags[!flags$.summary,]

# Export to ArcGIS Pro
write.csv(flagged_coordinates, "F:/Connectivity/outputs/01_cleanOcurrences/flaggedCoordinates.csv", row.names = FALSE)



# 6. Biomodelos and IUCN distributions ------------------------------------------------------------------


# Load sp distributions
Bio_IUCN <- readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN", layer = "ALL_sp", verbose = TRUE)

# Add a column named binomial
data_Bio_IUCN <- data %>% 
  
  mutate(BINOMIAL = `SCIENTIFIC NAME`) %>%
  
  filter(BINOMIAL %in% c("Henicorhina leucophrys", "Phaethornis malaris", "Grallaria hypoleuca",
                          "Picumnus squamulatus", "Picumnus lafresnayi", "Crypturellus soui", "Grallaria nuchalis", 
                         "Crypturellus kerriae", "Crypturellus erythropus", "Mitu salvini", "Mitu tomentosum", 
                         "Campephilus pollens", "Pyrrhura melanura", "Micrastur semitorquatus", "Spizaetus ornatus",
                         "Buteo albigula", "Ara militaris", "Ara macao")) %>%
  
  droplevels()

# Remove records that fall outside the distribution 
data_Bio_IUCN_clean <- cc_iucn(x = data_Bio_IUCN, 
                range = Bio_IUCN, 
                lon = "LONGITUDE", 
                lat = "LATITUDE",
                species = "BINOMIAL",
                buffer = 0.1,
                value = "clean",
                verbose = TRUE)

# check 
write.csv(data_Bio_IUCN, "F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/COMPLETE_RECORDS_EBIRD.csv", row.names = FALSE)
write.csv(data_Bio_IUCN_clean, "F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/CLEANED_RECORDS_EBIRD.csv", row.names = FALSE)



# 7. IUCN distribution and geographic outliers  --------------------------------------------------------------

# Load sp distributions

CB <- readOGR("F:/Connectivity/data/species_distributions/IUCN/low 5 km", layer = "Crypturellus berlepschi", verbose = TRUE)
TJ <- readOGR("F:/Connectivity/data/species_distributions/IUCN/low 5 km", layer = "Tangara johannae", verbose = TRUE)

ALL <- rbind(CB, TJ)
# Luego proyectarpor que no esta aplicando el buffer

try <- data %>% rename(BINOMIAL = `SCIENTIFIC NAME`) %>% filter(BINOMIAL %in% c("Crypturellus berlepschi", "Tangara johannae")) %>% droplevels()


# With a buffer of ~10 km . buffer = 0.1,
try2 <- cc_iucn(x = try, 
                range = ALL, 
                lon = "LONGITUDE", 
                lat = "LATITUDE",
                species = "BINOMIAL",
                verbose = TRUE)

# check 
write.csv(try, "F:/Connectivity/data/CHECK_COMPLETE.csv", row.names = FALSE)
write.csv(try2, "F:/Connectivity/data/CHECK_BUFFER.csv", row.names = FALSE)

# 8. Divide records prior and after 2000 ------------------------------------------------------------------


