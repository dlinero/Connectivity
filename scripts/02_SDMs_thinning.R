
library("tidyverse") # for data wrangling
library("Hmisc") # for %nin% 
library("spThin") # for spatial thinning 
library("raster") # For handling environmental variables
library("rgeos") # for buffers 
library("adehabitatHR") # For minimum convex polygons 
library("sf") # For turning dataframes into shapefiles 
library("janitor") # For cleanning dataset column names
library("rgdal")
library("foreach")
library("doParallel")
library("ranger")
library("kableExtra")
library("doSNOW") # To see progress bar in paralellized loops 


# detach("package:wallace", unload = TRUE) # Unload Wallace
# remotes::install_github("https://github.com/wallaceEcoMod/wallace/tree/multiSp", lib = "F:/R/R-4.1.1/library") # Install demo testing version 
# Go back to the original version because this one didn't worked

# install.packages("wallace", lib = "F:/R/R-4.1.1/library")

library("wallace")


run_wallace()

# --------------------------------------------------------------------------------------------

# Load cleaned data 
data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Create a variable with the name of the columns 
columns <- colnames(data) 
columns <- columns[columns  %nin% c("SCIENTIFIC.NAME", "LONGITUDE", "LATITUDE")]

# Change the names of some columns to comply with Wallace requirements 
data <- data %>%
  
  # Remove Spizaetus isidori, as we do not have yet the Colombia points 
  filter(SCIENTIFIC.NAME != "Spizaetus isidori") %>% 
  
  # Rename variables 
  rename(name = SCIENTIFIC.NAME, longitude = LONGITUDE, latitude = LATITUDE) %>%
  
  # The renamed variables should go first 
  dplyr::select(name, longitude, latitude, columns) %>%
  
  droplevels()

# Create a table with the number of records per species 
n_recors_species <- data %>% group_by(name) %>% summarise(records = n())


# 1. Spatial thinning - 1km  ---------------------------------------------

# Select species that have less of 1500 records (for running time)
sp_names <-  n_recors_species %>% filter(records < 1500) %>% pull(name) %>% as.character() 
 
# Este paso lo hago para eliminar las especies que se corrieron ayer, después se tiene que borrar
sp_names <- subset(sp_names, sp_names %nin% n_thinned_records_1km$name)

n_thinned_records_1km <- matrix(nrow = 1, ncol = 2) %>% data.frame()

colnames(n_thinned_records_1km) <- c("name", "n_records")

for (i in 1:length(sp_names)){ 
  
  print(i)
  # Spatial thinning 
  
  x <- data %>% filter(name == sp_names[i]) %>% droplevels()
  
  # Run 100 iterations of spatial thinning, with a distance of 1 km
  data_thinned_1km <- spThin::thin(loc.data = x, lat.col = "latitude", long.col = "longitude", spec.col = "name", 
                                   thin.par = 1, reps = 100, locs.thinned.list.return = TRUE, write.files = FALSE, verbose = TRUE)
  
  # find the iteration that returns the max number of occurrences
  maxThin <- which(sapply(data_thinned_1km, nrow) == max(sapply(data_thinned_1km, nrow)))
  # if there's more than one max, pick the first one
  maxThin <- data_thinned_1km[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
  # subset dataset to match only thinned occurrences with maximum number of localities
  x <- x[as.numeric(rownames(maxThin)),]  
  
  # Save number of thinned points 
  n_records <- x %>% group_by(name) %>% summarise(n_records = n())
  n_thinned_records_1km <- rbind(n_thinned_records_1km, n_records)
  
  # Define name of the file that is going to be exported 
  path_points <- paste0("F:/Connectivity/outputs/02_SDMs/Intentos/", x$name[1], "_thinned_1km.shp")
  
  # Clean names of columns
  x <- x %>%
    clean_names()
  
  # Turn to shapefile 
  x_shp <- sf::st_as_sf(x, coords = c("longitude", "latitude"), crs = 4326)
  
  # Export
  sf::st_write(x_shp, path_points)
  
  # Minimum convex polygon with a buffer of 0.1°
  
  # Extract only the coordinates
  x_mcp <- x[c('longitude', 'latitude')]
  # Assign the columns with the coordinates
  sp::coordinates(x_mcp) <- ~ longitude + latitude
  # Calculate minimum convex polygon
  calibration_Area <- mcp(x_mcp, percent = 100)
  # Apply buffer 
  calibration_Area <- rgeos::gBuffer(calibration_Area, width = 0.1)
  # Assign coordinates
  crs(calibration_Area) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  # Define name of the file that is going to be exported 
  path_mcp <- paste0("F:/Connectivity/outputs/02_SDMs/Intentos/", x$name[1], "_mcp_10km.shp")
  
  # Export 
  raster::shapefile(calibration_Area, filename = path_mcp)
  
  # Remove temporary variables
  rm(x, data_thinned_1km, maxThin, n_records, path_points, x_shp, x_mcp, calibration_Area, path_mcp)
  
}

# Run it again for the species that have less than 16 records when using the 5 km distance


#create the cluster
my.cluster <- parallel::makeCluster(
  3, 
  type = "PSOCK"
)

# Register the cluster
registerDoSNOW(my.cluster)

# Define the number of iterations 
iterations <- 3

# Define progress bar (code extracted from https://www.py4u.net/discuss/859837)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Select species with less than 16 records when using the 5 km distance 
sp <-  c("Mitu salvini", "Grallaria bangsi", "Myioborus flavivertex")

results <- foreach(i = 1:iterations, .combine = rbind, .packages = "tidyverse", .options.snow = opts) %dopar% {
  
  x <- data %>% dplyr::filter(name == sp[i]) %>% droplevels()
  
  # Run 100 iterations of spatial thinning, with a distance of 1 km
  data_thinned_1km <- spThin::thin(loc.data = x, lat.col = "latitude", long.col = "longitude", spec.col = "name", 
                                   thin.par = 1, reps = 100, locs.thinned.list.return = TRUE, write.files = FALSE, verbose = FALSE)
  
  # find the iteration that returns the max number of occurrences
  maxThin <- which(sapply(data_thinned_1km, nrow) == max(sapply(data_thinned_1km, nrow)))
  # if there's more than one max, pick the first one
  maxThin <- data_thinned_1km[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
  # subset dataset to match only thinned occurrences with maximum number of localities
  x1 <- x[as.numeric(rownames(maxThin)),]  
  
  x1
  
}

close(pb)
stopCluster(my.cluster) 

# Export
write.csv(results, "F:/Connectivity/outputs/02_SDMs/Intentos/final/thinned_data_1km.csv", row.names = FALSE)

# 2. Spatial thinning - 5km  ---------------------------------------------


# Load cleaned data 
data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE, stringsAsFactors = FALSE)

# Create a variable with the name of the columns 
columns <- colnames(data) 
columns <- columns[columns  %nin% c("SCIENTIFIC.NAME", "LONGITUDE", "LATITUDE")]

# Change the names of some columns to comply with Wallace requirements 
data <- data %>%
  
  # Remove Spizaetus isidori, as we do not have yet the Colombia points 
  filter(SCIENTIFIC.NAME != "Spizaetus isidori") %>% 
  
  # Rename variables 
  rename(name = SCIENTIFIC.NAME, longitude = LONGITUDE, latitude = LATITUDE) %>%
  
  # The renamed variables should go first 
  dplyr::select(name, longitude, latitude, columns) %>%
  
  droplevels()

# Create a table with the number of records per species 
n_recors_species <- data %>% group_by(name) %>% summarise(records = n())


#create the cluster
my.cluster <- parallel::makeCluster(
  6, 
  type = "PSOCK"
)

# Register the cluster
registerDoSNOW(my.cluster)

# Define the number of iterations 
iterations <- length(which(n_recors_species$records < 1500)) %>% as.numeric()

# Define progress bar (code extracted from https://www.py4u.net/discuss/859837)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Select species that have less than 1500 records 
sp <-  n_recors_species %>% filter(n_recors_species$records < 1500) %>% pull(name) %>% as.character() 

results <- foreach(i = 1:iterations, .combine = rbind, .packages = "tidyverse", .options.snow = opts) %dopar% {
  
  x <- data %>% filter(name == sp[i]) %>% droplevels()
  
  # Run 100 iterations of spatial thinning, with a distance of 5 km
  data_thinned_5km <- spThin::thin(loc.data = x, lat.col = "latitude", long.col = "longitude", spec.col = "name", 
                                   thin.par = 5, reps = 100, locs.thinned.list.return = TRUE, write.files = FALSE, verbose = FALSE)
  
  # find the iteration that returns the max number of occurrences
  maxThin <- which(sapply(data_thinned_5km, nrow) == max(sapply(data_thinned_5km, nrow)))
  # if there's more than one max, pick the first one
  maxThin <- data_thinned_5km[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
  # subset dataset to match only thinned occurrences with maximum number of localities
  x1 <- x[as.numeric(rownames(maxThin)),]  
  
  x1
  
}

# Export file 
write.csv(results, "F:/Connectivity/outputs/02_SDMs/Intentos/thinned_data_5km.csv", row.names = FALSE)

# ESTE LO HICE EN MI COMPUTADOR POR QUE POR ALGUNA RAZÓN ERA MÁS RÁPIDO 
results <- read.csv("F:/Connectivity/outputs/02_SDMs/Intentos/final/LARGE/thinned_largedata_5km_1.csv", header = TRUE)

# Calculate number of remaining records per species 
n_records_5km <- results %>% group_by(name) %>% summarise(records = n()) 

sp <- results %>% pull(name) %>% unique()

# For each species, calculate the minimum convex polygon and export shapefiles 

for (i in 1:length(sp)){ 
  
  print(i)
  
  # Filter species 
  x <- results %>% filter(name == sp[i]) %>% droplevels()
  
  # Define name of the file that is going to be exported 
  path_points <- paste0("F:/Connectivity/outputs/02_SDMs/Intentos/final/LARGE/", x$name[1], "_thinned_5km.shp")
  
  # Clean names of columns
  x <- x %>%
    clean_names()
  
  # Turn to shapefile 
  x_shp <- sf::st_as_sf(x, coords = c("longitude", "latitude"), crs = 4326)
  
  # Export
  sf::st_write(x_shp, path_points)
  
  # Minimum convex polygon with a buffer of 0.1°
  
  # Extract only the coordinates
  x_mcp <- x[c('longitude', 'latitude')]
  # Assign the columns with the coordinates
  sp::coordinates(x_mcp) <- ~ longitude + latitude
  # Calculate minimum convex polygon
  calibration_Area <- mcp(x_mcp, percent = 100)
  # Apply buffer 
  calibration_Area <- rgeos::gBuffer(calibration_Area, width = 0.1)
  # Assign coordinates
  crs(calibration_Area) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  # Define name of the file that is going to be exported 
  path_mcp <- paste0("F:/Connectivity/outputs/02_SDMs/Intentos/final/LARGE/", x$name[1], "_mcp_10km_thin5km.shp")
  
  # Export 
  raster::shapefile(calibration_Area, filename = path_mcp)
  
  # Remove temporary variables
  rm(x,  path_points, x_shp, x_mcp, calibration_Area, path_mcp)
  
}


close(pb)
stopCluster(my.cluster) 


# 3. Spatial thinning - 5 km with a polygon grid ----------------------------------------------

# Extract species with more than 5000 records 
sp <- n_recors_species %>% arrange(records) %>% filter(records > 5000) %>% pull(name) %>% unique()

# Create a table to collect thinned data for each species
large_thinned_data <- data[1, ] %>%
  clean_names()

for (i in 1:length(sp)){ 
  
  print(i)
  
  # Filter species 
  x <- data %>% filter(name == sp[i]) %>% droplevels()
  
  # Turn into vector 
  large_shp <- sf::st_as_sf(x, coords = c("longitude", "latitude"), crs = 4326)
  
  # Project to Web Mercator to intersect with a grid 
  dataset_transform <- st_transform(large_shp, 3857)
  
  # Define grid size 
  gridCell_Size <- 5000
  
  # Create grid 
  grid <- 
    
    # Find bounding box of the points
    st_as_sfc(st_bbox(dataset_transform)) %>% 
    
    # Make a grid within the bounding box, with the cell size in meters 
    st_make_grid(cellsize = gridCell_Size) %>% 
    
    # Add an ID to each cell 
    cbind(data.frame(ID = sprintf(paste("GID%0",nchar(length(.)),"d",sep=""), 1:length(.))))  %>% 
    
    st_sf()
  
  
  # Intersect points with the grid, so that each point has the ID of the cell it is in
  large_shp <- dataset_transform %>% st_join(grid, join = st_intersects )
  
  # transform the coordinates again to WGS84
  grid <- st_transform(grid, 4326)
  large_shp <- st_transform(large_shp, 4326)
  
  
  large_shp <- large_shp %>% 
    
    # Retain coordinates columns
    mutate(latitude = unlist(map(large_shp$geometry,2)),
           longitude = unlist(map(large_shp$geometry,1))) %>%
    
    # Clean names of the columns 
    clean_names()
  
  
  # Randomly select only one point per grid 
  thinned_data <- large_shp %>% group_by(id) %>% sample_n(1)
  
  
  # Define name of the file that is going to be exported 
  path <- paste0("F:/Connectivity/outputs/02_SDMs/Intentos/final/LARGE/TRY/", thinned_data$name[1], "_thinned_5km.shp")
  
  # Export 
  sf::st_write(thinned_data, path)
  
  # Remove extra columns 
  thinned_data <- thinned_data %>% st_drop_geometry() %>% 
    ungroup() %>% dplyr::select(-id)
  
  # Compile all into a single table 
  large_thinned_data <- rbind(large_thinned_data, thinned_data)

  # Minimum convex polygon with a buffer of 0.1°
  
  # Extract only the coordinates
  x_mcp <- thinned_data[c('longitude', 'latitude')]
  # Assign the columns with the coordinates
  sp::coordinates(x_mcp) <- ~ longitude + latitude
  # Calculate minimum convex polygon
  calibration_Area <- mcp(x_mcp, percent = 100)
  # Apply buffer 
  calibration_Area <- rgeos::gBuffer(calibration_Area, width = 0.1)
  # Assign coordinates
  crs(calibration_Area) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  # Define name of the file that is going to be exported 
  path_mcp <- paste0("F:/Connectivity/outputs/02_SDMs/Intentos/final/LARGE/TRY/", x$name[1], "_mcp_10km_thin5km.shp")
  
  # Export 
  raster::shapefile(calibration_Area, filename = path_mcp)
  
  
}

# Remove the Grallaria bangsi record
large_thinned_data <- large_thinned_data %>% filter(name != "Grallaria bangsi") %>% droplevels()


# Export file 
write.csv(large_thinned_data, "F:/Connectivity/outputs/02_SDMs/Intentos/final/LARGE/TRY/thinned_largedata_5km_2.csv", row.names = FALSE)


# 4. Create one consolidated excel with all thinned data ------------------------------------------------

# Load 5km thinned data for species with less than 1500 records
data_1 <- read.csv("F:/Connectivity/outputs/02_SDMs/Intentos/final/thinned_data_5km.csv") %>%
  
  # Remove the species that will have a 1 km thinning 
  filter(name %nin%  c("Mitu salvini", "Grallaria bangsi", "Myioborus flavivertex")) %>% 
  
  # Clean names of the columns 
  clean_names()

# Load 5km thinned data for species with between 1500 to 3000 records
data_2 <- read.csv("F:/Connectivity/outputs/02_SDMs/Intentos/final/LARGE/thinned_largedata_5km_1.csv") %>%
  
  # Clean names of the columns 
  clean_names()

# Load 1 km thinned data for species that had less than 16 records when using the 5km distance
data_3 <- read.csv("F:/Connectivity/outputs/02_SDMs/Intentos/final/thinned_data_1km.csv") %>%
  
  # Clean names of the columns 
  clean_names()

# Load 5km thinned data for species with more than 5000 records
data_4 <-  read.csv("F:/Connectivity/outputs/02_SDMs/Intentos/final/LARGE/TRY/thinned_largedata_5km_2.csv")

# Merge everything 
data_total <- rbind(data_1, data_2)
data_total <- rbind(data_total, data_3)
data_total <- rbind(data_total, data_4)

data_total <- data_total %>% dplyr::select(name, longitude, latitude, global_unique_identifier, 
                                    common_name, subspecies_common_name, subspecies_scientific_name, 
                                    observation_count, country, country_code, state, county, county_code, 
                                    iba_code, locality, locality_id, locality_type, observation_date, time_observations_started, 
                                    observer_id, sampling_event_identifier, protocol_type, protocol_code, duration_minutes, effort_distance_km, 
                                    group_identifier, year_observation)

# Remove Crypturellus kerriae because in both thin distances it only retains 6 records
data_total <- data_total %>% filter(name != "Crypturellus kerriae") %>% droplevels()

# Number of records per species 
n_records <- data_total %>% group_by(name) %>% summarise(records = n())

# Export consolidated data 
write.csv(data_total, "F:/Connectivity/outputs/02_SDMs/Intentos/final/FINAL_THINNED_DATA.csv", row.names = FALSE)



# 5. Sensitive data ------------------------------------------------


# Load sensitive data
s_data <- read.delim("F:/Connectivity/outputs/02_SDMs/Intentos/sensitive data/ebd_sensitive_relAug-2021_JorgeVelasquez.txt")

# Filter Spizaetus isidori 
s_data <- s_data %>% filter(SCIENTIFIC.NAME == "Spizaetus isidori") %>% droplevels()

# Let's clean that dataset ------------------------------------------------
# Select travelling, stationary and area protocols  
s_data <- s_data %>% 
  filter(PROTOCOL.TYPE %in% c("Traveling", "Stationary", "Area")) %>% 
  droplevels()

#Delete records with routes longer than 3 km and periods longer than 4 hours 
s_data <- s_data %>% 
  filter(!is.na(DURATION.MINUTES) & DURATION.MINUTES <= 240)%>% 
  filter(!is.na(EFFORT.DISTANCE.KM) & EFFORT.DISTANCE.KM <= 3) %>% 
  droplevels()


# Remove duplicates
s_data <- s_data %>%
  
  distinct(LONGITUDE,LATITUDE, SCIENTIFIC.NAME, OBSERVATION.DATE, TIME.OBSERVATIONS.STARTED, DURATION.MINUTES, OBSERVATION.COUNT, .keep_all = TRUE)


# Coordinate cleaner
#convert country code from ISO2c to ISO3c
s_data$ISO3 <-  countrycode::countrycode(s_data$COUNTRY.CODE, origin =  'iso2c', destination = 'iso3c')


#flag suspicious records using the most common tests
flags <-CoordinateCleaner::clean_coordinates(x = s_data, 
                           lon = "LONGITUDE", 
                           lat = "LATITUDE",
                           countries = "ISO3",
                           species = "SCIENTIFIC.NAME",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros", "countries", "seas")) 

# Summary of the number of flagged records
summary(flags)

# Load the merged IUCN and ayerbe distributions for all sp with a 0.1° buffer

Aye_IUCN_buffer <- readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN", layer = "ALL_sp_buffer", verbose = TRUE)

# Add a column named binomial to the ebird datase, since is the column that the 
# shapefile has 
data_Aye_IUCN <- s_data %>% 
  
  mutate(BINOMIAL = SCIENTIFIC.NAME) 

# Remove records that fall outside the distribution with the 0.1° buffer 
data_Aye_IUCN_clean <- CoordinateCleaner::cc_iucn(x = data_Aye_IUCN, 
                               range = Aye_IUCN_buffer, 
                               lon = "LONGITUDE", 
                               lat = "LATITUDE",
                               species = "BINOMIAL",
                               verbose = TRUE)


# Join records within and outside Colombia
# Load cleaned data 
data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Set same columns names
data <- data %>% filter(SCIENTIFIC.NAME == "Spizaetus isidori") %>% droplevels() %>% clean_names()
colnames(data) <- gsub(x = colnames(data), pattern = "\\_", replacement = "\\.") %>% toupper()

# Join
s_i_data<- full_join(data, data_Aye_IUCN_clean)

# Create year field
s_i_data <- s_i_data %>% mutate("YEAR.OBSERVATION" = format(OBSERVATION.DATE, format = "%Y")) %>%
  mutate(YEAR.OBSERVATION = format.Date(YEAR.OBSERVATION, format = "%Y"))

# Remove records before 2000
s_i_data <- s_i_data %>% filter(YEAR.OBSERVATION >= 2001) %>% droplevels()


write.csv(s_i_data, "F:/Connectivity/outputs/02_SDMs/Intentos/sensitive data/cleanedData_phase2_SpizaetusIsidori.csv", row.names = FALSE)


s_i_data <- read.csv("F:/Connectivity/outputs/02_SDMs/Intentos/sensitive data/cleanedData_phase2_SpizaetusIsidori.csv")

# 2. Spatial thinning - 5km  ---------------------------------------------

# Create a variable with the name of the columns 
columns <- colnames(s_i_data) 
columns <- columns[columns  %nin% c("SCIENTIFIC.NAME", "LONGITUDE", "LATITUDE")]

# Change the names of some columns to comply with Wallace requirements 
s_i_data <- s_i_data %>%
  
  # Rename variables 
  rename(name = SCIENTIFIC.NAME, longitude = LONGITUDE, latitude = LATITUDE) %>%
  
  # The renamed variables should go first 
  dplyr::select(name, longitude, latitude, columns) %>%
  
  droplevels()

# Run 100 iterations of spatial thinning, with a distance of 5 km
data_thinned_5km <- spThin::thin(loc.data = s_i_data, lat.col = "latitude", long.col = "longitude", spec.col = "name", 
                                   thin.par = 5, reps = 100, locs.thinned.list.return = TRUE, write.files = FALSE, verbose = TRUE)
  
# find the iteration that returns the max number of occurrences
maxThin <- which(sapply(data_thinned_5km, nrow) == max(sapply(data_thinned_5km, nrow)))

# if there's more than one max, pick the first one
maxThin <- data_thinned_5km[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  

# subset dataset to match only thinned occurrences with maximum number of localities
x1 <- s_i_data[as.numeric(rownames(maxThin)),]  

# Export file 
write.csv(x1, "F:/Connectivity/outputs/02_SDMs/Intentos/sensitive data/FINAL_THINNED_DATA_SpizaetusIsidori.csv", row.names = FALSE)

# Read csv
x1 <- read.csv("F:/Connectivity/outputs/02_SDMs/Intentos/sensitive data/FINAL_THINNED_DATA_SpizaetusIsidori.csv")


x1 <- subset(x1, select = -c(SUBSPECIES.COMMON.NAME))
x1 <- x1 %>%  clean_names()

# Define name of the file that is going to be exported 
path_points <- paste0("F:/Connectivity/outputs/02_SDMs/Intentos/sensitive data/", x1$name[1], "_thinned_5km.shp")

# Turn to shapefile 
sp::coordinates(x1) <- ~longitude+latitude
sp::proj4string(x1) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Export
rgdal::writeOGR(obj=x1, dsn="F:/Connectivity/outputs/02_SDMs/Intentos/sensitive data/H", layer="Spizaetus isidori_thinned_5km", driver="ESRI Shapefile", overwrite_layer = TRUE)

  
# Minimum convex polygon with a buffer of 0.1°
  
# Extract only the coordinates
x_mcp <- x1[c('longitude', 'latitude')]
# Assign the columns with the coordinates
sp::coordinates(x_mcp) <- ~ longitude + latitude
# Calculate minimum convex polygon
calibration_Area <- mcp(x_mcp, percent = 100)
# Apply buffer 
calibration_Area <- rgeos::gBuffer(calibration_Area, width = 0.1)
# Assign coordinates
crs(calibration_Area) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
# Define name of the file that is going to be exported 
path_mcp <- paste0("F:/Connectivity/outputs/02_SDMs/Intentos/sensitive data/", x1$name[1], "_mcp_10km_thin5km.shp")
  
# Export 
raster::shapefile(calibration_Area, filename = path_mcp)

  


















# Turn into shapefile 
s_data_p <- sf::st_as_sf(s_data, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
sf::st_write(s_data_p, "F:/Connectivity/outputs/02_SDMs/Intentos/sensitive data/SpizaetusIsidori_notThinned.shp")



