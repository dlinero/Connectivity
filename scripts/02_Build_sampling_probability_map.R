#### CREATION OF SAMPLING PROBABILITY MAP ###

library(tidyverse)
library(sf)
library(raster)
library(reticulate)
library(RPyGeo)
library(auk)
library(arcgisbinding)

# Edit PATH environment variable, so that it includes the ArcGIS bin directory
usethis::edit_r_environ(scope = "project")

# Inside of the file specify: PATH=F:/R/R-4.1.1/bin/x64;C:/Miniconda;C:/Miniconda/Library/mingw-w64/bin;C:/Miniconda/Library/usr/bin;C:/Miniconda/Library/bin;C:/Miniconda/Scripts;C:/Miniconda/bin;C:/Miniconda/condabin;C:/Windows/system32;C:/Windows;C:/Windows/System32/Wbem;C:/Windows/System32/WindowsPowerShell/v1.0;C:/Windows/System32/OpenSSH;C:/Windows/system32/config/systemprofile/AppData/Local/Microsoft/WindowsApps;C:/Users/builduser/AppData/Local/Microsoft/WindowsApps;.;C:/Program Files/Azure Data Studio/bin;C:/Users/builduser/AppData/Roaming/npm;C:/Users/builduser/.dotnet/tools;C:/Program Files/Docker;C:/Program Files/Microsoft MPI/Bin;C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.1/bin;C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.1/libnvvp;C:/Program Files (x86)/Microsoft SDKs/Azure/CLI2/wbin;C:/Program Files/Microsoft/jdk-11.0.11.9-hotspot/bin;C:/ProgramData/Chocolatey/bin;C:/Program Files/dotnet;C:/Program Files/CMake/bin;C:/Program Files (x86)/Graphviz2.38/bin;C:/Julia/bin;C:/Program Files/Pandoc;C:/Program Files/Microsoft VS Code/bin;C:/Program Files/nodejs;C:/Program Files (x86)/Microsoft SQL Server/150/DTS/Binn;C:/Program Files (x86)/dotnet;C:/Program Files/NVIDIA Corporation/Nsight Compute 2020.2.0;C:/Program Files (x86)/NVIDIA Corporation/PhysX/Common;C:/NVIDIA/Cuda;C:/Program Files/NVIDIA Corporation/NVSMI;C:/hadoop/bin;C:/dsvm/tools/spark-3.1.1-bin-hadoop2.7/bin;C:/Program Files/VowpalWabbit;C:/dsvm/tools/DataMovement/ADL;C:/Program Files/Microsoft SQL Server/130/Tools/Binn;C:/Program Files/Microsoft SQL Server/Client SDK/ODBC/170/Tools/Binn;C:/Program Files (x86)/Microsoft SQL Server/150/Tools/Binn;C:/Program Files/Microsoft SQL Server/150/Tools/Binn;C:/Program Files/Microsoft SQL Server/150/DTS/Binn;C:/Program Files/R/R-4.1.0/bin;C:/Program Files/Git/bin;C:/Users/dlinero/AppData/Roaming/Python/Scripts;C:/Program Files (x86)/dotnet/;C:/Users/dlinero/AppData/Local/Microsoft/WindowsApps;C:/Users/dlinero/AppData/Roaming/TinyTeX/bin/win32;C:/Program Files/ArcGIS/Pro/bin
# save and open R again

# Connect to an arcgis licence 
arc.check_product()

# Initialize the Python ArcPy site-package in R
arcpy <- rpygeo_build_env(path = "C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/python.exe", 
                          overwrite = TRUE,
                          extensions = c("3d", "Spatial", "na"),
                          pro = TRUE, workspace = "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model")

# Progress bar for raster operations 
rasterOptions(progress = 'text',timer=TRUE)

# Birds occurrences  -----------------------------------------
# Input data - Costa Rica 
occ_path <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_CR_relOct-2021/ebd_CR_relOct-2021.txt"

# Ouput path 
f_out <- "Costa_Rica.txt"

ebird_data <- occ_path %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Costa_Rica <- ebird_data %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Costa_Rica <- Costa_Rica %>%
  distinct(longitude, latitude, .keep_all = TRUE)

# Input data - Panama 
occ_path2 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_PA_relJun-2020/ebd_PA_relJun-2020.txt"

# Ouput path 
f_out2 <- "Panama.txt"

ebird_data2 <- occ_path2 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out2) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Panama <- ebird_data2 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Panama <- Panama %>%
  distinct(longitude, latitude, .keep_all = TRUE)

# Input data - Colombia 
occ_path3 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_CO_relJun-2020/ebd_CO_relJun-2020.txt"

# Ouput path 
f_out3 <- "Colombia2.txt"

ebird_data3 <- occ_path3 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out3) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Colombia <- ebird_data3 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Colombia <- Colombia %>%
  distinct(longitude, latitude, .keep_all = TRUE)

# Input data - Venezuela  
occ_path4 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_VE_relOct-2021/ebd_VE_relOct-2021.txt"

# Ouput path 
f_out4 <- "Venezuela.txt"

ebird_data4 <- occ_path4 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out4) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Venezuela <- ebird_data4 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Venezuela <- Venezuela %>%
  distinct(longitude, latitude, .keep_all = TRUE)


# Input data - Ecuador  
occ_path5 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_EC_relOct-2021/ebd_EC_relOct-2021.txt"

# Ouput path 
f_out5 <- "Ecuador.txt"

ebird_data5 <- occ_path5 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out5) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Ecuador <- ebird_data5 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Ecuador <- Ecuador %>%
  distinct(longitude, latitude, .keep_all = TRUE)

# Input data - Peru 
occ_path6 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_PE_relMay-2021/ebd_PE_relMay-2021.txt"

# Ouput path 
f_out6 <- "Peru.txt"

ebird_data6 <- occ_path6 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out6) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Peru <- ebird_data6 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Peru <- Peru %>%
  distinct(longitude, latitude, .keep_all = TRUE)


# Input data - Brazil
occ_path7 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_BR_relOct-2021/ebd_BR_relOct-2021.txt"

# Ouput path 
f_out7 <- "Brazil.txt"

ebird_data7 <- occ_path7 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out7) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Brazil <- ebird_data7 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Brazil <- Brazil %>%
  distinct(longitude, latitude, .keep_all = TRUE)

# Input data - Guyana
occ_path8 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_GY_relOct-2021/ebd_GY_relOct-2021.txt"

# Ouput path 
f_out8 <- "Guyana.txt"

ebird_data8 <- occ_path8 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out8) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Guyana <- ebird_data8 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Guyana <- Guyana %>%
  distinct(longitude, latitude, .keep_all = TRUE)

# Input data - Suriname
occ_path9 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_SR_relOct-2021/ebd_SR_relOct-2021.txt"

# Ouput path 
f_out9 <- "Suriname.txt"

ebird_data9 <- occ_path9 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out9) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Suriname <- ebird_data9 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Suriname <- Suriname %>%
  distinct(longitude, latitude, .keep_all = TRUE)

# Input data - French Guiana
occ_path10 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_GF_relOct-2021/ebd_GF_relOct-2021.txt"

# Ouput path 
f_out10 <- "FGuiana.txt"

ebird_data10 <- occ_path10 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out10) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
FGuiana <- ebird_data10 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
FGuiana <- FGuiana %>%
  distinct(longitude, latitude, .keep_all = TRUE)

# Input data - Bolivia
occ_path11 <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/ebd_BO_relOct-2021/ebd_BO_relOct-2021.txt"

# Ouput path 
f_out11 <- "Bolivia.txt"

ebird_data11 <- occ_path11 %>% 
  # 1. reference file
  auk_ebd() %>% 
  # 2. define filters
  auk_date(date = c("2000-01-01", "2021-12-31")) %>% 
  auk_protocol(protocol = c("Traveling", "Stationary", "Area")) %>%
  # 3. run filtering
  auk_filter(file = f_out11) %>% 
  # 4. read text file into r data frame
  read_ebd()

# Delete records with routes longer than 3 km and periods longer than 4 hours
Bolivia <- ebird_data11 %>% 
  
  filter((!is.na(duration_minutes) & duration_minutes <= 240) | (!is.na(effort_distance_km) & effort_distance_km <= 3))

# Leave only one record per coordinate combination
Bolivia <- Bolivia %>%
  distinct(longitude, latitude, .keep_all = TRUE)

## Compile all files 
Sampling_sites <- rbind(Bolivia, Brazil, Costa_Rica, Ecuador, FGuiana, Guyana, Peru, Suriname, Venezuela) %>%
  select(scientific_name, longitude, latitude, country_code, protocol_type, duration_minutes, effort_distance_km)

Sampling_sites2 <- rbind(Colombia, Panama) %>% 
  select(scientific_name, longitude, latitude, country_code, protocol_type, duration_minutes, effort_distance_km)

Sampling_sites <- rbind(Sampling_sites, Sampling_sites2)

# Final filtering 
Sampling_sites <- Sampling_sites %>% filter(duration_minutes <= 240) %>%
  filter(effort_distance_km <= 3)

# Export 
write.csv(Sampling_sites, "F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv", row.names = FALSE)

# Import
Sampling_sites <- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")


# Using human footprint map -------------------------------------------------------

## Prepare predictor layers -----------------------------------------------------------
# Crop human footprint to study extent:

# Load polygon to crop all data (goes from Costa Rica to Peru)
st_area <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/stArea_Nicaragua_to_Peru.shp")

# Load raster
hfootprint <- raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/doi_10.5061_dryad.052q5__v2/Dryadv3/Maps/HFP2009_WGS84.tif")

# Ensure both are in the same coordinate system 
st_area <- st_transform(st_area, 4326)

plot(hfootprint)
plot(st_area, add = TRUE)

# crop
hfootprint_crop <- mask(hfootprint, st_area) %>% crop(extent(st_area))

# Export as acii
writeRaster(hfootprint_crop, filename="F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/predictors/Human_footprint2009_crop.asc", format="ascii", overwrite = TRUE)

# Rasterize countries: 

# Load shapefile: 
countries <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/doi_10.5061_dryad.052q5__v2/Dryadv3/Maps/Latin_America_countries.shp")

# Create ID for each country
countries$ID <- 1:nrow(countries)

# Rasterize 
countries_raster <- rasterize(countries, hfootprint_crop, field="ID")

# Defining raster layer as categorical 
countries_raster_classes <- ratify(countries_raster)
classes <- data.frame(countries) %>% dplyr::select(ID, PAÍS) %>% filter(ID %in% unique(countries_raster$layer)) %>%
  rename(country = PAÍS )
levels(countries_raster_classes) <- classes

# plot
rasterVis::levelplot(countries_raster_classes)

# Export
writeRaster(countries_raster_classes, filename="F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/predictors/countries_crop.asc", format="ascii", overwrite = TRUE)

## Prepare sampling points -----------------------------------------------------------

# Load data
Sampling_sites <- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to vector
Sampling_sites <- sf::st_as_sf(Sampling_sites, coords = c("longitude", "latitude"), crs = 4326)


# Load cropped countries and HF rasters (turned into polygon in ArcGIS for processing speed)
countries <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/countries_cropped.shp")
HF <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/HF_polygon_cropped.shp")

# Don't use spherical geometry to avoid warning errors durin intersection
sf_use_s2(FALSE)

# Remove points that fall outside the study extent 
points_HF_model <- st_intersection(Sampling_sites, countries)
points_HF_model <- st_intersection(points_HF_model, HF)

# Organize dataset
coordinates <- as.data.frame(st_coordinates(points_HF_model)) %>% rename(longitude = X, latitude = Y)
points_HF_model <- as.data.frame(points_HF_model) %>% rename(species = scientific_name) %>% dplyr::select(species) %>% cbind(coordinates) %>%
  mutate(species = "aves")

# Export
write.csv(points_HF_model, "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/aves_occurrences_HF_model.csv", row.names = FALSE)

## Run model  -----------------------------------------------------------


setwd("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model")
wd <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/predictors"
maxent.location <- "F:/R/R-4.1.1/library/dismo/java/maxent.jar"
model.output='output'

all.models.cv <-" nowarnings noprefixes responsecurves jackknife outputformat=raw noaskoverwrite removeduplicates=false -a -z replicates=5 nothreshold nohinge noautofeature replicatetype=bootstrap randomtestpoints=50 betamultiplier=0 "

#Modelo de muestreo total
system(paste0("java -mx2000m -jar ", maxent.location, all.models.cv,"outputdirectory=F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/output",
              " environmentallayers=",wd," -t countries_crop", ' samplesfile=F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/aves_occurrences_HF_model.csv',
              " maximumbackground=",nrow(read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/aves_occurrences_HF_model.csv"))))



# Distance based on human footprint map -----------------------------------------

# Data on human pressures in 1993 and 2009 were collected or developed for: 
# 1) the extent of built environments, 2) population density,
# 3) electric infrastructure, 4) crop lands, 5) pasture lands, 6) roads, 
# 7) railways, and 8) navigable waterways (https://www.nature.com/articles/sdata201667). 

# Load polygon to crop all data (goes from Costa Rica to Peru)
st_area <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/stArea_Nicaragua_to_Peru.shp")

# Re-project to World Equidistant Cylindrical crs 
st_area_equidistant <- st_transform(st_area, crs = 4087)

# Create a sample raster, so all predictors have the same coordinates, extent, resolution, etc. 
sample_raster <- raster("F:/Connectivity/outputs/02_SDMs/wc/30s_bio1_LAC.tif")

# Crop 
sample_raster <- mask(sample_raster, st_area) %>% crop(extent(st_area))

# Export
writeRaster(sample_raster, "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Sample_Raster.tif")



## Extent of built environments --------------------------------------------------------

# Built environments are human produced areas that provide the setting for 
# human activity. In the context of the human footprint, we take these areas
# to be primarily urban settings, including buildings, paved land and urban 
# parks. Built environments do not provide viable habitats for many species
# of conservation concern, nor do they provide high levels of ecosystem 
# services24–27. 

# To map built environments, we used the Defence Meteorological Satellite 
# Program Operational Line Scanner (DMSP-OLS) composite images which gives
# the annual average brightness of 30 arc second (~1 km at the equator)
# pixels in units of digital numbers (DN)28. These data are provided for
# each year from 1992 to 2012. We extracted data for the years 1994 
# (1993 was excluded due to anomalies in the data), and 2009, and both 
# datasets were then inter-calibrated to facilitate comparison29.
# Using the DMSP-OLS datasets, we considered pixels to be ‘built’
# if they exhibited a calibrated DN greater than 20. We selected this 
# threshold based on a global analyses of the implications of a range of 
# thresholds for mapped extent of cities30, and visual validation against
# Landsat imagery for 10 cities spread globally.

# Load raster
built <- raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/doi_10.5061_dryad.052q5__v2/Dryadv3/Maps/Built2009.tif")

# Re-project
built_equidistant <- projectRaster(built,
                                       crs = crs(st_area_equidistant), method = "bilinear")

# Crop 
built_equidistant_crop <- mask(built_equidistant, st_area_equidistant) %>% crop(extent(st_area_equidistant))
  

# Set zeros to NAs
built_equidistant_crop[Which(built_equidistant_crop<10)]<-NA

# Export
writeRaster(built_equidistant_crop, "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/Built2009_crop.tif")

# Create distance raster (using python for processing speed)
rpygeo_search(search_term = "distance")

py_function_docs("arcpy$sa$EucDistance")

# Calculate distance raster
built_distance <- arcpy$sa$EucDistance(in_source_data = "data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/Built2009_crop.tif",  distance_method =  "PLANAR")

# Import result
built_distance_result <- rpygeo_load(str <- gsub("\\", "/", as.character(built_distance), fixed=TRUE))

# Turn to km
built_distance_kmresult <-built_distance_result/1000
                                     
# Snap
built_distance_kmresult <- projectRaster(built_distance_kmresult,
                                       crs = crs(sample_raster), method="bilinear")

built_distance_kmresult = resample(built_distance_kmresult, sample_raster, "bilinear")
   
built_distance_kmresult <- mask(built_distance_kmresult, sample_raster) %>% crop(extent(sample_raster))

# Export
writeRaster(built_distance_kmresult, filename="F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/predictors/DistanceTo_BuiltEnvironments.asc", format="ascii", overwrite = TRUE)

## Population density --------------------------------------------------------

# Many of the pressures humans impose on the environment are proximate to their location, 
# such as human disturbance, hunting and the persecution of non-desired species34.
# Moreover, even low-density human populations with limited technology and infrastructure
# developments can have significant impacts on biodiversity, as evidenced by the 
# widespread loss of various taxa, particularly mega fauna, following human colonization 
# of previously unpopulated areas35,36.

# Human population density was mapped using the Gridded Population of the World dataset 
# developed by the Centre for International Earth Science Information Network (CIESEN)37.
# The dataset provides a ~4 km×~4 km gridded summary of population census data for the 
# years 1990 and 2010, which we downscaled using bilinear sampling in ArcGIS 10.1 to 
# match the 1 km2 resolution of the other datasets.  For all locations with more than
# 1000 people·km−2, we assigned a pressure score of 10 (Table 2). For more sparsely
# populated areas with densities lower than 1000 people·km−2, we logarithmically 
# scaled the pressure score using, Human population density is scored in this way under the assumption
# that the pressures people induce on their local natural systems increase 
# logarithmically with increasing population density, and saturate at a level of 1000
# people per km2.
                                    
# Load raster
popdensity <- raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/doi_10.5061_dryad.052q5__v2/Dryadv3/Maps/Popdensity2010.tif")

# Snap
popdensity_result <- projectRaster(popdensity,
                                         crs = crs(sample_raster), method="ngb")

popdensity_result <- resample(popdensity_result, sample_raster, "ngb")

popdensity_result <- mask(popdensity_result, sample_raster) %>% crop(extent(sample_raster))

# Export
writeRaster(popdensity_result, filename="F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/predictors/PopDensity_2010.asc", format="ascii", overwrite = TRUE)

## Electric infrastructure --------------------------------------------------------

# The high sensitivity of the DMSP-OLS28 dataset provides a means for mapping the sparser
# electric infrastructure typical of more rural and suburban areas. In 2009, 79% 
# of the lights registered in the DMSP-OLS dataset had a Digital Number less than 20,
# and are therefore not included in our ‘built environments’ layers. However, these 
# lower DN values are often important human infrastructures, such as rural housing or 
# working landscapes, with associated pressures on natural environments.

# To include these pressures, we used the inter-calibrated DMSP-OLS layers28 used 
# for the built environments mapping. The equations for intercalibrating across years 
# are second order quadratics trained using data from Sicily, which was chosen as it
# had negligible infrastructure change over this period and where DN average roughly 
# 14 (ref. 28). For our purposes, DN values of six or less where excluded from 
# consideration prior to calibration of data, as the shape of the quadratic function 
# leads to severe distortion of very low DN values. The inter-calibrated DN data from 
# 1994 were then rescaled using an equal quintile approach into a 0–10 scale (Table 2). 
# To scale the data, we divided the calibrated night light data into 10 equal sample
# bins (each bin with a DN greater than 1 contains the same number of pixels) based on 
# the DN values and then assigned them scores of 1 through 10, starting with the lowest
# DN bin. DN values of 0 were assigned a score of 0. The thresholds used to bin the 1994 
# data where then used to convert the 2009 data into a comparable 0–10 scale.

# Load raster
nightlights <- raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/doi_10.5061_dryad.052q5__v2/Dryadv3/Maps/Lights2009.tif")

# Snap
nightlights_result <- projectRaster(nightlights,
                                   crs = crs(sample_raster), method="ngb")

nightlights_result <- resample(nightlights_result, sample_raster, "ngb")

nightlights_result <- mask(nightlights_result, sample_raster) %>% crop(extent(sample_raster))

# Export
writeRaster(nightlights_result, filename="F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/predictors/ElectricInfrastructure_2009.asc", format="ascii", overwrite = TRUE)

## Roads --------------------------------------------------------

# As one of humanity’s most prolific linear infrastructures, roads are an
# important direct driver of habitat conversion44. Beyond simply reducing the extent 
# of suitable habitat, roads can act as population sinks for many species through traffic
# induced mortality45. Roads also fragment otherwise contiguous blocks of habitat, 
# and create edge effects such as reduced humidity6 and increased fire frequency that 
# reach well beyond the roads immediate footprint46. Finally, roads provide conduits for
# humans to access nature, bringing hunters and nature users into otherwise wilderness 
# locations47."

# We acquired data on the distribution of roads from gROADS48, and excluded all trails
# and private roads, which were inconsistently mapped, with only a subset of countries
# mapping their linear infrastructure to this resolution. The dataset is the most 
# comprehensive publicly available database on roads, which compiles nationally mapped
# road data spanning the period 1980–2000 and has a spatial accuracy of around 500 m. 
# The gROADS data do not include all minor roads, and therefore should be viewed as a map
# of the major roadways. We mapped the direct and indirect influence of roads by 
# assigning an pressure score of 8 for 0.5 km out for either side of roads, and access 
# pressures were awarded a score of 4 at 0.5 km and decaying exponentially out to 15 km
# either side of the road (Table 2).

# Since the authors assigned a value of 8 to cells 0.5km out for either side of roads, I will 
# extract the roads by extracting cells with a score of 8. This buffer also allows to account for topography errors. 

# Load raster
roads <- raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/doi_10.5061_dryad.052q5__v2/Dryadv3/Maps/roads.tif")

# Set to NA all values less than 8
roads[Which(roads<8)]<-NA

# Re-project
roads_equidistant <- projectRaster(roads,
                                   crs = crs(st_area_equidistant))

# Crop 
roads_equidistant_crop <- mask(roads_equidistant, st_area_equidistant) %>% crop(extent(st_area_equidistant))

# Export
writeRaster(roads_equidistant_crop, "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/Roads2000_crop.tif")

# Create distance raster (using python for processing speed)
rpygeo_search(search_term = "distance")

py_function_docs("arcpy$sa$EucDistance")

# Calculate distance raster
roads_distance <- arcpy$sa$EucDistance(in_source_data = "data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/Roads2000_crop.tif",  distance_method =  "PLANAR")

# Import result
roads_distance_result <- rpygeo_load(str <- gsub("\\", "/", as.character(roads_distance), fixed=TRUE))

# Turn to km
roads_distance_kmresult <-roads_distance_result/1000

# Snap
roads_distance_kmresult <- projectRaster(built_distance_kmresult,
                                         crs = crs(sample_raster))

roads_distance_kmresult = resample(roads_distance_kmresult, sample_raster, "bilinear")

roads_distance_kmresult <- mask(roads_distance_kmresult, sample_raster) %>% crop(extent(sample_raster))

# Export
writeRaster(roads_distance_kmresult, filename="F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/predictors/DistanceTo_Roads.asc", format="ascii", overwrite = TRUE)

## Navigable waterways --------------------------------------------------------

# Like roads, coastlines and navigable rivers act as conduits for people to access
# nature. While all coastlines are theoretically navigable, for the purposes of the 
# human footprint we only considered coasts31 as navigable for 80 km either direction
# of signs of a human settlement, which were mapped as a night lights signal with a
# DN28 greater than 6 within 4 km of the coast. We chose 80 km as an approximation 
# of the distance a vessel can travel and return during daylight hours. 
# As new settlements can arise to make new sections of coast navigable, coastal layers
# were generated for the years 1994 and 2009.

# Large lakes can act essentially as inland seas, with their coasts frequently plied 
# by trade and harvest vessels. Based on their size and visually identified shipping
# traffic and shore side settlements, we treated the great lakes of North America,
# Lake Nicaragua, Lake Titicaca in South America, Lakes Onega and Peipus in Russia
# , Lakes Balkash and Issyk Kul in Kazakhstan, and Lakes Victoria, Tanganyika and Malawi
# in Africa as we did navigable marine coasts.

# Rivers were considered as navigable if their depth was greater than 2 m and
# there were signs of nighttime lights (DN>=6) within 4km of their banks, or if 
# contiguous with a navigable coast or large inland lake, and then for a distance 
# of 80 km or until stream depth is likely to prevent boat traffic (Table 2). To map 
# rivers and their depth we used the hydrosheds (hydrological data and maps based on
# shuttle elevation derivatives at multiple scales)49 dataset on stream discharge,
# and the following formulae from50,51:

# Navigable rivers layers were created for the years 1994 and 2009, and combined
# with the navigable coasts and inland seas layers to create the final navigable 
# waterways layers. The access pressure from navigable water bodies were awarded 
# a score of 4 adjacent to the water body, decaying exponentially out to 15 km.

# Since the authors assigned a value of 4 to  the navigable water bodies , I will 
# extract the rivers by extracting cells with a score of 4.  

# Load raster
waterways <- raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/doi_10.5061_dryad.052q5__v2/Dryadv3/Maps/NavWater2009.tif")

# Set to NA all values less than 8
waterways[Which(waterways<4)]<-NA

# Re-project
waterways_equidistant <- projectRaster(waterways,
                                   crs = crs(st_area_equidistant))

# Crop 
waterways_equidistant_crop <- mask(waterways_equidistant, st_area_equidistant) %>% crop(extent(st_area_equidistant))

# Export
writeRaster(waterways_equidistant_crop, "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/Waterways2009_crop.tif")

# Create distance raster (using python for processing speed)
rpygeo_search(search_term = "distance")

py_function_docs("arcpy$sa$EucDistance")

# Calculate distance raster
waterways_distance <- arcpy$sa$EucDistance(in_source_data = "data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/Waterways2009_crop.tif",  distance_method =  "PLANAR")

# Import result
waterways_distance_result <- rpygeo_load(str <- gsub("\\", "/", as.character(waterways_distance), fixed=TRUE))

# Turn to km
waterways_distance_kmresult <-waterways_distance_result/1000

# Snap
waterways_distance_kmresult <- projectRaster(waterways_distance_kmresult,
                                         crs = crs(sample_raster), method = "bilinear")

waterways_distance_kmresult = resample(waterways_distance_kmresult, sample_raster, "bilinear")

waterways_distance_kmresult <- mask(waterways_distance_kmresult, sample_raster) %>% crop(extent(sample_raster))

# Export
writeRaster(waterways_distance_kmresult, filename="F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/predictors/DistanceTo_Waterways.asc", format="ascii", overwrite = TRUE)

## Distance to protected areas  --------------------------------------------------------

# This one was not included in the study, so it is an additional layer (IUCN categories 1-4)

# Load shapefile: 
PAs <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/WDPA_1to4_crop.shp")

# Rasterize 
PAs_raster <- rasterize(PAs, sample_raster, field="ID")

# plot
plot(PAs_raster)

# Re-project
PAs_equidistant <- projectRaster(PAs_raster,
                                       crs = crs(st_area_equidistant))
# Export
writeRaster(PAs_equidistant, "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/PAs2021_raster_crop.tif")

# Create distance raster (using python for processing speed)
rpygeo_search(search_term = "distance")

py_function_docs("arcpy$sa$EucDistance")

# Calculate distance raster
PAs_distance <- arcpy$sa$EucDistance(in_source_data = "data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/PAs2021_raster_crop.tif",  distance_method =  "PLANAR")

# Import result
PAs_distance_result <- rpygeo_load(str <- gsub("\\", "/", as.character(PAs_distance), fixed=TRUE))

# Turn to km
PAs_distance_kmresult <-PAs_distance_result/1000

# Snap
PAs_distance_kmresult <- projectRaster(PAs_distance_kmresult,
                                             crs = crs(sample_raster), method = "bilinear")

PAs_distance_kmresult = resample(PAs_distance_kmresult, sample_raster, "bilinear")

PAs_distance_kmresult <- mask(PAs_distance_kmresult, sample_raster) %>% crop(extent(sample_raster))

# Export
writeRaster(PAs_distance_kmresult, filename="F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/predictors/DistanceTo_WDPA.asc", format="ascii", overwrite = TRUE)

## Countries  --------------------------------------------------------

# Load shapefile: 
countries <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/doi_10.5061_dryad.052q5__v2/Dryadv3/Maps/Latin_America_countries.shp")

# Create ID for each country
countries$ID <- 1:nrow(countries)

# Rasterize 
countries_raster <- rasterize(countries, sample_raster, field="ID")

# Export
writeRaster(countries_raster, filename="F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/predictors/countries_crop.asc", format="ascii", overwrite = TRUE)



## Prepare sampling points -----------------------------------------------------------

# Load data
Sampling_sites <- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to vector
Sampling_sites <- sf::st_as_sf(Sampling_sites, coords = c("longitude", "latitude"), crs = 4326)

# Load cropped countries and sample raster (turned into polygon in ArcGIS for processing speed)
countries <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/countries_cropped.shp")
sample_raster_shp <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Sample_Raster.shp")
electric_shp <- st_read("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Electric_infrastructure.shp")

# Don't use spherical geometry to avoid warning errors during intersection
sf_use_s2(FALSE)

# Remove points that fall outside the study extent 
points_HF_distance_model <- st_intersection(Sampling_sites, countries)
points_HF_distance_model <- st_intersection(points_HF_distance_model, sample_raster_shp)
points_HF_distance_model <- st_intersection(points_HF_distance_model, electric_shp)

# Organize dataset
coordinates <- as.data.frame(st_coordinates(points_HF_distance_model)) %>% rename(longitude = X, latitude = Y)
points_HF_distance_model <- as.data.frame(points_HF_distance_model) %>% rename(species = scientific_name) %>% dplyr::select(species) %>% cbind(coordinates) %>%
  mutate(species = "aves")

# Export
write.csv(points_HF_distance_model, "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/aves_occurrences_HF_distance_model.csv", row.names = FALSE)

## Run model  -----------------------------------------------------------


setwd("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model")
wd <- "F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/predictors"
maxent.location <- "F:/R/R-4.1.1/library/dismo/java/maxent.jar"
model.output='output'

all.models.cv <-" nowarnings noprefixes responsecurves jackknife outputformat=raw noaskoverwrite removeduplicates=false -a -z replicates=5 nothreshold nohinge noautofeature replicatetype=bootstrap randomtestpoints=50 betamultiplier=0 "

#Modelo de muestreo total
system(paste0("java -mx2000m -jar ", maxent.location, all.models.cv,"outputdirectory=F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output",
              " environmentallayers=",wd, " -t countries_crop", ' samplesfile=F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/aves_occurrences_HF_distance_model.csv',
              " maximumbackground=",nrow(read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/aves_occurrences_HF_distance_model.csv"))))


system(paste0("java -mx2000m -jar ", maxent.location, all.models.cv,"outputdirectory=F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/output",
              " environmentallayers=",wd," -t countries_crop", ' samplesfile=F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/aves_occurrences_HF_model.csv',
              " maximumbackground=",nrow(read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/aves_occurrences_HF_model.csv"))))



# Questions --------------------------------------------------------------------------


# Preguntarle a Jorge lo de dividir el mapa, cuando uno hace lo de muestrear 
# el sampling probability map para el background --> https://github.com/jamiemkass/ENMeval/issues/26 (comentario de bob muscarella)
# Si calcule bien el 10 percentil?