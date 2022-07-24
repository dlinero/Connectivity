library(sf)
library(raster)
library(tidyverse)
library(ncdf4)
library(rgeos)
library(lwgeom)

# Disable scientific notation
options(scipen=999)
rasterOptions(progress = 'text',timer=TRUE)

# Cropped selected best models to biogeographical boundaries drawn by Jorge Velásquez -------------------------------------------------

# Load list of final models
final_models <- list.files("./outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/FINAL_selected models", 
                           pattern = ".tif$",
                           full.names = TRUE, 
                           recursive =  FALSE)

# Load polygons that define model boundaries 
load_shapefile <- st_read("./outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Distribution_boundaries_velasquezWGS84.shp")


for (i in 1:length(final_models)){
  
  print(i)
  
  # Load raster
  raster <- raster::raster(final_models[i])
  
  # Extract species name 
  name <- str_replace_all(final_models[i], "[^[:alnum:]]", " ") %>% word(start=16,end=17,sep=fixed(" "))
  
  # Select respective shapefile 
  sp_boundary <-  load_shapefile[load_shapefile$sp_name == name, ] 
  
  
  # Clip
  crop <- raster::mask(x = raster, mask = sp_boundary) %>% raster::crop(extent(sp_boundary))
  
  # Export
  path1 <- "./outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/FINAL_selected models/CROPPED MODELS/"
  path2 <- paste0((str_replace_all(final_models[i], "[^[:alnum:]]", " ") %>% word(start=16,end=17,sep=fixed(" "))), "_", "cropped", "_")
  path3 <- str_replace_all(final_models[i], "[^[:alnum:]]", " ") %>% word(start=19,end=19,sep=fixed(" "))
  path4 <- ".tif"
  path5 <- paste0(path1, path2, path3, path4)
  raster::writeRaster(crop, path5, overwrite = TRUE)
}


# Remove cells with less that 50% forest cover - ESA  -----------------------------------------------------

# Load data
nc = nc_open("./data/03_connectivity/ESA/ESA_2020.nc")
# Get LCC variable
d_2020 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2020), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load Colombia polygon 
colombia <- sf::st_read("./data/Colombia_shapefile .shp")

# Clip raster
data_crop <- raster::mask(x = r, mask = colombia) %>% raster::crop(extent(colombia))

# Remove heavy variables
rm(d_2020, nc, r)

# Get unique land cover classes
unique(data_crop)

# Merge some land uses to decrease the number of classes
reclass_df <- c(50, 1, 
                60, 1, 
                90, 1, 
                160, 1, 
                10, 0, 
                11, 0, 
                12, 0, 
                30, 0,
                40, 0, 
                61, 1, 
                80, 0, 
                100, 0, 
                110, 0,
                120, 0, 
                122, 0, 
                130, 0, 
                150, 0, 
                153, 0,
                170, 0, 
                180, 0, 
                190, 0, 
                200, 0, 
                210, 0,
                220, 0)


# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
forest <-reclassify(data_crop, rcl=reclass_df)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(forest, fact=3)

# Calculate %forest per 1 km cell
forest_fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
  sum(vals==1, na.rm=na.rm)/length(vals)
})

# Set suitability to zero for percents lower than 50
forest_matrix <- cbind( from = c(0, 0.5), to = c(0.5, 1), becomes = c(0, 1))
forest_reclass <- reclassify(forest_fraction, forest_matrix)

# # Snap raster to bioclimatic variables
wc <- raster("./outputs/02_SDMs/wc/30s_bio1_LAC.tif")
wc_crop <- raster::mask(x = wc, mask = colombia) %>% raster::crop(extent(colombia))
forest_final <- raster::resample(forest_reclass, wc_crop, "ngb")

# Export
raster::writeRaster(forest_final, "./outputs/03_connectivity/ESA/percentsForest_larger50_2020.tif", overwrite = TRUE)

# Load list of final models
final_models <- list.files("./outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/FINAL_selected models/CROPPED MODELS", 
                           pattern = ".tif$",
                           full.names = TRUE, 
                           recursive =  FALSE)


# Load ESA forest larger than 50% per cell
forest_final<- raster::raster("./outputs/03_connectivity/ESA/percentsForest_larger50_2020.tif")


# Loop to crop HS models by % tree
for (i in 1:length(final_models)){
  
  print(i)
  
  x <- raster::raster(final_models[i])
  
  # Multiply
  raster_forest_50 <- x * forest_final
  
  # Export
  path1 <- "./outputs/03_connectivity/ESA/FINAL_CROPPED_MODELS_ESA/"
  path2 <- paste0((str_replace_all(final_models[i], "[^[:alnum:]]", " ") %>% word(start=18,end=19,sep=fixed(" "))), "_", "ESAforest50percent", "_")
  path3 <- str_replace_all(final_models[i], "[^[:alnum:]]", " ") %>% word(start=21,end=21,sep=fixed(" "))
  path4 <- ".tif"
  path5 <- paste0(path1, path2, path3, path4)
  raster::writeRaster(raster_forest_50, path5, overwrite = TRUE)
}

# Select core area per species  -----------------------------------------------------

# Load list of final cropped models
final_models <- list.files("./outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/FINAL_selected models/CROPPED MODELS", 
                           pattern = ".tif$",
                           full.names = TRUE, 
                           recursive =  FALSE)

# Load ESA forest larger than 50% per cell
forest_final<- raster::raster("./outputs/03_connectivity/ESA/percentsForest_larger50_2020.tif")


# Loop to prepare layers 
for (i in 1:length(final_models)){
  
  print(i)
  
  x <- raster::raster(final_models[i])
  
  min_hs <- cellStats(x, "min")
  
  # Reclassify all suitable values to 1
  m <- c((min_hs-1), 1, 1)
  m <- matrix(m, ncol=3, byrow=TRUE)
  climatically_suitable <- raster::reclassify(x, m)
  
  # Multiply
  climatically_habitat_suitable <- climatically_suitable * forest_final
  
  # Turn zeros to NAs 
  climatically_habitat_suitable[climatically_habitat_suitable == 0] <- NA 
  
  # Project
  final <- projectRaster(climatically_habitat_suitable, crs = "+proj=laea +lon_0=-72.6855469 +lat_0=0 +datum=WGS84 +units=m +no_defs", method = "ngb")
  
  # Export
  path1 <- "./outputs/03_connectivity/Core areas/Processing/"
  path2 <- paste0((str_replace_all(final_models[i], "[^[:alnum:]]", " ") %>% word(start=18,end=19,sep=fixed(" "))), "_", "only_climate_habitat_suitable")
  path3 <- ".tif"
  path4 <- paste0(path1, path2, path3)
  raster::writeRaster(climatically_habitat_suitable, path4, overwrite = TRUE)
}

# In ArcGIS I ran
# arcpy.Tmp.BatchInt(r"'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Ara macao_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Buteo albigula_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Campephilus pollens_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Crypturellus berlepschi_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Crypturellus erythropus_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Crypturellus soui_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Grallaria bangsi_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Grallaria flavotincta_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Grallaria hypoleuca_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Grallaria nuchalis_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Henicorhina leucophrys_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Micrastur semitorquatus_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Mitu tomentosum_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Myioborus flavivertex_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Odontophorus hyperythrus_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Ognorhynchus icterotis_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Patagioenas goodsoni_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Phaethornis malaris_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Picumnus cinnamomeus_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Picumnus squamulatus_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Pseudastur albicollis_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Pyrrhura melanura_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Ramphastos brevis_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Spizaetus isidori_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Spizaetus ornatus_only_climate_habitat_suitable.tif';'F:\Connectivity\outputs\03_connectivity\Core areas\Processing\Tangara johannae_only_climate_habitat_suitable.tif'", r"F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_%Name%")
# arcpy.Tmp.BatchProjectRaster(r"'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Ara_macao_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Buteo_albigula_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Campephilus_pollens_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Crypturellus_berlepschi_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Crypturellus_erythropus_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Crypturellus_soui_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Grallaria_bangsi_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Grallaria_flavotincta_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Grallaria_hypoleuca_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Grallaria_nuchalis_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Henicorhina_leucophrys_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Micrastur_semitorquatus_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Mitu_tomentosum_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Myioborus_flavivertex_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Odontophorus_hyperythrus_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Ognorhynchus_icterotis_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Patagioenas_goodsoni_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Phaethornis_malaris_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Picumnus_cinnamomeus_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Picumnus_squamulatus_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Pseudastur_albicollis_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Pyrrhura_melanura_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Ramphastos_brevis_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Spizaetus_isidori_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Spizaetus_ornatus_only_climate_habitat_suitable';'F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Int_Tangara_johannae_only_climate_habitat_suitable'", r"F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Pr_%Name%", 'PROJCS["WGS_1984_Lambert_Azimuthal_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-73.0371094],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]', "NEAREST", "924.435276945671 924.435276945671", None, None, 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]', "NO_VERTICAL")
# arcpy.management.Project("runap2_EqualAreaProjection", r"F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\runap_EqualArea", 'PROJCS["WGS_1984_Lambert_Azimuthal_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-73.0371094],PARAMETER["latitude_of_origin",0.0],UNIT["Meter",1.0]]', None, 'PROJCS["unknown",GEOGCS["GCS_unknown",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-72.6855469],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]', "PRESERVE_SHAPE", None, "NO_VERTICAL")
# arcpy.Tmp.BatchTabulateArea3("runap_EqualArea", "nombre", "Pr_Int_Ara_macao_only_climate_habitat_suitable;Pr_Int_Buteo_albigula_only_climate_habitat_suitable;Pr_Int_Campephilus_pollens_only_climate_habitat_suitable;Pr_Int_Crypturellus_berlepschi_only_climate_habitat_suitable;Pr_Int_Crypturellus_erythropus_only_climate_habitat_suitable;Pr_Int_Crypturellus_soui_only_climate_habitat_suitable;Pr_Int_Grallaria_bangsi_only_climate_habitat_suitable;Pr_Int_Grallaria_flavotincta_only_climate_habitat_suitable;Pr_Int_Grallaria_hypoleuca_only_climate_habitat_suitable;Pr_Int_Grallaria_nuchalis_only_climate_habitat_suitable;Pr_Int_Henicorhina_leucophrys_only_climate_habitat_suitable;Pr_Int_Micrastur_semitorquatus_only_climate_habitat_suitable;Pr_Int_Mitu_tomentosum_only_climate_habitat_suitable;Pr_Int_Myioborus_flavivertex_only_climate_habitat_suitable;Pr_Int_Odontophorus_hyperythrus_only_climate_habitat_suitable;Pr_Int_Ognorhynchus_icterotis_only_climate_habitat_suitable;Pr_Int_Patagioenas_goodsoni_only_climate_habitat_suitable;Pr_Int_Phaethornis_malaris_only_climate_habitat_suitable;Pr_Int_Picumnus_cinnamomeus_only_climate_habitat_suitable;Pr_Int_Picumnus_squamulatus_only_climate_habitat_suitable;Pr_Int_Pseudastur_albicollis_only_climate_habitat_suitable;Pr_Int_Pyrrhura_melanura_only_climate_habitat_suitable;Pr_Int_Ramphastos_brevis_only_climate_habitat_suitable;Pr_Int_Spizaetus_isidori_only_climate_habitat_suitable;Pr_Int_Spizaetus_ornatus_only_climate_habitat_suitable;Pr_Int_Tangara_johannae_only_climate_habitat_suitable", "Value", r"F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\TabArea_%Name%", r"F:\Connectivity\outputs\03_connectivity\Core areas\Default.gdb\Pr_Int_Ara_macao_only_climate_habitat_suitable", "CLASSES_AS_ROWS")

# Load  results - each species has a table showing the total suitable area per protected area
area_pas <- list.files("./outputs/03_connectivity/Core areas/Processing", 
                           pattern = ".csv$",
                           full.names = TRUE, 
                           recursive =  FALSE)

# Load polygons of protected areas 
pas <- st_read("./data/03_connectivity/runap2/runap2Polygon.shp")

# Load list of selected species with mean home range values 
home_range <- read.csv("./data/selected_species_FINAL.csv")

for (i in 1:length(area_pas)){
  
  print(i)
  
  data <- read.csv(area_pas[i])
  
  data <- data %>% mutate(Ha_suitable = Area/10000) %>% rename(objectid = objectid_1)
  
  complete_pas <- inner_join(pas, data, by = "objectid")
  
  species_name <- str_replace_all(area_pas[i], "[^[:alnum:]]", " ") %>% word(start=12,end=13,sep=fixed(" "))
  
  area_threshold <- home_range$Home.range..ha.[home_range$Species == species_name]
  
  complete_pas <- complete_pas %>% filter(Ha_suitable >= area_threshold) 
  
  path1 <- "F:/Connectivity/outputs/03_connectivity/Core areas/Core areas per species/"
  path2 <- paste0(species_name, "_core_areas", ".shp")
  path3 <- paste0(path1, path2)
  
  st_write(complete_pas, path3)
}


# Define study area per species ---------------------------------------------------------------------------

# Load list of final cropped models
final_cropped_models <- list.files("./outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/FINAL_selected models/CROPPED MODELS", 
                           pattern = ".tif$",
                           full.names = TRUE, 
                           recursive =  FALSE)


# Create table to store buffer length per species
buffers <- data.frame(matrix(NA,    # Create empty data frame
                                     nrow = 30,
                                     ncol = 6))

buffers <- buffers %>% rename(species = X1, Distribution_area_sqkm = X2, Distribution_perimeter_km = X3, 
                                   Buffer_distance_km = X4, Area_percent_increase = X6, Final_area_sqkm = X5)

for (i in 1:length(final_cropped_models)){
  
  print(i)
  
  # Load raster
  raster <- raster::raster(final_cropped_models[i])
  
  # Extract species name 
  buffers$species[i] <- str_replace_all(final_cropped_models[i], "[^[:alnum:]]", " ") %>% word(start=18,end=19,sep=fixed(" "))
  
  # Reclass all values to 1
  min_value <- cellStats(raster, "min")
  
  # Reclassify all suitable values to 1
  m <- c((min_value-0.1), 2, 1)
  m <- matrix(m, ncol=3, byrow=TRUE)
  raster <- raster::reclassify(raster, m)
  
  # Turn distribution into shapefile
  x <- rasterToPolygons(raster, na.rm = TRUE, dissolve = TRUE)
  
  # Reproject tp calculate distances and areas
  x <- as(x, "sf")
  x <- st_transform(x, 32618)
  
  # Calculate area of polygon
  buffers$Distribution_area_sqkm[i] <- st_area(x)/1000
  
  # Calculate perimeter of polygon
  buffers$Distribution_perimeter_km[i] <-  st_perimeter(x)/1000
  
  # Calculate buffer distance to achieve an increase of ~20% in polygon area
  buffers$Buffer_distance_km[i] <- (((buffers$Distribution_area_sqkm[i] * 1.4) - buffers$Distribution_area_sqkm[i]) / buffers$Distribution_perimeter_km[i])/1000
  
  # Apply buffer
  x <- as(x, Class = "Spatial")
  dataset_transform_buffer <- raster::buffer(x, width = (buffers$Buffer_distance_km[i]*1000), dissolve = TRUE)
  
  # Calculate actual % increase in area (due to irregular polygons)
  dataset_transform_buffer <- as(dataset_transform_buffer, "sf")
  buffers$Final_area_sqkm[i] <- st_area(dataset_transform_buffer)/1000
  buffers$Area_percent_increase[i] <- ((buffers$Final_area_sqkm[i] - buffers$Distribution_area_sqkm[i]) / buffers$Distribution_area_sqkm[i])*100

  
  # Export
  path1 <- "./outputs/03_connectivity/Species distributions with buffer/"
  path2 <- str_replace_all(final_cropped_models[i], "[^[:alnum:]]", " ") %>% word(start=18,end=19,sep=fixed(" "))
  path3 <- "_study_area_buffer.shp"
  path4 <- paste0(path1, path2, path3)
  
  sf::st_write(dataset_transform_buffer, path4)
}

# Dado que hay algunos en donde el aumento fue muy grande voy a volver a hacer esos 
# La idea es que el buffer este entre el 17-23%

species_revise_large <- buffers$species[which(buffers$Area_percent_increase >= 23)]


# manually erase the shp for these species
species_revise_large

# Indicate the species that we are going to re-run
final_cropped_models <- data.frame(final_cropped_models) %>% mutate(species_name = word(str_replace_all(final_cropped_models, "[^[:alnum:]]", " "), start = 18, end = 19, sep = fixed(" ")),
                                                                    new_buffer1 = case_when(species_name %in% species_revise ~ "Yes", 
                                                                                            TRUE ~ "No")) %>% filter(new_buffer1 == "Yes") %>% droplevels()


for (i in 20:nrow(final_cropped_models)){
  
  print(i)
  
  # Load raster
  raster <- raster::raster(final_cropped_models$final_cropped_models[i])
  
  # Extract species name 
  buffers_index <- which(buffers$species == final_cropped_models$species_name[i])
  
  # Reclass all values to 1
  min_value <- cellStats(raster, "min")
  
  # Reclassify all suitable values to 1
  m <- c((min_value-0.1), 2, 1)
  m <- matrix(m, ncol=3, byrow=TRUE)
  raster <- raster::reclassify(raster, m)
  
  # Turn distribution into shapefile
  x <- rasterToPolygons(raster, na.rm = TRUE, dissolve = TRUE)
  
  # Reproject tp calculate distances and areas
  x <- as(x, "sf")
  x <- st_transform(x, 32618)
  
  # Calculate area of polygon
  buffers$Distribution_area_sqkm[buffers_index] <- st_area(x)/1000
  
  # Calculate perimeter of polygon
  buffers$Distribution_perimeter_km[buffers_index] <-  st_perimeter(x)/1000
  
  # Calculate buffer distance to achieve an increase of ~20% in polygon area
  buffers$Buffer_distance_km[buffers_index] <- (((buffers$Distribution_area_sqkm[buffers_index] * 1.25) - buffers$Distribution_area_sqkm[buffers_index]) / buffers$Distribution_perimeter_km[buffers_index])/1000
  
  # Apply buffer
  x <- as(x, Class = "Spatial")
  dataset_transform_buffer <- raster::buffer(x, width = (buffers$Buffer_distance_km[buffers_index]*1000), dissolve = TRUE)
  
  # Calculate actual % increase in area (due to irregular polygons)
  dataset_transform_buffer <- as(dataset_transform_buffer, "sf")
  buffers$Final_area_sqkm[buffers_index] <- st_area(dataset_transform_buffer)/1000
  buffers$Area_percent_increase[buffers_index] <- ((buffers$Final_area_sqkm[buffers_index] - buffers$Distribution_area_sqkm[buffers_index]) / buffers$Distribution_area_sqkm[buffers_index])*100
  
  
  # Export
  path1 <- "./outputs/03_connectivity/Species distributions with buffer/"
  path2 <- buffers$species[buffers_index]
  path3 <- "_study_area_buffer.shp"
  path4 <- paste0(path1, path2, path3)
  
  sf::st_write(dataset_transform_buffer, path4, overwrite = TRUE)
}

# Dado que hay algunos en donde el aumento fue muy pequeño voy a volver a hacer esos 
# La idea es que el buffer este entre el 17-23%

buffers <- read.csv("F:/Connectivity/outputs/03_connectivity/Species distributions with buffer/buffer_per_species.csv")

species_revise_small <- buffers$species[which(buffers$Area_percent_increase <= 17)]


# manually erase the shp for these species
species_revise_small

# Indicate the species that we are going to re-run
final_cropped_models <- data.frame(final_cropped_models) %>% mutate(species_name = word(str_replace_all(final_cropped_models, "[^[:alnum:]]", " "), start = 18, end = 19, sep = fixed(" ")),
                                                                    new_buffer1 = case_when(species_name %in% species_revise_small ~ "Yes", 
                                                                                            TRUE ~ "No")) %>% filter(new_buffer1 == "Yes") %>% droplevels()


for (i in 1:nrow(final_cropped_models)){
  
  print(i)
  
  # Load raster
  raster <- raster::raster(final_cropped_models$final_cropped_models[i])
  
  # Extract species name 
  buffers_index <- which(buffers$species == final_cropped_models$species_name[i])
  
  # Reclass all values to 1
  min_value <- cellStats(raster, "min")
  
  # Reclassify all suitable values to 1
  m <- c((min_value-0.1), 2, 1)
  m <- matrix(m, ncol=3, byrow=TRUE)
  raster <- raster::reclassify(raster, m)
  
  # Turn distribution into shapefile
  x <- rasterToPolygons(raster, na.rm = TRUE, dissolve = TRUE)
  
  # Reproject tp calculate distances and areas
  x <- as(x, "sf")
  x <- st_transform(x, 32618)
  
  # Calculate area of polygon
  buffers$Distribution_area_sqkm[buffers_index] <- st_area(x)/1000
  
  # Calculate perimeter of polygon
  buffers$Distribution_perimeter_km[buffers_index] <-  st_perimeter(x)/1000
  
  # Calculate buffer distance to achieve an increase of ~20% in polygon area
  buffers$Buffer_distance_km[buffers_index] <- (((buffers$Distribution_area_sqkm[buffers_index] * 1.6) - buffers$Distribution_area_sqkm[buffers_index]) / buffers$Distribution_perimeter_km[buffers_index])/1000
  
  # Apply buffer
  x <- as(x, Class = "Spatial")
  dataset_transform_buffer <- raster::buffer(x, width = (buffers$Buffer_distance_km[buffers_index]*1000), dissolve = TRUE)
  
  # Calculate actual % increase in area (due to irregular polygons)
  dataset_transform_buffer <- as(dataset_transform_buffer, "sf")
  buffers$Final_area_sqkm[buffers_index] <- st_area(dataset_transform_buffer)/1000
  buffers$Area_percent_increase[buffers_index] <- ((buffers$Final_area_sqkm[buffers_index] - buffers$Distribution_area_sqkm[buffers_index]) / buffers$Distribution_area_sqkm[buffers_index])*100
  
  
  # Export
  path1 <- "./outputs/03_connectivity/Species distributions with buffer/"
  path2 <- buffers$species[buffers_index]
  path3 <- "_study_area_buffer.shp"
  path4 <- paste0(path1, path2, path3)
  
  sf::st_write(dataset_transform_buffer, path4, overwrite = TRUE)
}

# El único que sigue muy bajito es:

species_revise_small <- buffers$species[which(buffers$Area_percent_increase <= 17)]


# manually erase the shp for these species
species_revise_small

# Indicate the species that we are going to re-run
final_cropped_models <- data.frame(final_cropped_models) %>% mutate(species_name = word(str_replace_all(final_cropped_models, "[^[:alnum:]]", " "), start = 18, end = 19, sep = fixed(" ")),
                                                                    new_buffer2 = case_when(species_name %in% species_revise_small ~ "Yes", 
                                                                                            TRUE ~ "No")) %>% filter(new_buffer1 == "Yes") %>% droplevels()

# Load raster
raster <- raster::raster(final_cropped_models$final_cropped_models[final_cropped_models$new_buffer2 == "Yes"])
  
# Extract species name 
buffers_index <- which(buffers$species == final_cropped_models$species_name[final_cropped_models$new_buffer2 == "Yes"])
  
# Reclass all values to 1
min_value <- cellStats(raster, "min")
  
# Reclassify all suitable values to 1
m <- c((min_value-0.1), 2, 1)
m <- matrix(m, ncol=3, byrow=TRUE)
raster <- raster::reclassify(raster, m)
  
# Turn distribution into shapefile
x <- rasterToPolygons(raster, na.rm = TRUE, dissolve = TRUE)
  
# Reproject tp calculate distances and areas
x <- as(x, "sf")
x <- st_transform(x, 32618)
  
# Calculate area of polygon
buffers$Distribution_area_sqkm[buffers_index] <- st_area(x)/1000
  
# Calculate perimeter of polygon
buffers$Distribution_perimeter_km[buffers_index] <-  st_perimeter(x)/1000
  
# Calculate buffer distance to achieve an increase of ~20% in polygon area
buffers$Buffer_distance_km[buffers_index] <- (((buffers$Distribution_area_sqkm[buffers_index] * 1.8) - buffers$Distribution_area_sqkm[buffers_index]) / buffers$Distribution_perimeter_km[buffers_index])/1000
  
# Apply buffer
x <- as(x, Class = "Spatial")
dataset_transform_buffer <- raster::buffer(x, width = (buffers$Buffer_distance_km[buffers_index]*1000), dissolve = TRUE)
  
# Calculate actual % increase in area (due to irregular polygons)
dataset_transform_buffer <- as(dataset_transform_buffer, "sf")
buffers$Final_area_sqkm[buffers_index] <- st_area(dataset_transform_buffer)/1000
buffers$Area_percent_increase[buffers_index] <- ((buffers$Final_area_sqkm[buffers_index] - buffers$Distribution_area_sqkm[buffers_index]) / buffers$Distribution_area_sqkm[buffers_index])*100
  
  
# Export
path1 <- "./outputs/03_connectivity/Species distributions with buffer/"
path2 <- buffers$species[buffers_index]
path3 <- "_study_area_buffer.shp"
path4 <- paste0(path1, path2, path3)
  
sf::st_write(dataset_transform_buffer, path4, overwrite = TRUE)

write.table(buffers, "clipboard", sep="\t", row.names=FALSE)

# Revision of buffer sizes: ----------------------------------------------------

# Load list of final cropped models
final_cropped_models <- list.files("./outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/FINAL_selected models/CROPPED MODELS", 
                                   pattern = ".tif$",
                                   full.names = TRUE, 
                                   recursive =  FALSE)

buffers_per_species <- list.files("./outputs/03_connectivity/Species distributions with buffer", 
                                  pattern = ".shp$",
                                  full.names = TRUE, 
                                  recursive =  FALSE)


buffers_per_species <- data.frame(buffers_per_species) %>% mutate(species_name = word(str_replace_all(buffers_per_species, "[^[:alnum:]]", " "), start = 10, end = 11, sep = fixed(" ")))

# Create table to store buffer length per species
final_buffers <- data.frame(matrix(NA,    # Create empty data frame
                             nrow = 26,
                             ncol = 6))

final_buffers <- final_buffers %>% rename(species = X1, Distribution_area_sqkm = X2, Distribution_perimeter_km = X3, 
                              Buffer_distance_km = X4, Area_percent_increase = X6, Final_area_sqkm = X5)

for (i in 13:length(final_cropped_models)){
  
  print(i)
  
  # Load raster
  raster <- raster::raster(final_cropped_models[i])
  
  # Extract species name 
  species_name <- str_replace_all(final_cropped_models[i], "[^[:alnum:]]", " ") %>% word(start=18,end=19,sep=fixed(" "))
  final_buffers$species[i] <- species_name
  # Reclass all values to 1
  min_value <- cellStats(raster, "min")
  
  # Reclassify all suitable values to 1
  m <- c((min_value-0.1), 2, 1)
  m <- matrix(m, ncol=3, byrow=TRUE)
  raster <- raster::reclassify(raster, m)
  
  # Turn distribution into shapefile
  x <- rasterToPolygons(raster, na.rm = TRUE, dissolve = TRUE)
  
  # Reproject tp calculate distances and areas
  x <- as(x, "sf")
  x <- st_transform(x, 32618)
  
  # Calculate area of polygon
  final_buffers$Distribution_area_sqkm[i] <- st_area(x)/1000
  
  # Calculate perimeter of polygon
  final_buffers$Distribution_perimeter_km[i] <-  st_perimeter(x)/1000
  
  # Calculate buffer based on produced shapefile
  buffer_shp <- st_read(buffers_per_species$buffers_per_species[buffers_per_species$species_name == final_buffers$species[i]])
  
  # Calculate actual % increase in area 
  final_buffers$Final_area_sqkm[i] <- st_area(buffer_shp)/1000
  final_buffers$Area_percent_increase[i] <- ((final_buffers$Final_area_sqkm[i] - final_buffers$Distribution_area_sqkm[i]) / final_buffers$Distribution_area_sqkm[i])*100
  
}


write.table(final_buffers, "clipboard", sep="\t", row.names=FALSE)
