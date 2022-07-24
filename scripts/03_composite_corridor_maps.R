library(raster)
library(stringr)
library(dplyr)
library(arcgisbinding)


Equidistant_cylindrical <- crs("+proj=eqc +lon_0=-72.2460938 +lat_ts=0 +datum=WGS84 +units=m +no_defs")

# Progress bar
rasterOptions(progress = 'text',timer=TRUE)

# Check ArcGIS licence
arc.check_product()

# Upload corridor rasters 

# Load list of species
species <- read.csv("./data/Selected_species_FINAL.csv")
species <- species %>% 
  mutate(species_folder = str_replace_all(Species, " ", "_")) %>%
  dplyr::select(species_folder) %>%
  mutate(decil1 = NA, decil2 = NA, decil3 = NA)

# Establish path that holds the corridor raster
corridor_paths <- c("./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/")

# Calculate raster deciles ---------------------------------------------------------

# Load each one of the rasters and calculate decils 

for (i in 1:nrow(species)){
  print(i)
  
  species_path <- paste0(corridor_paths, species$species_folder[i], "/output/corridors.gdb/", species$species_folder[i], "_corridors")
  data <- as.raster(arc.raster(arc.open(species_path)))
  
  # Calculate deciles
  deciles <- raster::quantile(data, probs = seq(.1, .3, by = .1))
  
  # Save in table
  species$decil1[i] <- deciles[1]
  species$decil2[i] <- deciles[2]
  species$decil3[i] <- deciles[3]
  
  # Remove temporary variables
  rm(species_path, data, deciles)
  
}

# Export table 
write.csv(species, "./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/corridors_deciles_perSpecies.csv", row.names = FALSE)

# Reclassify rasters based on deciles ----------------------------------------------------------------

for (i in 1:nrow(species)){
  print(i)
  
  species_path <- paste0(corridor_paths, species$species_folder[i], "/output/corridors.gdb/", species$species_folder[i], "_corridors")
  data <- as.raster(arc.raster(arc.open(species_path)))
  
  # Calculate deciles
  deciles <- raster::quantile(data, probs = seq(.1, .9, by = .1))
  
  reclass_df <- c(-0.1, deciles[1], 10,
                  deciles[1], deciles[2], 9,
                  deciles[2], deciles[3], 8,
                  deciles[3], deciles[4], 7, 
                  deciles[4], deciles[5], 6,
                  deciles[5], deciles[6], 5,
                  deciles[6], deciles[7], 4, 
                  deciles[7], deciles[8], 3,
                  deciles[8], deciles[9], 2, 
                  deciles[9], Inf, 1)
  
  reclass_df <- matrix(reclass_df,
                      ncol = 3,
                      byrow = TRUE)
  
  data_classified <- reclassify(data,
                               reclass_df)
  
  # Set coordinate system
  crs(data_classified) <- Equidistant_cylindrical
  # Export
  writeRaster(data_classified, paste0(corridor_paths, species$species_folder[i], "/output/", species$species_folder[i], "_reclassified_1to10.tif"), overwrite = TRUE)
  
  # Remove temporary variables
  rm(species_path, data, deciles, reclass_df, data_classified)
  
}


# Mosaic results for each group ----------------------------------------------------------------

# This step is done in ArcGIS running the Mosaic function, this is the log 

#Mosaic To New Raster - Low
#=====================
  #Parameters

#Input Rasters     'Crypturellus berlepschi\Crypturellus_berlepschi_reclassified_1to10.tif';'Crypturellus erythropus\Crypturellus_erythropus_reclassified_1to10.tif';'Crypturellus soui\Crypturellus_soui_reclassified_1to10.tif';'Grallaria bangsi\Grallaria_bangsi_reclassified_1to10.tif';'Grallaria flavotincta\Grallaria_flavotincta_reclassified_1to10.tif';'Grallaria hypoleuca\Grallaria_hypoleuca_reclassified_1to10.tif';'Grallaria nuchalis\Grallaria_nuchalis_reclassified_1to10.tif';'Henicorhina leucophrys\Henicorhina_leucophrys_reclassified_1to10.tif';'Myioborus flavivertex\Myioborus_flavivertex_reclassified_1to10.tif';'Odontophorus hyperythrus\Odontophorus_hyperythrus_reclassified_1to10.tif';'Phaethornis malaris\Phaethornis_malaris_reclassified_1to10.tif';'Picumnus cinnamomeus\Picumnus_cinnamomeus_reclassified_1to10.tif';'Picumnus squamulatus\Picumnus_squamulatus_reclassified_1to10.tif';'Tangara johannae\Tangara_johannae_reclassified_1to10.tif'
#Output Location     F:\Connectivity\outputs\03_connectivity\Least Cost Corridors\Least Cost Corridors_corine\Results\Composite_maps.gdb
#Raster Dataset Name with Extension     TRY
#Spatial Reference for  Raster     PROJCS["Equirectangular",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Plate_Carree"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-72.2460938],UNIT["Meter",1.0]]
#Pixel Type     8_BIT_UNSIGNED
#Cellsize     
#Number of Bands     1
#Mosaic Operator     SUM
#Mosaic Colormap Mode     FIRST
#Output Raster Dataset     F:\Connectivity\outputs\03_connectivity\Least Cost Corridors\Least Cost Corridors_corine\Results\Composite_maps.gdb\TRY
#=====================
  #Messages

#Start Time: Friday, July 22, 2022 4:00:05 PM
#Crypturellus berlepschi\Crypturellus_berlepschi_reclassified_1to10.tif is loading...
#Crypturellus erythropus\Crypturellus_erythropus_reclassified_1to10.tif is loading...
#Crypturellus soui\Crypturellus_soui_reclassified_1to10.tif is loading...
#Grallaria bangsi\Grallaria_bangsi_reclassified_1to10.tif is loading...
#Grallaria flavotincta\Grallaria_flavotincta_reclassified_1to10.tif is loading...
#Grallaria hypoleuca\Grallaria_hypoleuca_reclassified_1to10.tif is loading...
#Grallaria nuchalis\Grallaria_nuchalis_reclassified_1to10.tif is loading...
#Henicorhina leucophrys\Henicorhina_leucophrys_reclassified_1to10.tif is loading...
#Myioborus flavivertex\Myioborus_flavivertex_reclassified_1to10.tif is loading...
#Odontophorus hyperythrus\Odontophorus_hyperythrus_reclassified_1to10.tif is loading...
#Phaethornis malaris\Phaethornis_malaris_reclassified_1to10.tif is loading...
#Picumnus cinnamomeus\Picumnus_cinnamomeus_reclassified_1to10.tif is loading...
#Picumnus squamulatus\Picumnus_squamulatus_reclassified_1to10.tif is loading...
#Tangara johannae\Tangara_johannae_reclassified_1to10.tif is loading...
#Succeeded at Friday, July 22, 2022 4:00:44 PM (Elapsed Time: 38.94 seconds)


#Mosaic Medium 
#=====================
  #Parameters

#Input Rasters     'Campephilus pollens\Campephilus_pollens_reclassified_1to10.tif';'Mitu tomentosum\Mitu_tomentosum_reclassified_1to10.tif';'Pseudastur albicollis\Pseudastur_albicollis_reclassified_1to10.tif';'Pyrrhura melanura\Pyrrhura_melanura_reclassified_1to10.tif';'Ramphastos brevis\Ramphastos_brevis_reclassified_1to10.tif';'Spizaetus isidori\Spizaetus_isidori_reclassified_1to10.tif'
#Target Raster     F:\Connectivity\outputs\03_connectivity\Least Cost Corridors\Least Cost Corridors_corine\Results\Composite_maps.gdb\Mosaic_medium_group
#Mosaic Operator     SUM
#Mosaic Colormap Mode     FIRST
#Ignore Background Value     
#NoData Value     -3.4E+38
#Convert 1 bit data to 8 bit     NONE
#Mosaicking Tolerance     0
#Updated Target Raster     F:\Connectivity\outputs\03_connectivity\Least Cost Corridors\Least Cost Corridors_corine\Results\Composite_maps.gdb\Mosaic_medium_group
#Color Matching Method     NONE
#=====================
  #Messages

#Start Time: Thursday, July 21, 2022 4:46:45 PM
#Campephilus pollens\Campephilus_pollens_reclassified_1to10.tif is loading...
#Mitu tomentosum\Mitu_tomentosum_reclassified_1to10.tif is loading...
#Pseudastur albicollis\Pseudastur_albicollis_reclassified_1to10.tif is loading...
#Pyrrhura melanura\Pyrrhura_melanura_reclassified_1to10.tif is loading...
#Ramphastos brevis\Ramphastos_brevis_reclassified_1to10.tif is loading...
#Spizaetus isidori\Spizaetus_isidori_reclassified_1to10.tif is loading...
#Succeeded at Thursday, July 21, 2022 4:47:04 PM (Elapsed Time: 18.28 seconds)

#Mosaic large
#=====================
  #Parameters

#Input Rasters     'Ara macao\Ara_macao_reclassified_1to10.tif';'Buteo albigula\Buteo_albigula_reclassified_1to10.tif';'Ognorhynchus icterotis\Ognorhynchus_icterotis_reclassified_1to10.tif';'Patagioenas goodsoni\Patagioenas_goodsoni_reclassified_1to10.tif'
#Target Raster     F:\Connectivity\outputs\03_connectivity\Least Cost Corridors\Least Cost Corridors_corine\Results\Composite_maps.gdb\Mosaic_large_group
#Mosaic Operator     SUM
#Mosaic Colormap Mode     FIRST
#Ignore Background Value     
#NoData Value     -3.4E+38
#Convert 1 bit data to 8 bit     NONE
#Mosaicking Tolerance     0
#Updated Target Raster     F:\Connectivity\outputs\03_connectivity\Least Cost Corridors\Least Cost Corridors_corine\Results\Composite_maps.gdb\Mosaic_large_group
#Color Matching Method     NONE
#=====================
  #Messages

#Start Time: Friday, July 22, 2022 11:32:29 AM
#Ara macao\Ara_macao_reclassified_1to10.tif is loading...
#Buteo albigula\Buteo_albigula_reclassified_1to10.tif is loading...
#Ognorhynchus icterotis\Ognorhynchus_icterotis_reclassified_1to10.tif is loading...
#Patagioenas goodsoni\Patagioenas_goodsoni_reclassified_1to10.tif is loading...
#Succeeded at Friday, July 22, 2022 11:32:40 AM (Elapsed Time: 10.98 seconds)


# Cut off rasters -----------------------------------------------------------------
# Mosaic deciles

low <- raster("./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/Mosaic_low_group.tif")
medium <- raster("./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/Mosaic_medium_group.tif")
large <- raster("./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/Mosaic_large_group.tif")
all <- raster("./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/Mosaic_all_species.tif")


# Remove zero values
low[low==0]<-NA
medium[medium==0]<-NA
large[large==0]<-NA
all[all==0]<-NA

# Calculate deciles
raster::quantile(low, probs = seq(.1, .9, by = .1))
raster::quantile(medium, probs = seq(.1, .9, by = .1))
raster::quantile(large, probs = seq(.1, .9, by = .1))
raster::quantile(all, probs = seq(.1, .9, by = .1))


# Intentar con el top 20% 
top_20_low <- low
top_20_low[top_20_low<22]<-NA

top_20_medium <- medium
top_20_medium[top_20_medium<23]<-NA

top_20_large <- large
top_20_large[top_20_large<15]<-NA

top_20_all <- all
top_20_all[top_20_all<58]<-NA


# Intentar con valores mayores a 20 (es decir corredores para al menos dos especies)
atleast_2_medium <- medium
atleast_2_medium[atleast_2_medium<20]<-NA

# Export
writeRaster(top_20_low, "./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/Mosaic_low_group_top20.tif")
writeRaster(top_20_medium, "./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/Mosaic_medium_group_top20.tif")
writeRaster(atleast_2_medium, "./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/Mosaic_medium_group_2species.tif")
writeRaster(top_20_large, "./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/Mosaic_large_group_top20.tif")
writeRaster(top_20_all, "./outputs/03_connectivity/Least Cost Corridors/Least Cost Corridors_corine/Results/Mosaic_all_species_top20.tif")
