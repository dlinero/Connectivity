library(tidyverse)
library(raster)
library(wallace)


remotes::install_github("https://github.com/wallaceEcoMod/wallace/tree/multiSp", force = TRUE)
# Open wallace
library(wallace)
run_wallace()


# Export one file of cleaned and thinned occurrences per species 
data <- read.csv("F:/Connectivity/outputs/02_SDMs/Intentos/final/FINAL_THINNED_DATA.csv")

# Empecemos por las que estoy segura que el MCP es el área de calibración 

# Crypturellus berlepschi - 23 records
data_CB <- data %>% filter(name == "Crypturellus berlepschi") %>% droplevels() %>%
  write.csv("F:/Connectivity/data/02_SDMs/Crypturellus_berlepschi_CleanedThinned.csv", row.names = FALSE)

# Ramphastos brevis - 128 records
data_RB <- data %>% filter(name == "Ramphastos brevis") %>% droplevels() %>%
  write.csv("F:/Connectivity/data/02_SDMs/Ramphastos_brevis_CleanedThinned.csv", row.names = FALSE)

# Grallaria hypoleuca - 71 records
data_RB <- data %>% filter(name == "Grallaria hypoleuca") %>% droplevels() %>%
  write.csv("F:/Connectivity/data/02_SDMs/Grallaria_hypoleuca_CleanedThinned.csv", row.names = FALSE)

# Grallaria flavotincta - 31 records & Picumnus lafresnayi - 215 records
data_GF_PL <- data %>% filter(name == "Grallaria flavotincta" | name == "Picumnus lafresnayi") %>% droplevels() %>%
  rename(scientific_name = name) %>% write.csv("F:/Connectivity/data/02_SDMs/GrallariaFlavotincta_PicumnusLafresnayi_CleanedThinned.csv", row.names = FALSE)

# Tangara johannae - 36 records &  Patagioenas goodsoni - 71 records
data_TJ_PG <- data %>% filter(name == "Tangara johannae" | name == "Patagioenas goodsoni") %>% droplevels() %>%
  rename(scientific_name = name)

write.csv(data_TJ_PG, "F:/Connectivity/data/02_SDMs/Tangarajohannae_PatagioenasGoodsoni_CleanedThinned.csv", row.names = FALSE)

# Picummus lafresnayi
data %>% filter(name == "Picumnus lafresnayi") %>% droplevels() %>%
  rename(scientific_name = name) %>% write.csv("F:/Connectivity/data/02_SDMs/PicumnusLafresnayi_CleanedThinned.csv", row.names = FALSE)

# Crypturellus soui - 2929 records
data_CS <- data %>% filter(name == "Crypturellus soui") %>% droplevels() %>%
  write.csv("F:/Connectivity/data/02_SDMs/Crypturellus_soui_CleanedThinned.csv", row.names = FALSE)

# Spizaetus ornatus - 471 records
data_SO <- data %>% filter(name == "Spizaetus ornatus") %>% droplevels() 
write.csv(data_SO, "F:/Connectivity/data/02_SDMs/Spizaetus_ornatus_CleanedThinned.csv", row.names = FALSE)

# Hypopyrrhus pyrohypogaster with 121 records
data_HP <- data %>% filter(name == "Hypopyrrhus pyrohypogaster") %>% droplevels() %>%
  rename(scientific_name = name) %>% write.csv("F:/Connectivity/data/02_SDMs/Hypopyrrhus_pyrohypogaster_CleanedThinned.csv", row.names = FALSE)

# Picumnus cinnamomeus - 45 records & Ognorhynchus icterotis - 39 records & Crypturellus erythropus - 48 records
data %>% filter(name == "Picumnus cinnamomeus" | name == "Ognorhynchus icterotis" | name == "Crypturellus erythropus") %>% droplevels() %>%
  rename(scientific_name = name) %>% write.csv("F:/Connectivity/data/02_SDMs/PicumnuCinnamomeus_OgnorhynchusIcterotis_CrypturellusErythropus_CleanedThinned.csv", row.names = FALSE)

# Myioborus flavivertex - 16 records & Grallaria bangsi - 18 records & Mitu salvini  - 22 records
data %>% filter(name == "Myioborus flavivertex" | name == "Grallaria bangsi" | name == "Mitu salvini") %>% droplevels() %>%
  rename(scientific_name = name) %>% write.csv("F:/Connectivity/data/02_SDMs/MyioborusFlavivertex_GrallariaBangsi_MituSalvini_CleanedThinned.csv", row.names = FALSE)

# Phaethornis malaris - 226 records

# Picumnus squamulatus - 332 records

# Odontophorus hyperythrus - 132 records

# Crypturellus soui - 2929 records

# Grallaria nuchalis - 94 records & Mitu tomentosum - 99 records
data %>% filter(name == "Grallaria nuchalis" | name == "Mitu tomentosum") %>% droplevels() %>%
  rename(scientific_name = name) %>% write.csv("F:/Connectivity/data/02_SDMs/GrallariaNuchalis_MituTomentosum_CleanedThinned.csv", row.names = FALSE)



# Spizaetus isidori

# 

# Campephilus pollens - 172 records

# Pyrrhura melanura - 229 records

# Pseudastur albicollis - 871 records

# Micrastur semitorquatus - 1989 records

# Buteo albigula - 179 records

# Ara militaris - 135 records

# Ara macao - 1020 records





# Project bioclim rasters 


r <- raster("F:/Connectivity/outputs/02_SDMs/wc/bio19_LAC.tif")
# an example for the common CRS called WGS84 (a latitude / longitude projection)
# the text format for these CRS's is called proj4
crs(r) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# save the raster with defined CRS to file
writeRaster(r, "F:/Connectivity/outputs/02_SDMs/wc/30s_bio19_LAC.tif")

gc()
memory.limit(99999999)


tempdir()
dir.create(tempdir())


# ------------------------------------------------------------------------------------

occs_Ce <- data %>% filter(name == "Grallaria nuchalis") %>% droplevels() %>%
  rename(scientific_name = name) 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Ce[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Ce <- occs_Ce[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Ce, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Intentos/final/Grallaria nuchalis_mcp_10km_thin5km.shp")


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Sample background points from the provided area
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 10000)

# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(scientific_name = paste0("bg_", "Crypturellus berlepschi"), bgSample_Cb,
                       occID = NA, year = NA, institution_code = NA, country = NA,
                       state_province = NA, locality = NA, elevation = NA,
                       record_type = NA, bgEnvsVals_Cb)


# R code to get partitioned data
groups_Ce <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Ce, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = FALSE,
  numCores = 7)


# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models/Grallaria_nuchalis/eval_table.csv", row.names = FALSE)


# 


# Select current model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)
#Get values of prediction
mapPredVals_Cb <- getRasterVals(predSel_Cb, "cloglog")
#Define colors and legend  
rasCols <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
legendPal <- colorNumeric(rev(rasCols), mapPredVals_Cb, na.color = 'transparent')
rasPal <- colorNumeric(rasCols, mapPredVals_Cb, na.color = 'transparent')
#Generate map
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) 
m  %>%
  addLegend("bottomright", pal = legendPal,
            title = "Predicted Suitability<br>(Training)",
            values = mapPredVals_Cb, layerId = "train",
            labFormat = reverseLabels(2, reverse_order = TRUE)) %>% 
  #add occurrence data
  addCircleMarkers(data = occs_Cb, lat = ~latitude, lng = ~longitude,
                   radius = 5, color = 'red', fill = TRUE, fillColor = "red",
                   fillOpacity = 0.2, weight = 2, popup = ~pop) %>% 
  ##Add model prediction
  addRasterImage(predSel_Cb, colors = rasPal, opacity = 0.7,
                 group = 'vis', layerId = 'mapPred', method = "ngb") %>%
  ##add background polygons
  addPolygons(data = bgExt_Cb,fill = FALSE,
              weight = 4, color = "blue", group = 'proj')


writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models/Grallaria_nuchalis/cloglog.tif")

# First must generate the projection area based on user provided files
##User must input the path to shapefile or csv file and the file name 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Crypturellus_berlepschi_Ayerbe_IUCN_buffer.shp")

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.0.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 
#store the cropped projection variables
projExt_Cb <- proj_area_Cb$projExt

###Make map of projection
bb_Cb <-  bgExt_Cb@bbox
bbZoom <- polyZoom(bb_Cb[1, 1], bb_Cb[2, 1], bb_Cb[1, 2], 
                   bb_Cb[2, 2], fraction = 0.05)
mapProjVals_Cb <- getRasterVals(proj_area_Cb$projArea,"cloglog")
rasCols_Cb <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Cb), mapProjVals_Cb, na.color = 'transparent')
rasPal_Cb <- colorNumeric(rasCols_Cb, mapProjVals_Cb, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) 
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  addLegend("bottomright", pal = legendPal,
            title = "Predicted Suitability<br>(Projected)",
            values = mapProjVals_Cb, layerId = 'proj',
            labFormat = reverseLabels(2, reverse_order = TRUE)) %>%
  # map model prediction raster and projection polygon
  clearMarkers() %>% clearShapes() %>% removeImage('projRas') %>%
  addRasterImage(proj_area_Cb$projArea, colors = rasPal_Cb, opacity = 0.7,
                 layerId = 'projRas', group = 'proj', method = "ngb") %>%
  ##add projection polygon (user provided area)
  addPolygons(data = proj_userExt_Cb, fill = FALSE,
              weight = 4, color = "blue", group = 'proj')


writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/thinned_projection_jacknife.tif")

