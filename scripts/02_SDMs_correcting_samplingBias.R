library(tidyverse)
library(wallace)
library(dismo)
library(rgeos)
library(ENMeval)
library(janitor) # For cleanning dataset column names
library(adehabitatHR) # For minimum convex polygons 
library(ggpubr)
library(rstatix)
library(nortest)


# Con los datos thinned + background random y TGB con evaluaciones AUC y AIC 
# (4 modelos). Y también datos not thinned + TGB con evaluaciones AUC y AIC (2 modelos )
# No tiene sentido hacer not thinned + random background por que no se controlaría
# el sampling bias. 

data_thinned <- read.csv("F:/Connectivity/outputs/02_SDMs/MCPs/FINAL_THINNED_DATA.csv")


# Myioborus flavivertex - 16 records --------------------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Myioborus flavivertex") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Myioborus flavivertex_mcp_10km_thin1km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 



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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Myioborus_flavivertex_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.L_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/projection_thinned_randomBackground_AIC.tif")


## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Myioborus flavivertex") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Myioborus flavivertex_mcp_10km_thin1km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 



# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Myioborus flavivertex") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/cloglog_MCP_thinned_TGB_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Myioborus_flavivertex_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.1",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.L_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/projection_thinned_TGB_AIC.tif")

## Not thinned + TGB ------------------------------------------------------------------


data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Myioborus flavivertex") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)



## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Myioborus flavivertex_mcp_10km_thin1km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 



# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Myioborus flavivertex") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 380 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/cloglog_MCP_notThinned_TGB_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Myioborus_flavivertex_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.L_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Myioborus flavivertex/projection_notThinned_TGB_AIC.tif")


# Grallaria flavotincta - 31 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Grallaria flavotincta") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria flavotincta_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Grallaria_flavotincta_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQH_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), "envs_Cb"))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Grallaria flavotincta") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria flavotincta_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Grallaria flavotincta") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/cloglog_MCP_thinned_TGB_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Grallaria_flavotincta_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.LQH_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), "envs_Cb"))

## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Grallaria flavotincta") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria flavotincta_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 



# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Grallaria flavotincta") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 531 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/cloglog_MCP_notThinned_TGB_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Grallaria_flavotincta_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.LQ_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.1",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria flavotincta/projection_notThinned_TGB_AIC.tif")

# Grallaria hypoleuca - 71 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Grallaria hypoleuca") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria hypoleuca_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Grallaria_hypoleuca_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.H_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), "envs_Cb"))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Grallaria hypoleuca") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria hypoleuca_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Grallaria hypoleuca") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/cloglog_MCP_thinned_TGB_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Grallaria_hypoleuca_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), "envs_Cb"))

## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Grallaria hypoleuca") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria hypoleuca_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Grallaria hypoleuca") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 557 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/cloglog_MCP_notThinned_TGB_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Grallaria_hypoleuca_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.1",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.0.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria hypoleuca/projection_notThinned_TGB_AIC.tif")



# Grallaria bangsi- 18 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Grallaria bangsi") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria bangsi_mcp_10km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/IUCN/Grallaria_bangsi_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.L_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

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

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), "envs_Cb"))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Grallaria bangsi") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria bangsi_mcp_10km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Grallaria bangsi") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/cloglog_MCP_thinned_TGB_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/IUCN/Grallaria_bangsi_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.LQ_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "proj_userExt_Cb")))

## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Grallaria bangsi") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria bangsi_mcp_10km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Grallaria bangsi") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)


# Run maxent model for the selected species ## Removed 313 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/cloglog_MCP_notThinned_TGB_AUC.tif")

# Generate projection

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.0.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.LQH_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria bangsi/projection_notThinned_TGB_AIC.tif")


# Picumnus cinnamomeus - 45 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Picumnus cinnamomeus") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Picumnus cinnamomeus_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Picumnus_cinnamomeus_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQ_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Picumnus cinnamomeus") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Picumnus cinnamomeus") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.LQ_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Picumnus cinnamomeus") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 



# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Picumnus cinnamomeus") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 531 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.LQH_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Picumnus cinnamomeus/projection_notThinned_TGB_AIC.tif")



# Odontophorus hyperythrus - 132 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Odontophorus hyperythrus") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Odontophorus hyperythrus_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 

# Run maxent model for the selected species
gc()
memory.limit(9999999999)
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Odontophorus_hyperythrus_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.0.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Odontophorus hyperythrus") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Odontophorus hyperythrus") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

gc()
# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.2.5 (also for AIC)

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/cloglog_MCP_thinned_TGB_AUC&AIC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.2.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/projection_thinned_TGB_AUC&AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))

## Not thinned + TGB  ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Odontophorus hyperythrus") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Odontophorus hyperythrus") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 778 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.LQH_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Odontophorus hyperythrus/projection_notThinned_TGB_AIC.tif")


# Grallaria nuchalis - 94 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Grallaria nuchalis") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria nuchalis_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Grallaria_nuchalis_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQH_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Grallaria nuchalis") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Grallaria nuchalis") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Grallaria nuchalis") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 



# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Grallaria nuchalis") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 564 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Grallaria nuchalis/projection_notThinned_TGB_AIC.tif")

# Crypturellus erythropus - 48 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Crypturellus erythropus") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Crypturellus erythropus_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Crypturellus_erythropus_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.0.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.L_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

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

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Crypturellus erythropus") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Crypturellus erythropus") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.LQH_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Crypturellus erythropus") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Crypturellus erythropus") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 564 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.LQ_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus erythropus/projection_notThinned_TGB_AIC.tif")

# Tangara johannae - 36 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Tangara johannae") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Tangara johannae_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Tangara_johannae_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.2.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = was the same selected for AUC 

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Tangara johannae") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Tangara johannae") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Tangara johannae") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Tangara johannae") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 74 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Tangara johannae/projection_notThinned_TGB_AIC.tif")


# Ramphastos brevis - 128 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Ramphastos brevis") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Ramphastos brevis_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 

gc()
memory.limit(99999999)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Ramphastos_brevis_Ayerbe_IUCN_buffer .shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQH_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Ramphastos brevis") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Ramphastos brevis") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Ramphastos brevis") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Ramphastos brevis") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 564 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ramphastos brevis/projection_notThinned_TGB_AIC.tif")



# Mitu salvini  - 22 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Mitu salvini") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Mitu salvini_mcp_10km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

gc()
memory.limit(99999999)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Mitu_salvini_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.H_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Mitu salvini") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Mitu salvini") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu salvini/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------
# The number of records once ENMeval removes occurrence localities that shared
# the same grid cell (25), is the same as the thinned number 



# Mitu tomentosum - 99 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Mitu tomentosum") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Mitu tomentosum_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Mitu_tomentosum_merge_10km.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQH_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Mitu tomentosum") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Mitu tomentosum") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/cloglog_MCP_thinned_TGB_AUC.tif", overwrite = TRUE)

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.L_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Mitu tomentosum") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Mitu tomentosum") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 130 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.L_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Mitu tomentosum/projection_notThinned_TGB_AIC.tif")


# Campephilus pollens  - 172 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Campephilus pollens") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Campephilus pollens_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 

gc()
memory.limit(99999999)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Campephilus_pollens_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQ_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Campephilus pollens") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Campephilus pollens") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Campephilus pollens") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Campephilus pollens") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 74 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Campephilus pollens/projection_notThinned_TGB_AIC.tif")

# Pyrrhura melanura  - 229 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Pyrrhura melanura") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Pyrrhura melanura_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 

gc()
memory.limit(99999999)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Pyrrhura_melanura_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQH_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Pyrrhura melanura") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Pyrrhura melanura") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Pyrrhura melanura") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Pyrrhura melanura") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 74 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Pyrrhura melanura/projection_notThinned_TGB_AIC.tif")


# Hypopyrrhus pyrohypogaster  - 121 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Hypopyrrhus pyrohypogaster") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Hypopyrrhus pyrohypogaster_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 

gc()
memory.limit(99999999)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Hypopyrrhus_pyrohypogaster_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.0.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Hypopyrrhus pyrohypogaster") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Hypopyrrhus pyrohypogaster") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.2.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Hypopyrrhus pyrohypogaster") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Hypopyrrhus pyrohypogaster") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 74 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.1",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Hypopyrrhus pyrohypogaster/projection_notThinned_TGB_AIC.tif")



# Patagioenas goodsoni  - 71 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Patagioenas goodsoni") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Patagioenas goodsoni_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 

gc()
memory.limit(99999999)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Patagioenas_goodsoni_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

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

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.0.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Patagioenas goodsoni") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Patagioenas goodsoni") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4 

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Patagioenas goodsoni") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Patagioenas goodsoni") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 299 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Patagioenas goodsoni/projection_notThinned_TGB_AIC.tif")


# Ognorhynchus icterotis  - 39 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Ognorhynchus icterotis") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Ognorhynchus icterotis_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 



# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Ognorhynchus_icterotis_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.2.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = Mismo que AUC

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Ognorhynchus icterotis") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Ognorhynchus icterotis") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Ognorhynchus icterotis") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Ognorhynchus icterotis") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 299 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.1.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Ognorhynchus icterotis/projection_notThinned_TGB_AIC.tif")

# Crypturellus berlepschi - 23 records ---------------------------------------------

## Thinned + random background -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Crypturellus berlepschi") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Crypturellus berlepschi_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 


# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/evalTable_thinned_randomBackground.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_thinned_randomBackground_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Crypturellus_berlepschi_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_thinned_randomBackground_AUC.tif")

# deltaAICc = fc.LQH_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_thinned_randomBackground_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_thinned_randomBackground_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Crypturellus berlepschi") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Crypturellus berlepschi") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/evalTable_thinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_thinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_thinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.4 

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_thinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_thinned_TGB_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Crypturellus berlepschi") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Crypturellus berlepschi") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/TGB_points.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 20 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/evalTable_notThinned_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_notThinned_TGB_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_notThinned_TGB_AUC.tif")

# deltaAICc = fc.H_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_notThinned_TGB_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_notThinned_TGB_AIC.tif")

## Thinned + random background + new MCP -----------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Crypturellus berlepschi") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/modified/Crypturellus berlepschi_mcp_modified.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 


# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/evalTable_thinned_randomBackground_newMCP.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_thinned_randomBackground_newMCP_AUC.tif")

# Generate projection
# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Crypturellus_berlepschi_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.0.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_thinned_randomBackground_newMCP_AUC.tif")

# deltaAICc = fc.L_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_thinned_randomBackground_newMCP_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.1",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_thinned_randomBackground_newMCP_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb", "data_thinned")))

## Thinned + TGB + new MCP -------------------------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Crypturellus berlepschi") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Crypturellus berlepschi") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/evalTable_thinned_TGB_newMCP.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.3.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_thinned_TGB_newMCP_AUC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.3.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_thinned_TGB_newMCP_AUC.tif")

# deltaAICc = fc.H_rm.4 

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_thinned_TGB_newMCP_AIC.tif")

# Generate projection

# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_thinned_TGB_newMCP_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "bgExt_Cb", "proj_userExt_Cb")))


## Not thinned + TGB ------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Crypturellus berlepschi") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Crypturellus berlepschi") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Export background points 
raster::shapefile(cleanned, filename = "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/TGB_points_newMCP.shp")


# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2) 


# Run maxent model for the selected species ## Removed 20 occurrence localities that shared the same grid cell
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/evalTable_notThinned_TGB_newMCP.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.4 (same for AIC)

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/cloglog_MCP_notThinned_TGB_newMCP_AUC&AIC.tif")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Crypturellus berlepschi/projection_notThinned_TGB_newMCP_AUC&AIC.tif")


# Crypturellus berlepschi - thinned -----------------------------------------------------------------------------

# NOTE: provide the folder path of the .csv file
occs_path <- "F:/Connectivity/data/02_SDMs"
occs_path <- file.path(occs_path, "Crypturellus_berlepschi_CleanedThinned.csv")
# get a list of species occurrence data
userOccs_Cb <- occs_userOccs(
  txtPath = occs_path, 
  txtName = "Crypturellus_berlepschi_CleanedThinned.csv", 
  txtSep = ",", 
  txtDec = ".")
occs_Cb <- userOccs_Cb$Crypturellus_berlepschi$cleaned


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Intentos/final/Crypturellus berlepschi_mcp_10km_thin5km.shp")

# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)


# Sample background points from the provided area

# Load cleaned, not thinned, data
cleanned <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

# Remove focal species points, and records that fall outside the calibration area
cleanned<- cleanned %>% 
  filter(SCIENTIFIC.NAME != "Crypturellus berlepschi") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(scientific_name = paste0("bg_", "Crypturellus berlepschi"), bgSample_Cb,
                       occID = NA, year = NA, institution_code = NA, country = NA,
                       state_province = NA, locality = NA, elevation = NA,
                       record_type = NA, bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = 23) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = FALSE,
  numCores = 7)

# https://githubmemory.com/repo/jamiemkass/ENMeval/issues/105
# Warnings: In cor(x, y) : the standard deviation is zero

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/thinned_TGB/eval_table_jacknife.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC

opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq # No funciona por que tiene filas con NA
# AUC = 0.59, ningún AUC pasa de 0.59

# Select current model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.4"]]
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


writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/thinned_cloglog_jacknife.tif")

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

# Try with a different partition
# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)


# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = FALSE,
  numCores = 7)

# https://githubmemory.com/repo/jamiemkass/ENMeval/issues/105
# Warnings: In cor(x, y) : the standard deviation is zero

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/thinned_TGB/eval_table_check1.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC

opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq # No funciona por que tiene filas con NA
# AUC = 0.6, ningún AUC pasa de 0.6

# Select current model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.1"]]
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


writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/thinned_cloglog_jacknife.tif")

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




# Crypturellus berlepschi - not thinned -----------------------------------------------------------------------------

data <- read.csv("F:/Connectivity/data/02_SDMs/cleanedData_phase2.csv", header = TRUE)

occs_Cb <- data %>% filter(SCIENTIFIC.NAME == "Crypturellus berlepschi") %>%
  
  clean_names() %>% droplevels() %>% dplyr::select(scientific_name, longitude, latitude)

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Create Minimum Convex Polygon with not thinned data (with a buffer of 0.1°)
# Turn to shapefile 
x_mcp <- occs_Cb[c('longitude', 'latitude')]
# Assign the columns with the coordinates
sp::coordinates(x_mcp) <- ~ longitude + latitude
# Calculate minimum convex polygon
bgExt_Cb <- mcp(x_mcp, percent = 100)
# Apply buffer 
bgExt_Cb <- rgeos::gBuffer(bgExt_Cb, width = 0.1)
# Assign coordinates
crs(bgExt_Cb) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)


# Sample background points from the provided area

# Remove focal species points, and records that fall outside the calibration area
cleanned<- data %>% 
  filter(SCIENTIFIC.NAME != "Crypturellus berlepschi") %>% droplevels() %>%
  dplyr::select(LONGITUDE, LATITUDE)

# Turn to spatial points dataframe 
sp::coordinates(cleanned) <- ~LONGITUDE+LATITUDE 
sp::proj4string(cleanned) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
cleanned <- gIntersection(bgExt_Cb, cleanned) 

# Turn back to dataframe 
bgSample_Cb <- as.data.frame(cleanned) %>%
  rename(longitude = x, latitude = y)


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)


##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(scientific_name = paste0("bg_", "Crypturellus berlepschi"), bgSample_Cb,
                       occID = NA, year = NA, institution_code = NA, country = NA,
                       state_province = NA, locality = NA, elevation = NA,
                       record_type = NA, bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = 52) 


# Run maxent model for the selected species
model_Cb2 <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)  # Less warnings, only 8 compared to 32 from before

# https://githubmemory.com/repo/jamiemkass/ENMeval/issues/105
# Warnings: In cor(x, y) : the standard deviation is zero

# Overall results
res2 <- eval.results(model_Cb2)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC

opt.seq <- res2 %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq # No funciona por que tiene filas con NA


# El siguiente omission rate más bajo es 0.1384615
opt.seq2 <- res2 %>% 
  filter(or.10p.avg >= 0.13 & or.10p.avg < 0.14) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq2 # AUC = 0.6, ningún AUC pasa de 0.60


# Select current model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.0.5"]]
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


writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/notthinned_cloglog.tif")

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


writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/notthinned_projection.tif")

# Crypturellus berlepschi - new MCP -----------------------------------------------------------------------------

# NOTE: provide the folder path of the .csv file
occs_path <- "F:/Connectivity/data/02_SDMs"
occs_path <- file.path(occs_path, "Crypturellus_berlepschi_CleanedThinned.csv")
# get a list of species occurrence data
userOccs_Cb <- occs_userOccs(
  txtPath = occs_path, 
  txtName = "Crypturellus_berlepschi_CleanedThinned.csv", 
  txtSep = ",", 
  txtDec = ".")
occs_Cb <- userOccs_Cb$Crypturellus_berlepschi$cleaned


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/new_MCP.shp")

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
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)


# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = FALSE,
  numCores = 7)

# Overall results
res <- eval.results(model_Cb)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC

opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq # No funciona por que tiene filas con NA


# El siguiente omission rate más bajo es 0.1384615
opt.seq2 <- res %>% 
  filter(or.10p.avg >= 0.13 & or.10p.avg < 0.14) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq2 # AUC = 0.6, ningún AUC pasa de 0.60

# Select current model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.0.5"]]
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


writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/thinned_cloglog.tif")

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


writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/thinned_projection.tif")



# Sampling bias map -----------------------------------------------------------------------------

# Aqui la idea es seleccionar los punto de background dependiendo de la probabilidad 
# que tengan de ser muestreados

# Inicialmente, voy a seguir un consejo de researchgate para el # de puntos de background:
# 20% of the cells which are neither NA nor holding a presence is a good starting point.

# Load raster
sampling_prob <- raster::raster("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/sampling_prob.tif")

# Count total number of cells 
ncell(sampling_prob)

# Check for NAs
cellStats(sampling_prob, "countNA")

# Calculate number of cells with species' presence points
npresence <- cellFromXY(sampling_prob,occs_xy_Cb)

# Calculate the number of pixels to be sampled
n_background <- (20 * (ncell(sampling_prob) - (cellStats(sampling_prob, "countNA") + length(unique(npresence)))))/100

# Replace NAs with zeros
sampling_prob[is.na(sampling_prob[])] <- 0 

# Should I not allow background to include the same points as presences?

# Sample points depending on survey probability
bg <- xyFromCell(sampling_prob, sample(ncell(sampling_prob), n_background, prob = values(sampling_prob)))

bg2 <- data.frame(bg)

# Turn into points 
phae_data <- sf::st_as_sf(bg2, coords = c("x", "y"), crs = 4326)

# Export
sf::st_write(phae_data, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/try.shp")


# Campephilus pollens  - 172 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Campephilus pollens") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Campephilus pollens_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load raster
sampling_prob <- raster::raster("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/sampling_prob.tif")

# Mask to MCP
spCrop <- raster::crop(sampling_prob, bgExt_Cb)
spMask <- raster::mask(spCrop, bgExt_Cb)

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), 10000, prob = values(spMask)))

bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Campephilus pollens/bgpoints_samplingprob.shp")


# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Campephilus pollens/evalTable_thinned_samplingProb_Bg.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Campephilus pollens/cloglog_MCP_thinned_samplingProb_AUC.tif")


# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Campephilus_pollens_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Campephilus pollens/projection_thinned_samplingProb_AUC.tif")

# deltaAICc = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Campephilus pollens/cloglog_MCP_thinned_samplingProb_AIC.tif")

# Generate projection

# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.0.5",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Campephilus pollens/projection_thinned_samplingProb_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))

# Odontophorus hyperythrus  - 132 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Odontophorus hyperythrus") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Odontophorus hyperythrus_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load raster
sampling_prob <- raster::raster("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/sampling_prob.tif")

# Mask to MCP
spCrop <- raster::crop(sampling_prob, bgExt_Cb)
spMask <- raster::mask(spCrop, bgExt_Cb)

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), 10000, prob = values(spMask)))

bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Odontophorus hyperythrus/bgpoints_samplingprob.shp")


# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Odontophorus hyperythrus/evalTable_thinned_samplingProb_Bg.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Odontophorus hyperythrus/cloglog_MCP_thinned_samplingProb_AUC.tif")


# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Odontophorus_hyperythrus_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.2.5",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Odontophorus hyperythrus/projection_thinned_samplingProb_AUC.tif")

# deltaAICc = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQ_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Odontophorus hyperythrus/cloglog_MCP_thinned_samplingProb_AIC.tif")

# Generate projection

# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQ_rm.0.5",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Odontophorus hyperythrus/projection_thinned_samplingProb_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))


# Myioborus flavivertex  - 16 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Myioborus flavivertex") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Myioborus flavivertex_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load raster
sampling_prob <- raster::raster("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/sampling_prob.tif")

# Mask to MCP
spCrop <- raster::crop(sampling_prob, bgExt_Cb)
spMask <- raster::mask(spCrop, bgExt_Cb)

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), 1000, prob = values(spMask)))

bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Myioborus flavivertex/bgpoints_samplingprob.shp")


# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Myioborus flavivertex/evalTable_thinned_samplingProb_Bg.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQH_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Myioborus flavivertex/cloglog_MCP_thinned_samplingProb_AUC.tif")


# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN/Myioborus_flavivertex_Ayerbe_IUCN_buffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Myioborus flavivertex/projection_thinned_samplingProb_AUC.tif")

# deltaAICc = fc.L_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Myioborus flavivertex/cloglog_MCP_thinned_samplingProb_AIC.tif")

# Generate projection

# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.1.5",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Myioborus flavivertex/projection_thinned_samplingProb_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))

# Mitu salvini  - 22 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Mitu salvini") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Mitu salvini_mcp_10km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load raster
sampling_prob <- raster::raster("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/sampling_prob.tif")

# Mask to MCP
spCrop <- raster::crop(sampling_prob, bgExt_Cb)
spMask <- raster::mask(spCrop, bgExt_Cb)

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), 10000, prob = values(spMask)))

bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Mitu salvini/bgpoints_samplingprob.shp")


# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 2, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 


# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Mitu salvini/evalTable_thinned_samplingProb_Bg.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.2

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Mitu salvini/cloglog_MCP_thinned_samplingProb_AUC.tif")


# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Biomodelos_and_IUCN/Mitu_salvini_merge_10kmbuffer.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Mitu salvini/projection_thinned_samplingProb_AUC.tif")

# deltaAICc = fc.H_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Mitu salvini/cloglog_MCP_thinned_samplingProb_AIC.tif")

# Generate projection

# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2.5",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/Bias/Mitu salvini/projection_thinned_samplingProb_AIC.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))


# FINAL SAMPLING PROBABILITY MAP -----------------------------------------------------------------------------

#### Campephilus pollens  - 172 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Campephilus pollens") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Campephilus pollens_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
####### Random background  ---------------------------------------
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 10000)

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, bgMask_Cb, n=1) 


####### Sampling Prob  ---------------------------------------

sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability

# Calculate number of sampling points within MCP
# Load Aves sampling sites
Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Sample points depending on survey probability 
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/bgpoints_SamplingProbNumberOfSurveyedSites.shp")

####### TGB ---------------------------------------

Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# EXPORT
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/bgpoints_TGB_Background.shp")


####### Model ---------------------------------------

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 

# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/evalTable_SamplingProb_Background.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Campephilus pollens/cloglog_MCP_SamplingProb_HF.tif")


# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.L_rm.1",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Campephilus pollens/AGOL/projection_HF_distance_random.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_Cb[ , c("longitude", "latitude")])

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]

# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA


# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/projection_p10_SamplingProb_Background.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))

#### Mitu salvini  - 22 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Mitu salvini") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Mitu salvini_mcp_10km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
# Load raster
sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/output/HF_model.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), 10000, prob = values(spMask)))

bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Mitu salvini/bgpoints_samplingprobHF_raw.shp")

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 2, color = "red") 

# TGB
TGB <- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(TGB) <- ~longitude+latitude
sp::proj4string(TGB) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
TGB <- gIntersection(bgExt_Cb, TGB) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(TGB, envs_Cb$X30s_bio1_LAC, n=1) 

# EXPORT
bgSample_Cb <- as.data.frame(bgSample_Cb) %>% rename(longitude = x, latitude = y)

sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Mitu salvini/bgpoints_TGB.shp")

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 




# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Mitu salvini/evalTable_TGB.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Mitu salvini/cloglog_MCP_SamplingProb_HF.tif")


# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.2.5",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Mitu salvini/projection_TGB.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_xy_Cb)

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]

# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA


# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Mitu salvini/projection_TGB_p10.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))


#### Myioborus flavivertex  - 16 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Myioborus flavivertex") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Myioborus flavivertex_mcp_10km_thin5km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
####### Random background  ---------------------------------------
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 1000)

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, bgMask_Cb, n=1) 


####### Sampling Prob  ---------------------------------------

sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability

# Calculate number of sampling points within MCP
# Load Aves sampling sites
Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Sample points depending on survey probability 
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Myioborus flavivertex/bgpoints_SamplingProbNumberOfSurveyedSites.shp")

####### TGB ---------------------------------------

Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# EXPORT
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/bgpoints_TGB_Background.shp")


####### Model ---------------------------------------

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 



# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Myioborus flavivertex/evalTable_SamplingProb_Background.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.1

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Myioborus flavivertex/cloglog_MCP_SamplingProb_HF.tif")


# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Myioborus flavivertex/projection_HF_distance.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_Cb[ , c("longitude", "latitude")])

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]

# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA


# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Myioborus flavivertex/projection_p10_SamplingProb_Background.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))

#### Odontophorus hyperythrus  - 132 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Odontophorus hyperythrus") %>% droplevels() 

## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Odontophorus hyperythrus_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
####### Random background  ---------------------------------------
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 10000)

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 


####### Sampling Prob  ---------------------------------------

sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability

# Calculate number of sampling points within MCP
# Load Aves sampling sites
Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Sample points depending on survey probability 
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Odontophorus hyperythrus/bgpoints_SamplingProbNumberOfSurveyedSites.shp")

####### TGB ---------------------------------------

Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# EXPORT
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/bgpoints_TGB_Background.shp")


####### Model ---------------------------------------

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Odontophorus hyperythrus/evalTable_SamplingProb_Background.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.LQ_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.3.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Odontophorus hyperythrus/cloglog_MCP_SamplingProb_HF.tif")


# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")

# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.3.5",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Odontophorus hyperythrus/projection_HF_distance.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_Cb[ , c("longitude", "latitude")])

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]

# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA

# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Odontophorus hyperythrus/projection_p10_SamplingProb_Background.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))


#### Crypturellus berlepschi - 23 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Crypturellus berlepschi") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Crypturellus berlepschi_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
####### Random background  ---------------------------------------
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 10000)

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, bgMask_Cb, n=1) 


####### Sampling Prob  ---------------------------------------

sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability
# PARA TIERRAS BAJAS (DADO QUE HAY POCAS ZONAS CON UNA ALTA PROBABILIDAD DE MUESTREO) EL NÚMERO DE PUNTOS A SER MUESTREADOS 
# LOS VAMOS A DEFINIR COMO EL NÚMERO DE PUNTOS DE MUESTREO QUE SE HAYAN REGISTRADO DENTRO DEL MCP

# Calculate number of sampling points within MCP
# Load Aves sampling sites
Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Sample points depending on survey probability 
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Crypturellus berlepschi/bgpoints_SamplingProbNumberOfSurveyedSites.shp")

####### TGB ---------------------------------------

Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(TGB, envs_Cb$X30s_bio1_LAC, n=1) 

# EXPORT
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Campephilus pollens/bgpoints_TGB.shp")



####### Model ---------------------------------------

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Crypturellus berlepschi/evalTable_SamplingProbNumberOfSurveyedSites_Background.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.1.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Crypturellus berlepschi/cloglog_MCP_SamplingProb_HF.tif")

# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")


# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Crypturellus berlepschi/projection_TGB2.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_Cb[ , c("longitude", "latitude")])

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]


# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA


# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Crypturellus berlepschi/projection_p10_SamplingProbNumberOfSurveyedSites_Background.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))


#### Hypopyrrhus pyrohypogaster  - 121 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Hypopyrrhus pyrohypogaster") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Hypopyrrhus pyrohypogaster_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
####### Random background  ---------------------------------------
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 10000)

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 


####### Sampling Prob  ---------------------------------------

sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability

# Calculate number of sampling points within MCP
# Load Aves sampling sites
Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Sample points depending on survey probability 
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Hypopyrrhus pyrohypogaster/bgpoints_SamplingProb.shp")

####### TGB ---------------------------------------

Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# EXPORT
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/bgpoints_TGB_Background.shp")


####### Model ---------------------------------------
# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb2",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Hypopyrrhus pyrohypogaster/evalTable_SamplingProb.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.2.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.H_rm.4"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Crypturellus berlepschi/cloglog_MCP_SamplingProb_HF.tif")

# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")


# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.H_rm.4",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Hypopyrrhus pyrohypogaster/projection_TGB.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_Cb[ , c("longitude", "latitude")])

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]

# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA


# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Hypopyrrhus pyrohypogaster/projection_p10_SamplingProb_Background.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))


#### Ognorhynchus icterotis  - 39 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Ognorhynchus icterotis") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Ognorhynchus icterotis_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
####### Random background  ---------------------------------------
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 10000)

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 


####### Sampling Prob  ---------------------------------------

sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability

# Calculate number of sampling points within MCP
# Load Aves sampling sites
Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Sample points depending on survey probability 
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Ognorhynchus icterotis/bgpoints_SamplingProb.shp")

####### TGB ---------------------------------------

Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# EXPORT
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/bgpoints_TGB_Background.shp")


####### Model ---------------------------------------
# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Ognorhynchus icterotis/evalTable_SamplingProb_Background.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.3

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.2.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Crypturellus berlepschi/cloglog_MCP_SamplingProb_HF.tif")

# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")


# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Create object of projection variables
projAreaEnvs_Cb <- envs_Cb
# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.2.5",
  envs = projAreaEnvs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Ognorhynchus icterotis/projection_HF_raw.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_Cb[ , c("longitude", "latitude")])

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]

# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA


# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Ognorhynchus icterotis/projection_p10_SamplingProb_Background.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))


#### Patagioenas goodsoni  - 71 records ---------------------------------------------


occs_Cb <- data_thinned %>% filter(name == "Patagioenas goodsoni") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Patagioenas goodsoni_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
####### Random background  ---------------------------------------
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 10000)

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 


####### Sampling Prob  ---------------------------------------

sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability

# Calculate number of sampling points within MCP
# Load Aves sampling sites
Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Sample points depending on survey probability 
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Patagioenas goodsoni/bgpoints_SamplingProb_Background.shp")

####### TGB ---------------------------------------

Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# EXPORT
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/bgpoints_TGB_Background.shp")


####### Model ---------------------------------------

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Patagioenas goodsoni/evalTable_SamplingProb_Background.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Crypturellus berlepschi/cloglog_MCP_SamplingProb_HF.tif")

# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")


# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Patagioenas goodsoni/projection_HF_raw.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_Cb[ , c("longitude", "latitude")])

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]

# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA


# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Patagioenas goodsoni/projection_p10_SamplingProb_Background.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))


#### Grallaria bangsi- 18 records ---------------------------------------------

occs_Cb <- data_thinned %>% filter(name == "Grallaria bangsi") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria bangsi_mcp_10km.shp")

# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
####### Random background  ---------------------------------------
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 1000)

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 


####### Sampling Prob  ---------------------------------------

sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability

# Calculate number of sampling points within MCP
# Load Aves sampling sites
Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Sample points depending on survey probability 
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Grallaria bangsi/bgpoints_SamplingProb_Background.shp")

####### TGB ---------------------------------------

Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# EXPORT
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/bgpoints_TGB_Background.shp")


####### Model ---------------------------------------

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "jack",
  kfolds = nrow(occs_Cb)) 

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Grallaria bangsi/evalTable_SamplingProb_Background.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.H_rm.4

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.LQH_rm.1.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Crypturellus berlepschi/cloglog_MCP_SamplingProb_HF.tif")

# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")


# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 



# Generate a projection of the model to the desired area
proj_area_Cb <- proj_area(
  evalOut = model_Cb,
  curModel = "fc.LQH_rm.1.5",
  envs = envs_Cb , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Cb) 

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Grallaria bangsi/projection_TGB.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_Cb[ , c('longitude', 'latitude')])

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]

# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA

# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Grallaria bangsi/projection_p10_SamplingProb_Background.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))

# Grallaria flavotincta - 31 records ---------------------------------------------


occs_Cb <- data_thinned %>% filter(name == "Grallaria flavotincta") %>% droplevels() 


## Specify the directory with the environmental variables
dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)
occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
# remove occurrence records with NA environmental values
occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Cb <- na.omit(occs_vals_Cb)
# add columns for env variable values for each occurrence record
occs_Cb <- cbind(occs_Cb, occs_vals_Cb)

# Load the user provided shapefile or csv file with the desired extent.
##User must input the path to shapefile or csv file and the file name 
bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Grallaria flavotincta_mcp_10km_thin5km.shp")


# Check that occurrence points and MCP intersect
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Mask environmental data to provided extent
bgMask_Cb <- penvs_bgMask(
  occs = occs_Cb,
  envs = envs_Cb,
  bgExt = bgExt_Cb)

# Define background points
####### Random background  ---------------------------------------
bgSample_Cb <- penvs_bgSample(
  occs = occs_Cb,
  bgMask =  bgMask_Cb,
  bgPtsNum = 10000)

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 


####### Sampling Prob  ---------------------------------------

sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")

# Mask to MCP
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample points depending on survey probability

# Calculate number of sampling points within MCP
# Load Aves sampling sites
Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Sample points depending on survey probability 
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# Turn into points and export
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Patagioenas goodsoni/bgpoints_SamplingProb_Background.shp")

####### TGB ---------------------------------------

Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")

# Turn to spatial points dataframe 
sp::coordinates(Sampling_sites) <- ~longitude+latitude
sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)

# Extract points that fall within calibration area
Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 

# Remove points that fall in the same cell
bgSample_Cb <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 

# Turn to dataframe 
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)

# EXPORT
sf::st_as_sf(bgSample_Cb, coords = c("longitude", "latitude"), crs = 4326) %>% 
  sf::st_write("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Campephilus pollens/bgpoints_TGB_Background.shp")


####### Model ---------------------------------------

# Plot
ggplot() + 
  geom_polygon(data = bgExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = bgSample_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


# Extract values of environmental layers for each background point
bgEnvsVals_Cb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb))

# Remove NA values 
bgSample_Cb <- bgSample_Cb[!(rowSums(is.na(bgEnvsVals_Cb)) > 1), ]
bgEnvsVals_Cb <- na.omit(bgEnvsVals_Cb)

##Add extracted values to background points table
bgEnvsVals_Cb <- cbind(bgSample_Cb,  bgEnvsVals_Cb)

# R code to get partitioned data
groups_Cb <- part_partitionOccs(
  occs = occs_Cb ,
  bg =  bgSample_Cb, 
  method = "cb1",
  bgMask = bgMask_Cb,
  aggFact = 2)

# Run maxent model for the selected species
model_Cb <- model_maxent(
  occs = occs_Cb,
  bg = bgEnvsVals_Cb,
  user.grp = groups_Cb, 
  bgMsk = bgMask_Cb,
  rms = c(0.5, 4), 
  rmsStep =  0.5,
  fcs = c('L', 'LQ', 'H', 'LQH'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 6)

# Overall results
res <- eval.results(model_Cb)

# export
write.csv(res, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Patagioenas goodsoni/evalTable_SamplingProb_Background.csv", row.names = FALSE)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC = fc.L_rm.0.5

# Select model and obtain raster prediction
m_Cb <- model_Cb@models[["fc.L_rm.0.5"]]
predSel_Cb <- predictMaxnet(m_Cb, bgMask_Cb,
                            type = "cloglog", 
                            clamp = TRUE)

writeRaster(predSel_Cb, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Crypturellus berlepschi/cloglog_MCP_SamplingProb_HF.tif")

# Load projection area 
proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")


# Check that points fall within the distribution
ggplot() + 
  geom_polygon(data = proj_userExt_Cb, aes(x = long, y = lat), colour='black', fill='white') +
  geom_point(data = occs_Cb, aes(x=longitude, y=latitude), size = 4, color = "red") 


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

writeRaster(proj_area_Cb$projArea, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Patagioenas goodsoni/projection_HF_raw.tif")

# Extract SDM predictions at occurrence points 
occPredVals <- raster::extract(predSel_Cb, occs_Cb[ , c("longitude", "latitude")])

# Calculate 10th percentile 
p10 <- ceiling(length(occPredVals) * 0.9)
thresh <- rev(sort(occPredVals))[p10]

# Apply threshold to prediction
p10_prediction <- proj_area_Cb$projArea
p10_prediction[Which(p10_prediction<thresh)]<-NA


# Export
writeRaster(p10_prediction, "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Patagioenas goodsoni/projection_p10_SamplingProb_Background.tif")

rm(list=setdiff(ls(), c("envs_Cb", "data_thinned")))



######## T TEST ####################### 

bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/MCPs/Crypturellus berlepschi_mcp_10km_thin5km.shp")
sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_model/output/HF_model.tif")
spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))

# Replace NAs with zeros
spMask[is.na(spMask[])] <- 0 

# Sample bg points
bgSample_Cb <- xyFromCell(spMask, sample(ncell(spMask), 10000, prob = values(spMask)))
# Remove points in the same cell
bgSample_Cb <- gridSample(bgSample_Cb, spMask)
# Turn to dataframe
bgSample_Cb <- as.data.frame(bgSample_Cb) %>%
  rename(longitude = x, latitude = y)



# Create a point per cell 
one <- rasterToPoints(spMask, spatial = TRUE)

# Clip points
two <- gIntersection(one, bgExt_Cb)  %>% 
  as.data.frame() %>% rename(longitude = x, latitude = y) %>% 
  dplyr::select(longitude, latitude)

# Remove sampled points
three <- rbind(two, bgSample_Cb)
# Round decimals
three <- three %>% mutate(longitude = round(longitude, digits = 4), latitude = round(latitude, digits = 4))
# Remove duplicates
dups <-which(duplicated(three)|duplicated(three, fromLast=TRUE))
four <- three[-(dups), ]

# Extract HF values per point
HF_vals_notS <- as.data.frame(raster::extract(spMask, four)) %>% cbind(four) %>% rename(Value = `raster::extract(spMask, four)`) %>% filter(Value != 0) %>% mutate(group = "Not surveyed")
HF_vals_S <- as.data.frame(raster::extract(spMask, bgSample_Cb)) %>% cbind(bgSample_Cb)  %>% rename(Value = `raster::extract(spMask, bgSample_Cb)`) %>% mutate(group = "Surveyed")

# Summary statistics and T test
Total_points <- rbind(HF_vals_notS, HF_vals_S)

Total_points %>%
  group_by(group) %>%
  summarise(mean = mean(Value), 
            sd = sd(Value))
  
ggboxplot(
  Total_points, x = "group", y = "Value", 
  ylab = "Value", xlab = "Groups", add = "jitter"
)

# Check assumptions 
# Normality
ad.test(Total_points$Value[Total_points$group == "Not surveyed"])
ad.test(Total_points$Value[Total_points$group == "Surveyed"])

ggqqplot(Total_points, x = "Value", facet.by = "group")

# Homocedasticity
Total_points %>% levene_test(Value ~ group)

# T test
stat.test <- Total_points %>% 
  t_test(Value ~ group) %>%
  add_significance()

test <- wilcox.test(Total_points$Value ~ Total_points$group,
                    alternative = "less")
test


# So, the sampled sites have a statistically higuer mean of HF values 





