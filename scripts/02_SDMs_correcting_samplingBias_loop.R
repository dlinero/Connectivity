library(tidyverse)
library(wallace)
library(dismo)
library(rgeos)
library(ENMeval)
library(janitor) # For cleanning dataset column names
library(adehabitatHR) # For minimum convex polygons 




data_thinned <- read.csv("F:/Connectivity/outputs/02_SDMs/MCPs/FINAL_THINNED_DATA.csv")

dir_envs_Cb <- "F:/Connectivity/outputs/02_SDMs/wc"
envs_path <- file.path(dir_envs_Cb, c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'))
# Create environmental object 
envs_Cb <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c('30s_bio1_LAC.tif', '30s_bio2_LAC.tif', '30s_bio3_LAC.tif', '30s_bio4_LAC.tif', '30s_bio5_LAC.tif', '30s_bio6_LAC.tif', '30s_bio7_LAC.tif', '30s_bio8_LAC.tif', '30s_bio9_LAC.tif', '30s_bio10_LAC.tif', '30s_bio11_LAC.tif', '30s_bio12_LAC.tif', '30s_bio13_LAC.tif', '30s_bio14_LAC.tif', '30s_bio15_LAC.tif', '30s_bio16_LAC.tif', '30s_bio17_LAC.tif', '30s_bio18_LAC.tif', '30s_bio19_LAC.tif'),
  doBrick = TRUE)

species_names <- c("Hypopyrrhus pyrohypogaster")



for (i in 1:length(species_names)){
  
  print(species_names[i])
  
  occs_Cb <- data_thinned %>% filter(name == species_names[i]) %>% droplevels() 
  
  
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
  filename1 <- "F:/Connectivity/outputs/02_SDMs/MCPs/"
  filename2 <- as.character(species_names[i])
  filename3 <- "_mcp_10km_thin5km.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  bgExt_Cb <- rgdal::readOGR(complete_filename)
  
  
  # Mask environmental data to provided extent
  bgMask_Cb <- penvs_bgMask(
    occs = occs_Cb,
    envs = envs_Cb,
    bgExt = bgExt_Cb)
  
  # Define background points
  ####### Random background  ---------------------------------------
  print("Random background")
  
  bgSample_Cb_rb <- penvs_bgSample(
    occs = occs_Cb,
    bgMask =  bgMask_Cb,
    bgPtsNum = 10000)
  
  # Remove points that fall in the same cell
  bgSample_Cb_rb <- gridSample(bgSample_Cb_rb, envs_Cb$X30s_bio1_LAC, n=1) 
  
  
  ####### Sampling Prob  ---------------------------------------
  
  print("Sampling prob")

  sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")
  
  # Mask to MCP
  spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))
  
  # Replace NAs with zeros
  spMask[is.na(spMask[])] <- 0 
  
  # Sample points depending on survey probability 
  bgSample_Cb_sp <- xyFromCell(spMask, sample(ncell(spMask), 10000, prob = values(spMask)))
  
  # Remove points that fall in the same cell
  bgSample_Cb_sp <- gridSample(bgSample_Cb_sp, envs_Cb$X30s_bio1_LAC, n=1) 
  
  # Turn to dataframe 
  bgSample_Cb_sp <- as.data.frame(bgSample_Cb_sp) %>%
    rename(longitude = x, latitude = y)
  
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/bgpoints_SamplingProb_Background.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  # Turn into points and export
  sf::st_as_sf(bgSample_Cb_sp, coords = c("longitude", "latitude"), crs = 4326) %>% 
    sf::st_write(complete_filename)
  
  ####### Sampling Prob - TGB number of points  ---------------------------------------
  
  print("Sampling prob - TGB points")
   
  # Calculate number of sampling points within MCP
  # Load Aves sampling sites
  Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")
  
  # Turn to spatial points dataframe 
  sp::coordinates(Sampling_sites) <- ~longitude+latitude
  sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)
  
  # Extract points that fall within calibration area
  print("intersection")
  Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 
  
  # Remove points that fall in the same cell
  TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 
  
  # Sample points depending on survey probability 
  bgSample_Cb_sp_tgb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))
  
  # Remove points that fall in the same cell
  bgSample_Cb_sp_tgb <- gridSample(bgSample_Cb_sp_tgb, envs_Cb$X30s_bio1_LAC, n=1) 
  
  # Turn to dataframe 
  bgSample_Cb_sp_tgb <- as.data.frame(bgSample_Cb_sp_tgb) %>%
    rename(longitude = x, latitude = y)
  
  # Turn into points and export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/bgpoints_SamplingProbNumberOfSurveyedSites_Background.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  sf::st_as_sf(bgSample_Cb_sp_tgb, coords = c("longitude", "latitude"), crs = 4326) %>% 
    sf::st_write(complete_filename)
  
  ####### TGB ---------------------------------------
  print("TGB")
  
  # Turn to dataframe 
  bgSample_Cb_TGB <- as.data.frame(TGB) %>%
    rename(longitude = x, latitude = y) %>%
    dplyr::select(longitude, latitude)
  
  # EXPORT
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/bgpoints_TGB_Background.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  sf::st_as_sf(bgSample_Cb_TGB, coords = c("longitude", "latitude"), crs = 4326) %>% 
    sf::st_write(complete_filename)
  
  
  ####### Model ---------------------------------------
  # Extract values of environmental layers for each background point
  print("extracting environmental values")
  bgEnvsVals_Cb_rb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_rb))
  bgEnvsVals_Cb_sp <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_sp))
  bgEnvsVals_Cb_sp_tgb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_sp_tgb))
  bgEnvsVals_Cb_TGB <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_TGB))
  
  # Remove NA values 
  bgSample_Cb_rb <- bgSample_Cb_rb[!(rowSums(is.na(bgEnvsVals_Cb_rb)) > 1), ]
  bgEnvsVals_Cb_rb <- na.omit(bgEnvsVals_Cb_rb)
  
  bgSample_Cb_sp <- bgSample_Cb_sp[!(rowSums(is.na(bgEnvsVals_Cb_sp)) > 1), ]
  bgEnvsVals_Cb_sp <- na.omit(bgEnvsVals_Cb_sp)
  
  bgSample_Cb_sp_tgb <- bgSample_Cb_sp_tgb[!(rowSums(is.na(bgEnvsVals_Cb_sp_tgb)) > 1), ]
  bgEnvsVals_Cb_sp_tgb <- na.omit(bgEnvsVals_Cb_sp_tgb)
  
  bgSample_Cb_TGB <- bgSample_Cb_TGB[!(rowSums(is.na(bgEnvsVals_Cb_TGB)) > 1), ]
  bgEnvsVals_Cb_TGB <- na.omit(bgEnvsVals_Cb_TGB)
  
  ##Add extracted values to background points table
  bgEnvsVals_Cb_rb <- cbind(bgSample_Cb_rb,  bgEnvsVals_Cb_rb)
  bgEnvsVals_Cb_sp <- cbind(bgSample_Cb_sp,  bgEnvsVals_Cb_sp)
  bgEnvsVals_Cb_sp_tgb <- cbind(bgSample_Cb_sp_tgb,  bgEnvsVals_Cb_sp_tgb)
  bgEnvsVals_Cb_TGB <- cbind(bgSample_Cb_TGB,  bgEnvsVals_Cb_TGB)
  
  # R code to get partitioned data
  groups_Cb_rb <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgSample_Cb_rb, 
    method = "cb2",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  groups_Cb_sp <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgSample_Cb_sp, 
    method = "cb2",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  groups_Cb_sp_tgb <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgSample_Cb_sp_tgb, 
    method = "cb2",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  groups_Cb_TGB <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgEnvsVals_Cb_TGB, 
    method = "cb2",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  # Run maxent model for the selected species
  print("first model")
  model_Cb_rb <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_rb,
    user.grp = groups_Cb_rb, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ', 'H', 'LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_rb <- eval.results(model_Cb_rb)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/evalTable_Random_Background.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_rb, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_rb <- min(res_rb$or.10p.avg)
  max_auc_rb <- max(res_rb$auc.val.avg[which(res_rb$or.10p.avg == min_or_rb)])
  best_model_rb <- as.character(res_rb$tune.args[res_rb$auc.val.avg == max_auc_rb & res_rb$or.10p.avg == min_or_rb])
  
  # Select model and obtain raster prediction
  m_Cb_rb <- model_Cb_rb@models[[best_model_rb]]
  predSel_Cb_rb <- predictMaxnet(m_Cb_rb, bgMask_Cb,
                              type = "cloglog", 
                              clamp = TRUE)
  
  # Load projection area 
  proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_rb <- proj_area(
    evalOut = model_Cb_rb,
    curModel = best_model_rb,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_rb <- raster::extract(predSel_Cb_rb, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_rb <- ceiling(length(occPredVals_rb) * 0.9)
  thresh_rb <- rev(sort(occPredVals_rb))[p10_rb]
  
  # Apply threshold to prediction
  p10_prediction_rb <- proj_area_Cb_rb$projArea
  p10_prediction_rb[Which(p10_prediction_rb<thresh_rb)]<-NA
  
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/projection_p10_Random_Background.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_rb, complete_filename)
  
  print("second model")
  model_Cb_sp <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_sp,
    user.grp = groups_Cb_sp, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ','H',  'LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_sp <- eval.results(model_Cb_sp)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/evalTable_SamplingProb_Background.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_sp, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_sp <- min(res_sp$or.10p.avg)
  max_auc_sp <- max(res_sp$auc.val.avg[which(res_sp$or.10p.avg == min_or_sp)])
  best_model_sp <- as.character(res_sp$tune.args[res_sp$auc.val.avg == max_auc_sp & res_sp$or.10p.avg == min_or_sp])
  
  # Select model and obtain raster prediction
  m_Cb_sp <- model_Cb_sp@models[[best_model_sp]]
  predSel_Cb_sp <- predictMaxnet(m_Cb_sp, bgMask_Cb,
                                 type = "cloglog", 
                                 clamp = TRUE)
  

  # Generate a projection of the model to the desired area
  proj_area_Cb_sp <- proj_area(
    evalOut = model_Cb_sp,
    curModel = best_model_sp,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_sp <- raster::extract(predSel_Cb_sp, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_sp <- ceiling(length(occPredVals_sp) * 0.9)
  thresh_sp <- rev(sort(occPredVals_sp))[p10_sp]
  
  # Apply threshold to prediction
  p10_prediction_sp <- proj_area_Cb_sp$projArea
  p10_prediction_sp[Which(p10_prediction_sp<thresh_sp)]<-NA
  
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/projection_p10_SamplingProb_Background.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_sp, complete_filename)
  
  print("third model")
  model_Cb_sp_tgb <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_sp_tgb,
    user.grp = groups_Cb_sp_tgb, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ',  'H','LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_sp_tgb <- eval.results(model_Cb_sp_tgb)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/evalTable_SamplingProbNumberOfSurveyedSites_Background.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_sp_tgb, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_sp_tgb <- min(res_sp_tgb$or.10p.avg)
  max_auc_sp_tgb <- max(res_sp_tgb$auc.val.avg[which(res_sp_tgb$or.10p.avg == min_or_sp_tgb)])
  best_model_sp_tgb <- as.character(res_sp_tgb$tune.args[res_sp_tgb$auc.val.avg == max_auc_sp_tgb & res_sp_tgb$or.10p.avg == min_or_sp_tgb])
  
  # Select model and obtain raster prediction
  m_Cb_sp_tgb <- model_Cb_sp_tgb@models[[best_model_sp_tgb]]
  predSel_Cb_sp_tgb <- predictMaxnet(m_Cb_sp_tgb, bgMask_Cb,
                                 type = "cloglog", 
                                 clamp = TRUE)
  
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_sp_tgb <- proj_area(
    evalOut = model_Cb_sp_tgb,
    curModel = best_model_sp_tgb,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_sp_tgb <- raster::extract(predSel_Cb_sp_tgb, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_sp_tgb <- ceiling(length(occPredVals_sp_tgb) * 0.9)
  thresh_sp_tgb <- rev(sort(occPredVals_sp_tgb))[p10_sp_tgb]
  
  # Apply threshold to prediction
  p10_prediction_sp_tgb <- proj_area_Cb_sp_tgb$projArea
  p10_prediction_sp_tgb[Which(p10_prediction_sp_tgb<thresh_sp_tgb)]<-NA
  
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/projection_p10_SamplingProbNumberOfSurveyedSites_Background.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_sp_tgb, complete_filename)
  
  
  print("fourth model")
  model_Cb_TGB <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_TGB,
    user.grp = groups_Cb_TGB, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ',  'H','LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_TGB <- eval.results(model_Cb_TGB)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/evalTable_TGB_Background.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_TGB, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_TGB <- min(res_TGB$or.10p.avg)
  max_auc_TGB <- max(res_TGB$auc.val.avg[which(res_TGB$or.10p.avg == min_or_TGB)])
  best_model_TGB <- as.character(res_TGB$tune.args[res_TGB$auc.val.avg == max_auc_TGB & res_TGB$or.10p.avg == min_or_TGB])

  # Select model and obtain raster prediction
  m_Cb_TGB <- model_Cb_TGB@models[[best_model_TGB]]
  predSel_Cb_TGB <- predictMaxnet(m_Cb_TGB, bgMask_Cb,
                                     type = "cloglog", 
                                     clamp = TRUE)
  
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_TGB <- proj_area(
    evalOut = model_Cb_TGB,
    curModel = best_model_TGB,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_TGB <- raster::extract(predSel_Cb_TGB, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_TGB <- ceiling(length(occPredVals_TGB) * 0.9)
  thresh_TGB <- rev(sort(occPredVals_TGB))[p10_TGB]
  
  # Apply threshold to prediction
  p10_prediction_TGB <- proj_area_Cb_TGB$projArea
  p10_prediction_TGB[Which(p10_prediction_TGB<thresh_TGB)]<-NA
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- as.character(species_names[i])
  filename3 <- "/projection_p10_TGB_Background_fc_H_rm_4.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_TGB, complete_filename)
  
  gc()
  
}


###### Pyrrhura melanura ----------------------------------------------------------------------

# Load cleanned, thinned data for Pyrrhura melanura, that differentiates among subspecies 
data <- read.csv("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/Pyrrhura melanura/correction/FINAL_THINNED_DATA_Pyrrhura_melanura_sbsp.csv")

pacifica <- data %>% filter(subspecie == "pacifica") %>% droplevels()
chapmani <- data %>% filter(subspecie == "chapmani") %>% droplevels()
souancei_melanura <- data %>% filter(subspecie == "souancei/melanura") %>% droplevels()


# Create MCP with a buffer of 0.1° for each subspecies °

# Extract only the coordinates
pacifica_MCP <- pacifica[c('longitude', 'latitude')]
chapmani_MCP <- chapmani[c('longitude', 'latitude')]
souancei_melanura_MCP <- souancei_melanura[c('longitude', 'latitude')]
# Assign the columns with the coordinates
sp::coordinates(pacifica_MCP) <- ~ longitude + latitude
sp::coordinates(chapmani_MCP) <- ~ longitude + latitude
sp::coordinates(souancei_melanura_MCP) <- ~ longitude + latitude
# Calculate minimum convex polygon
pacifica_MCP <- mcp(pacifica_MCP, percent = 100)
chapmani_MCP <- mcp(chapmani_MCP, percent = 100)
souancei_melanura_MCP <- mcp(souancei_melanura_MCP, percent = 100)
# Apply buffer 
pacifica_MCP <- rgeos::gBuffer(pacifica_MCP, width = 0.1)
chapmani_MCP <- rgeos::gBuffer(chapmani_MCP, width = 0.1)
souancei_melanura_MCP <- rgeos::gBuffer(souancei_melanura_MCP, width = 0.1)
# Assign coordinates
crs(pacifica_MCP) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(chapmani_MCP) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
crs(souancei_melanura_MCP) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

  
  occs_Cb <- souancei_melanura  
  
  
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
  bgExt_Cb <- souancei_melanura_MCP
  
  
  # Mask environmental data to provided extent
  bgMask_Cb <- penvs_bgMask(
    occs = occs_Cb,
    envs = envs_Cb,
    bgExt = bgExt_Cb)
  
  # Define background points
  ####### Random background  ---------------------------------------
  print("Random background")
  
  bgSample_Cb_rb <- penvs_bgSample(
    occs = occs_Cb,
    bgMask =  bgMask_Cb,
    bgPtsNum = 10000)
  
  # Remove points that fall in the same cell
  bgSample_Cb_rb <- gridSample(bgSample_Cb_rb, envs_Cb$X30s_bio1_LAC, n=1) 
  
  
  ####### Sampling Prob  ---------------------------------------
  
  print("Sampling prob")
  
  sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")
  
  # Mask to MCP
  spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))
  
  # Replace NAs with zeros
  spMask[is.na(spMask[])] <- 0 
  
  # Sample points depending on survey probability 
  bgSample_Cb_sp <- xyFromCell(spMask, sample(ncell(spMask), 10000, prob = values(spMask)))
  
  # Remove points that fall in the same cell
  bgSample_Cb_sp <- gridSample(bgSample_Cb_sp, envs_Cb$X30s_bio1_LAC, n=1) 
  
  # Turn to dataframe 
  bgSample_Cb_sp <- as.data.frame(bgSample_Cb_sp) %>%
    rename(longitude = x, latitude = y)
  
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/bgpoints_SamplingProb_Background_souancei_melanura.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  # Turn into points and export
  sf::st_as_sf(bgSample_Cb_sp, coords = c("longitude", "latitude"), crs = 4326) %>% 
    sf::st_write(complete_filename)
  
  ####### Sampling Prob - TGB number of points  ---------------------------------------
  
  print("Sampling prob - TGB points")
  
  # Calculate number of sampling points within MCP
  # Load Aves sampling sites
  Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")
  
  # Turn to spatial points dataframe 
  sp::coordinates(Sampling_sites) <- ~longitude+latitude
  sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)
  
  # Extract points that fall within calibration area
  print("intersection")
  Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 
  
  # Remove points that fall in the same cell
  TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 
  
  # Sample points depending on survey probability 
  bgSample_Cb_sp_tgb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))
  
  # Remove points that fall in the same cell
  bgSample_Cb_sp_tgb <- gridSample(bgSample_Cb_sp_tgb, envs_Cb$X30s_bio1_LAC, n=1) 
  
  # Turn to dataframe 
  bgSample_Cb_sp_tgb <- as.data.frame(bgSample_Cb_sp_tgb) %>%
    rename(longitude = x, latitude = y)
  
  # Turn into points and export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/bgpoints_SamplingProbNumberOfSurveyedSites_Background_souancei_melanura.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  sf::st_as_sf(bgSample_Cb_sp_tgb, coords = c("longitude", "latitude"), crs = 4326) %>% 
    sf::st_write(complete_filename)
  
  ####### TGB ---------------------------------------
  print("TGB")
  
  # Turn to dataframe 
  bgSample_Cb_TGB <- as.data.frame(TGB) %>%
    rename(longitude = x, latitude = y) %>%
    dplyr::select(longitude, latitude)
  
  # EXPORT
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/bgpoints_TGB_Background_souancei_melanura.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  sf::st_as_sf(bgSample_Cb_TGB, coords = c("longitude", "latitude"), crs = 4326) %>% 
    sf::st_write(complete_filename)
  
  
  ####### Model ---------------------------------------
  # Extract values of environmental layers for each background point
  print("extracting environmental values")
  bgEnvsVals_Cb_rb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_rb))
  bgEnvsVals_Cb_sp <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_sp))
  bgEnvsVals_Cb_sp_tgb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_sp_tgb))
  bgEnvsVals_Cb_TGB <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_TGB))
  
  # Remove NA values 
  bgSample_Cb_rb <- bgSample_Cb_rb[!(rowSums(is.na(bgEnvsVals_Cb_rb)) > 1), ]
  bgEnvsVals_Cb_rb <- na.omit(bgEnvsVals_Cb_rb)
  
  bgSample_Cb_sp <- bgSample_Cb_sp[!(rowSums(is.na(bgEnvsVals_Cb_sp)) > 1), ]
  bgEnvsVals_Cb_sp <- na.omit(bgEnvsVals_Cb_sp)
  
  bgSample_Cb_sp_tgb <- bgSample_Cb_sp_tgb[!(rowSums(is.na(bgEnvsVals_Cb_sp_tgb)) > 1), ]
  bgEnvsVals_Cb_sp_tgb <- na.omit(bgEnvsVals_Cb_sp_tgb)
  
  bgSample_Cb_TGB <- bgSample_Cb_TGB[!(rowSums(is.na(bgEnvsVals_Cb_TGB)) > 1), ]
  bgEnvsVals_Cb_TGB <- na.omit(bgEnvsVals_Cb_TGB)
  
  ##Add extracted values to background points table
  bgEnvsVals_Cb_rb <- cbind(bgSample_Cb_rb,  bgEnvsVals_Cb_rb)
  bgEnvsVals_Cb_sp <- cbind(bgSample_Cb_sp,  bgEnvsVals_Cb_sp)
  bgEnvsVals_Cb_sp_tgb <- cbind(bgSample_Cb_sp_tgb,  bgEnvsVals_Cb_sp_tgb)
  bgEnvsVals_Cb_TGB <- cbind(bgSample_Cb_TGB,  bgEnvsVals_Cb_TGB)
  
  # R code to get partitioned data
  groups_Cb_rb <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgSample_Cb_rb, 
    method = "cb2",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  groups_Cb_sp <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgSample_Cb_sp, 
    method = "cb2",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  groups_Cb_sp_tgb <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgSample_Cb_sp_tgb, 
    method = "cb2",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  groups_Cb_TGB <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgEnvsVals_Cb_TGB, 
    method = "cb2",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  # Run maxent model for the selected species
  print("first model")
  model_Cb_rb <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_rb,
    user.grp = groups_Cb_rb, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ', 'H', 'LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_rb <- eval.results(model_Cb_rb)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/evalTable_Random_Background_souancei_melanura.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_rb, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  res_rb <- res_rb %>% filter(or.10p.avg != 0) %>% droplevels()
  min_or_rb <- min(res_rb$or.10p.avg)
  max_auc_rb <- max(res_rb$auc.val.avg[which(res_rb$or.10p.avg == min_or_rb)])
  best_model_rb <- as.character(res_rb$tune.args[res_rb$auc.val.avg == max_auc_rb & res_rb$or.10p.avg == min_or_rb])
  
  # Select model and obtain raster prediction
  m_Cb_rb <- model_Cb_rb@models[[best_model_rb]]
  predSel_Cb_rb <- predictMaxnet(m_Cb_rb, bgMask_Cb,
                                 type = "cloglog", 
                                 clamp = TRUE)
  
  # Load projection area 
  proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_rb <- proj_area(
    evalOut = model_Cb_rb,
    curModel = best_model_rb,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_rb <- raster::extract(predSel_Cb_rb, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_rb <- ceiling(length(occPredVals_rb) * 0.9)
  thresh_rb <- rev(sort(occPredVals_rb))[p10_rb]
  
  # Apply threshold to prediction
  p10_prediction_rb <- proj_area_Cb_rb$projArea
  p10_prediction_rb[Which(p10_prediction_rb<thresh_rb)]<-NA
  
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/projection_p10_Random_Background_souancei_melanura.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_rb, complete_filename)
  
  print("second model")
  model_Cb_sp <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_sp,
    user.grp = groups_Cb_sp, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ','H',  'LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_sp <- eval.results(model_Cb_sp)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/evalTable_SamplingProb_Background_souancei_melanura.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_sp, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  res_sp <- res_sp %>% filter(or.10p.avg != 0) %>% droplevels()
  min_or_sp <- min(res_sp$or.10p.avg)
  max_auc_sp <- max(res_sp$auc.val.avg[which(res_sp$or.10p.avg == min_or_sp)])
  best_model_sp <- as.character(res_sp$tune.args[res_sp$auc.val.avg == max_auc_sp & res_sp$or.10p.avg == min_or_sp])
  
  # Select model and obtain raster prediction
  m_Cb_sp <- model_Cb_sp@models[[best_model_sp]]
  predSel_Cb_sp <- predictMaxnet(m_Cb_sp, bgMask_Cb,
                                 type = "cloglog", 
                                 clamp = TRUE)
  
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_sp <- proj_area(
    evalOut = model_Cb_sp,
    curModel = best_model_sp,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_sp <- raster::extract(predSel_Cb_sp, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_sp <- ceiling(length(occPredVals_sp) * 0.9)
  thresh_sp <- rev(sort(occPredVals_sp))[p10_sp]
  
  # Apply threshold to prediction
  p10_prediction_sp <- proj_area_Cb_sp$projArea
  p10_prediction_sp[Which(p10_prediction_sp<thresh_sp)]<-NA
  
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/projection_p10_SamplingProb_Background_souancei_melanura.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_sp, complete_filename)
  
  print("third model")
  model_Cb_sp_tgb <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_sp_tgb,
    user.grp = groups_Cb_sp_tgb, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ',  'H','LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_sp_tgb <- eval.results(model_Cb_sp_tgb)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/evalTable_SamplingProbNumberOfSurveyedSites_Background_souancei_melanura.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_sp_tgb, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_sp_tgb <- min(res_sp_tgb$or.10p.avg)
  max_auc_sp_tgb <- max(res_sp_tgb$auc.val.avg[which(res_sp_tgb$or.10p.avg == min_or_sp_tgb)])
  best_model_sp_tgb <- as.character(res_sp_tgb$tune.args[res_sp_tgb$auc.val.avg == max_auc_sp_tgb & res_sp_tgb$or.10p.avg == min_or_sp_tgb])
  
  # Select model and obtain raster prediction
  m_Cb_sp_tgb <- model_Cb_sp_tgb@models[[best_model_sp_tgb]]
  predSel_Cb_sp_tgb <- predictMaxnet(m_Cb_sp_tgb, bgMask_Cb,
                                     type = "cloglog", 
                                     clamp = TRUE)
  
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_sp_tgb <- proj_area(
    evalOut = model_Cb_sp_tgb,
    curModel = best_model_sp_tgb,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_sp_tgb <- raster::extract(predSel_Cb_sp_tgb, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_sp_tgb <- ceiling(length(occPredVals_sp_tgb) * 0.9)
  thresh_sp_tgb <- rev(sort(occPredVals_sp_tgb))[p10_sp_tgb]
  
  # Apply threshold to prediction
  p10_prediction_sp_tgb <- proj_area_Cb_sp_tgb$projArea
  p10_prediction_sp_tgb[Which(p10_prediction_sp_tgb<thresh_sp_tgb)]<-NA
  
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/projection_p10_SamplingProbNumberOfSurveyedSites_Background_souancei_melanura.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_sp_tgb, complete_filename)
  
  
  print("fourth model")
  model_Cb_TGB <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_TGB,
    user.grp = groups_Cb_TGB, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ',  'H','LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_TGB <- eval.results(model_Cb_TGB)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/evalTable_TGB_Background_souancei_melanura.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_TGB, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_TGB <- min(res_TGB$or.10p.avg)
  max_auc_TGB <- max(res_TGB$auc.val.avg[which(res_TGB$or.10p.avg == min_or_TGB)])
  best_model_TGB <- as.character(res_TGB$tune.args[res_TGB$auc.val.avg == max_auc_TGB & res_TGB$or.10p.avg == min_or_TGB])
  
  # Select model and obtain raster prediction
  m_Cb_TGB <- model_Cb_TGB@models[[best_model_TGB]]
  predSel_Cb_TGB <- predictMaxnet(m_Cb_TGB, bgMask_Cb,
                                  type = "cloglog", 
                                  clamp = TRUE)
  
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_TGB <- proj_area(
    evalOut = model_Cb_TGB,
    curModel = best_model_TGB,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_TGB <- raster::extract(predSel_Cb_TGB, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_TGB <- ceiling(length(occPredVals_TGB) * 0.9)
  thresh_TGB <- rev(sort(occPredVals_TGB))[p10_TGB]
  
  # Apply threshold to prediction
  p10_prediction_TGB <- proj_area_Cb_TGB$projArea
  p10_prediction_TGB[Which(p10_prediction_TGB<thresh_TGB)]<-NA
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Pyrrhura melanura"
  filename3 <- "/correction/projection_p10_TGB_Background_souancei_melanura.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_TGB, complete_filename)
  
  gc()
  
  ###### Ognorhynchus icterotis ----------------------------------------------------------------------
  

  species_names <- c("Ognorhynchus icterotis")
  
    
    occs_Cb <- data_thinned %>% filter(name == species_names[1]) %>% droplevels() 
    
    
    occs_xy_Cb <- occs_Cb[c('longitude', 'latitude')]
    occs_vals_Cb <- as.data.frame(raster::extract(envs_Cb, occs_xy_Cb))
    # remove occurrence records with NA environmental values
    occs_Cb <- occs_Cb[!(rowSums(is.na(occs_vals_Cb)) > 1), ]
    # also remove variable value rows with NA environmental values
    occs_vals_Cb <- na.omit(occs_vals_Cb)
    # add columns for env variable values for each occurrence record
    occs_Cb <- cbind(occs_Cb, occs_vals_Cb)
    
    # Load the Andes shapefile 
    bgExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/andes.shp")
    # Reproject
    bgExt_Cb <- spTransform(bgExt_Cb,
                                  crs(envs_Cb))
    
  
  # Mask environmental data to provided extent
  bgMask_Cb <- penvs_bgMask(
    occs = occs_Cb,
    envs = envs_Cb,
    bgExt = bgExt_Cb)
  
  # Define background points
  ####### Random background  ---------------------------------------
  print("Random background")
  
  bgSample_Cb_rb <- penvs_bgSample(
    occs = occs_Cb,
    bgMask =  bgMask_Cb,
    bgPtsNum = 10000)
  
  # Remove points that fall in the same cell
  bgSample_Cb_rb <- gridSample(bgSample_Cb_rb, envs_Cb$X30s_bio1_LAC, n=1) 
  
  
  ####### Sampling Prob  ---------------------------------------
  
  print("Sampling prob")
  
  sampling_prob <- raster::raster("F:/Connectivity/data/02_SDMs/Sampling_probability_map/human_footprint/Human_footprint_distance_model/output/corrected/HF_distance_model_corrected.tif")
  
  # Mask to MCP
  spMask <- mask(sampling_prob, bgExt_Cb) %>% crop(extent(bgExt_Cb))
  
  # Replace NAs with zeros
  spMask[is.na(spMask[])] <- 0 
  
  # Sample points depending on survey probability 
  bgSample_Cb_sp <- xyFromCell(spMask, sample(ncell(spMask), 10000, prob = values(spMask)))
  
  # Remove points that fall in the same cell
  bgSample_Cb_sp <- gridSample(bgSample_Cb_sp, envs_Cb$X30s_bio1_LAC, n=1) 
  
  # Turn to dataframe 
  bgSample_Cb_sp <- as.data.frame(bgSample_Cb_sp) %>%
    rename(longitude = x, latitude = y)
  
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/bgpoints_SamplingProb_Background.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  # Turn into points and export
  sf::st_as_sf(bgSample_Cb_sp, coords = c("longitude", "latitude"), crs = 4326) %>% 
    sf::st_write(complete_filename)
  
  ####### Sampling Prob - TGB number of points  ---------------------------------------
  
  print("Sampling prob - TGB points")
  
  # Calculate number of sampling points within MCP
  # Load Aves sampling sites
  Sampling_sites<- read.csv("F:/Connectivity/data/02_SDMs/Sampling_probability_map/Aves_occurrences/eBird/compiled/Total_sampling_sites.csv")
  
  # Turn to spatial points dataframe 
  sp::coordinates(Sampling_sites) <- ~longitude+latitude
  sp::proj4string(Sampling_sites) <- crs(bgExt_Cb)
  
  # Extract points that fall within calibration area
  print("intersection")
  Sampling_sites <- gIntersection(bgExt_Cb, Sampling_sites) 
  
  # Remove points that fall in the same cell
  TGB <- gridSample(Sampling_sites, envs_Cb$X30s_bio1_LAC, n=1) 
  
  # Sample points depending on survey probability 
  bgSample_Cb_sp_tgb <- xyFromCell(spMask, sample(ncell(spMask), nrow(TGB), prob = values(spMask)))
  
  # Remove points that fall in the same cell
  bgSample_Cb_sp_tgb <- gridSample(bgSample_Cb_sp_tgb, envs_Cb$X30s_bio1_LAC, n=1) 
  
  # Turn to dataframe 
  bgSample_Cb_sp_tgb <- as.data.frame(bgSample_Cb_sp_tgb) %>%
    rename(longitude = x, latitude = y)
  
  # Turn into points and export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/bgpoints_SamplingProbNumberOfSurveyedSites_Background.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  sf::st_as_sf(bgSample_Cb_sp_tgb, coords = c("longitude", "latitude"), crs = 4326) %>% 
    sf::st_write(complete_filename)
  
  ####### TGB ---------------------------------------
  print("TGB")
  
  # Turn to dataframe 
  bgSample_Cb_TGB <- as.data.frame(TGB) %>%
    rename(longitude = x, latitude = y) %>%
    dplyr::select(longitude, latitude)
  
  # EXPORT
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/bgpoints_TGB_Background.shp"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  sf::st_as_sf(bgSample_Cb_TGB, coords = c("longitude", "latitude"), crs = 4326) %>% 
    sf::st_write(complete_filename)
  
  
  ####### Model ---------------------------------------
  # Extract values of environmental layers for each background point
  print("extracting environmental values")
  bgEnvsVals_Cb_rb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_rb))
  bgEnvsVals_Cb_sp <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_sp))
  bgEnvsVals_Cb_sp_tgb <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_sp_tgb))
  bgEnvsVals_Cb_TGB <- as.data.frame(raster::extract(bgMask_Cb,  bgSample_Cb_TGB))
  
  # Remove NA values 
  bgSample_Cb_rb <- bgSample_Cb_rb[!(rowSums(is.na(bgEnvsVals_Cb_rb)) > 1), ]
  bgEnvsVals_Cb_rb <- na.omit(bgEnvsVals_Cb_rb)
  
  bgSample_Cb_sp <- bgSample_Cb_sp[!(rowSums(is.na(bgEnvsVals_Cb_sp)) > 1), ]
  bgEnvsVals_Cb_sp <- na.omit(bgEnvsVals_Cb_sp)
  
  bgSample_Cb_sp_tgb <- bgSample_Cb_sp_tgb[!(rowSums(is.na(bgEnvsVals_Cb_sp_tgb)) > 1), ]
  bgEnvsVals_Cb_sp_tgb <- na.omit(bgEnvsVals_Cb_sp_tgb)
  
  bgSample_Cb_TGB <- bgSample_Cb_TGB[!(rowSums(is.na(bgEnvsVals_Cb_TGB)) > 1), ]
  bgEnvsVals_Cb_TGB <- na.omit(bgEnvsVals_Cb_TGB)
  
  ##Add extracted values to background points table
  bgEnvsVals_Cb_rb <- cbind(bgSample_Cb_rb,  bgEnvsVals_Cb_rb)
  bgEnvsVals_Cb_sp <- cbind(bgSample_Cb_sp,  bgEnvsVals_Cb_sp)
  bgEnvsVals_Cb_sp_tgb <- cbind(bgSample_Cb_sp_tgb,  bgEnvsVals_Cb_sp_tgb)
  bgEnvsVals_Cb_TGB <- cbind(bgSample_Cb_TGB,  bgEnvsVals_Cb_TGB)
  
  # R code to get partitioned data
  groups_Cb_rb <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgSample_Cb_rb, 
    method = "cb1",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  groups_Cb_sp <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgSample_Cb_sp, 
    method = "cb1",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  groups_Cb_sp_tgb <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgSample_Cb_sp_tgb, 
    method = "cb1",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  groups_Cb_TGB <- part_partitionOccs(
    occs = occs_Cb ,
    bg =  bgEnvsVals_Cb_TGB, 
    method = "cb1",
    bgMask = bgMask_Cb,
    aggFact = 2)
  
  # Run maxent model for the selected species
  print("first model")
  model_Cb_rb <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_rb,
    user.grp = groups_Cb_rb, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ', 'H', 'LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_rb <- eval.results(model_Cb_rb)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/evalTable_Random_Background.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_rb, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_rb <- min(res_rb$or.10p.avg)
  max_auc_rb <- max(res_rb$auc.val.avg[which(res_rb$or.10p.avg == min_or_rb)])
  best_model_rb <- as.character(res_rb$tune.args[res_rb$auc.val.avg == max_auc_rb & res_rb$or.10p.avg == min_or_rb])
  
  # Select model and obtain raster prediction
  m_Cb_rb <- model_Cb_rb@models[[best_model_rb]]
  predSel_Cb_rb <- predictMaxnet(m_Cb_rb, bgMask_Cb,
                                 type = "cloglog", 
                                 clamp = TRUE)
  
  # Load projection area 
  proj_userExt_Cb <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/Colombia.shp")
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_rb <- proj_area(
    evalOut = model_Cb_rb,
    curModel = best_model_rb,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_rb <- raster::extract(predSel_Cb_rb, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_rb <- ceiling(length(occPredVals_rb) * 0.9)
  thresh_rb <- rev(sort(occPredVals_rb))[p10_rb]
  
  # Apply threshold to prediction
  p10_prediction_rb <- proj_area_Cb_rb$projArea
  p10_prediction_rb[Which(p10_prediction_rb<thresh_rb)]<-NA
  
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/projection_p10_Random_Background.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_rb, complete_filename)
  
  print("second model")
  model_Cb_sp <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_sp,
    user.grp = groups_Cb_sp, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ','H',  'LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_sp <- eval.results(model_Cb_sp)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/evalTable_SamplingProb_Background.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_sp, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_sp <- min(res_sp$or.10p.avg)
  max_auc_sp <- max(res_sp$auc.val.avg[which(res_sp$or.10p.avg == min_or_sp)])
  best_model_sp <- as.character(res_sp$tune.args[res_sp$auc.val.avg == max_auc_sp & res_sp$or.10p.avg == min_or_sp])
  
  # Select model and obtain raster prediction
  m_Cb_sp <- model_Cb_sp@models[[best_model_sp]]
  predSel_Cb_sp <- predictMaxnet(m_Cb_sp, bgMask_Cb,
                                 type = "cloglog", 
                                 clamp = TRUE)
  
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_sp <- proj_area(
    evalOut = model_Cb_sp,
    curModel = best_model_sp,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_sp <- raster::extract(predSel_Cb_sp, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_sp <- ceiling(length(occPredVals_sp) * 0.9)
  thresh_sp <- rev(sort(occPredVals_sp))[p10_sp]
  
  # Apply threshold to prediction
  p10_prediction_sp <- proj_area_Cb_sp$projArea
  p10_prediction_sp[Which(p10_prediction_sp<thresh_sp)]<-NA
  
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/projection_p10_SamplingProb_Background.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_sp, complete_filename)
  
  print("third model")
  model_Cb_sp_tgb <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_sp_tgb,
    user.grp = groups_Cb_sp_tgb, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ',  'H','LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_sp_tgb <- eval.results(model_Cb_sp_tgb)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/evalTable_SamplingProbNumberOfSurveyedSites_Background.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_sp_tgb, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_sp_tgb <- min(res_sp_tgb$or.10p.avg)
  max_auc_sp_tgb <- max(res_sp_tgb$auc.val.avg[which(res_sp_tgb$or.10p.avg == min_or_sp_tgb)])
  best_model_sp_tgb <- as.character(res_sp_tgb$tune.args[res_sp_tgb$auc.val.avg == max_auc_sp_tgb & res_sp_tgb$or.10p.avg == min_or_sp_tgb])
  
  # Select model and obtain raster prediction
  m_Cb_sp_tgb <- model_Cb_sp_tgb@models[[best_model_sp_tgb]]
  predSel_Cb_sp_tgb <- predictMaxnet(m_Cb_sp_tgb, bgMask_Cb,
                                     type = "cloglog", 
                                     clamp = TRUE)
  
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_sp_tgb <- proj_area(
    evalOut = model_Cb_sp_tgb,
    curModel = best_model_sp_tgb,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_sp_tgb <- raster::extract(predSel_Cb_sp_tgb, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_sp_tgb <- ceiling(length(occPredVals_sp_tgb) * 0.9)
  thresh_sp_tgb <- rev(sort(occPredVals_sp_tgb))[p10_sp_tgb]
  
  # Apply threshold to prediction
  p10_prediction_sp_tgb <- proj_area_Cb_sp_tgb$projArea
  p10_prediction_sp_tgb[Which(p10_prediction_sp_tgb<thresh_sp_tgb)]<-NA
  
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/projection_p10_SamplingProbNumberOfSurveyedSites_Background.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_sp_tgb, complete_filename)
  
  
  print("fourth model")
  model_Cb_TGB <- model_maxent(
    occs = occs_Cb,
    bg = bgEnvsVals_Cb_TGB,
    user.grp = groups_Cb_TGB, 
    bgMsk = bgMask_Cb,
    rms = c(0.5, 4), 
    rmsStep =  0.5,
    fcs = c('L', 'LQ',  'H','LQH'),
    clampSel = TRUE,
    algMaxent = "maxnet",
    parallel = TRUE,
    numCores = 6)
  
  # Overall results
  res_TGB <- eval.results(model_Cb_TGB)
  
  # export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/evalTable_TGB_Background.csv"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  
  write.csv(res_TGB, complete_filename, row.names = FALSE)
  
  # Selection based on the lowest average test omission rate (10th percentile), 
  # and to break ties, with the highest average validation AUC = fc.L_rm.0.5
  
  min_or_TGB <- min(res_TGB$or.10p.avg)
  max_auc_TGB <- max(res_TGB$auc.val.avg[which(res_TGB$or.10p.avg == min_or_TGB)])
  best_model_TGB <- as.character(res_TGB$tune.args[res_TGB$auc.val.avg == max_auc_TGB & res_TGB$or.10p.avg == min_or_TGB])
  
  # Select model and obtain raster prediction
  m_Cb_TGB <- model_Cb_TGB@models[[best_model_TGB]]
  predSel_Cb_TGB <- predictMaxnet(m_Cb_TGB, bgMask_Cb,
                                  type = "cloglog", 
                                  clamp = TRUE)
  
  
  # Generate a projection of the model to the desired area
  proj_area_Cb_TGB <- proj_area(
    evalOut = model_Cb_TGB,
    curModel = best_model_TGB,
    envs = envs_Cb , 
    outputType = "cloglog",
    alg = "maxnet",
    clamp = TRUE,
    pjExt = proj_userExt_Cb) 
  
  # Extract SDM predictions at occurrence points 
  occPredVals_TGB <- raster::extract(predSel_Cb_TGB, occs_Cb[ , c("longitude", "latitude")])
  
  # Calculate 10th percentile 
  p10_TGB <- ceiling(length(occPredVals_TGB) * 0.9)
  thresh_TGB <- rev(sort(occPredVals_TGB))[p10_TGB]
  
  # Apply threshold to prediction
  p10_prediction_TGB <- proj_area_Cb_TGB$projArea
  p10_prediction_TGB[Which(p10_prediction_TGB<thresh_TGB)]<-NA
  
  # Export
  filename1 <- "F:/Connectivity/outputs/02_SDMs/Models_samplingBias/BIAS/SAMPLING_PROB_MODELS/AGOL/"
  filename2 <- "Ognorhynchus icterotis"
  filename3 <- "/correction/projection_p10_TGB_Background.tif"
  complete_filename <- paste0(filename1, filename2)
  complete_filename <- paste0(complete_filename, filename3)
  
  writeRaster(p10_prediction_TGB, complete_filename)
  
  gc()

