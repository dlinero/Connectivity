library(tidyverse)
library(wallace)
library(dismo)
library(rgeos)
library(ENMeval)
library(foreach)
library(doParallel)

# PRIMER INTENTO -----------------------------------------------------------------------------

# Picumnus cinnamomeus - 45 records 
data <- read.csv("F:/Connectivity/outputs/02_SDMs/Intentos/final/FINAL_THINNED_DATA.csv")

occs_Pc <- data %>% 
  
  # Select the species
  filter(name == "Picumnus cinnamomeus") %>% 
  
  # Remove 2021 data, because ESA maps are available until 2020
  filter(year_observation != 2021) %>%
  
  # droplevels
  droplevels() %>% 
  
  # rename variable
  rename(scientific_name= name )

# Load environmental variables - Bioclimatic variables  --------------------------------------------

## Specify the directory with the environmental variables
dir_bio_Pc <- "F:/Connectivity/data/02_SDMs/Environmental_variables"
bio_path <- file.path(dir_bio_Pc, c('Bio1.tif', 'Bio2.tif',
                                      'Bio3.tif', 'Bio4.tif',
                                      'Bio5.tif', 'Bio6.tif',
                                      'Bio7.tif', 'Bio8.tif',
                                      'Bio9.tif', 'Bio10.tif',
                                      'Bio11.tif', 'Bio12.tif',
                                      'Bio13.tif', 'Bio14.tif',
                                      'Bio15.tif', 'Bio16.tif',
                                      'Bio17.tif', 'Bio18.tif', 'Bio19.tif'))
# Create environmental object 
bio_Pc <- envs_userEnvs(
  rasPath = bio_path,
  rasName = c('Bio1.tif', 'Bio2.tif',
              'Bio3.tif', 'Bio4.tif',
              'Bio5.tif', 'Bio6.tif',
              'Bio7.tif', 'Bio8.tif',
              'Bio9.tif', 'Bio10.tif',
              'Bio11.tif', 'Bio12.tif',
              'Bio13.tif', 'Bio14.tif',
              'Bio15.tif', 'Bio16.tif',
              'Bio17.tif', 'Bio18.tif', 'Bio19.tif'),
  doBrick = TRUE)

occs_xy_Pc_bio <- occs_Pc[c('longitude', 'latitude')]
occs_vals_Pc_bio <- as.data.frame(raster::extract(bio_Pc, occs_xy_Pc_bio))

# Load environmental variables - Land use -------------------------------------------------------
UseCores <- detectCores() - 2

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)

length(unique(occs_Pc$year_observation))
data_years <- sort(unique(occs_Pc$year_observation)) 

land_cover <- foreach (i = 1:11) %dopar%{ 
  
  library(dplyr)
  library(raster)
  
  pattern_files <- paste0("ESA.+", data_years[i])
  # Load ESA layers 
  ESA <- list.files("F:/Connectivity/data/02_SDMs/Environmental_variables",full.names = T, pattern = pattern_files) 
  ESA <- raster::stack(ESA)
  new_names <- paste0("ESA", c(1:8)) 
  names(ESA)<- new_names
  ESA <- raster::brick(ESA)
  
  occs_xy_Pc <- occs_Pc %>% mutate(longitude = case_when(year_observation == data_years[i] ~ longitude, 
                                                              TRUE ~ NA_real_)) %>% dplyr::select(longitude, latitude)
  
  values <- as.data.frame(raster::extract(ESA, occs_xy_Pc))

return(values)
  
}

stopCluster(cl)

complete_LC <- land_cover[[1]]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[2]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[3]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[4]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[5]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[6]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[7]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[8]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[9]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[10]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]
complete_LC[which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ] <- land_cover[[11]][which(rowSums(is.na(complete_LC)) == ncol(complete_LC)), ]

# Join occurrences with environmental variables --------------------------------------------------

# Join all environmental variables
occs_vals_Pc <- cbind(occs_vals_Pc_bio, complete_LC)

# remove occurrence records with NA environmental values
occs_Pc <- occs_Pc[!(rowSums(is.na(occs_vals_Pc)) > 1), ]

# also remove variable value rows with NA environmental values
occs_vals_Pc <- na.omit(occs_vals_Pc)

# add columns for env variable values for each occurrence record
occs_Pc <- cbind(occs_Pc, occs_vals_Pc)

# Process environmental data ---------------------------------------------

# Load MCP with a ~10km buffer as calibration area and background extent 
bgExt_Pc <- rgdal::readOGR("F:/Connectivity/outputs/02_SDMs/Intentos/final", layer= "Picumnus cinnamomeus_mcp_10km_thin5km", verbose = TRUE)

# Create a raster brick with Bioclimatic variables and ESA 2015 (that is the one we are going to use for background sampling)
# Load layers 
Env_layers_bck <- list.files("F:/Connectivity/data/02_SDMs/Environmental_variables",full.names = T, pattern = "Bio|ESA.+2015") 
Env_layers_bck <- raster::stack(Env_layers_bck)
new_names <- paste0("ESA", c(1:8)) 
names(Env_layers_bck)[20:27]<- new_names
Env_layers_bck <- raster::brick(Env_layers_bck)


# Mask environmental data to provided extent
bgMask_Pc <- penvs_bgMask(
  occs = occs_Pc,
  envs = Env_layers_bck,
  bgExt = bgExt_Pc)


# Sample background points from the provided area
bgSample_Pc <- penvs_bgSample(
  occs = occs_Pc,
  bgMask =  bgMask_Pc,
  bgPtsNum = 10000)


# Extract values of environmental layers for each background point
bgEnvsVals_Pc <- as.data.frame(raster::extract(bgMask_Pc,  bgSample_Pc))


##Add extracted values to background points table
bgEnvsVals_Pc <- cbind(scientific_name = paste0("bg_", "Picumnus cinnamomeus"), bgSample_Pc,
                       bgEnvsVals_Pc)


# Partition occurrence data ---------------------------------------------------

# Partition occurrences and background points for model training and
# validation using “hierarchical checkerboard”, a spatial partition method
# with an aggregation factor of 2.


# R code to get partitioned data
groups_Pc <- part_partitionOccs(
  occs = occs_Pc ,
  bg =  bgSample_Pc, 
  method = "cb2",
  bgMask = bgMask_Pc,
  aggFact = 2) 

# Build and Evaluate Niche Model -------------------------------------------------------

# Generating a species distribution model using the maxnet algorithm as
# implemented in ENMeval V2.0 (with clamping = TRUE). For tuning using L,
# LQ, H, LQH feature classes and regularization multipliers in the 0.5, 4
# range increasing by 0.5. Not using any categorical predictor variables.

# Remove columns from the occs_Pc variables, so that it has the same dimensions 
# of the background values dataframe 
cols <- colnames(bgEnvsVals_Pc)
cols <- cols[2:30]

occs_Pc_model <- occs_Pc %>% 
  
  dplyr::select(all_of(cols))

bgEnvsVals_Pc_model <- bgEnvsVals_Pc %>%
  
  dplyr::select(all_of(cols))

# Primer intento 

tune.args = list(fc = c('L'), rm = 0.5)

e <- ENMeval::ENMevaluate(occs = occs_Pc_model, 
                          bg = bgEnvsVals_Pc_model, partitions = "user", 
                          user.grp = groups_Pc, tune.args = tune.args, 
                          doClamp = TRUE, algorithm = "maxnet",  
                          parallel = FALSE, 
                          updateProgress = FALSE, quiet = FALSE)



# Overall results
res <- eval.results(e)

# Select the model with delta AICc equal to 0, or the one with the lowest AICc score.
# In practice, models with delta AICc scores less than 2 are usually considered 
# statistically equivalent.
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC

opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.
mod.seq <- eval.models(e)[[opt.seq$tune.args]]

# Here are the non-zero coefficients in our model.
mod.seq$betas

# And these are the marginal response curves for the predictor variables wit non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

pdf(file= "F:/Connectivity/outputs/02_SDMs/Models_Bio_ESA/plot1.pdf")
plot(mod.seq, type = "cloglog")
dev.off()

# Select best model and obtain raster prediction ---------------------------------------


# Generate projection area: species distribution + 10km buffer
proj_userExt_Pc <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN", layer = "Picumnus_cinnamomeus_Ayerbe_IUCN_buffer", verbose = TRUE)
crs(proj_userExt_Pc) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Create raster brick with 2020 ESA data 
Env_layers_bck_2020 <- list.files("F:/Connectivity/data/02_SDMs/Environmental_variables",full.names = T, pattern = "Bio|ESA.+2020") 
Env_layers_bck_2020 <- raster::stack(Env_layers_bck_2020)
new_names <- paste0("ESA", c(1:8)) 
names(Env_layers_bck_2020)[20:27]<- new_names
Env_layers_bck_2020 <- raster::brick(Env_layers_bck_2020)


# Generate a projection of the model to the desired area
proj_area_Pc <- proj_area(
  evalOut = e,
  curModel = "fc.L_rm.0.5",
  envs = Env_layers_bck_2020 , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Pc) 

# Save results
Projected_raster <- proj_area_Pc$projArea

writeRaster(Projected_raster, "F:/Connectivity/outputs/02_SDMs/Models_Bio_ESA/raster2.tif")

# Plot in R
#store the cropped projection variables
projExt_Pc <- proj_area_Pc$projExt

###Make map of projection
bb_Pc <-  bgExt_Pc@bbox
bbZoom <- polyZoom(bb_Pc[1, 1], bb_Pc[2, 1], bb_Pc[1, 2], 
                   bb_Pc[2, 2], fraction = 0.05)
mapProjVals_Pc <- getRasterVals(proj_area_Pc$projArea,"cloglog")
rasCols_Pc <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Pc), mapProjVals_Pc, na.color = 'transparent')
rasPal_Pc <- colorNumeric(rasCols_Pc, mapProjVals_Pc, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) 
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  addLegend("bottomright", pal = legendPal,
            title = "Predicted Suitability<br>(Projected)",
            values = mapProjVals_Pc, layerId = 'proj',
            labFormat = reverseLabels(2, reverse_order = TRUE)) %>%
  # map model prediction raster and projection polygon
  clearMarkers() %>% clearShapes() %>% removeImage('projRas') %>%
  addRasterImage(proj_area_Pc$projArea, colors = rasPal_Pc, opacity = 0.7,
                 layerId = 'projRas', group = 'proj', method = "ngb") %>%
  ##add projection polygon (user provided area)
  addPolygons(data = proj_userExt_Pc, fill = FALSE,
              weight = 4, color = "blue", group = 'proj')


# Run maxent model for the selected species ------------------------------------------------------------

rms.interval <- seq(0.5, 4, 0.5)
tune.args = list(fc = c('L', 'LQ', 'H', 'LQH'), rm = rms.interval)

e <- ENMeval::ENMevaluate(occs = occs_Pc_model, 
                          bg = bgEnvsVals_Pc_model, partitions = "user", 
                          user.grp = groups_Pc, tune.args = tune.args, 
                          doClamp = TRUE, algorithm = "maxnet",  
                          parallel = TRUE, numCores = 6, parallelType = "doSNOW", 
                          updateProgress = FALSE, quiet = FALSE)

# Overall results
res <- eval.results(e)

# Selection based on the lowest average test omission rate (10th percentile), 
# and to break ties, with the highest average validation AUC

opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq

# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.
mod.seq <- eval.models(e)[[opt.seq$tune.args]]

# Here are the non-zero coefficients in our model.
mod.seq$betas

# And these are the marginal response curves for the predictor variables wit non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).

pdf(file= "F:/Connectivity/outputs/02_SDMs/Models_Bio_ESA/Picumnus_cinnamomeus/response_curves.pdf")
plot(mod.seq, type = "cloglog")
dev.off()

# Select best model and obtain raster prediction ---------------------------------------


# Generate projection area: species distribution + 10km buffer
proj_userExt_Pc <- rgdal::readOGR("F:/Connectivity/data/species_distributions/Ayerbe_and_IUCN", layer = "Picumnus_cinnamomeus_Ayerbe_IUCN_buffer", verbose = TRUE)
crs(proj_userExt_Pc) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Create raster brick with 2020 ESA data 
Env_layers_bck_2020 <- list.files("F:/Connectivity/data/02_SDMs/Environmental_variables",full.names = T, pattern = "Bio|ESA.+2020") 
Env_layers_bck_2020 <- raster::stack(Env_layers_bck_2020)
new_names <- paste0("ESA", c(1:8)) 
names(Env_layers_bck_2020)[20:27]<- new_names
Env_layers_bck_2020 <- raster::brick(Env_layers_bck_2020)

# Generate a projection of the model to the desired area
proj_area_Pc <- proj_area(
  evalOut = e,
  curModel = opt.seq$tune.args,
  envs = Env_layers_bck_2020 , 
  outputType = "cloglog",
  alg = "maxnet",
  clamp = TRUE,
  pjExt = proj_userExt_Pc) 

# Save results
Projected_raster <- proj_area_Pc$projArea

writeRaster(Projected_raster, "F:/Connectivity/outputs/02_SDMs/Models_Bio_ESA/Picumnus_cinnamomeus/HS_Picumnus_cinnamomeus.tif")



