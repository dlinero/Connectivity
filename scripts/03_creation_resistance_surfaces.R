library(sf)
library(tidyverse)
library(plyr)
library(ggplot2)
library(ggridges)
library(raster)
library(ncdf4)

# Increase memory to allocate results
memory.limit(999999999999)

# Progress bar
rasterOptions(progress = 'text',timer=TRUE)

# Create density plots with expert answers ----------------------------------------------------

theme_set(theme_bw(16))

data <- st_read("F:/Connectivity/outputs/03_connectivity/Resistance surface/Expert_resistance_surface___Connectivity.shp")

data <- as.data.frame(data)

data <- data %>% filter(Nombre != "BORRAR")  %>% droplevels() %>% select(-c("Creator", "EditDate", "Editor", "geometry"))

dataset <- data %>% 
  pivot_longer(!c("globalid", "Fecha", "Nombre", "Institucio", "Cargo", "Species", "CreationDa"), names_to = "Land_use", values_to = "Resistance")

write.csv(dataset, "F:/Connectivity/outputs/03_connectivity/Resistance surface/Pivot_table_fiz_acronyms.csv", row.names = FALSE)

rm(list=ls())

dataset <- read.csv("F:/Connectivity/outputs/03_connectivity/Resistance surface/Pivot_table_fiz_acronyms.csv")

# Hay una especie que el acr?nimo es PM y otra es PME, voy a cambiar el PM por ZZ para que el siguiente c?digo funciones

dataset <- dataset %>% mutate(selected_species = case_when(str_detect(Land_use, "AM") == TRUE ~ "Ara_macao", 
                                                           str_detect(Land_use, "BA") == TRUE ~ "Buteo_albigula",
                                                           str_detect(Land_use, "CP") == TRUE ~ "Campephilus_pollens", 
                                                           str_detect(Land_use, "CB") == TRUE ~ "Crypturellus_berlepschi", 
                                                           str_detect(Land_use, "CE") == TRUE ~ "Crypturellus_erythropus", 
                                                           str_detect(Land_use, "CS") == TRUE ~ "Crypturellus_soui", 
                                                           str_detect(Land_use, "GB") == TRUE ~ "Grallaria_bangsi", 
                                                           str_detect(Land_use, "GF") == TRUE ~ "Grallaria_flavotincta", 
                                                           str_detect(Land_use, "GH") == TRUE ~ "Grallaria_hypoleuca", 
                                                           str_detect(Land_use, "GN") == TRUE ~ "Grallaria_nuchalis", 
                                                           str_detect(Land_use, "HL") == TRUE ~ "Henicorhina_leucophrys", 
                                                           str_detect(Land_use, "MS") == TRUE ~ "Micrastur_semitorquatus", 
                                                           str_detect(Land_use, "MT") == TRUE ~ "Mitu_tomentosum", 
                                                           str_detect(Land_use, "MF") == TRUE ~ "Myioborus_flavivertex", 
                                                           str_detect(Land_use, "OH") == TRUE ~ "Odontophorus_hyperythrus", 
                                                           str_detect(Land_use, "OI") == TRUE ~ "Ognorhynchus_icterotis", 
                                                           str_detect(Land_use, "PG") == TRUE ~ "Patagioenas_goodsoni", 
                                                           str_detect(Land_use, "ZZ") == TRUE ~ "Phaethornis_malaris", 
                                                           str_detect(Land_use, "PC") == TRUE ~ "Picumnus_cinnamomeus", 
                                                           str_detect(Land_use, "PS") == TRUE ~ "Picumnus_squamulatus", 
                                                           str_detect(Land_use, "PA") == TRUE ~ "Pseudastur_albicollis", 
                                                           str_detect(Land_use, "PME") == TRUE ~ "Pyrrhura_melanura", 
                                                           str_detect(Land_use, "RB") == TRUE ~ "Ramphastos_brevis", 
                                                           str_detect(Land_use, "SI") == TRUE ~ "Spizaetus_isidori", 
                                                           str_detect(Land_use, "SO") == TRUE ~ "Spizaetus_ornatus", 
                                                           str_detect(Land_use, "TJ") == TRUE ~ "Tangara_johannae", 
                                                           TRUE ~ "REVISAR")) %>%
     mutate(land_category = case_when(str_detect(Land_use, "Agri") == TRUE ~ "Agriculture",
                                      str_detect(Land_use, "Gss") == TRUE ~ "Grassland and Shrubland",
                                      str_detect(Land_use, "For") == TRUE ~ "Forest",
                                      str_detect(Land_use, "Spar") == TRUE ~ "Sparse vegetation",
                                      str_detect(Land_use, "Urb") == TRUE ~ "Artificial surface or urban area",
                                      str_detect(Land_use, "Bar") == TRUE ~ "Bare area",
                                      str_detect(Land_use, "Swa") == TRUE ~ "Swampy or Often Flooded Vegetation",
                                      str_detect(Land_use, "Wat") == TRUE ~ "Surface Water", 
                                      TRUE ~ "REVISAR"))

# A mano voy a eliminar las especies que no lleno cada persona, por que ArcGIS lo llen? con ceros
# autom?ticamente 

write.csv(dataset, "F:/Connectivity/outputs/03_connectivity/Resistance surface/Only_filled_species.csv", row.names = FALSE)

rm(list=ls())

dataset <- read.csv("F:/Connectivity/outputs/03_connectivity/Resistance surface/Only_filled_species.csv")

# Crear el land use en espa?ol 
dataset <- dataset %>% mutate(land_use_espanol = case_when(land_category == "Agriculture" ~ "Mayormente agricultura", 
                                                           land_category == "Grassland and Shrubland" ~ "Pastizales y matorrales", 
                                                           land_category == "Forest" ~ "Bosque", 
                                                           land_category == "Sparse vegetation" ~ "Vegetaci?n dispersa",
                                                           land_category == "Artificial surface or urban area" ~ "Urbano",
                                                           land_category == "Bare area" ~ "?reas desnudas",
                                                           land_category == "Swampy or Often Flooded Vegetation" ~ "Vegetaci?n pantanosa o a menudo inundada",
                                                           land_category == "Surface Water" ~ "Cuerpos de Agua",
                                                           TRUE ~ "REVISAR"))


# Change order of land classes

dataset <- dataset %>% mutate(land_use_espanol = factor(land_use_espanol, levels = c("Cuerpos de Agua", "Vegetaci?n pantanosa o a menudo inundada", 
                                                          "?reas desnudas","Urbano","Vegetaci?n dispersa", "Pastizales y matorrales", "Mayormente agricultura", "Bosque" )))
levels(dataset$land_use_espanol)
# Plot 

first <- dataset %>% filter(selected_species %in% c("Ara_macao", "Buteo_albigula", "Campephilus_pollens", "Crypturellus_berlepschi")) %>% droplevels()


pdf(file = "F:/Connectivity/outputs/03_connectivity/Resistance surface/First_plot.pdf",   # The directory you want to save the file in
    width = 29.5, # The width of the plot in inches
    height = 11.7) # The height of the plot in inches


first %>%
  ggplot(aes(x=Resistance,y=land_use_espanol, fill=land_use_espanol)) +
  ggtitle("Distribuci?n de resistencia al movimiento de cada especie a trav?s de diferentes coberturas") +
  geom_density_ridges()+
  theme(legend.position = "none")+
  scale_fill_cyclical(values = c("navy", "paleturquoise4", "navajowhite3", "darkorange3", "olivedrab1", "darkolivegreen3", "brown3", "chartreuse4"))+
  facet_wrap(~selected_species, ncol = 2) + 
  scale_x_continuous(name="Resistencia al movimiento", breaks=c(1,50,100), label = c("Sin resistencia", "Resistencia moderada", "Resistencia absoluta"))+
  ylab("Cobetura del suelo")

dev.off()

second <- dataset %>% filter(selected_species %in% c("Crypturellus_erythropus", "Crypturellus_soui", "Grallaria_bangsi", "Grallaria_flavotincta")) %>% droplevels()


pdf(file = "F:/Connectivity/outputs/03_connectivity/Resistance surface/Second_plot.pdf",   # The directory you want to save the file in
    width = 25.3, # The width of the plot in inches
    height = 11.7) # The height of the plot in inches


second %>%
  ggplot(aes(x=Resistance,y=land_use_espanol, fill=land_use_espanol)) +
  ggtitle("Distribuci?n de resistencia al movimiento de cada especie a trav?s de diferentes coberturas") +
  geom_density_ridges()+
  theme(legend.position = "none")+
  scale_fill_cyclical(values = c("navy", "paleturquoise4", "navajowhite3", "darkorange3", "olivedrab1", "darkolivegreen3", "brown3", "chartreuse4"))+
  facet_wrap(~selected_species, ncol = 2) + 
  scale_x_continuous(name="Resistencia al movimiento", breaks=c(1,50,100), label = c("Sin resistencia", "Resistencia moderada", "Resistencia absoluta"))+
  ylab("Cobetura del suelo")

dev.off()


Third <- dataset %>% filter(selected_species %in% c("Grallaria_hypoleuca", "Grallaria_nuchalis", "Henicorhina_leucophrys", "Micrastur_semitorquatus")) %>% droplevels()


pdf(file = "F:/Connectivity/outputs/03_connectivity/Resistance surface/Third_plot.pdf",   # The directory you want to save the file in
    width = 29.5, # The width of the plot in inches
    height = 11.7) # The height of the plot in inches


Third %>%
  ggplot(aes(x=Resistance,y=land_use_espanol, fill=land_use_espanol)) +
  ggtitle("Distribuci?n de resistencia al movimiento de cada especie a trav?s de diferentes coberturas") +
  geom_density_ridges()+
  theme(legend.position = "none")+
  scale_fill_cyclical(values = c("navy", "paleturquoise4", "navajowhite3", "darkorange3", "olivedrab1", "darkolivegreen3", "brown3", "chartreuse4"))+
  facet_wrap(~selected_species, ncol = 2) + 
  scale_x_continuous(name="Resistencia al movimiento", breaks=c(1,50,100), label = c("Sin resistencia", "Resistencia moderada", "Resistencia absoluta"))+
  ylab("Cobetura del suelo")

dev.off()

Fourth <- dataset %>% filter(selected_species %in% c("Mitu_tomentosum", "Myioborus_flavivertex", "Odontophorus_hyperythrus", "Ognorhynchus_icterotis")) %>% droplevels()


pdf(file = "F:/Connectivity/outputs/03_connectivity/Resistance surface/Fourth_plot.pdf",   # The directory you want to save the file in
    width = 25.3, # The width of the plot in inches
    height = 11.7) # The height of the plot in inches


Fourth %>%
  ggplot(aes(x=Resistance,y=land_use_espanol, fill=land_use_espanol)) +
  ggtitle("Distribuci?n de resistencia al movimiento de cada especie a trav?s de diferentes coberturas") +
  geom_density_ridges()+
  theme(legend.position = "none")+
  scale_fill_cyclical(values = c("navy", "paleturquoise4", "navajowhite3", "darkorange3", "olivedrab1", "darkolivegreen3", "brown3", "chartreuse4"))+
  facet_wrap(~selected_species, ncol = 2) + 
  scale_x_continuous(name="Resistencia al movimiento", breaks=c(1,50,100), label = c("Sin resistencia", "Resistencia moderada", "Resistencia absoluta"))+
  ylab("Cobetura del suelo")

dev.off()

Fifth <- dataset %>% filter(selected_species %in% c("Picumnus_cinnamomeus", "Picumnus_squamulatus", "Pyrrhura_melanura", "Ramphastos_brevis")) %>% droplevels()


pdf(file = "F:/Connectivity/outputs/03_connectivity/Resistance surface/Fifth_plot.pdf",   # The directory you want to save the file in
    width = 25.3, # The width of the plot in inches
    height = 11.7) # The height of the plot in inches


Fifth %>%
  ggplot(aes(x=Resistance,y=land_use_espanol, fill=land_use_espanol)) +
  ggtitle("Distribuci?n de resistencia al movimiento de cada especie a trav?s de diferentes coberturas") +
  geom_density_ridges()+
  theme(legend.position = "none")+
  scale_fill_cyclical(values = c("navy", "paleturquoise4", "navajowhite3", "darkorange3", "olivedrab1", "darkolivegreen3", "brown3", "chartreuse4"))+
  facet_wrap(~selected_species, ncol = 2) + 
  scale_x_continuous(name="Resistencia al movimiento", breaks=c(1,50,100), label = c("Sin resistencia", "Resistencia moderada", "Resistencia absoluta"))+
  ylab("Cobetura del suelo")

dev.off()


six <- dataset %>% filter(selected_species %in% c("Spizaetus_isidori","Spizaetus_ornatus","Tangara_johannae")) %>% droplevels()


pdf(file = "F:/Connectivity/outputs/03_connectivity/Resistance surface/Six_plot.pdf",   # The directory you want to save the file in
    width = 25.3, # The width of the plot in inches
    height = 11.7) # The height of the plot in inches


six %>%
  ggplot(aes(x=Resistance,y=land_use_espanol, fill=land_use_espanol)) +
  ggtitle("Distribuci?n de resistencia al movimiento de cada especie a trav?s de diferentes coberturas") +
  geom_density_ridges()+
  theme(legend.position = "none")+
  scale_fill_cyclical(values = c("navy", "paleturquoise4", "navajowhite3", "darkorange3", "olivedrab1", "darkolivegreen3", "brown3", "chartreuse4"))+
  facet_wrap(~selected_species, ncol = 2) + 
  scale_x_continuous(name="Resistencia al movimiento", breaks=c(1,50,100), label = c("Sin resistencia", "Resistencia moderada", "Resistencia absoluta"))+
  ylab("Cobetura del suelo")

dev.off()


# Create resistance surfaces for each species - ESA  ----------------------------------------------------

# Load data
nc_ESA2020 <- nc_open('./data/03_connectivity/Resistance surface/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.nc')

# Check dimension of coordinates
lon <- ncvar_get(nc_ESA2020, "lon")
lat <- ncvar_get(nc_ESA2020, "lat")

# Get LCC variable
d_2020 <- ncvar_get(nc_ESA2020, "lccs_class")

# Save data in a raster
r <- raster(t(d_2020), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Close nc connection
nc_close(nc_ESA2020) 

# Empty memory
rm(d_2020)
gc()

# Load Colombia polygon 
crop_extent <- rgdal::readOGR("./data/03_connectivity/Resistance surface/col-administrative-divisions-shapefiles/col_admbnda_adm0_mgn_itos_20200416.shp")

# Aply a ~100km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 1, dissolve = TRUE)

# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Project
Equidistant_cylindrical <- crs("+proj=eqc +lon_0=-72.2460938 +lat_ts=0 +datum=WGS84 +units=m +no_defs")

data_project <- projectRaster(data_crop, crs = Equidistant_cylindrical, method = "ngb")

# Check resolution
res(data_project)

# Export
writeRaster(data_project, "./data/03_connectivity/Resistance surface/ESA2020_Colombia_100km_buffer.tif")

# Read
ESA2020 <- raster("./data/03_connectivity/Resistance surface/ESA2020_Colombia_100km_buffer.tif")

unique(ESA2020)

# Reclassify to generalized land-cover
reclass_df <- c(10, 1,
                11, 1,
                12, 1, 
                20, 1, 
                30, 1, 
                40, 2, 
                100, 2, 
                110, 2, 
                120, 2, 
                122, 2, 
                130, 2, 
                50, 3, 
                60, 3, 
                61, 3, 
                62, 3,
                80, 3, 
                90, 3, 
                150, 4, 
                153, 4, 
                160, 5, 
                170, 5, 
                180, 5, 
                190, 6, 
                200, 7, 
                210, 8, 
                220, 9)

# Turn to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
ESA2020_generalized <-reclassify(ESA2020, rcl=reclass_df)

# Export
writeRaster(ESA2020_generalized, "./outputs/03_connectivity/Resistance surface/ESA2020_Colombia_100kmBuffer_generalized.tif")

### Reclassify values per species --------------------------------------------------------

ESA2020_generalized <- raster("./outputs/03_connectivity/Resistance surface/ESA2020_Colombia_100kmBuffer_generalized.tif")

expert_responses <- read.csv("./outputs/03_connectivity/Resistance surface/Only_filled_species.csv", stringsAsFactors = FALSE)

# Calculate number of responses per species
number_responses <- expert_responses %>% group_by(selected_species) %>% dplyr::summarise(number_responses = length(unique(Nombre)))

# Calculate average resistance value per land cover/ per species
expert_responses_average <- expert_responses %>%
  dplyr::group_by(selected_species, land_category) %>%
  dplyr::summarise(mean_resistance = mean(Resistance), sd_resistance = sd(Resistance), min_resistance = min(Resistance), max_resistance = max(Resistance))


# Export mean values
write.csv(expert_responses_average, "F:/Connectivity/outputs/03_connectivity/Resistance surface/Mean_resistance_per_species.csv", row.names = FALSE)


species <- unique(expert_responses_average$selected_species)

# List of study areas 
st_areas <- list.files("./outputs/03_connectivity/Species distributions with buffer", 
                                  pattern = ".shp$",
                                  full.names = TRUE, 
                                  recursive =  FALSE)

st_areas <- data.frame(st_areas) %>% mutate(species_name = word(str_replace_all(st_areas, "[^[:alnum:]]", " "), start = 10, end = 11, sep = fixed(" ")))


for (i in 1:length(species)){
  
  print(i)
  
  data <- expert_responses_average %>% filter(selected_species == species[i]) %>% droplevels()
  
  reclass_data <- c(1, round(data$mean_resistance[data$land_category == "Agriculture"]),
                  2, round(data$mean_resistance[data$land_category == "Grassland and Shrubland"]),
                  3, round(data$mean_resistance[data$land_category == "Forest"]),
                  4, round(data$mean_resistance[data$land_category == "Sparse vegetation"]),
                  5, round(data$mean_resistance[data$land_category == "Swampy or Often Flooded Vegetation"]),
                  6, round(data$mean_resistance[data$land_category == "Artificial surface or urban area"]),
                  7, round(data$mean_resistance[data$land_category == "Bare area"]),
                  8, round(data$mean_resistance[data$land_category == "Surface Water"]),
                  9, NA)
  
  # Turn to matrix
  reclass_data <- matrix(reclass_data,
                       ncol = 2,
                       byrow = TRUE)
  
  # Reclassify land covers with resistance value
  data_resistance <-reclassify(ESA2020_generalized, rcl=reclass_data)
  
  # Cut to species study area
  species_index <- species[i] %>% str_replace_all("_", " ")
  species_st_area <- st_read(st_areas$st_areas[st_areas$species_name == species_index])
  species_st_area <- st_transform(species_st_area, crs =  crs(ESA2020_generalized))
  
  resistance <- mask(data_resistance, species_st_area) %>% crop(extent(species_st_area))
  
  # Export
  path1 <- "./outputs/03_connectivity/Resistance surface/rasters/"
  path2 <- species_index
  path3 <- "_resistance_surface.tif"
  path4 <- paste0(path1, path2, path3)
  
  writeRaster(resistance, path4, datatype='INT4S', overwrite=TRUE)
  
}
  
  

# Create resistance surfaces for each species - Corine 2018  ----------------------------------------------------

# Load  corine data
corine <- sf::st_read("F:/Connectivity/data/03_connectivity/Resistance surface/CORINE_2018/COBERTURAS CORINE 2018/ECOSISTEMAS_062021.gdb", layer = "cobertura_tierra_clc_2018")

corine_dissolve <- sf::st_read("F:/Connectivity/data/03_connectivity/Resistance surface/CORINE_2018/CORINE_2018/Default.gdb", layer = "CORINE_2018_generalized_Diss")

# Project
Equidistant_cylindrical <- crs("+proj=eqc +lon_0=-72.2460938 +lat_ts=0 +datum=WGS84 +units=m +no_defs")

data_project <- st_transform(corine_dissolve, crs = Equidistant_cylindrical)

# Add field of new generalized land-cover

data_project <- data_project %>% mutate(Land_cover_generalized = case_when(nivel_1 == "1" ~ "6", 
                                                                           nivel_1 == "2" ~ "1", 
                                                                           nivel_2 == "31" ~ "3", 
                                                                           nivel_2 == "32" ~ "2",
                                                                           nivel_2 == "33" ~ "7", 
                                                                           nivel_1 == "4" ~ "5", 
                                                                           nivel_1 == "5" ~ "8", 
                                                                           TRUE ~ "CHECK"
)) %>%
  mutate(Land_cover_generalized = case_when(nivel_5 == "31112" ~ "5", 
                                            nivel_5 == "31122" ~ "5", 
                                            nivel_5 == "31212" ~ "5", 
                                            nivel_5 == "31222" ~ "5", 
                                            nivel_5 == "32112" ~ "5",
                                            TRUE ~ Land_cover_generalized)) %>%
  dplyr::select(leyenda, nivel_1, nivel_2, nivel_3, nivel_4, nivel_5, nivel_6, Land_cover_generalized) %>%
  mutate(Land_cover_generalized = as.integer(Land_cover_generalized))


# Turn into raster 
# Create a generic raster, set the extent to the same as corine map
r.raster <- raster() 
raster::extent(r.raster) <- raster::extent(data_project)
res(r.raster) <- 300 # set cell size to 30 meters

# Identify unsupported geomtry types (namely multisurface)
as.data.frame(table(st_geometry_type(data_project)))

# Cast to multipolygon to remove unsupported geometry types
data_project = st_cast(data_project, "MULTIPOLYGON")

# Make a raster of the polygon:
corine_raster <- rasterize(x = data_project, y = r.raster, field = "Land_cover_generalized", background = -9999)



### Reclassify values per species --------------------------------------------------------

# Load result (all the previous steps were conducted in ArcGIS for greater speed)
corine_generalized <- raster("./data/03_connectivity/Resistance surface/CORINE_generalized/Corine_2018_generalized_projected.tif")

# Remove zero values
corine_generalized[corine_generalized==0]<-NA

expert_responses <- read.csv("./outputs/03_connectivity/Resistance surface/Only_filled_species.csv", stringsAsFactors = FALSE)

# Calculate number of responses per species
number_responses <- expert_responses %>% group_by(selected_species) %>% dplyr::summarise(number_responses = length(unique(Nombre)))

# Calculate average resistance value per land cover/ per species
expert_responses_average <- expert_responses %>%
  dplyr::group_by(selected_species, land_category) %>%
  dplyr::summarise(mean_resistance = mean(Resistance), sd_resistance = sd(Resistance), min_resistance = min(Resistance), max_resistance = max(Resistance))

species <- unique(expert_responses_average$selected_species)

# List of study areas 
st_areas <- list.files("./outputs/03_connectivity/Species distributions with buffer", 
                       pattern = ".shp$",
                       full.names = TRUE, 
                       recursive =  FALSE)

st_areas <- data.frame(st_areas) %>% mutate(species_name = word(str_replace_all(st_areas, "[^[:alnum:]]", " "), start = 10, end = 11, sep = fixed(" ")))


for (i in 1:length(species)){
  
  print(i)
  
  data <- expert_responses_average %>% filter(selected_species == species[i]) %>% droplevels()
  
  reclass_data <- c(1, round(data$mean_resistance[data$land_category == "Agriculture"]),
                    2, round(data$mean_resistance[data$land_category == "Grassland and Shrubland"]),
                    3, round(data$mean_resistance[data$land_category == "Forest"]),
                    5, round(data$mean_resistance[data$land_category == "Swampy or Often Flooded Vegetation"]),
                    6, round(data$mean_resistance[data$land_category == "Artificial surface or urban area"]),
                    7, round(data$mean_resistance[data$land_category == "Bare area"]),
                    8, round(data$mean_resistance[data$land_category == "Surface Water"]))
  
  # Turn to matrix
  reclass_data <- matrix(reclass_data,
                         ncol = 2,
                         byrow = TRUE)
  
  # Reclassify land covers with resistance value
  data_resistance <-reclassify(corine_generalized, rcl=reclass_data)
  
  # Remove zero values
  data_resistance[data_resistance==0]<-NA
  
  # Cut to species study area
  species_index <- species[i] %>% str_replace_all("_", " ")
  species_st_area <- st_read(st_areas$st_areas[st_areas$species_name == species_index])
  species_st_area <- st_transform(species_st_area, crs =  crs(corine_generalized))
  
  resistance <- mask(data_resistance, species_st_area) %>% crop(extent(species_st_area))
  
  # Export
  path1 <- "./outputs/03_connectivity/Resistance surface/rasters_corine/"
  path2 <- species_index
  path3 <- "_resistance_surface_corine.tif"
  path4 <- paste0(path1, path2, path3)
  
  writeRaster(resistance, path4, datatype='INT4S', overwrite=TRUE)
  
}


