# Load libraries 
library(raster)
library(ncdf4)
library(foreach)
library(doParallel)
library(rgdal)


# Set working directory
setwd("E:/VM/VM/Connectivity/ESA/")

# 2006 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2006.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nlon <- dim(lon)
nlat <- dim(lat)
print(c(nlon, nlat))

# Get LCC variable
d_2006 <- ncvar_get(nc, "lccs_class")
# Check dimensions
dim(d_2006)

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2006), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Save to GeoTiff file 
writeRaster(data_crop, "ESA_data_2006_crop.tif", "GTiff", overwrite=TRUE)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
              121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
              201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                    ncol = 2,
                    byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

landuse <- raster("ESA_data_2006_reclassify.tif")

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)

# Aggregate raster to calculate fractional cover of each land use
unique(landuse)

c_1 <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
  sum(vals==1, na.rm=na.rm)/length(vals)
})


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 4

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command
foreach(i=2:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2006_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  

# http://zevross.com/blog/2015/03/30/map-and-analyze-raster-data-in-r/
# https://gis.stackexchange.com/questions/297908/calculating-block-statistics-from-land-cover-raster-using-r
# https://gis.stackexchange.com/questions/410802/calculate-percentage-of-a-class-in-a-raster-within-a-buffer-in-moving-windows-in/410811#410811

plot(crop_extent_buffer)
plot(crop_extent, add = TRUE)  


# 2010 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2010.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nlon <- dim(lon)
nlat <- dim(lat)
print(c(nlon, nlat))

# Get LCC variable
d_2010 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2010), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2010, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2010_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  

# 2011 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2011.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2011 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2011), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2011, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2011_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  


# 2012 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2012.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2012 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2012), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2012, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2012_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  

# 2013 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2013.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2013 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2013), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2013, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2013_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  


# 2014 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2014.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2014 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2014), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2014, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2014_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  


# 2015 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2015.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2015 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2015), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2015, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2015_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  

# 2016 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2016.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2016 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2016), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2016, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command ------------------------------------------------------
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2016_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(data_crop, "Crop_2016.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  


# 2017 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2017.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2017 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2017), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2017, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command ------------------------------------------------------
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2017_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(data_crop, "Crop_2016.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  


# 2018 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2018.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2018 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2018), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2018, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command ------------------------------------------------------
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2018_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(data_crop, "Crop_2016.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  


# 2019 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2019.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2019 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2019), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2019, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command ------------------------------------------------------
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2019_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(data_crop, "Crop_2016.tif", "GTiff", overwrite=TRUE)


# Snap rater to WC variables  


# 2020 data -----------------------------------------------------------------------
# Load data
nc = nc_open("ESA_2020.nc")
# Check dimension of coordinates
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")

# Get LCC variable
d_2020 <- ncvar_get(nc, "lccs_class")

# Close nc connection
nc_close(nc) 

# Save data in a raster
r <- raster(t(d_2020), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

# Load LAC polygon 
crop_extent <- rgdal::readOGR("ESA_data/Latin_America.shp")
# Aplly a ~50km buffer
crop_extent_buffer <- raster::buffer(crop_extent, width = 0.5, dissolve = TRUE)
# Clip raster
data_crop <- crop(r, crop_extent_buffer)

# Remove heavy variables
rm(d_2020, nc, r, lat, lon)

# Merge some land uses to decrease the number of classes
reclass_df <- c(10, 1, 11,
                1, 12, 1,
                20, 1, 30,
                1, 40, 1, 
                50, 2, 60, 2, 
                61, 2, 62, 2, 
                70, 2, 71, 2, 72, 
                2, 80, 2, 81, 2, 
                82, 2, 90, 2, 100, 2, 
                160, 2, 170, 2, 110, 3, 
                130, 3, 180, 4, 120, 5, 
                121, 5, 122, 5, 140, 6, 150, 6, 152, 6, 153, 6, 200, 6, 
                201, 6, 202, 6, 190, 7, 210, 8, 220, 9)

# Turn vector to matrix
reclass_df <- matrix(reclass_df,
                     ncol = 2,
                     byrow = TRUE)

# Reclassify
landuse<-reclassify(data_crop, rcl=reclass_df)

unique(landuse)

# Disaggregate raster to change resolution to ~100 m
landuse_disa <- raster::disaggregate(landuse, fact=3)


# Paralelized for other land uses 
#Define how many cores you want to use
UseCores <- detectCores() - 3

#Register CoreCluster
cl  <- makeCluster(UseCores)

registerDoParallel(cl)


#Use foreach loop and %dopar% command ------------------------------------------------------
foreach(i=1:9) %dopar% {
  library(raster)
  
  fraction <- aggregate(landuse_disa, fact=10, fun=function(vals, na.rm) {
    sum(vals==i, na.rm=na.rm)/length(vals)
  })
  
  outname <- paste0("Aggregated_maps/ESA_data_2020_reclassify_c", i, ".tif")
  
  writeRaster(fraction, 
              filename  = outname,
              "GTiff",
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)


# Export
writeRaster(landuse, "ESA_data_2006_reclassify.tif", "GTiff", overwrite=TRUE)

writeRaster(landuse_disa, "ESA_data_2006_reclassify_disa.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(c_1, "Aggregated_maps/ESA_data_2006_reclassify_c1.tif", "GTiff", overwrite=TRUE)

writeRaster(data_crop, "Crop_2016.tif", "GTiff", overwrite=TRUE)



### Align Worldclim variables -------------------------------------------------------

# Load ESA raster
ESA_2006_1 <-  raster::raster("F:/Connectivity/data/02_SDMs/ESA/Aggregated_maps/ESA_data_2006_reclassify_c1.tif")

# List of bio rasters 
bio_files <- list.files("F:/Connectivity/outputs/02_SDMs/wc",full.names = T, pattern = "30s.+.tif$") 

for(i in 1:length(bio_files)){
  
  print(i)

  # load bio variable
  temp <- raster::raster(bio_files[i])
  
  # Resample to ESA grid
  temp_r = raster::resample(temp, ESA_2006_1, "bilinear")
  
  # If the raster are aligned, export object
  if (raster::compareRaster(temp_r, ESA_2006_1) == TRUE){
    print("rasters aligned")
    outdir <- paste0("F:/Connectivity/data/02_SDMs/Environmental_variables/Bio", i, ".tif")
    writeRaster(temp_r, outdir, "GTiff", overwrite=TRUE)
    
  } else {
    
    # Otherwise stop the loop
    print("rasters are not aligned")
    break
  }
}
  
  