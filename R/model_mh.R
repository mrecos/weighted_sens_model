## packages
library("raster")
library('tidyverse')
library("sf")
library("sp")
library("mapview")
library("FedData")
library("elevatr")
library("doParallel")    # Parallel
library("foreach")    # Parallel
library("velox")


## Data locations
data_loc <- "U:/z_OldServer/Projects/NY/Oakdale_Fayette_model/GIS"
SHP_loc  <- file.path(data_loc,"SHP","Client_files")

output_res_ft = 50     # US_survey feet
use_parallel  = FALSE

if(isTRUE(use_parallel)){
  # http://www.gis-blog.com/increasing-the-speed-of-raster-processing-with-r-part-23-parallelisation/
  #Define how many cores you want to use
  UseCores <- detectCores() -1
  cl       <- makeCluster(UseCores)
  registerDoParallel(cl)
}

if(!dir.exists("RASTERS")){
  dir.create("RASTERS")
}
if(!dir.exists("GRD")){
  dir.create("GRD")
}

## projections details
# unprojected WGS84
WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# New York Central SP NAD 83
SPNYCentral83 <- "+proj=tmerc +lat_0=40 +lon_0=-76.58333333333333 +k=0.9999375 +x_0=250000 +y_0=0 +datum=NAD83 +units=us-ft +no_defs"

## Read SHPs
sitesSP    <- read_sf(file.path(SHP_loc, "ArchaeoSitesPt.shp"))
sitesSP_sp <- as(sitesSP, "Spatial")
clSP       <- read_sf(file.path(SHP_loc, "ProposedFinalRouteCL_10-05-2017.shp"))
clSP_sp    <- as(clSP, "Spatial")
temp_bboxSP<- read_sf(file.path(data_loc,"SHP","Temp","temp_bbox.shp"))


start.time <- Sys.time()

### buffering
cl_5mi_buffSP  <- st_buffer(st_combine(clSP), 26400) %>%
  st_sf(.) %>%
  mutate(value = 1)
# cl_5mi_buffSP  <- st_intersection(cl_5mi_buffSP, temp_bboxSP)  #### <- CROPPED FOR TESTING!!!!!!!
# convert to WGS84 for raster work & cast to sp object
cl_5mi_buffWGS <- st_transform(cl_5mi_buffSP, WGS84)
cl_5mi_buffWGS_sp <- as(cl_5mi_buffWGS, "Spatial")

### buffering bigger!
cl_525mi_buffSP  <- st_buffer(st_combine(clSP), 27720) %>%
  st_sf(.) %>%
  mutate(value = 1)
# cl_525mi_buffSP    <- st_intersection(cl_525mi_buffSP, temp_bboxSP)  #### <- CROPPED FOR TESTING!!!!!!!
cl_525mi_buffWGS <- st_transform(cl_525mi_buffSP, WGS84)
cl_525mi_buffWGS_sp <- as(cl_525mi_buffWGS, "Spatial")

#### RASTERS!!!! FTW!
### get Elevation data & crop
# elevation_full <- get_elev_raster(cl_5mi_buffWGS_sp, z = 12, src = "aws")
# elevation      <- mask(elevation_full, cl_5mi_buffWGS_sp)
# elevationSP    <- projectRaster(elevation, crs = CRS(SPNYCentral83),
#                              res = output_res_ft)
# writeRaster(elevationSP,file.path("RASTERS","model_elevationSP.tiff"), overwrite=TRUE)
# writeRaster(elevationSP,file.path("GRD","model_elevationSP.grd"),overwrite=TRUE)
elevationSP <- raster(file.path("GRD","model_elevationSP.grd"))

####### MAKE VELOX TEMPLATE!!!!!!
template <- elevationSP      # <- pre masked!
template[template >= 0] <- 0 # <- zero out values
temp_vx <- velox(template)   # <- cast to velox class
###############

# slope
# slope   <- raster::terrain(elevation, out = "slope", unit='radians')
# slope   <- tan(slope)*100
# slopeSP <- projectRaster(slope, elevationSP)
# writeRaster(slopeSP,file.path("RASTERS","model_slope_degree.tiff"), overwrite=TRUE)
# writeRaster(slopeSP,file.path("GRD","model_slope_degree.GRD"), overwrite=TRUE)
slopeSP <- raster(file.path("GRD","model_slope_degree.GRD"))

### NED 
NHD <- FedData::get_nhd(template=cl_525mi_buffWGS_sp, label='Oakdale_Fayette_model') ### <- change for full run
# flowlines
flowlines     <- NHD[[grep("Flowline",names(NHD))]]
flowlinesSP   <- spTransform(flowlines, SPNYCentral83)
flowlinesSP$raster <- 1
flowline_vx   <- temp_vx$copy() # make copy of template
flowline_vx$rasterize(flowlinesSP, field="raster", band=1) # rasterize from "raster" field
flowlinesSP_r <- flowline_vx$as.RasterLayer(band = 1) # cast back to raster object
writeRaster(flowlinesSP_r,file.path("RASTERS","flowlinesSP.tiff"), overwrite=TRUE)

## AREA (large streams)
areas     <-  NHD[[grep("Area",names(NHD))]]
areasSP   <- spTransform(areas, SPNYCentral83)
areasSP$raster <- 1
areas_vx  <- temp_vx$copy() # make copy of template
areas_vx$rasterize(areasSP, field="raster", band=1) # rasterize from "raster" field
areasSP_r <- areas_vx$as.RasterLayer(band = 1) # cast back to raster object
writeRaster(areasSP_r,file.path("RASTERS","areasSP.tiff"), overwrite=TRUE)

## Waterbodies
waterbody     <-  NHD[[grep("Waterbody",names(NHD))]]
waterbodySP   <- spTransform(waterbody, SPNYCentral83)
waterbodySP$raster <- 1
waterbody_vx  <- temp_vx$copy() # make copy of template
waterbody_vx$rasterize(waterbodySP, field="raster", band=1) # rasterize from "raster" field
waterbodySP_r <- waterbody_vx$as.RasterLayer(band = 1) # cast back to raster object
writeRaster(waterbodySP_r,file.path("RASTERS","waterbodySP.tiff"), overwrite=TRUE)

### max of digitized rasters then distance
h20_stack <- stack(flowlinesSP_r, areasSP_r, waterbodySP_r )
names(h20_stack) <- c("flowlines","areas","waterbodies")
min_h20SP_r <- max(h20_stack)
min_h20SP_r[min_h20SP_r != 1] <- NA
writeRaster(min_h20SP_r,file.path("RASTERS","min_h20SP_r.tiff"), overwrite=TRUE)

# min_h20SP_dist <- distance(min_h20SP_r, progress = "text")
min_h20SP_dist <- raster(file.path("RASTERS","model_min_h20SP_dist_esri.tif")) ### <- ESRI!!!
min_h20SP_dist <- mask(min_h20SP_dist, elevationSP, progress = "text")
writeRaster(min_h20SP_dist,file.path("RASTERS","model_min_h20SP_dist.tiff"), overwrite=TRUE)
writeRaster(min_h20SP_dist,file.path("GRD","model_min_h20SP_dist.grd"), overwrite=TRUE)

gc()
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

