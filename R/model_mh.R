## packages
library("raster")
library('tidyverse')
library("sf")
library("sp")
library("mapview")
library("FedData")
library("elevatr")
library("snow")    # Parallel

## Data locations
data_loc <- "U:/z_OldServer/Projects/NY/Oakdale_Fayette_model/GIS"
SHP_loc  <- file.path(data_loc,"SHP","Client_files")

## ignored folders
if(!dir.exists("RASTERS")){
  dir.create("RASTERS")
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

m1 <- mapview(list(sitesSP, clSP, temp_bboxSP))


### buffering
cl_5mi_buffSP  <- st_buffer(st_combine(clSP), 26400) %>%
  st_sf(.) %>%
  mutate(value = 1)
cl_5mi_buffSP  <- st_intersection(cl_5mi_buffSP, temp_bboxSP)  #### <- CROPPED FOR TESTING!!!!!!!
# convert to WGS84 for raster work & cast to sp object
cl_5mi_buffWGS <- st_transform(cl_5mi_buffSP, WGS84)
cl_5mi_buffWGS_sp <- as(cl_5mi_buffWGS, "Spatial")

### buffering bigger!
cl_525mi_buffSP  <- st_buffer(st_combine(clSP), 27720) %>%
  st_sf(.) %>%
  mutate(value = 1)
cl_525mi_buffSP  <- st_intersection(cl_525mi_buffSP, temp_bboxSP)  #### <- CROPPED FOR TESTING!!!!!!!
cl_525mi_buffSP_sp <- as(cl_5mi_buffSP, "Spatial")

beginCluster(8) ###### <- PARALLEL n=cores!!

#### RASTERS!!!! FTW!
### get Elevation data & crop
elevation_full <- get_elev_raster(cl_5mi_buffWGS_sp, z = 12, src = "aws")
elevation <- mask(elevation_full, cl_5mi_buffWGS_sp)
elevationSP <- projectRaster(elevation, crs = CRS(SPNYCentral83), res = 50)
# writeRaster(elevationSP,file.path("RASTERS","elevationSP.tiff"))
m1 <- mapview(elevationSP)
m2 <- mapview(sitesSP)
# m1 + m2

# slope
slope <- raster::terrain(elevation, out = "slope", unit='radians')
slope <- tan(slope)*100
slopeSP <- projectRaster(slope, elevationSP)
# writeRaster(slopeSP,file.path("RASTERS","slop_degree.tiff"))
m3 <- mapview(slopeSP)
# m2 + m3

### NED 
NHD <- FedData::get_nhd(template=cl_525mi_buffWGS_sp, label='test')
# flowlines
flowlines <- NHD[["_Flowline"]]
flowlinesSP <- spTransform(flowlines, SPNYCentral83)
flowlinesSP_r <- rasterize(flowlinesSP, elevationSP)
flowlinesSP_r[flowlinesSP_r > 0] <- 1
# writeRaster(flowlinesSP_r,file.path("RASTERS","flowlines.tiff"), overwrite=TRUE)
m4 <- mapview(flowlinesSP_r)
flowlinesSP_dist <- raster::distance(flowlinesSP_r)
flowlinesSP_dist <- projectRaster(flowlinesSP_dist, elevationSP)
flowlinesSP_dist <- mask(flowlinesSP_dist, elevationSP)
# writeRaster(flowlinesSP_dist,file.path("RASTERS","fflowlinesSP_dist.tiff"), overwrite=TRUE)
m5 <- mapview(flowlinesSP_dist)

## AREA (large streams)
areas <- NHD[["_Area"]]
m6_a <- mapview(areas)
areasSP <- spTransform(areas, SPNYCentral83)
areasSP_r <- rasterize(areasSP, elevationSP)
areasSP_r[areasSP_r > 0] <- 1
# writeRaster(areasSP_r,file.path("RASTERS","areasSP.tiff"), overwrite=TRUE)
m6 <- mapview(areasSP_r)
areasSP_dist <- raster::distance(areasSP_r)
areasSP_dist <- projectRaster(areasSP_dist, elevationSP)
areasSP_dist <- mask(areasSP_dist, elevationSP)
# writeRaster(areasSP_dist,file.path("RASTERS","areasSP_dist.tiff"), overwrite=TRUE)
m7 <- mapview(areasSP_dist)

## Waterbodies
waterbody <- NHD[["_Waterbody"]]
m8_a <- mapview(waterbody)
waterbodySP <- spTransform(waterbody, SPNYCentral83)
waterbodySP_r <- rasterize(waterbodySP, elevationSP)
waterbodySP_r[waterbodySP_r > 0] <- 1
writeRaster(waterbodySP_r,file.path("RASTERS","waterbodySP.tiff"), overwrite=TRUE)
m8 <- mapview(waterbodySP_r)
waterbodySP_dist <- raster::distance(waterbodySP_r)
waterbodySP_dist <- projectRaster(waterbodySP_dist, elevationSP)
waterbodySP_dist <- mask(waterbodySP_dist, elevationSP)
# writeRaster(waterbodySP_dist,file.path("RASTERS","waterbodySP_dist.tiff"), overwrite=TRUE)
m9 <- mapview(waterbodySP_dist)

### take min of H20 dists for min to H20


