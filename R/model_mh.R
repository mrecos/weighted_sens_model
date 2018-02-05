## packages
library("raster")
library('tidyverse')
library("sf")
library("sp")
library("mapview")
library("FedData")
library("elevatr")

## Data locations
data_loc <- "U:/z_OldServer/Projects/NY/Oakdale_Fayette_model/GIS"
SHP_loc  <- file.path(data_loc,"SHP","Client_files")

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


### RASTERS!!!! FTW!
### get Elevation data & crop
elevation_full <- get_elev_raster(cl_5mi_buffWGS_sp, z = 12, src = "aws")
elevation <- mask(elevation_full, cl_5mi_buffWGS_sp)
elevationSP <- projectRaster(elevation, crs = CRS(SPNYCentral83), res = 50)
m1 <- mapview(elevationSP)
m2 <- mapview(sitesSP)
# m1 + m2

# slope
slope <- raster::terrain(elevation, out = "slope", unit='radians')
slope <- tan(slope)*100
slopeSP <- projectRaster(slope, elevationSP)
writeRaster(slopeSP,file.path("RASTERS","slop_degree.tiff"))

m3 <- mapview(slopeSP)
# m2 + m3

### NED 
NHD <- FedData::get_nhd(template=cl_5mi_buffWGS_sp, label='test')
flowlines <- NHD[["_Flowline"]]
flowlinesSP <- spTransform(flowlines, SPNYCentral83)
flowlinesSP_r <- rasterize(flowlinesSP, elevationSP)
# writeRaster(flowlinesSP_r,file.path("RASTERS","flowlines.tiff"))

m4 <- mapview(flowlinesSP_r)



##### TESTING BELOW HERE ##############
## Flowlines and distance to
flowlines <- NHD[["_Flowline"]]
flowlines_UTM18N83 <- spTransform(flowlines, CRS(proj4string(PASS_bound_UTM18N83_buff_0_r)))
flowlines_UTM18N83_r <- rasterize(flowlines_UTM18N83, PASS_bound_UTM18N83_buff_0_r)
flowlines_UTM18N83_r[flowlines_UTM18N83_r > 0] <- 1
flowlines_dist_UTM18N83 <- raster::distance(flowlines_UTM18N83_r)
flowlines_dist <- projectRaster(flowlines_dist_UTM18N83, crs = proj_dd)
flowlines_dist <- crop(flowlines_dist, PASS_bound_r)
flow_scale_factor <- ppside * (floor(dim(flowlines_dist)[1]/ppside))
flow_rescale <- flowlines_dist
dim(flow_rescale) <- c(flow_scale_factor, flow_scale_factor)
flowlines_dist <- resample(flowlines_dist, elevation)
flowlines_dist <- raster::scale(flowlines_dist)
# plot
plot(flowlines_dist)
plot(PASS_bound_sp, add = TRUE)
plot(PASS_sp, add = TRUE)