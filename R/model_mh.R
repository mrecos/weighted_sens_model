## packages
library("raster")
library('tidyverse')
library("sf")
library("sp")
library("mapview")

## Data locations
data_loc <- "U:/z_OldServer/Projects/NY/Oakdale_Fayette_model/GIS"
SHP_loc  <- file.path(data_loc,"SHP","Client_files")

## projections details
# unprojected WGS84
WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# New York Central SP NAD 83
SPNYCentral83 <- "+proj=tmerc +lat_0=40 +lon_0=-76.58333333333333 +k=0.9999375 +x_0=250000 +y_0=0 +datum=NAD83 +units=us-ft +no_defs"

## Read SHPs
sites    <- read_sf(file.path(SHP_loc, "ArchaeoSitesPt.shp"))
sites_sp <- as(sites, "Spatial")
cl       <- read_sf(file.path(SHP_loc, "ProposedFinalRouteCL_10-05-2017.shp"))
cl_sp    <- as(cl, "Spatial")

m1 <- mapview(sites)
m2 <- mapview(cl)
m1 + m2


##### testubg


#### reproject and buffer boundary for use with distance rasters.
PASS_bound_UTM18N83 <- st_transform(PASS_bound, crs = UTM18N83)
PASS_bound_UTM18N83_buff <- st_buffer(PASS_bound_UTM18N83, 1000)
PASS_bound_UTM18N83_buff_0_r <- raster(as(PASS_bound_UTM18N83_buff, "Spatial"), res = 10)
# PASS_bound_WGS84_buff <- st_transform(PASS_bound_UTM18N83_buff, proj_dd)
# plot
plot(PASS_bound_UTM18N83_buff)
plot(PASS_bound_UTM18N83, add = T)
