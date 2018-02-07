# require(rgeos)
# r <- template
# Sl <- gUnion(flowlinesSP, flowlinesSP)
# dd = gDistance(Sl, as(r,"SpatialPoints"), byid=TRUE)
# r[] = apply(dd,1,min)
# plot(r)


library("raster")
library("tidyverse")
library("doParallel")

GRD_loc <- "C:/R_local/weighted_sens_model"

SLOPE <- raster(file.path(GRD_loc, "GRD", "model_slope_degree.grd"))
H20 <- raster(file.path(GRD_loc, "GRD", "model_min_h20SP_dist.grd"))

### Slope weighting Models###
slp_rcl1 <- c(0, 3, 5, 8, 15)
slp_rcl2 <- c(3, 5, 8, 15, 99999)
slp_1 <- c(10, 5, 3, 2, 99)
slp_2 <- c(9, 8, 2, 1, 99)
slp_3 <- c(8, 6, 4, 2, 99)
slp_4 <- c(8, 9, 2, 1, 99)
slp_5 <- c(7, 8, 3, 2, 99)
slp_rcl_m <- cbind(slp_rcl1, slp_rcl2, slp_1, slp_2, slp_3, slp_4, slp_5) 

### Dist to h20 weighting Models###
h20_rcl1 <- c(0, 100, 200, 400, 800)
h20_rcl2 <- c(100,200, 400, 800, 9999)
h20_1 <- c(15, 8, 4, 2, 1)
h20_2 <- c(12, 10, 5, 2, 1)
h20_3 <- c(10, 8, 6, 4, 2)
h20_4 <- c(9, 14, 4, 2, 1)
h20_5 <- c(8, 12, 6, 3, 1)
h20_rcl_m <- cbind(h20_rcl1, h20_rcl2, h20_1, h20_2, h20_3, h20_4, h20_5) 

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

start.time <- Sys.time()

slp_stack <- stack()
### Reclassify slope to all of the weighting models ###
for (i in 1:as.integer(ncol(slp_rcl_m)-2)){
  print(i)
  rclmat <- cbind(slp_rcl_m[,1], slp_rcl_m[,2], slp_rcl_m[,i+2])
  print(rclmat)
  rc_slope <- reclassify(SLOPE, rclmat, progress='text')
  print("reclassified")
  fname <- file.path(GRD_loc, "GRD", "RECLASS", paste0("slp_rcl_", i, ".grd"))
  rname <- paste(c("slp_rcl_", i), collapse = '')
  print(fname)
  writeRaster(rc_slope, filename = fname,  progress='text', overwrite=TRUE)
  print("Wrote Raster")
  assign(rname, rc_slope)
}

slp_stack <- stack(slp_rcl_1, slp_rcl_2, slp_rcl_3, slp_rcl_4, slp_rcl_5)  

hyd_stack <- stack()
### Reclassify H20 Dist to all of the weighting models ###
for (i in 1:as.integer(ncol(h20_rcl_m)-2)){
  print(i)
  rclmat <- cbind(h20_rcl_m[,1], h20_rcl_m[,2], h20_rcl_m[,i+2])
  print(rclmat)
  rc_h20 <- reclassify(H20, rclmat,  progress='text')
  print("reclassified")
  fname <- file.path(GRD_loc, "GRD", "RECLASS", paste0("h20_rcl_", i, ".grd"))
  rname <- paste(c("h20_rcl_", i), collapse = '')
  print(fname)
  writeRaster(rc_h20, filename = fname,  progress='text', overwrite=TRUE)
  assign(rname, rc_h20)
}

hyd_stack <- stack(h20_rcl_1, h20_rcl_2, h20_rcl_3, h20_rcl_4, h20_rcl_5)  

sum_stack <- stack()
system.time({
  for (i in 1:as.integer(ncol(slp_rcl_m)-2)){       
    slp_tmp = raster(slp_stack,layer=i)
    print(paste(c("slp_rcl_", i, " summed with:"), collapse = ''))
    for (j in 1:as.integer(ncol(h20_rcl_m)-2)){
      h20_tmp = raster(hyd_stack, layer=j)
      z = slp_tmp + h20_tmp
      sum_stack <- stack(sum_stack, z)
      print(paste(c("h20_rcl_", j, " = layer", nlayers(sum_stack)), collapse = ''))	
    }
  }
})

#########cycle through sum_stack and save each raster######
system.time({
  for (i in 1:nlayers(sum_stack)) {
    out_fname <- file.path(GRD_loc, "GRD", "SUMMED", paste0("sum_", i, ".grd"))
    out_rast = raster(sum_stack, layer=i)
    print(paste(c("Acquired summed raster: ", i) , collapse = ''))
    writeRaster(out_rast, filename= out_fname, progress='text')
    print(paste(c("Saved summed raster: ", out_fname) , collapse = ''))
  }
})


gc()
removeTmpFiles(h=1)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

stopCluster(cl)



