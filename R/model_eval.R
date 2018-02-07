library("raster")
library("tidyverse")
library("sf")
library("sp")

data_loc <- "U:/z_OldServer/Projects/NY/Oakdale_Fayette_model/GIS"
SHP_loc  <- file.path(data_loc,"SHP","Client_files")
GRD_loc <- "C:/R_local/weighted_sens_model"

sites <- read_sf(file.path(SHP_loc, "ArchaeoSitesPt.shp")) %>%
  mutate(presence = 1)
sites_buff <- st_buffer(sites, 50) # 50 ft radius buffer
sites_buff_sp <- as(sites_buff, "Spatial")

system.time({
sum_stack <- stack()
for (i in 1:25) {
  in_rast <- raster(file.path(GRD_loc, "GRD", "SUMMED", paste0("sum_", i, ".grd")))
  sum_stack <- stack(sum_stack, in_rast)
}
})
names(sum_stack) <- seq(1:25)

extract_values <- NULL         
system.time({
  for (i in 1:nlayers(sum_stack)) {
    # PASS_out_fname <- paste(c(output_location, "Model_1_Output/", zone, "/PASS_sum_", i, ".grd"), collapse = '')
    sum_in_fname <- paste0("sum_", i, ".grd")
    sum_in = sum_stack[[i]]
    cat("Acquired summed raster: ", sum_in_fname, "/n")
    temp_extract <- raster::extract(sum_in, sites_buff_sp, progress='text',
                                    weights = TRUE, small = TRUE, df = TRUE) %>%
      rename_at(vars(starts_with("X")), funs(sub(pattern = "X\\d+", "sens", .))) %>%
      mutate(model_num = i)
    cat("Masked summed raster: ", i, "/n")
    extract_values <- rbind(extract_values, temp_extract)
    rm(sum_in)
    gc()
  }
})

write.csv(extract_values, file.path(GRD_loc, "RESULTS", "extracted_sens_values.csv" ))

extract_values %>%
  group_by(model_num) %>%
  summarise(median = median(sens)) %>%
  arrange(desc(median))
