library("raster")
library("tidyverse")
library("sf")
library("sp")

data_loc <- "U:/z_OldServer/Projects/NY/Oakdale_Fayette_model/GIS"
SHP_loc  <- file.path(data_loc,"SHP","Client_files")
GRD_loc <- "C:/R_local/weighted_sens_model"

sites <- read_sf(file.path(SHP_loc, "ArchaeoSitesPt_PreHist.shp")) %>%
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
  summarise(median = median(sens),
            n = n()) %>%
  arrange(desc(median))

get_random_sample <- function(model_num, stack, sample_size){
  raster <- stack[[model_num]]
  cat("Sampling raster:", model_num, "\n")
  rand_smpl <- sampleRandom(raster, sample_size)
}

tabulate_sample <- function(smpl, steep = 99){
  smpl_freq <- as.data.frame(table(smpl), stringsAsFactors = FALSE)  %>%
    mutate_if(is.character, as.numeric)
  if(max(smpl_freq[,1]) > 98){
    over <- smpl_freq[which(smpl_freq[,1] >= 98), ]
    over98_sum <- sum(over[,2])
    smpl_freq[1,2] <- smpl_freq[1,2] + over98_sum
    smpl_freq <- smpl_freq[-which(smpl_freq[,1] >= 98), ]
  }
  colnames(smpl_freq) <- c("value", "count")	
  return(smpl_freq)
}
cumulative_sum <- function(rand_smpl_freq,site_smpl_freq){
  rand_smpl_CumFreq = cumsum(rand_smpl_freq[,2])
  site_smpl_CumFreq = cumsum(site_smpl_freq[,2])
  rand_smpl_revCumFreq <- rev(cumsum(rev(rand_smpl_freq[,2])))
  site_smpl_revCumFreq <- rev(cumsum(rev(site_smpl_freq[,2])))
  rand_smpl_revCumPcnt <- (rand_smpl_revCumFreq/sum(rand_smpl_freq[,2]))*100
  site_smpl_revCumPcnt <- (site_smpl_revCumFreq/sum(site_smpl_freq[,2]))*100
  rand_smpl_CumPcnt <- (rand_smpl_CumFreq/sum(rand_smpl_freq[,2]))*100
  site_smpl_CumPcnt <- (site_smpl_CumFreq/sum(site_smpl_freq[,2]))*100
  rand_smpl_freq <- cbind(rand_smpl_freq, rand_smpl_revCumFreq, rand_smpl_revCumPcnt)
  site_smpl_freq <- cbind(site_smpl_freq, site_smpl_revCumFreq, site_smpl_revCumPcnt)
  freq_merge <- merge(rand_smpl_freq, site_smpl_freq, all = TRUE, by = 'value')
  freq_merge[is.na(freq_merge)] <- 0
  freq_merge$kg <- (1 - (freq_merge$rand_smpl_revCumPcnt/freq_merge$site_smpl_revCumPcnt))
  return(freq_merge)
}

model_sens1 <- extract_values %>%
  as_tibble() %>%
  tidyr::nest(-model_num) %>%
  mutate(rand_smpl = map(model_num, get_random_sample, sum_stack, 50000))

######### min(asb()) needs to be closest without going over (price is right rules) intead of absolute closest
model_sens2 <- model_sens1 %>%
  mutate(rand_smpl_freq = map(rand_smpl, tabulate_sample),
         site_smpl_freq = map(data, ~tabulate_sample(.[["sens"]])),
         freq_merge     = map2(rand_smpl_freq, site_smpl_freq, cumulative_sum),
         # max_kg         = map_dbl(freq_merge, ~ max(.[["kg"]])),
         # kg_sites_pcnt  = map2_dbl(freq_merge, max_kg, ~ .x[.x[,"kg"] == .y, "site_smpl_revCumPcnt"]),
         # kg_backg_pcnt  = map2_dbl(freq_merge, max_kg, ~ .x[.x[,"kg"] == .y, "rand_smpl_revCumPcnt"]),
         # kg_threshold   = map2_dbl(freq_merge, max_kg, ~ .x[.x[,"kg"] == .y, "value"]),
         bkg_sites_high    = map_dbl(freq_merge,  ~ .[which.min(abs(.[,"site_smpl_revCumPcnt"]-75)), "rand_smpl_revCumPcnt"]),
         sites_high_thold  = map2_dbl(freq_merge, bkg_sites_high, ~ .x[.x[,"rand_smpl_revCumPcnt"] == .y, "value"]),
         bkg_sites_mod     = map_dbl(freq_merge,  ~ .[which.min(abs(.[,"site_smpl_revCumPcnt"]-92)), "rand_smpl_revCumPcnt"]),
         sites_mod_thold   = map2_dbl(freq_merge, bkg_sites_mod, ~ .x[.x[,"rand_smpl_revCumPcnt"] == .y, "value"])) 

max_kg_results <- model_sens2 %>%
  dplyr::select(model_num, bkg_sites_high, sites_high_thold, bkg_sites_mod, sites_mod_thold)

# view selected model's freq table
model_sens2 %>% filter(model_num == 15) %>% select(freq_merge) %>% unnest()





