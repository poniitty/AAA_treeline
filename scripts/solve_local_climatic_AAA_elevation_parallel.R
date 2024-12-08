library(sf)
library(tidyverse)
library(terra)
library(mgcv)
library(scales)
library(foreach)
library(doParallel)


setwd("/projappl/project_2007415/repos/AAA_treeline")
dl_dir <- "/scratch/project_2007415/AAA_treelines/"

workers <- min(c(length(future::availableWorkers()),
                 future::availableCores()))
print(workers)

g <- st_read("/scratch/project_2007415/AAA_treelines/AAA_treeline_initial_grid_020.gpkg")
dem <- rast(paste0(dl_dir, "CHELSA_dem.tif"))
tmin <- rast(paste0(dl_dir, "CHELSA_bio6_1981-2010_V.2.1.tif"))
warmt <- rast(paste0(dl_dir, "CHELSA_bio10_1981-2010_V.2.1.tif"))
mat <- rast(paste0(dl_dir, "CHELSA_bio1_1981-2010_V.2.1.tif"))

utmall <- st_read("data/utm_zones.gpkg") %>% 
  filter(ZONE != 0)

g <- bind_cols(g,
               g %>% st_centroid() %>% st_coordinates())

polids <- g %>% 
  filter(X > 16,
         X < 32,
         Y > 62,
         Y < 72) %>% 
  pull(polid)

polids <- polids[!polids %in% results2$polid]
polids <- sample(g$polid, 2000) %>% sort

polids <- polids[!polids %in% results2$polid]
polids <- polids[!polids %in% results$polid]

polids <- g$polid

polids <- sample(g$polid, 10000)
length(polids)
st <- Sys.time()
results <- mclapply(polids, 
                    mc.cores = workers,  
                    FUN = function(i){
                      # i <- 101467
                      # i <- 106991
                      print(i)
                      aoi <- g %>% filter(polid == i)
                      
                      # WGS84 UTM zones to set the correct projection
                      utm <- utmall[aoi %>% st_centroid(),] # Which zone the study points falls in
                      
                      if(nrow(utm) == 0){
                        utm <- utmall[st_nearest_feature(aoi, utmall),]
                      }
                      
                      lat <- st_coordinates(aoi %>% st_centroid())[,"Y"]
                      utm$ZONE <- ifelse(nchar(utm$ZONE) == 1, paste0("0",utm$ZONE), utm$ZONE)
                      epsg <- as.numeric(ifelse(lat > 0, paste0(326, utm$ZONE), paste0(327, utm$ZONE)))
                      
                      aoi_pr <- aoi %>% st_transform(epsg) %>% 
                        st_buffer(5000)
                      
                      aoi <- aoi_pr %>% 
                        st_transform(4326)
                      
                      r <- c(dem %>% crop(aoi),
                             tmin %>% crop(aoi),
                             warmt %>% crop(aoi),
                             mat %>% crop(aoi))
                      names(r) <- c("dem","tmin","warmt","mat")
                      df <- as.data.frame(r, na.rm = TRUE)
                      
                      if(nrow(df) > 100 & diff(range(df$dem)) > 50){
                        # summary(df)
                        # gg <- df %>% ggplot(aes(x = dem, y = mat)) +
                        #   geom_point()+
                        #   geom_smooth(method = 'lm', formula = 'y ~ x')
                        # print(gg)
                        
                        tminm <- summary(lm(tmin ~ dem, data = df))
                        warmtm <- summary(lm(warmt ~ dem, data = df))
                        matm <- summary(lm(mat ~ dem, data = df))
                        
                        results <- tibble(polid = as.integer(i),
                                          lat = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(Y),
                                          lon = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
                                          area = nrow(df),
                                          minele = min(df$dem),
                                          maxele = max(df$dem),
                                          tmin = mean(df$tmin),
                                          warmt = mean(df$warmt),
                                          mat = mean(df$mat),
                                          tmin_r2 = tminm$r.squared,
                                          tmin_interc = tminm$coefficients[1,1],
                                          tmin_slope = tminm$coefficients[2,1],
                                          warmt_r2 = warmtm$r.squared,
                                          warmt_interc = warmtm$coefficients[1,1],
                                          warmt_slope = warmtm$coefficients[2,1],
                                          mat_r2 = matm$r.squared,
                                          mat_interc = matm$coefficients[1,1],
                                          mat_slope = matm$coefficients[2,1])
                        return(results)
                      } else {
                        return(NULL)
                      }
                    })
Sys.time() - st

results2 <- results[lapply(results, function(x) class(x)[1]) %>% unlist != "try-error"] %>% bind_rows()
results <- results %>% bind_rows()
results2 <- bind_rows(results2,
                      results[lapply(results, function(x) class(x)[1]) %>% unlist != "try-error"] %>% bind_rows())

results %>% bind_rows() %>% write_csv("output/temp_slopes.csv")
# lapply(results, function(x) class(x)[1]) %>% unlist

po <- results %>% 
  arrange(polid) %>% 
  mutate(rangeele = maxele - minele) %>% 
  filter(rangeele > 100,
         mat_r2 > 0.6) %>% 
  st_as_sf(coords = c("lon","lat"), crs = 4326)
results %>% view
plot(po[,"mat_slope"], pch = 20)

po <- results %>% 
  arrange(polid) %>% 
  st_as_sf(coords = c("lon","lat"), crs = 4326)
po %>% st_write("output/temp_slopes.gpkg", append = FALSE)
