library(sf)
library(tidyverse)
library(terra)
library(rstac)
library(mgcv)
library(scales)
library(foreach)
library(doParallel)

make_vsicurl_url_dem <- function(base_url) {
  paste0(
    "/vsicurl", 
    "?pc_url_signing=yes",
    "&pc_collection=alos-dem",
    "&url=",
    base_url
  )
}
make_vsicurl_url_esa <- function(base_url) {
  paste0(
    "/vsicurl", 
    "?pc_url_signing=yes",
    "&pc_collection=esa-worldcover",
    "&url=",
    base_url
  )
}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

setwd("/projappl/project_2007415/repos/AAA_treeline")

workers <- min(c(length(future::availableWorkers()),
                 future::availableCores()))
print(workers)

g <- st_read("/scratch/project_2007415/AAA_treelines/AAA_treeline_initial_grid_020.gpkg")

predictor_dir <- tempdir()

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

polids <- sample(g$polid, 2000) %>% sort

polids <- polids[!polids %in% results2$polid]
polids <- polids[!polids %in% results$polid]

polids <- g$polid[1:50000]


length(polids)
st <- Sys.time()
results <- mclapply(polids, 
                    mc.cores = workers-1,  
                    FUN = function(i){
  # i <- 101467
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
  
  e <- try({
    
    # ESA WORLDCOVER
    
    s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
    
    it_obj <- s_obj %>% 
      stac_search(collections = "esa-worldcover",
                  bbox = st_bbox(aoi_pr %>% st_transform(4326)),
                  datetime = "2021-01-01/2021-12-31",
                  limit = 1000) %>%
      get_request()
    
    juuh <- lapply(it_obj$features, function(ft){
      # print(ft$id)
      # ft <- it_obj$features[[1]]
      full_url <- make_vsicurl_url_esa(assets_url(ft) %>% sort)
      full_url <- full_url[endsWith(full_url, "_Map.tif")]
      
      juuh <- lapply(seq_len(length(full_url)), function(nr){
        # nr <- 1
        file_name <- paste0(tempfile(), ".tif")
        
        e <- try({
          gdal_utils(
            util = "warp",
            source = full_url[[nr]],
            destination = file_name,
            options = c(
              "-t_srs", sf::st_crs(aoi_pr)$wkt,
              "-te", sf::st_bbox(aoi_pr),
              "-tr", c(10, 10)
            )
          )
        }, silent = TRUE)
        if(class(e)[[1]] == "try-error"){
          return(NULL)
        } else {
          return(file_name)
        }
      }) %>% unlist
      return(juuh)
    }) %>% unlist
    
    if(!is.null(juuh)){
      esas <- juuh
      
      esa <- lapply(esas, function(x){
        esa <- rast(x)
        esa[esa == 0] <- NA
        # dem <- resample(dem, r)
        return(esa)
      })
      
      esa <- sprc(esa)
      esa <- mosaic(esa, fun = "max")
      esa[esa < 0] <- 0
      esa[is.na(esa)] <- 0
      
      # plot(esa)
      
      # DEM
      
      s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
      
      it_obj <- s_obj %>% 
        stac_search(collections = "alos-dem",
                    bbox = st_bbox(aoi_pr %>% st_transform(4326)),
                    limit = 1000) %>%
        get_request()
      
      juuh <- lapply(it_obj$features, function(ft){
        # ft <- it_obj$features[[1]]
        
        # print(ft$id)
        full_url <- make_vsicurl_url_dem(assets_url(ft) %>% sort)
        full_url <- full_url[endsWith(full_url, "_DSM.tif")]
        
        juuh <- lapply(seq_len(length(full_url)), function(nr){
          file_name <- paste0(tempfile(), ".tif")
          
          e <- try({
            gdal_utils(
              "warp",
              source = full_url[[nr]],
              destination = file_name,
              options = c(
                "-t_srs", sf::st_crs(aoi_pr)$wkt,
                "-te", sf::st_bbox(aoi_pr),
                "-tr", c(30, 30)
              )
            )
          }, silent = TRUE)
          if(class(e)[[1]] == "try-error"){
            return(NULL)
          } else {
            return(file_name)
          }
        }) %>% unlist
        return(juuh)
      }) %>% unlist
      
      dems <- juuh
      
      dem <- lapply(dems, function(x){
        dem <- rast(x)
        dem[dem == 0] <- NA
        # dem <- resample(dem, r)
        return(dem)
      })
      
      dem <- sprc(dem)
      dem <- mosaic(dem)
      dem[dem < 0] <- 0
      dem[is.na(dem)] <- 0
      
      # plot(dem)
      
      unlink(esas)
      unlink(dems)
      
      esa <- aggregate(esa, 3, getmode)
      esa <- project(esa, dem, method="near")
      # plot(esa)
      
      esa[esa %in% c(0,40,50,80,90)] <- NA
      rmask <- focal(esa, 7, min, na.rm = FALSE, expand = TRUE, fillvalue = 100)
      esa <- mask(esa, rmask)
      esa <- ifel(esa == 10, 1, 0)
      # plot(esa)
      r <- c(dem, esa)
      names(r) <- c("dem","esa")
      
      df <- as.data.frame(r, na.rm = TRUE)
      if(nrow(df) > 1000){
        # summary(df)
        if(mean(df$esa) > 0.001 & mean(df$esa) < 0.999 & diff(range(df$dem)) >= 20){
          # df %>%
          #   sample_n(size = 10000) %>%
          #   ggplot(aes(y = esa, x = dem)) +
          #   geom_smooth(method = "gam",
          #               method.args=list(family="binomial"))
          
          m <- gam(esa ~ s(dem), data = df %>% sample_n(size = ifelse(nrow(df) > 20000, 20000, nrow(df))),
                   family = "binomial")
          s <- summary(m)
          
          pr <- tibble(dem = seq(from = min(df$dem), to = max(df$dem), by = 1)) %>% 
            mutate(fprob = predict(m, ., type = "response")) %>% 
            mutate(fprob_scaled = rescale(fprob, to = c(0,1)))
          
          if(max(pr$fprob) >= 0.1 & min(pr$fprob) <= 0.9){
            results <- tibble(polid = as.integer(i),
                              lat = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(Y),
                              lon = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
                              minele = min(df$dem),
                              maxele = max(df$dem),
                              landprop = nrow(df)/ncell(dem),
                              forestprop = mean(df$esa),
                              forestq95 = quantile(df %>% filter(esa == 1) %>% pull(dem), 0.95),
                              forestq99 = quantile(df %>% filter(esa == 1) %>% pull(dem), 0.99),
                              maxforestele = quantile(df %>% filter(esa == 1) %>% pull(dem), 1),
                              gamR2 = s$r.sq,
                              maxfprob = max(pr$fprob),
                              minfprob = min(pr$fprob),
                              forestprob99ele = pr %>% filter(fprob > 0.99) %>% pull(dem) %>% max(),
                              forestprob95ele = pr %>% filter(fprob > 0.95) %>% pull(dem) %>% max(),
                              forestprob90ele = pr %>% filter(fprob > 0.9) %>% pull(dem) %>% max(),
                              forestprob75ele = pr %>% filter(fprob > 0.75) %>% pull(dem) %>% max(),
                              forestprob50ele = pr %>% filter(fprob > 0.5) %>% pull(dem) %>% max(),
                              forestprob25ele = pr %>% filter(fprob > 0.25) %>% pull(dem) %>% max(),
                              forestprob10ele = pr %>% filter(fprob > 0.1) %>% pull(dem) %>% max(),
                              forestprob05ele = pr %>% filter(fprob > 0.05) %>% pull(dem) %>% max(),
                              forestprob01ele = pr %>% filter(fprob > 0.01) %>% pull(dem) %>% max(),
                              forestprobscaled99ele = pr %>% filter(fprob_scaled > 0.99) %>% pull(dem) %>% max(),
                              forestprobscaled95ele = pr %>% filter(fprob_scaled > 0.95) %>% pull(dem) %>% max(),
                              forestprobscaled90ele = pr %>% filter(fprob_scaled > 0.9) %>% pull(dem) %>% max(),
                              forestprobscaled75ele = pr %>% filter(fprob_scaled > 0.75) %>% pull(dem) %>% max(),
                              forestprobscaled50ele = pr %>% filter(fprob_scaled > 0.5) %>% pull(dem) %>% max(),
                              forestprobscaled25ele = pr %>% filter(fprob_scaled > 0.25) %>% pull(dem) %>% max(),
                              forestprobscaled10ele = pr %>% filter(fprob_scaled > 0.1) %>% pull(dem) %>% max(),
                              forestprobscaled05ele = pr %>% filter(fprob_scaled > 0.05) %>% pull(dem) %>% max(),
                              forestprobscaled01ele = pr %>% filter(fprob_scaled > 0.01) %>% pull(dem) %>% max()) %>% 
              mutate(across(landprop:forestprobscaled01ele, ~round(.x, 3))) %>% 
              mutate(across(where(is.double), ~ifelse(is.infinite(.x), NA, .x)))
          } else {
            results <- tibble(polid = as.integer(i),
                              lat = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(Y),
                              lon = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
                              minele = min(df$dem),
                              maxele = max(df$dem),
                              landprop = nrow(df)/ncell(dem),
                              forestprop = mean(df$esa)) %>% 
              mutate(across(landprop:forestprop, ~round(.x, 3))) %>% 
              mutate(across(where(is.double), ~ifelse(is.infinite(.x), NA, .x)))
          }
          return(results)
        } else {
          results <- tibble(polid = as.integer(i),
                            lat = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(Y),
                            lon = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
                            minele = min(df$dem),
                            maxele = max(df$dem),
                            landprop = nrow(df)/ncell(dem),
                            forestprop = mean(df$esa)) %>% 
            mutate(across(landprop:forestprop, ~round(.x, 3))) %>% 
            mutate(across(where(is.double), ~ifelse(is.infinite(.x), NA, .x)))
          return(results)
        }
      } else {
        mm <- minmax(dem)
        results <- tibble(polid = as.integer(i),
                          lat = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(Y),
                          lon = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
                          minele = mm[1,1],
                          maxele = mm[2,1])
        return(results)
      }
    } else {
      results <- tibble(polid = as.integer(i),
                        lat = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(Y),
                        lon = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(X))
      return(results)
    }
  }, silent = TRUE)
  
  if(class(e)[[1]] == "try-error"){
    results <- tibble(polid = as.integer(i),
                      lat = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(Y),
                      lon = st_centroid(aoi) %>% st_coordinates() %>% as.data.frame() %>% pull(X),
                      error = TRUE)
    return(results)
  }
})
Sys.time() - st

results <- results %>% bind_rows()
results2 <- bind_rows(results2,
                      results[lapply(results, function(x) class(x)[1]) %>% unlist != "try-error"] %>% bind_rows())

results %>% bind_rows() %>% write_csv("output/sample1-50k_020.csv")
# lapply(results, function(x) class(x)[1]) %>% unlist

po <- results %>% 
  arrange(polid) %>% 
  mutate(rangeele = maxele - minele) %>% 
  filter(rangeele > 100) %>% 
  filter(maxele != forestprob50ele) %>% 
  st_as_sf(coords = c("lon","lat"), crs = 4326)
results %>% view
plot(po[,"forestprob50ele"], pch = 20)

po <- results %>% 
  arrange(polid) %>% 
  st_as_sf(coords = c("lon","lat"), crs = 4326)
po %>% st_write("output/random3.gpkg")
