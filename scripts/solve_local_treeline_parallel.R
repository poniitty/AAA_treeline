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

workers <- min(c(length(future::availableWorkers()),
                 future::availableCores()))
print(workers)

g <- st_read("/scratch/project_2007415/AAA_treelines/AAA_treeline_initial_grid.gpkg")

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

results <- mclapply(polids[1:10], mc.cores = 5, 
                    FUN = function(i){
  # i <- 29053
  
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
  
  aoi_pr <- aoi %>% st_transform(epsg)
  
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
    file_names <- gsub("TIF$","tif",basename(full_url))
    
    juuh <- lapply(seq_len(length(full_url)), function(nr){
      # nr <- 1
      e <- try({
        gdal_utils(
          util = "warp",
          source = full_url[[nr]],
          destination = paste0(predictor_dir,"/",file_names[[nr]]),
          options = c(
            "-t_srs", sf::st_crs(aoi_pr)$wkt,
            "-te", sf::st_bbox(aoi_pr),
            "-tr", c(10, 10)
          )
        )
      }, silent = TRUE)
      if(class(e)[[1]] == "try-error"){
        return(FALSE)
      } else {
        return(TRUE)
      }
    })
  })
  
  esas <- list.files(predictor_dir, pattern = "_Map.tif", full.names = TRUE)
  
  esas <- lapply(esas, function(x){
    esa <- rast(x)
    esa[esa == 0] <- NA
    # dem <- resample(dem, r)
    return(esa)
  })
  
  esa <- sprc(esas)
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
    file_names <- gsub("TIF$","tif",basename(full_url))
    
    juuh <- lapply(seq_len(length(full_url)), function(nr){
      e <- try({
        gdal_utils(
          "warp",
          source = full_url[[nr]],
          destination = paste0(predictor_dir,"/",file_names[[nr]]),
          options = c(
            "-t_srs", sf::st_crs(aoi_pr)$wkt,
            "-te", sf::st_bbox(aoi_pr),
            "-tr", c(30, 30)
          )
        )
      }, silent 