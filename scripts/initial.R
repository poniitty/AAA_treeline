library(terra)
library(sf)
library(tidyverse)


dl_dir <- "/scratch/project_2007415/AAA_treelines/"

tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/input/dem_latlong.sdat")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))
tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/input/dem_latlong.mgrd")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))
tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/input/dem_latlong.prj")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))
tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/input/dem_latlong.sgrd")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))


tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/input/landseamask_buffer.sdat")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))
tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/input/landseamask_buffer.mgrd")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))
tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/input/landseamask_buffer.prj")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))
tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/input/landseamask_buffer.sgrd")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))


tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",
                          10,"_1981-2010_V.2.1.tif")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))
tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",
                          11,"_1981-2010_V.2.1.tif")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))
tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",
                          6,"_1981-2010_V.2.1.tif")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))
tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",
                          1,"_1981-2010_V.2.1.tif")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))

tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_gst_1981-2010_V.2.1.tif")
download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))

r <- rast(paste0(dl_dir, "CHELSA_gst_1981-2010_V.2.1.tif"))
dem <- rast(paste0(dl_dir, "dem_latlong.sdat"))
plot(dem)
oce <- rast(paste0(dl_dir, "landseamask_buffer.sdat"))
plot(oce)
dem[is.na(oce)] <- NA
dem[oce == 1] <- NA
plot(dem)

dem %>% round %>% writeRaster(paste0(dl_dir, "CHELSA_dem.tif"), datatype = "INT2S")

r <- rast(paste0(dl_dir, "CHELSA_gst_1981-2010_V.2.1.tif"))
r <- ifel(r > 8, 1, 0)
plot(r)
writeRaster(r, "/scratch/project_2007415/AAA_treelines/8dec_gst.tif")

list.files(dl_dir)

r <- rast(paste0(dl_dir, "CHELSA_gst_1981-2010_V.2.1.tif"))
r[is.na(r)] <- 0
r <- ifel(r <= 9, 1, 0)
plot(r)

warmt <- rast(paste0(dl_dir, "CHELSA_bio10_1981-2010_V.2.1.tif"))
r[warmt > 12] <- 0

gst <- rast(paste0(dl_dir, "CHELSA_gst_1981-2010_V.2.1.tif"))
r[gst < 6.7] <- 1

tmin <- rast(paste0(dl_dir, "CHELSA_bio6_1981-2010_V.2.1.tif"))
r[tmin > 0] <- 0


# Sea mask

oce <- rast(paste0(dl_dir, "landseamask_buffer.sdat"))

plot(oce)
r[is.na(oce)] <- NA
plot(r)
gc()

writeRaster(r, "/scratch/project_2007415/AAA_treelines/AAA_initial.tif", datatype = "INT1U", overwrite = TRUE)

r[r != 1] <- NA

p <- as.polygons(r) %>% st_as_sf() %>% st_cast("POLYGON") %>% st_transform("ESRI:53009")
pb <- p %>% st_buffer(20000)
pb <- st_crop(pb, st_bbox(p) %>% st_as_sfc %>% st_as_sf %>% st_transform("ESRI:53009"))
pb <- pb %>% st_cast("POLYGON")
pb$xdim <- lapply(1:nrow(pb), function(i){
  # i <- 1
  bb <- pb %>% slice(i) %>% st_bbox()
  return(diff(c(bb$xmin,bb$xmax)))
}) %>% unlist

pb %>% arrange(desc(xdim)) %>% slice(1:100) %>% st_geometry() %>% plot
pb %>% arrange(desc(xdim)) %>% slice(-1) %>% st_geometry() %>% plot
pb %>% arrange(desc(xdim)) %>% slice(-1) %>% st_union() %>% st_geometry() %>% plot

pb2 <- st_crop(pb, st_bbox(p) %>% st_as_sfc %>% st_as_sf %>% st_transform("ESRI:53009")) %>% st_union()

plot(st_geometry(p %>% slice(1:10) %>% st_buffer(20000)))
pb <- st_buffer(p %>% slice(1:10), 20000) %>% st_union() %>% st_cast("POLYGON")

p %>% slice(1:20) %>% st_buffer(20000) %>% st_write("/scratch/project_2007415/AAA_treelines/TEST1.gpkg")
p %>% st_buffer(20000) %>% st_union() %>% st_cast("POLYGON") %>% st_write("/scratch/project_2007415/AAA_treelines/AAA_maxpoly.gpkg", append = FALSE)
pb2 %>% st_cast("POLYGON") %>% st_write("/scratch/project_2007415/AAA_treelines/AAA_maxpoly2.gpkg", append = FALSE)
pb %>% arrange(desc(xdim)) %>% slice(-1) %>% st_union() %>% st_write("/scratch/project_2007415/AAA_treelines/AAA_maxpoly3.gpkg", append = FALSE)

##############################################################################

sf_use_s2(FALSE)
p <- st_read("/scratch/project_2007415/AAA_treelines/AAA_maxpoly_comb.gpkg") %>% 
  st_transform(4326) %>% 
  st_union()
sf_use_s2(TRUE)
# Glaciers
glc <- bind_rows(st_read("/scratch/project_2007415/glaciers/glims_download_58903/glims_polygons.shp"),
                 st_read("/scratch/project_2007415/glaciers/glims_download_73332/glims_polygons.shp")) %>% 
  st_zm(drop = T) %>% 
  st_cast("POLYGON")

ant <- st_read("/scratch/project_2007415/glaciers/antarctica-icesheet-polygons-3857/icesheet_polygons.shp") %>% 
  st_transform(4326) %>% 
  st_simplify(dTolerance = 0.001)

gre <- st_read("/scratch/project_2007415/glaciers/Greenland_Hydrologic_Sub_Basins.shp") %>% 
  st_transform(4326) %>% 
  st_cast("POLYGON")

sf_use_s2(FALSE)
ant <- ant %>% st_union
gre <- gre %>% st_union

p2 <- st_crop(p, ant)
p2 <- st_difference(p2, ant)
plot(st_geometry(p2))
plot(st_geometry(ant), add = TRUE, col = "red")

p <- st_difference(p, gre)
p <- st_difference(p, glc %>% 
                     arrange(desc(area)) %>% 
                     slice_head(n = 100) %>% 
                     st_union())

g <- st_make_grid(p, cellsize = 0.3333) %>% 
  st_as_sf()

g1 <- st_difference(g, ant %>% st_bbox() %>% st_as_sfc() %>% st_as_sf)
g2 <- st_intersection(g, ant %>% st_bbox() %>% st_as_sfc() %>% st_as_sf)

# plot(g1)
# plot(g2, add = T)

g1 <- g1[world %>% filter(continent != "Antarctica"),]
g2 <- g2[world %>% filter(continent == "Antarctica"),]

p <- st_difference(p, ant)
gc()


gt <- g %>% 
  mutate(id = rep(1:ceiling(nrow(g)/10000), each = 10000, length.out = nrow(g)))

gtt <- split(gt, gt$id)
length(gtt)

gtt <- lapply(gtt, function(x){
  print(x$id[[1]])
  gttt <- st_intersects(x, p)
  return(x[!is.na(as.numeric(gttt)),])
}) %>% bind_rows()

g <- gtt

sf_use_s2(FALSE)
g <- g[p,]
plot(g)

library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires, lib.loc = "/projappl/project_2003061/Rpackages/")

world <- ne_countries(scale = "large", returnclass = "sf") %>% 
  st_transform(crs = 4326)

g <- g[world,]

g <- g %>%
  mutate(polid = 1:nrow(.))

g %>% st_write("/scratch/project_2007415/AAA_treelines/AAA_treeline_initial_grid_033.gpkg", append = FALSE)
