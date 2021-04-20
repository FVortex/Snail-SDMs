  library(dplyr)
  library(sdm) 
  library(sdmpredictors)
  library(dismo)
  library(spThin)
  library(raster)
  library(maps)
  #custom data
  tb <- readxl::read_xlsx(path = '/home/mikhail/Documents/slug_sdm/invasive_snails_SDM.xlsx')
  tb %>% filter(species == 'B_cylindrica') -> BCtb
  tb %>% filter(species == 'X_derbentina') -> XDtb
  #gbif data 
  ##Brephulospis cylindrica
  BCgb <- gbif(genus = 'Brephulopsis', 
               species = "cylindrica", geo = T
  )
  !is.na(BCgb$lon) ->BCind
  
  BCgb[BCind,] -> BCgbif
  
  BCall <- data.frame(lon = as.numeric(c(BCgbif$lon, BCtb$longitude)),
                      lat = as.numeric(c(BCgbif$lat, BCtb$latitude)),
                      species = 1)
  # BCall %>% ##View
  #no europe
  BCall %>% filter(lon >20) ->BCall
  BCall %>% plot
  
  ##Xeropicta derbentina
  XDgb <- gbif(genus = 'Xeropicta', 
               species = "derbentina", geo = T
  )
  !is.na(XDgb$lon) ->XDind
  XDind %>% table()
  XDgb[XDind,] -> XDgbif
  
  # XDgbif %>% ##View
  
  XDall <- data.frame(lon = as.numeric(c(XDgbif$lon, XDtb$longitude)),
                      lat = as.numeric(c(XDgbif$lat, XDtb$latitude)),
                      species = 1)
  # BCall %>% ##View
  XDall %>% filter(lon > 20) ->XDall
  # # # # #XDall %>% filter(between(lon, 30, 55)) -> XDall
  # # # ##XDall %>% filter(between(lat, 44, 55)) -> XDall
  XDall[,1:2] %>% plot
  
  # boundaries are set
  
  prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  # 
  #  
  # ru_adm <- getData("GADM", country=c("Russia", "Ukraine"), level=1)
  # ru_adm %>% plot
  xs <- c(0, 80)
  ys <- c(30, 67)
  e = extent( xs , ys)
  
  #predictors 
  
  # first to select layers
  lays <- list_layers(datasets="WorldClim", 
                      terrestrial = T,
                      freshwater = FALSE,
                      marine = FALSE)
  # lays %>% ##View
  ##no minors montlhy variabls
  dat <- load_layers(lays[1:20], datadir = tempdir())
  # save(dat, file = '/home/mikhail/Documents/slug_sdm/worldclim_europe.rda')
  # load( '/home/mikhail/Documents/slug_sdm/worldclim_europe.rda')
  roi <- crop(dat, e)
  names(roi) <- lays$name
  # roi %>% plot
  
  raster::scale(roi) -> roi
  
  
  #occurences on top of predictors
  par(mfrow = c(1,2))
  plot(roi[[1]], main = 'Species occurences on scaled elevation map', asp = 1)
  points(x = BCall$lon, y = BCall$lat, cex = 0.5, pch = 4, col = 'red', asp = 1)
  points(x = XDall$lon, y = XDall$lat, cex = 0.5, pch = 4, col = 'blue', asp = 1)
  #close
  plot(roi[[1]], main = 'Species occurences on scaled elevation map', asp = 1, 
       ylim = c(37, 53), 
       xlim = c(23, 47))
  points(x = BCall$lon, y = BCall$lat, cex = 0.5, pch = 4, col = 'red', asp = 1)
  points(x = XDall$lon, y = XDall$lat, cex = 0.5, pch = 4, col = 'blue', asp = 1)
  #how many observations ono ur hands
  nrow(BCall); nrow(XDall)
  
  #
  # Generate 100 random points for no-presence data
  set.seed(90)
  bg_spdf <- sampleRandom(roi, size = 150, sp = T)
  ###B. cylindrca comes first
  
  #     
  # #coords are to be extracted
  BC_df <- data.frame(decimalLatitude = BCall$lat,
                      decimalLongitude = BCall$lon)
  # #and crs came from ppredictor object whiteSea
  BC_spdf <- SpatialPointsDataFrame(coords = BC_df,
                                    data = BC_df,
                                    proj4string = roi@crs)
  # 
  # # Combine presence/no-presence data
  BC_and_bg <- SpatialPointsDataFrame(coords = rbind(bg_spdf@coords, BC_spdf@coords), 
                                      data = data.frame(Occurence = c(rep(0, nrow(bg_spdf)), rep(1, nrow(BC_spdf)))),
                                      proj4string = roi@crs)
  #     
  #     
  # # next as it goes in sdm vignette
  #     
  # # PART ML 
  
  d_BC <- sdmData(formula = Occurence ~ ., train = BC_and_bg, predictors= roi)
  set.seed(90)
  m1_BC <-sdm(Occurence ~ ., data = d_BC, methods = 'rf', replication='sub',  test.percent = 30, n = 10)
  #     
  # #model diagnostics
  roc(m1_BC)
  m1_BC
  #selecting the best-performing model
  eval_BC <- getEvaluation(m1_BC)
  #which are the best model here ?
  BC_which <- which.max(eval_BC$AUC)
  
  #     
  p1_BC <-predict(m1_BC[[BC_which]],newdata= roi,filename='p1.img', overwrite = T)# many commonly used raster format is supported (throughplot(p1)
  plot(p1_BC, main = 'Predicted distribution for B. cylindrica')
  
  # ###D. xeropicta comes second
  
  #     
  # #coords are to be extracted
  XD_df <- data.frame(decimalLatitude = XDall$lat,
                      decimalLongitude = XDall$lon)
  # #and crs came from ppredictor object whiteSea
  XD_spdf <- SpatialPointsDataFrame(coords = XD_df,
                                    data = XD_df,
                                    proj4string = roi@crs)
  # 
  # # Combine presence/no-presence data
  XD_and_bg <- SpatialPointsDataFrame(coords = rbind(bg_spdf@coords, XD_spdf@coords), 
                                      data = data.frame(Occurence = c(rep(0, nrow(bg_spdf)), rep(1, nrow(XD_spdf)))),
                                      proj4string = roi@crs)
  #     
  # # PART ML 
  
  d_XD <- sdmData(formula = Occurence ~ ., train = XD_and_bg, predictors= roi)
  set.seed(90)
  m1_XD <-sdm(Occurence ~ ., data = d_XD, methods = 'rf', replication='sub',  test.percent = 30, n = 10)
  #     
  # #model diagnostics
  roc(m1_XD)
  m1_XD
  #selecting the best-performing model
  eval_XD <- getEvaluation(m1_XD)
  #which are the best model here ?
  XD_which <- which.max(eval_XD$AUC)
  
  #     
  p1_XD <-predict(m1_XD[[XD_which]],newdata= roi,filename='p2.img', overwrite = T)# many commonly used raster format is supported (throughplot(p1)
  plot(p1_XD, main = 'Predicted distribution for X. derbentina')
  
  ##al pics
  #
  par(mfrow = c(1,2))
  plot(p1_BC, main = 'Predicted distribution for B. cylindrica', asp = 1)
  points(x = BCall$lon, y = BCall$lat, cex = 0.1, pch = 4, col = 'red', asp = 1)
  points(x = XDall$lon, y = XDall$lat, cex = 0.1, pch = 4, col = 'blue', asp = 1)
  legend('topleft', legend = c('B. cylindrica', 'X. derbetnita'), col = c('red', 'blue'), pch = 4, bty = 'n')
  #
  plot(p1_XD, main = 'Predicted distribution for X. derbentina', asp = 1)
  points(x = BCall$lon, y = BCall$lat, cex = 0.1, pch = 4, col = 'red', asp = 1)
  points(x = XDall$lon, y = XDall$lat, cex = 0.1, pch = 4, col = 'blue', asp = 1)
  legend('topleft', legend = c('B. cylindrica', 'X. derbentina'), col = c('red', 'blue'), pch = 4, bty = 'n')
  
  ####multicollinearity porblem
  
  usdm::Variogram(roi)
  
    # pear <- raster::layerStats(x = roi, stat = 'pearson', na.rm = T)
  # save(dat, file = '/home/mikhail/Documents/slug_sdm/predictors_europe_slice.rda')
  # save(dat, file = '/home/mikhail/Documents/slug_sdm/predictors_europe_slice.rda')
  #custom data
  tb <- readxl::read_xlsx(path = '/home/mikhail/Documents/slug_sdm/invasive_snails_SDM.xlsx')
  # tb %>% ##View()
  tb %>% filter(species == 'B_cylindrica') -> BC_cus
  tb %>% filter(species == 'X_derbentina') -> XD_cus
  #GBIF data
  # #reading in occurences
  #next data retrieval
  BC_gbif <- gbif("Brephulopsis", species="cylindrica", ext = e, geo = T )
  XD_gbif <- gbif("Xeropicta", species="derbentina", ext = e, geo = T )
  
  
  
  #raster sonversion
  
  BC$species <- 1
  sp_CL <- dismo_CL[,c('lon', 'lat', 'species')]
  head(sp_CL)
  coordinates(sp_CL) <- ~lon + lat
  class(sp_CL)
  head(sp_CL)
  
  sp_BC <- tibble(cbind(longitude = as.numeric(BC$longitude),
                        latitude = as.numeric(BC$latitude),
                        species = BC$species
                        ))
  head(sp_BC)
  coordinates(sp_BC) <- ~longitude + latitude
  class(sp_CL)
  head(sp_CL)
  
  BCspdf <- SpatialPointsDataFrame(coords = )
  
  
  sp::SpatialMultiPointsDataFrame(coords = )
  coordinates(BC_cus) <- ~"latitute" + "longitude"
  
  vif <-usdm::vifcor(roi[[1:20]])
  vif
                  
  niche(x = scale(roi), h = p1_BC, n = names(roi)[c(1,4)])
  niche(x = scale(roi), h = p1_XD, n = names(roi)[c(1,4)])
  
  
  df <- sapply(roi, function(x){})
  roi %>% values() ->df
  df[,1:20] %>% cor(method = 'pear', use = 'compl') ->cor
  colnames(cor)
  corrplot::corrplot(cor)
  
  FactoMineR::PCA(df[,1:20]) ->pcad
  
  
  library(sen2r)
  # 
  
  r1 <- raster('/home/mikhail/Documents/travelog/LC08_L1TP_178029_20200516_20200516_01_RT(1)/LC08_L1TP_178029_20200516_20200516_01_RT_B4.TIF')
  r2 <-   raster('/home/mikhail/Documents/travelog/LC08_L1TP_178029_20200516_20200516_01_RT(1)/LC08_L1TP_178029_20200516_20200516_01_RT_B3.TIF')
  r3 <-   raster('/home/mikhail/Documents/travelog/LC08_L1TP_178029_20200516_20200516_01_RT(1)/LC08_L1TP_178029_20200516_20200516_01_RT_B2.TIF')
  rgb = raster::stack(r1, r2, r3)
  
  
  crs(rgb); roi %>% crs
  extent(roi); extent(rgb)
  class(rgb); class(roi)
  raster::projectRaster(rgb, crs, crs = crs(roi))
  
  
  crop(rgb, roi)
  
#distant probing parrt
rgb %>% extent

Xs <- c(30, 55)
Ys <- c(44, 55)


#modistools
library(MODISTools)
bands <- mt_bands(product = "MOD13Q1")
head(bands)

arcachon_lc <- mt_subset(product = "MCD12Q1",
                         #lat = 44.656286,
                         #lon =  -1.174748,
                         band = "LC_Type1",
                         start = "2004-01-01",
                         end = "2004-01-01",
                         #km_lr = 200,
                         #km_ab = 20,
                         site_name = "europe",
                         internal = TRUE,
                         progress = FALSE)


LC_r <- mt_to_raster(df = arcachon_lc, reproject = TRUE)
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"):
#> Discarded ellps unknown in CRS definition: +proj=sinu +lon_0=0 +x_0=0 +y_0=0
#> +R=6371007.181 +units=m +no_defs
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"): Discarded
#> datum unknown in CRS definition
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"):
#> Discarded ellps unknown in CRS definition: +proj=sinu +lon_0=0 +x_0=0 +y_0=0
#> +R=6371007.181 +units=m +no_defs
#> Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"): Discarded
#> datum unknown in CRS definition

# plot the raster data as a map
raster::plot(LC_r)


df <- data.frame("site_name" = paste("test",1:2), stringsAsFactors = FALSE)
df$lat <- 40
df$lon <- -110
subsets <- mt_batch_subset(df = df,
                           product = "MOD11A2",
                           band = "LST_Day_1km",
                           internal = TRUE,
                           start = "2004-01-01",
                           end = "2004-02-28",
                           out_dir = "~")
LC_r <- mt_to_raster(df = subsets, reproject = TRUE)


#modis

library(gdalUtils)
# Get a list of sds names
sds <- get_subdatasets('/home/mikhail/Documents/slug_sdm/MOD13C1.A2020209.006.2020226060035.hdf')
# Any sds can then be read directly using the raster function
ndvi <- raster(sds[1])
plot(r)

crop(ndvi, roi) ->ndvi_roi
ndvi_roi %>% plot(main = "MODIS NDVI")

#Sasha's files
setwd('/home/mikhail/Documents/slug_sdm/')
dir <- dir('/home/mikhail/Documents/slug_sdm/modis_txts/')
for (txt in dir){
  dir.create(paste0('/home/mikhail/Documents/slug_sdm/', 'retrieved_', txt))
    https <- readLines(paste0('/home/mikhail/Documents/slug_sdm/modis_txts/',txt))
  for (n in seq_along(https)){
    download.file(https[n],
              #destfile= paste0('/home/mikhail/Documents/modis_files/', i, '.hcf'),
              destfile= paste0('/home/mikhail/Documents/slug_sdm/', 'retrieved_', txt, '/', n, '.hdf'),
              method="curl")
  }
}

wget(https[n])
#read in hdfs and merge
library(gdalUtils)
# Get a list of sds names
sds <- get_subdatasets('full/path/filename.hdf')
# Any sds can then be read directly using the raster function
r <- raster(sds[1])


retrieved_MCD12Q1_2018 <- dir('/home/mikhail/Documents/slug_sdm/retrieved_MCD12Q1-2018.txt', full.names = T)
get_subdatasets(retrieved_MCD12Q1_2018[1])
gdalinfo(retrieved_MCD12Q1_2018[1])



#autumn att

library(MODIS)
ti <- getTile()
ti
 ###View(getProduct())
 #MOD13A3	
library(MODIS)
Austria <- extent(9.2, 17.47, 46.12, 49.3)
b2 <- getHdf(product = "MOD13A3", 
             begin = "2020001",
             #end = "2020002",
             extent = Austria)
  