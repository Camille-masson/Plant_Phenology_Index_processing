#
#                    Double Logistic Models for phenological data fitting
#               Estimation of parameters using Copernicus 10m PPI-based phenology
#                                    Philippe Choler
#                               last update : 15.05.2025
#

# Bibliography ----

# 2018.	Jonsson, P., Cai, Z., Melaas, E., Friedl, M. A. & Eklundh, L. A Method for Robust Estimation of Vegetation Seasonality from Landsat and Sentinel-2 Time Series Data. Remote Sensing 10, doi:10.3390/rs10040635 (2018). -> method implemented for Copernicus (producer is VITO, https://remotesensing.vito.be/)

# 2012. Elmore & al. Global Change Biol. 18, 656-674, doi:10.1111/j.1365-2486.2011.02521.x (2012). -> develop a 7-parameter DL that better simulates the senescence. The 6 param model is the one retained in Copernicus (data cubes produced by VITO)

# 2007. Julien, Y. & Sobrino, J. A. Changes in the global vegetal cover through a phenological analysis of GIMMS data. 2007 International Workshop on the Analysis of Multi-Temporal Remote Sensing Images, 228-232 290 (2007). -> used the DL function to derive global land-surface phenology trends from the GIMMS database. These studies have proven the effectiveness and practicability of the DL function in filtering NDVI TS.

# 2007. Fisher, J. I. & Mustard, J. F. Cross-scalar satellite phenology from ground, Landsat, and MODIS data. Remote Sens. Environ. 109, 261-273 (2007). -> used the DL function (from Zhang)  to extract MODIS-based seasonal NDVI in USA forests and compare with ground-based long-term records from Harvard Forest (Massachusetts) and Hubbard Brook (New Hampshire) 

# 2007. Richardson, A. D. et al. Use of digital webcam images to track spring green-up in a deciduous broadleaf forest. Oecologia 152, 323-334, doi:10.1007/s00442-006-0657-z (2007). -> use a Four-parameter model (simple DL) of green-up and compare with time-lapse

# 2006. Beck, P. S. A., C. Atzberger, K. A. Hogda, B. Johansen, and A. K. Skidmore. 2006. Improved monitoring of vegetation dynamics at very high latitudes: A new method using MODIS NDVI. Remote Sensing of Environment 100:321-334. -> show that the DL function out-performed the Fourier series and the asymmetric Gaussian function in describing NDVI time series.

# 2003. Zhang, X. Y. et al. Monitoring vegetation phenology using MODIS. Remote Sens. Environ. 84, 471-475 (2003). -> show that one may capture the seasonal time course of VI using  piecewise logistic function of time. show aloso the analytic formula of curvature


rm(list=ls())

library(terra)
library(hdar)      # to access Copernicus data using wekeo services
library(jsonlite)  # to modify the query with json

# tempdir()
# [1] "C:\Users\XYZ~1\AppData\Local\Temp\Rtmp86bEoJ\Rtxt32dcef24de2"
# dir.create(tempdir())

# I. 6-parameters Double Logistic (DL) Function ----
  # I.A. Define functions ----
# Double Logistic Function -> DL.f
DL.f <- function(DOY,NDVImin,NDVImax,ONSET,OFFSET,GROWTH,SENESC){
  NDVI = NDVImin + (NDVImax - NDVImin)*(1/(1+exp(-GROWTH*(DOY-ONSET)))+1/(1+exp(SENESC*(DOY-OFFSET)))-1)
  # NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}

# Growth part of the DL Function -> DLG.f
DLG.f  <-  function(DOY,NDVImin,NDVImax,ONSET,GROWTH){
  NDVI = NDVImin + (NDVImax - NDVImin)*(1/(1+exp(-GROWTH*(DOY-ONSET))))
  # NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}

# Inverse of the growth part -> DIG.f
DIG.f <- function(NDVI,NDVImin,NDVImax,ONSET,GROWTH){
  DOY = (GROWTH*ONSET-log((NDVImax-NDVImin)/(NDVI - NDVImin)-1))/GROWTH
  # NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(DOY)
}
# DLG.f(150,NDVImin,NDVImax,ONSET,GROWTH)
# DIG.f(0.5,NDVImin,NDVImax,ONSET,GROWTH)


# Senescence part of the DL function -> DLS.f
DLS.f  <-  function(DOY,NDVImin,NDVImax,OFFSET,SENESC){
  NDVI = NDVImin + (NDVImax - NDVImin)*(1/(1+exp(SENESC*(DOY-OFFSET))))
  return(NDVI)
}

# Inverse of the senescence part -> DIS.f
DIS.f <- function(NDVI,NDVImin,NDVImax,OFFSET,SENESC){
  DOY = (SENESC*OFFSET+log((NDVImax-NDVImin)/(NDVI - NDVImin)-1))/SENESC
  # NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(DOY)
}
# DLS.f(300,NDVImin,NDVImax,OFFSET,SENESC)
# DIS.f(0.2582529,NDVImin,NDVImax,OFFSET,SENESC)

# IRG formule analytique de la dérivée de la DL -> DLd.f
DLd.f <- function(DOY,NDVImin,NDVImax,ONSET,OFFSET,GROWTH,SENESC){
  IRG = (NDVImax - NDVImin)*(GROWTH*exp(-GROWTH*(DOY-ONSET))/(1+exp(-GROWTH*(DOY-ONSET)))^2 - SENESC*exp(-SENESC*(DOY-OFFSET))/(1+exp(-SENESC*(DOY-OFFSET)))^2)
  # NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # if rescaling
  return(IRG)
}

# IRG formule analytique de la dérivée pour la phase de growth -> DLGd.f
DLGd.f <-  function(DOY,NDVImin,NDVImax,ONSET,GROWTH){
  NDVId = (NDVImax - NDVImin)*GROWTH*exp(-GROWTH*(DOY-ONSET))/(1+exp(-GROWTH*(DOY-ONSET)))^2
  # NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # if rescaling
  return(NDVId)
}

# IRG formule analytique de la dérivée pour la phase de senescence -> DLSd.f
DLSd.f <-  function(DOY,NDVImin,NDVImax,OFFSET,SENESC){
  NDVId = -(NDVImax - NDVImin)*SENESC*exp(-SENESC*(DOY-OFFSET))/(1+exp(-SENESC*(DOY-OFFSET)))^2
  # NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVId)
}

  # I.B. Define a set of parameters ----
DOY.seq     = seq(1,365,1)     # Day of (non bissextile) year
NDVImin     = 0.1              # annual minimum NDVI
NDVImax     = 0.9              # annual maximum NDVI
NDVIamp     = NDVImax - NDVImin # amplitude
OFFSET      = 280              # date of the autumn inflection point
ONSET       = 150              # date of the spring inflection point
GROWTH      = 0.15             # rate of NDVI increase at ONSET
SENESC      = 0.07             # rate of NDVI decrease at OFFSET

SOSV        = NDVImin + 0.25*NDVIamp # VI at SOSD, i.e. 25% of NDVIamp is reached
SOSD        = DIG.f(SOSV,NDVImin,NDVImax,ONSET,GROWTH) # Start of Season Date (Copernicus)
LSLOP       = DLGd.f(SOSD,NDVImin,NDVImax,ONSET,GROWTH)

EOSV        = NDVImin + 0.15*NDVIamp # VI at EOSD, i.e. 15% of NDVIamp is reached
EOSD        = DIS.f(EOSV,NDVImin,NDVImax,OFFSET,SENESC) # End of Season sensu Copernicus - 
RSLOP       = DLSd.f(EOSD,NDVImin,NDVImax,OFFSET,SENESC)

# 3. Simulate yearly NDVI time series
NDVIyts  <- DL.f(DOY=DOY.seq,NDVImin=NDVImin,NDVImax=NDVImax,ONSET=ONSET,OFFSET=OFFSET,GROWTH=GROWTH,SENESC=SENESC)

  # I.C. Reference figure ----
graphics.off()
windows(10,8)
par(mfrow=c(2,1), oma=c(0,0,0,0),mar=c(1.5,3,0.5,1),mgp=c(1.5,0.5,0))

  # Panel a - seasonal trajectory of VI
plot(DOY.seq,NDVIyts,type="l",xlab="", ylab="Vegetation Indice",col="yellow",lwd=10,ylim=c(0,1))
abline(h=c(NDVImin,NDVImax),lty=c(2,2))
# for(i in 1:length(ONSET.seq)) lines(DOY.seq,NDVIs[,i],col=rainbow(length(ONSET.seq))[i])

# SOS copernicus 25% du seasonal amplitude durant le greenup
abline(v=SOSD,lty=2)
points(SOSD,SOSV,bg="green",pch=21,cex=2)
# EOS copernicus 15% du seasonal amplitude durant le greendown
abline(v=EOSD,lty=2)
points(EOSD,EOSV,bg="blue",pch=21,cex=2)

# Left Inflection point DLL
abline(v=ONSET,lty=2)
points(ONSET,DLG.f(ONSET,NDVImin,NDVImax,ONSET,GROWTH),bg="green",pch=23,cex=2)
# Right Inflexion point DLL
abline(v=OFFSET,lty=2)
points(OFFSET,DLS.f(OFFSET,NDVImin,NDVImax,OFFSET,SENESC),bg="blue",pch=23,cex=2)

lines(DOY.seq,DLG.f(DOY.seq,NDVImin,NDVImax,ONSET,GROWTH),col="green",lwd=2,lty=1)
lines(DOY.seq,DLS.f(DOY.seq,NDVImin,NDVImax,OFFSET,SENESC),col="blue",lwd=2,lty=1)

  # Panel b - curvature = IGP 
plot(DOY.seq,c(NA,diff(NDVIyts)),col="yellow",type="l",xlab="Day Of Year", ylab="Diff VI",lwd=10,ylim=c(-0.03,0.03))
# lines(DOY.seq,DLd.f(DOY.seq,NDVImin,NDVImax,ONSET.seq[1],OFFSET,GROWTH,SENESC),col="gray",type="l",xlab="Day Of Year", ylab="Diff VI",lwd=10,ylim=c(-0.03,0.03))
abline(h=0,lty=1)

# Left Inflection point DLL
abline(v=c(ONSET),lty=2)
points(ONSET,DLGd.f(ONSET,NDVImin,NDVImax,ONSET,GROWTH),bg="green",pch=23,cex=2)
# Right Inflection point DLL
abline(v=OFFSET,lty=2)
points(OFFSET,DLSd.f(OFFSET,NDVImin,NDVImax,OFFSET,SENESC),bg="blue",pch=23,cex=2)
# SOS copernicus 25% du seasonal amplitude durant le greenup
abline(v=SOSD,lty=2)
points(SOSD,LSLOP,bg="green",pch=21,cex=2)
# EOS copernicus 15% du seasonal amplitude durant le greendown
abline(v=EOSD,lty=2)
points(EOSD,RSLOP,bg="blue",pch=21,cex=2)

lines(DOY.seq,DLGd.f(DOY.seq,NDVImin,NDVImax,ONSET,GROWTH),col="green",lwd=2,lty=1)
lines(DOY.seq,DLSd.f(DOY.seq,NDVImin,NDVImax,OFFSET,SENESC),col="blue",lwd=2,lty=1)

# II. Estimate DL parametes using Copernicus 20m res, PPI-based land phenology  ----
  # II.A. Downloading CLMS_HRVPP_VPP products ----
# A tutorial can be found here :
# https://help.wekeo.eu/en/articles/7035318-how-to-use-the-hdar-package-for-accessing-the-wekeo-hda-api-in-r

# wekeo credentials
username <- "pcholer"
password <- "Cardamine@2021"
client   <- Client$new(username, password, save_credentials = TRUE)

query_template      <- client$generate_query_template("EO:EEA:DAT:CLMS_HRVPP_VPP")
query_templateINI   <- fromJSON(query_template, flatten=F)
query_templateINI



getwd
output_directory










# Copy an example of API from wekeo
query_BEL <- paste0('{ 
    "dataset_id":"EO:EEA:DAT:CLMS_HRVPP_VPP-LAEA",
    "httpAccept":"application%2Fgeo%2Bjson",
    "recordSchema":"geojson",
    "productType": "EOSD",
    "productGroupId": "s1",
    "resolution":"10",
    "bbox": [
    6.670773080321873,
    44.2089367691059,
    6.849267152601807,
    44.296048517023806],
    "startdate": "2022-01-01T00:00:00.000Z",
    "enddate"  : "2023-12-31T23:59:59.999Z",
    "itemsPerPage": 200,
    "startIndex": 0
}')




#### 0. LIBRARY and PARAMETERS ####
#---------------------------------#

site <- "Rouanette"









#### 1. Downloading Plant Phenology Index (PPI) ####
#--------------------------------------------------#

if (TRUE) {
  
  
  # Définition des chemins 
  root_dir <- getwd()  # Récupère le chemin du projet automatiquement
  input_case <- file.path(root_dir, "input")
  if (!dir.exists(input_case)) {
    dir.create(input_case, recursive = TRUE)
  }
  
  
  data_case <- file.path(input_case, paste0("data_",site))
  if (!dir.exists(data_case)) {
    dir.create(data_case, recursive = TRUE)
  }
  
  
  output_directory <- file.path(data_case, "downloads")
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  
  # download 2020:2023 products
  VAR <- c("AMPL","EOSD","EOSV","LSLOPE","MAXD","MAXV","MINV","RSLOPE","SOSD","SOSV")
  for (i in 1:length(VAR)){
    print(i)
    query2_BEL <- paste0("{ \n    \"dataset_id\":\"EO:EEA:DAT:CLMS_HRVPP_VPP-LAEA\",\n    \"httpAccept\":\"application%2Fgeo%2Bjson\",\n    \"recordSchema\":\"geojson\",\n    \"productType\": \"",VAR[i],"\",\n    \"productGroupId\": \"s1\",\n    \"resolution\":\"10\",\n    \"bbox\":[6.207312871270429,
    44.57325963145546,
    6.4282245342827595,
    44.71470256740122],\n    \"startdate\": \"2023-01-01T00:00:00.000Z\",\n    \"enddate\": \"2023-12-31T23:59:59.999Z\",\n    \"itemsPerPage\": 200,\n    \"startIndex\": 0\n}")
    matches <- client$search(query2_BEL)
    # Assuming 'matches' is an instance of SearchResults obtained from the search 
    matches$download(output_directory,force=F,prompt=F)
  }
  








}






#### 2. Estimate the missing parameters  (ONSET, GROWTH, OFFSET, SENESC) ####
#---------------------------------------------------------------------------#

if (TRUE){
  
  # 0. Paramètres généraux
  site <- "Rouanette"
  YEAR <- 2023
  VAR  <- c("AMPL","EOSD","EOSV","LSLOPE","MAXD","MAXV","MINV","RSLOPE","SOSD","SOSV")
  
  # 1. Création des répertoires
  root_dir         <- getwd()
  input_case       <- file.path(root_dir, "input")
  data_case        <- file.path(input_case, paste0("data_", site))
  extent_case      <- file.path(data_case, "extent")
  output_directory <- file.path(data_case, "downloads")
  
  for (d in c(input_case, data_case, extent_case, output_directory)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  # 2. Chargement du masque L93 et création de la grille REF10 (L93, 10 m)
  mask_shp <- file.path(extent_case, paste0("Masque_", site, ".shp"))
  mask_vec <- vect(mask_shp)        # doit être en EPSG:2154
  AOIext   <- ext(mask_vec)
  REF10    <- rast(
    xmin = 1000 * floor(AOIext$xmin / 1000),
    xmax = 1000 * ceiling(AOIext$xmax / 1000),
    ymin = 1000 * floor(AOIext$ymin / 1000),
    ymax = 1000 * ceiling(AOIext$ymax / 1000),
    crs  = crs(mask_vec),           # EPSG:2154
    res  = 10
  )
  
  # 3. Pour chaque variable VAR, reprojeter → mosaïque (si besoin) → crop & mask → rescale si PPI → sauver
  # 3. Assemble & rééchelle comme dans le "script de base II.C"
  for (prod in VAR) {
    
    # a) Liste des tuiles brutes (LAEA) pour l’année et le produit
    files <- list.files(
      output_directory,
      pattern    = paste0(YEAR, ".*", prod, ".*\\.tif$"),
      full.names = TRUE
    )
    if (length(files) == 0) {
      warning("Pas de fichier trouvé pour ", prod, " ", YEAR)
      next
    }
    
    # b) Lecture de toutes les tuiles en LAEA
    r_la_list <- lapply(files, terra::rast)
    
    # c) Mosaïque en LAEA (pas de reproj par tuile)
    sprc_tiles <- terra::sprc(r_la_list)
    mos_laea   <- terra::mosaic(sprc_tiles)
    
    # d) Projection de la mosaïque LAEA → REF10 (L93 10 m)
    mos_l93 <- terra::project(mos_laea, REF10)
    
    # e) Crop + mask en L93
    mos_crop <- terra::crop(mos_l93, mask_vec)
    mos_mask <- terra::mask(mos_crop, mask_vec)
    
    
    
    # g) Sauvegarde avec le nom <prod>_<site>_<YEAR>.tif
    out_file <- file.path(
      output_directory,
      paste0(prod, "_", site, "_", YEAR, ".tif")
    )
    terra::writeRaster(mos_mask, out_file, overwrite = TRUE)
    message("→ Écrit : ", out_file)
  }
  
  # example for
  
  # on liste tous les 10 rasters <VAR>_<site>_<YEAR>.tif
  LF <- list.files(output_directory,
                   pattern    = paste0("_", site, "_", YEAR, "\\.tif$"),
                   full.names = TRUE)
  
  # une seule itération (grep("") renvoie tous les fichiers)
  TILES <- ""
  for (j in TILES){
    print(j)
    FILES <- LF[grep(j,LF)] # should return 10 files
    
    lag    <- 1000*as.numeric(substr(YEAR,3,4))
    AMPL   <- terra::rast(FILES[grep("AMPL",FILES)])
    MINV   <- terra::rast(FILES[grep("MINV",FILES)])
    SOSV   <- terra::rast(FILES[grep("SOSV",FILES)])
    EOSV   <- terra::rast(FILES[grep("EOSV",FILES)])
    LSLOPE <- terra::rast(FILES[grep("LSLOPE",FILES)])
    RSLOPE <- terra::rast(FILES[grep("RSLOPE",FILES)])
    SOSD   <- terra::rast(FILES[grep("SOSD",FILES)])-lag ; SOSD[SOSD<0 | SOSD>365] <- NA
    EOSD   <- terra::rast(FILES[grep("EOSD",FILES)])-lag ; EOSD[EOSD<0 | EOSD>365] <- NA
    
    # let A = exp(-GROWTH*(SOSD-ONSET))
    # Compute A using DLG.f(SOSV)
    A       <- (AMPL)/(SOSV-MINV)-1
    A[A<=0] <- NA  
    
    # Compute GROWTH using DLGd.f(SOSV,LSLOPE)
    GROWTH              <- round(LSLOPE/(A/(1+A)^2*(AMPL)),4)
    GROWTH[GROWTH[]<0]  <-NA
    GROWTH[GROWTH[]>0.5]<-NA
    hist(GROWTH[],breaks=seq(0,0.5,0.01),xlim=c(0,0.5))
    terra::writeRaster(GROWTH,gsub("AMPL","GROWTH",FILES[1]),overwrite=T)
    
    # Compute ONSET using DIG.f(SOSD)
    ONSET <- round((log(A) + GROWTH*SOSD)/GROWTH)
    ONSET[(ONSET-SOSD)<0] <- NA; ONSET[ONSET>365]<-NA
    hist(ONSET[],breaks=seq(0,365,1),xlim=c(0,365))
    terra::writeRaster(ONSET,gsub("AMPL","ONSET",FILES[1]),overwrite=T)
    
    # let B = exp(SENESC*(EOSD-OFFSET))
    # Compute B using DLS.f(EOSV) # to be done    
    B      <- (AMPL)/(EOSV-MINV)-1
    B[B<=0] <- NA  
    
    # Compute SENESC using DLSd.f(EOSV,RSLOPE)
    SENESC <- round((1+1/B)^2*B*RSLOPE/AMPL,4)
    SENESC[SENESC[]<=0]  <-NA
    SENESC[SENESC[]>0.5] <-NA
    hist(SENESC[],breaks=seq(0.0,0.5,0.01),xlim=c(0,0.5))
    terra::writeRaster(SENESC,gsub("AMPL","SENESC",FILES[1]),overwrite=T)
    
    # Compute OFFSET using DIS.f(EOSD)
    OFFSET <- round(EOSD - log(B)/SENESC)
    OFFSET[OFFSET>365] <- NA
    OFFSET[OFFSET<0]   <- NA
    hist(OFFSET[],breaks=seq(0,365,1),xlim=c(0,365))
    terra::writeRaster(OFFSET,gsub("AMPL","OFFSET",FILES[1]),overwrite=T)
  }
  
  
  
  
  







}




#### 3. Calcul of double logistic in PPI ####
#-------------------------------------------#

if (TRUE){
  library(terra)
  
  # 1. Paramètres
  site      <- "Rouanette"
  YEAR      <- 2023
  output_directory <- file.path(getwd(), "input", paste0("data_", site), "downloads")
  
  # 2. On importe directement les rasters dérivés produits en partie 2
  params <- c("MINV","MAXV","ONSET","OFFSET","GROWTH","SENESC")
  
  # Construire une liste nommée de SpatRaster
  rasters <- setNames(
    lapply(params, function(p) {
      fp <- file.path(output_directory, paste0(p, "_", site, "_", YEAR, ".tif"))
      if (!file.exists(fp)) stop("Fichier introuvable : ", fp)
      rast(fp)
    }),
    params
  )
  
  NDVImin <- rasters[["MINV"]]
  NDVImax <- rasters[["MAXV"]]
  ONSET   <- rasters[["ONSET"]]
  OFFSET  <- rasters[["OFFSET"]]
  GROWTH  <- rasters[["GROWTH"]]
  SENESC  <- rasters[["SENESC"]]
  
  # 3. Définition du modèle DL.f
  DL.f <- function(DOY, NDVImin, NDVImax, ONSET, OFFSET, GROWTH, SENESC) {
    NDVImin + (NDVImax - NDVImin) * (
      1/(1 + exp(-GROWTH * (DOY - ONSET))) +
        1/(1 + exp( SENESC * (DOY - OFFSET))) -
        1
    )
  }
  
  # 4. Application journalière & empilement
  DOY_range <- 0:365  # du 1er mai (121) au 30 nov (334)
  
  ndvi_list <- lapply(DOY_range, function(doy) {
    lapp(
      c(NDVImin, NDVImax, ONSET, OFFSET, GROWTH, SENESC),
      fun      = DL.f,
      filename = "",      # pas d’écriture immédiate
      DOY      = doy      # passe 'doy' au premier argument de DL.f
    )
  })
  
  ndvi_stack <- rast(ndvi_list)
  names(ndvi_stack) <- paste0("NDVI_DOY", DOY_range)
  
  # 5. Sauvegarde finale
  writeRaster(
    ndvi_stack,
    filename  = file.path(output_directory, paste0("NDVI_season_", YEAR, ".tif")),
    overwrite = TRUE
  )
  
  message("▶ NDVI saisonnier ", YEAR, " (DOY ", min(DOY_range),
          "–", max(DOY_range), ") empilé et sauvé dans :\n", 
          file.path(output_directory, paste0("NDVI_season_", YEAR, ".tif")))
  




}









## NDVI max 

library(terra)

# 1. Chemin du stack NDVI saisonnier
site      <- "Rouanette"
YEAR      <- 2023
output_directory <- file.path(getwd(), "input", paste0("data_", site), "downloads")
season_stack_fp  <- file.path(output_directory, paste0("NDVI_season_", YEAR, ".tif"))

# 2. Lecture du SpatRaster multi-couches
season_stack <- rast(season_stack_fp)

# 3. Calcul du NDVI max par pixel
# Méthode 1 : avec terra::app()
ndvi_max <- app(season_stack,
                fun = function(...) max(..., na.rm = TRUE),
                filename = "",    # pas d’écriture immédiate
                overwrite = FALSE)

# OU plus simplement :
# ndvi_max <- terra::max(season_stack, na.rm = TRUE)

# 4. Sauvegarde
out_fp <- file.path(output_directory,
                    paste0("NDVI_max_", site, "_", YEAR, ".tif"))
writeRaster(ndvi_max, out_fp, overwrite = TRUE)

message("⇒ Raster NDVI max saisonnier généré et sauvé sous :\n", out_fp)








#### 4. Calcul et empilement de l’IRG (Instantaneous Rate of Green‐up) ####
#------------------------------------------------------------------------#


if (TRUE){
  
  library(terra)
  
  # 1. Paramètres
  site      <- "Rouanette"
  YEAR      <- 2023
  output_directory <- file.path(getwd(), "input", paste0("data_", site), "downloads")
  
  # 2. Chargement des rasters dérivés (déjà produits en partie 2)
  params <- c("MINV","MAXV","ONSET","OFFSET","GROWTH","SENESC")
  rasters <- setNames(
    lapply(params, function(p) {
      fp <- file.path(output_directory, paste0(p, "_", site, "_", YEAR, ".tif"))
      if (!file.exists(fp)) stop("Fichier introuvable : ", fp)
      rast(fp)
    }),
    params
  )
  NDVImin <- rasters[["MINV"]]
  NDVImax <- rasters[["MAXV"]]
  ONSET   <- rasters[["ONSET"]]
  OFFSET  <- rasters[["OFFSET"]]
  GROWTH  <- rasters[["GROWTH"]]
  SENESC  <- rasters[["SENESC"]]
  
  # 3. Définition de la dérivée de la double logistique (IRG)
  DLd.f <- function(DOY, NDVImin, NDVImax, ONSET, OFFSET, GROWTH, SENESC) {
    (NDVImax - NDVImin) * (
      GROWTH * exp(-GROWTH * (DOY - ONSET)) / (1 + exp(-GROWTH * (DOY - ONSET)))^2 -
        SENESC * exp(-SENESC * (DOY - OFFSET)) / (1 + exp(-SENESC * (DOY - OFFSET)))^2
    )
  }
  
  # 4. Empilement IRG sur la saison (DOY 121–334)
  DOY_range <- 121:334
  irg_list  <- lapply(DOY_range, function(doy) {
    lapp(
      c(NDVImin, NDVImax, ONSET, OFFSET, GROWTH, SENESC),
      fun      = DLd.f,
      filename = "",    # pas d’écriture immédiate
      DOY      = doy
    )
  })
  
  irg_stack <- rast(irg_list)
  names(irg_stack) <- paste0("IRG_DOY", DOY_range)
  
  # 5. Sauvegarde de l’IRG saisonnier
  out_season_irg <- file.path(
    output_directory,
    paste0("IRG_season_", site, "_", YEAR, ".tif")
  )
  writeRaster(irg_stack, out_season_irg, overwrite = TRUE)
  message("▶ IRG saisonnier empilé et sauvé sous :\n", out_season_irg)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #’IRG max et son DOY
  
  IRG_max <- app(irg_stack, fun = function(...) max(..., na.rm=TRUE))
   vals <- IRG_max[]; vals <- vals[!is.na(vals) & is.finite(vals)]
   hist(vals, breaks=seq(floor(min(vals)*100)/100, ceiling(max(vals)*100)/100, 0.01),
        xlim=c(floor(min(vals)*100)/100, ceiling(max(vals)*100)/100),
        main="Distribution de l'IRG max (2023)", xlab="IRG max", ylab="Nombre de pixels")
   
   out_max <- file.path(output_directory, paste0("IRG_max_", site, "_", YEAR, ".tif"))
   writeRaster(IRG_max, out_max, overwrite=TRUE)
   message("▶ IRG max sauvé sous :\n", out_max)
   
   
   
   
   
   
   
  
  
   
   
   
   
   
   
   # ─── repères déjà prêts en DOY (MAXD et EOSD lag-corrigés) ───────────────────
   MAXD   <- rast(file.path(output_directory, paste0("MAXD_"  , site, "_", YEAR, ".tif")))
   OFFSET <- rast(file.path(output_directory, paste0("OFFSET_", site, "_", YEAR, ".tif")))
   EOSD   <- rast(file.path(output_directory, paste0("EOSD_"  , site, "_", YEAR, ".tif")))
   
   lag <- 1000 * as.numeric(substr(YEAR, 3, 4))      # 2023 → 23000
   MAXD <- MAXD - lag ; MAXD[MAXD < 1 | MAXD > 365] <- NA
   EOSD <- EOSD - lag ; EOSD[EOSD < 1 | EOSD > 365] <- NA
   
   OM10 <- OFFSET - 10                                # plateau = MAXD+1 … OFFSET−10
   
   print(MAXD)
   print(OFFSET)
   print(EOSD)
   print(OM10)
   
   
   # ─── phases journalières : 1 pousse / 2 plateau / 3 dépérissement / 4 sénescence
   DOY_range   <- 121:334
   mask_valid  <- !is.na(MAXD) & !is.na(OM10) & !is.na(EOSD)  # pixels où tous les repères existent
   template    <- MAXD                                        # support vierge
   
   phase_stack <- rast(lapply(DOY_range, function(d){
     
     # un raster constant rempli du DOY courant
     d_rast <- setValues(template, d)
     
     # phase provisoire (1-4) sans gérer les NA
     phase <- ifel(d_rast <= MAXD,                 1,
                   ifel(d_rast <= OM10,                  2,
                        ifel(d_rast <= EOSD,                  3, 4)))
     
     # on remet NA là où au moins un repère manquait
     phase <- mask(phase, mask_valid, maskvalue = 0, updatevalue = NA)
     
     phase
   }))
   
   names(phase_stack) <- paste0("Phase_DOY", DOY_range)
   
  
   
   # ─── sauvegarde ──────────────────────────────────────────────────────────────
   out_fp <- file.path(output_directory,
                       paste0("PhenologyPhase_", site, "_", YEAR, ".tif"))
   writeRaster(
     phase_stack,
     out_fp,
     overwrite = TRUE,
     datatype  = "INT1U"   # 0-255, parfait pour 1-4
   )
   
   message("✓ Phases phénologiques écrites : ", out_fp,
           "\n   1=pousse, 2=plateau, 3=déperissement, 4=sénescence")
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   






}























library(terra)
## rasters déjà chargés ------------------------------------------------------
MAXD   <- rast(file.path(output_directory, paste0("MAXD_"  , site, "_", YEAR, ".tif")))
OFFSET <- rast(file.path(output_directory, paste0("OFFSET_", site, "_", YEAR, ".tif")))
EOSD   <- rast(file.path(output_directory, paste0("EOSD_"  , site, "_", YEAR, ".tif")))

## même traitement “DOY” que pour SOSD / EOSD dans votre script --------------
lag <- 1000 * as.numeric(substr(YEAR, 3, 4))   # ex. 2023  →  23000

MAXD <- MAXD - lag ; MAXD[MAXD < 0 | MAXD > 365] <- NA
EOSD <- EOSD - lag ; EOSD[EOSD < 0 | EOSD > 365] <- NA


print(MAXD)
print(OFFSET)
print(EOSD)


DOY_range <- 0:365                         # 1er mai → 30 nov.
offset_m10 <- OFFSET - 10                    # seuil OFFSET − 10 j

# Une couche « phase » par jour
phase_list <- lapply(DOY_range, function(doy) {
  lapp(c(MAXD, offset_m10, EOSD), fun = function(maxd, om10, eosd) {
    ifelse(doy <= maxd,                1,   # pousse
           ifelse(doy <= om10,              2,   # plateau
                  ifelse(doy <= eosd,            3,   # déperissement
                         4)))     # sénescence
  })
})

phase_stack <- rast(phase_list)
names(phase_stack) <- paste0("Phase_DOY", DOY_range)

# sauvegarde
out_fp <- file.path(output_directory,
                    paste0("PhenologyPhase_", site, "_", YEAR, ".tif"))
writeRaster(phase_stack, out_fp, overwrite = TRUE)

message(
  "▶ Raster des phases phénologiques écrit : ", out_fp, "\n",
  "   1 = pousse, 2 = plateau, 3 = déperissement, 4 = sénescence"
)








# rasters repères déjà chargés ----------------------------------------------
MAXD   <- rast(file.path(output_directory, paste0("MAXD_"  , site, "_", YEAR, ".tif")))
OFFSET <- rast(file.path(output_directory, paste0("OFFSET_", site, "_", YEAR, ".tif")))
EOSD   <- rast(file.path(output_directory, paste0("EOSD_"  , site, "_", YEAR, ".tif")))


EOSD   <- terra::rast(FILES[grep("EOSD",FILES)])-lag ; EOSD[EOSD<0 | EOSD>365] <- NA

EOSD   <- terra::rast(FILES[grep("EOSD",FILES)])-lag ; EOSD[EOSD<0 | EOSD>365] <- NA





print(MAXD)
print(OFFSET)
print(EOSD)





om10   <- OFFSET - 10                        # seuil OFFSET – 10 j

DOY_range <- 121:334                         # 1er mai → 30 nov.

# fonction qui renvoie la classe (1-4) ou NA si un repère manque --------------
phase_fun <- function(maxd, om10, eosd, doy) {
  out <- rep(NA_integer_, length(maxd))
  ok  <- is.finite(maxd) & is.finite(om10) & is.finite(eosd)
  out[ ok & doy <= maxd             ] <- 1   # pousse
  out[ ok & doy > maxd & doy <= om10] <- 2   # plateau
  out[ ok & doy > om10 & doy <= eosd] <- 3   # déperissement
  out[ ok & doy > eosd              ] <- 4   # sénescence
  out
}

# empilement -----------------------------------------------------------------
phase_stack <- rast(lapply(DOY_range, function(doy) {
  lapp(c(MAXD, om10, EOSD),
       fun      = phase_fun,
       dtype    = "INT1U",      # écrit directement un entier 1 octet (0–255)
       doy      = doy)
}))

names(phase_stack) <- paste0("Phase_DOY", DOY_range)

# sauvegarde -----------------------------------------------------------------
out_fp <- file.path(output_directory,
                    paste0("PhenologyPhase_", site, "_", YEAR, ".tif"))
writeRaster(phase_stack, out_fp, overwrite = TRUE)

message(
  "▶ Raster multi-couches des phases phénologiques écrit :\n", out_fp,
  "\n     1 = pousse, 2 = plateau, 3 = déperissement, 4 = sénescence"
)


















library(terra)

# rasters repères déjà prêts (0–365)
MAXD   <- rast(file.path(output_directory, paste0("MAXD_"  , site, "_", YEAR, ".tif")))
OFFSET <- rast(file.path(output_directory, paste0("OFFSET_", site, "_", YEAR, ".tif")))
EOSD   <- rast(file.path(output_directory, paste0("EOSD_"  , site, "_", YEAR, ".tif")))



lag <- 1000 * as.numeric(substr(YEAR, 3, 4))   # ex. 2023  →  23000

MAXD <- MAXD - lag ; MAXD[MAXD < 0 | MAXD > 365] <- NA
EOSD <- EOSD - lag ; EOSD[EOSD < 0 | EOSD > 365] <- NA


print(MAXD)
print(OFFSET)
print(EOSD)



om10 <- OFFSET - 10                         # seuil plateau
DOY_range <- 121:334                        # 1er mai → 30 nov.

phase_fun <- function(maxd, om10, eosd, d) {
  out <- rep(NA_integer_, length(maxd))
  ok  <- is.finite(maxd) & is.finite(om10) & is.finite(eosd)
  
  out[ ok & d <= maxd                      ] <- 1
  out[ ok & d >  maxd & d <= om10          ] <- 2
  out[ ok & d >  om10 & d <= eosd          ] <- 3
  out[ ok & d >  eosd                      ] <- 4
  out
}

phase_stack <- rast(lapply(DOY_range, function(d)
  lapp(c(MAXD, om10, EOSD), phase_fun, dtype = "INT1U", d = d)
))
names(phase_stack) <- paste0("Phase_DOY", DOY_range)

# sauvegarde
out_fp <- file.path(output_directory,
                    paste0("PhenologyPhase_", site, "_", YEAR, ".tif"))
writeRaster(phase_stack, out_fp, overwrite = TRUE)

message("✔ Phases enregistrées dans ", out_fp,
        "\n   1=pousse, 2=plateau, 3=déperissement, 4=sénescence")























































library(terra)

# 1. Chemin vers le dossier de sorties
output_directory <- file.path(getwd(), "input", paste0("data_", site), "downloads")

# 2. Lister les tif MAXV pour l'année 2022
maxv_files <- list.files(
  output_directory,
  pattern    = "2022.*MAXV.*\\.tif$",
  full.names = TRUE
)
if (length(maxv_files) == 0) stop("Aucun fichier MAXV trouvé.")

# 3. Lecture et reprojection
#    On prend la première tuile (ou toutes, si vous voulez mosaïquer par la suite)
r_laea <- rast(maxv_files[1])
# Forcer la CRS source si manquante
if (is.na(crs(r_laea))) crs(r_laea) <- "EPSG:3035"

# 4. Projeté sur EPSG:2154
r_l93 <- project(r_laea, "EPSG:2154")

# 5. Sauvegarde
out_file <- file.path(output_directory, "test.tif")
writeRaster(r_l93, out_file, overwrite = TRUE)

message("✓ Reprojection terminée : ", out_file)




























DL.f <- function(DOY,NDVImin,NDVImax,ONSET,OFFSET,GROWTH,SENESC){
  NDVI = NDVImin + (NDVImax - NDVImin)*(1/(1+exp(-GROWTH*(DOY-ONSET)))+1/(1+exp(SENESC*(DOY-OFFSET)))-1)
  # NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}




DL.f(DOY = )



























# IRG formule analytique de la dérivée de la DL -> DLd.f
DLd.f <- function(DOY,NDVImin,NDVImax,ONSET,OFFSET,GROWTH,SENESC){
  IRG = (NDVImax - NDVImin)*(GROWTH*exp(-GROWTH*(DOY-ONSET))/(1+exp(-GROWTH*(DOY-ONSET)))^2 - SENESC*exp(-SENESC*(DOY-OFFSET))/(1+exp(-SENESC*(DOY-OFFSET)))^2)
  # NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # if rescaling
  return(IRG)
}


















































# PRELIMINARY TRIALS 
# client$get_token()
# client$terms_and_conditions(term_id = 'all')
# all_datasets <- client$datasets()
# 
# filtered_datasets <- client$datasets("EO:EEA:DAT:CLMS_HRVPP_VPP")
# length(filtered_datasets)
# filtered_datasets[[1]] # LAEA projection
# filtered_datasets[[2]] # UTM projection

# transform the query in readable/ediable format
# myquery             <- fromJSON(query_template, flatten = FALSE)
# myquery$productType <- "MINV"  # choose parameter
# myquery$bbox        <- c(5.8, 45.05,6.29,45.53)    # choose parameter
# myquery$startdate   <- "2023-01-01T00:00:00"  # choose parameter
# myquery$enddate     <- "2023-12-31T00:00:00"  # choose parameter
# myquery$recordSchema<- "geojson"
# myquery$productionStatus<- "ARCHIVED"
# myquery$productGroupId <-"s1"
# 
# # Convert back to JSON format
# myqueryJSON         <- toJSON(myquery, auto_unbox = TRUE, digits = 17)
# matches <- client$search(myqueryJSON)

# SUPPLEMENTARY REQUEST
# Persistent snow area 20m resolution
query_PSA <- '{
    "dataset_id": "EO:CRYO:DAT:HRSI:PSA",
    "startdate" : "2016-09-01T00:00:00.000Z",
    "enddate"   : "2025-05-15T23:59:59.999Z",
    "itemsPerPage": 200,
    "startIndex": 0
}'










  # II.C. Assemble for Belledonne ----
AOIext <- terra::ext(terra::vect("~/PROJETSmc/RESALP/DATA_ANALYSIS/BouBel_emprise.shp"))
mybbox <- terra::project(AOIext,from="epsg:2154",to="epsg:4326")
REF10  <- terra::rast(xmin=1000*floor(AOIext$xmin/1000),
                     xmax=1000*ceiling(AOIext$xmax/1000),
                     ymin=1000*floor(AOIext$ymin/1000),
                     ymax=1000*ceiling(AOIext$ymax/1000),
                     crs=terra::crs("epsg:2154"),
                     res=10)

YEAR  <- 2022
VAR   <- c("AMPL","EOSD","EOSV","LSLOPE","MAXD","MAXV","MINV","RSLOPE","SOSD","SOSV")
LF    <- list.files(output_directory,pattern="ONSET")
tmp        <- lapply(LF,function(x) terra::rast(x)) #
rsrc       <- terra::sprc(tmp)
m          <- terra::mosaic(rsrc)
tmp        <- terra::project(m,REF10)

# ZOOM for test
ZOOM <- terra::ext(931000,936000,6454000,6459000)

graphics.off()
windows(20,20)
terra::plot(terra::crop(tmp,ZOOM),colNA="black",range=c(100,250))

# SUPP 1 : Elmore DL6 and DL7 parameters ----
# Note that this is similar to the 6 DL parameter in Elmore (eq. 4) 
# Elmore & al. Global Change Biol. 18, 656-674, doi:10.1111/j.1365-2486.2011.02521.x (2012).
DL6.f <- function(PAR,DOY){
  NDVIs  <- PAR[1] + (PAR[2] - PAR[1])*(1/(1+exp(-PAR[4]*(DOY-PAR[3])))-1/(1+exp(-PAR[6]*(DOY-PAR[5]))))
  return(NDVIs)
}
# where PAR[3] = PAR[4]*PAR[3], and PAR[5] = PAR[5]*PAR[6]

# Elmore & al. DL with 7 parameters
DL7.f <- function(PAR,DOY){
  NDVIs  <- PAR[1] + (PAR[2] - PAR[7]*DOY)*(1/(1+exp((PAR[3]-DOY)*PAR[4]))-1/(1+exp((PAR[5]-DOY)*PAR[6])))
  return(NDVIs)
}
# where m3p = PAR[3]   # m3/m4,  m4p = 1/PAR[4] # 1/m4, m5p = PAR[5] # m5/m6, m6p = 1/PAR[6] # 1/m6
# any analytical derivative of it ? 
START7 <- c(0,0.75,140,0.075,300,0.054,0.001)
plot(1:366,DL7.f(START7,1:366),ylim=c(0,1),type="l")

# SUPP 2 : Sampling the DL for Ecography -> DO.f  ----
# see bias.r (script for Ecography)
# simulates dates of observation as a function of the revisit interval (in days)
DO.d <- function(REVISIT){
  SAMP <- list()
  for (i in 1:(REVISIT-1)) SAMP[[i]] <- seq(i,365,REVISIT)
  return(SAMP)  # a list of length REVISIT-1 with vectors of observation dates in DOY
}

DO.d(16)

# Sample NDVIyts under varying revisit intervals and estimate NDVImax using a quantile function 
# NDVImax.q is a list of length REVISIT.seq

NSAMP       = 10000          # number of sampling
REVISIT.seq = seq(4,32,8)   # revisit interval in days

NDVImax.q  <- list()

QUANT <- quantile(x, probs = seq(0, 1, 0.25), na.rm = FALSE,
                  names = TRUE, type = 7, ...)

for (k in 1:length(REVISIT.seq)){

  print(k)
  DO.tmp   <- DO.d(REVISIT.seq[k])
  NOBSmin  <- min(unlist(lapply(DO.tmp,function(x) length(x)))) 
  # NOBSmin is the minimum number of dates for a given revisit interval
  
  NDVImax.tmp   <- matrix(NA,NOBSmin,ncol(NDVIyts))
  # initialize the matrix to store data
  
for(j in 1:ncol(NDVIyts)){

  print(j)
  
  for (i in 2:NOBSmin){

  tmp    <- sapply(1:NSAMP, function(x) NDVIyts[sample(DO.tmp[[sample(1:length(DO.tmp),1)]],i),j])
  
  # retrieve the NDVImax value using a quantile function
  # MED[i,j] <- median(apply(tmp,2,max))
  NDVImax.tmp[i,j]     <- median(apply(tmp,2,quantile,probs=0.75))
  
    }
}  

  NDVImax.q[[k]] <- NDVImax.tmp
}


# Figure 2
graphics.off()
windows(10,10)
K   <- 8
N4K <- nrow(NDVImax.q[[K]])
plot(spline(x=1:N4K,y=NDVImax.q[[K]][,1],n=2*N4K),type="l",ylim=c(0.2,0.8),xlab="Number of observation per year")
for(i in 2:length(ONSET.seq)) lines(spline(x=1:N4K,y=NDVImax.q[[K]][,i],n=2*N4K),col=rainbow(length(ONSET.seq))[i])
abline(h=0.75)












































