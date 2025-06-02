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






#### 0. LIBRARY and PARAMETERS ####
#---------------------------------#
 
 # Source config
 source("config.R")
 
 # Parameters
 site      <- "OBJ"

 








#### 1. Downloading Plant Phenology Index (PPI) ####
#--------------------------------------------------#
if (TRUE) {
  
  ## DESCRIPTION
  #
  # Downloading PPI data for this, modify the parameters below:
  # - AREA: extent of your study zone
  # - YEARS: the desired time range
  #
  # The download of the tiles for all parameters will be automatically saved
  # in the folder input/data_"siste"/downloads
  
  ## FUNCTION 
  source(file.path(functions_case, "Function_processing.R"))
  
  
  ## PARAMETERS 
  AREA  <- c(6.2073128713, 44.5732596315, 6.4282245343, 44.7147025674)  # xmin,ymin,xmax,ymax
  YEARS <- 2023:2024                                                    # Load year
  
  VAR <- c("AMPL","EOSD","EOSV","LSLOPE","MAXD","MAXV","MINV","RSLOPE","SOSD","SOSV")
  
  
  # CONNECTION WEkEO 
  username <- "pcholer"
  password <- "Cardamine@2021"
  client   <- Client$new(username, password, save_credentials = TRUE)
  
  
  ## INPUT
  root_dir   <- getwd()
  data_case  <- file.path(input_case, paste0("data_", site))
  
  ## OUTPUT
  out_dir    <- file.path(data_case, "downloads")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  ## CODE
  
  data_download(out_dir, AREA, client, YEARS, VAR)


}

#### 2. Estimate the missing parameters  (ONSET, GROWTH, OFFSET, SENESC) ####
#---------------------------------------------------------------------------#
if (TRUE){
  ## DESCRIPTION
  #
  # Retrieves the 10 parameters, merges the different tiles if the study area requires multiple,
  # then reprojects to EPSG:2154 before cropping the mosaic tiles to the provided study area,
  # then extracts the 4 missing parameters.
  #
  # Requires:
  # - The 10 parameters downloaded in part 1 in the folder input/data_"siste"/downloads
  # - The study area boundary in .SHP format named Masque_"site" and placed in the folder
  #   created automatically below: input/data_"siste"/extent/Masque_"site"
  #
  # Generates:
  # - The ten parameters reprojected and cropped to the study area in format "VAR"_"site"_"Year".tif
  # - The 4 missing parameters: ONSET, GROWTH, OFFSET, SENESC in format "VAR"_"site"_"Year".tif
  #
  # All these parameters are saved in the following path: output/data_"site"/raw_data/"VAR"_"site"_"Year".tif
  
  
  ## FUNCTION 
  source(file.path(functions_case, "Function_estimate_parameters.R"))
  
  
  
  
  ## PARAMETERS 
  VAR  <- c("AMPL","EOSD","EOSV","LSLOPE","MAXD","MAXV","MINV","RSLOPE","SOSD","SOSV")
  resolution = 10
  YEAR = 2023
  
  
  
  ## INPOUT
  data_case        <- file.path(input_case, paste0("data_", site))
  downloads_case    <- file.path(data_case, "downloads")
  extent_case      <- file.path(data_case, "extent")
    
  # Un .SHP de la zone d'étude a DEPOSER dans le dossier extent
  mask_shp <- file.path(extent_case, paste0("Masque_", site, ".shp"))
  
  ## OUTPUT
  output_data_case  <- file.path(output_case, paste0("data_", site))
  brut_data_case <- file.path(output_data_case, "raw_data")
  
  
  for (d in c(input_case, data_case, extent_case, output_directory, output_data_case, brut_data_case)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  ## CODE
  
  
  estimate_parameters(VAR, mask_shp, resolution, downloads_case, YEAR, brut_data_case) 
  
  
  
}

#### 3. Calcul et empilement de l’IRG (Instantaneous Rate of Green‐up) ####
#------------------------------------------------------------------------#

if (TRUE){
  ## DESCRIPTION
  #  IRG processing based on PPI parameters previously computed in parts 1 & 2
  #
  #  Function 1– Compute the daily IRG stack.
  #
  #  Function 2– Derive phenological phases:
  #    1. Growth      (IRG positive up to MAXD)
  #    2. Maturity    (from MAXD to the threshold where |IRG| reaches 10 % of
  #                    its minimum — the “plateau”, bounded between MAXD and EOSD
  #                    so the zeros at the very beginning and end of season are
  #                    ignored)
  #    3. Decline     (from the 10 % threshold to EOSD)
  #    4. Senescence  (from EOSD to the end of the season)
  #
  #  Function 3– Extract IRG-max, i.e. the maximum IRG reached (normally at MAXD).
  #
  #  Requirements
  #    • PPI parameters already computed and re-projected for the study area
  #      in part 2.
  #
  #  Outputs
  #    • IRG_season_<site>_<YEAR>.tif   — multi-band stack (one band per DOY)
  #    • Phenology_Phase_<site>_<YEAR>.tif  — multi-band phenological phases
  #    • IRG_max_<site>_<YEAR>.tif      — IRG-max raster
  #
  #  All outputs are saved to:   output/data_<site>/IRG/
  
  
  
  ## FUNCTION 
  source(file.path(functions_case, "Function_calcul_indicators.R"))
  
  
  ## PARAMETERS
  site      <- "LALA"
  YEAR      <- 2023
  DOY_range <- 121:334
  
  
  
  
  ## INPUT
  input_data_case  <- file.path(output_case, paste0("data_", site))
  input_brut_data_case <- file.path(input_data_case, "raw_data")
  
  
  ## OUTPUT
  output_data_case  <- file.path(output_case, paste0("data_", site))
  IRG_data_case <- file.path(output_data_case, "IRG")
  
  for (d in c(IRG_data_case, input_brut_data_case)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  ## CODE
  
  # Function 1 : calcul IRG and stack in one rast by day (DOY)
  IRG_processing(YEAR, DOY, input_brut_data_case, IRG_data_case, site)
  
  # Function 2 : calcul the phenology phase with IRG
  calcul_phenology_phase(YEAR, input_brut_data_case, IRG_data_case, site)
  
  
  
  
  
 
  
  
  
  
  
  
  
  
  
  
  #’IRG max et son DOY
  
  IRG_max <- app(irg_stack, fun = function(...) max(..., na.rm=TRUE))
   vals <- IRG_max[]; vals <- vals[!is.na(vals) & is.finite(vals)]
   hist(vals, breaks=seq(floor(min(vals)*100)/100, ceiling(max(vals)*100)/100, 0.01),
        xlim=c(floor(min(vals)*100)/100, ceiling(max(vals)*100)/100),
        main="Distribution de l'IRG max (2023)", xlab="IRG max", ylab="Nombre de pixels")
   
   out_max <- file.path(output_directory, paste0("IRG_max_", site, "_", YEAR, ".tif"))
   writeRaster(IRG_max, out_max, overwrite=TRUE)
   message("▶ IRG max sauvé sous :\n", out_max)
   
   
   
   
   
   
   
   
   ## V1
   
   
   
   library(terra)
   library(glue)
   
   # ── 1. Repères (en DOY) ──────────────────────────────────────────────────────
   MAXD   <- rast(file.path(output_directory, glue("MAXD_{site}_{YEAR}.tif")))
   OFFSET <- rast(file.path(output_directory, glue("OFFSET_{site}_{YEAR}.tif")))
   EOSD   <- rast(file.path(output_directory, glue("EOSD_{site}_{YEAR}.tif")))
   
   lag  <- 1000 * as.numeric(substr(YEAR, 3, 4))   # 2023 → 23000
   MAXD <- MAXD - lag ; MAXD[MAXD < 1 | MAXD > 365] <- NA
   EOSD <- EOSD - lag ; EOSD[EOSD < 1 | EOSD > 365] <- NA
   valid <- is.finite(MAXD) & is.finite(EOSD)
   
   
   
   print(MAXD)
   print(OFFSET)
   print(EOSD)
   print(OM10)
   
   
   
   
   
   # ── 2. IRG empilé 121–334 ────────────────────────────────────────────────────
   irg_stack <- rast(file.path(output_directory, glue("IRG_season_{site}_{YEAR}.tif")))
   DOY_range <- 121:334                              # même ordre que les bandes IRG
   
   # seuil 10 % de |IRG_min| pour chaque pixel
   IRG_min  <- app(irg_stack, min, na.rm = TRUE)
   thr10    <- abs(IRG_min) * 0.10
   
   template <- MAXD                                  # support vierge
   
   # ── 3. Phase par jour ────────────────────────────────────────────────────────
   phase_stack <- rast(lapply(seq_along(DOY_range), function(i){
     
     d    <- DOY_range[i]
     irg  <- irg_stack[[i]]
     d_r  <- setValues(template, d)
     phase <- setValues(template, NA_integer_)
     
     ## 1 • pousse (avant MAXD) ou IRG >= 0
     sel <- valid & (d_r <= MAXD | irg >= 0)
     phase[sel] <- 1
     
     ## 2 • plateau : MAXD < d ≤ EOSD  &  |IRG| < 10 % |IRG_min|
     sel <- valid & d_r > MAXD & d_r <= EOSD & abs(irg) < thr10
     phase[sel] <- 2
     
     ## 3 • déperissement : dès que |IRG| ≥ 10 % ET d ≤ EOSD
     sel <- valid & d_r > MAXD & d_r <= EOSD & abs(irg) >= thr10
     phase[sel] <- 3
     
     ## 4 • sénescence : après EOSD
     sel <- valid & d_r > EOSD
     phase[sel] <- 4
     
     phase
   }))
   
   names(phase_stack) <- paste0("Phase_DOY", DOY_range)
   
   # ── 4. Sauvegarde (entiers 1-4) ───────────────────────────────────────────────
   out_fp <- file.path(output_directory, glue("Phenology_Phase_v3_{site}_{YEAR}.tif"))
   writeRaster(phase_stack, out_fp, overwrite = TRUE, datatype = "INT1U")
   
   message(
     "✓ Phases phénologiques écrites dans ", out_fp, "\n",
     "   1 = pousse, 2 = plateau (|IRG| < 10 %), 3 = dépérissement, 4 = sénescence"
   )
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   






}

























#### 3. Calcul of double logistic in PPI ####
#-------------------------------------------#

if (TRUE){
  library(terra)
  
  # 1. Paramètres
  
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















#### Bonus : partie 1 de philippe ####
#------------------------------------#





if (TRUE){
  
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
  
  
  
  
  
}




























