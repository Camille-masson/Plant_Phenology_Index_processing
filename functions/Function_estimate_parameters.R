
estimate_parameters <- function(VAR, mask_shp, resolution, downloads_case, YEAR, brut_data_case) {
  
  # 1. Chargement du masque L93 et création de la grille REF10 (L93, 10 m)
  mask_vec <- try(vect(mask_shp), silent = TRUE)  # doit être en EPSG:2154
  if (inherits(mask_vec, "try-error") || length(mask_vec) == 0) {
    warning("Warning !! : the shapefile is missing or could not be read: ", mask_shp)
    next
  }
  
  AOIext   <- ext(mask_vec)
  REF10    <- rast(
    xmin = 1000 * floor(AOIext$xmin / 1000),
    xmax = 1000 * ceiling(AOIext$xmax / 1000),
    ymin = 1000 * floor(AOIext$ymin / 1000),
    ymax = 1000 * ceiling(AOIext$ymax / 1000),
    crs  = crs(mask_vec),  # EPSG:2154
    res  = resolution
  )
  
  # 2. Assemble & rééchelle
  for (prod in VAR) {
    # a) Liste des tuiles brutes (LAEA) pour l’année et le produit
    files <- list.files(
      downloads_case,
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
      brut_data_case,
      paste0(prod, "_", site, "_", YEAR, ".tif")
    )
    terra::writeRaster(mos_mask, out_file, overwrite = TRUE)
    message("→ Écrit : ", out_file)
  }
  
  # on liste tous les rasters <VAR>_<site>_<YEAR>.tif
  LF <- list.files(
    brut_data_case,
    pattern    = paste0("_", site, "_", YEAR, "\\.tif$"),
    full.names = TRUE
  )
  
  # une seule itération (grep("") renvoie tous les fichiers)
  TILES <- ""
  for (j in TILES) {
    print(j)
    FILES <- LF[grep(j, LF)] # devrait renvoyer 10 fichiers
    
    lag    <- 1000 * as.numeric(substr(YEAR, 3, 4))
    AMPL   <- terra::rast(FILES[grep("AMPL", FILES)])
    MINV   <- terra::rast(FILES[grep("MINV", FILES)])
    SOSV   <- terra::rast(FILES[grep("SOSV", FILES)])
    EOSV   <- terra::rast(FILES[grep("EOSV", FILES)])
    LSLOPE <- terra::rast(FILES[grep("LSLOPE", FILES)])
    RSLOPE <- terra::rast(FILES[grep("RSLOPE", FILES)])
    SOSD   <- terra::rast(FILES[grep("SOSD", FILES)]) - lag
    SOSD[SOSD < 0 | SOSD > 365] <- NA
    EOSD   <- terra::rast(FILES[grep("EOSD", FILES)]) - lag
    EOSD[EOSD < 0 | EOSD > 365] <- NA
    
    # let A = exp(-GROWTH*(SOSD-ONSET))
    # Compute A using DLG.f(SOSV)
    A       <- (AMPL) / (SOSV - MINV) - 1
    A[A <= 0] <- NA  
    
    # Compute GROWTH using DLGd.f(SOSV,LSLOPE)
    GROWTH             <- round(LSLOPE / (A / (1 + A)^2 * (AMPL)), 4)
    GROWTH[GROWTH[] < 0]   <- NA
    GROWTH[GROWTH[] > 0.5] <- NA
    hist(GROWTH[], breaks = seq(0, 0.5, 0.01), xlim = c(0, 0.5), plot = FALSE)
    terra::writeRaster(GROWTH, gsub("AMPL", "GROWTH", FILES[1]), overwrite = TRUE)
    
    # Compute ONSET using DIG.f(SOSD)
    ONSET <- round((log(A) + GROWTH * SOSD) / GROWTH)
    ONSET[(ONSET - SOSD) < 0] <- NA
    ONSET[ONSET > 365]       <- NA
    hist(ONSET[], breaks = seq(0, 365, 1), xlim = c(0, 365), plot = FALSE)
    terra::writeRaster(ONSET, gsub("AMPL", "ONSET", FILES[1]), overwrite = TRUE)
    
    # let B = exp(SENESC*(EOSD-OFFSET))
    # Compute B using DLS.f(EOSV)
    B       <- (AMPL) / (EOSV - MINV) - 1
    B[B <= 0] <- NA  
    
    # Compute SENESC using DLSd.f(EOSV,RSLOPE)
    SENESC <- round((1 + 1 / B)^2 * B * RSLOPE / AMPL, 4)
    SENESC[SENESC[] <= 0]   <- NA
    SENESC[SENESC[] > 0.5]  <- NA
    hist(SENESC[], breaks = seq(0, 0.5, 0.01), xlim = c(0, 0.5), plot = FALSE)
    terra::writeRaster(SENESC, gsub("AMPL", "SENESC", FILES[1]), overwrite = TRUE)
    
    # Compute OFFSET using DIS.f(EOSD)
    OFFSET <- round(EOSD - log(B) / SENESC)
    OFFSET[OFFSET > 365] <- NA
    OFFSET[OFFSET < 0]   <- NA
    hist(OFFSET[], breaks = seq(0, 365, 1), xlim = c(0, 365), plot = FALSE)
    terra::writeRaster(OFFSET, gsub("AMPL", "OFFSET", FILES[1]), overwrite = TRUE)
  }
  
}






















###############################################################################
#  estimate_parameters()  – version corrigée
#   • pas d’arrondi pré‑mature  (GROWTH / SENESC)
#   • fenêtre de détection LSLOPE / RSLOPE plus robuste
#   • vérifications post‑écriture conservées
###############################################################################
library(terra)

estimate_parameters <- function(VAR,
                                mask_shp,
                                resolution      = 10,
                                downloads_case,
                                YEAR,
                                brut_data_case,
                                site) {
  ## ------------------------------------------------------------------ ##
  ## 1. Zone d’étude : masque EPSG 2154 et grille REF10 (10 m)          ##
  ## ------------------------------------------------------------------ ##
  mask_vec <- vect(mask_shp)
  AOIext   <- ext(mask_vec)
  REF10    <- rast(xmin = 1000*floor(AOIext$xmin/1000),
                   xmax = 1000*ceiling(AOIext$xmax/1000),
                   ymin = 1000*floor(AOIext$ymin/1000),
                   ymax = 1000*ceiling(AOIext$ymax/1000),
                   crs  = crs(mask_vec), res = resolution)
  
  ## ------------------------------------------------------------------ ##
  ## 2. Pré‑traitement LAEA → L93, normalisations automatiques          ##
  ## ------------------------------------------------------------------ ##
  ndvi_like  <- c("MINV","MAXV","SOSV","EOSV","AMPL")
  slope_like <- c("LSLOPE","RSLOPE")
  doy_like   <- c("SOSD","EOSD","MAXD","SOST","MAXT","EOST")
  
  for (prod in VAR) {
    files <- list.files(downloads_case,
                        pattern = paste0(YEAR, ".*", prod, ".*\\.tif$"),
                        full.names = TRUE)
    if (!length(files)) { warning("Manque ", prod); next }
    
    mos_lae <- mosaic(sprc(lapply(files, rast)))
    mos_l93 <- project(mos_lae, REF10)
    
    ## -- dates YYYYDDD -------------------------------------------------------
    if (prod %in% doy_like) {
      mos_l93 <- mos_l93 %% 1000
      mos_l93[mos_l93 < 1 | mos_l93 > 365] <- NA
    }
    
    ## -- valeurs NDVI/PPI ----------------------------------------------------
    if (prod %in% ndvi_like) {
      vmax <- global(mos_l93, "max", na.rm = TRUE)[[1]]
      if (!is.na(vmax) && vmax > 3.5) {
        sf <- if (vmax > 30000) 30000 else 10000
        mos_l93 <- mos_l93 / sf
        message("→ ", prod, " normalisé (/ ", sf, ")")
      }
    }
    
    ## -- slopes : détection sur la moyenne (plus fiable que max) ------------
    if (prod %in% slope_like) {
      vmed <- global(mos_l93, "mean", na.rm = TRUE)[[1]]
      if (!is.na(vmed) && vmed > 0.5) {        # encore codé 0‑30000
        mos_l93 <- mos_l93 / 10000
        message("→ ", prod, " normalisé (/10 000)")
      } else {
        message("→ ", prod, " semble déjà normalisé (mean = ",
                round(vmed, 3), ")")
      }
    }
    
    mos_mask <- mask(crop(mos_l93, mask_vec), mask_vec)
    writeRaster(mos_mask,
                file.path(brut_data_case,
                          paste0(prod, "_", site, "_", YEAR, ".tif")),
                overwrite = TRUE)
  }
  
  ## ------------------------------------------------------------------ ##
  ## 3. Calcul des paramètres dynamiques                                ##
  ## ------------------------------------------------------------------ ##
  lf <- list.files(brut_data_case,
                   pattern = paste0("_", site, "_", YEAR, "\\.tif$"),
                   full.names = TRUE)
  rget <- \(pat) rast(lf[grep(pat, lf)])
  
  AMPL   <- rget("AMPL"); MINV <- rget("MINV")
  SOSV   <- rget("SOSV");  EOSV <- rget("EOSV")
  LSLOPE <- rget("LSLOPE"); RSLOPE <- rget("RSLOPE")
  SOSD   <- rget("SOSD");  EOSD <- rget("EOSD")
  
  ## -- phase de verdissement -----------------------------------------------
  A <- (AMPL)/(SOSV - MINV) - 1;  A[A <= 0] <- NA
  
  GROWTH <- LSLOPE / (A/(1+A)^2 * AMPL)          # PAS d’arrondi ici
  GROWTH[GROWTH <= 0 | GROWTH > 2] <- NA
  writeRaster(round(GROWTH, 6),
              gsub("AMPL","GROWTH", lf[grep("AMPL",lf)][1]), overwrite = TRUE)
  
  ONSET <- (log(A) + GROWTH*SOSD)/GROWTH
  ONSET[ONSET < 1 | ONSET > 365 | (ONSET - SOSD) < 0] <- NA
  writeRaster(round(ONSET), gsub("AMPL","ONSET",lf[grep("AMPL",lf)][1]),
              overwrite = TRUE)
  
  ## -- phase de sénescence ---------------------------------------------------
  B <- (AMPL)/(EOSV - MINV) - 1;  B[B <= 0] <- NA
  
  SENESC <- (1 + 1/B)^2 * B * RSLOPE / AMPL
  SENESC[SENESC <= 0 | SENESC > 2] <- NA
  writeRaster(round(SENESC, 6),
              gsub("AMPL","SENESC",lf[grep("AMPL",lf)][1]), overwrite = TRUE)
  
  OFFSET <- EOSD - log(B)/SENESC
  OFFSET[OFFSET < 1 | OFFSET > 365] <- NA
  writeRaster(round(OFFSET),
              gsub("AMPL","OFFSET",lf[grep("AMPL",lf)][1]), overwrite = TRUE)
  
  message("✓ Paramètres dynamiques calculés pour ", YEAR)
  
  ## ------------------------------------------------------------------ ##
  ## 4. Résumé min‑max (optionnel)                                       ##
  ## ------------------------------------------------------------------ ##
  rng <- sapply(list.files(brut_data_case,
                           pattern = paste0("_", site, "_", YEAR, ".tif$"),
                           full.names = TRUE),
                \(f) c(min = global(rast(f),"min",na.rm=TRUE)[[1]],
                       max = global(rast(f),"max",na.rm=TRUE)[[1]]))
  print(round(rng, 4))
}






















































