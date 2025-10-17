
estimate_parameters <- function(VAR, mask_shp, resolution, downloads_case, YEAR, brut_data_case, site) {
  
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















































