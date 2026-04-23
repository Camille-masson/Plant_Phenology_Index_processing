
estimate_parameters <- function(VAR,
                                template,
                               downloads_case,
                               YEAR,
                               brut_data_case,
                               site) {
  ## ------------------------------------------------------------------ ##
  ## 1. Zone d’étude : masque EPSG 2154 et grille REF10 (10 m)          ##
  ## ------------------------------------------------------------------ ##
  
  REF10 <- rast(template)
  
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










library(terra)

estimate_parameters <- function(VAR,
                                template,        
                                downloads_case,
                                YEAR,
                                brut_data_case,
                                site,
                                cycle = "s2") {  
  
  ## 1) Raster de référence (grille parfaite)
  REF10 <- terra::rast(template)
  
  ## 2) dossier de sortie séparé par cycle (évite mélange s1/s2)
  brut_cycle <- file.path(brut_data_case, cycle)
  if (!dir.exists(brut_cycle)) dir.create(brut_cycle, recursive = TRUE)
  
  ## 3) catégories de produits
  ndvi_like  <- c("MINV","MAXV","SOSV","EOSV","AMPL")
  slope_like <- c("LSLOPE","RSLOPE")
  doy_like   <- c("SOSD","EOSD","MAXD","SOST","MAXT","EOST")
  
  for (prod in VAR) {
    
    # -> IMPORTANT : on filtre aussi sur le cycle
    files <- list.files(
      downloads_case,
      pattern = paste0("^VPP_", YEAR, ".*_", cycle, "_", prod, "\\.tif$"),
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    if (!length(files)) {
      warning("Manque ", prod, " pour ", YEAR, " ", cycle)
      next
    }
    
    mos_lae <- terra::mosaic(terra::sprc(lapply(files, terra::rast)))
    
    # reprojection + resampling exactement sur la grille REF10
    mos_l93 <- terra::project(mos_lae, REF10)
    
    ## -- dates YYYYDDD
    if (prod %in% doy_like) {
      mos_l93 <- mos_l93 %% 1000
      mos_l93[mos_l93 < 1 | mos_l93 > 365] <- NA
    }
    
    ## -- valeurs NDVI/PPI
    if (prod %in% ndvi_like) {
      vmax <- terra::global(mos_l93, "max", na.rm = TRUE)[[1]]
      if (!is.na(vmax) && vmax > 3.5) {
        sf <- if (vmax > 30000) 30000 else 10000
        mos_l93 <- mos_l93 / sf
        message("→ ", prod, " normalisé (/ ", sf, ")")
      }
    }
    
    ## -- slopes
    if (prod %in% slope_like) {
      vmed <- terra::global(mos_l93, "mean", na.rm = TRUE)[[1]]
      if (!is.na(vmed) && vmed > 0.5) {
        mos_l93 <- mos_l93 / 10000
        message("→ ", prod, " normalisé (/10 000)")
      } else {
        message("→ ", prod, " déjà normalisé (mean = ", round(vmed, 3), ")")
      }
    }
    
    # crop/mask sur l’emprise du template (simple, fiable)
    mos_crop <- terra::crop(mos_l93, REF10)
    mos_mask <- terra::mask(mos_crop, REF10)
    
    out_file <- file.path(brut_cycle, paste0(prod, "_", site, "_", YEAR, "_", cycle, ".tif"))
    terra::writeRaster(mos_mask, out_file, overwrite = TRUE)
  }
  
  ## 4) calcul des paramètres dynamiques à partir des fichiers sortis
  lf <- list.files(brut_cycle,
                   pattern = paste0("_", site, "_", YEAR, "_", cycle, "\\.tif$"),
                   full.names = TRUE)
  
  rget <- function(pat) terra::rast(lf[grep(pat, lf)])
  
  AMPL   <- rget("AMPL");   MINV <- rget("MINV")
  SOSV   <- rget("SOSV");    EOSV <- rget("EOSV")
  LSLOPE <- rget("LSLOPE"); RSLOPE <- rget("RSLOPE")
  SOSD   <- rget("SOSD");    EOSD <- rget("EOSD")
  
  # verdissement
  A <- (AMPL)/(SOSV - MINV) - 1
  A[A <= 0] <- NA
  
  GROWTH <- LSLOPE / (A/(1+A)^2 * AMPL)
  GROWTH[GROWTH <= 0 | GROWTH > 2] <- NA
  terra::writeRaster(round(GROWTH, 6),
                     file.path(brut_cycle, paste0("GROWTH_", site, "_", YEAR, "_", cycle, ".tif")),
                     overwrite = TRUE)
  
  ONSET <- (log(A) + GROWTH*SOSD)/GROWTH
  ONSET[ONSET < 1 | ONSET > 365 | (ONSET - SOSD) < 0] <- NA
  terra::writeRaster(round(ONSET),
                     file.path(brut_cycle, paste0("ONSET_", site, "_", YEAR, "_", cycle, ".tif")),
                     overwrite = TRUE)
  
  # sénescence
  B <- (AMPL)/(EOSV - MINV) - 1
  B[B <= 0] <- NA
  
  SENESC <- (1 + 1/B)^2 * B * RSLOPE / AMPL
  SENESC[SENESC <= 0 | SENESC > 2] <- NA
  terra::writeRaster(round(SENESC, 6),
                     file.path(brut_cycle, paste0("SENESC_", site, "_", YEAR, "_", cycle, ".tif")),
                     overwrite = TRUE)
  
  OFFSET <- EOSD - log(B)/SENESC
  OFFSET[OFFSET < 1 | OFFSET > 365] <- NA
  terra::writeRaster(round(OFFSET),
                     file.path(brut_cycle, paste0("OFFSET_", site, "_", YEAR, "_", cycle, ".tif")),
                     overwrite = TRUE)
  
  message("✓ Paramètres dynamiques calculés pour ", YEAR, " (", cycle, ")")
}
