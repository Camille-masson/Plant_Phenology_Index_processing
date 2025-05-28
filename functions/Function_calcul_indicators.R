IRG_processing <- function(YEAR, DOY, input_brut_data_case, IRG_data_case, site){
  
  
  
  
  # 1. Chargement des rasters dérivés (déjà produits en partie 2)
  params <- c("MINV","MAXV","ONSET","OFFSET","GROWTH","SENESC")
  rasters <- setNames(
    lapply(params, function(p) {
      fp <- file.path(input_brut_data_case, paste0(p, "_", site, "_", YEAR, ".tif"))
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
  
  # 2. Définition de la dérivée de la double logistique (IRG)
  DLd.f <- function(DOY, NDVImin, NDVImax, ONSET, OFFSET, GROWTH, SENESC) {
    (NDVImax - NDVImin) * (
      GROWTH * exp(-GROWTH * (DOY - ONSET)) / (1 + exp(-GROWTH * (DOY - ONSET)))^2 -
        SENESC * exp(-SENESC * (DOY - OFFSET)) / (1 + exp(-SENESC * (DOY - OFFSET)))^2
    )
  }
  
  # 4. Empilement IRG sur la saison (DOY)
  
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
    IRG_data_case,
    paste0("IRG_season_", site, "_", YEAR, ".tif")
  )
  writeRaster(irg_stack, out_season_irg, overwrite = TRUE)
  message("▶ IRG saisonnier empilé et sauvé sous :\n", out_season_irg)
  
  
  
  
  
  
  
  
  
}
  


calcul_phenology_phase <- function (YEAR, input_brut_data_case, IRG_data_case, site){
  
  
  library(terra)
  library(glue)
  
  # ── 1. Repères (en DOY) ──────────────────────────────────────────────────────
  MAXD   <- rast(file.path(input_brut_data_case, glue("MAXD_{site}_{YEAR}.tif")))
  OFFSET <- rast(file.path(input_brut_data_case, glue("OFFSET_{site}_{YEAR}.tif")))
  EOSD   <- rast(file.path(input_brut_data_case, glue("EOSD_{site}_{YEAR}.tif")))
  
  lag  <- 1000 * as.numeric(substr(YEAR, 3, 4))   # 2023 → 23000
  MAXD <- MAXD - lag ; MAXD[MAXD < 1 | MAXD > 365] <- NA
  EOSD <- EOSD - lag ; EOSD[EOSD < 1 | EOSD > 365] <- NA
  valid <- is.finite(MAXD) & is.finite(EOSD)
  
  
  
  print(MAXD)
  print(OFFSET)
  print(EOSD)
  
  
  
  
  
  
  # ── 2. IRG empilé 121–334 ────────────────────────────────────────────────────
  irg_stack <- rast(file.path(IRG_data_case, glue("IRG_season_{site}_{YEAR}.tif")))
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
  out_fp <- file.path(IRG_data_case, glue("Phenology_Phase_v3_{site}_{YEAR}.tif"))
  writeRaster(phase_stack, out_fp, overwrite = TRUE, datatype = "INT1U")
  
  message(
    "✓ Phases phénologiques écrites dans ", out_fp, "\n",
    "   1 = pousse, 2 = plateau (|IRG| < 10 %), 3 = dépérissement, 4 = sénescence"
  )
  
  
  
  
  
  
  
  
}






