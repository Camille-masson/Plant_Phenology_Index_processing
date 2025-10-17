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
  






calcul_phenology_phase <- function (YEAR, input_brut_data_case, IRG_data_case, site,threshold ){
  
  
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
  thr10    <- abs(IRG_min) * threshold
  
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
  out_fp <- file.path(IRG_data_case, glue("Phenology_Phase_{site}_{YEAR}.tif"))
  writeRaster(phase_stack, out_fp, overwrite = TRUE, datatype = "INT1U")
  
  message(
    "✓ Phases phénologiques écrites dans ", out_fp, "\n",
    "   1 = pousse, 2 = plateau (|IRG| < 10 %), 3 = dépérissement, 4 = sénescence"
  )
  
  
  
  
  
  
  
  
}







#’IRG max et son DOY

IRG_MAX_processing <- function (input_IRG, site, YEAR, IRG_data_case) {
  IRG  <- rast(input_IRG)
  IRG_max <- app(IRG, fun = function(...) max(..., na.rm=TRUE))
  vals <- IRG_max[]; vals <- vals[!is.na(vals) & is.finite(vals)]
  hist(vals, breaks=seq(floor(min(vals)*100)/100, ceiling(max(vals)*100)/100, 0.01),
       xlim=c(floor(min(vals)*100)/100, ceiling(max(vals)*100)/100),
       main="Distribution de l'IRG max ", xlab="IRG max", ylab="Nombre de pixels")
  
  out_max <- file.path(IRG_data_case, paste0("IRG_max_", site, "_", YEAR, ".tif"))
  writeRaster(IRG_max, out_max, overwrite=TRUE)
  message("▶ IRG max save at :\n", out_max)
}















## ──────────────────────────────────────────────────────────────────────────────
##  Function 4 – Delta (jours) par rapport au maximum de verdissement (MAXD)
## ──────────────────────────────────────────────────────────────────────────────
#  Pour chaque jour de l’année (ou sous-période choisie), on calcule :
#        Δd = DOYcourant − MAXD
#  →  négatif avant MAXD, 0 le jour MAXD, positif après MAXD
#
#  • YEAR  : année numérique (ex. 2021)
#  • site  : nom du massif / alpage (ex. "Cayolle")
#  • input_brut_data_case : dossier contenant MAXD_<site>_<YEAR>.tif
#  • output_dir          : dossier où écrire Delta_day_of_IRGmax_*.tif
#  • DOY_range           : vecteur des DOY à traiter (défaut = 1:365)
#  • datatype            : format de sortie (INT2S = entier signé 16 bits ; couvre ±32 768)
#
#  Ex. d’appel :
#        delta_IRGmax_processing(2021, "Cayolle", input_brut_data_case,
#                                IRG_data_case, DOY_range = 121:334)
# ──────────────────────────────────────────────────────────────────────────────
delta_IRGmax_processing <- function(YEAR,
                                    site,
                                    input_brut_data_case,
                                    output_dir,
                                    DOY_range = 1:365,
                                    datatype  = "INT2S")
{
  library(terra)
  library(glue)
  
  ## 1. Lecture du MAXD et mise au propre --------------------------------------
  maxd_fp <- file.path(input_brut_data_case,
                       glue("MAXD_{site}_{YEAR}.tif"))
  if (!file.exists(maxd_fp))
    stop("Raster MAXD manquant : ", maxd_fp)
  
  MAXD <- rast(maxd_fp)
  
  # correction du « lag » (ex. 2021 → 21000) utilisé plus tôt dans le workflow
  lag  <- 1000 * as.numeric(substr(YEAR, 3, 4))
  MAXD <- MAXD - lag
  MAXD[MAXD < 1 | MAXD > 365] <- NA               # on garde uniquement 1–365
  valid <- is.finite(MAXD)                        # pixels d’intérêt
  
  ## 2. Calcul vectorisé du Δ jour par jour ------------------------------------
  #  On ré-utilise la géométrie de MAXD et on évite les boucles coûteuses :
  delta_stack <- rast(lapply(DOY_range, function(d) {
    # raster constant = d, puis on soustrait MAXD
    out <- setValues(MAXD, d) - MAXD
    out[!valid] <- NA
    out
  }))
  
  names(delta_stack) <- sprintf("Delta_DOY%03d", DOY_range)
  
  ## 3. Sauvegarde -------------------------------------------------------------
  out_fp <- file.path(output_dir,
                      glue("Delta_day_of_IRGmax_{site}_{YEAR}.tif"))
  writeRaster(delta_stack, out_fp,
              overwrite = TRUE,
              datatype  = datatype)
  
  message("✓ Raster Δ jours écrit : ", out_fp)
}



















###############################################################################
#  NDVI_processing  – Calcul et empilement saisonnier du NDVI
###############################################################################
NDVI_processing <- function(YEAR,
                            DOY_range,
                            input_brut_data_case,
                            NDVI_data_case,
                            site) {
  library(terra)
  
  ##--------------------------------------------------------------------------##
  ## 1. Chargement des rasters PARAMÈTRES (phase 2)
  ##--------------------------------------------------------------------------##
  params  <- c("MINV", "MAXV", "ONSET", "OFFSET", "GROWTH", "SENESC")
  
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
  
  ##--------------------------------------------------------------------------##
  ## 2. NORMALISATIONS
  ##    – a) Onset / Offset : YYYYDDD  →  DDD
  ##--------------------------------------------------------------------------##
  ONSET  <- ONSET  %% 1000
  OFFSET <- OFFSET %% 1000
  
  ONSET [ONSET  < 1 | ONSET  > 365] <- NA
  OFFSET[OFFSET < 1 | OFFSET > 365] <- NA
  
  ## Variante 1 - compatible toutes versions terra
  rng <- as.numeric( global(NDVImax, "max", na.rm = TRUE)[[1]] )   # -> nombre
  # ou, encore plus robuste :
  # rng <- as.numeric( unlist(global(NDVImax, "max", na.rm = TRUE)) )
  
  if (!is.na(rng) && is.finite(rng) && rng > 1.5) {
    NDVImin <- NDVImin / 10000
    NDVImax <- NDVImax / 10000
  }
  
  ##--------------------------------------------------------------------------##
  ## 3. Double logistique (niveau NDVI)
  ##--------------------------------------------------------------------------##
  DL.f <- function(DOY, NDVImin, NDVImax, ONSET, OFFSET, GROWTH, SENESC) {
    NDVImin +
      (NDVImax - NDVImin) *
      (  1 / (1 + exp(-GROWTH * (DOY - ONSET))) +
           1 / (1 + exp( SENESC * (DOY - OFFSET))) -
           1 )
  }
  
  ##--------------------------------------------------------------------------##
  ## 4. Boucle journalière & empilement
  ##--------------------------------------------------------------------------##
  ndvi_list <- lapply(DOY_range, function(doy) {
    lapp(
      c(NDVImin, NDVImax, ONSET, OFFSET, GROWTH, SENESC),
      fun      = DL.f,
      filename = "",
      DOY      = doy
    )
  })
  
  ndvi_stack <- rast(ndvi_list)
  names(ndvi_stack) <- paste0("NDVI_DOY", DOY_range)
  
  ## – borne stricte à [-1 ; 1] -----------------------------------------------
  ndvi_stack <- clamp(ndvi_stack, -100000, 1000000, values = TRUE)
  
  ##--------------------------------------------------------------------------##
  ## 5. Sauvegarde
  ##--------------------------------------------------------------------------##
  if (!dir.exists(NDVI_data_case))
    dir.create(NDVI_data_case, recursive = TRUE)
  
  out_season_ndvi <- file.path(
    NDVI_data_case,
    paste0("NDVI_season_", site, "_", YEAR, ".tif")
  )
  writeRaster(ndvi_stack, out_season_ndvi, overwrite = TRUE)
  
  message("▶ NDVI saisonnier empilé et sauvé sous :\n", out_season_ndvi)
}



###############################################################################
#  PPI_processing  – Calcul et empilement saisonnier du PPI
###############################################################################
PPI_processing <- function(YEAR,
                           DOY_range,
                           input_brut_data_case,
                           PPI_data_case,
                           site) {
  suppressPackageStartupMessages(library(terra))
  
  ##--------------------------------------------------------------------------##
  ## 1. Charger les rasters PARAMÈTRES (déjà normalisés à l’unité PPI)
  ##--------------------------------------------------------------------------##
  params  <- c("MINV", "MAXV", "ONSET", "OFFSET", "GROWTH", "SENESC")
  rasters <- setNames(
    lapply(params, function(p) {
      fp <- file.path(input_brut_data_case,
                      paste0(p, "_", site, "_", YEAR, ".tif"))
      if (!file.exists(fp))
        stop("Fichier introuvable : ", fp)
      rast(fp)
    }),
    params
  )
  
  PPImin <- rasters$MINV
  PPImax <- rasters$MAXV
  ONSET  <- rasters$ONSET        # jour 1‑365
  OFFSET <- rasters$OFFSET
  GROWTH <- rasters$GROWTH       # jour‑1
  SENESC <- rasters$SENESC
  
  ##--------------------------------------------------------------------------##
  ## 2. Contrôles rapides (pas de modification automatique)
  ##--------------------------------------------------------------------------##
  vmax <- global(PPImax, "max", na.rm = TRUE)[[1]]
  vmin <- global(PPImin, "min", na.rm = TRUE)[[1]]
  if (vmax > 3.5 || vmin < -0.05)
    warning("PPImin/max hors [0 ; 3] (min = ",
            round(vmin, 3), ", max = ", round(vmax, 3), ").")
  
  if (any(values(ONSET)  < 1 | values(ONSET)  > 365, na.rm = TRUE) ||
      any(values(OFFSET) < 1 | values(OFFSET) > 365, na.rm = TRUE))
    warning("ONSET ou OFFSET hors [1‑365] détectés (non modifiés).")
  
  ##--------------------------------------------------------------------------##
  ## 3. Double logistique (niveau PPI)
  ##--------------------------------------------------------------------------##
  DL.f <- function(DOY, PPImin, PPImax, ONSET, OFFSET, GROWTH, SENESC) {
    PPImin +
      (PPImax - PPImin) *
      (  1 / (1 + exp(-GROWTH * (DOY - ONSET))) +
           1 / (1 + exp( SENESC * (DOY - OFFSET))) -
           1 )
  }
  
  ##--------------------------------------------------------------------------##
  ## 4. Boucle journalière & empilement
  ##--------------------------------------------------------------------------##
  ppi_list <- lapply(DOY_range, function(doy) {
    lapp(
      c(PPImin, PPImax, ONSET, OFFSET, GROWTH, SENESC),
      fun      = DL.f,
      filename = "",
      DOY      = doy
    )
  })
  
  ppi_stack <- rast(ppi_list)
  names(ppi_stack) <- paste0("PPI_DOY", DOY_range)
  
  ## – borne stricte à [0 ; 3] ------------------------------------------------
  ppi_stack <- clamp(ppi_stack, 0, 3, values = FALSE)
  
  ##--------------------------------------------------------------------------##
  ## 5. Sauvegarde
  ##--------------------------------------------------------------------------##
  if (!dir.exists(PPI_data_case))
    dir.create(PPI_data_case, recursive = TRUE)
  
  out_season_ppi <- file.path(
    PPI_data_case,
    paste0("PPI_season_", site, "_", YEAR, ".tif")
  )
  writeRaster(ppi_stack, out_season_ppi, overwrite = TRUE)
  
  message("▶ PPI saisonnier empilé et sauvé sous :\n", out_season_ppi)
}



###############################################################################
# Fonction check_ndvi ---------------------------------------------------------
###############################################################################
check_ndvi <- function(ndvi_path, site, YEAR, n_sample = 1e5, seed = 123) {
  # 1. Charger les librairies (si pas déjà fait)
  suppressPackageStartupMessages({
    library(terra)
    library(ggplot2)
  })
  
  # 2. Lire le raster NDVI
  ndvi_stack <- rast(ndvi_path)
  
  # 3. Échantillonnage aléatoire jusqu’à n_sample valeurs
  set.seed(seed)
  vals <- spatSample(ndvi_stack,
                     size      = n_sample,
                     method    = "random",
                     as.raster = FALSE,
                     na.rm     = TRUE)[ , -c(1:2)]
  vals_vec <- as.numeric(unlist(vals))
  
  # 4. Création de l’histogramme avec ggplot2
  ggplot(data.frame(NDVI = vals_vec), aes(NDVI)) +
    geom_histogram(bins   = 100,
                   colour = "grey30",
                   fill   = "forestgreen",
                   alpha  = .6) +
    labs(
      title    = paste0("Distribution des valeurs NDVI/PPI – saison ", YEAR),
      subtitle = paste0("Site : ", site,
                        "  •  échantillon aléatoire de ",
                        format(n_sample, big.mark = " "), " pixels × jours"),
      x        = "Valeur NDVI (ou PPI rescalé)",
      y        = "Nombre d’observations"
    ) +
    theme_bw(base_size = 11)}




















###############################################################################
# Fonction check_irg ----------------------------------------------------------
###############################################################################
check_irg <- function(irg_season_path, site, YEAR,
                      n_sample = 1e5, seed = 123) {
  # 1. Charger les librairies
  suppressPackageStartupMessages({
    library(terra)
    library(ggplot2)
  })
  
  # 2. Lire le raster IRG saisonnier
  if (!file.exists(irg_season_path)) {
    stop("Fichier introuvable : ", irg_season_path)
  }
  irg_stack <- rast(irg_season_path)
  
  # 3. Échantillonnage aléatoire
  set.seed(seed)
  samp_df <- spatSample(
    irg_stack,
    size      = n_sample,
    method    = "random",
    as.raster = FALSE,
    na.rm     = TRUE
  )
  # samp_df a les colonnes x,y puis IRG_DOY...
  vals_mat <- as.data.frame(samp_df)[ , -c(1,2)]
  vals_vec <- as.numeric(unlist(vals_mat))
  
  # 4. Histogramme de la distribution IRG
  p1 <- ggplot(data.frame(IRG = vals_vec), aes(IRG)) +
    geom_histogram(
      bins   = 100,
      colour = "grey30",
      fill   = "steelblue",
      alpha  = .6
    ) +
    labs(
      title    = paste0("Distribution des valeurs IRG – saison ", YEAR),
      subtitle = paste0(
        "Site : ", site,
        "  • échantillon aléatoire de ",
        format(n_sample, big.mark = " "), " pixels × jours"
      ),
      x = "Valeur IRG (dérivée double logistique)",
      y = "Nombre d’observations"
    ) +
    theme_bw(base_size = 11)
  
  print(p1)
  
  # 5. Évolution moyenne IRG par DOY
  #    on calcule la moyenne de chaque colonne (chaque couche = un DOY)
  mean_irg <- colMeans(vals_mat)
  # extraire les DOY depuis les noms de colonnes
  DOYs <- as.integer(gsub("^IRG_DOY", "", colnames(vals_mat)))
  evo_df <- data.frame(
    DOY      = DOYs,
    mean_IRG = mean_irg
  )
  evo_df <- evo_df[order(evo_df$DOY), ]
  
  p2 <- ggplot(evo_df, aes(x = DOY, y = mean_IRG)) +
    geom_line(size = 1) +
    labs(
      title    = paste0("Évolution moyenne de l’IRG – saison ", YEAR),
      subtitle = paste0(
        "Site : ", site,
        "  • moyenne sur un échantillon aléatoire de ",
        format(n_sample, big.mark = " "), " pixels"
      ),
      x = "Jour de l’année (DOY)",
      y = "IRG moyen (dPPI/dt)"
    ) +
    theme_bw(base_size = 11)
  
  print(p2)
  
  invisible(list(histogram = p1, evolution = p2))
}











###############################################################################
# extra_indicators  – Stack RAW + indicateurs dérivés (17 bandes)
#    noms courts : MINV, MAXV, AMPL, SOSD, …, AsymSlope
###############################################################################
library(terra)

extra_indicators <- function(YEAR,
                             site,
                             input_brut_data_case,
                             extra_data_case) {
  
  ## 1. Rasters RAW à empiler (ONSET seulement pour le calcul de ONSET10)
  raw_keep <- c("MINV","MAXV","AMPL",
                "SOSD","EOSD","MAXD",
                "OFFSET",
                "GROWTH","SENESC",
                "LSLOPE","RSLOPE")
  
  read_r <- function(nm) rast(file.path(input_brut_data_case,
                                        paste0(nm, "_", site, "_", YEAR, ".tif")))
  
  raw_list <- lapply(raw_keep, read_r)
  names(raw_list) <- raw_keep
  
  ONSET <- read_r("ONSET")                     # non conservé dans le stack
  
  ## 2. Calcul des indicateurs dérivés
  MINV   <- raw_list$MINV ;  MAXV <- raw_list$MAXV ;  AMPL <- raw_list$AMPL
  SOSD   <- raw_list$SOSD ;  EOSD <- raw_list$EOSD ;  MAXD <- raw_list$MAXD
  OFFSET <- raw_list$OFFSET
  GROWTH <- raw_list$GROWTH
  LSLOPE <- raw_list$LSLOPE ; RSLOPE <- raw_list$RSLOPE
  
  LENGTH        <- EOSD - SOSD
  ONSET10       <- round(ONSET - log(9) / GROWTH)
  MaxSlope      <- AMPL * GROWTH / 4
  GreenUpDur    <- MAXD   - ONSET10
  GreenDownDur  <- OFFSET - MAXD
  AsymSlope     <- RSLOPE / LSLOPE
  
  ## Nettoyage
  ONSET10[ONSET10 < 1 | ONSET10 > 365] <- NA
  GreenUpDur [GreenUpDur  < 0]         <- NA
  GreenDownDur[GreenDownDur < 0]       <- NA
  LENGTH[LENGTH < 0]                   <- NA
  
  derived_list <- list(
    LENGTH       = LENGTH,
    ONSET10      = ONSET10,
    MaxSlope     = MaxSlope,
    GreenUpDur   = GreenUpDur,
    GreenDownDur = GreenDownDur,
    AsymSlope    = AsymSlope
  )
  
  ## 3. Empilement — Reduce pour être sûr d’avoir un SpatRaster
  full_stack <- Reduce(function(x, y) c(x, y), c(raw_list, derived_list))
  
  ## 4. Renommer les bandes de façon unique
  band_names <- c(names(raw_list), names(derived_list))  # 17 noms courts
  names(full_stack) <- band_names
  
  ## 5. Écriture
  if (!dir.exists(extra_data_case))
    dir.create(extra_data_case, recursive = TRUE)
  
  out_tif <- file.path(extra_data_case,
                       paste0("PPI_extra_", site, "_", YEAR, ".tif"))
  writeRaster(full_stack, out_tif,
              overwrite = TRUE, datatype = "FLT4S")
  
  message("▶ Stack RAW + indicateurs écrit (17 bandes) :\n", out_tif)
}



















 # ╔══════════════════════════════════════════════════════════════════════╗
  # ║  EOS10pct_processing()                                              ║
  # ║  → produit 1 raster : "EOS10pct_<site>_<YEAR>.tif"                  ║
  # ║    (DOY où PPI ≤ 0.10 × MAXV après le MAXD)                         ║
  # ╚══════════════════════════════════════════════════════════════════════╝
  EOS10pct_processing <- function(YEAR,
                                  DOY_range,
                                  input_brut_data_case,   # dossier raw_data
                                  PPI_data_case,          # dossier PPI
                                  site,
                                  out_dir = file.path(input_brut_data_case, "..", "extra"))
  {
    suppressPackageStartupMessages(library(terra))
    
    ## 1.  Rasters indispensables -----------------------------------------------
    ras_files <- function(x) file.path(input_brut_data_case,
                                       sprintf("%s_%s_%d.tif", x, site, YEAR))
    PPI_stack <- rast(file.path(PPI_data_case,
                                sprintf("PPI_season_%s_%d.tif", site, YEAR)))
    MAXV <- rast(ras_files("MAXV"))
    MAXD <- rast(ras_files("MAXD"))
    
    if (nlyr(PPI_stack) != length(DOY_range))
      stop("Le stack PPI ne correspond pas exactement au DOY_range fourni.")
    
    ## 2.  Préparation des valeurs en matrices -----------------------------------
    matPPI <- as.matrix(PPI_stack)      # n_pix × n_DOY
    thr10  <- 0.05 * values(MAXV)       # seuil 10 %
    dMax   <- values(MAXD)
    DOYvec <- DOY_range                 # vecteur des DOY
    
    ## 3.  Fonction pixel → DOY_EOS10 -------------------------------------------
    find_EOS10 <- function(pixPPI, thr, maxD) {
      
      if (is.na(thr) || is.na(maxD))
        return(NA_real_)
      
      idx <- which(DOYvec > maxD & pixPPI <= thr)   # après le pic & < 10 %
      if (length(idx) == 0) NA_real_ else DOYvec[idx[1]]
    }
    
    eos10 <- mapply(find_EOS10,
                    split(matPPI, row(matPPI)),
                    thr10, dMax,
                    SIMPLIFY = TRUE)
    
    ## 4.  Conversion en raster + export ----------------------------------------
    EOS10_r <- setValues(MAXD, eos10)               # même géométrie
    
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_file <- file.path(out_dir,
                          sprintf("EOS10pct_%s_%d.tif", site, YEAR))
    writeRaster(EOS10_r, out_file, overwrite = TRUE)
    cat("✓ EOS10pct enregistré :", out_file, "\n")
    
    invisible(out_file)
  }
  
  
  
  
  
  
  extra_indicators <- function(YEAR,
                               site,
                               input_brut_data_case,
                               extra_data_case) {
    library(terra)
    
    ## 1) Rasters RAW à empiler (ONSET seulement pour le calcul de ONSET10 & AUC)
    raw_keep <- c("MINV","MAXV","AMPL",
                  "SOSD","EOSD","MAXD",
                  "OFFSET",
                  "GROWTH","SENESC",
                  "LSLOPE","RSLOPE")
    
    read_r <- function(nm)
      rast(file.path(input_brut_data_case,
                     paste0(nm, "_", site, "_", YEAR, ".tif")))
    
    raw_list <- lapply(raw_keep, read_r)
    names(raw_list) <- raw_keep
    
    ONSET <- read_r("ONSET")  # utilisé mais pas empilé dans RAW
    
    ## 2) Raccourcis
    MINV   <- raw_list$MINV ;  MAXV <- raw_list$MAXV ;  AMPL <- raw_list$AMPL
    SOSD   <- raw_list$SOSD ;  EOSD <- raw_list$EOSD ;  MAXD <- raw_list$MAXD
    OFFSET <- raw_list$OFFSET
    GROWTH <- raw_list$GROWTH
    LSLOPE <- raw_list$LSLOPE ; RSLOPE <- raw_list$RSLOPE
    
    ## 3) Indicateurs dérivés “classiques”
    LENGTH        <- EOSD - SOSD
    ONSET10       <- round(ONSET - log(9) / GROWTH)   # ~10 % de l’amplitude
    MaxSlope      <- AMPL * GROWTH / 4
    GreenUpDur    <- MAXD   - ONSET10
    GreenDownDur  <- OFFSET - MAXD
    AsymSlope     <- RSLOPE / LSLOPE
    
    ## 4) AUC croissance entre ONSET10 et MAXD (logistique de croissance seule)
    AUC_growth_fun <- function(MINV, MAXV, ONSET, GROWTH, t1, t2){
      bad <- is.na(t1) | is.na(t2) | (t2 <= t1) | (GROWTH == 0)
      out <- MINV * (t2 - t1) +
        (MAXV - MINV) / GROWTH * (
          log1p(exp(GROWTH * (t2 - ONSET))) -
            log1p(exp(GROWTH * (t1 - ONSET)))
        )
      out[bad] <- NA
      out
    }
    
    AUCg_ONSET10_MAXD <- lapp(
      c(MINV, MAXV, ONSET, GROWTH, ONSET10, MAXD),
      fun = AUC_growth_fun
    )
    
    ## 5) Nettoyage
    ONSET10[ONSET10 < 1 | ONSET10 > 365] <- NA
    GreenUpDur [GreenUpDur  < 0]         <- NA
    GreenDownDur[GreenDownDur < 0]       <- NA
    LENGTH[LENGTH < 0]                   <- NA
    
    ## 6) Empilement
    derived_list <- list(
      LENGTH            = LENGTH,
      ONSET10           = ONSET10,
      MaxSlope          = MaxSlope,
      GreenUpDur        = GreenUpDur,
      GreenDownDur      = GreenDownDur,
      AsymSlope         = AsymSlope,
      AUCg_ONSET10_MAXD = AUCg_ONSET10_MAXD
    )
    
    full_stack <- Reduce(function(x, y) c(x, y), c(raw_list, derived_list))
    
    ## 7) Noms & écriture
    band_names <- c(names(raw_list), names(derived_list))  # 18 bandes
    names(full_stack) <- band_names
    
    if (!dir.exists(extra_data_case))
      dir.create(extra_data_case, recursive = TRUE)
    
    out_tif <- file.path(extra_data_case,
                         paste0("PPI_extra_", site, "_", YEAR, ".tif"))
    writeRaster(full_stack, out_tif,
                overwrite = TRUE, datatype = "FLT4S")
    
    message("▶ Stack RAW + indicateurs (incl. AUCg_ONSET10_MAXD) écrit (18 bandes) :\n",
            out_tif)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #' DeltaDay_MAXD_processing
  #' 
  #' Génère un empilement (stack) journalier où chaque couche représente, pour
  #' chaque pixel, le delta (DOY - MAXD), avec MAXD = jour julien du maximum
  #' (1–365). Les valeurs négatives indiquent des jours avant le maximum,
  #' positives après le maximum.
  #'
  #' @param YEAR integer. Année de travail (utilisée pour les noms de fichiers).
  #' @param DOY_range integer vector. Jours julien à traiter (ex. 60:365).
  #' @param input_brut_data_case character. Dossier où se trouve MAXD_*.tif.
  #' @param delta_data_case character. Dossier de sortie pour le stack delta.
  #' @param site character. Nom du site (utilisé pour les noms de fichiers).
  #'
  #' @return (invisible) chemin du fichier TIF écrit.
  #'
  #' @examples
  #' # Exemple d'appel, cohérent avec votre pipeline :
  #' # delta_data_case <- file.path(output_data_case, "delta_MAXD")
  #' # DeltaDay_MAXD_processing(YEAR, DOY_range, input_brut_data_case, delta_data_case, site)
  DeltaDay_MAXD_processing <- function(YEAR,
                                       DOY_range,
                                       input_brut_data_case,
                                       delta_data_case,
                                       site) {
    suppressPackageStartupMessages(library(terra))
    
    ##--------------------------------------------------------------------------##
    ## 1) Charger MAXD (jour julien du maximum)
    ##--------------------------------------------------------------------------##
    fp_maxd <- file.path(input_brut_data_case, paste0("MAXD_", site, "_", YEAR, ".tif"))
    if (!file.exists(fp_maxd))
      stop("Fichier introuvable : ", fp_maxd)
    MAXD <- rast(fp_maxd)
    
    ## Contrôles rapides
    if (!is.integer(DOY_range)) DOY_range <- as.integer(DOY_range)
    if (length(DOY_range) == 0L) stop("DOY_range est vide.")
    if (any(DOY_range < 1L | DOY_range > 365L))
      warning("Des DOY hors [1–365] ont été fournis (non modifiés).")
    
    gmin <- global(MAXD, "min", na.rm = TRUE)[[1]]
    gmax <- global(MAXD, "max", na.rm = TRUE)[[1]]
    if (!is.na(gmin) && (gmin < 1 || gmax > 365))
      warning("MAXD comporte des valeurs hors [1–365] (non modifiées).")
    
    ##--------------------------------------------------------------------------##
    ## 2) Calcul : pour chaque jour DOY, couche = DOY - MAXD
    ##--------------------------------------------------------------------------##
    delta_list <- lapply(DOY_range, function(doy) {
      doy - MAXD  # terra gère l'opération scalaire - raster
    })
    
    delta_stack <- rast(delta_list)
    names(delta_stack) <- paste0("DELTA_DOY", DOY_range)
    
    ## NB: Les deltas sont des entiers dans [-364 ; +364]. On écrit en INT2S.
    
    ##--------------------------------------------------------------------------##
    ## 3) Sauvegarde
    ##--------------------------------------------------------------------------##
    if (!dir.exists(delta_data_case)) dir.create(delta_data_case, recursive = TRUE)
    out_file <- file.path(delta_data_case, paste0("DeltaDay_from_MAXD_", site, "_", YEAR, ".tif"))
    
    writeRaster(delta_stack, out_file, overwrite = TRUE, datatype = "INT2S", NAflag = -32768)
    
    message("\u25B6 Empilement DeltaDay sauvegardé sous :\n", out_file)
    invisible(out_file)
  }
  
  
  # --------------------------------------------------------------------------- #
  # Mini-check visuel optionnel (facultatif) : retourne un raster unique pour un
  # DOY donné, pratique en debug.
  # --------------------------------------------------------------------------- #
  DeltaDay_for_one_DOY <- function(DOY, input_brut_data_case, site, YEAR) {
    suppressPackageStartupMessages(library(terra))
    fp_maxd <- file.path(input_brut_data_case, paste0("MAXD_", site, "_", YEAR, ".tif"))
    if (!file.exists(fp_maxd)) stop("Fichier introuvable : ", fp_maxd)
    MAXD <- rast(fp_maxd)
    out <- DOY - MAXD
    names(out) <- paste0("DELTA_DOY", DOY)
    out
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # =====================================================================
  #  Extraire la date du pic (PMAX_DOY) et la valeur au pic (PMAX_VAL)
  #  depuis le stack PPI_season_<site>_<year>.tif
  # =====================================================================
  PPI_extract_peak <- function(site,
                               YEAR,
                               PPI_data_case,
                               amp_min = 0.02,   # seuil amplitude pour filtrer non-végétation
                               write_out = TRUE) {
    
    suppressPackageStartupMessages(library(terra))
    
    # -- fichier d'entrée (créé par ta fonction PPI_processing)
    f_in <- file.path(PPI_data_case, paste0("PPI_season_", site, "_", YEAR, ".tif"))
    if (!file.exists(f_in)) stop("Introuvable : ", f_in)
    PPI <- rast(f_in)
    
    # -- vecteur DOY (lu dans les noms : "PPI_DOY60", "PPI_DOY61", ...)
    nm <- names(PPI)
    DOY <- suppressWarnings(as.integer(sub(".*DOY", "", nm)))
    if (any(is.na(DOY))) {
      message("Noms sans DOY explicite : j'assume DOY = 1..nlayers.")
      DOY <- seq_len(nlyr(PPI))
    }
    
    # -- fonction cellule : renvoie c(PMAX_DOY, PMAX_VAL) ou NA si non valide
    f_peak <- function(v, DOY, thr) {
      if (all(is.na(v))) return(c(NA_integer_, NA_real_))
      amp <- max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
      if (!is.finite(amp) || amp <= thr) return(c(NA_integer_, NA_real_))  # masque non-vég.
      i <- which.max(v)                        # en cas de plateau, prend l'EARLIEST pic
      c(DOY[i], v[i])
    }
    
    res <- app(PPI, f_peak, DOY = DOY, thr = amp_min)
    names(res) <- c("PMAX_DOY", "PMAX_VAL")
    
    if (write_out) {
      out_doy <- file.path(PPI_data_case, paste0("PMAX_DOY_", site, "_", YEAR, ".tif"))
      out_val <- file.path(PPI_data_case, paste0("PMAX_VAL_", site, "_", YEAR, ".tif"))
      
      writeRaster(res[[1]], out_doy, overwrite = TRUE,
                  wopt = list(datatype = "INT2U",
                              gdal = c("COMPRESS=DEFLATE","ZLEVEL=9")))
      writeRaster(res[[2]], out_val, overwrite = TRUE,
                  wopt = list(datatype = "FLT4S",
                              gdal = c("COMPRESS=DEFLATE","ZLEVEL=9")))
      
      message("▶ PMAX_DOY écrit : ", out_doy,
              "\n▶ PMAX_VAL écrit : ", out_val)
    }
    
    invisible(res)   # retourne aussi les rasters en mémoire
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # =====================================================================
  #  DOY à 10% du MAXV (premier passage, rampe montante)
  #  - lit PPI_season_<site>_<YEAR>.tif
  #  - renvoie un raster DOY_10pctMAX  (+ deux couches informatives : MAX_val, thr_val)
  # =====================================================================
  PPI_extract_thrdate <- function(site,
                                  YEAR,
                                  PPI_data_case,
                                  pct = 0.10,            # 10% du max (modifiable)
                                  amp_min = 0.02,        # masque zones peu végétalisées
                                  side = c("rising","falling","any"),
                                  relative_to = c("max","amp"),
                                  write_out = TRUE) {
    
    suppressPackageStartupMessages(library(terra))
    side <- match.arg(side)
    relative_to <- match.arg(relative_to)
    
    # -- entrée
    f_in <- file.path(PPI_data_case, paste0("PPI_season_", site, "_", YEAR, ".tif"))
    if (!file.exists(f_in)) stop("Introuvable : ", f_in)
    PPI <- rast(f_in)
    
    # -- DOY depuis les noms ("PPI_DOY60", ...)
    nm  <- names(PPI)
    DOY <- suppressWarnings(as.integer(sub(".*DOY", "", nm)))
    if (any(is.na(DOY))) DOY <- seq_len(nlyr(PPI))
    
    # -- fonction cellule
    f_thr <- function(v, DOY, pct, thr_amp, side, mode_rel){
      if (all(is.na(v))) return(c(NA_integer_, NA_real_, NA_real_))
      vmax <- max(v, na.rm=TRUE); vmin <- min(v, na.rm=TRUE); amp <- vmax - vmin
      if (!is.finite(amp) || amp <= thr_amp) return(c(NA_integer_, vmax, NA_real_))
      
      tval <- if (mode_rel == "max") pct * vmax else vmin + pct * amp
      
      idx <- which(v >= tval)
      if (length(idx) == 0) return(c(NA_integer_, vmax, tval))
      
      i <- switch(side,
                  rising  = min(idx),                 # premier passage
                  falling = max(idx),                 # dernier passage
                  any     = round(mean(idx)))         # milieu du plateau éventuel
      
      c(DOY[i], vmax, tval)
    }
    
    res <- app(PPI, f_thr,
               DOY = DOY, pct = pct, thr_amp = amp_min,
               side = side, mode_rel = relative_to)
    
    names(res) <- c(sprintf("DOY_%dpct_%s",
                            round(pct*100),
                            ifelse(relative_to=="max","MAX","AMP")),
                    "MAX_val","thr_val")
    
    if (write_out) {
      out_doy <- file.path(
        PPI_data_case,
        sprintf("DOY_%dpct%s_%s_%d.tif",
                round(pct*100),
                ifelse(relative_to=="max","MAX","AMP"),
                site, YEAR)
      )
      writeRaster(res[[1]], out_doy, overwrite=TRUE,
                  wopt=list(datatype="INT2U",
                            gdal=c("COMPRESS=DEFLATE","ZLEVEL=9")))
      message("▶ Écrit : ", out_doy)
    }
    
    invisible(res)
  }