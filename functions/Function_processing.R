


data_download <- function (out_dir, AREA, client, YEARS, VAR) {
  for (yr in YEARS) {
    for (v in VAR) {
      
      ## si le tiff existe déjà, on passe
      already <- length(list.files(out_dir,
                                   pattern = paste0(yr, ".*", v, ".*\\.tif$"))) > 0
      if (already) next
      
      ## construire la requête JSON (sans httpAccept)
      query_json <- jsonlite::toJSON(list(
        dataset_id     = "EO:EEA:DAT:CLMS_HRVPP_VPP-LAEA",
        recordSchema   = "geojson",
        productType    = v,
        productGroupId = "s1",
        resolution     = "10",
        bbox           = AREA,
        startdate      = sprintf("%d-01-01T00:00:00.000Z", yr),
        enddate        = sprintf("%d-12-31T23:59:59.999Z", yr),
        itemsPerPage   = 200,
        startIndex     = 0
      ), auto_unbox = TRUE)
      
      cat("→", v, yr, "\n")
      res <- client$search(query_json)
      res$download(out_dir, force = FALSE, prompt = FALSE)
    }
  }
  
  
  }



# ---- robust_wekeo_download.R -----------------------------------------
library(jsonlite)

ensure_tc <- function(client) {
  id <- "Copernicus_Land_Monitoring_Service_Data_Policy"
  tc <- client$terms_and_conditions()
  if (!any(tc$term_id == id & tc$accepted)) client$terms_and_conditions(id)
}

search_safe <- function(client, body, retries = 5, wait = 2) {
  for (i in 0:retries) {
    res <- try(client$search(body), silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    if (!grepl("HTTP 500", as.character(res))) stop(res)
    Sys.sleep(wait * 2^i)
    try(client$get_token(), silent = TRUE)
  }
  stop("HDA indisponible après plusieurs tentatives (HTTP 500).")
}

# grille régulière dans la BBOX (deg)
tile_bbox <- function(AREA, dx = 0.25, dy = 0.25) {
  xmin <- AREA[1]; ymin <- AREA[2]; xmax <- AREA[3]; ymax <- AREA[4]
  xs <- seq(xmin, xmax, by = dx); if (tail(xs,1) < xmax) xs <- c(xs, xmax)
  ys <- seq(ymin, ymax, by = dy); if (tail(ys,1) < ymax) ys <- c(ys, ymax)
  tiles <- list()
  k <- 1
  for (i in seq_len(length(xs)-1)) {
    for (j in seq_len(length(ys)-1)) {
      tiles[[k]] <- c(xs[i], ys[j], xs[i+1], ys[j+1]); k <- k+1
    }
  }
  tiles
}

# tranches temporelles (par trimestre)
quarter_slices <- function(year) {
  list(
    c(sprintf("%d-01-01T00:00:00Z", year), sprintf("%d-04-01T00:00:00Z", year)),
    c(sprintf("%d-04-01T00:00:00Z", year), sprintf("%d-07-01T00:00:00Z", year)),
    c(sprintf("%d-07-01T00:00:00Z", year), sprintf("%d-10-01T00:00:00Z", year)),
    c(sprintf("%d-10-01T00:00:00Z", year), sprintf("%d-12-31T23:59:59Z", year))
  )
}

download_ppi_qflag_chunk <- function(out_dir, client, bbox, start_iso, end_iso,
                                     items = 50) {
  for (ptype in c("PPI","QFLAG")) {
    body <- toJSON(list(
      dataset_id   = "EO:EEA:DAT:CLMS_HRVPP_ST-LAEA",
      productType  = ptype,
      bbox         = bbox,
      start        = start_iso,
      end          = end_iso,
      itemsPerPage = items,
      startIndex   = 0
    ), auto_unbox = TRUE)
    res <- search_safe(client, body)
    if (!is.null(res) && isTRUE(res$total_count > 0)) {
      res$download(out_dir, force = FALSE, prompt = FALSE)
    }
  }
}

download_ppi_robust <- function(out_dir, AREA, client, year,
                                dx = 0.25, dy = 0.25, items = 50) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ensure_tc(client)
  tiles <- tile_bbox(AREA, dx, dy)
  slices <- quarter_slices(year)
  for (b in tiles) {
    for (s in slices) {
      message(sprintf("→ tile [%0.4f,%0.4f,%0.4f,%0.4f], %s → %s",
                      b[1],b[2],b[3],b[4], s[1], s[2]))
      try(download_ppi_qflag_chunk(out_dir, client, b, s[1], s[2], items),
          silent = TRUE)
    }
  }
  message("✓ Téléchargement PPI/QFLAG terminé (mode robuste).")
}
















