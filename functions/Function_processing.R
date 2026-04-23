
data_download <- function(out_dir, AREA, client, YEARS, VAR, S = c("s2")) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (yr in YEARS) {
    for (v in VAR) {
      for (s in S) {
        
        # si le tif existe déjà pour CE cycle, on passe
        already <- length(list.files(
          out_dir,
          pattern = paste0(yr, ".*", v, ".*", s, ".*\\.tif$"),
          ignore.case = TRUE
        )) > 0
        if (already) {
          cat("✓ déjà présent :", v, yr, s, "\n")
          next
        }
        
        query_json <- jsonlite::toJSON(list(
          dataset_id     = "EO:EEA:DAT:CLMS_HRVPP_VPP-LAEA",
          recordSchema   = "geojson",
          productType    = v,
          productGroupId = s,        # <- "s1" ou "s2"
          resolution     = "10",
          bbox           = AREA,
          startdate      = sprintf("%d-01-01T00:00:00.000Z", yr),
          enddate        = sprintf("%d-12-31T23:59:59.999Z", yr),
          itemsPerPage   = 200,
          startIndex     = 0
        ), auto_unbox = TRUE)
        
        cat("→", v, yr, s, "\n")
        res <- client$search(query_json)
        res$download(out_dir, force = FALSE, prompt = FALSE)
      }
    }
  }
}
