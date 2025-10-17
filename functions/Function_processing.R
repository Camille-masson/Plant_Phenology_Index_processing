


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