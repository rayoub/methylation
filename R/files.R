
suppressWarnings(suppressPackageStartupMessages({
    library(here)
}))

loadGeoData <- function (batch_id, suffix) {

    obj <- readRDS(file=here::here("results", "geo", paste0(batch_id, "_", suffix, ".rds")))
    return(obj) 
}

loadLabData <- function (batch_id, suffix) {

    obj <- readRDS(file=here::here("results", "lab", paste0(batch_id, "_", suffix, ".rds")))
    return(obj) 
}

saveGeoData <- function (batch_id, suffix, obj) {

    saveRDS(obj, file=here::here("results", "geo", paste0(batch_id, "_", suffix, ".rds")))
}

saveLabData <- function (batch_id, suffix, obj) {

    saveRDS(obj, file=here::here("results", "lab", paste0(batch_id, "_", suffix, ".rds")))
}