
suppressWarnings(suppressPackageStartupMessages({
    library(here)
}))

loadSavedAnno <- function (id) {

    anno <- readRDS(file=here::here("results", paste0(id, "_anno.rds")))
    return(anno) 
}

loadSavedMset <- function (id) {

    mset <- readRDS(file=here::here("results", paste0(id, "_mset.rds")))
    return(mset) 
}

loadSavedCoefs <- function (id) {

    coefs <- readRDS(file=here::here("results", paste0(id, "_coefs.rds")))
    return(coefs) 
}

loadSavedBetas <- function (id) {

    betas <- readRDS(file=here::here("results", paste0(id, "_betas.rds")))
    return(betas) 
}

loadSavedUmapBetas <- function (id) {

    umap_betas <- readRDS(file=here::here("results", paste0(id, "_umap_betas.rds")))
    return(umap_betas) 
}

loadSavedModel <- function (id) {

    rf_model <- readRDS(file=here::here("results", paste0(id, "_model.rds")))
    return(rf_model) 
}

loadSavedCalfit <- function (id) {

    calfit <- readRDS(file=here::here("results", paste0(id, "_calfit.rds")))
    return(calfit) 
}

loadSavedCnvRef <- function () {

    cnv_ref <- readRDS(file=here::here("results", "CNV_REF_EPICv2.rds"))
    return(cnv_ref) 
}

loadSavedCnvAnno <- function () {

    cnv_anno <- readRDS(file=here::here("results", "CNV_ANNO_EPICv2.rds"))
    return(cnv_anno) 
}