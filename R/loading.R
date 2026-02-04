
loadSavedAnno <- function (id) {

    anno <- readRDS(file=file.path("results", paste0(id, "_anno.rds")))
    return(anno) 
}

loadSavedMset <- function (id) {

    mset <- readRDS(file=file.path("results", paste0(id, "_mset.rds")))
    return(mset) 
}

loadSavedCoefs <- function (id) {

    coefs <- readRDS(file=file.path("results", paste0(id, "_coefs.rds")))
    return(coefs) 
}

loadSavedBetas <- function (id) {

    betas <- readRDS(file=file.path("results", paste0(id, "_betas.rds")))
    return(betas) 
}

loadSavedFilteredBetas <- function (id) {

    filtered_betas <- readRDS(file=file.path("results", paste0(id, "_filtered_betas.rds")))
    return(filtered_betas) 
}

loadSavedModel <- function (id) {

    rf_model <- readRDS(file=file.path("results", paste0(id, "_model.rds")))
    return(rf_model) 
}

loadSavedCalfit <- function (id) {

    calfit <- readRDS(file=file.path("results", paste0(id, "_calfit.rds")))
    return(calfit) 
}

loadSavedCnvRef <- function () {

    cnv_ref <- readRDS(file=file.path("results", "CNV_REF_EPICv2.rds"))
    return(cnv_ref) 
}

loadSavedCnvAnno <- function () {

    cnv_anno <- readRDS(file=file.path("results", "CNV_ANNO_EPICv2.rds"))
    return(cnv_anno) 
}