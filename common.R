
REF_GSE_ID <- "GSE90496"
VAL_GSE_ID <- "GSE109379"

loadSavedAnno <- function (gse_id, save_tag = "") {

    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID)) 
    anno <- readRDS(file=file.path("results", paste0(save_tag, gse_id, "_anno.rds")))
    return(anno) 
}

loadSavedMset <- function (gse_id, save_tag = "") {

    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID)) 
    mset <- readRDS(file=file.path("results", paste0(save_tag, gse_id, "_mset.rds")))
    return(mset) 
}

loadSavedCoefs <- function (gse_id, save_tag = "") {

    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID)) 
    coefs <- readRDS(file=file.path("results", paste0(save_tag, gse_id, "_coefs.rds")))
    return(coefs) 
}

loadSavedBetas <- function (gse_id, save_tag = "") {

    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID)) 
    betas <- readRDS(file=file.path("results", paste0(save_tag, gse_id, "_betas.rds")))
    return(betas) 
}

loadSavedFilteredBetas <- function (gse_id, save_tag = "") {

    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID)) 
    betas <- readRDS(file=file.path("results", paste0(save_tag, gse_id, "_filtered_betas.rds")))
    return(betas) 
}

loadSavedModel <- function (gse_id, save_tag = "") {

    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID)) 
    rf.model <- readRDS(file=file.path("results", paste0(save_tag, gse_id, "_model.rds")))
    return(rf.model) 
}
