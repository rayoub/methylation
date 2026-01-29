
library(ranger)

source("common.R")
options(ranger.num.threads = 8)
    
getFilteredBetasBasedOnImportance <- function (betas, y) {

    N_CORES <- 8
    N_TREES <- 10000
    N_TOP_BETAS <- 100000 # top variably methylated probes
    N_TOP_FEATURES <- 30000 # top feature importance
    
    seed <- 180314
    set.seed(seed,kind ="L'Ecuyer-CMRG")

    betas <- betas[,order(-apply(betas,2,sd))[1:N_TOP_BETAS]]
    
    rf_varsel <- ranger(
            x=betas,
            y=y,
            num.trees=N_TREES,
            num.threads=N_CORES,
            importance="permutation")

    # get permutation variable importance
    or <- order(rf_varsel$variable.importance, decreasing=TRUE)

    # reduce data matrix
    betas <- betas[,or[1:N_TOP_FEATURES]]

    return(betas)
}

getRandomForestModel <- function (betas, y) {
    
    N_CORES <- 8
    N_TREES <- 10000
    
    seed <- 180314
    set.seed(seed,kind ="L'Ecuyer-CMRG")
    
    rf_model <- ranger(
            x=betas,
            y=y,
            probability=TRUE,
            num.trees=N_TREES,
            num.threads=N_CORES,
            keep.inbag=TRUE,
            verbose=TRUE)

    return(rf_model)
}

do.getFilteredBetasBasedOnImportance <- function() {

    gse_id <- REF_GSE_ID

    betas <- loadSavedBetas(gse_id)
    anno <- loadSavedAnno(gse_id)
    
    y <- as.factor(anno$`methylation class:ch1`)
    
    betas <- getFilteredBetasBasedOnImportance(betas, y)

    saveRDS(betas, file=file.path("results", paste0(gse_id,"_filtered_betas.rds")))
}

do.getRandomForestModel <- function() {

    gse_id <- REF_GSE_ID

    betas <- loadSavedFilteredBetas(gse_id)
    anno <- loadSavedAnno(gse_id)
    
    y <- as.factor(anno$`methylation class:ch1`)

    rf_model <- getRandomForestModel(betas, y)

    saveRDS(rf_model, file=file.path("results", paste0(gse_id,"_model.rds")))
}