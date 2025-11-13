
library(ranger)

source("common.R")
options(ranger.num.threads = 8)
    
getFilteredBetasBasedOnImportanceParallel <- function (betas, anno) {

    N_CORES <- 8
    N_TREES <- 10000
    N_TOP_BETAS <- 100000 # top variably methylated probes
    N_TOP_FEATURES <- 10000 # top feature importance

    y <- as.factor(anno$`methylation class:ch1`)
    betas <- betas[,order(-apply(betas,2,sd))[1:N_TOP_BETAS]]
    
    rf.varsel <- ranger(
            x=betas,
            y=y,
            num.trees=N_TREES,
            num.threads=N_CORES,
            importance="permutation")

    # get permutation variable importance
    or <- order(rf.varsel$variable.importance, decreasing=TRUE)

    # reduce data matrix
    betas <- betas[,or[1:N_TOP_FEATURES]]

    return(betas)
}

getRandomForestModel <- function (betas, anno) {
    
    N_CORES <- 8
    N_TREES <- 10000

    y <- as.factor(anno$`methylation class:ch1`)
    
    rf.model <- ranger(
            x=betas,
            y=y,
            probability=TRUE,
            num.trees=N_TREES,
            num.threads=N_CORES,
            keep.inbag=TRUE,
            verbose=TRUE)

    return(rf.model)
}

do.getFilteredBetasBasedOnImportance <- function(gse_id, save_tag = "") {

    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID))

    betas <- loadSavedBetas(gse_id, save_tag)
    anno <- loadSavedAnno(gse_id, save_tag)

    seed <- 180314
    set.seed(seed,kind ="L'Ecuyer-CMRG")

    betas <- getFilteredBetasBasedOnImportanceParallel(betas, anno)

    saveRDS(betas, file=file.path("results", paste0(save_tag, gse_id,"_filtered_betas.rds")))
}

do.getRandomForestModel <- function(gse_id, save_tag = "") {
    
    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID))

    betas <- loadSavedFilteredBetas(gse_id, save_tag)
    anno <- loadSavedAnno(gse_id, save_tag)

    seed <- 180314
    set.seed(seed,kind ="L'Ecuyer-CMRG")

    rf.model <- getRandomForestModel(betas, anno)

    saveRDS(rf.model, file=file.path("results", paste0(save_tag, gse_id,"_model.rds")))
}