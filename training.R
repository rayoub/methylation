
library(doParallel)
library(foreach)
library(randomForest)

source("common.R")
    
N_CORES <- 10
N_TREES <- 10000
N_TOP_BETAS <- 30000 # top variably methylated probes
N_TOP_FEATURES <- 10000 # top feature importance

getFilteredBetasBasedOnImportance <- function (betas, anno) {

    y <- as.factor(anno$`methylation class:ch1`)
    betas <- betas[,order(-apply(betas,2,sd))[1:N_TOP_BETAS]]

    t1 <- Sys.time()

    cl <- makeCluster(N_CORES)
    registerDoParallel(cl)

    rf.varsel <- foreach(i = 1:N_CORES, .packages="randomForest", .combine=randomForest::combine) %dopar% {

        rf <- randomForest(
                betas,
                y=y,
                ntree=ceiling(N_TREES/N_CORES),
                sampsize=rep(min(table(y)),length(table(y))),
                norm.votes=FALSE,
                importance=TRUE)
        
        return(rf)
    }

    stopCluster(cl)

    t2 <- Sys.time()

    td <- t2 - t1
    
    print(td)

    # get permutation variable importance
    imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]

    # reduce data matrix
    or <- order(imp.meandecrease, decreasing=T)
    betas <- betas[,or[1:N_TOP_FEATURES]]

    return(betas)
}

getRandomForestModel <- function (betas, anno) {

    y <- as.factor(anno$`methylation class:ch1`)

    rf.model <- randomForest(
                    betas,
                    y,
                    ntree=N_TREES,
                    sampsize=rep(min(table(y)),length(table(y))),
                    proximity=TRUE,
                    oob.prox=TRUE,
                    importance=TRUE,
                    keep.inbag=TRUE,
                    do.trace=TRUE)

    return(rf.model)
}

do.getFilteredBetasBasedOnImportance <- function(gse_id, save_tag = "") {

    betas <- loadSavedBetas(gse_id, save_tag)
    anno <- loadSavedAnno(gse_id, save_tag)

    seed <- 180314
    set.seed(seed,kind ="L'Ecuyer-CMRG")

    betas <- getFilteredBetasBasedOnImportanceParallel(betas, anno)

    saveRDS(betas, file=file.path("results", paste0(save_tag, gse_id,"_filtered_betas.rds")))
}

do.getRandomForestModel <- function(gse_id, save_tag = "") {

    betas <- loadSavedFilteredBetas(gse_id, save_tag)
    anno <- loadSavedAnno(gse_id, save_tag)

    seed <- 180314
    set.seed(seed,kind ="L'Ecuyer-CMRG")

    rf.model <- getRandomForestModel(betas, anno)

    saveRDS(rf.model, file=file.path("results", paste0(save_tag, gse_id,"_model.rds")))
}