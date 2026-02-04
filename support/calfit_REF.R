
library(minfi)
library(limma)
library(glmnet)
library(here)

source(here("R","loading.R"))
source(here("R","preprocessing.R"))
source(here("R","prediction.R"))
source(here("R","training.R"))

FOLDS <- 3

make.folds <- function(y, cv.fold = 3) {
	n <- length(y)
	nlvl <- table(y)
	idx <- numeric(n)
	folds <- list()
	for (i in 1:length(nlvl)) {
		idx[which(y == levels(y)[i])] <- sample(rep(
			1:cv.fold,
			length = nlvl[i]
		))
	}
	for (i in 1:cv.fold) {
		folds[[i]] <- list(train = which(idx != i), test = which(idx == i))
	}
	return(folds)
}

make.nestedfolds <- function(y, cv.fold = 3) {

    seed <- 180314
    set.seed(seed,kind ="L'Ecuyer-CMRG")

	nfolds <- list()
	folds <- make.folds(y, cv.fold)
	names(folds) <- paste0("outer", 1:length(folds))
	for (k in 1:length(folds)) {
		inner = make.folds(y[folds[[k]]$train], cv.fold)
		names(inner) <- paste0("inner", 1:length(folds))
		for (i in 1:length(inner)) {
			inner[[i]]$train <- folds[[k]]$train[inner[[i]]$train]
			inner[[i]]$test <- folds[[k]]$train[inner[[i]]$test]
		}
		nfolds[[k]] <- list(folds[k], inner)
	}
	names(nfolds) <- paste0("outer", 1:length(nfolds))
	return(nfolds)
}

calculate.fold <- function(mset, y, material, fold) {

	# get train batch adjusted methys
	methy <- getMeth(mset[,fold$train])
	unmethy <- getUnmeth(mset[,fold$train])
	ba <- getBatchAdjustedMethyls(methy, unmethy, material[fold$train])
    
	# recalculate betas using the Illumina Genome Studio offset
	illumina_offset <- 100
	betas_train <- ba$methy.ba / (ba$methy.ba + ba$unmethy.ba + illumina_offset)
	betas_train <- as.data.frame(t(betas_train))

	# calculate batch adjust coefs to apply to the test data
	coefs <- getBatchEffectCoefs(methy, unmethy, ba$methy.ba, ba$unmethy.ba, material[fold$train])

	# get test batch adjusted methys
	methy <- getMeth(mset[,fold$test])
	unmethy <- getUnmeth(mset[,fold$test])
    ba <- applyBatchEffectCoefs(methy, unmethy, material[fold$test], coefs)

	# recalculate betas using the Illumina Genome Studio offset
	illumina_offset <- 100
	betas_test <- ba$methy.ba / (ba$methy.ba + ba$unmethy.ba + illumina_offset)
	betas_test <- as.data.frame(t(betas_test))

	# train random forest based on training data
	betas_train_filtered <- getFilteredBetasBasedOnImportance(betas_train, y[fold$train])
	rf_model <- getRandomForestModel(betas_train_filtered, y[fold$train])

	# now to calculate the scores
	scores <- predictSampleScores(rf_model, betas_test)

	return(scores)
}

calculate.folds <- function(gse_id) {

	mset <- loadSavedMset(gse_id)
	anno <- loadSavedAnno(gse_id)

	y <- as.factor(anno$`methylation class:ch1`)
	material <- as.factor(anno$`material:ch1`)

	nfolds <- make.nestedfolds(y, FOLDS)
  	save(nfolds,file=here("temp","nfolds.RData"))

	for (K in 1:FOLDS) {

		# calculate scores for the outer fold
		fold <- nfolds[[K]][[1]][[1]] 
		scores <- calculate.fold(mset, y, material, fold) 
		fname <- paste("CVfold", K, 0, "RData", sep = ".")
		save(scores, file = here("temp", fname))
	}
}

gse_id <- REF_GSE_ID
	
mset <- loadSavedMset(gse_id)
anno <- loadSavedAnno(gse_id)

y <- as.factor(anno$`methylation class:ch1`)
material <- as.factor(anno$`material:ch1`)

nfolds <- make.nestedfolds(y, FOLDS)

all_scores <- list()
for (K in 1:FOLDS) {

	# calculate scores for the outer folds
	fold <- nfolds[[K]][[1]][[1]] 
	all_scores[[K]] <- calculate.fold(mset, y, material, fold) 
}

# load all the outer test scores
idx <- list()
for (i in 1:length(nfolds)) {
	idx[[i]] <- nfolds[[i]][[1]][[1]]$test
}
idx <- unlist(idx) # collapse the list of indices used for testing

y <- anno$`methylation class:ch1`[idx] # these y's have to match the scores

# fit calibration model and save
suppressWarnings(
	calfit <- cv.glmnet(
		y = y,
		x = do.call(rbind, all_scores),
		family = "multinomial",
		type.measure = "mse",
		alpha = 0,
		nlambda = 100,
		lambda.min.ratio = 10^-6,
		parallel = TRUE
	)
)
saveRDS(calfit, file = file.path("results", paste0(gse_id, "_calfit.rds")))