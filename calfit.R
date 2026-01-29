
library(minfi)
library(limma)
library(glmnet)

source("common.R")
source("preprocess.R")
source("train.R")
source("predict.R")

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
	scores <- do.predict(rf_model, betas_test)

	return(scores)
}

calculate.folds <- function(gse_id) {

	FOLDS <- 3

	mset <- loadSavedMset(gse_id)
	anno <- loadSavedAnno(gse_id)

	y <- as.factor(anno$`methylation class:ch1`)
	material <- as.factor(anno$`material:ch1`)

	dir.create("temp", showWarnings = FALSE)

	nfolds <- make.nestedfolds(y, FOLDS)
  	save(nfolds,file=file.path("temp","nfolds.RData"))

	for (K in 1:FOLDS) {

		# calculate scores for the outer fold
		fold <- nfolds[[K]][[1]][[1]] 
		scores <- calculate.fold(mset, y, material, fold) 
		fname <- paste("CVfold", K, 0, "RData", sep = ".")
		save(scores, file = file.path("temp", fname))

# THIS IS NOT UNUSED TO CREATE FINAL CALFIT MODEL
#		# calculate scores for the inner fold
#		for (k in 1:FOLDS) {
#		
#			fold <- nfolds[[K]][[2]][[k]] 
#			scores <- calculate.fold(mset, y, material, fold) 
#			fname <- paste("CVfold", K, k, "RData", sep = ".")
#			save(scores, file = file.path("temp", fname))
#		}
	}
}

fit.calibration <- function (gse_id) {

	# calculate scores to build calibration model
	calculate.folds(gse_id)

	# load the fold definitions - nfolds variable
	load(file.path("temp", "nfolds.RData"))

	# load all the outer test scores
	s <- list()
	idx <- list()
	for (i in 1:length(nfolds)) {
		fname <- paste0("CVfold.", i, ".", 0, ".RData")
		load(file.path("temp", fname)) # this is where the scores variable is coming from
		s[[i]] <- scores
		idx[[i]] <- nfolds[[i]][[1]][[1]]$test
	}
	s <- do.call(rbind, s)
	idx <- unlist(idx) # collapse the list of indices used for testing

	anno <- loadSavedAnno(gse_id)
	y <- anno$`methylation class:ch1`[idx] # these y's have to match the scores

	# fit calibration model and save
	suppressWarnings(
		calfit <- cv.glmnet(
			y = y,
			x = s,
			family = "multinomial",
			type.measure = "mse",
			alpha = 0,
			nlambda = 100,
			lambda.min.ratio = 10^-6,
			parallel = TRUE
		)
	)
	saveRDS(calfit, file = file.path("results", paste0(gse_id, "_calfit.rds")))
}