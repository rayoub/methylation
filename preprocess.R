
library(minfi)
library(GEOquery)
library(limma)
library(dplyr)

source("common.R")
source("MNPpreprocess.R")

getSampleAnnotations <- function (gse_id) {

    # get sample annotation from GEO
    gse <- getGEO(gse_id, GSEMatrix=TRUE, getGPL=FALSE)
    mat <- paste0(gse_id, "_series_matrix.txt.gz")
    anno <- pData(gse[[mat]])

    return(anno)
}

getFilteredMethylSet <- function (mset) {

    # only consider probes present in 450k, epic, and epicv2 bead chips
    common.probes <- read.table(file.path("probes","common_450k_epic_epicv2.txt"), header=FALSE)
    mset <- mset[na.omit(match(common.probes[[1]],rownames(mset))),]

    # CpG probes 
    amb.filter <- read.table(file.path("probes","amb_3965probes.vh20151030.txt"), header=FALSE)
    epic.filter <- read.table(file.path("probes","epicV1B2_32260probes.vh20160325.txt"), header=FALSE)
    snp.filter <- read.table(file.path("probes","snp_7998probes.vh20151030.txt"), header=FALSE)
    xy.filter <- read.table(file.path("probes","xy_11551probes.vh20151030.txt"), header=FALSE)
    rs.filter <- grep("rs",rownames(mset))
    ch.filter <- grep("ch",rownames(mset))

    # additional filter CpG probes
    remove <- unique(c(match(amb.filter[,1], rownames(mset)),
                    match(epic.filter[,1], rownames(mset)),
                    match(snp.filter[,1], rownames(mset)),
                    match(xy.filter[,1], rownames(mset)),
                    rs.filter,
                    ch.filter))

    mset <- mset[-na.omit(remove),]

    return(mset)
}

getBatchAdjustedMethyls <- function (methy, unmethy, material) {

    # batch based on material
    batch <- ifelse(material == "FFPE", 2, 1)

    # remove batch effects by linear model (the exponential reverses the log2)
    methy.ba <- 2^removeBatchEffect(log2(methy + 1), batch)
    unmethy.ba <- 2^removeBatchEffect(log2(unmethy + 1), batch)

    return(list(methy.ba = methy.ba, unmethy.ba = unmethy.ba))
}

getBatchEffectCoefs <- function (methy, unmethy, methy.ba, unmethy.ba, material) {
    
    # batch based on material
    batch <- ifelse(material == "FFPE", 2, 1)

    # this can probably be improved upon: 
    # instead of picking the first column of Frozen or FFPE, try picking the medoid

    # extract effects to adjust diagnostic samples
    s.frozen <- min(which(batch == 1)) 
    s.ffpe <- min(which(batch == 2))
    methy.coef <- unmethy.coef <- list()
    methy.coef[["Frozen"]] <- log2(methy.ba[, s.frozen]) - log2(methy[, s.frozen] +1)
    methy.coef[["FFPE"]] <- log2(methy.ba[, s.ffpe]) - log2(methy[, s.ffpe] +1)
    unmethy.coef[["Frozen"]] <- log2(unmethy.ba[, s.frozen]) - log2(unmethy[, s.frozen] +1)
    unmethy.coef[["FFPE"]] <- log2(unmethy.ba[, s.ffpe]) - log2(unmethy[, s.ffpe] +1)

    return(list(methy.coef = methy.coef, unmethy.coef = unmethy.coef))
}

applyBatchEffectCoefs <- function (methy, unmethy, material, coefs) {

    # repeat material if length 1 the number of samples
    if (length(material) == 1) {
        material = rep(material, ncol(methy))
    }

    # batch adjust 
    # WARNING: this code depends on Frozen is the first element in the *.coef list
    methy.b <- log2(methy + 1) + matrix(unlist(coefs$methy.coef[match(material,names(coefs$methy.coef))]),ncol=length(material)) 
    unmethy.b <- log2(unmethy + 1) + matrix(unlist(coefs$unmethy.coef[match(material,names(coefs$unmethy.coef))]),ncol=length(material))
    methy.b[methy.b < 0] <- 0
    unmethy.b[unmethy.b < 0] <- 0
    methy.ba <- 2^methy.b
    unmethy.ba <- 2^unmethy.b

    return(list(methy.ba = methy.ba, unmethy.ba = unmethy.ba))
}

do.preprocessGEO <- function (gse_id) {

    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID, NM_GSE_ID))

    # annotations
    message("getting sample annotations ... ", Sys.time())
    anno <- getSampleAnnotations(gse_id)
    saveRDS(anno, file=file.path("results", paste0(gse_id,"_anno.rds")))

    # get material and normalize
    material <- if ("material:ch1" %in% names(anno)) anno$`material:ch1` else anno$`sample type:ch1`
    material <- case_when(
        material == "DNA_FFPE" ~ "FFPE",
        material == "DNA_KRYO" ~ "Frozen",
        TRUE ~ material
    )

    # get basenames for idat files
    data_dir = file.path("data", gse_id)
    basenames <- unique(file.path(data_dir, gsub("_Grn.*", "", gsub("_Red.*", "", list.files(path = data_dir, pattern = "*.idat")))))

    message("reading meth arrays ... ", Sys.time())

    # read idat files into an RGSet MAYBE THE PROBE NAMES ARE HERE ALREADY AT THIS POINT
    rgset <- read.metharray(basenames, verbose=TRUE)

    message("preprocessing meth arrays ... ", Sys.time())

    # preprocess RGSet and convert to a MethylSet
    mset <- MNPpreprocessIllumina(rgset)
    
    message("applying probe filters ... ", Sys.time())

    # apply probe filtering to the MethylSet
    mset <- getFilteredMethylSet(mset)
    saveRDS(mset, file=file.path("results", paste0(gse_id,"_mset.rds")))

    methy <- getMeth(mset)
    unmethy <- getUnmeth(mset)
        
    message("performing batch adjustments ... ", Sys.time())

    if (gse_id == REF_GSE_ID) {
    
        # get batch adjusted methyls
        ba <- getBatchAdjustedMethyls(methy, unmethy, material)

        # save batch effect coefs for future diagnostic samples
        coefs <- getBatchEffectCoefs(methy, unmethy, ba$methy.ba, ba$unmethy.ba, material)
        saveRDS(coefs, file=file.path("results", paste0(gse_id,"_coefs.rds")))
    }
    else { # gse_id == VAL_GSE_ID || NM_GSE_ID

        # apply batch effects from saved coefs
        coefs <- loadSavedCoefs(REF_GSE_ID)
        ba <- applyBatchEffectCoefs(methy, unmethy, material, coefs)
    }

    # recalculate betas using the Illumina Genome Studio offset
    illumina_offset <- 100
    betas <- ba$methy.ba / (ba$methy.ba + ba$unmethy.ba + illumina_offset)
    betas <- as.data.frame(t(betas))
    saveRDS(betas, file=file.path("results", paste0(gse_id, "_betas.rds")))

    message("preprocessing finished ... ", Sys.time())
}

do.preprocessDIAG <- function (diag_id, material) {

    # get basenames for idat files
    data_dir = file.path("data", "diagnostic", diag_id)
    basenames <- unique(file.path(data_dir, gsub("_Grn.*", "", gsub("_Red.*", "", list.files(path = data_dir, pattern = "*.idat")))))

    message("reading meth arrays ... ", Sys.time())

    # read idat files into an RGSet
    rgset <- read.metharray(basenames, verbose=TRUE)

    message("preprocessing meth arrays ... ", Sys.time())

    # preprocess RGSet and convert to a MethylSet
    mset <- MNPpreprocessIllumina(rgset)
    rownames(mset) <- sub("_.*", "", rownames(mset))

    message("applying probe filters ... ", Sys.time())

    # apply probe filtering to the MethylSet
    mset <- getFilteredMethylSet(mset)

    methy <- getMeth(mset)
    unmethy <- getUnmeth(mset)
        
    message("performing batch adjustments ... ", Sys.time())

    # apply batch effects from saved coefs
    coefs <- loadSavedCoefs(REF_GSE_ID)
    ba <- applyBatchEffectCoefs(methy, unmethy, material, coefs)

    # recalculate betas using the Illumina Genome Studio offset
    illumina_offset <- 100
    betas <- ba$methy.ba / (ba$methy.ba + ba$unmethy.ba + illumina_offset)
    betas <- as.data.frame(t(betas))
    saveRDS(betas, file=file.path("results", paste0(diag_id, "_betas.rds")))

    message("preprocessing finished ... ", Sys.time())
}
