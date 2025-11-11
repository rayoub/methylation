
library(minfi)
library(GEOquery)
library(limma)

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

    # probe filtering
    amb.filter <- read.table(file.path("probes","amb_3965probes.vh20151030.txt"), header=F)
    epic.filter <- read.table(file.path("probes","epicV1B2_32260probes.vh20160325.txt"), header=F)
    snp.filter <- read.table(file.path("probes","snp_7998probes.vh20151030.txt"), header=F)
    xy.filter <- read.table(file.path("probes","xy_11551probes.vh20151030.txt"), header=F)
    rs.filter <- grep("rs",rownames(mset))
    ch.filter <- grep("ch",rownames(mset))

    # filter CpG probes
    remove <- unique(c(match(amb.filter[,1], rownames(mset)),
                    match(epic.filter[,1], rownames(mset)),
                    match(snp.filter[,1], rownames(mset)),
                    match(xy.filter[,1], rownames(mset)),
                    rs.filter,
                    ch.filter))

    mset <- mset[-remove,]

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

preprocessSamples <- function (gse_id, save_tag = "") {

    gse_id <- match.arg(gse_id, c(REF_GSE_ID, VAL_GSE_ID))

    # annotations
    message("getting sample annotations ... ", Sys.time())
    anno <- getSampleAnnotations(gse_id)
    material <- anno$`material:ch1`
    saveRDS(anno, file=file.path("results", paste0(save_tag, gse_id,"_anno.rds")))

    # when testing apply the filter below when getting the anno
    #anno <- anno[grep("ATRT", anno$`methylation class:ch1`),]

    # get basenames for idat files
    data_dir = file.path("data", gse_id)
    basenames <- file.path(data_dir, gsub("_Grn.*", "", gsub(".*suppl/", "", anno$supplementary_file)))

    # read idat files into an RGSet
    rgset <- read.metharray(basenames, verbose=TRUE)

    # preprocess RGSet and convert to a MethylSet
    mset <- MNPpreprocessIllumina(rgset)

    # apply probe filtering to the MethylSet
    mset <- getFilteredMethylSet(mset)
    saveRDS(mset, file=file.path("results", paste0(save_tag, gse_id,"_mset.rds")))

    methy <- getMeth(mset)
    unmethy <- getUnmeth(mset)

    if (gse_id == REF_GSE_ID) {
    
        # get batch adjusted methyls
        ba <- getBatchAdjustedMethyls(methy, unmethy, material)

        # save batch effect coefs for future diagnostic samples
        coefs <- getBatchEffectCoefs(methy, unmethy, ba$methy.ba, ba$unmethy.ba, material)
        saveRDS(coefs, file=file.path("results", paste0(save_tag, gse_id,"_coefs.rds")))
    }
    else { # gse_id == VAL_GSE_ID

        # TODO: use coefs to make these batch adjustments
        ba <- list(methy.ba = methy, unmethy.ba = unmethy)
    }

    # recalculate betas using the Illumina Genome Studio offset
    illumina_offset <- 100
    betas <- ba$methy.ba / (ba$methy.ba + ba$unmethy.ba + illumina_offset)
    betas <- as.data.frame(t(betas))
    saveRDS(betas, file=file.path("results", paste0(save_tag, gse_id, "_betas.rds")))

    message("preprocessing finished ... ", Sys.time())
}
