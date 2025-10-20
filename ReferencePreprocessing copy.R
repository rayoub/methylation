
library(GEOquery)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# create directories
data_dir = "data_raw"
data_ref = "data_raw/idats_reference"
data_val = "data_raw/idats_validation"
data_tmp <- "data_raw/tmp"

dir.create(data_dir, showWarnings = FALSE)
dir.create(data_ref, showWarnings = FALSE)
dir.create(data_val, showWarnings = FALSE)
dir.create(data_tmp, showWarnings = FALSE)

# ******************************************************************************
# *** DOWNLOAD IDAT FILES FOR GSE90496 REFERENCE DATA
# ******************************************************************************

# get series and sample list
gse_id <- "GSE90496"
gse <- getGEO(gse_id, GSEMatrix = FALSE) 
gsm_list <- GSMList(gse)

# get the IDAT supplementary file URLs from the samples
supp_urls <- list()
for (gsm_name in names(gsm_list)) {
  gsm <- gsm_list[[gsm_name]]
  sfiles <- Meta(gsm)$supplementary_file
  if (length(sfiles) > 0) supp_urls[[gsm_name]] <- sfiles
}

# get rid of dupes
unique_supp <- gsub("^ftp://", "https://", unique(unlist(supp_urls)))

# increase timeout so big IDATs don't fail
options(timeout = 600)   # 10 minutes per file; increase if needed

# function for downloading a batch
download_batch <- function(urls, batch_num, folder) {
  message("Starting batch ", batch_num, " with ", length(urls), " files...")
  for (url in urls) {
    fname <- file.path(folder, basename(url))
    if (!file.exists(fname)) {
      message("Downloading: ", basename(url))
      try(
        download.file(url, destfile = fname, mode = "wb", method = "curl"),
        silent = TRUE
      )
    }
    # If it's an archive, unpack
    if (grepl("\\.tar|\\.tgz|\\.zip", fname)) {
      try({
        if (grepl("\\.zip$", fname)) unzip(fname, exdir = "data_raw/idats")
        else untar(fname, exdir = "data_raw/idats")
      }, silent = TRUE)
    }
  }
  message("Batch ", batch_num, " done.")
}

# split into batches of 10
batch_size <- 10
num_batches <- 2 #ceiling(length(unique_supp) / batch_size)

# download batches - this will take a few hours
for (i in seq_len(num_batches)) {
  start <- (i - 1) * batch_size + 1
  end <- min(i * batch_size, length(unique_supp))
  batch_urls <- unique_supp[start:end]
  
  download_batch(batch_urls, i, data_ref)
  
  # pause between batches to avoid hammering NCBI
  Sys.sleep(5)
}

# unzip IDAT files using gzip utility
files <- list.files(data_ref, pattern = "idat.gz$", full.names = TRUE, recursive = TRUE)
for (f in files) {
  R.utils::gunzip(f, overwrite = TRUE)  # removes .gz, leaves .idat
}

# ******************************************************************************
# *** GET BETA FILE FROM IDAT FILES
# ******************************************************************************

# get a list of the download IDAT files
idat_files <- list.files(data_ref, pattern = "\\.(idat|IDAT)$", full.names = TRUE, recursive = TRUE)

# samples must have both green and red
samples <- sub("_(Grn|Red)\\.idat$", "", basename(idat_files))
counts <- table(samples)
missing <- names(counts[counts != 2])
if(length(missing) > 0){
  warning("Dropping samples missing Red/Grn pair: ", paste(missing, collapse=", "))
}

# get unique samples and exclude ones without both green and red
samples <- setdiff(unique(samples), missing)

# preprocess beta files 
batch_size <- 10    # 900
for (i in seq(1, length(samples), by = batch_size)) {
  idx <- i:min(i + batch_size - 1, length(samples))
  message("Processing batch ", i, " to ", idx[length(idx)])
  
  batch_samples <- samples[idx]
  batch_files <- idat_files[
    sub("_(Grn|Red)\\.idat$", "", basename(idat_files)) %in% batch_samples
  ]
  
  file.copy(batch_files, data_tmp, overwrite = TRUE)
  
  # read and normalize
  rgSet <- read.metharray.exp(data_tmp, recursive = TRUE)
  mset <- preprocessFunnorm(rgSet)
  beta <- getBeta(mset)
  
  # save batch
  saveRDS(beta, file = paste0(file.path(data_dir, "beta_batch_", i, ".rds"))
  
  # clean tmp
  unlink(data_tmp, recursive = TRUE)
  gc()  # free memory
}

# merge beta batches
batch_files <- list.files(data_dir, pattern="^beta_batch_.*\\.rds$", full.names=TRUE)
beta_full <- do.call(cbind, lapply(batch_files, readRDS))

# save the merged beta file
beta_full_file <- "beta_normalized_GSE90496_full.rds"
saveRDS(beta_full, file.path(data_dir, beta_full_file))
message("Full beta matrix saved: ", beta_full_file)

# ******************************************************************************
# *** PROCESS THE BETA FILE ACCORDING TO PAPER 
# ******************************************************************************

# load 450k annotation
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# align probes
common <- intersect(rownames(beta_full), rownames(ann))
beta_full <- beta_full[common, ]
ann <- ann[common, ]

# remove sex chromosome probes
keep <- !(ann$chr %in% c("chrX", "chrY"))

# remove probes with specific SNP (CpG or within 5bp)
if ("SNPs.137CommonSingle" %in% colnames(ann)) {
  snp_mask <- ann$SNPs.137CommonSingle != ""
} else {
  snp_mask <- rep(FALSE, nrow(ann))  # fallback
}
keep <- keep & !snp_mask

# remove cross-reactive probes (Zhou 2016 / Chen 2013)
xreactive_file <- "48639-non-specific-probes-Illumina450k.csv"
xreactive <- read.csv(xreactive_file, header = FALSE, stringsAsFactors = FALSE)[,1]
keep <- keep & !(rownames(ann) %in% xreactive)

# keep only probes common to 450K and EPIC
ann_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
probes_450k <- rownames(ann)
probes_EPIC <- rownames(ann_epic)
common_probes_list <- intersect(probes_450k, probes_EPIC)

keep <- keep & (rownames(ann) %in% common_probes_list)

# apply filter
beta_filtered <- beta_full[keep, ]

# save final filtered matrix
saveRDS(beta_filtered, "data_raw/beta_filtered_GSE90496.rds")
message("Filtering complete. Final beta matrix: ", nrow(beta_filtered), " probes x ", ncol(beta_filtered), " samples")
