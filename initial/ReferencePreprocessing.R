# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")

# Set repositories over HTTP (This step may not be necessary if you have contacted IT)
options(repos = c(
  CRAN    = "http://cran.us.r-project.org",
  BioCsoft = "http://bioconductor.org/packages/3.21/bioc",
  BioCann  = "http://bioconductor.org/packages/3.21/data/annotation",
  BioCexp  = "http://bioconductor.org/packages/3.21/data/experiment"
))

# Install CRAN and Bioconductor packages (I had to contact IT for a firewall bypass. These only need to be installed
# once, then you can just load the libraries).
install.packages("minfi")
install.packages("GEOquery")
install.packages("IlluminaHumanMethylation450kmanifest")
install.packages("IlluminaHumanMethylationEPICmanifest")
install.packages("IlluminaHumanMethylation450kanno.ilmn12.hg19")
install.packages("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
install.packages("matrixStats")
install.packages("RSpectra")
install.packages("Rtsne")
install.packages("devtools")
install.packages("abind")
install.packages("SummarizedExperiment")
install.packages("statmod")
BiocManager::install("IlluminaHumanMethylation450kmanifest")

## Start here after initial setup
# Load packages
library(GEOquery)
library(minfi)
library(matrixStats)
library(data.table)
library(RSpectra)
library(Rtsne)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)



# === START: 01 Download and Preprocess ===
# Purpose: download GSE90496 (reference set) where possible, load IDATs, normalize, and filter probes

# 1) Download GEO series metadata (GSE90496 - reference/training set)
gse_id <- "GSE90496"
gse <- getGEO(gse_id, GSEMatrix = FALSE)  # metadata;
gsm_list <- GSMList(gse)

# 2) Inspect supplementary files for raw IDATs
supp_urls <- list()
for (gsm_name in names(gsm_list)) {
  gsm <- gsm_list[[gsm_name]]
  sfiles <- Meta(gsm)$supplementary_file
  if (length(sfiles) > 0) supp_urls[[gsm_name]] <- sfiles
}
length(supp_urls)

# 3) Try to download series matrix / raw files fallback
# (A) If IDAT bundles exist, download and untar into local directory - data_raw/idats/
dir.create("data_raw", showWarnings = FALSE)
dir.create("data_raw/idats", showWarnings = FALSE)

# Increase timeout so big IDATs don't fail
options(timeout = 600)   # 10 minutes per file; increase if needed

# Collapse to unique URLs
unique_supp <- gsub("^ftp://", "https://", unique(unlist(supp_urls)))

# Function to download a batch of files
download_batch <- function(urls, batch_num) {
  message("Starting batch ", batch_num, " with ", length(urls), " files...")
  for (url in urls) {
    fname <- file.path("data_raw/idats", basename(url))
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

# Split into batches of 10
batch_size <- 10
num_batches <- ceiling(length(unique_supp) / batch_size)

for (i in seq_len(num_batches)) {
  start <- (i - 1) * batch_size + 1
  end <- min(i * batch_size, length(unique_supp))
  batch_urls <- unique_supp[start:end]
  
  download_batch(batch_urls, i)
  
  # Optional: pause between batches to avoid hammering NCBI
  Sys.sleep(5)
}
### This is going to take a few hours per dataset

# (B) Only if no raw IDATs downloaded, attempt to get series matrix processed file
if (length(list.files("data_raw/idats")) == 0) {
  gse_matrix <- try(getGEO(gse_id, GSEMatrix = TRUE), silent = TRUE)
  if (!inherits(gse_matrix, "try-error")) {
    saveRDS(gse_matrix, file = file.path("data_raw", paste0(gse_id, "_matrix.RDS")))
    message("Series matrix downloaded and saved. If you need raw IDATs, consider contacting data submitters.")
  } else {
    message("No series matrix available via getGEO(). May need manual download.")
  }
}

# 4) If IDAT files are present, direct R to folder and unzip
idat_dir <- "data_raw/idats"
files <- list.files(idat_dir, pattern = "idat.gz$", full.names = TRUE, recursive = TRUE)
for (f in files) {
  R.utils::gunzip(f, overwrite = TRUE)  # removes .gz, leaves .idat
}

#5 Process IDATs with minfi, noramlize, and extract Beta values
library(minfi)

# List all .idat files
idat_files <- list.files(idat_dir, pattern = "\\.(idat|IDAT)$", full.names = TRUE, recursive = TRUE)

# Check samples have both Grn + Red
samples <- sub("_(Grn|Red)\\.idat$", "", basename(idat_files))
counts <- table(samples)
missing <- names(counts[counts != 2])
if(length(missing) > 0){
  warning("Dropping samples missing Red/Grn pair: ", paste(missing, collapse=", "))
}
samples <- setdiff(unique(samples), missing)

batch_size <- 900

for (i in seq(1, length(samples), by = batch_size)) {
  idx <- i:min(i + batch_size - 1, length(samples))
  message("Processing batch ", i, " to ", idx[length(idx)])
  
  batch_samples <- samples[idx]
  
  # Copy only batch IDATs into a temp folder
  tmp_dir <- "data_raw/tmp"
  dir.create(tmp_dir, showWarnings = FALSE)
  
  batch_files <- idat_files[
    sub("_(Grn|Red)\\.idat$", "", basename(idat_files)) %in% batch_samples
  ]
  
  file.copy(batch_files, tmp_dir, overwrite = TRUE)
  
  # Read and normalize
  rgSet <- read.metharray.exp(tmp_dir, recursive = TRUE)
  mset <- preprocessFunnorm(rgSet)
  beta <- getBeta(mset)
  
  # Save batch
  saveRDS(beta, file = paste0("data_raw/beta_batch_", i, ".RDS"))
  
  # Clean tmp
  unlink(tmp_dir, recursive = TRUE)
  gc()  # free memory
}

#6 Merge batches 
batch_files <- list.files("data_raw", pattern="^beta_batch_.*\\.RDS$", full.names=TRUE)
beta_full <- do.call(cbind, lapply(batch_files, readRDS))
dim(beta_full)  # probes x all samples
saveRDS(beta_full, "data_raw/beta_normalized_GSE90496_full.RDS")
message("Full beta matrix saved: beta_normalized_GSE90496_full.RDS")


#7 Probe filtering
# Load beta matrix
beta_full <- readRDS("data_raw/beta_normalized_GSE90496_full.RDS")

# Load 450k annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Align probes
common <- intersect(rownames(beta_full), rownames(ann))
beta_full <- beta_full[common, ]
ann <- ann[common, ]

# 1) Remove sex chromosome probes
keep <- !(ann$chr %in% c("chrX", "chrY"))

# 2) Remove probes with SNPs (CpG or within 5bp)
if ("SNPs.137CommonSingle" %in% colnames(ann)) {
  snp_mask <- ann$SNPs.137CommonSingle != ""
} else {
  snp_mask <- rep(FALSE, nrow(ann))  # fallback
}
keep <- keep & !snp_mask

# 3) Remove cross-reactive probes (Zhou 2016 / Chen 2013)
# Download cross-reactive probe list if not already present - place locally
xreactive_file <- "C:/Users/a27668/Downloads/48639-non-specific-probes-Illumina450k.csv"
xreactive <- read.csv(xreactive_file, header = FALSE, stringsAsFactors = FALSE)[,1]
keep <- keep & !(rownames(ann) %in% xreactive)

# 4) Keep only probes common to 450K and EPIC
# Load EPIC annotation
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly=TRUE)) {
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
}
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

probes_450k <- rownames(ann)
probes_EPIC <- rownames(ann_epic)
common_probes_list <- intersect(probes_450k, probes_EPIC)

keep <- keep & (rownames(ann) %in% common_probes_list)

# Apply filter
beta_filtered <- beta_full[keep, ]

# Save final filtered matrix
saveRDS(beta_filtered, "data_raw/beta_filtered_GSE90496.RDS")
message("Filtering complete. Final beta matrix: ", 
        nrow(beta_filtered), " probes x ", ncol(beta_filtered), " samples")

# === END: 01_download_and_preprocess.R ===