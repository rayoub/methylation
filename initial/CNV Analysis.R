# --- CNV Analysis (Uses GSE109379 validation set, already downloaded locally on my computer.  
# Use GetGEO if needed for others or direct file path accordingly ---
library(minfi)
library(conumee)

# Load reference data - Download from suplement of Capper paper or Github - https://github.com/mwsill/mnp_training
load("CNV_data/CNanalysis4_conumee_ANNO.vh20150715.RData")
load("CNV_data/CNanalysis4_conumee_REF-M.vh20150715.RData")
load("CNV_data/CNanalysis4_conumee_REF-F.vh20150715.RData")

# Load your IDATs
idat_dir <- "data_raw/validation_idats"
RGset <- read.metharray.exp(idat_dir, verbose = TRUE)
RGset <- RGset[, 1] # Uses first sample only, change as needed/wanted

# Preprocess (no normalization needed for CNV analysis)
Mset <- preprocessRaw(RGset)

# Run CNV analysis for each sample
  cndata <- CNV.load(Mset)
  
  # Choose reference (adjust refM.data vs refF.data for sex)
  x <- CNV.fit(cndata, refF.data, annoXY)
  # Double check sex
  
  x <- CNV.bin(x)
  x <- CNV.detail(x)
  x <- CNV.segment(x)
  
  # Plot results
  CNV.genomeplot(x, chrY=TRUE)
  
  #If you want PDF generated
  pdf("CNV_data/CNV Plots/GSM2940725.pdf", width = 12, height = 6)
  conumee::CNV.genomeplot(x, chrY = TRUE)
  dev.off()
  