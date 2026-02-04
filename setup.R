# **** this script should be run as administrator ****

# IMPORTANT NOTE: you may need to run this setup script twice.

# IMPORTANT NOTE: install rtools for Windows before running this script.

# ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

# BiocConductor version 3.22 is compatible with version 4.5.2+
BiocManager::install(version = "3.22")

# install all packages using BiocManager
BiocManager::install(c(
  "ade4"
  "GEOquery",
  "ggtext",
  "glmnet",
  "here",
  "lumi",
  "minfi",
  "openxlsx",
  "plotly",
  "ranger",
  "sesame",
  "splitstackshape",
  "tidyverse",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICv2manifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
))

# install conumee2 from source
install.packages("devtools")
devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")
devtools::install_github("badozor/mgmtstp27")

# to use sesame
library("sesame")
sesameDataCacheAll()


