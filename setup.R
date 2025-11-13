# **** this script should be run as administrator ****

# IMPORTANT NOTE: you may need to run this setup script twice.

# IMPORTANT NOTE: install rtools for Windows before running this script.

# ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

# BiocConductor version 3.21 is compatible with version 4.5.*
BiocManager::install(version = "3.21")

# install all packages using BiocManager
BiocManager::install(c(
  "minfi",
  "GEOquery",
  "randomForest",
  "dplyr",
  "doParallel",
  "foreach",
  "glmnet",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
))

# create the results directory in the project root
dir.create("results", showWarnings = FALSE)