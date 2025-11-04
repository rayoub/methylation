# **** this script should be run as administrator ****

# IMPORTANT NOTE: you may need to run this setup script twice.
# This is a known issue for some Bioc packages

# IMPORTANT NOTE: rtools and git-sdk conflict
# If you have git-sdk for Windows installed, the rtools installed with it
# was built with a incompatible tool chain for C/C++ compilation of
# BiocParallel and Rsamtools. To resolve, install rtools for Windows
# before running this script.

# ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

# BiocConductor version 3.21 is compatible with version 4.5.*
BiocManager::install(version = "3.21")

# install Bioconductor packages
BiocManager::install(c(
  # non Bioconductor packages
  #    "RSpectra",
  #    "Rtsne",
  #    "devtools",
  #    "abind",
  #    "statmod",

  # Bioconductor packages
  "minfi",
  "GEOquery",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "ggplot2",
  "randomForest",
  "doParallel"
))
