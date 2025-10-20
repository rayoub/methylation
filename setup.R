
# **** this script should be run as administrator ****

# ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")

# BiocConductor version 3.21 is compatible with version 4.5.0 
# These versions must match. 
# If you have to, install an older or newer version of R and update your PATH environment variable
BiocManager::install(version = "3.21")

# install Bioconductor packages
# Since I had git-sdk for Windows on my machine the rtools installed with it was built with a different tool chain 
# and was resulting in installation problems with BiocParallel and Rsamtools. I had to install rtools for Windows and then it worked
BiocManager::install(c(

    # non Bioconductor packages
#    "matrixStats",
#    "RSpectra",
#    "Rtsne",
#    "devtools",
#    "abind",
#    "statmod",

    # Bioconductor packages
    "minfi",
    "GEOquery",
    "SummarizedExperiment",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
))
