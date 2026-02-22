
library(here)

source(here::here("R","files.R"))
source(here::here("R","training.R"))
    
betas <- loadGeoData(REF_GSE_ID, "betas")
anno <- loadGeoData(REF_GSE_ID, "anno")

y <- as.factor(anno$`methylation class:ch1`)

betas <- getFilteredBetasBasedOnImportance(betas, y)

rf_model <- getRandomForestModel(betas, y)

saveGeoData(REF_GSE_ID, "model", rf_model)