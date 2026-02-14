
library(here)

source(here("R","loading.R"))
source(here("R","training.R"))
    
betas <- loadSavedBetas(REF_GSE_ID)
anno <- loadSavedAnno(REF_GSE_ID)

y <- as.factor(anno$`methylation class:ch1`)

betas <- getFilteredBetasBasedOnImportance(betas, y)

rf_model <- getRandomForestModel(betas, y)

saveRDS(rf_model, file=file.path("results", paste0(REF_GSE_ID,"_model.rds")))