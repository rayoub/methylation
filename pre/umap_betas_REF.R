
library(here)

source(here::here("R", "constants.R"))
source(here::here("R", "loading.R"))

N_TOP_BETAS <- 32000 # top variably methylated probes

# load betas and filter top 32,000 variably methylated probes
betas <- loadSavedBetas(REF_GSE_ID)
filtered_betas <- betas[, order(-apply(betas, 2, sd))[1:N_TOP_BETAS]]

saveRDS(filtered_betas, file=here::here("results", paste0(REF_GSE_ID, "_umap_betas.rds")))


