
library(conumee2)
library(here)

load(here("cnv","CNanalysis4_conumee_ANNO.vh20150715.RData"))

cnv_anno <- CNV.create_anno(array_type = "EPICv2", detail_regions = anno@detail, exclude_regions = anno@exclude)

saveRDS(cnv_anno, file=here("results", "CNV_ANNO_EPICv2.rds"))
