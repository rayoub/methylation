
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(here)

anno_450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_epicv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
probes_450k <- rownames(anno_450k)
probes_epic <- rownames(anno_epic)
probes_epicv2 <- rownames(anno_epicv2)

probes_epicv2 <- sub("_.*", "", probes_epicv2)

valid_probes1 <- intersect(probes_450k, probes_epic)
valid_probes2 <- intersect(valid_probes1, probes_epicv2)

writeLines(valid_probes2, con=here("probes", "common_450k_epic_epicv2.txt"))


