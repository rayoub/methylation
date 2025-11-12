

source("common.R")

library(dplyr)
library(minfi)

  mset <- loadSavedMset(VAL_GSE_ID) 

  methy <- getMeth(mset)
  unmethy <- getUnmeth(mset)

  anno <- loadSavedAnno(VAL_GSE_ID)
  material <- anno$`material:ch1`
  material <- case_when(
    material == "DNA_FFPE" ~ "FFPE",
    material == "DNA_KRYO" ~ "Frozen",
    TRUE ~ material
  )


# FUNCTION BEGINS HERE ARGS: methy unmethy material
  

  coefs <- loadSavedCoefs(REF_GSE_ID)

  # coefs is a nest list of vectors - the vectors are 
  # coefs$methy.coef$Frozen and FFPE
  # coefs$unmethy.coef$Frozen or FFPE

  # this code depends on Frozen is the first element in the *.coef list
  methy.b <- log2(methy + 1) + matrix(unlist(coefs$methy.coef[match(material,names(coefs$methy.coef))]),ncol=length(material)) 
  unmethy.b <- log2(unmethy + 1) + matrix(unlist(coefs$unmethy.coef[match(material,names(coefs$unmethy.coef))]),ncol=length(material))
  methy.b[methy.b < 0] <- 0
  unmethy.b[unmethy.b < 0] <- 0
  methy.ba <- 2^methy.b
  unmethy.ba <- 2^unmethy.b

  # return methy.ba and unmethy.ba