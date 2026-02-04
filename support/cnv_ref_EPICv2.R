
library(GEOquery)
library(ntile)
library(splitstackshape)
library(sesame)
library(conumee2)
library(here)

source(here("R","constants.R"))
source(here("R","preprocessing.R"))

num_bins <- 8
samples_per_bin <- 3

anno <- getSampleAnnotations(CNV_GSE_ID)

m <- anno[anno$`Sex:ch1` == "M",c("geo_accession", "description", "age:ch1", "Sex:ch1")] 
f <- anno[anno$`Sex:ch1` == "F",c("geo_accession", "description", "age:ch1", "Sex:ch1")] 

m_binned <- m %>% mutate(bin = ntile(`age:ch1`, n = num_bins))
f_binned <- f %>% mutate(bin = ntile(`age:ch1`, n = num_bins))

m_sample <- stratified(
  m_binned,          			# The binned dataframe
  group = "bin",      			# The column to stratify by
  size = samples_per_bin, 		# Number of samples per bin
  replace = FALSE     			# Sample without replacement
)
f_sample <- stratified(
  f_binned,          			# The binned dataframe
  group = "bin",      			# The column to stratify by
  size = samples_per_bin, 		# Number of samples per bin
  replace = FALSE     			# Sample without replacement
)

m_files <- paste(m_sample$geo_accession, m_sample$description, sep = "_")
f_files <- paste(f_sample$geo_accession, f_sample$description, sep = "_")

files <- c(m_files, f_files)

# copy files to ref directory
src <- here("data","GSE246337")
des <- here("data","CNV_REF_EPICv2")
for (file in files) {

	file_name <- paste0(file, "_Grn.idat")
	source_path <- file.path(src, file_name)	
	destination_path <- file.path(des, file_name)
	file.copy(from = source_path, to = destination_path)
	
	file_name <- paste0(file, "_Red.idat")
	source_path <- file.path(src, file_name)	
	destination_path <- file.path(des, file_name)
	file.copy(from = source_path, to = destination_path)
}

# read in files using sesame
sdfs_c <- openSesame(des, prep = "QCDPB", func = NULL)
cnv_ref <- CNV.load(do.call(cbind, lapply(sdfs_c, totalIntensities)))

# save ref data
saveRDS(cnv_ref, file=here("results", "CNV_REF_EPICv2.rds"))