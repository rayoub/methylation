
library(Biobase)

MCF <- list(
	MTGF_IDH_GLM = 
		c("A IDH","A IDH, HG","O IDH"),
	MTGF_MB_G3G4 = 
		c("MB, G3","MB, G4"),
	MTGF_MB_SHH = 
		c("MB, SHH CHL AD","MB, SHH INF"),
	MTGF_ENB_nC = 
		c("ENB, A","ENB, B"),
	MTGF_GB = 
		c("GBM, MES","GBM, RTK I","GBM, RTK II","GBM, MYCN","GBM, MID","GBM, RTK III"),
	MTGF_ATRT =
		c("ATRT, TYR","ATRT, SHH","ATRT, MYC"),
	MTGF_PA = 
		c("LGG, PA/GG ST", "LGG, PA PF","LGG, PA MID"),
	MTGF_PLEX_T = 
		c("PLEX, PED A", "PLEX, AD","PLEX, PED B")
)

MCF_LOOKUP <- reverseSplit(MCF)






