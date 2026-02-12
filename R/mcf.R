
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

mcf_lookup <- function (mc) {
	case_when(
		mc %in% c("A IDH","A IDH, HG","O IDH") ~ "MTGF_IDH_GLM",
		mc %in% c("MB, G3","MB, G4") ~ "MTGF_MB_G3G4",
		mc %in% c("MB, SHH CHL AD","MB, SHH INF") ~ "MTGF_MB_SHH",
		mc %in% c("ENB, A","ENB, B") ~ "MTGF_ENB_nC",
		mc %in% c("GBM, MES","GBM, RTK I","GBM, RTK II","GBM, MYCN","GBM, MID","GBM, RTK III") ~ "MTGF_GB",
		mc %in% c("ATRT, TYR","ATRT, SHH","ATRT, MYC") ~ "MTGF_ATRT",
		mc %in% c("LGG, PA/GG ST", "LGG, PA PF","LGG, PA MID") ~ "MTGF_PA",
		mc %in% c("PLEX, PED A", "PLEX, AD","PLEX, PED B") ~ "MTGF_PLEX_T",
		TRUE ~ mc
	)
}




