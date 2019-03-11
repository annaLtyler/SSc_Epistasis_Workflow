#This function 
# data.obj <- readRDS('~/Documents/Data/Scleroderma/Results/Dominant_lung_pheno_2ET/cross.RData')
# path <- "~/Documents/Data/Scleroderma/Results/Dominant_lung_pheno_2ET/SSc_snpnexus_11708"

write.snp.data <- function(data.obj, p.or.q = 0.05, path = "."){
	
	cur.dir <- getwd()
	var.inf <- writeVariantInfluences(data.obj, p.or.q, include.main.effects = FALSE, write.file = FALSE)
	u_snps <- unique(c(var.inf[,1], var.inf[,4]))
	just.snps <- sapply(strsplit(u_snps, "_"), function(x) x[1])
	
	setwd(path)
	for(i in 1:length(just.snps)){
		system(paste0("grep ", just.snps[i], " *.txt > compiled.", just.snps[i], ".txt"))
		}	
	setwd(cur.dir)	
	
	
	
}