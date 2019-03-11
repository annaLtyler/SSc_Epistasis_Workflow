#this function recodes the SSc snp matrix
#so that homozygous for the minor allele is 1 
#and homozygous for the major allele is 0

source('~/Documents/git_repositories/useful_r_code/usefulScripts/report.progress.R', chdir = FALSE)

geno.dir <- "~/Documents/Grants/R21/Scleroderma/GWAS_documents/geno/"
setwd(geno.dir)
geno.mat <-readRDS("GenotypeMatrix.RData") 

recode.snp <- function(genotypes){
	zero.allele.freq <- length(which(genotypes == 0))
	het.allele.freq <- length(which(genotypes == 0.5))
	one.allele.freq <- length(which(genotypes == 1))
	
	new.genotypes <- genotypes
	if(one.allele.freq > zero.allele.freq){
		new.genotypes[which(genotypes == 0)] <- 1
		new.genotypes[which(genotypes == 1)] <- 0
		}
	return(new.genotypes)
	}

new.geno.mat <- matrix(NA, nrow = nrow(geno.mat), ncol = ncol(geno.mat))

for(i in 209024:dim(geno.mat)[2]){
	report.progress(i, dim(geno.mat)[2])
	new.geno.mat[,i] <- recode.snp(geno.mat[,i])
	}
	
saveRDS(new.geno.mat, "GenotypeMatrixRecoded.RData")