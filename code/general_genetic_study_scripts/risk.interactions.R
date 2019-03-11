#This function calculates the interaction
#between exposures and an outcome
# snp1 <- "rs204990"; snp2 <- "rs10508168" #protective motif
# snp1 <- "rs3130933"; snp2 <- "rs10501351" #risk motif
# snp1 <- sample(colnames(pop$geno.for.pairscan), 1); snp2 <- sample(colnames(pop$geno.for.pairscan), 1) #random
#exposure.vector1 <- pop$geno.for.pairscan[,which(colnames(pop$geno.for.pairscan) == snp1)]
#exposure.vector1[which(exposure.vector1 == 0.5)] <- 1
#exposure.vector2 <- pop$geno.for.pairscan[,which(colnames(pop$geno.for.pairscan) == snp2)]
#exposure.vector2[which(exposure.vector2 == 0.5)] <- 1
#outcome.vector <- pop$raw.pheno[,2] #topo


risk.interactions <- function(exposure.vector1, exposure.vector2, outcome.vector, outcome.threshold = 0.5, plot.results = FALSE){
	
	no.exp1 <- which(exposure.vector1 == 0); exp1 <- which(exposure.vector1 == 1)
	no.exp2 <- which(exposure.vector2 == 0); exp2 <- which(exposure.vector2 >= 1)
	with.outcome <- which(outcome.vector >= outcome.threshold)
	
	#probability of outcome if exposed to both factors
	r11 <- length(Reduce("intersect", list(exp1, exp2, with.outcome)))/length(intersect(exp1, exp2))
	r01 <- length(Reduce("intersect", list(no.exp1, exp2, with.outcome)))/length(intersect(no.exp1, exp2))
	r10 <- length(Reduce("intersect", list(exp1, no.exp2, with.outcome)))/length(intersect(exp1, no.exp2))
	r00 <- length(Reduce("intersect", list(no.exp1, no.exp2, with.outcome)))/length(intersect(no.exp1, no.exp2))
	
	IC <- r11 - r10 - r01 + r00
	
	result <- c("0/0" = r00, "1/0" = r10, "0/1" = r01, "1/1" = r11, "IC" = IC)
	
	if(plot.results){
		barplot(result)
		}
	
	return(result)
	
}