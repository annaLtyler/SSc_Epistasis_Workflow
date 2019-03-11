#This function calculates the risk ratio
#for an exposure and outcome
#This is the risk of the outcome in the exposed
#group relative to the risk of outcome in the
#unexposed group. If RR > 1, the exposure increase
#risk. If RR < 1, the exposure decreases risk
# snp1 <- "rs204990"; snp2 <- "rs10508168" #protective motif
# snp1 <- "rs3130933"; snp2 <- "rs10501351" #risk motif
# snp1 <- sample(colnames(pop$geno.for.pairscan), 1); snp2 <- sample(colnames(pop$geno.for.pairscan), 1) #random
#exposure.vector <- pop$geno.for.pairscan[,which(colnames(pop$geno.for.pairscan) == snp1)]
#exposure.vector[which(exposure.vector == 0.5)] <- 1
#exposure.vector <- pop$geno.for.pairscan[,which(colnames(pop$geno.for.pairscan) == snp2)]
#exposure.vector[which(exposure.vector == 0.5)] <- 1
#outcome.vector <- pop$raw.pheno[,2] #topo
#those with an outcome greater than the threshold have the outcome
#those with an outcome less than the threshold do not have the outcome

risk.ratio <- function(exposure.vector, outcome.vector, outcome.threshold = 0.5){
	
	with.outcome <- which(outcome.vector >= outcome.threshold)
	without.outcome <- which(outcome.vector < outcome.threshold)
	exposed <- which(exposure.vector == max(exposure.vector, na.rm = TRUE))
	unexposed <- which(exposure.vector == min(exposure.vector, na.rm = TRUE))
		
	exposed.with.outcome <- intersect(exposed, with.outcome)
	exposed.without.outcome <- intersect(exposed, without.outcome)
	unexposed.with.outcome <- intersect(unexposed, with.outcome)
	unexposed.without.outcome <- intersect(unexposed, without.outcome)
		
	outcome.with.exp.ci <- length(exposed.with.outcome)/length(exposed)
	outcome.without.exp.ci <- length(unexposed.with.outcome)/length(unexposed)

	rr <- outcome.with.exp.ci/outcome.without.exp.ci
	return(rr)
}