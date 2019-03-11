#This function calculates risk ratios for pairs of SNPs
#as well as a number of interaction terms as defined in
#Kalilani and Atashili 2006
#the outcome threshold is for when we have corrected the
#phenotype for a covariate and the outcome is no longer
#0 or 1 

get.snp.risk.ratio <- function(data.obj, geno.obj = NULL, snps, phenotype, allele.which = 2, covar = NULL, trait.type = c("ET", "Normalized", "Raw"), geno.coding = c("Additive", "Dominant", "Recessive"), outcome.threshold = 0.5, plot.results = TRUE, las = 2){
	
	if(is.null(geno.obj)){
		geno <- data.obj$geno.for.pairscan
		snp.locale <- match(snps, colnames(geno))
		geno.table <- geno[,snp.locale]
		}else{
		geno <- get.geno(data.obj, geno.obj)
		geno.table <- geno[,allele.which,snps,drop=FALSE]
		}
		
	cor.pheno <- get.pheno(data.obj, trait.type[1], covar)	
	pheno.locale <- which(colnames(data.obj$pheno) == phenotype)
	pheno <- cor.pheno[,pheno.locale, drop=FALSE]
	
	if(geno.coding == "Dominant"){
		geno.table[which(geno.table >= 0.5)] <- 1
		}
	if(geno.coding == "Recessive"){
		geno.table[which(geno.table <= 0.5)] <- 0
		}
					
	
	#calculate risk ratio for each SNP
	snp.rr.by.pheno <- matrix(NA, nrow = ncol(geno.table), ncol = ncol(pheno))
	colnames(snp.rr.by.pheno) <- colnames(pheno)
	rownames(snp.rr.by.pheno) <- dimnames(geno.table)[[3]]
	for(ph in 1:length(phenotype)){
	snp.rr.by.pheno[,ph] <- apply(geno.table, 2, function(x) risk.ratio(exposure.vector = x, outcome.vector = pheno[,ph], outcome.threshold = outcome.threshold))
	}
	
	if(plot.results){
	for(ph in 1:ncol(geno.table)){
		barplot(snp.rr.by.pheno[,ph], names = snps, ylab = "Risk Relative to Reference Allele", main = phenotype[ph], ylim = c(0, max(c(snp.rr.by.pheno[,ph], 1), na.rm = TRUE)), las = las)
		abline(h = 1)
		}
	}

	return(snp.rr)	
	}