#This function tests the effects of a marker
#conditioned on another marker. For example,
#if you want to see whether marker A has
#a different effecct in males and females
#this function compares the models of the
#marker in the two populations.

test.marker.environment <- function(data.obj, test.marker, cond.marker, error.bars = FALSE){
	
	all.phenotypes <- cbind(data.obj$pheno, data.obj$ET)
	pheno.names <- colnames(all.phenotypes)

		
	test.locale <- get.col.num(data.obj$geno, test.marker, warn = TRUE)
	cond.locale <- get.col.num(data.obj$geno, cond.marker, warn = TRUE)


	plot.effects(data.obj, marker = test.marker, marker2 = cond.marker, error.bars = error.bars)

	#take out any missing values so we can compare the models
	non.missing.vals <- sort(unique(intersect(which(!is.na(data.obj$geno[,test.locale])), which(!is.na(data.obj$geno[,cond.locale])))))

	all.p <- matrix(NA, ncol = 1, nrow = length(pheno.names))
	rownames(all.p) <- pheno.names
	colnames(all.p) <- "p.value"
	for(ph in 1:length(pheno.names)){
		bare.model <- lm(all.phenotypes[non.missing.vals,ph]~data.obj$geno[non.missing.vals,test.locale]+data.obj$geno[non.missing.vals,cond.locale])
		# summary(bare.model)
		full.model <- lm(all.phenotypes[non.missing.vals,ph]~data.obj$geno[non.missing.vals,cond.locale]*data.obj$geno[non.missing.vals,test.locale])
		# summary(full.model)
		
		comparison <- anova(full.model, bare.model)
		all.p[ph,1] <- signif(comparison$"Pr(>F)"[2], 3)
		}

	return(all.p)
}