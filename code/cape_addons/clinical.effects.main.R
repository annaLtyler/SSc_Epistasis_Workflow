#This function is related to clinical.benefits, but focuses
#on markers that have significant main effects. Main effects
#in cape are all conditioned on a second marker.
#This function looks at whether the marker with the main effect
#has consistent clinical effects across the traits when conditioned
#on the conditioning marker.
#It also looks at whether the interaction between the main effect
#marker and its conditioning marker are consistently clinically
#relevant.
# Clinidal relevance is based on whether high/low values of each
#phenotype is good/bad, the function decides which main effects are 
#beneficial for a patient, and which are harmful
#clinical.benefit is a vector that says whether each phenotype is
#beneficial if it is high (h), low (l), or doesn't matter (*)

clinical.benefits.main <- function(data.obj, clinical.benefit, p.or.q = 0.05, covar = NULL, scan.what = c("normalized.traits", "raw.traits"), geno.coding = c("Additive", "Dominant", "Recessive")){
	
	geno.coding <- geno.coding[1]
	pheno <- get.pheno(data.obj, scan.what)
	geno <- data.obj$geno.for.pairscan
	covar.info <- get.covar(cross)
	
	if(geno.coding == "Dominant"){
		geno[which(geno >= 0.5)] <- 1
		}
	if(geno.coding == "Recessive"){
		geno[which(geno <= 0.5)] <- 0
		}

	var.inf <- writeVariantInfluences(data.obj, p.or.q, include.main.effects = TRUE, write.file = FALSE)	
	
	#restrict to only main effects. clinical.benefits.R looks at interactions
	var.inf <- var.inf[which(!is.na(var.inf[,"Conditioning.Marker"])),]
	
	if(nrow(var.inf) == 0){
		stop("There are no main effects.")
		}
		
	
	
	
	#=========================================================
	# internal functions
	#=========================================================
		
	get.genotype <- function(marker.name){
		marker.locale <- which(colnames(geno) == marker.name)
		if(length(marker.locale) == 1){
			return(geno[,marker.locale])
			}else{
			marker.locale <- which(data.obj$p.covar == marker.name)	
			if(length(marker.locale) == 1){
				return(data.obj$p.covar.table[,marker.locale])
				}
			}
		}	
		
		
	get.pheno.vals <- function(marker1, marker2, phenotype){		
		marker1.geno <- get.genotype(marker1)
		marker2.geno <- get.genotype(marker2)
	
		max1 <- which(marker1.geno == max(marker1.geno))
		max2 <- which(marker2.geno == max(marker2.geno))
		min1 <- which(marker1.geno == min(marker1.geno))
		min2 <- which(marker2.geno == min(marker2.geno))

		baseline <- mean(phenotype[intersect(min1, min2)])
		just1.pheno <- mean(phenotype[intersect(max1, min2)]) - baseline
		just2.pheno <- mean(phenotype[intersect(min1, max2)]) - baseline
		both.pheno <- mean(phenotype[intersect(max1, max2)]) - baseline
		add.exp <- just1.pheno + just2.pheno
		
		if(is.finite(add.exp)){
			if(add.exp < 0){
				int.expect <- "expect.negative"
				}else{
				int.expect <- "expect.positive"
				}
			
			if(abs(both.pheno) < abs(add.exp)){
				int.effect <- "alleviating"
				}else{
				int.effect <- "aggravating"
				}
			}else{
				int.expect <- "none"
				int.effect <- "none"
				}
		
		pheno.result <- c(just1.pheno, just2.pheno, add.exp, both.pheno, int.expect, int.effect)
		return(pheno.result)
		}


	#find which effects on phenotypes are beneficial
	#and which are harmful
	evaluate.int <- function(result.mat){
		clin.effects <- rep(NA, nrow(result.mat))
		
		high.locale <- which(clinical.benefit == "h")
		low.locale <- which(clinical.benefit == "l")

		diff.expect <- as.numeric(result.mat[,"actual"]) - as.numeric(result.mat[,"expected.additive"])
		higher.than.expected <- which(diff.expect > 0)
		lower.than.expected <- which(diff.expect < 0)		
				
		#==============
		#beneficial
		#==============
		#for traits that we want high, interactions are
		#beneficial if the actual value is higher than 
		#the expected value
		ben.high <- intersect(high.locale, higher.than.expected)
		if(length(ben.high) > 0){
			clin.effects[ben.high] <- "beneficial"
			}
	
		#for traits that we want low, interactions
		#are beneficial if the actual value is lower
		#than the expected value
		ben.low <- intersect(low.locale, lower.than.expected)
		if(length(ben.low) > 0){
			clin.effects[ben.low] <- "beneficial"
			}
		
		#==============
		#harmful
		#==============
		#for traits that we want high, interactions are
		#harmful if the actual value is lower than 
		#the expected value
		harm.high <- intersect(high.locale, lower.than.expected)
		if(length(harm.high) > 0){
			clin.effects[harm.high] <- "harmful"
			}
	
		#for traits that we want low, interactions
		#are harmful if the actual value is higher than
		#expected
		harm.low <- intersect(low.locale, higher.than.expected)
		if(length(harm.low) > 0){
			clin.effects[harm.low] <- "harmful"
			}
		
		return(clin.effects)
		
		}
	
	evaluate.main <- function(pheno.effects){
		main.sign <- sign(as.numeric(pheno.effects[,1]))
		main.clin <- rep(NA, length(clinical.benefit))
		clin.sign <- rep(NA, length(clinical.benefit))
		clin.sign[which(clinical.benefit == "l")] <- -1
		clin.sign[which(clinical.benefit == "h")] <- 1
		clin.sign[which(clinical.benefit == "*")] <- 0
		
		for(i in 1:length(clinical.benefit)){
			if(main.sign[i] == clin.sign[i]){main.clin[i] <- "beneficial"}
			if(main.sign[i] == clin.sign[i]*-1){main.clin[i] <- "harmful"}
			if(clinical.benefit[i] == "*"){main.clin[i] <- NA}
			}
		return(main.clin)
	}
	
	
	assign.final <- function(effects){
		no.na <- effects[which(!is.na(effects))]
		u_effects <- unique(no.na)
		if(length(u_effects) == 1){
			return(u_effects)
			}else{
			return("mixed")
			}
	}
	
	scan.effects <- function(full.effects.mat){
		all.main <- full.effects.mat[,"main.effect"]
		all.int <- full.effects.mat[,"interaction.effect"]
		final.main <- assign.final(all.main)
		final.int <- assign.final(all.int)
		results <- c("combined.main" = final.main, "combined.int" = final.int)
		return(results)		
		}
	#=========================================================

	all.effects <- vector(mode = "list", length = nrow(var.inf))
	# snp.names <- apply(var.inf, 1, function(x) paste(x[1], x[7], sep = "_"))
	names(all.effects) <- var.inf[,1]

	
	for(i in 1:nrow(var.inf)){
		marker <- var.inf[i,1] #get the main effect marker
		c.marker <- var.inf[i,7] #and the conditioning marker

		# plot.effects(data.obj, marker, c.marker, plot.type = "b", error.bars = "se", covar = covar)

		#find phenotype values for individuals with different genotype combinations
		#and 
		all.pheno.effects <- matrix(NA, nrow = ncol(pheno), ncol = 6)
		for(ph in 1:ncol(pheno)){
			all.pheno.effects[ph,] <- get.pheno.vals(marker, c.marker, phenotype = pheno[,ph])	
			}
		rownames(all.pheno.effects) <- colnames(pheno)
		colnames(all.pheno.effects) <- c("marker", "conditioning.marker", "expected.additive", "actual", "expected.effect", "interaction.effect")

		int.effect <- evaluate.int(all.pheno.effects)
		main.effect <- evaluate.main(all.pheno.effects)
		
		all.pheno.effects <- cbind(all.pheno.effects, main.effect, int.effect)
		colnames(all.pheno.effects)[c(7,8)] <- c("main.effect", "interaction.effect")
		
		all.effects[[i]] <- all.pheno.effects
	
		}	

	clin.int <- Reduce("rbind", lapply(all.int.effects, scan.effects))
	final.clin <- cbind(var.inf[,1], var.inf[,4], clin.int)
	colnames(final.clin) <- c("Marker", "Conditioning.Marker", "Clinical_Main_Effect", "Clinical_Interaction_Effect")
	rownames(final.clin) <- NULL

	final.results <- list("clinical_summary" = final.clin, "details" = all.effects)
	return(final.results)
	}