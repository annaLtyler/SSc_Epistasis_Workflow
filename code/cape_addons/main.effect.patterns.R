#This function finds specific main effects
#in a cape object. 
#you can look for combinations of 
#positive "+"
#negative "-" 
#no significant main effects "0"
#or irrelevant "*"

main.effect.patterns <- function(data.obj, effect.pattern = "+", collapsed.net = TRUE){
	
	if(collapsed.net){
		full.net <- data.obj$collapsed.net
		}else{
		full.net <- data.obj$full.net	
		}
	
	just.main <- full.net[, (nrow(full.net)+1):ncol(full.net)]
	n.pheno <- ncol(just.main)
	if(length(effect.pattern) == 1){
		effect.pattern <- rep(effect.pattern, n.pheno)
		}
		
	if(length(effect.pattern) != n.pheno){stop("You must specify an effect pattern for each phenotype")}
	
	match.effect <- function(phenoV, effect){
		if(effect == "+"){pheno.ind <- which(phenoV > 0)}
		if(effect == "-"){pheno.ind <- which(phenoV < 0)}
		if(effect == "0"){pheno.ind <- which(phenoV == 0)}
		if(effect == "*"){pheno.ind <- 1:length(phenoV)}
		return(pheno.ind)
		}
	
	all.pheno.ind <- vector(mode = "list", length = n.pheno)
	for(i in 1:n.pheno){
		all.pheno.ind[[i]] <- match.effect(just.main[,i], effect.pattern[i])
		}
	
	match.pattern.ind <- Reduce('intersect', all.pheno.ind)

	if(length(match.pattern.ind) == 1){
		cat("One locus found that match pattern.\n")
		}else{
		cat(length(match.pattern.ind), "loci found.\n")
		}
	result <- just.main[match.pattern.ind,,drop=FALSE]
	return(result)
	
	
}