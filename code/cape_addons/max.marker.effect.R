#This function finds the maximum marker effect when conditioned
#on a secondary marker. It then plots the interaction between 
#the markers
#code from the shiny capeDO interaction viewer needs to be loaded
#to run this function


max.marker.effect <- function(data.obj, perm.data, marker.name, allele = NULL, standardized = TRUE, covar = NULL){
	
	
	library(RColorBrewer)
	cols <- brewer.pal(9, "Set1")
	
	get.marker.name <- function(marker.name){
		if(length(grep("Chr", marker.name)) > 0){
			marker.locale <- which(names(data.obj$linkage.blocks.full) == marker.name)
			marker.label <- data.obj$linkage.blocks.full[[marker.locale]]
			return(marker.label)
			}else{
			if(is.null(allele)){
				full.name <- marker.name
				}else{
				full.name <- paste0(marker.name, "_", allele)
				}
			return(marker.name)
			}
		}
		
	true.marker.name <- get.marker.name(marker.name)
	
	pheno <- get.pheno(data.obj, "normalized.trait", covar)
	geno <- data.obj$geno.for.pairscan

	var.to.pheno.influence <- perm.data[[1]]
	pheno.names <- names(var.to.pheno.influence)

	all.inf <- vector(mode = "list", length = length(pheno.names))
	names(all.inf) <- pheno.names	
	layout.mat <- get.layout.mat(length(pheno.names))
		
	layout(layout.mat)
	for(i in 1:length(var.to.pheno.influence)){	
		marker.locale1 <- which(var.to.pheno.influence[[i]][,1] == true.marker.name)
		marker.locale2 <- which(var.to.pheno.influence[[i]][,2] == true.marker.name)

		marker.effects <- c(as.numeric(var.to.pheno.influence[[i]][marker.locale1,3]), as.numeric(var.to.pheno.influence[[i]][marker.locale2,4]))
		marker.se <- c(as.numeric(var.to.pheno.influence[[i]][marker.locale1,5]), as.numeric(var.to.pheno.influence[[i]][marker.locale2,6]))
	
		if(standardized){
			max.locale <- which.max(abs(marker.effects/marker.se))
			}else{
			max.locale <- which.max(abs(marker.effects))	
			}
		
		marker.pair <- var.to.pheno.influence[[i]][c(marker.locale1, marker.locale2)[max.locale],,drop=FALSE]
		
		marker2 <- setdiff(marker.pair[1,1:2], true.marker.name)
		
		marker1.geno <- geno[,which(colnames(geno) == true.marker.name)]
		marker2.geno <- geno[,which(colnames(geno) == marker2)]
		
		cor.pheno <- rep(NA, length(marker2.geno))
		not.na.locale <- which(!is.na(marker2.geno))
		cor.pheno[not.na.locale] <- resid(lm(pheno[,i]~marker2.geno))
		
		genotypes <- sort(unique(marker1.geno))
		pheno.means <- unlist(lapply(genotypes, function(x) mean(cor.pheno[which(marker1.geno == x)])))
		pheno.se <- unlist(lapply(genotypes, function(x) sd(cor.pheno[which(marker1.geno == x)])/sqrt(length(which(marker1.geno == x)))))
		
		max.y <- max((pheno.means + pheno.se), na.rm = TRUE)
		min.y <- min((pheno.means - pheno.se), na.rm = TRUE)
				
		plot(genotypes, pheno.means, type = "b", main = paste("Effect of", marker.name, "on", pheno.names[i], "\nconditioned on", marker2), axes = FALSE, ylab = pheno.names[i], xlab = paste(marker.name, "Genotype"), ylim = c(min.y, max.y))
		segments(x0 = genotypes, y0 = (pheno.means - pheno.se), y1 = (pheno.means + pheno.se))
		axis(1, at = genotypes)
		axis(2)
		
		}


	
	
	
	
}
