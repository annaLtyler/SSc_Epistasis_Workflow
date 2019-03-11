#This function finds markers with main effects revealed by cape
#some markers don't have a significant main effect in the main
#scan, but do in cape.


main.effects.revealed <- function(data.obj, singlescan.obj, p.or.q = 0.05){
	
	
	cape.results <- data.obj$max.var.to.pheno.influence
	pheno.names <- names(cape.results)
	
	singlescan.results <- singlescan.obj$singlescan.results
	singlescan.names <- names(singlescan.results)
	
	if(!identical(sort(pheno.names), sort(singlescan.names))){
		stop("The phenotype names from the singlescan do not match the phenotype names in the final cape analysis.")
		}
		
	revealed.markers <- vector(mode = "list", length = length(singlescan.names))
	names(revealed.markers) <- singlescan.names
	
	for(i in 1:length(pheno.names)){
		cape.marker.locale <- which(cape.results[[i]][,7] <= p.or.q)
		cape.markers <- cape.results[[i]][cape.marker.locale,1]
		
		#get the p values for these markers in the singlescan.obj
		single.marker.locale <- match(cape.markers, rownames(singlescan.results[[i]]))
		single.cape.overlap <- singlescan.results[[i]][single.marker.locale,]
		
		non.sig <- which(single.cape.overlap[,"p.val"] > p.or.q)
		cape.sig.single.non.sig <- single.cape.overlap[non.sig,]
		
		revealed.markers[[i]] <- sort(as.numeric(rownames(cape.sig.single.non.sig)))
		}
	
	all.markers <- sort(unique(unlist(revealed.markers)))
	num.markers <- length(all.markers)
	results.table <- matrix(0, nrow = num.markers, ncol = length(singlescan.names))
	colnames(results.table) <- singlescan.names
	rownames(results.table) <- 	all.markers
	
	for(i in 1:length(revealed.markers)){
		if(length(revealed.markers[[i]]) > 0){
			marker.locale <- match(revealed.markers[[i]], all.markers)
			results.table[marker.locale,i] <- 1
			}
		}
	
	results <- list(revealed.markers, results.table)
	
	return(results)
	
	
	
	
}