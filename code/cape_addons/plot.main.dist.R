#This function plots main effects distributions for
#individual markers conditioned on all other markers
#the gen.coding argument is only used if marker2 is 
#not null

plot.main.dist <- function(data.obj, perm.data, marker.name, marker2 = NULL, allele = NULL, standardized = FALSE, plot.results = TRUE, show.means = FALSE, geno.coding = c("Additive", "Dominant", "Recessive")){
	
	library(RColorBrewer)
	cols <- brewer.pal(9, "Set1")
	
	#============================================================================
	# internal functions
	#============================================================================
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
	#============================================================================
	
	marker2.mat <- matrix(0);colnames(marker2.mat) <- "double.effect"
	if(!is.null(marker2)){
		marker2.effects <- lapply(colnames(data.obj$pheno), function(x) get.pheno.effects(data.obj, marker.name, marker2, x, scan.what = "norm", collapsed.net = FALSE, geno.coding = geno.coding))
		marker2.mat <- Reduce("rbind", marker2.effects)
		rownames(marker2.mat) <- colnames(data.obj$pheno)
		}	
	
	true.marker.name <- get.marker.name(marker.name)

	var.to.pheno.influence <- perm.data[[1]]
	pheno.names <- names(var.to.pheno.influence)

	all.inf <- vector(mode = "list", length = length(pheno.names))
	names(all.inf) <- pheno.names
	
	for(i in 1:length(var.to.pheno.influence)){	
		marker.locale1 <- which(var.to.pheno.influence[[i]][,1] == true.marker.name)
		marker.locale2 <- which(var.to.pheno.influence[[i]][,2] == true.marker.name)
		marker.effects <- c(as.numeric(var.to.pheno.influence[[i]][marker.locale1,3]), as.numeric(var.to.pheno.influence[[i]][marker.locale2,4]))
		marker.se <- c(as.numeric(var.to.pheno.influence[[i]][marker.locale1,5]), as.numeric(var.to.pheno.influence[[i]][marker.locale2,6]))
		
		
		if(standardized){
			all.inf[[i]] <- marker.effects/marker.se
			x.text <- "Standardized Effects"
			}else{
			all.inf[[i]] <- marker.effects
			x.text <- "Effects"
			}
		}


	all.dens <- lapply(all.inf, density)
	min.x <- min(c(unlist(lapply(all.dens, function(x) min(x$x))), as.numeric(marker2.mat[,"double.effect"])))
	max.x <- max(c(unlist(lapply(all.dens, function(x) min(x$x))), as.numeric(marker2.mat[,"double.effect"])))
	min.y <- min(unlist(lapply(all.dens, function(x) min(x$y))))
	max.y <- max(unlist(lapply(all.dens, function(x) max(x$y))))


	if(plot.results){
		plot.new()
		plot.window(xlim = c(min.x, max.x), ylim = c(min.y,max.y))
		for(i in 1:length(all.inf)){
			points(all.dens[[i]], type = "l", col = cols[i], lwd = 3)
			if(show.means){
				arrows(x0 = mean(all.inf[[i]]), y0 = max.y*0.05, y1 = 0, col = cols[i], lwd = 3, length = 0.1)
				}
			}
		if(!is.null(marker2)){
			par(xpd = TRUE)
			points(x = as.numeric(marker2.mat[,"double.effect"]), y = rep(max.y*1.05, length(pheno.names)), col = cols[1:length(pheno.names)], pch = "*", cex = 2)
			par(xpd = FALSE)
			mtext(paste("*interaction with", marker2), side = 3, line = 1.5)
			}
		legend("topright", lty = 1, col = cols[1:length(pheno.names)], legend = pheno.names, lwd = 3)
		axis(1);axis(2)	
		mtext(x.text, side = 1, line = 2.5)
		mtext("Density", side = 2, line = 2.5)
		mtext(marker.name, side = 3, line = 2.5)
		
		}
		
	invisible(all.inf)
		
		
}