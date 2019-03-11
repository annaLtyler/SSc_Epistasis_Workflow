#This function takes in the data object, a marker name, and a phenotype
#and plots phenotype means and standard deviations for each factor in
#the genotype
#if a second marker is added, the function plots an interaction plot
#for the two markers.

plot.effects.2 <- function(data.obj, geno.obj = NULL, marker, marker2 = NULL, error.bars = TRUE, ymin = NULL, ymax = NULL, bin.geno = c(0, 0.5, 1), box.plot = FALSE){
	

	error.type = "se"
	
	geno <- get.geno(data.obj, geno.obj)
	
	marker.names <- c(marker, marker2)
	
	all.pheno.mat <- cbind(data.obj$pheno, data.obj$ET)
	all.pheno <- colnames(all.pheno.mat)

	if(is.null(marker2)){
		layout.mat <- matrix(1:length(all.pheno), nrow = 1)
		}else{
		layout.mat <- matrix(c(1:(length(all.pheno)*3)), nrow = 3, byrow = TRUE)	
		}
		
	dev.new(width = dim(layout.mat)[2]*3, height = dim(layout.mat)[1]*4)
	layout(layout.mat)
	
	for(m in 1:length(markers)){	
		
		marker.mat <- geno[,markers[m]]
			if(!is.null(bin.geno)){
			marker.mat <- bin.vector(marker.mat, bins = bin.geno)
			}
		
		for(ph in 1:length(all.pheno)){
		
			genotypes <- levels(as.factor(marker.mat))
			
	
			val.list <- lapply(genotypes, function(x) all.pheno.mat[which(marker.mat == x),ph])
			names(val.list) <- genotypes
			all.mean <- unlist(lapply(val.list, function(x) mean(x, na.rm = TRUE)))
			all.sd <- unlist(lapply(val.list, function(x) sqrt(var(x, na.rm = TRUE))))
		
			val.lims <- c((all.mean+all.sd), (all.mean-all.sd))

			if(error.bars){
				if(is.null(ymin)){ymin.cur <- floor(min(val.lims))}
				if(is.null(ymax)){ymax.cur <- ceiling(max(val.lims))}
				}else{
				if(is.null(ymin)){ymin.cur <- min(all.mean)}
				if(is.null(ymax)){ymax.cur <- max(all.mean)}
				}
				
			if(box.plot){				
				boxplot(val.list, notch = TRUE, main = paste(marker.names[m], all.pheno[ph], sep = "\n"))
				}else{
					stripchart(val.list, vertical = TRUE, pch = 16, method = "jitter", axes = FALSE, main = paste(marker.names[m], all.pheno[ph], sep = "\n"))
					axis(1, at = 1:length(all.mean), labels = genotypes)
					axis(2)
					#add lines for means
					segments(x0 = 1:length(all.mean)-0.1, y0 = all.mean, x1 = 1:length(all.mean)+0.1, y1 = all.mean, lwd = 3, col = "red")
					#add lines for se
					if(error.bars){
						segments(x0 = 1:length(all.sd)-0.1, y0 = all.mean+all.sd, x1 = 1:length(all.sd)+0.1, y1 = all.mean+all.sd, lwd = 1, col = "red", lty = 2)
						segments(x0 = 1:length(all.sd)-0.1, y0 = all.mean-all.sd, x1 = 1:length(all.sd)+0.1, y1 = all.mean-all.sd, lwd = 1, col = "red", lty = 2)
						}
				}
			
			}
		}
			if(!is.null(marker2)){
				for(ph in 1:length(all.pheno)){
					if(error.bars){
						marker1.geno <- geno[,markers[1]]; if(!is.null(bin.geno)){marker1.geno <- bin.vector(marker1.geno, bin.geno)}
						marker2.geno <- geno[,markers[2]]; if(!is.null(bin.geno)){marker2.geno <- bin.vector(marker2.geno, bin.geno)}
						errors <- get.interaction.error(marker1.geno, marker2.geno, all.pheno.mat[,ph], error.type = error.type)
						ylim <- c(min((errors$means - errors[[2]]), na.rm = TRUE), max((errors$means + errors[[2]]), na.rm = TRUE))
						interaction.plot(marker1.geno, marker2.geno, all.pheno.mat[,ph], xlab = marker, trace.label = marker2, ylab = all.pheno[ph], lwd = 3, cex.lab = 1.7, main = all.pheno[ph], ylim = ylim)
						}else{
						interaction.plot(marker1.geno, marker2.geno, all.pheno.mat[,ph], xlab = marker, trace.label = marker2, ylab = all.pheno[ph], lwd = 3, cex.lab = 1.7, main = all.pheno[ph])	
							}
					if(error.bars){
						for(i in 1:length(errors$means[,1])){
							segments(1:length(colnames(errors$means)), (errors$means[i,]+errors[[2]][i,]), 1:length(colnames(errors$means)), (errors$means[i,]-errors[[2]][i,]))
							}
						}

					}
				}
		
	
	
}