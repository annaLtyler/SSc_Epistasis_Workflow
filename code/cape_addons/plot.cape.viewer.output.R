#This function takes in the data object, a marker name, and a phenotype
#and plots phenotype means and standard deviations for each factor in
#the genotype
#if a second marker is added, the function plots an interaction plot
#for the two markers.

plot.cape.viewer.output <- function(viewer.output, num.rows = 1){

 	
 	for(i in 1:length(viewer.output)){
		assign(names(viewer.output)[i], viewer.output[[i]])
		}

 	marker <- colnames(viewer.output[[3]])[1]
 	marker2 <- colnames(viewer.output[[3]])[2]

	geno.col = c("purple", "green", "blue", "red", "black", "orange", "gold")
	
	if(plot.type == "Interaction Plot"){plot.type <- "l"}
	if(plot.type == "Scatter Plot"){plot.type <- "p"}
	
	upper.plot.buffer = 0.5
	lower.plot.buffer = 0.2
	mean.bar.width = 0.15
	jitter.factor = 0.5

	type <- grep("l", plot.type)
	if(length(type) == 1){
		plot.type = "l"
		}else{
		error.bars = FALSE	
		}
	
	if(is.null(marker) && is.null(marker2)){
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		text(0.5, 0.5, "No Markers Selected")
		}else{

		error.type = "se"
	
	
		marker.names <- c(marker, marker2)
	
		if(is.null(marker) || is.null(marker2)){ #if one of the markers is null
			if(marker1.label == "" && marker2.label == ""){ #and both of the labels are blank
				if(marker1.label == ""){ #automatically assign a markername
					marker1.label = marker
					}
				if(marker2.label == ""){
					marker2.label = marker2
					}										
				}
			marker.label = c(marker1.label, marker2.label)
			}else{	#if both markers are assigned
			if(marker1.label == ""){ #automatically assign a marker name
				marker1.label = marker
				}
			if(marker2.label == ""){
				marker2.label = marker2
				}										
			}


		all.pheno.mat <- viewer.output$phenotypes
		all.pheno <- colnames(all.pheno.mat)
	
		num.col <- ceiling(length(all.pheno)/num.rows)
		plot.cells <- num.col*num.rows
		na.padding <- rep(0, plot.cells - length(all.pheno))
		layout.mat <- matrix(c(1:length(all.pheno), na.padding), nrow = num.rows)
	
		layout(layout.mat)
		par(mar = c(5,5,5,5))
		
		for(m in 1:length(marker.names)){
			
			marker.mat <- viewer.output$marker.genotypes[,m]
			geno.bins <- sort(unique(marker.mat[which(!is.na(marker.mat))]))
			if(length(geno.bins) > 3){
				geno.bins <- c(0, 0.5, 1)
				}
			geno.col <- geno.col[1:length(geno.bins)]
				if(!is.null(geno.bins)){
				marker.mat <- bin.vector(marker.mat, bins = geno.bins)
				}
			
			
			ind.cols <- rep(NA, dim(all.pheno.mat)[1])
			for(cl in 1:length(geno.bins)){
				ind.cols[which(marker.mat == geno.bins[cl])] <- geno.col[cl]
				}
			
			for(ph in 1:length(all.pheno)){
				
				if(length(covar) > 0){
					covar.mat <- covar
					covar.names <- colnames(covar.mat)
					covar.means <- apply(covar.mat, 2, function(x) mean(x, na.rm = TRUE))
					for(i in 1:length(covar.names)){
						covar.mat[,i] <- covar.mat[,i] - covar.means[i]
						}
					model <- lm(all.pheno.mat[,ph]~covar.mat)
					resids <- residuals(model)
							
					na.locale <- unique(which(is.na(cbind(all.pheno.mat[,ph], covar.mat)), arr.ind = TRUE)[,1])
					not.na <- c(1:dim(all.pheno.mat)[1])
					if(length(na.locale) > 0){
						not.na <- not.na[-na.locale]
						}
	
					all.pheno.mat[not.na,ph] <- resids
					}
					model <- lm(all.pheno.mat[,ph]~marker.mat)
			
				genotypes <- levels(as.factor(marker.mat))
		
				all.sd <- NULL
				all.mean <- NULL
			
				for(i in 1:length(genotypes)){
					all.mean <- c(all.mean, mean(all.pheno.mat[which(marker.mat == genotypes[i]),ph], na.rm = TRUE))
					if(error.type == "sd"){
						all.sd <- c(all.sd, sd(all.pheno.mat[which(marker.mat == genotypes[i]),ph], na.rm = TRUE))
						}else{
						all.sd <- c(all.sd, sqrt(var(all.pheno.mat[which(marker.mat == genotypes[i]),ph], na.rm = TRUE)))
						}
					}
				val.lims <- c((all.mean+all.sd), (all.mean-all.sd))
	
				ymin.cur <- ymin
				ymax.cur <- ymax
	
				if(error.bars > 0){
					if(is.null(ymin)){ymin.cur <- floor(min(val.lims))}
					if(is.null(ymax)){ymax.cur <- ceiling(max(val.lims))}
					}else{
					if(is.null(ymin)){ymin.cur <- min(all.mean)}
					if(is.null(ymax)){ymax.cur <- max(all.mean)}
					}
					
				if(length(marker.names) == 1){
					if(plot.type == "l"){
						plot(all.mean, type = "b", ylim = c(ymin.cur, ymax.cur), xlim = c(1,length(all.sd)), axes = FALSE, ylab = "", xlab = marker.label, main = paste(all.pheno[ph], "\np =", signif(anova(model)$"Pr(>F)", 2)[1]), lwd = 3, cex.lab = 2.5, cex.main = 2.5)
					axis(1, at = 1:length(all.mean), labels = genotypes, cex.axis = 2)
					axis(2, at = geno.bins, cex.axis = 2)
						}else{
							plot(jitter(marker.mat, factor = jitter.factor), all.pheno.mat[,ph], col = ind.cols, xlim = c(min(as.numeric(genotypes))-lower.plot.buffer, max(as.numeric(genotypes))+upper.plot.buffer), axes = FALSE, xlab = marker.label, ylab = all.pheno[ph])
							axis(1, at = as.numeric(genotypes), labels = genotypes)
							axis(2)
							segments(x0 = (as.numeric(genotypes) - mean.bar.width), y0 = all.mean, x1 = (as.numeric(genotypes) + mean.bar.width), col = "black", lwd = 3)
							}
	
					#add error bars
					if(error.bars){
						segments(x0 = 1:length(all.mean), y0 = (all.mean-all.sd), x1 = 1:length(all.mean), y1 = (all.mean+all.sd), lwd = 3)
						} 
					} #end case for if there is one marker
				} #end looping through phenotypes
			}#end looping through markers

			# if(!is.null(marker2)){
			if(length(marker.names) == 2){
				for(ph in 1:length(all.pheno)){

						marker1.geno <- viewer.output$marker.genotypes[,1]; if(!is.null(geno.bins)){marker1.geno <- bin.vector(marker1.geno, geno.bins)}
						marker2.geno <- viewer.output$marker.genotypes[,2]; if(!is.null(geno.bins)){marker2.geno <- bin.vector(marker2.geno, geno.bins)}
						errors <- get.interaction.error(marker1.geno, marker2.geno, all.pheno.mat[,ph], error.type = error.type)
					if(error.bars){
						if(is.null(ymin)){
							ylim <- c(min((errors$means - errors[[2]]), na.rm = TRUE), max((errors$means + errors[[2]]), na.rm = TRUE))
							}else{
								ylim <- c(ymin, ymax)
								}
							interaction.plot(marker1.geno, marker2.geno, all.pheno.mat[,ph], xlab = marker1.label, trace.label = marker2.label, ylab = "", lwd = 3, cex.lab = 2.5, cex.main = 2.5, main = all.pheno[ph], ylim = ylim, axes = FALSE, fixed = TRUE)
							plot.lim <- par("usr")
							min.area = plot.lim[1]*1.17; max.area = plot.lim[2]*0.85
							axis(1, at = c(seq(min.area, max.area, (max.area-min.area)/(length(geno.bins)-1))), labels = geno.bins, cex.axis = 2); axis(2, cex.axis = 2)							
						}else{
						if(is.null(ymin)){
							ylim <- c(min(errors$means, na.rm = TRUE), max(errors$means, na.rm = TRUE))
							# ylim <- c(min(all.pheno.mat[,ph], na.rm = TRUE), max(all.pheno.mat[,ph], na.rm = TRUE))
							}else{
								ylim <- c(ymin, ymax)
								}

						if(plot.type == "l"){
						interaction.plot(marker1.geno, marker2.geno, all.pheno.mat[,ph], xlab = marker1.label, trace.label = marker2.label, ylab = "", lwd = 3, cex.lab = 2.5, main = all.pheno[ph], axes = FALSE, ylim = ylim, fixed = TRUE, cex.main = 2.5)
						plot.lim <- par("usr")
						min.area = plot.lim[1]*1.17; max.area = plot.lim[2]*0.85
						axis(1, at = c(seq(min.area, max.area, (max.area-min.area)/(length(geno.bins)-1))), labels = geno.bins, cex.axis = 2); axis(2, cex.axis = 2)
						}else{
							ind.cols <- rep(NA, length(marker2.geno))
							for(cl in 1:length(geno.bins)){
								ind.cols[which(marker2.geno == geno.bins[cl])] <- geno.col[cl]
								}
							plot(jitter(marker1.geno, factor = jitter.factor), all.pheno.mat[,ph], col = ind.cols, xlim = c(min(as.numeric(genotypes))-lower.plot.buffer, max(as.numeric(genotypes))+upper.plot.buffer), axes = FALSE, xlab = marker1.label, ylab = all.pheno[ph], pch = 1)	
							axis(1, at = as.numeric(genotypes), labels = genotypes)
							axis(2)
							for(g in 1:length(genotypes)){
								segments(x0 = (as.numeric(genotypes) - mean.bar.width), y0 = errors$means[g,], x1 = (as.numeric(genotypes) + mean.bar.width), col = geno.col[g], lwd = 3)
								}
							legend("topright", legend = geno.bins, col = geno.col, pch = 16, title = marker2.label)

							}

							}
					if(error.bars){
						for(i in 1:length(errors$means[,1])){
							segments(1:length(colnames(errors$means)), (errors$means[i,]+errors[[2]][i,]), 1:length(colnames(errors$means)), (errors$means[i,]-errors[[2]][i,]))
							}
						}

					}
				}#end case for if we have selected two markers
		}#end of case for if no markers are selected			

	
}