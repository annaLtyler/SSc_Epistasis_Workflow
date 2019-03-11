
plot.points <- function(marker.geno, marker.name, phenotype.vals, phenotype.name, ymin = NULL, ymax = NULL){
	
	geno.col = c("purple", "green", "blue", "red", "black", "orange", "gold")
	mean.bar.width = 0.15
	jitter.factor = 0.1
	upper.plot.buffer = 0.5
	lower.plot.buffer = 0.2


	if(is.null(ymin)){
		ymin <- min(phenotype.vals, na.rm = TRUE)	
		}
	if(is.null(ymax)){
		ymax <- max(phenotype.vals, na.rm = TRUE)
		}

	if(dim(marker.geno)[2] == 1){

		#get the genotype values
		genotypes <- sort(unique(marker.geno[which(!is.na(marker.geno))]))
			
		stripchart(phenotype.vals~as.factor(marker.geno), vertical = TRUE, method = "jitter", jitter = jitter.factor, pch = 1, col = geno.col[1:length(genotypes)], xlab = marker.name, ylim = c(ymin, ymax))
		xlim = c(min(as.numeric(genotypes))-lower.plot.buffer, max(as.numeric(genotypes))+upper.plot.buffer)
		
		mean.test <- boxplot(phenotype.vals~as.factor(marker.geno), plot = FALSE)
		
		segments(x0 = (1:length(genotypes) - mean.bar.width), y0 = mean.test[[1]][3,], x1 = (1:length(genotypes) + mean.bar.width), col = "black", lwd = 3)
		}
	
	
	if(dim(marker.geno)[2] == 2){

		#get the genotype values for marker 2
		genotypes <- sort(unique(marker.geno[which(!is.na(marker.geno[,2])),2]))


		ind.cols <- rep(NA, length(phenotype.vals))
		for(g in 1:length(genotypes)){
			ind.cols[which(marker.geno[,2] == genotypes[g])] <- geno.col[g]
			}

		plot(jitter(marker.geno[,1], factor = jitter.factor*5), phenotype.vals, col = ind.cols, xlim = c(min(as.numeric(genotypes))-lower.plot.buffer, max(as.numeric(genotypes))+upper.plot.buffer), axes = FALSE, xlab = marker.name[1], ylab = phenotype.name, pch = 1, ylim = c(ymin, ymax), main = phenotype.name)
		axis(1, at = as.numeric(genotypes), labels = genotypes)
		axis(2)
		for(g in 1:length(genotypes)){
			errors <- get.interaction.error(marker.geno[,1], marker.geno[,2], phenotype.vals, error.type = "se")
			segments(x0 = (as.numeric(genotypes) - mean.bar.width), y0 = errors$means[g,], x1 = (as.numeric(genotypes) + mean.bar.width), col = geno.col[g], lwd = 3)
			}
		legend("topright", legend = genotypes, col = geno.col, pch = 1, title = marker.name[2])

		
		} #end case for if there are two markers
		
	} #end function