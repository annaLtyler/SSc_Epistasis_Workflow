plot.permutations <-
function(data.obj, min.per.genotype = 6, num.pairs = 10, num.perm = 10, percentile = 95){

	verbose = TRUE
	geno <- data.obj$geno.for.pairscan
	covar <- data.obj$covar.for.pairscan
	
	pairs <- get.pairs.for.pairscan(geno, min.per.genotype = min.per.genotype, verbose = verbose)

	rnd.pairs <- pairs[sample(1:dim(pairs)[1], num.pairs),]
	
	
	ph <- 1
	#the num.perm we put in as an argument is the number of permutations for the total
	#calculate the number of total number of permutations to give each pair the number
	#permutations specified in the argument for the outer function
	total.perm <- num.pairs*num.perm
	results <- one.pairscan(data.obj$ET[,ph], geno, covar[,ph], rnd.pairs, n.perm = total.perm, verbose = verbose)
	total.stats <- dim(results$pairscan.perm[[1]])[2]

	# overall.min <- min(as.numeric(results$pairscan.perm[[1]][,3:total.stats])/as.numeric(results$pairscan.perm[[2]][,3:total.stats]), na.rm = TRUE)
	overall.min = 0
	overall.max <- max(as.numeric(results$pairscan.perm[[1]][,3:total.stats])/as.numeric(results$pairscan.perm[[2]][,3:total.stats]), na.rm = TRUE)
	
	get.percentile <- function(distX, percentile){
		ranked <- rank(distX)
		locale <- which(ranked == round(length(distX)*(percentile/100)))
		return(distX[locale])
		}

	
	#plot the histograms of the permutations for each pair
	plot.hist <- function(pair){
		pair.locale <- intersect(which(results$pairscan.perm[[1]][,1] == pair[1]), which(results$pairscan.perm[[1]][,2] == pair[2]))

		marker1.null <- abs(as.numeric(results$pairscan.perm[[1]][pair.locale,(total.stats-2)]))/(as.numeric(results$pairscan.perm[[2]][pair.locale,(total.stats-2)]))
		marker1.perc <- get.percentile(marker1.null, percentile)
		hist(marker1.null, main = "Marker1", xlab = "Marker1", xlim = c(overall.min, overall.max))
		abline(v = marker1.perc, col = "red")

		marker2.null <- abs(as.numeric(results$pairscan.perm[[1]][pair.locale,(total.stats-1)]))/(as.numeric(results$pairscan.perm[[2]][pair.locale,(total.stats-1)]))
		marker2.perc <- get.percentile(marker2.null, percentile)
		hist(marker2.null, main = "Marker2", xlab = "Marker2", xlim = c(overall.min, overall.max))
		abline(v = marker2.perc, col = "red")
		
		interaction.null <- abs(as.numeric(results$pairscan.perm[[1]][pair.locale,(total.stats)]))/(as.numeric(results$pairscan.perm[[2]][pair.locale,(total.stats)]))
		interaction.perc <- get.percentile(interaction.null, percentile)
		hist(interaction.null, main = "Interaction", xlab = "Interaction", xlim = c(overall.min, overall.max))
		abline(v = interaction.perc, col = "red")
		
		results <- c(marker1.perc, marker2.perc, interaction.perc)
		return(results)
		}
		
	
	pdf(paste("Null.Distributions.for.", num.perm, ".Permutations.pdf", sep = ""), width = 9, height = 10)
	layout(matrix(c(1:12), ncol = 3, byrow = TRUE))
	all.percentiles <- apply(rnd.pairs, 1, plot.hist)
	dev.off()

	#now plot the distribution of the percentiles found
	pdf(paste("Distribution.of.", percentile, ".percentile.for.", num.perm, ".Permutations.pdf", sep = ""), width = 9, height = 3)
	layout(matrix(c(1:3), ncol = 3))
	for(i in 1:dim(all.percentiles)[1]){
		hist(all.percentiles[i,], xlim = c(overall.min, overall.max), main = rownames(all.percentiles)[i], xlab = paste(percentile, "percentile distribution"))
	}
	dev.off()
}
