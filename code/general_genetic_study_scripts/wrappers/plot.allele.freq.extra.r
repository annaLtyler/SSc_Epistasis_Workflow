#extra code from plot.allele.freq.R

		#This function plots frequencies at individual locus pairs		
		plot.pair.freq <- function(freq.mat, locus1, locus2){
			data.mat <- freq.mat[[1]]
			cols <- c("blue", "purple", "green")
			plot.new()
			plot.window(c(1,3), c(0,0.5))
			mtext("Frequency", 2, line = 3)
			mtext(locus1, 1, line = 3)
			axis(1, at = c(1:3), labels = genotypes)
			axis(2)
			for(i in 1:3){
				points(1:3, data.mat[i,], col = cols[i], pch = 16, type = "b", cex = 1.5, lwd = 2)
				}
			legend("topright", legend = genotypes, col = cols, lty = 1, lwd = 2)
			mtext(paste(locus1, "and", locus2), outer = FALSE, line = 0)
			}


		pdf("Locus_Pair_Counts_per_Genotype_Pair.pdf", width = 15, height = 5)
		plot.new()
		plot.window(ylim = c(0, max(pair.freq)), xlim = c(1,length(pair.freq[,1])))
		for(i in 1:3){
			points(pair.freq[,i], type = "l", col = i)
			}
		mtext("Number of Animals", side = 2, line = 2)
		mtext("Locus Pair", side = 1, line = 2)
		axis(1)
		axis(2)
		dev.off()

		
		pdf("Locus_Pair_Genotype_Frequencies.pdf", width = 10, height = 10)
		par(mfrow = c(3,3))
		for(i in 1:length(pair.freq)){
			plot.pair.freq(pair.freq[[i]], colnames(geno.table)[pairs.mat[i,1]], colnames(geno.table)[pairs.mat[i,2]])
			}
		dev.off()
