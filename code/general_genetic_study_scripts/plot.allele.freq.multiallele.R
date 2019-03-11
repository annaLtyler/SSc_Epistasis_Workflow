
plot.allele.freq.multiallele <- function(geno.table, label = NULL, pairs = FALSE, percent.in.rare.geno = 5){


	loci <- colnames(geno.table)
	alleles <- dimnames(geno.table)[[3]]
	genotypes <- pairs.matrix(alleles, self.pairs = TRUE)
	cols = c("yellow", "black", "salmon", "blue", "lightblue", "green", "red", "purple")
	
	
	#given a 0, 0.5, 1 set of probabilities
	#for allele A at a given locus, this 
	#function returns the frequency of each
	#allele over all individuals
	allele.freq <- function(prob.mat){
		allele.totals <- colSums(prob.mat)
		allele.frequency <- allele.totals/(sum(allele.totals))
		return(allele.frequency)
		}
	
	#This function counts the frequency of all genotypes
	#at each locus
	diploid.freq <- function(prob.mat){
		diploid.matrix <- matrix(0, ncol = length(alleles), nrow = length(alleles))
		dips <- apply(prob.mat, 1, function(x) which(x > 0.4))
		for(i in 1:length(dips)){
			diploid.call <- dips[[i]]
			diploid.matrix[diploid.call[1], diploid.call[length(diploid.call)]] <- (diploid.matrix[diploid.call[1], diploid.call[length(diploid.call)]] + 1)	
			}
		dip.freqs <- diploid.matrix/(sum(diploid.matrix))
		colnames(dip.freqs) <- row.names(dip.freqs) <- alleles
		return(dip.freqs[upper.tri(dip.freqs, diag = TRUE)])
		}
	
	

	
	allele.frequencies <- apply(geno.table, 2, allele.freq)
	geno.freq <- apply(geno.table, 2, diploid.freq)
	rownames(geno.freq) <- apply(genotypes, 1, function(x) paste(x, sep = "", collapse = "_"))
	
	#reverse each column, so the A/J frequencies are at the top
	#The colors also need to be reversed
	allele.frequencies <- apply(allele.frequencies, 2, function(x) x[1:length(x)] <- x[length(x):1])
	
	pdf(paste("Allele.Frequencies.Color", label, ".pdf", sep = ""), width = 20, height = 5)
	# quartz(width = 10, height = 5)
	par(mar = c(4, 2, 4, 0))
	bar.mids <- barplot(allele.frequencies, border = NA, space = 0, axes = FALSE, ylim = c(0, 1), names.arg = rep("", length(loci)), main = "Allele Frequencies", xpd = TRUE, col = rev(cols))
	legend(-(length(loci)*0.04),1, legend = alleles, fill = cols)
	dev.off()

	#reverse the columns back, so we plot the individual
	#alleles in the right order
	allele.frequencies <- apply(allele.frequencies, 2, function(x) x[1:length(x)] <- x[length(x):1])

	layout.mat <- matrix(1:length(alleles), nrow = 8)
	pdf(paste("Allele.Frequencies.Individual", label, ".pdf", sep = ""), width = 10, height = 12)
	layout(layout.mat)
	par(mar = c(2,3,2,0))
	for(i in 1:length(allele.frequencies[,1])){
		plot(allele.frequencies[i,], type = "h", main = rownames(allele.frequencies)[i], col = cols[i], axes = FALSE, xlab = "", ylab = "", ylim = c(0, max(allele.frequencies)))
		axis(2)
		mtext("Frequency", side = 2, line = 2.5)
		abline(h = (percent.in.rare.geno/100), col = "red")
		}
	
	for(i in 1:length(geno.freq[,1])){
		plot(geno.freq[i,], type = "h", main = rownames(geno.freq)[i], axes = FALSE, xlab = "", ylab = "", ylim = c(0, max(geno.freq)))
		axis(2)
		mtext("Frequency", side = 2, line = 2.5)
		abline(h = (sqrt(percent.in.rare.geno)/100), col = "red")
		}
		
	dev.off()

	
	
	#if pairs is true, plot the frequency of all pairs
	if(pairs){
		
		#This function gets the frequencies of all pairs of genotypes
		count.pairs <- function(locus.pair){
			probs1 <- geno.table[,locus.pair[1]]
			probs2 <- geno.table[,locus.pair[2]]
			poss.probs <- c(0, 0.5, 1)
			freq.mat <- matrix(NA, 3, 3)
			colnames(freq.mat) <- rownames(freq.mat) <- genotypes
			for(i in 1:length(poss.probs)){
				for(j in 1:length(poss.probs)){
					freq.mat[i,j] <- length(intersect(which(probs1 == poss.probs[i]), which(probs2 == poss.probs[j])))
					}
				}
			# result <- list(freq.mat/sum(freq.mat))
			# names(result) <- paste(colnames(geno.table)[locus.pair[1]], colnames(geno.table)[locus.pair[2]], sep = "_")
			# return(result)

			result <- as.vector(freq.mat)
			return(freq.mat)
			}

		
		#do a check to see how many pairs this will be
		num.pairs <- choose(length(loci), 2)
		
		if(num.pairs > 1000){
			
		}
		pairs.mat <- pairs.matrix(1:length(loci))	#all locus pairs (by index)
		two.way.genos <- pairs.matrix(genotypes, TRUE, TRUE)	#all genotype pairs
		
		#count the frequency of each genotype pair for each locus pair
		pair.freq <- t(apply(pairs.mat, 1, count.pairs))
		cutoff <- dim(geno.table)[1]*(percent.in.rare.geno/100)
		pair.freq.names <- rep(NA, dim(pair.freq)[2])
		
		#plot the freqeuncies by genotype combination
		pdf("Locus_Pair_Counts_per_Genotype_Pair.pdf", width = 15, height = 5)
		for(i in 1:length(pair.freq[1,])){
			plot(pair.freq[,i], type = "l", main = paste(two.way.genos[i,1], "and", two.way.genos[i,2]), ylab = "Number of Animals", xlab = "Locus Pair", ylim = c(0, max(pair.freq)))
			pair.freq.names[i] <- paste(two.way.genos[i,1], two.way.genos[i,2], sep = "_", collapse = "_")
			abline(h = cutoff, col = "red")
			}
		dev.off()
		
		
		all.locus.pairs <- pairs.matrix(loci)
		final.pair.freq <- cbind(all.locus.pairs, pair.freq)
		colnames(final.pair.freq) <- c("marker1", "marker2", pair.freq.names)
		write.table(final.pair.freq, paste("Genotype.Frequency.Locus.Pairs.", label, ".txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")

		
		}

}