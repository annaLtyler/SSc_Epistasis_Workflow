
plot.allele.freq.biallele <- function(geno.table, label = NULL, pairs = FALSE, percent.in.rare.geno = 5){

	loci <- colnames(geno.table)
	genotypes <- c("AA", "AB", "BB")
	alleles <- c("A", "B")
	
	#given a 0, 0.5, 1 set of probabilities
	#for allele A at a given locus, this 
	#function returns the frequency of each
	#allele over all individuals
	count.alleles <- function(prob.list){
		num.A <- length(which(prob.list == 0))*2 + length(which(prob.list == 0))*0.5
		num.B <- length(which(prob.list == 1))*2 + length(which(prob.list == 0))*0.5
		total.alleles <- num.A + num.B
		freq.a.b <- c((num.A/total.alleles), (num.B/total.alleles))
		return(freq.a.b)
		}
	
	#This function counts the frequency of all genotypes
	#at each locus
	get.geno.freq <- function(prob.list){
		AA.freq <- length(which(prob.list == 0))
		AB.freq <- length(which(prob.list == 0.5))
		BB.freq <- length(which(prob.list == 1))
		total.geno <- AA.freq + AB.freq + BB.freq
		freq.geno <- c(AA.freq/total.geno, AB.freq/total.geno, BB.freq/total.geno)
		return(freq.geno)
		}
	
	allele.freq <- apply(geno.table, 2, count.alleles)
	geno.freq <- apply(geno.table, 2, get.geno.freq)
	
	
	pdf(paste("Allele.Frequencies.", label, ".pdf", sep = ""), width = 20, height = 5)
	# quartz(width = 20, height = 5)
	par(mar = c(4, 2, 4, 0))
	cols <- gray.colors(2)
	bar.mids <- barplot(allele.freq, axes = FALSE, ylim = c(0, 1), names.arg = rep("", length(loci)), main = "Allele Frequencies", xpd = TRUE)
	legend(-4,1, legend = alleles, fill = rev(cols))
	axis(1, at = bar.mids, labels = FALSE)
	text(bar.mids, par("usr")[3] - 0.05, srt = 45, adj = 1, labels = loci, xpd = TRUE, cex = 0.7)

	# quartz(width = 20, height = 5)
	par(mar = c(4, 2, 4, 0))
	bar.mids <- barplot(geno.freq, axes = FALSE, ylim = c(0, 1), names.arg = rep("", length(loci)), main = "Genotype Frequencies", xpd = TRUE)
	cols <- gray.colors(3)
	legend(-5, 1, legend = genotypes, fill = rev(cols))
	axis(1, at = bar.mids, labels = FALSE)
	text(bar.mids, par("usr")[3] - 0.05, srt = 45, adj = 1, labels = loci, xpd = TRUE, cex = 0.7)

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