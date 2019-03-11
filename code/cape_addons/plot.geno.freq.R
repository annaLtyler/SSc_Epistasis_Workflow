#This function is a cape add-on that plots genotype frequencies
#using the geno table in the data object
#For the pairwise genotype frequencies, this script plots only
#Those frequencies above the given threshold

plot.geno.freq <- function(data.obj, min.thresh = 0.1){
	
	library(RColorBrewer)
	library(lattice)
	mypal <- brewer.pal(3, "Blues")
	# barplot(rep(1, length(mypal)), col = mypal)

	geno <- data.obj$geno
	
	#first look at genotype frequencies at single loci
	
	#======================================================
	#functions
	#======================================================
	
	get.allele.freq <- function(genotypes){
		count.table <- table(genotypes)
		result <- list(count.table/length(genotypes), names(count.table))
		names(result) <- c("freq", "alleles")
		return(result)
		}
		
	plot.allele.freq <- function(freq.mat, main = "Genotype Frequencies"){
		alleles <- rownames(freq.mat)
		bar.col <- mypal[1:3]
		bar.mids <- barplot(freq.mat, names.arg = rep("", dim(geno)[2]), main = main, col = bar.col)
		axis(2)
		axis(1, at = bar.mids, labels = NA)
		par(xpd = TRUE)
		text(x = bar.mids, y = par("usr")[3]-(max(freq.mat)*0.12), labels = colnames(geno), srt = 90)
		legend("topright", legend = alleles, fill = bar.col[1:length(alleles)])	
		}

	get.pair.freq <- function(pair){
		geno1 <- geno[,pair[1]]
		geno2 <- geno[,pair[2]]
		count.mat <- table(geno1, geno2)/dim(geno)[1]
		return(list(count.mat))
		}
		
	get.genotype.count <- function(count.mat, genotype){
		return(count.mat[as.character(genotype)[1], as.character(genotype[2])])
		}
		
	get.pair.counts <- function(count.table, allele.pair){
		loc1 <- match(allele.pair, rownames(count.table[[1]]))
		loc2 <- match(allele.pair, colnames(count.table[[1]]))
		if(length(which(is.na(c(loc1, loc2)))) == 0){
			return(count.table[[1]][as.vector(allele.pair[1]), as.vector(allele.pair[2])])
			}else{
			return(NA)
			}
		}
	
	plot.pair.counts <- function(count.mat, label){
		layout(matrix(1:2, ncol = 2), widths = 5,1)
		par(mar = c(3,2.5,5,2))
		image(count.mat, axes = FALSE, main = paste("Genotype:", paste(label, collapse = ", ")), col=ColorRamp)
		image(1, ColorLevels,matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp,xlab="",ylab="",xaxt="n")
		}
	#======================================================
	#======================================================
	
	
	allele.freq <- apply(geno, 2, get.allele.freq)
	all.alleles <- unique(unlist(lapply(allele.freq, function(x) x$alleles)))
	allele.mat <- matrix(NA, nrow = length(all.alleles), ncol = length(allele.freq))
	rownames(allele.mat) <- all.alleles
	for(i in 1:length(allele.freq)){
		spec.alleles <- allele.freq[[i]]$alleles
		allele.mat[spec.alleles, i] <- allele.freq[[i]][[1]]
		}
	allele.mat[is.na(allele.mat)] <- 0
	pdf("Genotype.Frequencies.pdf", width = 14, height = 6)
	plot.allele.freq(allele.mat)
	
	norm.allele.freq <- apply(allele.mat, 2, function(x) x/sum(x))
	plot.allele.freq(norm.allele.freq, main = "Normalized Genotype Frequencies")
	dev.off()

	#Now look at all the pair frequencies
	pair.mat <- pair.matrix(1:dim(geno)[2])


	#set some plotting parameters
	ColorRamp <- rgb( seq(0,1,length=256),  # Red
	                   seq(0,1,length=256),  # Green
	                   seq(1,0,length=256))  # Blue
	
	
	pairwise.counts <- apply(pair.mat, 1, get.pair.freq)
	ColorLevels <- seq(min(unlist(pairwise.counts), na.rm = TRUE), max(unlist(pairwise.counts), na.rm = TRUE), length=length(ColorRamp))
	
	#make a table out of the counts
	all.alleles <- unique(unlist(lapply(pairwise.counts, function(x) colnames(x[[1]]))))




	allele.pairs <- pair.matrix(all.alleles, ordered = TRUE, self.pairs = TRUE)	

	pdf("Genotype.Pairs.Frequency.pdf")
	for(i in 1:length(allele.pairs[,1])){ 
		#for each of the genotype pairs, make a matrix that has the counts
		#for this genotype over all pairs of loci
		all.pairs <- matrix(NA, ncol = dim(geno)[2], nrow = dim(geno)[2])
		#get the counts for each pair
		all.counts <- lapply(pairwise.counts, function(x) get.pair.counts(x, allele.pairs[i,]))
		#put these counts into the appropriate places in the matrix
		for(j in 1:length(pair.mat[,1])){
			all.pairs[pair.mat[j,1], pair.mat[j,2]] <- all.counts[[j]]
			all.pairs[pair.mat[j,2], pair.mat[j,1]] <- all.counts[[j]]
			}
			all.pairs[all.pairs <= min.thresh] <- NA
		if(length(which(!is.na(all.pairs))) > 0){
			plot.pair.counts(all.pairs/dim(geno)[1], allele.pairs[i,])
			}
		}
	dev.off()

	




}