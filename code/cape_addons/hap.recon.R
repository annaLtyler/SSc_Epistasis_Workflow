#This script takes in a data object and does haplotype
#reconstruction for each individual using maximum 
#parsimony


hap.recon <- function(data.obj){
	
	chromosomes <- data.obj$chromosome
	u_chr <- unique(chromosomes)
		
	chr.mat <- consec.pairs(chromosomes)
	break.point.mat <- matrix(0, nrow = (length(chromosomes)-1), ncol = dim(data.obj$geno)[1])
	
	plot.width = max((dim(data.obj$geno)[1]*0.1), 5)
	plot.height <- max((length(u_chr)*0.25), 5)
	quartz(width = plot.width, height = plot.height)
	plot.new()
	plot.window(xlim = c(1,dim(data.obj$geno)[1]), ylim = c(1,dim(data.obj$geno)[2]))
	par(mar = c(0,0,0,0))
	
	for(i in 1:dim(data.obj$geno)[1]){
		ind <- data.obj$geno[i,]
	
		cons.mat <- consec.pairs(ind)
		
		#find all the breakpoints for each chromosome
		for(ch in 1:length(u_chr)){
			chr.locale <- intersect(which(chr.mat[,1] == u_chr[ch]), which(chr.mat[,2] == u_chr[ch]))
			na.points <- intersect(which(is.na(cons.mat[chr.locale,1])), which(is.na(cons.mat[chr.locale,2])))
			break.points <- which(apply(cons.mat[chr.locale,,drop=FALSE], 1, function(x) x[1] != x[2]))
			if(length(na.points) > 0){
				break.points <- setdiff(break.points, na.points)
				}
			break.point.mat[chr.locale[break.points],i] <- 1
			}
	
		
		hap1 <- rep(NA, length(ind))
		hap2 <- rep(NA, length(ind))
		hap1.col <- hap1
		hap2.col <- hap2
		
		#find all the homozygous locations
		aa.locale <- which(ind == 0)
		hap1[aa.locale] <- "A"
		hap2[aa.locale] <- "A"	
		bb.locale <- which(ind == 1)
		hap1[bb.locale] <- "B"
		hap2[bb.locale] <- "B"	
		#fill in the heterozygous locations
		ab.locale <- which(ind == 0.5)
		hap1[ab.locale] <- "A"
		hap2[ab.locale] <- "B"
		
		hap1.col[which(hap1 == "A")] <- "red"
		hap1.col[hap1 == "B"] <- "purple"
		hap2.col[hap2 == "A"] <- "red"
		hap2.col[hap2 == "B"] <- "purple"
	
	
		points(rep(i-0.2, length(ind)), 1:length(ind), col = hap1.col, type = "p", lwd = 3, pch = 16, cex = 0.5)
		points(rep(i+0.2, length(ind)), 1:length(ind), col = hap2.col, type = "p", lwd = 3, pch = 16, cex = 0.5)
	
		}	
	
	par(xpd = TRUE)
	text(x = 1:dim(data.obj$geno)[1], y = rep(-4, dim(data.obj$geno)[2]), labels = 1:dim(data.obj$geno)[1], cex = 0.6, srt = 90)

	#mark the chromosome break points
	for(ch in u_chr){
		max.chr <- max(which(chromosomes == ch))
		min.chr <- min(which(chromosomes == ch))
		abline(h = max.chr+0.01)
		text(-2, mean(c(max.chr, min.chr)), ch, cex = 0.5)
		}

	# rf <- rowSums(break.point.mat)/dim(data.obj$geno)[2]
	# plot(rf, type = "h")
	
}