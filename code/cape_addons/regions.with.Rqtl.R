#This function uses R/qtl to impute genotypes and find
#usable regions for markers. We can use these regions to
#look for candidate genes


regions.with.Rqtl <- function(filename, chr = NULL, scan.pheno = NULL, covar = NULL){
	
	require(qtl)
	
	data.obj <- read.cross(format = "csv", file = filename)
	num.pheno <- dim(data.obj$pheno)[2]
	
	if(is.null(scan.pheno)){
		pheno.to.scan <- 1:num.pheno
		}else{
		pheno.to.scan <- which(dimnames(data.obj$pheno)[[2]] %in% scan.pheno)
		if(length(pheno.to.scan) < length(scan.pheno)){
			didnt.find <- which(scan.pheno %in% dimnames(data.obj$pheno)[[2]] == FALSE)
			message("I couldn't find the following phenotypes:")
			cat(scan.pheno[didnt.find], sep = "\n")
			}
		}
	
	if(is.null(chr)){
		chr <- names(data.obj$geno)
		}
		
	
	data.obj <- calc.genoprob(data.obj, step = 1, off.end = 10)
	
	if(is.null(covar)){
		onescan <- scanone(data.obj, chr = chr, pheno.col = pheno.to.scan, model = "normal", method = "hk")
		}else{
		covar.locale <- which(names(data.obj$pheno) %in% covar)
		covar.col <- data.obj$pheno[,covar.locale]
		onescan <- scanone(data.obj, chr = chr, pheno.col = pheno.to.scan, model = "normal", method = "hk", addcovar = covar.col)
		}
	
	pdf(file = "qtl.scan.pdf", width = 10, height = 3*length(scan.pheno))
	par(mfrow = c(length(scan.pheno), 1))
	for(i in 1:length(scan.pheno)){
		plot(onescan, lodcol = i)
		}
		dev.off()


	#find the intervals for each chromosome and each phenotype
	interval.list.lods <- vector(mode = "list", length = length(scan.pheno))
	interval.list.pos <- vector(mode = "list", length = length(scan.pheno))
	names(interval.list.lods) <- names(interval.list.pos) <- scan.pheno
	start <- 1
	for(j in 1:length(scan.pheno)){
		interval.lod.table <- matrix(NA, ncol = 3, nrow = length(chr))
		interval.pos.table <- matrix(NA, ncol = 3, nrow = length(chr))
		rownames(interval.table) <- chr; colnames(interval.table) <- c("min", "peak", "max")
		for(i in 1:length(chr)){
			#put the interval for this chromosome and this phenotype into the matrix
			interval <- lodint(onescan, chr[i], drop = 1.5, lodcolumn = j)
			interval.lod <- interval[,scan.pheno[j]]
			interval.pos <- interval[,"pos"]
			if(length(interval.pos) == 3){
				interval.lod.table[i,] <- interval.lod
				interval.pos.table[i,] <- interval.pos
				}
			}
		#put the interval matrix for this phenotype into the list
		interval.list.lods[[j]] <- interval.lod.table
		interval.list.pos[[j]] <- interval.pos.table
		}
	
	#for each chromosome, find the intervals of the maximum peak across all phenotypes
	get.peak.range <- function(chr.num){
		all.lod <- sapply(interval.list.lods, function(x) x[chr.num,])
		all.pos <- sapply(interval.list.pos, function(x) x[chr.num,])
		peak.lod.pheno.locale <- which(all.lod == max(all.lod, na.rm = TRUE), arr.ind = TRUE)[2]
		peak.range <- all.pos[,peak.lod.pheno.locale]
		return(peak.range)
		}
		
	all.peaks <- t(apply(matrix(1:length(chr), nrow = 1), 2, get.peak.range))
	colnames(all.peaks) <- c("min", "peak", "max")
	rownames(all.peaks) <- chr

	neg.vals <- which(all.peaks < 0)
	all.peaks[neg.vals] <- 0

	#put the table into a form that we can load directly into the mouse mapper
	block.names <- paste("block", 1:length(chr), sep = "")
	final.min <- cbind(chr, all.peaks[,1])
	final.max <- cbind(chr, all.peaks[,3])
	rownames(final.min) <- paste(block.names, ".min", sep = "")
	rownames(final.max) <- paste(block.names, ".max", sep = "")

	final.table <- rbind(final.min, final.max)

	write.table(final.table, file = "block.coord.cM.txt", sep = "\t", quote = FALSE, col.names = FALSE)

			
}