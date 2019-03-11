#This function plots allele frequencies around a SNP


plot.allele.freq <- function(geno.obj, snp.name, organism = c("human", "mouse"), upstream.buffer = 5e5, downstream.buffer = 5e5){
	
	
	gene.info <- gene.near.snp(snp.name, organism, upstream.buffer, downstream.buffer)	
	snp.info <- gene.info[[1]]
	
	genes <- gene.info[[2]]
	gene.locale <- which(!is.na(genes[,1]))
	
	if(length(gene.locale) > 0){
		genes <- genes[gene.locale,,drop=FALSE]
		genes <- condense.table(genes, condense.by = 1, col.to.collapse = c(2,3,4,5), col.to.concat = c(6,7))
		gene.regions <- apply(genes, 1, function(x) paste0(x[3], ":", x[4], ":", x[5]))
		names(gene.regions) <- genes[,2]
		}else{
		gene.regions <- NULL
		}
	
	
	get.gene.plot.pts <- function(one.gene.region){
		gene.pts <- as.numeric(unlist(strsplit(one.gene.region, ":")))
		gene.start <- gene.pts[2]
		gene.end <- gene.pts[3]
		rel.start.x <- ((gene.start - region.start)/region.width)*length(region.locale)
		rel.end.x <- ((gene.end - region.start)/region.width)*length(region.locale)
		
		return(c(rel.start.x, rel.end.x))
		}


	
	snp.pos <- paste0(snp.info[1,3], ":", snp.info[1,4])
	chr.region <- paste0(snp.info[1,3], ":", as.numeric(snp.info[1,4])-upstream.buffer, ":", as.numeric(snp.info[1,4])+downstream.buffer) 

	chr.pts <- as.numeric(unlist(strsplit(chr.region, ":")))
	chr <- chr.pts[1]
	region.start <- chr.pts[2]
	region.end <- chr.pts[3]

	region.width <- region.end - region.start

	chr.locale <- which(geno.obj$chromosome == chr)
	after.start <- which(geno.obj$marker.location >= region.start)
	before.end <- which(geno.obj$marker.location <= region.end)
	
	region.locale <- Reduce(intersect, list(chr.locale, after.start, before.end))

	allele.freq <- apply(geno.obj$geno[,region.locale], 2, mean)
	maj.allele <- which(allele.freq > 0.5)
	if(length(maj.allele) > 0){
		allele.freq[maj.allele] <- 1 - allele.freq[maj.allele]
		}
	
	plot(allele.freq, type = "h", axes = FALSE, ylim = c(0, 0.5), ylab = "Minor Allele Frequency", xlab = "")
		if(!is.null(gene.regions)){
		
		gene.x <- lapply(gene.regions, get.gene.plot.pts)
		par(xpd = TRUE)
		for(i in 1:length(gene.x)){
			segments(x0 = gene.x[[i]][1], y0 = -0.01, x1 = gene.x[[i]][2], y1 = -0.01, lwd = 1)
			label.x <- mean(gene.x[[i]][c(1,2)])
			label.y <- -0.02
			text(label.x-(length(region.locale)*0.01), label.y, names(gene.regions)[i], adj = 1, cex = 0.5, srt = 90)
			}
		par(xpd = FALSE)
		}

	if(!is.null(snp.pos)){
		par(xpd = TRUE)
		snp.info <- as.numeric(strsplit(snp.pos, ":")[[1]])[2]
		snp.x <- (snp.info - region.start)/region.width*length(region.locale)
		snp.y <- -0.02
		points(snp.x, snp.y, pch = "*", col = "red")
		text(snp.x, snp.y, labels = snp.name, adj = 1, cex = 0.5, col = "red", srt = 90)
		par(xpd = FALSE)
		}
	axis(2)
	
}