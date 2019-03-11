plot.LD.around.SNP <- function(geno.obj, snp.name, mart, snp.db, upstream.buffer = 1e6, downstream.buffer = 1e6, ask.before.plotting = TRUE, test.type = c("correlation", "shared.geno")){


	gene.info <- gene.near.snp(snp.name, mart, snp.db, upstream.buffer, downstream.buffer)	
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
	
	snp.pos <- paste0(snp.info[1,3], ":", snp.info[1,4])
	chr.region <- paste0(snp.info[1,3], ":", as.numeric(snp.info[1,4])-upstream.buffer, ":", as.numeric(snp.info[1,4])+downstream.buffer) 
	
	plot.LD(geno.obj, chr.region = chr.region, gene.regions = gene.regions, snp.pos = snp.pos, snp.name = snp.name, ask.before.plotting, test.type)
	
	}