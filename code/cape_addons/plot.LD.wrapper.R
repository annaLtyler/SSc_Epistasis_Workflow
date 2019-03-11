#This function plots the proportion
#of pts with shared genotypes at 
# two different SNPs
#chr.region = "6:32084741:32324067"
# gene.regions = c("PBX2"= "6:32184741:32190186", "GPSM3"= "6:32190766:32195523", "NOTCH4" = "6:32194843:32224067")
#snp: rs2451280 -- snp.pos = "2:163111659"; chr.region = c("2:163106659:163116659"); gene.regions = c()
#snp: rs3130933, snp.pos = "6:31132335"; chr.region = "6:31107335:31157335"; gene.regions = c("TCF19" = "6:31126303:31131992", "POU5F1" = "6:31132114:31138451")
# chr.region <- "11:55231406:57588649"; snp.pos = "6:56298726"; plot.genes = TRUE; snp.name = "rs10501351";ask.before.plotting = FALSE; test.type = "correlation"; organism = "human"


plot.LD.wrapper <- function(geno.obj, chr.region, organism = c("mouse", "human"), plot.genes = TRUE, snp.pos = NULL, snp.name = NULL, ask.before.plotting = FALSE, test.type = c("correlation", "r2")){

	require(igraph)
	require(biomaRt)

	organism <- organism[1]
	test.type <- test.type[1]

	if(organism == "mouse"){		
	# lib <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
		lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2017.archive.ensembl.org")
		snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="mmusculus_snp", host = "may2017.archive.ensembl.org")	
	}else{
	# lib <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")	
	lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2017.archive.ensembl.org")
	snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp", host = "may2017.archive.ensembl.org")
	}


	#first find all the genes in the region
	gene.table <- getBM(c("external_gene_name", "entrezgene", "chromosome_name", "start_position", "end_position"), "chromosomal_region", chr.region, lib)
	#filter to genes with entrez IDs
	gene.table <- gene.table[which(!is.na(gene.table[,"entrezgene"])),]
	
	genes <- gene.table[,1]
	gene.regions <- apply(gene.table, 1, function(x) paste0(x[3], ":", x[4], ":", x[5]))
	names(gene.regions) <- genes

	chr.pts <- as.numeric(unlist(strsplit(chr.region, ":")))
	chr <- chr.pts[1]
	region.start <- chr.pts[2]
	region.end <- chr.pts[3]

	region.width <- region.end - region.start

	chr.locale <- which(geno.obj$chromosome == chr)
	after.start <- which(geno.obj$marker.location >= region.start)
	before.end <- which(geno.obj$marker.location <= region.end)
	
	region.locale <- Reduce(intersect, list(chr.locale, after.start, before.end))
	geno.pairs <- pair.matrix(colnames(geno.obj$geno)[region.locale])

	if(ask.before.plotting){
	choice <- readline(paste("The pair matrix has", nrow(geno.pairs), "pairs.\nWould you like to continue (y/n)?\n"))
	if(choice == "n"){stop()}
	}
	
	#=========================================================
	#internal functions
	#=========================================================
		
	get.gene.plot.pts <- function(one.gene.region){
		gene.pts <- as.numeric(unlist(strsplit(one.gene.region, ":")))
		gene.start <- gene.pts[2]
		gene.end <- gene.pts[3]
		rel.start.x <- ((gene.start - region.start)/region.width)*length(region.locale)
		rel.start.y <- length(region.locale) - ((gene.start - region.start)/region.width)*length(region.locale)
		
		rel.end.x <- ((gene.end - region.start)/region.width)*length(region.locale)
		rel.end.y <- length(region.locale) - ((gene.end - region.start)/region.width)*length(region.locale)
		
		return(c(rel.start.x, rel.start.y, rel.end.x, rel.end.y))
		}
	
	get.r2 <- function(pair.num){
		snp1 <- geno.pairs[pair.num,1]
		snp2 <- geno.pairs[pair.num,2]
		snp1.geno <- geno.obj$geno[,snp1]
		snp2.geno <- geno.obj$geno[,snp2]
		
		combo <- table(bin.vector(snp1.geno), bin.vector(snp2.geno))
		
		maf1 <- mean(snp1.geno)
		if(maf1 > 0.5){snp1.geno <- 1 - snp1.geno}

		p1 <- mean(snp1.geno)
		q1 <- 1 - p1
		
		maf2 <- mean(snp2.geno)
		if(maf2 > 0.5){snp2.geno <- 1 - snp2.geno}
		
		p2 <- mean(snp2.geno)
		q2 <- 1 - p2
		
		p12 <- length(intersect(which(snp1.geno > 0), which(snp2.geno > 0)))/length(snp1.geno)
		
		D = p12 - (p1*p2)
		
		if(D < 0){
			Dmin <- max(c(p1*p2*-1, q1*q2*-1))
			}else{
			Dmin <- min(c(p1*q2, q1*p2))
			}
		
		Dprime <- D/Dmin
		
		denom <- sqrt(p1*q1*p2*q2)
		r <- D/denom
		
		return(r^2)
		
		}
	#=========================================================

	if(test.type == "correlation"){
	geno.equal <- unlist(lapply_pb(1:nrow(geno.pairs), function(x) cor(geno.obj$geno[,geno.pairs[x,1]], geno.obj$geno[,geno.pairs[x,2]], use = "complete")))
	}else{
	geno.equal <- unlist(lapply_pb(1:nrow(geno.pairs), function(x) get.r2(x)))
	}
	net <- graph_from_edgelist(geno.pairs, directed = FALSE)
	E(net)$weight <- geno.equal
	
	adj <- as.matrix(as_adjacency_matrix(net, attr = "weight", type = "upper"))
	adj[which(adj == 0)] <- NA

	# quartz(width = 8, height = 5)
	layout(matrix(c(1,2), nrow = 1), widths = c(1, 0.2))
	par(mar = c(2,2,2,2))
	imageWithText(abs(adj), show.text = FALSE, split.at.vals = TRUE, split.points = c(0.25, 0.5, 0.75), col.scale = c("blue","green","orange","red"), global.color.scale = TRUE, grad.dir = "high", main = paste0("Chr ", chr, ": ", region.start, "-", region.end))
	
	if(!is.null(gene.regions)){
		gene.xy <- lapply(gene.regions, get.gene.plot.pts)

		for(i in 1:length(gene.xy)){
			segments(x0 = gene.xy[[i]][1], y0 = gene.xy[[i]][2], x1 = gene.xy[[i]][3], y1 = gene.xy[[i]][4], lwd = 1)
			label.x <- mean(gene.xy[[i]][c(1,3)])
			label.y <- mean(gene.xy[[i]][c(2,4)])
			text(label.x-(length(region.locale)*0.01), label.y, names(gene.regions)[i], adj = 1, cex = 0.5)
			}
		}
		
	if(!is.null(snp.pos)){
		snp.info <- as.numeric(strsplit(snp.pos, ":")[[1]])[2]
		snp.x <- (snp.info - region.start)/region.width*length(region.locale)
		snp.y <- length(region.locale) - ((snp.info - region.start)/region.width)*length(region.locale)
		points(snp.x, snp.y, pch = "*")
		text(snp.x, snp.y, labels = snp.name, adj = 1, cex = 0.5, col = "red")
		}
		
	par(mar = c(2,2,2,2))
	imageWithTextColorbar(abs(adj), split.at.vals = TRUE, split.points = c(0.25, 0.5, 0.75), col.scale = c("blue","green","orange","red"), global.color.scale = TRUE, grad.dir = "high")

	
}