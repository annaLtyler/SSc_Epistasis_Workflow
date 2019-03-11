#look for endocrine enrichment in blocks interacting with a given block
# library(evd); library(Matrix); library(corpcor); library(qpcR); library(fdrtool); library(nlme); library(RColorBrewer); library(shape); library(regress); library(doParallel); library(igraph)
# library(gProfileR)
# #source all the functions we'll need
# source('~/Documents/git_repositories/useful_r_code/source.fun.R')
# source('~/Documents/git_repositories/Scleroderma/filter.snps.netwas.R')
# source.fun(c("cape_addons", "capeRel", "useful_r_code"))
# results.dir <- "~/Documents/Data/Little_Cross/little_cross_data_cape/Body_Comp_IGF_New_P" 
# setwd(results.dir)
# cross <- readRDS("cross.RData")
# cross <- get.network(cross, p.or.q = 0.0005)
# fun.terms <- c("sex", "steroid", "estrogen", "testosterone", "androgen")
# data.obj  <- cross
# block.name = "sex"
# calc.enrichment = TRUE
# n.iter = 100
# block.name = "Final.IGF.1"; fun.terms = c("bone", "osteoblast", "osteoclast", "growth", "ossif")

enrichment.in.blocks <- function(data.obj, block.name, fun.terms, calc.enrichment = TRUE, combine.blocks = TRUE, n.iter = 100){

	#================================================================
	#read in data files
	#downloaded all data from: ftp://ftp.informatics.jax.org/pub/reports/index.html#marker
	#gene locations: MGI_MRK_Coord.rpt 
	#gene-associated GO terms: gene_association.mgi (deleted all ' characters)
	#all GO terms: go_terms.mgi (deleted all ' characters)
	#================================================================	
	all.genes <- read.table("~/Documents/Data/Mice/mouseGenes.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
	pc.gene.locale <- which(all.genes[,2] == "protein coding gene")
	all.genes <- all.genes[pc.gene.locale,]
	
	go.genes <- read.table("~/Documents/Data/Mice/mouseGOterms_Genes.txt", sep = "\t", stringsAsFactors = FALSE, skip = 6, fill = TRUE)
	chr.sizes <- read.table("~/Documents/Data/Mice/mouse.chr.sizes.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
	all.go <- read.table("~/Documents/Data/Mice/MouseGOterms.txt", sep = "\t", stringsAsFactors = FALSE)
	#================================================================

	collapsed.net <- cross$collapsed.net
	just.interactions <- collapsed.net[1:min(dim(collapsed.net)), 1:min(dim(collapsed.net))]
	int.block.locale <- which(rownames(just.interactions) == block.name)
	targets.block <- names(which(cross$collapsed.net[,int.block.locale] != 0))
	block.targets <- names(which(cross$collapsed.net[int.block.locale,] != 0))

	#================================================================
	#internal functions
	#================================================================
	translate.go.terms <- function(the.terms){
		term.locale <- match(the.terms, all.go[,2])
		term.descr <- all.go[term.locale,3]
		return(term.descr)
		}


	get.block.functions <- function(go.terms, gene.list){
	
		block.function.list <- vector(mode = "list", length = length(gene.list))
		names(block.function.list) <- names(gene.list)
		
		for(i in 1:length(gene.list)){
			block.genes <- gene.list[[i]]
			gene.locale <- which(go.genes[,3] %in% block.genes)
			gene.block <- go.genes[gene.locale,]
			gene.go.terms <- gene.block[,5]
			go.descr <- translate.go.terms(gene.go.terms)
			
			#find the terms of interest among the GO terms
			int.term.locale <- lapply(go.terms, function(x) grep(x, go.descr))
			int.term.genes <- lapply(int.term.locale, function(x) gene.block[x,3])
			int.term.fun <- lapply(int.term.locale, function(x) go.descr[x])
			names(int.term.genes) <- names(int.term.fun) <- go.terms
			
			term.mat <- NULL
			for(j in 1:length(go.terms)){
				#get all related terms
				term.mat <- rbind(term.mat, unique(cbind(int.term.genes[[j]], int.term.fun[[j]])))
				} #end looping through term counts
			block.function.list[[i]] <- term.mat
			}
		return(block.function.list)
		}	
		
		
		
		find.overlapping.genes <- function(chr, start.pos, end.pos){
			chr.locale <- which(all.genes[,5] == chr)
			chr.table <- all.genes[chr.locale,]
			gene.start.in.region <- intersect(which(as.numeric(chr.table[,6]) >= start.pos), which(as.numeric(chr.table[,6]) <= end.pos))
			gene.end.in.region <- intersect(which(as.numeric(chr.table[,7]) >= start.pos), which(as.numeric(chr.table[,7]) <= end.pos))
			#genes in the region have either start or end in the region
			inc.gene.locale <- unique(c(gene.start.in.region, gene.end.in.region))
			inc.genes <- chr.table[inc.gene.locale,3]
			return(unique(inc.genes))
			}
		
	get.block.location <- function(block.id){
		block.locale <- which(names(data.obj$linkage.blocks.collapsed) == block.id)
		block.markers <- data.obj$linkage.blocks.collapsed[[block.locale]]
		block.pos <- get.marker.location(data.obj, block.markers)
		block.chr <- as.numeric(get.marker.chr(data.obj, block.markers)[1])
		results <- c(block.chr, min(block.pos), max(block.pos))
		return(results)
		}
		
		
	get.enrichment <- function(block.names, filename){
		block.enrichment <- NULL
		just.terms <- matrix(NA, ncol = 2, nrow = length(block.names))
		just.terms[,1] <- block.names
		for(i in 1:length(block.names)){
			report.progress(i, length(block.names))
			block.pos <- get.block.location(block.names[i])
			block.genes <- find.overlapping.genes(block.pos[1], block.pos[2], block.pos[3])
			one.block.enrichment <- gprofiler(block.genes, organism = "mmusculus")
			just.terms[i,2] <- paste(one.block.enrichment[,"term.name"], collapse = ", ")
			one.block.enrichment <- cbind(rep(block.names[i], dim(one.block.enrichment)[1]), one.block.enrichment)
			block.enrichment <- rbind(block.enrichment, one.block.enrichment)
			}
		colnames(block.enrichment)[1] <- "block.name"
		cat("\nWriting table", filename, "\n")
		write.table(block.enrichment, filename, quote = FALSE, sep = "\t", row.names = FALSE)
		cat("Writing table", paste("Just.Terms.", filename, sep = ""), "\n")
		write.table(just.terms, paste("Just.Terms.", filename, sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
		}
	
	block.p <- function(true.num, rnd.dist){
		p.vals <- rep(NA, length(true.num))
		for(i in 1:length(true.num)){
			# hist(rnd.dist[,i]);abline(v = true.num[i], col = "red")
			p.vals[i] <- length(which(rnd.dist[,i] >= true.num[i]))/dim(rnd.dist)[1]
			}
		return(p.vals)
		}

	#================================================================
	#end internal functions
	#================================================================
	

	if(calc.enrichment){
		cat("Checking block enrichment...\n")
		get.enrichment(targets.block, filename = paste("Enrichment.", block.name, ".sources.txt", sep = ""))
		get.enrichment(block.targets, filename = paste("Enrichment.", block.name, ".targets.txt", sep = ""))
		}

	#for each block
	#cound genes with sex-related terms
	#then select random genomic location
	#find overlapping genes
	#count genes with sex-related terms

	all.block.interactors <- unique(c(targets.block, block.targets))
	non.interactors <- setdiff(rownames(data.obj$collapsed.net), all.block.interactors)
	
	true.gene.list <- vector(mode = "list", length = length(all.block.interactors))
	names(true.gene.list) <- all.block.interactors
	for(i in 1:length(all.block.interactors)){
		block.pos <- get.block.location(all.block.interactors[i])
		true.gene.list[[i]] <- find.overlapping.genes(block.pos[1], block.pos[2], block.pos[3])
		}

	alt.gene.list <- vector(mode = "list", length = length(non.interactors))
	names(alt.gene.list) <- non.interactors	
	for(i in 1:length(non.interactors)){
		block.pos <- get.block.location(non.interactors[i])
		alt.gene.list[[i]] <- find.overlapping.genes(block.pos[1], block.pos[2], block.pos[3])
		}

	if(combine.blocks){
		combined.gene.list <- list(unique(unlist(true.gene.list)))
		true.gene.list <- combined.gene.list
		names(true.gene.list) <- "all.blocks"
		combined.alt.list <- list(unique(unlist(alt.gene.list)))
		names(combined.alt.list) <- "non.interactors"
		alt.gene.list <- combined.alt.list
		}


	true.gene.function <- get.block.functions(go.terms = fun.terms, gene.list = true.gene.list)
	true.gene.table <- list2Matrix(true.gene.function, preserve.dim = TRUE)
	cat("Writing table Genes.Related.to.Terms.Block.",block.name, ".txt", sep = "", "\n")
	write.table(true.gene.table, paste("Genes.Related.to.Terms.Block.",block.name, ".txt", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE)
	alt.gene.function <- get.block.functions(go.terms = fun.terms, gene.list = alt.gene.list)	

	num.related.genes <- unlist(lapply(true.gene.function, function(x) length(unique(x[,1]))))
	just.gene.table <- matrix(NA, ncol = 3, nrow = length(num.related.genes))
	just.gene.table[,1] <- names(num.related.genes)
	just.gene.table[,2] <- num.related.genes
	for(i in 1:length(num.related.genes)){
		block.locale <- grep(names(num.related.genes)[i], rownames(true.gene.table))
		just.gene.table[i,3] <- paste(unique(true.gene.table[block.locale,1]), collapse = ", ")
		}
	cat("Writing table Gene.Num.Related.to.Terms.Block.",block.name, ".txt", sep = "", "\n")
	write.table(just.gene.table, paste("Gene.Num.Related.to.Terms.Block.",block.name, ".txt", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
	
	
	all.block.stats <- matrix(unlist(lapply(all.block.interactors, get.block.location)), ncol = 3, byrow = TRUE)
	rownames(all.block.stats) <- all.block.interactors
	all.block.size <- all.block.stats[,3] - all.block.stats[,2]

	num.related.genes.table <- matrix(NA, ncol = length(all.block.interactors), nrow = n.iter)
	#select random regions on random chromosomes corresponding to each block
	for(n in 1:n.iter){
		report.progress(n, n.iter)
		rnd.gene.list <- vector(mode = "list", length = length(all.block.interactors))
		names(rnd.gene.list) <- 1:length(all.block.interactors)

		#sample genes in a comparable genomic region
		for(i in 1:length(all.block.interactors)){	
			block.size <- all.block.size[i]
			rnd.chr <- sample(19, 1)
			chr.size <- chr.sizes[which(chr.sizes[,"chromosome"] == rnd.chr),2]*10^6
			max.start <- chr.size - block.size
			rnd.start <- sample(max.start, 1)
			rnd.gene.list[[i]] <- find.overlapping.genes(chr = rnd.chr, start.pos = rnd.start, end.pos = rnd.start+block.size)
			}
		
		#check these genes for enrichment in the selected terms
		rnd.gene.function <- get.block.functions(go.terms = fun.terms, gene.list = rnd.gene.list)
		num.related.genes.table[n,] <- unlist(lapply(rnd.gene.function, function(x) length(unique(x[,1]))))		
		}
	
	
		p.table <- block.p(num.related.genes, num.related.genes.table)
		names(p.table) <- all.block.interactors
	
		cat("\nP values of enrichment in each block are the following:\n")
		print(p.table)
	}