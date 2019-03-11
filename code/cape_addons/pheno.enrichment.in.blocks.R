#This function is similar to enrichment.in.blocks, but uses phenotype
#terms instead of GO processes
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

pheno.enrichment.in.blocks <- function(data.obj, block.name, fun.terms, calc.enrichment = TRUE, n.iter = 100){

	#================================================================
	#read in data files
	#downloaded all data from: ftp://ftp.informatics.jax.org/pub/reports/index.html#marker
	#gene locations: MGI_MRK_Coord.rpt 
	#gene-associated GO terms: gene_association.mgi (deleted all ' characters)
	#all GO terms: go_terms.mgi (deleted all ' characters)
	#================================================================	
	all.genes <- read.table("~/Documents/Data/Mice/mouseGenes.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
	pheno.go <- read.table("~/Documents/Data/Mice/mousePhenoGO.txt", sep = "\n", stringsAsFactors = FALSE, skip = 6, fill = TRUE)
	chr.sizes <- read.table("~/Documents/Data/Mice/mouse.chr.sizes.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
	pheno.genes <- read.table("~/Documents/Data/Mice/mousePhenoGenes.txt", sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
	#================================================================

	#================================================================
	#process pheno.go to get separate terms
	#================================================================
	term.locale <- which(pheno.go == "[Term]")
	pheno.go.list <- vector(mode = "list", length = length(term.locale))
	for(i in 1:(length(term.locale)-1)){
		pheno.go.list[[i]] <- pheno.go[(term.locale[i]+1):(term.locale[(i+1)]-1),1]
		}
	pheno.go.list[[length(pheno.go.list)]] <- pheno.go[(term.locale[length(term.locale)]+1):dim(pheno.go)[1],1]
	just.mp.terms <- unlist(lapply(lapply(lapply(pheno.go.list, function(x) x[1]), function(x) strsplit(x, "id: ")), function(x) x[[1]][2]))
	#================================================================
	
	
	collapsed.net <- data.obj$collapsed.net
	just.interactions <- collapsed.net[1:min(dim(collapsed.net)), 1:min(dim(collapsed.net))]
	int.block.locale <- which(rownames(just.interactions) == block.name)
	targets.block <- names(which(just.interactions[,int.block.locale] != 0))
	block.targets <- names(which(just.interactions[int.block.locale,] != 0))


	#================================================================
	#internal functions
	#================================================================
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

	get.gene.ids <- function(gene.names){
		gene.locale <- match(gene.names, all.genes[,3])
		gene.ids <- rownames(all.genes)[gene.locale]
		result <- cbind(gene.ids, gene.names)
		return(result)	
		}
	
	get.gene.names <- function(gene.ids){
		gene.locale <- match(gene.ids, rownames(all.genes))
		gene.names <- all.genes[gene.locale,3]
		result <- cbind(gene.ids, gene.names)
		return(result)	
		}	
	
	get.block.functions <- function(gene.list){
	
		block.function.list <- vector(mode = "list", length = length(gene.list))
		names(block.function.list) <- names(gene.list)
		
		for(i in 1:length(gene.list)){
			block.genes <- gene.list[[i]]
			#find all the MP terms associated with these genes
			block.gene.ids <- get.gene.ids(block.genes)
			#find the genes in the phenotype table
			gene.locale <- which(pheno.genes[,6] %in% block.gene.ids[,1])
			mp.terms <- unique(pheno.genes[gene.locale,4])
			mp.table <- translate.mp.terms(mp.terms)
			block.function.list[[i]] <- mp.table
			}
			return(block.function.list)
		}	

	translate.mp.terms <- function(mp.terms){
		all.go <- NULL
		term.locale <- match(mp.terms, just.mp.terms)
		term.id.list <- lapply(pheno.go.list[term.locale], function(x) c(x[1], x[2], x[3]))
		term.ids <- unlist(lapply(term.id.list, function(x) strsplit(x[1], "id: ")[[1]][2]))
		term.names <- unlist(lapply(strsplit(unlist(lapply(term.id.list, function(x) x[2])), "name: "), function(x) x[2]))
		term.def <- unlist(lapply(strsplit(unlist(lapply(term.id.list, function(x) x[3])), "def: "), function(x) x[2]))
		all.go <- rbind(all.go, cbind(term.ids, term.names, term.def))
		return(all.go)
		}
	
	filter.function.table <- function(go.terms, fun.table){
		term.locale <- NULL
		for(i in 1:length(go.terms)){
			locale1 <- grep(go.terms[i], fun.table[,2], ignore.case = TRUE)
			locale2 <- grep(go.terms[i], fun.table[,3], ignore.case = TRUE)
			term.locale <- unique(c(term.locale, locale1, locale2))
			}
		filtered.table <- fun.table[term.locale,]
		return(filtered.table)
		}
		
	get.assoc.genes <- function(go.table, gene.block.table){
		row.locale <- which(pheno.genes[,4] %in% go.table[,1])
		term.gene.ids <- unique(unlist(strsplit(pheno.genes[row.locale,6], ",")))
		gene.names <- get.gene.names(term.gene.ids)
		#filter out the NAs, these tend to be alleles of genes and not genes themselves
		not.na.locale <- which(!is.na(gene.names[,2]))
		assoc.genes <- gene.names[not.na.locale,]
		assoc.block.genes <- intersect(assoc.genes[,1], gene.block.table)
		return(assoc.block.genes)
		}

	#================================================================

	#get the stats for the interacting blocks

	all.block.interactors <- unique(c(targets.block, block.targets))
	true.gene.list <- vector(mode = "list", length = length(all.block.interactors))
	names(true.gene.list) <- all.block.interactors
	for(i in 1:length(all.block.interactors)){
		block.pos <- get.block.location(all.block.interactors[i])
		true.gene.list[[i]] <- find.overlapping.genes(block.pos[1], block.pos[2], block.pos[3])
		}

	#get all the functions/phenotypes in the block
	true.gene.function <- get.block.functions(gene.list = true.gene.list)
	#filter the functions down to those relevant to our terms
	filtered.functions <- lapply(true.gene.function, function(x) filter.function.table(fun.terms, x))
	#get all the genes in the region associated with these MP terms
	mp.associated.genes <- vector(mode = "list", length = length(all.block.interactors))
	names(mp.associated.genes) <- all.block.interactors
	for(i in 1:length(filtered.functions)){
		mp.associated.genes[[i]] <- get.assoc.genes(go.table = filtered.functions[[i]], gene.block.table = true.gene.list[[i]])	
		}
	
	
	true.gene.table <- list2Matrix(true.gene.function, preserve.dim = TRUE)
	cat("Writing table Genes.Related.to.Terms.Block.",block.name, ".txt", sep = "", "\n")
	write.table(true.gene.table, paste("Genes.Related.to.Terms.Block.",block.name, ".txt", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE)
	
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