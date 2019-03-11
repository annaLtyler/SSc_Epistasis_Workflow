#code for merging cape and capeDO
#using SSc data
#make sure that if you want to run the singlescan
#there are no singlescan objects in the results directory
#otherwise, this code reads that in, and skips the singlescan
#!!! Remember to place the gene list in the results directory!!! 
testing = FALSE
build.cape.data <- FALSE
exp.name <- "single.perm.for.real"
geno.coding <- "Dominant"
parameter.file <- "cape.parameters.txt"
#===============================================================
# check whether we are on the cluster or home machine
#===============================================================

	cur.dir <- getwd()
	test.dir <- strsplit(cur.dir, "home")
	on.cluster = FALSE
	if(length(test.dir[[1]]) > 1){on.cluster = TRUE}
	
	if(on.cluster){		
		base.dir <- "/home/atyler/SSc"
		r.dir <- "/home/atyler/cape/"
		results.base <- "/home/atyler/"
		cape.do.dir <- "/home/atyler/code/cape"
		cape.addons.dir <- "/home/atyler/code/cape_addons"
		useful.scripts.dir <- "/home/atyler/code/useful_r_code"
		results.dir <- paste0("/home/atyler/SSc/Results/", exp.name)
		init.dir <- "/home/atyler/SSc/data/"
		n.cores = 20
		}else{
		genotype.dir <- "~/Documents/Data/Scleroderma/GWAS_documents/geno/"
		phenotype.dir <- "~/Documents/Data/Scleroderma/GWAS_documents/pheno/phenotypes_temp"
		cape.do.dir <- "~/Documents/git_repositories/cape-r/capempp"
		cape.addons.dir <- "~/Documents/git_repositories/cape_addons"
		useful.scripts.dir <- "~/Documents/git_repositories/useful_r_code"
		results.dir <- paste0("~/Documents/Data/Scleroderma/Results/", exp.name)
		init.dir <- "~/Documents/Data/Scleroderma/project_data/"
		n.cores = 4		
		}

if(!file.exists(results.dir)){system(paste("mkdir", results.dir))}

#===============================================================
# install and load all the necessary libraries
needed.packages <- c("evd", "Matrix", "fdrtool", "shape", "corpcor", "RColorBrewer", "doParallel", "foreach", "caTools", "stringr", "abind", "propagate", "biomaRt")

if(on.cluster){
	all.packages <- installed.packages(lib.loc = "/home/atyler/R/x86_64-pc-linux-gnu-library/3.4")
	
	package.locale <- match(needed.packages, all.packages[,1])
	missing.packages <- needed.packages[which(is.na(package.locale))]
	
	if(length(missing.packages) > 0){install.packages(missing.packages, lib = "/home/atyler/R/x86_64-pc-linux-gnu-library/3.4", repos = "https://cloud.r-project.org")}
	}

for(i in 1:length(needed.packages)){library(needed.packages[i], character.only = TRUE)}
#===============================================================


#===============================================================
#source all the code from needed packages
#===============================================================
all.fun <- list.files(pattern = ".R", path = cape.do.dir, full.names = TRUE)
for(i in 1:length(all.fun)){source(all.fun[i])}


all.fun <- list.files(pattern = ".R", path = cape.addons.dir, full.names = TRUE)
for(i in 1:length(all.fun)){source(all.fun[i])}


all.fun <- list.files(pattern = ".R", path = useful.scripts.dir, full.names = TRUE)
if(length(all.fun) > 0){for(i in 1:length(all.fun)){source(all.fun[i])}}

# all.shiny.fun <- list.files(path = "~/Documents/R/shiny/capeDO_interactions/code", pattern = ".R", full.names = TRUE)
# for(i in 1:length(all.shiny.fun)){source(all.shiny.fun[i])}
#===============================================================
#specify parameters for the analysis
#some parameters are read in from exp_params.csv
#===============================================================

#===============================================================

#===============================================================
# read in the data
# This only needs to be done once
#===============================================================
if(build.cape.data){
	build.ssc.cape.data()
	}else{
	cross <- readRDS(paste0(init.dir, "cross.RData"))
	geno <- readRDS(paste0(init.dir, "geno.RData"))
	}

	if(testing){
		setwd(results.dir)
		chosen.snps <- read.table("fileredSNPstesting.txt", stringsAsFactors = FALSE)
		sampled.idx <- sort(match(chosen.snps[,1], cross$geno.names[[3]]))
		cross$geno.names[[3]] <- cross$geno.names[[3]][sampled.idx]
		cross$marker.num <- cross$marker.num[sampled.idx]
		cross$chromosome <- cross$chromosome[sampled.idx]
		cross$marker.location <- cross$marker.location[sampled.idx]
		parameter.file <- "cape.parameters.testing.txt"
		}
	
# add.maf <- colMeans(geno[,2,], na.rm = TRUE)	
	
if(geno.coding == "Dominant"){
	geno[which(geno >= 0.5)] <- 1
	}


# dom.maf <- colMeans(geno[,2,], na.rm = TRUE)
# plot(add.maf, dom.maf)


#===============================================================
#keep track of the experimental parameters
#===============================================================
if(!file.exists(results.dir)){system(paste("mkdir", results.dir))}
setwd(results.dir)

cross <- run.cape(cross, geno, parameter.file, p.or.q = 0.05, n.cores = n.cores)


# pairscan.obj <- readRDS("cross.pairscan.RData")
# plot.null.dist(cross, pairscan.obj)

		
# snp.motifs <- pheno.effects.DO(data.obj = cross, covar = cross$p.covar, scan.what = "normalized.traits", geno.coding = "Dominant")
		
		# #get genes asssociated with each motif snp
		# library(biomaRt)
		# hum = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
		# snp.db = useEnsembl(biomart="snp", dataset="hsapiens_snp")
		# motif.snps <- gsub("_B", "", unique(unlist(lapply(snp.motifs, function(x) unique(c(x[,1], x[,2]))))))

		# full.motif.table <- unique(Reduce("rbind", snp.motifs))
		# write.table(full.motif.table, "Motifs.All.txt", sep = '\t', quote = FALSE)

		# snp.table <- getBM(c("refsnp_id", "chr_name", "chrom_start", "minor_allele", "associated_gene"), filters = "snp_filter", values = motif.snps, mart = snp.db)
		# num.chr <- which(!is.na(as.numeric(snp.table[,2])))
		# snp.table <- snp.table[num.chr,]
		# snp.table <- condense.table(snp.table, 1, c(2,3,4), 5)
		# regions.for.viewing <- apply(snp.table, 1, function(x) paste0("Chr", x[2],":", as.numeric(x[3])-1000, "-", as.numeric(x[3])+1000))
		# cat(regions.for.viewing, sep = "\n")
		# cat(snp.table[,1], sep = "\n")

		# gene.table <- as.matrix(read.table("Motif.snp.gene.table.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, fill = TRUE))

		# full.snp.table <- matrix(NA, nrow = nrow(full.motif.table), ncol = 12)
		# for(i in 1:nrow(full.motif.table)){
			# snp1.locale <- which(gene.table[,1] == gsub("_B", "", full.motif.table[i,1]))	
			# snp1.info <- gene.table[snp1.locale,]

			# snp2.locale <- which(gene.table[,1] == gsub("_B", "", full.motif.table[i,2]))	
			# snp2.info <- gene.table[snp2.locale,]
			
			# full.snp.table[i,] <- c(full.motif.table[i,1], snp1.info[2:3], full.motif.table[i,2], snp2.info[2:3], full.motif.table[i,3:ncol(full.motif.table)])
			# }
			# colnames(full.snp.table) <- c("Source", "Pos", "Gene", "Target", "Pos", "Gene", colnames(full.motif.table)[3:ncol(full.motif.table)])
		# write.table(full.snp.table, "Motifs.All.txt", sep = '\t', quote = FALSE, row.names = FALSE)


		# # write.data.for.shiny(cross, "SSc.marker.selection", "~/Documents/R/shiny/capeDO_interactions/data")


		# # test.source <- "rs1344655_B"; test.target <- "rs660895_B"
		# # quartz();plot.effects(cross, test.source, test.target, "Normalized", plot.type = "b", error.bars = "se", num.rows = 2)

		# # quartz();plot.main.dist(cross, perm.stats, test.source, standardized = FALSE)
		# # quartz();max.marker.effect(cross, perm.stats, test.source, standardized = FALSE, covar = cross$p.covar)

		# # quartz();plot.main.dist(cross, perm.stats, test.target, standardized = FALSE)
		# # quartz();max.marker.effect(cross, perm.stats, test.target, standardized = FALSE, covar = cross$p.covar)

		# tissue.name = "skin"
		# # tissue.name = "skin_fibroblast"
		# library(gProfileR)
		# setwd("~/Documents/Data/Scleroderma/Results/Dominant_2ET_filtered_SNPs")
		# raven.fun <- list.files("~/Documents/git_repositories/raven", full.names = TRUE)
		# for(i in 1:length(raven.fun)){source(raven.fun[i])}
		# useful.fun <- list.files("~/Documents/git_repositories/useful_r_code", full.names = TRUE) 
		# for(i in 1:length(useful.fun)){source(useful.fun[i])}
		# skin.net <- readRDS(paste0("~/Documents/Data/FGN/human/", tissue.name, "_top.RData"))
		# library(biomaRt)
		# hum <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
		# # hum <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
		# snp.db = useEnsembl(biomart="snp", dataset="hsapiens_snp")
		# source('~/Documents/git_repositories/expression_plotting/plot.Hinchcliff.coexpression.R')
		# source('~/Documents/git_repositories/expression_plotting/plot.Hinchcliff.expression.R')
		
		# var.inf <- writeVariantInfluences(cross, 0.05, TRUE, write.file = FALSE)
		# int.locale <- grep("rs", var.inf[,"Target"])
		# just.int <- var.inf[int.locale,]
		# int.snps <- unique(gsub("_B", "", c(just.int[,1], just.int[,4])))

		# snp.table <- getBM(c("refsnp_id", "chr_name", "chrom_start", "minor_allele", "associated_gene"), filters = "snp_filter", values = int.snps, mart = snp.db)
		# num.chr <- which(!is.na(as.numeric(snp.table[,2])))
		# snp.table <- snp.table[num.chr,]
		# snp.info <- condense.table(snp.table, condense.by = 1, col.to.collapse = c(2,3,4), col.to.concat = 5)
		# write.table(snp.info, "Interaction.SNPs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
		
		# genes.near.snps <- get.nearest.gene2(snp.names = int.snps, "human", restrict.to.entrez = FALSE)
		# #the list of genes that are near the SNPs and whose expression is 
		# #influenced by these SNPs in the skin
		# coding.and.qtl <- read.table("coding.and.eQTL.txt", sep = "\t", stringsAsFactors = FALSE)

		# pdf("Gene.Expr.in.SSc.pdf")
		# plot.Hinchcliff.expression(unique(coding.and.qtl[,1]))
		# dev.off()
			
		# test.gene.info <- getBM(c("external_gene_name", "entrezgene"), "external_gene_name", coding.and.qtl[,1], hum)
		
		# test.gene.ids <- test.gene.info[which(!is.na(test.gene.info[,2])),]
		# gene.text <- "Int.SNPs"

		# snp.ints <- just.int[,c(1,4)]
		# snp.ints <- apply(snp.ints, 2, function(x) gsub("_B", "", x))
	
		# gene.ints <- snp.ints
		# gene.ints[,1] <- unlist(lapply(gene.ints[,1], function(x) genes.near.snps[which(genes.near.snps[,1] == x),"entrez_id"]))
		# gene.ints[,2] <- unlist(lapply(gene.ints[,2], function(x) genes.near.snps[which(genes.near.snps[,1] == x),"entrez_id"]))
		
		# #we also need to add to these interactions the interactions
		# #between eQTL genes, with each SNP representing its eQTL
		# #targets
		# get.gene.ids <- function(gene.names){
			# gene.ids <- rep(NA, length(gene.names))
				# for(i in 1:length(gene.names)){
					# gene.locale <- which(test.gene.ids[,1] == gene.names[i])[1]
					# if(length(gene.locale) > 0){
						# gene.ids[i] <- test.gene.ids[gene.locale,2]
						# }
					# }
				# return(gene.ids)
				# }

		
		# snp.eqtl <- read.table("snps.and.eqtl.txt", sep = "\t", stringsAsFactors = FALSE)
		# for(i in 1:nrow(snp.ints)){
			# snp1.locale <- which(snp.eqtl[,1] == snp.ints[i,1])
			# snp1.eqtl <- str_trim(unlist(strsplit(snp.eqtl[snp1.locale,5], ",")))
			# if(snp1.eqtl[1] == "none"){snp1.eqtl <- snp.eqtl[snp1.locale,3]}
			
			# snp2.locale <- which(snp.eqtl[,1] == snp.ints[i,2])
			# snp2.eqtl <- str_trim(unlist(strsplit(snp.eqtl[snp2.locale,5], ",")))
			# if(snp2.eqtl[1] == "none"){snp2.eqtl <- snp.eqtl[snp2.locale,3]}

			# gene.pairs <- cbind(rep(snp1.eqtl, length(snp2.eqtl)), rep(snp2.eqtl, each = length(snp1.eqtl)))
			# id1 <- get.gene.ids(gene.pairs[,1])
			# id2 <- get.gene.ids(gene.pairs[,2])
			# entrez.chunk <- cbind(id1, id2)
			# gene.ints <- rbind(gene.ints, entrez.chunk)
			# }
		# na.locale <- which(is.na(gene.ints), arr.ind = TRUE)
		# gene.ints <- gene.ints[-unique(na.locale[,1]),]


		# net.file <- paste0("GIANT.", gene.text, ".genes.RData")
		# if(!file.exists(net.file)){
			# # test.net <- get.network.context(net.genes = test.gene.ids[,2], tissue.net = skin.net, mart = hum, min.edge.weight = 0.15, add.genes.per.subnet = 20)
			# #get the network that is trimmed to match epistasis patterns
			# # test.net <- get.network.context2(net.genes = unique(c(gene.ints[,1], gene.ints[,2])), skin.net, hum, min.edge.weight = 0.15)
			# test.net <- get.network.context3(gene.edges = gene.ints, tissue.net = skin.net, mart = hum, start.thresh = 0.1, step.size = 0.01)
			# saveRDS(test.net, net.file)
			# }else{
			# test.net <- readRDS(net.file)
			# }

		
		# # pdf(paste0("GIANT.", tissue.name, ".net_",gene.text, ".pdf"))
		# # net.full <- plot.fntm.net(fntm.net = test.net, highlight.nodes = test.gene.ids[,2], min.edge.weight = min.edge.weight, vertex.size = 5, stretch.layout = 1, layout.x.adjust = 0, label.size = 0.5, layout.call = layout_with_kk)
		# # dev.off()
		
		# net.full <- test.net
		# vcount(net.full)
		
	# #find clusters in the network either with the fastgreedy network clustering algorithm
	# #or hierarchical clustering
	
	# get.mod <- function(gene.name){
		# gene.mod <- names(unlist(lapply(comm.genes, function(x) which(x == gene.name))))
		# return(gene.mod)
		# }
		
	# comm.gene.obj <- module.enrichments(fntm.net = test.net, num.terms = 30, set.names.manually = TRUE, module.fun = "fast_greedy", organism = "hsapiens", order.by = "gprofiler")
	# # comm.gene.obj <- module.enrichments(net.full, set.names.manually = TRUE, module.fun = "fast_greedy", organism = "hsapiens")
	# comm.genes <- comm.gene.obj[[1]]
	# gene.order <- comm.gene.obj[[2]]

	# gene.mods <- cbind(test.gene.ids[,1], rep(NA, nrow(test.gene.ids)))
	# for(i in 1:length(test.gene.ids[,1])){
		# gene.mod <- get.mod(test.gene.ids[i,1])
		# if(length(gene.mod) > 0){
			# gene.mods[i,2] <- gene.mod
			# }else{
			# gene.mods[i,2] <- "not present"	
			# }	
		# }

	
	# saveRDS(comm.genes, paste0("comm.genes.",gene.text,".RData"))
	
	# #look for previously reported SSc genes in the network
	# ssc.genes <- read.table("~/Documents/Data/Scleroderma/Previous_Results/Genopedia_SSc_related_genes.txt", stringsAsFactors = FALSE)
	# ssc.genes2 <- read.table("~/Documents/Data/Scleroderma/Previous_Results/SSc_published_SNP_associations.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
	# ssc.genes.names <- unlist(strsplit(ssc.genes2[,"MAPPED_GENE"], ", "))
	# all.ssc.genes <- unique(c(ssc.genes[,1], ssc.genes.names))
	# ssc.genes.in.net <- intersect(all.ssc.genes, V(net.full)$name)
	# ssc.gene.locale <- match(ssc.genes.in.net, V(net.full)$name)
	# SSc.genes <- rep("-", vcount(net.full))
	# SSc.genes[ssc.gene.locale]  <- "SSc"
	# SSc.genes <- as.data.frame(SSc.genes)
	# net.names <- V(net.full)$name
	# na.locale <- which(is.na(net.names))
	# net.names[na.locale] <- paste0("NA", 1:length(na.locale))
	# SSc.genes <- cbind(net.names, SSc.genes)
	# write.table(SSc.genes, paste0("SSc.genes.", gene.text, ".txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
	# length(which(SSc.genes[,2] == "SSc"))
	# nrow(SSc.genes)
	
	# ssc.genes.in.net <- as.vector(SSc.genes[which(SSc.genes[,2] == "SSc"),1])
	# ssc.mods <- unlist(lapply(ssc.genes.in.net, get.mod))
	# ssc.mod <- cbind(ssc.genes.in.net, ssc.mods)
	
	# write.table(ssc.mod, "SSc.Gene.Module.Assignment.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

	# full.adj <- as.matrix(as_adjacency_matrix(net.full, attr = c("weight")))
	# mod.locale <- lapply(comm.genes, function(x) match(x, rownames(full.adj)))

	# #plot the clustered adjacency matrix with labeled modules
	# mod.mem <- rep(names(comm.genes[1]), nrow(full.adj))
	# for(m in 2:length(mod.locale)){
		# mod.mem[mod.locale[[m]]] <- names(comm.genes)[m]
		# }
	# mod.mem <- data.frame(mod.mem)
	# rownames(mod.mem) <- net.names
	
	# write.table(mod.mem, "Network.Gene.Module.Assignments.txt", sep = "\t", quote = FALSE, col.names = FALSE)
	# # pdf(paste0("Module_Membership_", tissue.name, "_", gene.text, "_", min.edge.weight, ".pdf"), width = 13, height = 7)
	# # # pheatmap(full.adj[gene.order, gene.order], cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = mod.mem)
	# # pheatmap(full.adj[gene.order, gene.order], cluster_rows = FALSE, cluster_cols = FALSE)
	# # dev.off()

	# net.layout <- layout_nicely(net.full)
	# mods <- cluster_fast_greedy(net.full, weights = 1-E(net.full)$weight)$membership
	# plot(net.layout[,1], net.layout[,2] , col = mods)

	# pdf("Network2.pdf", width = 20, height = 20)
	# plot(net.full, vertex.color = mods, vertex.size = 5, vertex.label = NA)
	# dev.off()
	
	# shift.x <- 25; shift.y = 25; vertex.size = 2;mod.width = 15;mod.height = 15

	# pdf(paste0("Modular_Network_", tissue.name, "_", gene.text, ".grid.pdf"), width = mod.width, height = mod.height)
	# if(length(comm.genes) <= 3){
		# cluster.layout.matrix <- matrix(c(0,1,0,2,0,3), nrow = 2, byrow = TRUE)
		# }else{
			
		# cluster.layout.matrix <- layout_with_kk(erdos.renyi.game(length(comm.genes), p = 1))
		# # cluster.layout.matrix <- layout_on_sphere(erdos.renyi.game(length(comm.genes), p = 1))
		# # cluster.layout.matrix <- layout_in_circle(erdos.renyi.game(length(comm.genes), p = 1))
		# # cluster.layout.matrix <- layout_on_grid(erdos.renyi.game(length(comm.genes), p = 1))
		# # plot(cluster.layout.matrix[,1], cluster.layout.matrix[,2])
		# }
	# mod.net <- plot.modular.net(net = net.full, modules = as.vector(mod.mem[,1]), cluster.layout.matrix = cluster.layout.matrix, shiftx = shift.x, shifty = shift.y, vertex.size = vertex.size, vertex.names = NA, layout.fun = layout_on_grid)
	# dev.off()
	
	# #plot the same network with differential expression information
	# # diff.exp <- as.matrix(read.csv("~/Documents/Data/Scleroderma/data_files/Expression/Hinchcliff.Differential.Gene.expression.Probes.Averaged.csv", stringsAsFactors = FALSE))
	# diff.exp <- as.matrix(read.csv("~/Documents/Data/Scleroderma/data_files/Expression/Hinchcliff.Differential.Gene.expression.Individual.Probes.csv", stringsAsFactors = FALSE))

	# sig.val <- 0.01
	# sig.locale <- which(as.numeric(diff.exp[,"p.value"]) <= sig.val)
	
	# sig.genes <- diff.exp[sig.locale,1]

	# diff.exp.net <- intersect(V(mod.net)$name, sig.genes)
	# diff.exp.idx <- match(diff.exp.net, diff.exp[,1])
	# diff.exp.table <- diff.exp[diff.exp.idx,]

	# diff.val <- as.numeric(diff.exp.table[,"SSc.Mean"]) - as.numeric(diff.exp.table[,"Normal.Mean"])
	# all.diff.val <- rep(0, vcount(mod.net))
	# for(i in 1:nrow(diff.exp.table)){
		# v.locale <- which(V(mod.net)$name == diff.exp.table[i,1])
		# all.diff.val[v.locale] <- diff.val[i]
		# }
		
	# v.col <- colors.from.values(all.diff.val, split.at.vals = TRUE, split.points = 0.001, grad.dir = "ends", light.dark = "full", col.scale = c("blue", "brown"))

	# pdf(paste0("Modular_Network_Expression_", tissue.name, "_", gene.text, ".pdf"), width = mod.width, height = mod.height)
	# layout(matrix(c(1,2), ncol = 1), heights = c(1, 0.2))
	# par(mar = c(1,1,1,1))
	# plot(mod.net, vertex.color = v.col, edge.color = E(mod.net)$edge.color, vertex.size = vertex.size, vertex.label = NA)
	# par(mar = c(2,4,2,2))
	# imageWithTextColorbar(matrix(all.diff.val, ncol = 1), split.at.vals = TRUE, split.points = 0.001, col.scale = c("blue", "brown"), light.dark = "full", grad.dir = "ends", cex = 1, axis.line = 1, orientation = "h")
	# dev.off()
	
	# pdf(paste0("Modular_Network_Expression_", tissue.name, "_", gene.text, ".quantified.pdf"), width = 10, height = 5)
	# mod.expr <- lapply(comm.genes, function(x) all.diff.val[match(x, V(mod.net)$name)])
	# mod.expr <- lapply(mod.expr, function(x) x[which(x != 0)])
	# boxplot(mod.expr, las = 2);abline(h = 0)
	# plot.new()
	# plot.window(ylim = c(min(all.diff.val), max(all.diff.val)), xlim = c(1, length(comm.genes)))
	# par(bty = "n")
	# par(xpd = TRUE)
	# for(i in 1:length(mod.expr)){
		# if(length(mod.expr[[i]]) > 2){
			# vioplot(mod.expr[[i]], at = i, add = TRUE, col = "lightgray", box = FALSE)
			# }
		# }
	# stripchart(mod.expr, las = 2, vertical = TRUE, pch = 16, method = "jitter", col = rgb(31/256,120/256,180/256, alpha = 0.5), add = TRUE, cex = 0.5)
	# text(1:length(comm.genes), y = rep(min(all.diff.val), length(comm.genes)), labels = names(comm.genes), srt = 45, adj = 1, cex = 0.7)
	# par(xpd = FALSE)
	# abline(h =0);axis(2)
	# dev.off()


	# #plot the network highlighting genes previously associated with SSc
	# pdf("Modular.Network.SSc.genes2.pdf", width = mod.width, height = mod.height)
	# ssc.col <- rep("lightgray", vcount(mod.net))
	# ssc.locale <- match(ssc.genes.in.net, V(mod.net)$name)
	# ssc.col[ssc.locale] <- "red"
	
	# epi.locale <- match(as.vector(gene.ints), V(mod.net)$entrez.name)
	# ssc.col[epi.locale] <- "blue"
	
	# both.locale <- intersect(ssc.locale, epi.locale)
	# if(length(both.locale) > 0){
		# ssc.col[both.locale] <- "purple"
		# }
	# plot(mod.net, vertex.color = ssc.col, edge.color = E(mod.net)$edge.color, vertex.size = 3, vertex.label = NA)
	# dev.off()
	

	
	# plot.Hinchcliff.expression(test.genes[1])
	# plot.Hinchcliff.expression(test.genes[2])
	# plot.Hinchcliff.coexpression(test.genes)

	# # #plot a finer-scale enrichment for each module
	# # enrichment.sources <- c("GO", "KEGG","REAC","TF", "MI","CORUM", "HP", "HPA", "OMIM")
	# # pdf(paste0("GIANT.enrichment_", gene.text, "_", add.genes.per.subnet, ".pdf"))
	# # for(m in 1:length(comm.genes)){
			# # enrichment <- gprofiler(comm.genes[[m]], "hsapiens", max_p_value = 0.05)
			# # plot.enrichment(enrichment, 30, plot.label = paste("Module", m, "Overall"), text.size = 1)	
		# # for(e in 1:length(enrichment.sources)){
			# # print(enrichment.sources[e])
			# # enrichment <- gprofiler(comm.genes[[m]], "hsapiens", src_filter = enrichment.sources[e], max_p_value = 0.05)
			# # plot.enrichment(enrichment, 30, plot.label = paste("Module", m, enrichment.sources[e]), text.size = 1)
			# # }
		# # }
		# # dev.off()
	
		
	# expression.fun <- list.files(path = "~/Documents/git_repositories/expression_plotting", pattern = ".R", full.names = TRUE)
	# for(i in 1:length(expression.fun)){source(expression.fun[i])}

	# net.genes <- unlist(comm.genes)
	# net.genes <- net.genes[which(!is.na(net.genes))]

	# pdf(paste0("Expression.", gene.text, ".pdf")	)
	# for(i in 1:length(net.genes)){
		# cat(net.genes[i], "\n")
		# plot.Hinchcliff.expression(net.genes[i], average.multiple.probes = TRUE)
		# }		
	# dev.off()	
		
	# #look at gene expression between all pairs of genes in the network
	# #are there positive correlations within a module and negative between?

	# gene.pairs <- pair.matrix(net.genes)
	# pair.data <- vector(mode = "list", length = nrow(gene.pairs))
	# pdf(paste0("Expression.Net.Gene.Pairs.", gene.text, ".pdf"), width = 8, height = 4)
	# pair.data <- plot.Hinchcliff.coexpression(gene.pairs, rank.z.normalize = TRUE, average.multiple.probes = TRUE, color.by.grp = FALSE)
	# dev.off()
		
	# sig.val <- 0.05
	# adj.p <- p.adjust(as.numeric(pair.data[,"p"]), "fdr"); hist(adj.p, breaks = 100)
	# cor.net <- graph_from_edgelist(pair.data[,1:2], directed = FALSE)
	# edge.weights <- as.numeric(pair.data[,"r"])
	# edge.weights[which(adj.p > sig.val)] <- NA
	# edge.weights[which(is.na(edge.weights))] <- 0
	# E(cor.net)$weight <- edge.weights
	# cor.adj <- as_adjacency_matrix(cor.net, attr = "weight")
	# pheatmap(cor.adj, cluster_rows = FALSE, cluster_cols = FALSE)

	# cor.vals <- lapply(comm.genes, function(x) as.vector(cor.adj[x,x]))	
	# cross.vals <- as.vector(cor.adj[comm.genes[[1]], comm.genes[[2]]])
	# boxplot(cor.vals, main = "Correlation Values Between Genes")
	# abline(h = 0)
	# #correlation of expression between genes is positive in the GPSM3 module
	# #but all over the place in the FHIT module and between modules

	
	
	# #look at risk ratios for AA
	# snp.rr <- matrix(NA, ncol = 2, nrow = ncol(cross$pheno))
	# rownames(snp.rr) <- colnames(cross$pheno)
	# colnames(snp.rr) <- c("rs2301156", "rs204991")
	# for(ph in 1:ncol(cross$pheno)){
		# snp.rr[ph,] <- get.snp.risk.ratio(data.obj = cross, snps = c("rs2301156_B", "rs204991_B"), phenotype = colnames(cross$pheno)[ph], covar = "sex", trait.type = "norm", outcome.threshold = 0, geno.coding = "Dominant", plot.results = FALSE)
		# }
		
	# for(sn in 1:ncol(snp.rr)){
		# quartz()
		# barplot(snp.rr[,sn], ylim = c(0, max(1, snp.rr[,sn])), main = colnames(snp.rr)[sn])
		# abline(h = 1)
		# }
	
	# int.mat <- matrix(NA, ncol = ncol(cross$pheno), nrow = 10)
	# colnames(int.mat) <- colnames(cross$pheno)
	# for(ph in 1:ncol(cross$pheno)){
		# result <- get.snp.int.risk.ratio(cross, snp1 = "rs2301156_B", snp2 = "rs204991_B", allele.which = 2, phenotype = colnames(cross$pheno)[ph], covar = "sex", trait.type = "norm", outcome.threshold = 0, geno.coding = "Dominant", plot.results = FALSE)
		# int.mat[,ph] <- result
		# }
	# rownames(int.mat) <- rownames(result)
	
	# barplot(int.mat[3,], ylim = c(0, max(1, int.mat[3,])), main = )
	
	# pheno <- get.pheno(cross, "norm", "sex")
	# snp <- cross$geno.for.pairscan[,"rs204991_B"]
	# # snp <- cross$geno.for.pairscan[,"rs2301156_B"]
	# u_snp <- sort(unique(snp))
	# ph = 1
	# phenoV <- pheno[,ph]
	# phenoV[which(phenoV < 0)] <- 0
	# phenoV[which(phenoV > 0)] <- 1
	# table(phenoV, snp)

	# pheno.means <- lapply(u_snp, function(x) pheno[which(snp == x),ph])
	# boxplot(pheno.means)
	
	
	# total.aa <- apply(cross$pheno, 1, function(x) if(any(x == 1)){return(1)}else{return(0)})
	# risk.ratio(cross$geno.for.pairscan[,"rs2301156_B"], total.aa)
	# risk.ratio(cross$geno.for.pairscan[,"rs204991_B"], total.aa)
	# risk.interactions(cross$geno.for.pairscan[,"rs2301156_B"], cross$geno.for.pairscan[,"rs204991_B"], total.aa, plot.results = TRUE)


	# risk.ratio(cross$geno.for.pairscan[,"rs2301156_B"], cross$pheno[,"centromere"])
	# risk.ratio(cross$geno.for.pairscan[,"rs204991_B"], cross$pheno[,"centromere"])
	# risk.interactions(cross$geno.for.pairscan[,"rs2301156_B"], cross$geno.for.pairscan[,"rs204991_B"], cross$pheno[,"centromere"], plot.results = TRUE)
	
	# pheno.cor <- cor(cross$pheno)
	# diag(pheno.cor) <- NA
	# imageWithText(pheno.cor, cex = 1, split.at.vals = TRUE, col.scale = c("blue", "brown"), light.dark = "l", col.names = colnames(cross$pheno), row.name = colnames(cross$pheno), grad.dir = "ends")


	# #look for papers written about the community genes and SSc
	# comm.genes <- readRDS("comm.genes.RData")
	# library(easyPubMed)
	# all.genes <- c(unlist(comm.genes), "Wnt")
	# # all.pubmed <- lapply_pb(as.vector(all.genes), function(x) get_pubmed_ids(paste("Scleroderma AND", x)))
	# all.pubmed <- lapply_pb(as.vector(all.genes), function(x) get_pubmed_ids(paste("Systemic sclerosis AND", x)))

	# num.papers <- unlist(lapply(all.pubmed, function(x) x$Count))
	# with.papers <- which(num.papers != 0)
	# cbind(all.genes[with.papers], num.papers[with.papers])
	
	# has.data <- all.pubmed[with.papers]
	# all.abstracts <- lapply(has.data, function(x) fetch_pubmed_data(x, format = "abstract"))
	# for(i in 1:length(all.abstracts)){
		# write.table(all.abstracts[[i]], paste0("PubMed_Query_", all.genes[with.papers[i]], ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
		# }
	
	# }
	

