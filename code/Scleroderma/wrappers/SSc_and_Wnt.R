#SSc genes and Wnt signaling
#some objects in this script are set in 
#SSc_cape

	min.edge.weight <- 0.13

	setwd("~/Documents/Data/Scleroderma/Results/Dominant_max.pair_2ET")
	source('~/Documents/git_repositories/useful_r_code/pair.matrix.R')
	source('~/Documents/git_repositories/useful_r_code/draw.rectangle.R')
	source('~/Documents/git_repositories/useful_r_code/report.progress.R')
	expr.fun <- list.files("~/Documents/git_repositories/expression_plotting", pattern = ".R", full.names = TRUE)
	for(i in 1:length(expr.fun)){source(expr.fun[i])}
	
	library(KEGGREST)
	library(igraph)
	library(pheatmap)
	
	comm.genes <- paste0("comm.genes.", min.edge.weight, ".RData")
	#look for overlap with Wnt signaling
	path.obj <- keggGet("hsa04310")
	wnt.gene.list <- path.obj[[1]]$GENE
	just.genes <- wnt.gene.list[seq(2, length(wnt.gene.list), 2)]
	wnt.genes <- unlist(lapply(strsplit(just.genes, ";"), function(x) x[1]))
	
	wnt.genes[grep("CSNK1E", wnt.genes)] <- "CSNK1E"

	pdf("Wnt.Gene.Expr.in.SSc.pdf")
	all.wnt.expr <- plot.Hinchcliff.expression(wnt.genes, average.multiple.probes = TRUE)
	dev.off()
	
	gene.symbol <- hinchcliff.probe.table[,7]
	col.genes <- gene.symbol[grep("COL", gene.symbol)]

	pdf("Col.Gene.Expr.in.SSc.pdf")
	all.col.expr <- plot.Hinchcliff.expression(col.genes, average.multiple.probes = TRUE)
	dev.off()
	
	tgf.genes <- unique(gene.symbol[grep("TGFB", gene.symbol)])
	pdf("TGF.Gene.Expr.in.SSc.pdf")
	all.tgf.expr <- plot.Hinchcliff.expression(tgf.genes, average.multiple.probes = TRUE)
	dev.off()
	
	#===================================================
	#internal functions
	#===================================================
	
	plot.sig.diff <- function(all.expr.table, sig.val, label = "WNT"){
		ssc.n.med <- apply(all.expr.table[,9:10], 2, as.numeric)
		rownames(ssc.n.med) <- all.expr.table[,1]
		ssc.n.p <- as.numeric(all.expr.table[,11])
		adj.p <- p.adjust(ssc.n.p, "none")
		sig.locale <- which(adj.p < sig.val)
		sig.diff <- ssc.n.med[sig.locale,]
		ssc.diff <- sig.diff[,2] - sig.diff[,1]
	
		pdf(paste0(label, ".Gene.Expr.in.SSc.p.", sig.val, ".pdf"), width = 12, height = 5)
		barplot(sort(ssc.diff), las = 2)
		barplot(ssc.diff[order(names(ssc.diff))], las = 2)
		dev.off()
		}
		
	plot.heatmap <- function(pair.data, sig.val, cex = 1){
		sig.locale <- which(as.numeric(pair.data[,"adj.p"]) <= sig.val)	
		cor.net <- graph_from_edgelist(pair.data[sig.locale,1:2], directed = FALSE)
		edge.weights <- as.numeric(pair.data[sig.locale,"r"])
		E(cor.net)$weight <- edge.weights
		cor.adj <- as_adjacency_matrix(cor.net, attr = "weight")
		pheatmap(cor.adj, fontsize = cex, border_color = NA)
		return(list(cor.net, cor.adj))
		}
	
	get.pairwise.cor <- function(gene.list, label = "WNT", fresh.start = FALSE){
		data.file <- paste0(label, ".Gene.Pairwise.Correlations.RData")
		if(fresh.start){unlink(data.file)}
		if(!file.exists(data.file)){
			gene.pairs <- pair.matrix(gene.list)
			pair.data <- vector(mode = "list", length = nrow(gene.pairs))
			pdf(paste0(label, "Expression.Net.Gene.Pairs.pdf"), width = 8, height = 4)
			pair.data <- plot.Hinchcliff.coexpression(gene.pairs, rank.z.normalize = TRUE, average.multiple.probes = TRUE, color.by.grp = FALSE, verbose = TRUE)
			dev.off()
			not.na.locale <- which(!is.na(pair.data[,"p"]))
			pair.data <- pair.data[not.na.locale,]	
			adj.p <- p.adjust(as.numeric(pair.data[,"p"]), "fdr")
			pair.data <- cbind(pair.data, adj.p)
			saveRDS(pair.data, data.file)	
			}else{
			pair.data <- readRDS(data.file)
			}
		return(pair.data)
		}
		
	plot.gene.cor.to.list <- function(genes.to.show = c("FHIT", "GPSM3"), cor.adj, label = "WNT"){
	
		gene.locale <- match(genes.to.show, colnames(cor.adj))
		stripped.mat <- as.matrix(cor.adj[gene.locale,,drop=FALSE])
		
		remove.locale <- match(genes.to.show, colnames(stripped.mat))
		stripped.mat <- stripped.mat[,-remove.locale]
				
		pdf(paste0(label, ".Individual.Correlations.pdf"), width = 30, height = 7)	
		for(i in 1:nrow(stripped.mat)){
			barplot(sort(stripped.mat[i,]), main = rownames(stripped.mat)[i], las = 2)
			}
		dev.off()
	
		gene.to.all.cor <- as.matrix(cor.adj[,gene.locale])
		not.both.zero <- which(apply(gene.to.all.cor, 1, function(x) !all(x == 0)))
		just.sig <- gene.to.all.cor[not.both.zero,]
		gene.locale <- match(genes.to.show, rownames(just.sig))
		just.sig <- just.sig[-gene.locale,,drop=FALSE]
					
		pdf(paste0(label, ".Correlation.to.GPSM3.FHIT.p.", sig.val, ".pdf"), height = 10)
		pheatmap(just.sig)
		dev.off()	
		}
	#===================================================

	plot.sig.diff(all.wnt.expr, 0.05, "WNT")
	plot.sig.diff(all.col.expr, 0.05, "COL")
	plot.sig.diff(all.tgf.expr, 0.05, "TGF")
	
	lapply(comm.genes, function(x) intersect(x, wnt.genes))
	
	#look at correlations between module genes and WNT genes
	net.wnt.genes <- c("GPSM3", "FHIT", wnt.genes)
	net.col.genes <- c("GPSM3", "FHIT", col.genes)
	net.tgf.genes <- c("GPSM3", "FHIT", tgf.genes)

	wnt.pair.data <- get.pairwise.cor(gene.list = net.wnt.genes, "WNT", fresh.start = TRUE)
	col.pair.data <- get.pairwise.cor(net.col.genes, "COL", fresh.start = TRUE)
	tgf.pair.data <- get.pairwise.cor(net.tgf.genes, "TGF", fresh.start = TRUE)	

	sig.val <- 0.05

	pdf(paste0("WNT.Gene.Net.p.", sig.val, ".pdf"))
	cor.obj <- plot.heatmap(wnt.pair.data, sig.val, 4)
	wnt.cor.adj <- cor.obj[[2]]	
	dev.off()

	pdf(paste0("COL.Gene.Net.p.", sig.val, ".pdf"))
	cor.obj <- plot.heatmap(col.pair.data, sig.val, 4)
	col.cor.adj <- cor.obj[[2]]	
	dev.off()

	pdf(paste0("TGF.Gene.Net.p.", sig.val, ".pdf"))
	cor.obj <- plot.heatmap(tgf.pair.data, sig.val, 4)
	tgf.cor.adj <- cor.obj[[2]]	
	dev.off()


	plot.gene.cor.to.list(genes.to.show = c("GPSM3", "FHIT"), cor.adj = as.matrix(wnt.cor.adj), "WNT")
	plot.gene.cor.to.list(genes.to.show = c("GPSM3", "FHIT"), as.matrix(col.cor.adj), "COL")	
	plot.gene.cor.to.list(genes.to.show = c("GPSM3", "FHIT"), as.matrix(tgf.cor.adj), "TGF")	


	plot.Hinchcliff.coexpression(c("GPSM3", "ACTA2"))
	plot.Hinchcliff.coexpression(c("FHIT", "ACTA2"))