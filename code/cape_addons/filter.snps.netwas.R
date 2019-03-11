#This function take in a data.obj, a geno.obj
#the snp-gene table and NetWAS results and 
#prioritizes SNPs in the data object for running
#in the pairscan


filter.snps.netwas <- function(data.obj, geno.obj = NULL, snp.info.file = "SNP.Gene.Info.txt", sig.val, tissue.type, top.n.genes = NULL, select.by = c("mean", "median", "max"), verbose = TRUE, plot.diagnostics = FALSE){
	
	if(length(select.by) > 1){
		select.by <- "mean"
		}
	
	snp.info <- read.table(snp.info.file, stringsAsFactors = FALSE, sep = "\t")
	
	netwas.filenames <- get.files(want = c("NetWAS", sig.val, tissue.type), dont.want = "netwas", ignore.case = FALSE)
	
	netwas.genes <- NULL
	if(verbose){cat("Reading in NetWAS results...\n")}
	for(n in 1:length(netwas.filenames)){
		results <- read.table(netwas.filenames[n], stringsAsFactors = FALSE, header = TRUE)
		#take only genes in the dataset that are above the hyperplane
		gene.locale <- intersect(which(results[,2] != 0), which(results[,3] > 0)) 
		netwas.genes <- rbind(netwas.genes, results[gene.locale,])
		}
	
	get.all.vals <- function(gene){
		gene.locale <- which(netwas.genes[,1] == gene)
		gene.table <- netwas.genes[gene.locale,]
		return(gene.table)
		}
	

	u.genes <- unique(netwas.genes[,1])
	all.gene.vals <- lapply(u.genes, get.all.vals)

	gene.stat <- unlist(lapply(all.gene.vals, function(x) eval(call(select.by, x[,3]))))
	stat.order <- order(gene.stat, decreasing = TRUE)

	num.sig.genes <- length(u.genes)

	if(is.null(top.n.genes)){
		num.sel.genes <- num.sig.genes
		}else{
		num.sel.genes <- top.n.genes
		}

	gene.col <- rep("black", num.sig.genes)
	gene.col[1:num.sel.genes] <- "red"

	if(verbose){cat("There are", num.sig.genes, "genes above the hyperplane for all traits.\n")}
	if(verbose){cat("Selecting the top", num.sel.genes, "genes.\n")}
	
	selected.genes <- u.genes[stat.order[1:num.sel.genes]]


	if(plot.diagnostics){
		if(verbose){cat("Plotting Gene Statistics...\n")}
		pdf(paste("Gene.vals.ordered.by.", select.by, ".pdf", sep = ""), width = 20, height = 5)
		plot.new()
		plot.window(xlim = c(1,length(all.gene.vals)), ylim = c(0, 1))
		na <- lapply_pb(1:length(gene.stat), function(x) boxplot(all.gene.vals[[stat.order[x]]][,3], add = TRUE, at = x, axes = FALSE, col = gene.col[x]))
		axis(1); axis(2)
		dev.off()
		}
		
	
			
	if(plot.diagnostics){
		pdf(paste("Selected.genes.", select.by, ".distance.from.hyperplane.pdf", sep = ""))
		par(mar = c(5, 5, 4, 2))
		hist(gene.stat, xlab = paste(select.by, "distance from hyperplane"), ylab = "number of genes", main = "All Positive Genes from NetWAS", cex.lab = 2, cex.axis = 1.5)
		hist(gene.stat[stat.order[1:num.sel.genes]], add = TRUE, col = "red")
		legend("topright", fill = "red", legend = "Selected", cex = 2)
		dev.off()
		}

	
	
	#find the snps that are in these genes
	gene.snp.locale <- which(snp.info[,2] %in% selected.genes)
	sel.snps <- snp.info[gene.snp.locale,]
	
	#find these snps in the data.obj
	snp.locale <- match(sel.snps[,1], data.obj$marker.names)
	snp.locale <- sort(snp.locale[which(!is.na(snp.locale))])
	num.snps <- length(snp.locale)
	
	if(verbose){cat(num.snps, "have been selected for the pairscan based on NetWAS results.\n")}
	
	
	if(verbose){cat("Selecting markers and calculating linear independence...\n")}
	geno <- get.geno(data.obj, geno.obj)
	
	pair.geno <- geno[,snp.locale]
	geno.ind <- get.linearly.independent(pair.geno)
	rejected.markers <- geno.ind[[2]]
	geno.for.pairscan <- geno.ind[[1]]
	
		
	if(length(rejected.markers) > 0){
		message("\n", length(rejected.markers), " marker(s) rejected due to linear non-independence.\n For more information see markers.removed.for.non.independence.txt")
		write.table(colnames(pair.geno)[sort(rejected.markers)], "markers.removed.for.non.independence.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
		}

	data.obj$geno.for.pairscan <- geno.for.pairscan
	
	return(data.obj)

}