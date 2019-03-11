
plot.pval.dist <- function(data.obj, perm.obj){

	emp.p.table <- perm.obj$var.to.pheno.test.stat
	markers <- unique(emp.p.table[[1]][,1])
	marker.names <- data.obj$marker.names[match(markers, colnames(data.obj$geno))]


	pdf("marker.pvals.pdf")
	par(mfrow = c(3,3))
	
	for(i in unique(emp.p.table[[1]][,1])){
		for(j in 1:length(emp.p.table)){
			marker.p <- emp.p.table[[j]][which(emp.p.table[[j]][,1] == i),"emp.p"]
			fdr.out <- suppressWarnings(fdrtool(marker.p, statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = "pct0"))
			adj.p <- fdr.out$qval
			hist(marker.p, xlim = c(0,1), main = paste(marker.names[i], names(emp.p.table)[j], "\nadj.p =", signif(min(adj.p), 2)), breaks = seq(0,1,0.01))
			}
		}
		dev.off()


	# pdf("marker.pvals.box.pdf")
	# par(mfrow = c(3,3))
	
	# for(i in unique(emp.p.table[[1]][,1])){
		# for(j in 1:length(emp.p.table)){
			# marker.p <- emp.p.table[[j]][which(emp.p.table[[j]][,1] == i),"emp.p"]
			# boxplot(marker.p, ylim = c(0,1), main = paste("marker", i, names(emp.p.table)[j], "\n", signif(median(marker.p), 2)))
			# }
		# }
		# dev.off()
	
	
	
	pdf("marker.effects.pdf")
	par(mfrow = c(3,3))
	
	for(i in unique(emp.p.table[[1]][,1])){
		for(j in 1:length(emp.p.table)){
			marker.t <- emp.p.table[[j]][which(emp.p.table[[j]][,1] == i),"|t.stat|"]
			hist(marker.t, main = paste("marker", i, names(emp.p.table)[j]), xlim = c(0,6))
			}
		}
		dev.off()
	
}