#This function takes in an expression data set that has been processed by
#GEO2R, as in Download_SSc_lung_dataset.R
#It expects that the group labels are in the column names, separated from the
#sample names by "_"
#The gene IDs are in rownames. The gene id table is also created by the 
#GEO2R data set. It is the top table printed for all genes. This table is used
#to get gene probe IDs for genes
#For now, if this function finds multiple probes for a gene, they are averaged
#together.

DEG.in.GEO <- function(gene.names = c("WNT5A", "RBMS3"), gene.id.table, expr, plot.results = TRUE, las = 1, avg.probes = TRUE, return.data = c("means", "all")){

	groups <- sapply(strsplit(colnames(expr), "_"), function(x) x[2])
	u_groups <- unique(groups)

	all.expr <- vector(mode = "list", length = length(gene.names))
	names(all.expr) <- gene.names
	all.num.probes <- rep(NA, length(gene.names))
	all.p <- vector(mode = "list", length = length(gene.names))
	all.means <- vector(mode = "list", length = length(gene.names))
	names(all.means) <- gene.names
	symbol.col <- grep("symbol", colnames(gene.id.table), ignore.case = TRUE)
	genes.with.expr <- 0
	
	for(g in 1:length(gene.names)){	
		gene.locale <- which(gene.id.table[,symbol.col] == gene.names[g])
		if(length(gene.locale) == 0){
			next()
			}else{

			gene.id <- gene.id.table[gene.locale,1]
			id.locale <- match(gene.id, rownames(expr))
			gene.expr <- expr[id.locale,,drop=FALSE]
			num.probes <- nrow(gene.expr)
			all.num.probes[g] <- num.probes
			
			if(avg.probes){
				avg.expr <- matrix(colMeans(gene.expr), nrow = 1)
				colnames(avg.expr) <- colnames(gene.expr)
				gene.expr <- avg.expr
				}
		
			test <- apply(gene.expr, 1, function(x) anova(aov(as.numeric(x)~as.factor(groups))))
			# test <- anova(aov(as.numeric(gene.expr)~as.factor(groups)))
			pval <- sapply(test, function(x) signif(x$"Pr(>F)"[1], 2))
			all.p[[g]] <- pval
			
			group.vals <- apply(gene.expr, 1, function(y) lapply(u_groups, function(x) as.numeric(y[which(groups == x)])))
			for(gr in 1:length(group.vals)){
				names(group.vals[[gr]]) <- u_groups
				}
			group.means <- t(sapply(group.vals, function(x) sapply(x, mean)))
			
			# group.vals <- lapply(u_groups, function(x) as.numeric(gene.expr[which(groups == x)]))
			# names(group.vals) <- u_groups
			# group.means <- sapply(group.vals, mean)
			all.means[[g]] <- group.means
			all.expr[[g]] <- group.vals
			genes.with.expr <- genes.with.expr + 1
			}
		}


		if(plot.results){
			layout.mat <- get.layout.mat(min(genes.with.expr, 9))
			layout(layout.mat)
			for(g in 1:length(all.expr)){
				if(!is.null(all.expr[[g]])){
				
				group.mats <- vector(mode = "list", length = length(u_groups))
				names(group.mats) <- u_groups
				for(gr in 1:length(u_groups)){
					group.mats[[gr]] <- sapply(all.expr[[g]], function(x) x[[gr]])	
					}

				plot.grouped.boxes(group.mats, group.labels = u_groups, type = "matrix", main = paste(gene.names[g], "\np =", min(all.p[[g]])))
				# stripchart(all.expr[[g]], vertical = TRUE, add = TRUE, pch = 16, col = "darkgray", method = "jitter")
				}
				}
			}

		if(return.data[1] == "all"){
			invisible(all.expr)
			}else{

			final.results <- sapply(1:length(all.means), function(x) cbind(all.means[[x]], all.p[[x]]))
			
			if(class(final.results) == "matrix"){
				colnames(final.results) <- c(u_groups, "pval")
				rownames(final.results) <- gene.names
				}else{
				not.null <- which(!sapply(final.results, is.null))
				final.results <- final.results[not.null]
				for(i in 1:length(final.results)){
					colnames(final.results[[i]])[3]  <- "pval"
					}
				names(final.results) <- gene.names[not.null]
				}
			
			invisible(final.results)
			}

}