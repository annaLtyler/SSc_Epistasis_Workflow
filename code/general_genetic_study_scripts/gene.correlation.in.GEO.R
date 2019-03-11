#This function takes in an expression data set that has been processed by
#GEO2R, as in Download_SSc_lung_dataset.R
#It expects that the group labels are in the column names, separated from the
#sample names by "_"
#The gene IDs are in rownames. The gene id table is also created by the 
#GEO2R data set. It is the top table printed for all genes. This table is used
#to get gene probe IDs for genes
#For now, if this function finds multiple probes for a gene, they are averaged
#together.

gene.correlation.in.GEO <- function(gene.pairs, gene.id.table, expr, adjust.by.group = TRUE){

	require(RColorBrewer)
	group.colors <- brewer.pal(8, "Set1")
	
	groups <- sapply(strsplit(colnames(expr), "_"), function(x) x[2])
	u_groups <- unique(groups)
	group.cols <- groups
	for(i in 1:length(u_groups)){
		group.cols[which(groups == u_groups[i])] <- group.colors[i]
		}
	
	adjust.na <- function(out.var, adj.var){
		na.locale <- unique(which(is.na(cbind(out.var, adj.var)), arr.ind = TRUE)[,1])
		not.na.locale <- setdiff(1:length(out.var), na.locale)
		resid.var <- rep(NA, length(out.var))
		resids <- residuals(lm(out.var~adj.var))
		resid.var[not.na.locale] <- resids
		names(resid.var) <- names(out.var)
		return(resid.var)
		}

	all.expr <- list()
	adj.expr <- list()
	all.gene.pairs <- NULL
	all.num.probes <- NULL
	all.p <- NULL
	all.r <- NULL
	with.expr.count <- 1
	symbol.col <- grep("symbol", colnames(gene.id.table), ignore.case = TRUE)
		
	for(g in 1:nrow(gene.pairs)){	
		gene1.locale <- which(gene.id.table[,symbol.col] == gene.pairs[g,1])
		gene2.locale <- which(gene.id.table[,symbol.col] == gene.pairs[g,2])
		if(length(gene1.locale) == 0 || length(gene2.locale) == 0){
			next()
			}else{

			gene1.id <- gene.id.table[gene1.locale,1]
			gene2.id <- gene.id.table[gene2.locale,1]
			
			id1.locale <- match(gene1.id, rownames(expr))
			id2.locale <- match(gene2.id, rownames(expr))			
			
			gene1.expr <- expr[id1.locale,,drop=FALSE]
			gene2.expr <- expr[id2.locale,,drop=FALSE]			
			
			all.gene.pairs <- rbind(all.gene.pairs, gene.pairs[g,])
			num.probes <- c(nrow(gene1.expr), nrow(gene2.expr))
			all.num.probes <- rbind(all.num.probes, num.probes)
			
			#average across all probes
			gene1.expr <- colMeans(gene1.expr)
			gene2.expr <- colMeans(gene2.expr)
					
			#mean center expression for each gene
			gene1.expr <- gene1.expr - mean(gene1.expr, na.rm = TRUE)
			gene2.expr <- gene2.expr - mean(gene2.expr, na.rm = TRUE)
			
			all.expr[[with.expr.count]] <- rbind(gene1.expr, gene2.expr)		
			
			#adjust for group status
			if(adjust.by.group){
			adj.expr <- apply(all.expr[[with.expr.count]], 1, function(x) adjust.na(x, as.factor(groups)))
			all.expr[[with.expr.count]] <- adj.expr
			}else{
			all.expr[[with.expr.count]] <- t(all.expr[[with.expr.count]])	
			}

				
			test <- cor.test(all.expr[[with.expr.count]][,1], all.expr[[with.expr.count]][,2])
			r <- signif(test$estimate, 2)
			pval <- signif(test$p.value, 2)
			all.p <- c(all.p, pval)
			all.r <- c(all.r, r)

			with.expr.count = with.expr.count + 1
			}
		}
		

		layout.mat <- get.layout.mat((length(all.expr)+1))
		layout(layout.mat)
		for(g in 1:length(all.expr)){
			plot(all.expr[[g]][,1], all.expr[[g]][,2], col = group.cols, pch = 16, xlab = paste(all.gene.pairs[g,1], "Expression"), ylab = paste(all.gene.pairs[g,2], "Expression"), main = paste0("r = ", all.r[g], ": p = ", all.p[g]))
			model <- lm(all.expr[[g]][,1]~all.expr[[g]][,2])
			abline(model)
			}
		par(xpd = TRUE)
		plot.new()
		plot.window(xlim = c(0,1), ylim = c(0,1))
		legend(0,0.9, legend = u_groups, col = group.colors[1:length(u_groups)], pch = 16, cex = 2)
		par(xpd = FALSE)

}