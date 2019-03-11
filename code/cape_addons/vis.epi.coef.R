#visualize the epistasis coefficients between a marker and a chromosome

vis.epi.coef <- function(data.obj, pairscan.obj, chr = 6, marker = "sex", pval.correction = c("none","holm", "fdr")){

chr.marker.locale <- which(data.obj$chromosome == chr)
chr.markers <- colnames(data.obj$geno)[chr.marker.locale]
sex.locale <- which(data.obj$marker.names == marker)
sex <- colnames(data.obj$geno)[sex.locale]

	if(missing(pairscan.obj)){
	pair.scan <- data.obj$pairscan.results	
	}else{
	pair.scan <- pairscan.obj$pairscan.results
	}

#find epistasis coefficients for chr 6 markers and sex for each phenotype
coef.table <- NULL
pval.table <- NULL
for(i in 1:length(pair.scan)){
		last.col <- dim(pair.scan[[i]][[1]])[2]
		m1 <- pair.scan[[i]][[1]][,1]
		m2 <- pair.scan[[i]][[1]][,2]
		
		chr.locale <- c(which(m1 %in% as.numeric(chr.markers)), which(m2 %in% as.numeric(chr.markers)))
		sex.locale <- c(which(m1 %in% sex), which(m2 %in% sex))
		
		sex.chr.pair <- intersect(chr.locale, sex.locale)
		epi.coeff <- pair.scan[[i]][[1]][sex.chr.pair,last.col]/pair.scan[[i]][[2]][sex.chr.pair,last.col]
		coef.table <- cbind(coef.table, epi.coeff)
		
		all.pvals <- pair.scan[[i]]$model.pvalues[,3]
		if(pval.correction == "fdr"){
			cpvals <- fdrtool(all.pvals, statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = "fndr")$qval
			}
		if(pval.correction == "holm"){
			cpvals <- p.adjust(all.pvals, "holm")
			}
		if(pval.correction == "none"){
			cpvals <- all.pvals
			}
		pvals <- cpvals[sex.chr.pair]
		pval.table <- cbind(pval.table, pvals)
		}
		
	pairs <- matrix(pair.scan[[1]][[1]][sex.chr.pair,1:2], ncol = 2)
	pair.names <- apply(pairs, 2, function(x) data.obj$marker.names[match(x, colnames(data.obj$geno))])
	
	colnames(coef.table) <- colnames(pval.table) <- names(pair.scan)
	rownames(coef.table) <- rownames(pval.table) <- apply(pair.names, 1, function(x) paste(x, collapse = "-"))
	pdf(paste("Epistasis.Coefficients.for.Chr", ".and.", marker, ".", pval.correction, ".pdf", sep = ""), width = 7, height = 12)
	par(mar = c(4,7,2,2), mfrow = c(2,1))
	imageWithText(coef.table, col.names = colnames(coef.table), col.text.rotation = 90, col.text.adj = 1, row.text.adj = 0, row.text.shift = 1.2, split.at.vals = TRUE, row.names = rownames(coef.table), main = "Standardized Epistasis Coefficients", grad.dir = "ends")

	imageWithText(pval.table, col.names = colnames(coef.table), col.text.rotation = 90, col.text.adj = 1, row.text.adj = 0, row.text.shift = 1.2, split.at.vals = FALSE, row.names = rownames(coef.table), main = "P Values", grad.dir = "low")
	dev.off()
	
	invisible(list("coef.table" = coef.table, "pval.table" = pval.table))
	}