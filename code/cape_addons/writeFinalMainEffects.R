writeFinalMainEffects <- function(data.obj){
	
	D1.results <- data.obj$max.var.to.pheno.influence
	#order the results by marker order
	D1.results <- lapply(D1.results, function(x) x[order(x[,1]),])
	marker.names <- data.obj$marker.names
	ind.markers <- data.obj$geno.for.pairscan
		
	if(is.null(D1.results)){
		stop("direct.influence() must be run before plotting the results")
		}

	chr <- unique(data.obj$chromosome[which(colnames(data.obj$geno) %in% D1.results[[1]][,1])])
	traits <- names(D1.results)
	
	
	chr.locale <- which(data.obj$chromosome %in% chr)
	markers.used <- which(colnames(data.obj$geno) %in% D1.results[[1]][,1])
	markers.which <- intersect(chr.locale, markers.used)
	results.rows <- which(D1.results[[1]][,1] %in% colnames(data.obj$geno)[markers.which])	
	
	
	results.el <- which(names(D1.results) %in% traits)
	results.to.plot <- NULL
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(D1.results[[1]]))))
	tval.col <- which(colnames(D1.results[[1]]) == "|t.stat|")
	coef.col <- which(colnames(D1.results[[1]]) == "coef")
	for(r in results.el){
		if(standardized){
			results.table <- D1.results[[r]][results.rows,c(tval.col, var.sig.col)]
			colnames(results.table) <- paste(traits[r], colnames(results.table), sep = "_")
			results.to.plot <- cbind(results.to.plot, results.table)
			}else{
			results.table <- D1.results[[r]][results.rows,c(coef.col, var.sig.col)]
			colnames(results.table) <- paste(traits[r], colnames(results.table), sep = "_")
			results.to.plot <- cbind(results.to.plot, results.table)
			}	
		}

	final.table <- cbind(data.obj$marker.names[markers.which], data.obj$chromosome[markers.which], data.obj$marker.location[markers.which], results.to.plot)
	colnames(final.table)[1:3] <- c("marker", "chromosome", "location")
	write.table(final.table, "Variant.Influences.Main.csv", sep = ",", quote = FALSE, row.names = FALSE)

	
}