#This function writes out the significant influences
#as an edge list

writeInteractions <- function(data.obj, source.block = NULL, target.block = NULL, p.or.q = 0.05, filename = "Variant.Interactions.csv", delim = ",", mark.covar = FALSE, write.file = TRUE){
	
	var.influences <- data.obj$var.to.var.p.val
	
	if(is.null(var.influences)){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}
	

	get.block <- function(marker.name){
		return(names(unlist(lapply(data.obj$linkage.blocks.collapsed, function(x) which(x == marker.name)))))	
		}

	source.locale <- NULL
	target.locale <- NULL
	subset.locale <- 1:dim(var.influences)[1]
	if(!is.null(source.block)){
		cat("Finding source linkage blocks...\n")
		source.blocks <- t(apply(matrix(var.influences[,1], ncol = 1), 1, get.block))
		source.locale <- which(source.blocks %in% source.block)
		}

	if(!is.null(target.block)){
		cat("Finding target linkage blocks...\n")
		target.blocks <- t(apply(matrix(var.influences[,2], ncol = 1), 1, get.block))
		target.locale <- which(target.blocks %in% target.block)
		}

	if(!is.null(source.locale) && !is.null(target.locale)){
		subset.locale <- intersect(source.locale, target.locale)
		}
	if(!is.null(source.locale) || !is.null(target.locale)){
		subset.locale <- c(source.locale, target.locale)	
		}
	
	if(length(subset.locale) > 0){
		var.influences <- var.influences[subset.locale,,drop=FALSE]
		}
		
	if(data.obj$transform.to.phenospace){
		pheno.names <- colnames(data.obj$pheno)
		}else{
			pheno.names <- names(data.obj$pairscan.results)
			}
	num.pheno <- length(pheno.names)
	
	
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences))))
	sig.var <- which(as.numeric(var.influences[, var.sig.col]) <= p.or.q)
	
	
	if(length(sig.var) > 0){
		var.table <- var.influences[sig.var,,drop=FALSE]
		}else{
			stop("There are no significant interactions at this p.or.q value.")
			}

	for(j in 1:2){
		marker.locale <- match(var.table[,j], colnames(data.obj$geno))
		marker.names <- data.obj$marker.names[marker.locale]
		marker.block <- apply(matrix(var.table[,j], ncol = 1), 1, get.block)
		full.names <- paste(marker.block, marker.names, sep = ":")
		var.table[,j] <- full.names
		}

	
	final.table <- var.table[order(var.table[,"|Effect|/SE"], decreasing = TRUE),]

	if(mark.covar){
		covar.flags <- data.obj$covar.for.pairscan
		is.covar <- rownames(covar.flags)[unique(as.vector(unlist(apply(covar.flags, 2, function(x) which(x == 1)))))]
		covar.names <- data.obj$marker.names[match(is.covar, colnames(data.obj$geno))]
		covar.source.locale <- which(final.table[,1] %in% covar.names)
		covar.target.locale <- which(final.table[,2] %in% covar.names)
		if(length(covar.source.locale) > 0){
			final.table[covar.source.locale,1] <- paste(final.table[covar.source.locale,1], "*", sep = "")
			}
		if(length(covar.target.locale) > 0){
			final.table[covar.target.locale,2] <- paste(final.table[covar.target.locale,2], "*", sep = "")
			}
		}
	
	if(write.file){
		write.table(final.table, file = filename, quote = FALSE, sep = delim, row.names = FALSE)	
		invisible(final.table)
		}else{
			return(final.table)
			}
}