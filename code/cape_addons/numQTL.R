numQTL <- function(data.obj, p.or.q, separate.covar = TRUE){
	
	data.obj <- get.network(data.obj, p.or.q = p.or.q)
	
	net <- data.obj$collapsed.net
	
	geno.ind <- 1:min(dim(net))
	pheno.ind <- (min(dim(net))+1):max(dim(net))
	just.geno <- net[geno.ind, geno.ind]
	just.pheno <- net[,pheno.ind]

	if(!separate.covar){
		covar.names <- data.obj$marker.names[which(data.obj$chromosome == 0)]
		covar.locale <- match(covar.names, colnames(net))
		just.covar.source <- just.geno[covar.locale,]
		just.covar.target <- just.geno[,covar.locale]
		covar.pheno <- net[covar.locale, pheno.ind]
		just.geno <- just.geno[-covar.locale,]
		just.geno <- just.geno[,-covar.locale]
		
		covar.as.source <- colnames(just.covar.source)[which(apply(just.covar.source, 2, sum) != 0)]
		covar.as.target <- rownames(just.covar.target)[which(apply(just.covar.target, 1, sum) != 0)]
		unique.covar.qtl <- length(unique(c(covar.as.source, covar.as.target)))
		
		}
	
	num.int <- length(which(just.geno != 0))
	neg.int <- length(which(just.geno < 0))
	pos.int <- length(which(just.geno > 0))

	sourceQTL <- colnames(just.geno)[which(apply(just.geno, 1, sum) != 0)]
	targetQTL <- colnames(just.geno)[which(apply(just.geno, 2, sum) != 0)]

	uniqueQTL <- length(unique(c(sourceQTL, targetQTL)))
	
	phenoQTL <- apply(just.pheno, 2, function(x) length(which(x != 0)))


	if(!separate.covar){
		table.report <- vector(mode = "list", length = 3)
		names(table.report) <- c("Interactions", "QTL", "Phenotoype QTL")

		table.report[[1]] <- matrix(c(num.int, neg.int, pos.int), ncol = 1)
		rownames(table.report[[1]]) <- c("Interactions", "Interactions Negative", "Interactions Positive")
		table.report[[2]] <- matrix(c(length(sourceQTL), length(targetQTL), uniqueQTL), ncol = 1)
		rownames(table.report[[2]]) <- c("QTL Source", "QTL Target", "QTL Unique")
		table.report[[3]] <- matrix(phenoQTL, ncol = 1)
		rownames(table.report[[3]]) <- paste("QTL", names(phenoQTL))

		}else{
		table.report <- vector(mode = "list", length = 4)
		names(table.report)  <- c("Genetic Interactions", "QTL", "Phenotoype QTL", "Covariate Interactions")

		table.report[[1]] <- matrix(c(num.int, neg.int, pos.int), ncol = 1)
		rownames(table.report[[1]]) <- c("Interactions", "Interactions Negative", "Interactions Positive")
		table.report[[2]] <- matrix(c(length(sourceQTL), length(targetQTL), uniqueQTL), ncol = 1)
		rownames(table.report[[2]]) <- c("QTL Source", "QTL Target", "QTL Unique")
		table.report[[3]] <- matrix(phenoQTL, ncol = 1)
		rownames(table.report[[3]]) <- paste("QTL", names(phenoQTL))

		table.report[[4]] <- matrix(c(length(covar.as.source), length(covar.as.target), unique.covar.qtl), ncol = 1)
		rownames(table.report[[4]]) <- c("covar -> QTL", "QTL -> covar", "unique QTL")

		}
	
	

	return(table.report)
	}