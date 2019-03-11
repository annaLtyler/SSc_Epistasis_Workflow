#This function pulls out the top n interactors
#from a cape data object


top.interactor.info <- function(data.obj, num.pairs = 25, r2.thresh = 0.8){
	
	var.influences <- data.obj$var.to.var.p.val
	previous.r2 <- data.obj$r2.thresh
	if(is.null(previous.r2) || previous.r2 != r2.thresh){
		cat("Calculating linkage blocks...\n")
		data.obj <- all.linkage.blocks(data.obj, r2.thresh = r2.thresh)
		linkage.blocks <- data.obj$linkage.blocks
		}

	interactions.of.interest <- var.influences[1:num.pairs,]
	
	get.marker.info <- function(marker.num){
		marker.mat <- matrix(NA, ncol = 3, nrow = length(marker.num))
		colnames(marker.mat) <- c("marker", "chromosome", "linkage.block")
		marker.locale <- match(marker.num, colnames(data.obj$geno))
		marker.mat[,1] <- data.obj$marker.names[marker.locale]
		marker.mat[,2] <- data.obj$chromosome[marker.locale]
		marker.mat[,3] <- names(linkage.blocks)[apply(matrix(marker.num, ncol = 1), 1, function(x) min(grep(x, linkage.blocks)))]
		return(marker.mat)
		}
	

	marker.list <- vector(mode = "list", length = 2)
	for(i in 1:2){
		marker.list[[i]] <- get.marker.info(interactions.of.interest[,i])
		}
	

	final.table <- cbind(marker.list[[1]], marker.list[[2]], interactions.of.interest[,3:dim(interactions.of.interest)[2]])
	write.table(final.table, paste("Top.", num.pairs, ".Interactions.r2.thresh.", r2.thresh, ".txt", sep = ""), row.names = FALSE, sep = "\t", quote = FALSE)

}