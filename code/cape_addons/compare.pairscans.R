compare.pairscans <- function(pairscan1, pairscan2, label1 = "", label2 = ""){
	
	pairscan1.results <- pairscan1$pairscan.results
	pairscan2.results <- pairscan2$pairscan.results	
	
	layout.mat <- get.layout.mat(3*length(pairscan1.results))
	
	# quartz(width = ncol(layout.mat)*3, height = nrow(layout.mat)*3)
	layout(layout.mat)
	for(i in 1:length(pairscan1.results)){
		marker.pairs1 <- apply(pairscan1.results[[i]][[1]][,1:2], 1, function(x) paste(x, collapse = ","))
		marker.pairs2 <- apply(pairscan2.results[[i]][[1]][,1:2], 1, function(x) paste(x, collapse = ","))
		
		common.pairs <- intersect(marker.pairs1, marker.pairs2)
		common.locale1 <- match(common.pairs, marker.pairs1)
		common.locale2 <- match(common.pairs, marker.pairs2)

		for(j in 3:5){
		std.stat1 <- abs(pairscan1.results[[i]][[1]][common.locale1,j]/pairscan1.results[[i]][[2]][common.locale1,j])
		std.stat2 <- abs(pairscan2.results[[i]][[1]][common.locale2,j]/pairscan2.results[[i]][[2]][common.locale2,j])
		plot(std.stat1, std.stat2, main = paste(names(pairscan1.results)[i], "\n", colnames(pairscan1.results[[i]][[1]])[j]), xlab = label1, ylab = label2)
		abline(0,1, col = "red")
		}
	
	}

}