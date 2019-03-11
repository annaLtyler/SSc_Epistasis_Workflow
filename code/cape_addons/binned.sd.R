#This function gets the standard deviation of
#Y values along bins of X

binned.sd <- function(X, Y, n.bins = 10){
	
	bin.placement <- segment.region(min(X, na.rm = TRUE), max(X, na.rm = TRUE), num.points = n.bins, alignment = "ends")

	binned.V <- vector(mode = "list", length = n.bins-1)
	for(i in 1:(length(bin.placement)-1)){
		bin.locale <- intersect(which(X >= bin.placement[i]), which(X < bin.placement[i+1]))
		binned.V[[i]] <- Y[bin.locale]
		}

	bin.mean <- lapply(binned.V, mean)
	bin.sd <- lapply(binned.V, sd)
	bin.range <- consec.pairs(bin.placement)
	bin.mids <- apply(bin.range, 1, mean)
	final.table <- cbind(bin.range[,1], bin.mids, bin.range[,2], bin.mean, bin.sd)
	colnames(final.table) <- c("minX", "midX", "maxX", "meanY", "sdY")
	
	return(final.table)
	}