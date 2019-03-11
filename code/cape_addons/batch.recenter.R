#This function recenters batches to match the 
#overall mean, here batch.size should be one
#value determined by examining the results from
#batch.id.

batch.recenter <- function(v, batch.size, standardize = FALSE){
	
	dev.new(height = 8, width = 8)
	par(mfrow = c(2,2))
	segs <- round(segment.region(1,length(v), round(length(v)/batch.size), alignment = "ends"))
	# plot(segs, rep(0.5, length(segs)), xlim = c(0, length(v)))
	cons.seg <- consec.pairs(segs)
	batches <- apply(cons.seg, 1, function(x) list(v[x[1]:(x[2]-1)]))
	batch.list <- lapply(batches, function(x) x[[1]])
	#add on the last point
	batch.list[[length(batch.list)]] <- c(batch.list[[length(batch.list)]], v[length(v)])
	
	boxplot(batch.list, main = "Before Centering", notch = TRUE)
	plot(v, main = "Uncentered Data By Individual")

	#recenter each batch independently
	mean.centered <- batch.list
	if(standardize){
		mean.centered <- lapply(mean.centered, function(x) x/sd(x, na.rm = TRUE))
		}
	mean.centered <- lapply(mean.centered, function(x) x - mean(x, na.rm = TRUE))
	boxplot(mean.centered, main = "After Centering", notch = TRUE)
	
	new.v <- unlist(mean.centered)
	plot(new.v, main = "Recentered Data By Individual")

	return(new.v)	
}