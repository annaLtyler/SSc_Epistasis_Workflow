#This function looks for batch effects in a vector of quantities

batch.id <- function(v, batch.size = seq(100, 500, 10), pdf.label = "Batch"){
	
	if(length(batch.size) < 9){
		layout.mat <- get.layout.mat(length(batch.size), "upright")
		}else{
		layout.mat <- matrix(1:9, 3,3, byrow = TRUE)
		}
	
	pdf(paste(pdf.label, ".Boxes.pdf", sep = ""))
	layout(layout.mat)
	#get window means (non-overlapping)
	for(b in 1:length(batch.size)){
		segs <- round(segment.region(1,length(v), round(length(v)/batch.size[b]), alignment = "ends"))
		# plot(segs, rep(0.5, length(segs)), xlim = c(0, length(v)))
		cons.seg <- consec.pairs(segs)
		batches <- apply(cons.seg, 1, function(x) list(v[x[1]:(x[2]-1)]))
		batch.list <- lapply(batches, function(x) x[[1]])
		boxplot(batch.list, main = paste("batch size:", batch.size[b]), notch = TRUE)
		}
	dev.off()
	
	
	pdf(paste(pdf.label, ".images.pdf", sep = ""))
	for(b in 1:length(batch.size)){
		segs <- round(segment.region(1,length(v), round(length(v)/batch.size[b]), alignment = "ends"))
		# plot(segs, rep(0.5, length(segs)), xlim = c(0, length(v)))
		cons.seg <- consec.pairs(segs)
		batches <- apply(cons.seg, 1, function(x) list(v[x[1]:(x[2]-1)]))
		batch.list <- lapply(batches, function(x) x[[1]])
		num.row <- max(sapply(batch.list, length))
		batch.mat <- matrix(NA, nrow = num.row, ncol = length(batch.list))
		for(i in 1:length(batch.list)){
			batch.mat[1:length(batch.list[[i]]), i] <- batch.list[[i]]
			}
		image(batch.mat, main = paste("batch size:", batch.size[b]))
		}
	dev.off()
}