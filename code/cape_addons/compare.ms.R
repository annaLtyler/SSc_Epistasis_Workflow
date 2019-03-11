#This function compares the significant cape interactions
#between two runs of the same data


compare.ms <- function(data.obj1, data.obj2, label1 = "data.obj1", label2 = "data.obj2"){
	
	net1 <- data.obj1$full.net
	net2 <- data.obj2$full.net
	
	plot(net1, net2, xlab = label1, ylab = label2)
	
	# motifs1 <- find.motifs(data.obj1, collapsed.net = FALSE)
	# motifs2 <- find.motifs(data.obj2, collapsed.net = FALSE)
	
	
}