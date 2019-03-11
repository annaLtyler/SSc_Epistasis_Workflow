#This plots the main effects vs. the interaction effects for DO data

main.v.int <- function(data.obj, singlescan.obj){
	
	int.data <- data.obj$var.to.var.influences
	main.data <- singlescan.obj$singlescan.t.stats
	
	u_markers <- unique(c(int.data[,1], int.data[,2]))
	overall.int <- c(as.numeric(int.data[,3])/as.numeric(int.data[,4]), as.numeric(int.data[,5])/as.numeric(int.data[,6]))
	
	split.markers <- strsplit(u_markers, "_")
	marker.names <- unlist(lapply(split.markers, function(x) x[1]))
	allele.names <- unlist(lapply(split.markers, function(x) x[2]))
	allele.names[which(is.na(allele.names))] <- "A"
	
pdf("Main.Effects.V.Interaction.Effects.pdf")
	for(j in 1:dim(main.data)[2]){
		cat("Working on ET", j, "\n")
		plot.new()
		plot.window(xlim = c(min(main.data), max(main.data)), ylim = c(min(overall.int), max(overall.int)))
		axis(1);axis(2)

		for(i in 1:length(u_markers)){
			report.progress(i, length(u_markers))
			source.locale <- which(int.data[,1] == u_markers[i])
			target.locale <- which(int.data[,2] == u_markers[i])
			source.int <- as.numeric(int.data[source.locale,"m12"])/as.numeric(int.data[source.locale,"m12.std.dev"])
			target.int <- as.numeric(int.data[source.locale,"m21"])/as.numeric(int.data[source.locale,"m21.std.dev"])
			# plot(density(source.int))
			# points(density(target.int), type = "l", col = "red")
			all.int <- c(source.int, target.int)
			
			main.locale <- which(dimnames(main.data)[[1]] == marker.names[i])
			allele.locale <- which(dimnames(main.data)[[3]] == allele.names[i])
			main.eff <- matrix(main.data[main.locale,,allele.locale], ncol = 3, nrow = length(all.int), byrow = TRUE)
			
			points(main.eff[,j], all.int)
			
			}
			
		}	

dev.off()	
}