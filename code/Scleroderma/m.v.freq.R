#correlate m12/m21 values with allele frequency

m.v.freq <- function(data.obj, singlescan.obj){
	library(rgl)
	var.to.var <- data.obj$var.to.var.p.val
	m.table <- data.obj$var.to.var.influences
	
	main.table <- data.obj$max.var.to.pheno.influence
	
	get.allele.freq <- function(marker.name){
		marker.locale <- which(colnames(data.obj$geno.for.pairscan) == marker.name)
		allele.freq <- mean(data.obj$geno.for.pairscan[,marker.locale], na.rm = TRUE)
		if(allele.freq < 0.5){
			return(allele.freq)
			}else{
			return(1-allele.freq)	
			}
		}
		
	source.freq <- apply(matrix(var.to.var[,1], ncol = 1), 1, get.allele.freq)
	target.freq <- apply(matrix(var.to.var[,2], ncol = 1), 1, get.allele.freq)

	quartz(width = 10, height = 5)
	par(mfrow = c(1,2))
	plot(source.freq, var.to.var[,5], xlab = "Minor Allele Frequency of Source Marker", ylab = "Standardized Effect Size of Interaction")
	plot(target.freq, var.to.var[,5], xlab = "Minor Allele Frequency of Target Marker", ylab = "Standardized Effect Size of Interaction")
	
	quartz()
	layout.mat <- get.layout.mat(length(main.table), "landscape")
	layout(layout.mat)
	for(i in 1:length(main.table)){
		marker.freq <- apply(matrix(main.table[[i]][,1], ncol = 1), 1, get.allele.freq)
		plot(marker.freq, main.table[[i]][,5], xlab = "Allele Frequency", ylab = "|Standardized Main Effect|", main = names(main.table)[i])
		}

	quartz()
	layout.mat <- get.layout.mat(length(singlescan.obj$singlescan.results), "landscape")
	layout(layout.mat)
	if(!missing(singlescan.obj)){
		for(i in 1:length(singlescan.obj$singlescan.results)){
			marker.freq <- apply(matrix(rownames(singlescan.obj$singlescan.results[[i]]), ncol = 1), 1, get.allele.freq)
			plot(marker.freq, singlescan.obj$singlescan.results[[i]][,"t.stat"], ylab = "Standardized Main Effect", xlab = "Allele Frequency")
			}
		}

}

