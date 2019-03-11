#This function plots significant interactions


plot.interaction.heatmaps <- function(data.obj, p.or.q = 0.05, pdf.label = "Interaction.Heatmaps.pdf", raw.pheno.normalized = TRUE){
	
	var.influences <- data.obj$var.to.var.p.val
	
	marker.names <- matrix(NA, ncol = 2, nrow = dim(var.influences)[1])
	marker.names[,1] <- data.obj$marker.names[match(var.influences[,1], colnames(data.obj$geno))]
	marker.names[,2] <- data.obj$marker.names[match(var.influences[,2], colnames(data.obj$geno))]	

	pheno.names <- colnames(data.obj$pheno)
	num.pheno <- length(pheno.names)
	num.et <- dim(data.obj$ET)[2]
	et.names <- colnames(data.obj$ET)

	if(raw.pheno.normalized){
		all.pheno <- cbind(data.obj$pheno, data.obj$ET)
		}else{
		#use the uncentered, unnormalized data
		all.pheno <- cbind(data.obj$raw.pheno, data.obj$ET)
		}
		
	all.pheno.names <- c(pheno.names, et.names)
	total.response <- num.et + num.pheno
	
	if(is.null(var.influences)){
		stop("variant-to-variant influences must be calculated before plotting.")
		}

	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences))))
	sig.int.locale <- which(var.influences[,var.sig.col] <= p.or.q)
	
	if(length(sig.int.locale) == 0){
		
		pdf(pdf.label)
		plot.new()
		plot.window(xlim = c(0,1), ylim = c(0,1))
		text(0.5, 0.5, "No Significant Interactions")
		dev.off()
		}else{
		
		num.int <- length(sig.int.locale)
		if(num.int > 3){
			num.per.page <- 3
			}else{
				num.per.page <- num.int
				}
	
		mean.table <- function(marker1, marker2, pheno){
			geno.table <- table(marker1, marker2)
			geno.pairs <- pair.matrix(colnames(geno.table), ordered = TRUE, self.pairs = TRUE)
			mean.vals <- apply(geno.pairs, 1, function(x) mean(pheno[intersect(which(marker1 == x[1]), which(marker2 == x[2]))], na.rm = TRUE))
			for(i in 1:length(geno.pairs[,1])){
				geno.table[as.character(geno.pairs[i,1]), as.character(geno.pairs[i,2])] <- mean.vals[i]
				}
			return(geno.table)
			}
	
	
		pdf(pdf.label)
		for(i in sig.int.locale){
			marker1 <- data.obj$geno[,as.character(var.influences[i,1])]
			marker2 <- data.obj$geno[,as.character(var.influences[i,2])]

			for(ph in 1:total.response){
				val.table <- mean.table(marker1, marker2, pheno = all.pheno[,ph])
				myImagePlot(val.table, min.x = (max(abs(val.table))*-1), max.x = max(abs(val.table)), main = paste(marker.names[i,1], ", ", marker.names[i,2], "\n", all.pheno.names[ph], sep = ""))
				}	
			}

		dev.off()
		}
}

