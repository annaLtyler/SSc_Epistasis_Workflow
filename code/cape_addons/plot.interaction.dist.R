#This script plots the histogram for var to var effects
#shade.over indicates where you would like to shade the 
#empirical distribution.
#You can shade it over a given percentile of the null 
#distribution. The default is to shade the portion of
#the empirical distribution over the 95th percentil of
#the null distribution.
#You can also shade based on a Holm's adjusted p value
#To do this, set shade.over to c.pval and set per.or.p
#to the adjusted p value. This will show any part of
#the distribution that is above the effect level corresponding
#to the adjusted p value. 
#markers <- c("D2Mit148", "D18Mit140")
#markers <- c("D6Mit275", "D17Mit240")
#markers <- c("D1Mit76", "D12Mit132")
plot.interaction.dist <- function(data.obj, markers, shade.over = c("percentile", "c.pval"), perc.or.p = 95){
	
	
	var.inf <- data.obj$var.to.var.influences
	var.inf.perm <- data.obj$var.to.var.influences.perm
	m12.m21 <- c(abs(as.numeric(var.inf[,"m12"])/as.numeric(var.inf[,"m12.std.dev"])), abs(as.numeric(var.inf[,"m21"])/as.numeric(var.inf[,"m21.std.dev"])))
	
	marker.locale <- match(data.obj$marker.names, markers)
	marker.order <- marker.locale[!is.na(marker.locale)]
	markers.num <- colnames(data.obj$geno)[which(!is.na(marker.locale))[marker.order]]

	if(length(grep("perc", shade.over)) > 0){
		shade.above.percentile <- TRUE
		}else{
		shade.above.percentile <- FALSE
		}


	marker1.position1 <- which(var.inf[,1] %in% markers.num[1])
	marker1.position2 <- which(var.inf[,2] %in% markers.num[1])
	if(length(marker1.position1) == 0 && length(marker1.position2) == 0){
		stop(paste(markers[1], "was not tested."))
		}
	marker2.position1 <- which(var.inf[,1] %in% markers.num[2])
	marker2.position2 <- which(var.inf[,2] %in% markers.num[2])
	if(length(marker1.position1) == 0 && length(marker1.position2) == 0){
		stop(paste(markers[2], "was not tested."))
		}
	
	pair.locale <- c(intersect(marker1.position1, marker2.position2), intersect(marker2.position1, marker1.position2))
			
			
		if(length(pair.locale) == 0){
			stop("The markers were not tested together.")
			table(data.obj$geno[,markers.num[1]], data.obj$geno[,markers.num[2]])
			}
			
		var.inf[pair.locale,]
	
		layout.mat <- get.layout.mat(2, "landscape")



		null.dist <- c(abs(as.numeric(var.inf.perm[,"m12"])/as.numeric(var.inf.perm[,"m12.std.dev"])), abs(as.numeric(var.inf.perm[,"m21"])/as.numeric(var.inf.perm[,"m21.std.dev"])))
	
		m12 <- abs(as.numeric(var.inf[pair.locale,"m12"])/as.numeric(var.inf[pair.locale,"m12.std.dev"]))
		m21 <- abs(as.numeric(var.inf[pair.locale,"m21"])/as.numeric(var.inf[pair.locale,"m21.std.dev"]))
			
		xmin <- min(density(null.dist)$x); xmax <- max(c(density(null.dist)$x, m12, m21))
		ymax <- max(c(density(null.dist)$y)*1.15)
		
				
		
		plot.null.with.arrows <- function(null, arrow.locale, main.label){
			plot(density(null), xlim  = c(xmin, xmax), xlab = "Direct Influence t Statistic", main = main.label, ylim = c(0, ymax))
			#add arrows for m12 and m21
			arrows(arrow.locale, ymax*0.05, arrow.locale, 0 , col = "red", lwd = 3, length = 0.1)

			#add a shaded region to the emp dist above the percentile in the null dist
			if(shade.above.percentile){
				perc <- as.vector(quantile(null, perc.or.p/100))
				
				}else{
				#get the p values of the distribution
				all.p <- apply(matrix(arrow.locale, ncol = 1), 1, function(x) length(which(null.dist >= x))/length(null.dist))
				c.p <- p.adjust(all.p, method = "holm")
				sig.locale <- which(c.p <= perc.or.p)
				if(length(sig.locale) > 0){
					perc <- arrow.locale[min(sig.locale)]
					}else{
						perc <- max(null.dist)
						}
				}
				
			
			null.dist.greater <- which(density(null)$x >= perc)
			
			if(length(null.dist.greater) > 0){
				x1 <- min(density(null)$x[null.dist.greater])
				x2 <- max(density(null)$x[null.dist.greater])
				polygon(x = c(x1, density(null)$x[null.dist.greater], x2), y = c(0, density(null)$y[null.dist.greater], 0), col = "blue")
				# arrows(perc, ymax*0.05, perc, 0 , col = "red", lwd = 3, length = 0.1)
				}
			
			if(shade.above.percentile){
				legend.label <- paste("Above", perc.or.p, "percentile")
				}else{
				legend.label <- paste("Below p =", perc.or.p)
				}
			legend("topleft", legend = c("Observed Interaction", legend.label), fill = c("red", "blue"), cex = 0.6)
			
			}
	
	dev.new(width = 10, height = 5)
	layout(layout.mat)

	main.label <- paste(markers[1], "<<", markers[2])
	plot.null.with.arrows(null.dist, m12, main.label)
	
	main.label <- paste(markers[1], ">>", markers[2])	
	plot.null.with.arrows(null.dist, m21, main.label)	

	
	
		
	
}