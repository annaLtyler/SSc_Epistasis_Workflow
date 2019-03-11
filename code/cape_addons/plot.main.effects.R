#This script plots significant main effects of a cross, even if there is no interaction effect

plot.main.effects <- function(data.obj, p.or.q = 0.05, pdf.label = "Main.Effect.Plots.pdf", pheno.normalized = TRUE){
	
	
	pheno.names <- colnames(data.obj$pheno)
	
	var.influences <- data.obj$max.var.to.pheno.influence
	pheno.names <- colnames(data.obj$pheno)

	
	if(is.null(var.influences)){
		stop("variant-to-phenotype influences must be calculated before plotting.")
		}


	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences[[1]]))))
	
	sig.table <- NULL
	for(i in 1:length(var.influences)){
		sig.locale <- which(var.influences[[i]][,var.sig.col] <= p.or.q)
		table.section <- cbind(rep(pheno.names[i], length(sig.locale)), var.influences[[i]][sig.locale,])
		sig.table <- rbind(sig.table, table.section)
		}
	
	if(is.null(sig.table)){
		stop("There are no significant main effects.")
		}
	
	u_markers <- as.character(sort(as.numeric(unique(sig.table[,"marker"]))))
	marker.names <- data.obj$marker.names[match(u_markers, colnames(data.obj$geno))]
	
	num.total <- length(u_markers)*length(pheno.names)
	if(num.total <= 4){
		layout.mat <- matrix(1:num.total, ncol = length(pheno.names), byrow = TRUE)
		}else{
		layout.mat <- matrix(1:4,  ncol = length(pheno.names), byrow = TRUE)
		}
	
	pdf(pdf.label, width = dim(layout.mat)[2]*3, height = dim(layout.mat)[1]*3)
	layout(layout.mat)
	#for each of the markers with significant main effects plot effects on all phenotypes
	for(i in 1:length(u_markers)){
		geno <- data.obj$geno[,u_markers[i]]
		u_geno <- unique(geno); u_geno <- sort(u_geno[which(!is.na(u_geno))])
		geno.locale <- apply(matrix(u_geno, ncol = 1), 1, function(x) which(geno == x))

		if(class(geno.locale) == "list"){
			for(p in 1:length(pheno.names)){
				if(pheno.normalized){
					pheno.vals <- lapply(geno.locale, function(x) data.obj$pheno[x,p])
					}else{
					pheno.vals <- lapply(geno.locale, function(x) data.obj$raw.pheno[x,p])			
					}
					pheno.means <- unlist(lapply(pheno.vals, function(x) mean(x, na.rm = TRUE)))
					pheno.se <- unlist(lapply(pheno.vals, function(x) sd(x, na.rm = TRUE)/sqrt(length(which(!is.na(x))))))
					ylim = c(min(pheno.means-pheno.se), max(pheno.means+pheno.se))
					x.vals <- segment.region(0.1, 0.9, length(u_geno), alignment = "ends")
					plot(x = x.vals, y = lapply(pheno.vals, mean), xlim = c(0,1), ylim = ylim, type = "l", axes = FALSE, xlab = "Genotype", ylab = pheno.names[p], main = paste(marker.names[i], pheno.names[p]))
					axis(1, at = x.vals, labels = u_geno)
					axis(2)
					segments(x0 = x.vals, y0 = pheno.means - pheno.se, x1 = x.vals, y1 = pheno.means + pheno.se)
					}
				}
			}
	dev.off()

	
	
	
}