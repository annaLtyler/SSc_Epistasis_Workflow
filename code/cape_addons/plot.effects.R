#This function takes in the data object, a marker name, and a phenotype
#and plots phenotype means and standard deviations for each factor in
#the genotype
#if a second marker is added, the function plots an interaction plot
#for the two markers.
#marker = colnames(data.obj$geno)[1]; marker2 = colnames(data.obj$geno)[2]; pheno.transform = "Normalized"; plot.type = "l"; error.bars = "none"; ymin = NULL; ymax = NULL; covar = data.obj$p.covar; marker1.label = ""; marker2.label = ""; num.rows = 1; ref.centered = TRUE; sort.by = "genotype"; pheno = NULL; allele.subset = "none";pos.col = "brown"; neg.col = "blue"

plot.effects <- function(data.obj, marker = NULL, marker2 = NULL, pheno.transform = "Normalized", plot.type = c("l", "p", "b", "a"), error.bars = "none", ymin = NULL, ymax = NULL, covar = NULL, marker1.label = "", marker2.label = "", num.rows = 1, ref.centered = TRUE, sort.by = "genotype", pheno = NULL, allele.subset = "none", pos.col = "brown", neg.col = "blue", gen.model = "Additive", report.p = c("empirical", "adjusted"), pheno.labels = NULL, genotype.labels = NULL, cex.axis = 1.5, cex.text = 1.5){




	if(length(grep("m", report.p[1], ignore.case = TRUE)) == 1){pval.col = 7}
	if(length(grep("j", report.p[1], ignore.case = TRUE)) == 1){pval.col = 8}

	type <- grep("l", plot.type)
	if(length(type) == 1){
		plot.type = "l"
		}
		
	if(plot.type == "8" || plot.type == "36" || plot.type == "2"){
		pheno.boxes(data.obj, markers = c(marker, marker2), sort.by = sort.by, pheno = pheno, allele.subset = allele.subset, covar = covar, states = plot.type)
		}
	
	if(plot.type == "h"){
		inter.heat(data.obj, marker1, marker2, pheno, allele.subset, covar, pos.col, neg.col)
		}
	
	centering <- grep("u", plot.type)
	if(length(centering) == 1){
		plot.type <- "b"
		ref.centered = FALSE
		}
	
	if(is.null(marker) && is.null(marker2)){
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		text(0.5, 0.5, "No Markers Selected")
		}else{

	
		#============================================================
		#get the marker names
		#============================================================

		# markers <- colnames(data.obj$geno)[match(c(marker, marker2), data.obj$marker.names)]
		markers <- c(marker, marker2)
		if(is.null(marker)){marker1.label = ""}
		if(!is.null(marker) && marker1.label == ""){marker1.label = marker}
		if(is.null(marker2)){marker2.label = ""}
		if(!is.null(marker2) && marker2.label == ""){marker2.label = marker2}
		
		marker.names <- c(marker1.label, marker2.label)
		#============================================================
		#transform the phenotype based on user input
		#============================================================


		
		all.pheno.mat <- get.pheno(data.obj, scan.what = pheno.transform, covar = covar)
		
		if(is.null(pheno.labels)){
			pheno.names <- colnames(all.pheno.mat)
			}else{	
			pheno.names <- pheno.labels
			}
		

		#============================================================
		#bin the genotypes if necessary
		#============================================================
		if(gen.model == "Additive"){geno.bins <- c(0, 0.5, 1)}
		if(gen.model == "Dominant"){geno.bins <- c(0, 1)}
		if(gen.model == "Recessive"){geno.bins <- c(0, 1)}
		
		marker.locale <- match(markers, colnames(data.obj$geno.for.pairscan))
	
		marker.mat <- data.obj$geno.for.pairscan[,marker.locale,drop=FALSE]
		
		if(gen.model == "Dominant"){
			marker.mat[which(marker.mat >= 0.5)] <- 1
			marker.mat[which(marker.mat < 0.5)] <- 0
			}
		if(gen.model == "Recessive"){
			marker.mat[which(marker.mat <= 0.5)] <- 0			
			}
		
		if(length(which(is.na(marker.locale))) > 0){ #if there are covariates selected
			marker.covar.locale <- match(markers, data.obj$p.covar)
			marker.mat.c <- data.obj$p.covar.table[,marker.covar.locale[which(!is.na(marker.covar.locale))],drop=FALSE]
			colnames(marker.mat.c) <- data.obj$p.covar[marker.locale[which(!is.na(marker.locale))]]
			marker.mat[,which(is.na(marker.locale))] <- marker.mat.c
			}
			
		
		for(m in 1:length(markers)){
			marker.geno.bins <- sort(unique(marker.mat[which(!is.na(marker.mat[,m])),m]))
			if(length(marker.geno.bins) > 3){
				marker.mat[,m] <- bin.vector(marker.mat[,m], geno.bins)
				}
			}
		#============================================================		
		
		
		#============================================================
		#figure out the layout for the plots
		#============================================================
		num.col <- ceiling(length(pheno.names)/num.rows)
		plot.cells <- num.col*num.rows
		na.padding <- rep(0, plot.cells - length(pheno.names))
		layout.mat <- matrix(c(1:length(pheno.names), na.padding), nrow = num.rows)
		layout(layout.mat)
		#============================================================
		
		
		for(ph in 1:dim(all.pheno.mat)[2]){
			pheno.influence <- data.obj$max.var.to.pheno.influence[[ph]]
			
			marker.pvals <- rep(NA, length(markers))
			for(m in 1:length(markers)){
				marker.main.locale <- which(pheno.influence[,1] == markers[m])
				marker.pvals[m] <- as.numeric(pheno.influence[marker.main.locale,pval.col])
				}
			
			if(plot.type == "l"){
				plot.lines(data.obj, marker.geno = marker.mat, marker.name = marker.names, phenotype.vals = all.pheno.mat[,ph], phenotype.name = pheno.names[ph], ymin = ymin, ymax = ymax, error.bars = error.bars, marker.pvals = marker.pvals, ref.centered = ref.centered)
				}
			if(plot.type == "p"){
				plot.points(marker.mat, marker.names, all.pheno.mat[,ph], pheno.names[ph], ymin = ymin, ymax = ymax)
				}
			if(plot.type == "b"){
				plot.bars(marker.geno = marker.mat,  marker.name = marker.names, phenotype.vals = all.pheno.mat[,ph], pheno.name = pheno.names[ph], ref.centered = ref.centered, error.bars = error.bars, ymin = ymin, ymax = ymax, marker.pvals = marker.pvals, genotype.symbols = genotype.labels, cex.axis = cex.axis, text.cex = cex.text)
			}
			
		} #end looping through phenotypes

	} #end case for when there are no markers being input

} #end function