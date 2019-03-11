#This function plot phenotype correlations in 
#panels.

plot.pheno.by.ind <- function(data.obj, geno.obj = NULL, pheno.which = NULL, color.by = NULL, group.labels = NULL, text.cex = 1, pheno.labels = NULL){

	library(RColorBrewer)
	
	all.pheno <- data.obj$pheno
	all.geno <- get.geno(data.obj, geno.obj)
	
	if(length(dim(all.geno)) == 3){is.do = TRUE}else{is.do = FALSE}
	
	if(is.null(pheno.which)){
		pheno.names <- colnames(all.pheno)
		}else{
			if(is.numeric(pheno.which)[1]){
				pheno.names <- colnames(all.pheno)[pheno.which]
				}else{
				pheno.names <- pheno.which	
				}
			}
	
	if(is.null(pheno.labels)){
		pheno.labels <- pheno.names
		}
		

	all.cols <- brewer.pal(8, "Accent")


	if(!is.null(color.by)){
		cols <- rep(NA, dim(all.pheno)[1])
		if(is.do){
			group.col <- which(colnames(data.obj$covar.table) %in% color.by)
			}else{
			group.col <- which(data.obj$p.covar %in% color.by)
			}
			
		if(length(group.col) == 0){
			stop(paste("I couldn't find the", color.by, "column. Please check the case and the spelling."))
			}
			
		if(is.do){
			group.mem <- data.obj$covar.table[,group.col]
			}else{
			group.mem <- data.obj$p.covar.table[,group.col]
			}
		groups <- sort(unique(group.mem))
		num.groups <- length(groups)
		if(num.groups > 8){
			stop("There cannot be more than 8 groups")
			}
		for(i in 1:num.groups){
			cols[which(group.mem == groups[i])] <- all.cols[i]
			}
		}else{
		cols <- rep("black", dim(all.geno)[1])
		groups <- NULL
		}

	if(is.null(group.labels)){
		group.labels <- groups
		}

		layout.mat <- get.layout.mat(dim(data.obj$pheno)[2])
		layout(layout.mat)
		for(p in 1:length(pheno.labels)){
			plot(data.obj$pheno[,p], col = cols, pch = 16, xlab = "Individual", ylab = pheno.labels[p], main = pheno.labels[p])
			legend("topleft", pch = 16, col = all.cols[1:length(groups)], legend = group.labels)			
			}
		
	}
