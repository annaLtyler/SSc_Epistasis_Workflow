

plot.pheno.cor <- function(data.obj, pdf.label = "Phenotype.Correlations.pdf", pheno.which = NULL, color.by = NULL, group.labels = NULL){

	library(RColorBrewer)
	
	all.pheno <- data.obj$pheno
	all.geno <- data.obj$geno
	
	if(is.null(pheno.which)){
		pheno.names <- colnames(all.pheno)
		}else{
			if(is.numeric(pheno.which)[1]){
				pheno.names <- colnames(all.pheno)[pheno.which]
				}else{
				pheno.names <- pheno.which	
				}
			}
			
	pairs.matrix <- pair.matrix(pheno.names)

	all.cols <- brewer.pal(8, "Accent")


	if(!is.null(color.by)){
		cols <- rep(NA, dim(all.pheno)[1])
		group.col <- which(data.obj$p.covar %in% color.by)
		group.mem <- data.obj$p.covar.table[,group.col]
		groups <- sort(unique(group.mem))
		num.groups <- length(groups)
		if(num.groups > 8){
			stop("There cannot be more than 8 groups")
			}
		for(i in 1:num.groups){
			cols[which(all.geno[,group.col] == groups[i])] <- all.cols[i]
			}
		}else{
		cols <- rep("black", dim(all.geno)[1])
		groups <- NULL
		}

	if(is.null(group.labels)){
		group.labels <- groups
		}

	plot.cor <- function(pair){
		if(length(which(!is.na(all.pheno[,pair[1]]+all.pheno[,pair[2]]))) > 0){
			if(!is.null(groups)){
				group.cor <- try(apply(matrix(groups, ncol = 1), 1, function(x) cor(all.pheno[which(group.mem == x),pair[1]], all.pheno[which(group.mem == x),pair[2]], use = "complete")), silent = TRUE)
				}
				pheno.cor <- try(cor(all.pheno[,pair[1]], all.pheno[,pair[2]], use = "complete"), silent = TRUE)

			
			cor.label <- paste("r =", round(pheno.cor, 2))
			pheno.label <- paste(pair[1], pair[2], "\n", sep = " ")
				
			if(!is.null(groups)){
				cor.group.label <- paste(c(group.labels, "total"), "r =", round(c(group.cor, pheno.cor), 2), collapse = "; ")
				if(class(group.cor) != "try-error" && class(pheno.cor) != "try.error"){
					plot(all.pheno[,pair[1]], all.pheno[,pair[2]], col = cols, pch = 16, xlab = pair[1], ylab = pair[2], main = paste(cor.group.label, pheno.label, sep = "\n"))
					legend("topleft", legend = group.labels, fill = all.cols[1:num.groups])
					}
				}else{
					plot(all.pheno[,pair[1]], all.pheno[,pair[2]], col = cols, pch = 16, xlab = pair[1], ylab = pair[2], main = paste(cor.label, pheno.label, sep = "\n"))
					}
			}

		}

	if(length(pairs.matrix[,1]) > 9){
		layout.mat <- matrix(1:9, nrow = 3, ncol = 3)
		}else{
		layout.mat <- get.layout.mat(length(pairs.matrix[,1]), type = "landscape")	
		}

	pdf(pdf.label, width = max(dim(layout.mat)[2]*3, 5), height = max(dim(layout.mat)[1]*3, 5))
	layout(layout.mat)
	apply(pairs.matrix, 1, plot.cor)
	dev.off()
	
	}
