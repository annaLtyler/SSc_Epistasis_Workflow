#This function plots boxplots for phenotype distributions
#split by the designated group, such as sex or diet

pheno.boxplots <- function(data.obj, split.by, group.labels = NULL, pdf.label = "Phenotype.Boxplots"){

	library(RColorBrewer)
	cols <- brewer.pal(8, "Accent")
	
	all.pheno <- data.obj$pheno
	all.geno <- data.obj$geno
	split.col <- match(split.by, data.obj$marker.names)
	all.groups <- all.geno[,split.col]
	
	if(length(split.by) > 1){
		all.groups <- apply(all.groups, 1, function(x) paste(x, collapse = ","))
		groups <- sort(unique(all.groups))
		if(length(group.labels) != length(groups)){
			group.labels <- groups
			}
		}else{
		groups <- unique(all.groups)	
		}


	if(length(group.labels) == 0){
		group.labels <- groups
		}
	
		
	plot.box <- function(pheno.num){
		if(length(groups) > 2){
			pval <- kruskal.test(all.pheno[,pheno.num]~as.factor(all.groups))$p.value
			
			if(length(unique(all.groups[which(!is.na(all.pheno[,pheno.num]))])) == length(group.labels)){
				if(pval <= 0.05){
					boxplot(all.pheno[,pheno.num]~all.groups, notch = TRUE, main = colnames(all.pheno)[pheno.num], col = cols[2], names = group.labels)
					}else{
					boxplot(all.pheno[,pheno.num]~all.groups, notch = TRUE, main = colnames(all.pheno)[pheno.num], names = group.labels)	
					}
				}
			}
	
		if(length(groups) == 2){
			group1.num <- length(which(!is.na(all.pheno[which(all.groups == groups[1]),pheno.num])))
			group2.num <- length(which(!is.na(all.pheno[which(all.groups == groups[2]),pheno.num])))
			if(group1.num > 2 && group2.num > 2){
				pval <- wilcox.test(all.pheno[which(all.groups == groups[1]),pheno.num], all.pheno[which(all.groups == groups[2]),pheno.num])$p.value
				mean1 <- median(all.pheno[which(all.geno[,split.col] == groups[1]), pheno.num], na.rm = TRUE)
				mean2 <- median(all.pheno[which(all.geno[,split.col] == groups[2]), pheno.num], na.rm = TRUE)
			
				if(pval[1] <= 0.05 && mean1 < mean2){
					boxplot(all.pheno[,pheno.num]~all.geno[,split.col], notch = TRUE, main = colnames(all.pheno)[pheno.num], col = cols[1], names = group.labels, cex.axis = 2)
					}
				if(pval[1] <= 0.05 && mean2 < mean1){
					boxplot(all.pheno[,pheno.num]~all.geno[,split.col], notch = TRUE, main = colnames(all.pheno)[pheno.num], col = cols[2], names = group.labels, cex.axis = 2)
					}
				if(pval[1] > 0.05){
					boxplot(all.pheno[,pheno.num]~all.geno[,split.col], notch = TRUE, main = colnames(all.pheno)[pheno.num], names = group.labels, cex.axis = 2)		
					}
				}
			}
		}
	
	if(dim(all.pheno)[2] > 9){
		layout.matrix <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
		}else{
		layout.matrix <- get.layout.mat(dim(all.pheno)[2])
		}
	
	split.by.text <- paste(split.by, collapse = "_")
	pdf(paste(pdf.label, ".Split.By.", split.by.text, ".pdf", sep = ""), width = 9, height = 9)
	layout(layout.matrix)
	apply(matrix(1:dim(all.pheno)[2], ncol = 1), 1, plot.box)
	dev.off()
	
}