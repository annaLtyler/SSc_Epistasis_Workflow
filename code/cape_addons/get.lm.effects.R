#This function uses a linear model to calculate the
#main and combined effects of markers for three
#different genotype codings


get.lm.effects <- function(data.obj, marker1.name, marker2.name, pheno.name, scan.what = c("normalized.traits", "raw.traits", "eigentraits"), collapsed.net = FALSE, geno.coding = c("Additive", "Dominant", "Recessive")){

	#=========================================================
	#internal functions
	#=========================================================
	get.marker.name <- function(marker.name){
		if(length(grep("Chr", marker.name)) > 0){
			marker.locale <- which(names(data.obj$linkage.blocks.full) == marker.name)
			marker.label <- data.obj$linkage.blocks.full[[marker.locale]]
			return(marker.label)
			}else{
			if(is.null(allele)){
				full.name <- marker.name
				}else{
				full.name <- paste0(marker.name, "_", allele)
				}
			return(marker.name)
			}
		}
	#=========================================================

	pheno <- get.pheno(data.obj, scan.what, covar)
	
	if(collapsed.net){
		m1.markers <- data.obj$linkage.blocks.collapsed[[which(names(data.obj$linkage.blocks.collapsed) == m1)]]
		m2.markers <- data.obj$linkage.blocks.collapsed[[which(names(data.obj$linkage.blocks.collapsed) == m2)]]
		m1.ind <- match(m1.markers, colnames(geno))
		m2.ind <- match(m2.markers, colnames(geno))
		all.int <- data.obj$full.net[m1.ind, m2.ind, drop=FALSE]
		max.int <- which(abs(all.int) == max(abs(all.int)), arr.ind = TRUE)
		m1.geno <- geno[not.na.locale,m1.ind[max.int[1,1]]]
		m2.geno <- geno[not.na.locale,m2.ind[max.int[1,2]]]
		m1.name <- 	m1.markers[max.int[1,1]]
		m2.name <- 	m2.markers[max.int[1,2]]
		}else{
		m1.name <- get.marker.name(m1)
		m2.name <- get.marker.name(m2)
		m1.geno <- geno[not.na.locale,which(colnames(geno) == m1.name)]
		m2.geno <- geno[not.na.locale,which(colnames(geno) == m2.name)]
		}

	geno.bins <- c(0, 0.5, 1)
	if(geno.coding == "Dominant"){
		geno.bins <- c(0, 1)
		m1.geno[which(m1.geno >= 0.5)] <- 1
		m2.geno[which(m2.geno >= 0.5)] <- 1
		}
		
	if(geno.coding == "Recessive"){
		geno.bins <- c(0, 1)
		m1.geno[which(m1.geno <= 0.5)] <- 0
		m2.geno[which(m2.geno <= 0.5)] <- 0
		}

	pheno.locale <- which(colnames(pheno) == pheno.name)
	pheno.vals <- pheno[,pheno.locale]

	full.model <- lm(pheno[,pheno.locale]~m1.geno*m2.geno)
	
	m1.coef <- coefficients(full.model)[2] - coefficients(full.model)[1]
	m2.coef <- coefficients(full.model)[3] - coefficients(full.model)[1]
	int.term <- coefficients(full.model)[4] - coefficients(full.model)[1]
	add.prediction <- m1.coef + m2.coef
	actual.int <- m1.coef + m2.coef + int.term
			
	return(c("m1.name" = m1.name, "m2.name" = m2.name, "m1.effect" = m1.coef, "m2.effect" = m2.coef, "double.effect" = actual.int))
	}