#This function gets the mean phenotype effects for
#two markers and their combined effect using three
#different patterns of genotype coding


get.pheno.effects <- function(data.obj, marker1.name, marker2.name, pheno.name, covar = NULL, scan.what = c("normalized.traits", "raw.traits", "eigentraits"), collapsed.net = FALSE, geno.coding = c("Additive", "Dominant", "Recessive")){

	geno <- data.obj$geno.for.pairscan

	covar.info <- get.covar(data.obj)
		if(!is.null(covar.info$covar.table)){
			not.na.locale <- which(!is.na(rowSums(covar.info$covar.table)))
			}else{
			not.na.locale <- 1:nrow(pheno)	
			}

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
		m1.markers <- data.obj$linkage.blocks.collapsed[[which(names(data.obj$linkage.blocks.collapsed) == marker1.name)]]
		m2.markers <- data.obj$linkage.blocks.collapsed[[which(names(data.obj$linkage.blocks.collapsed) == marker2.name)]]
		m1.ind <- match(m1.markers, colnames(geno))
		m2.ind <- match(m2.markers, colnames(geno))
		all.int <- data.obj$full.net[m1.ind, m2.ind, drop=FALSE]
		max.int <- which(abs(all.int) == max(abs(all.int)), arr.ind = TRUE)
		m1.geno <- geno[not.na.locale,m1.ind[max.int[1,1]]]
		m2.geno <- geno[not.na.locale,m2.ind[max.int[1,2]]]
		m1.name <- 	m1.markers[max.int[1,1]]
		m2.name <- 	m2.markers[max.int[1,2]]
		}else{
		m1.name <- get.marker.name(marker1.name)
		m2.name <- get.marker.name(marker2.name)
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

	m1.min.locale <- which(m1.geno == min(geno.bins))
	m1.max.locale <- which(m1.geno == max(geno.bins))
	m2.min.locale <- which(m2.geno == min(geno.bins))
	m2.max.locale <- which(m2.geno == max(geno.bins))

	no.m1.no.m2 <- intersect(m1.min.locale, m2.min.locale)
	baseline.val <- mean(pheno.vals[no.m1.no.m2])

	just.m1 <- intersect(m1.max.locale, m2.min.locale)
	m1.coef <- mean(pheno.vals[just.m1]) - baseline.val

	just.m2 <- intersect(m1.min.locale, m2.max.locale)
	m2.coef <- mean(pheno.vals[just.m2]) - baseline.val
	
	m1.m2 <- intersect(m1.max.locale, m2.max.locale)
	double.coef <- mean(pheno.vals[m1.m2]) - baseline.val
		
	return(c("m1.name" = m1.name, "m2.name" = m2.name, "m1.effect" = m1.coef, "m2.effect" = m2.coef, "double.effect" = double.coef))
	}