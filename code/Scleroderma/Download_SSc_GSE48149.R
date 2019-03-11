# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Mon Aug 13 12:06:30 EDT 2018

################################################################
#   Differential expression analysis with limma

Download_SSc_GSE48149 <- function(groupings = NULL, path = "~/Documents/Data/Scleroderma/project_data/"){
	
	library(Biobase)
	library(GEOquery)
	library(limma)
	
	# load series and platform data from GEO
	
	gset <- getGEO("GSE48149", GSEMatrix =TRUE, AnnotGPL=FALSE)
	if (length(gset) > 1) idx <- grep("GPL16221", attr(gset, "names")) else idx <- 1
	gset <- gset[[idx]]
	
	# make proper column names to match toptable 
	fvarLabels(gset) <- make.names(fvarLabels(gset))
	
	# group names for all samples
	pheno.data <- phenoData(gset)
	sample.table <- pheno.data@data

	sample.info <- as.character(sample.table$source_name_ch1)

	sml <- rep("X", length(sample.info))

	if(is.null(groupings)){	#if we are plotting all types separately
		sample.types <- unique(sample.info)
		for(i in 1:length(sample.types)){
			group.locale <- which(sample.info %in% sample.types[i])
			sml[group.locale] <- (i - 1)
			}
		sample.map <- cbind(sample.types, paste0("G", 0:(i-1)))
		if(is.null(group.names)){sample.labels <- sample.types} #assign the group names from data
		}else{ #if we have grouped types
		sample.labels <- names(groupings)
		for(i in 1:length(sample.labels)){
			group.locale <- which(sample.info %in% groupings[[i]])
			# cat(groupings[[i]], length(sample.info[group.locale]), "\n")
			sml[group.locale] <- (i-1)
			}
		sample.map <- cbind(sample.labels, paste0("G", 0:(i-1)))
		}
	
	# eliminate samples marked as "X"
	sel <- which(sml != "X")
	sml <- sml[sel]
	gset <- gset[ ,sel]
	

	#Check to see if we need to log2 transform
	#and perform the transform if yes.
	#it is not performed in this data set because
	#the transform was done prior to uploading
	ex <- exprs(gset)
	qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
	LogC <- (qx[5] > 100) ||
	          (qx[6]-qx[1] > 50 && qx[2] > 0) ||
	          (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	if (LogC) { ex[which(ex <= 0)] <- NaN
	  exprs(gset) <- log2(ex) }
	
	# set up the data and proceed with analysis
	sml <- paste("G", sml, sep="")    # set group names

	#look for group differences
	fl <- as.factor(sml)
	gset$description <- fl
	design <- model.matrix(~ description + 0, gset)
	colnames(design) <- levels(fl)
	fit <- lmFit(gset, design)
	
	group.pairs <- pair.matrix(unique(sml))
	group.diff <- apply(group.pairs, 1, function(x) paste0(x[1], "-", x[2], collapse = ""))
	cont.matrix <- makeContrasts(contrasts = group.diff, levels = design)
	
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2, 0.01)
	tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(ex))
	
	tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","Symbol","Synonyms"))
	
	cat("writing to",paste0(path, "/Top_table_GSE48149.txt"))
	write.table(tT, file= paste0(path, "/Top_table_GSE48149.txt"), row.names=F, sep="\t", quote = FALSE)
	
	u_groups <- sort(unique(sml))
	for(i in 1:length(u_groups)){
		sml[which(sml == u_groups[i])] <- sample.labels[i]
		}
	
	# sapply(unique(sml), function(x) print(length(which(sml == x))))
	
	#write out the expression table
	ex <- exprs(gset)
	colnames(ex) <- paste(colnames(ex), sml, sep = "_")	

	write.table(ex, file = paste0(path, "/Expression_Table_GSE48149.txt"), sep = "\t", quote = FALSE)

	}