#This function retrieves variants from the Sanger SNP database
#using the function tabix in Samtools (http://www.htslib.org).
#The name of the ftp file may need to be updated periodically.
#to find the file, go to https://www.sanger.ac.uk/sanger/Mouse_SnpViewer/rel-1505
#in the upper right corner, click on "Mouse FTP site" and allow the finder to 
#open the page.
#This will load a disk image on the Desktop. Click on current_snps, which will 
#take you to the directory with the most recent SNPs and find the file that
#follows this pattern: mgp.v5.merged.snps_all.dbSNP142.vcf.gz
#If the name has changed, change this in the arguments below.
#chr = 4; start.pos = 143023930; end.pos = 148864661
#DBA/2J indels: ftp.file = "ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/strain_specific_vcfs/DBA_2J.mgp.v5.indels.dbSNP142.normed.vcf.gz"
#DBA/2J snps: ftp.file = "ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/strain_specific_vcfs/DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz"
# ftp.file = "ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
# chr <- qtl.info[[1]][,2]; start.pos = qtl.info[[1]][,3]; end.pos = qtl.info[[1]][,4]
#strains = "DBA_2J"
# strains = c("129S1_SvImJ","A_J","DBA_2J","CAST_EiJ","PWK_PhJ","WSB_EiJ","NZO_HlLtJ","NOD_ShiLtJ")
# chr = 12; start.pos = 102406588; end.pos = 102407421
#chr = 8; start.pos = 34965259; end.pos = 34965765
#chr = 10; start.pos = 20148796; end.pos = 20149051

get.variants <- function(chr = 10, start.pos = 10000000, end.pos = 11000000, strains = NULL, intergenic.variants = TRUE, ftp.file = "ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"){
	
	require("Rsamtools")
	
	con = TabixFile(ftp.file)
    hdr = headerTabix(con)
    hdr = strsplit(hdr$header, split = "\t")[[length(hdr$header)]]
    hdr = sub("^#", "", hdr)
    
    gr <- GRanges(chr, IRanges(start.pos, end.pos))
    
    variants = scanTabix(con, param = gr)
    variants = lapply(variants, strsplit, split = "\t")
    num.variants = sapply(variants, length)

    if(any(num.variants == 0)) {
        warning(paste("Some regions returned no variants:", paste(names(num.variants)[num.variants == 
            0], collapse = ",")))
    		}

    rng = which(num.variants > 0)
    if(length(rng) == 0){return("No Variants")}
       
     #=========================================================
     # internal functions
     #=========================================================      
   	filter.genotype <- function(genotype.mat){
   		if(is.null(dim(genotype.mat))){
   			genotype.mat <- matrix(genotype.mat, ncol = 1)
   			}
   		split.geno <- apply(genotype.mat, 2, function(x) strsplit(x, ":"))
   		geno.mat <- matrix(NA, ncol = length(split.geno), nrow = length(split.geno[[1]]))
   		for(i in 1:length(split.geno)){
	   		geno.mat[,i] <- sapply(split.geno[[i]], function(x) x[1])
   			}
   		hom.geno.locale <- apply(geno.mat, 2, function(x) c(which(x == "1/1"), which(x == "2/2")))
   		hom.geno.idx <- Reduce("union", hom.geno.locale)   		
   		return(list(geno.mat, hom.geno.idx))
   		}
   
     parse.snp.info <- function(info){
     	split.info <- unlist(strsplit(info, "|", fixed = TRUE))
     	gene.locale <- grep("ENSMUSG", split.info)
     	genes <- split.info[gene.locale]
     	snp.csq <- split.info[gene.locale+3]
     	# strand.info <- split.info[gene.locale+11]
     	# csq.table <- unique(cbind(genes, snp.csq, strand.info))
  	    csq.table <- unique(cbind(genes, snp.csq))
  	    csq.v <- paste(as.vector(t(csq.table)), collapse = ":")
     	return(csq.v)
     	}
     #=========================================================       
     
     
     
    if(length(rng) > 0){
	    	strain.columns = which(hdr %in% strains)
	    	keep = c(c(1:5,8), strain.columns)
	        hdr = hdr[keep]
	        for(i in rng) {
	            variants[[i]] = lapply(variants[[i]], function(a) a[keep])
	            # variants[[i]][variants[[i]] == "N"] = NA
	            variants[[i]] = matrix(unlist(variants[[i]]), nrow = length(variants[[i]]), ncol = length(variants[[i]][[1]]), byrow = T)
	            colnames(variants[[i]]) = hdr
	        	}
		}

	var.calls <- lapply(variants, function(x) x[,7:ncol(x),drop=FALSE])
	good.calls <- lapply(var.calls, filter.genotype)
	
	
	sub.var <- vector(mode = "list", length = length(variants))

	
	for(i in 1:length(variants)){
		var.mat <- variants[[i]]
		var.mat[,7:ncol(var.mat)] <- good.calls[[i]][[1]]
		sub.mat <- var.mat[good.calls[[i]][[2]],,drop=FALSE]
		sub.mat[which(sub.mat == "0/0")] <- "-"
		alt.ind <- which(sub.mat == "1/1", arr.ind = TRUE)
		u_col <- unique(alt.ind[,2])
		for(j in 1:length(u_col)){
			col.locale <- which(alt.ind[,2,drop=FALSE] == u_col[j])
			row.idx <- alt.ind[col.locale,1,drop=FALSE]
			sub.mat[row.idx,u_col[j]]  <- sub.mat[row.idx,"ALT"]
			}
		sub.var[[i]] <- sub.mat
		}
	
	num.check <- sum(sapply(sub.var, nrow))
	if(num.check == 0){return("No Variants")}
	
	var.info <- lapply(sub.var, function(x) strsplit(x[,"INFO"], ";"))
	info.det <- lapply(var.info, function(x) sapply(x, function(y) parse.snp.info(y)))

	final.variants <- vector(mode = "list", length = length(variants))
	for(i in 1:length(variants)){
		trimmed.var <- sub.var[[i]][,-6,drop=FALSE]
		info.var <- cbind(trimmed.var, info.det[[i]])
		colnames(info.var)[ncol(info.var)] <- "Csq"
		if(!intergenic.variants){
			info.var <- info.var[which(info.var[,"Csq"] != ""),,drop=FALSE]
			}
		final.variants[[i]] <- info.var
		}

	num.check <- sapply(final.variants, nrow)
	if(sum(num.check) == 0){return("No Variants")}

    if(length(variants) == 1) {
        final.variants = final.variants[[1]]
    		}
    		
    	

    return(final.variants)


	
}