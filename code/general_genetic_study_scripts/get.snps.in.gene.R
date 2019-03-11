# pop <- readRDS("~/Documents/Grants/R21/Scleroderma/Results/capeRel_testing/pop.RData")
# geno <- readRDS("~/Documents/Grants/R21/Scleroderma/Results/capeRel_testing/popGeno.RData")
# lig4.snps <- get.snps.in.gene("LIG4", "human", 5e5, 5e5)
# snp.locale <- which(geno$marker.names %in% lig4.snps[,1])



get.snps.in.gene <- function(gene.name, organism = c("mouse", "human"), upstream.buffer = 0, downstream.buffer = 0){
	require(biomaRt)
		if(organism == "mouse"){		
		# lib <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
			lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2017.archive.ensembl.org")
			snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="mmusculus_snp", host = "may2017.archive.ensembl.org")	
		}else{
		# lib <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")	
		lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2017.archive.ensembl.org")
		snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp", host = "may2017.archive.ensembl.org")
		}

	gene.coord <- getBM(c("external_gene_name", "chromosome_name", "start_position", "end_position"), "external_gene_name", values = gene.name, mart = lib)
	gene.snps <- getBM("refsnp_id", "chromosomal_region", paste(gene.coord[1,2], as.numeric(gene.coord[1,3] - upstream.buffer), as.numeric(gene.coord[1,4]+downstream.buffer), sep = ":"), snp.db)
	return(gene.snps)
	
	
}