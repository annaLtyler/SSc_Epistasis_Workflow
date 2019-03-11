#SSc analyses for paper
setwd("~/Documents/Data/Scleroderma/Results/Dominant_2ET_filtered_SNPs")

library(biomaRt)
hum = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
snp.db = useEnsembl(biomart="snp", dataset="hsapiens_snp")

hum <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2017.archive.ensembl.org")
snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp", host = "may2017.archive.ensembl.org")


cross <- readRDS("cross.RData")
geno <- readRDS('~/Documents/Data/Scleroderma/GWAS_documents/geno/geno.RData')

useful.scripts <- list.files("~/Documents/git_repositories/useful_r_code", pattern = ".R", full.names = TRUE)
for(i in 1:length(useful.scripts)){source(useful.scripts[i])}

cape.fun <- list.files("~/Documents/git_repositories/capempp", pattern = ".R", full.names = TRUE)
for(i in 1:length(cape.fun)){source(cape.fun[i])}

#=====================================================================
#replication of previously known SNPs in singlescan
#=====================================================================
	gene.fun <- list.files("~/Documents/git_repositories/general_genetic_study_scripts", ".R", full.names = TRUE)
	for(i in 1:length(gene.fun)){source(gene.fun[i])}
	
	cross.singlescan <- readRDS('cross.singlescan.raw.RData')
	p.thresh <- cross.singlescan$alpha.thresh
	lower.thresh <- p.thresh$"0.05"
	sig.locale <- apply(cross.singlescan$singlescan.t.stats[,,1], 2, function(x) which(abs(x) > lower.thresh))
	sig.snps <- unique(unlist(lapply(sig.locale, names)))
	
	prev.results <- read.table("~/Documents/Data/Scleroderma/Previous_Results/SSc_published_SNP_associations.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
	
	rr.results <- sig.locale
	for(i in 1:length(sig.locale)){
		if(length(sig.locale[[i]]) > 0){
			phenotype <- names(sig.locale[i])
			outcomeV <- cross$pheno[,phenotype]
			for(j in 1:length(sig.locale[[i]])){
				exposureV <- geno[,2,sig.locale[[i]][j]]
				rr.results[[i]][j] <- risk.ratio(exposureV, outcomeV, 0.5)
				}
			}
		}
	
	for(i in 1:length(rr.results)){
		if(length(rr.results[[i]]) > 0){
			quartz();barplot(rr.results[[i]], las = 2, main = names(rr.results)[i]);abline(h = 1)
			}
		}
	
	cbind(unlist(sig.locale), unlist(rr.results))
	rep.snps <- intersect(prev.results[,"SNPS"], sig.snps)
	rep.locale <- match(rep.snps, prev.results[,"SNPS"])
	prev.results[rep.locale,]
	
	snp.table <- getBM(c("refsnp_id", "chr_name", "chrom_start", "minor_allele", "associated_gene"), filters = "snp_filter", values = unlist(lapply(sig.locale, names)), mart = snp.db)
	num.chr <- which(!is.na(as.numeric(snp.table[,2])))
	snp.table <- snp.table[num.chr,]
	snp.info <- condense.table(snp.table, condense.by = 1, col.to.collapse = c(2,3,4), col.to.concat = 5)
	
	#add the phenotype information
	snp.rr <- rep(NA, length(sig.snps))
	snp.pheno <- rep(NA, length(sig.snps))
	for(i in 1:length(rr.results)){
		for(j in 1:length(rr.results[[i]])){
			snp.locale <- which(snp.info[,1] == names(rr.results[[i]])[j])
			snp.rr[snp.locale] <- rr.results[[i]][j]
			snp.pheno[snp.locale] <- names(rr.results)[i]
			}
		}
	
	
	rep.locale <- match(rep.snps, snp.info[,1])
	snp.info[rep.locale,1] <- paste0(snp.info[rep.locale,1], "*")
	
	final.table <- cbind(snp.info, snp.rr, snp.pheno)
	ordered.table <- final.table[order(as.numeric(final.table[,"chrom_start"])),]
	ordered.table[,"snp.rr"] <- signif(as.numeric(ordered.table[,"snp.rr"]), 2)
	write.table(ordered.table, "Singlescan.SNPs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
	
#=====================================================================
# replication of previously known SNPs from CAPE pipeline
# We did not replicate any SNPs with the additional cape pipeline
# but none of the previously identified SNPs were in the pairscan
# there are 22 SNPs on the chip that have previously been associated
# with SSc, but none of them made it into the pair scan. This is because
# I was averaging over too many SNPs to sample SNPS within peaks
# I fixed this in the code, and now I'm re-running everything. Ugh.
#=====================================================================
	
	var.inf <- writeVariantInfluences(cross, 0.05, TRUE, write.file = FALSE)
	not.main.locale <- grep("rs", var.inf[,"Target"])
	just.main <- var.inf[-not.main.locale,]
	just.main.snps <- unique(gsub("_B", "", just.main[,1]))
	rep.snps <- intersect(prev.results[,"SNPS"], just.main.snps)
	
	et.singlescan <- readRDS("cross.singlescan.et.RData")
	et.t.stats <- et.singlescan$singlescan.t.stats[,,1]
	raw.t.stats <- cross.singlescan$singlescan.t.stats[,,1] 
	# cbind(rownames(et.t.stats), rownames(raw.t.stats))
	
	jpeg("ET.v.Traits.jpg", height = 8, width = 12, units = "in", res = 300)
	par(mfrow = c(2,4))
	for(i in 1:ncol(et.t.stats)){
		for(j in 1:ncol(raw.t.stats)){
		plot(abs(et.t.stats[,i]), abs(raw.t.stats[,j]), main = paste(colnames(et.t.stats)[i], colnames(raw.t.stats)[j]), xlab = colnames(et.t.stats)[i], ylab = colnames(raw.t.stats)[j])
		}
	}
	dev.off()
	
	#The effects of ET1 correlate positively with ANA
	#The effects of ET2 correlate positively with centromere
	#This is definitely what we expect
	

	

#=====================================================================
# plot the effects of an interaction on a phenotype
#=====================================================================
	shiny.code <- list.files("~/Documents/R/shiny/capeDO_interactions/code", pattern = ".R", full.names = TRUE)
	for(i in 1:length(shiny.code)){source(shiny.code[i])}
	
	snp1 <- "rs2301156"
	snp2 <- "rs204991"
	phenotype <- "centromere"
	
	snp1.geno <- cross$geno.for.pairscan[,paste0(snp1, "_B")]
	snp2.geno <- cross$geno.for.pairscan[,paste0(snp2, "_B")]
	pheno.vals <- cross$pheno[,phenotype]
	
	pdf("Centromere.Effect.pdf")
	plot.bars(cbind(snp1.geno, snp2.geno), marker.name = c(snp1, snp2), pheno.name = phenotype, pheno.vals, ref.centered = TRUE, error.bars = "se")
	dev.off()


#=====================================================================
# get risk ratios for SNPs and interactions
#=====================================================================
	gen.code <- list.files("~/Documents/git_repositories/general_genetic_study_scripts", pattern = ".R", full.names = TRUE)
	for(i in 1:length(gen.code)){source(gen.code[i])}

	snp1 <- "rs2301156"
	snp2 <- "rs204991"
	phenotype <- "centromere"
	
	snp.rr <- get.snp.risk.ratio(cross, snps = paste0(c(snp1, snp2), "_B"), phenotype = phenotype, trait.type = "raw", geno.coding = "Dominant")

	snp.rr.int <- get.snp.int.risk.ratio(data.obj = cross, snp1 = paste0(snp1, "_B"), snp2 = paste0(snp2, "_B"), phenotype = phenotype, covar = NULL, trait.type = "raw", geno.coding = "Additive", outcome.threshold = 0.5)
	
	#proportion of the population with ACA for each genotype
	prop.aca <- risk.interactions(cross$pheno[,phenotype], cross$geno.for.pairscan[,paste0(snp1, "_B")], cross$geno.for.pairscan[,paste0(snp2, "_B")], plot.results = TRUE)
	barplot(prop.aca[1:4] - prop.aca[1])