colnames(cross$pheno)

limited <- main.effect.patterns(cross, effect.pattern = c("+", "-", "0", "0"), collapsed.net = FALSE)
diffuse <- main.effect.patterns(cross, effect.pattern = c("-", "+", "0", "0"), collapsed.net = FALSE)

limited.snps <- gsub("_B", "", unlist(cross$linkage.blocks.full[rownames(limited)]))
diffuse.snps <- gsub("_B", "", unlist(cross$linkage.blocks.full[rownames(diffuse)]))

limited.rr.topo <- get.snp.risk.ratio(cross, geno, limited.snps, "topo", covar = cross$p.covar, trait.type = "Normalized")
quartz();limited.rr.cent <- get.snp.risk.ratio(cross, geno, limited.snps, "centromere", covar = cross$p.covar, trait.type = "Normalized")
max.val <- max(c(limited.rr.topo, limited.rr.cent))
plot(limited.rr.topo, limited.rr.cent, xlim = c(0,max.val), ylim = c(0,max.val), xlab = "Risk Ratio for Topo", ylab = "Risk Ratio for Centromere")
abline(h = 1); abline(v = 1)

diffuse.rr.topo <- get.snp.risk.ratio(cross, geno, diffuse.snps, "topo", covar = cross$p.covar, trait.type = "Normalized")
diffuse.rr.cent <- get.snp.risk.ratio(cross, geno, diffuse.snps, "centromere", covar = cross$p.covar, trait.type = "Normalized")

max.val <- max(c(diffuse.rr.topo, diffuse.rr.cent))
plot(diffuse.rr.topo, diffuse.rr.cent, xlim = c(0,max.val), ylim = c(0,max.val), xlab = "Risk Ratio for Topo", ylab = "Risk Ratio for Centromere")
abline(h = 1); abline(v = 1)


max.val <- max(c(limited.rr.topo, diffuse.rr.cent))
plot(diffuse.rr.topo, diffuse.rr.cent, xlim = c(0,max.val), ylim = c(0,max.val), xlab = "Risk Ratio for Topo", ylab = "Risk Ratio for Centromere", col = "blue", pch = 16)
points(limited.rr.topo, limited.rr.cent, col = "red", pch = 16)
abline(h = 1); abline(v = 1)
legend("topright", legend = c("Diffuse", "Limited"), pch = 16, col = c("blue", "red"))


