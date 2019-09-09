---
title: "Genetic interactions in Systemic Sclerosis"
author: Anna L Tyler
date: 8/Aug/2018
output: 
  html_document:
    keep_md: true
    code_folding: hide
    collapsed: no
    toc: yes
    toc_float: yes
bibliography: SSc.bib
---





## Introduction
This workflow accompanies the paper ``Genetic interactions affect lung function in 
patients with systemic sclerosis."

This purpose of this project was to investigate genetic interactions influencing the 
presence of autoantibodies as well as lung function in a cohort of patients with systemic 
sclerosis (SSc). We use the Combines Analysis of Pleiotropy and Epistasis (CAPE) \cite{cape}, 
which increases power to detect interactions by combining information across multiple traits.

We identify three high-confidence candidates that influence lung function in SSc patients.

This workflow is a companion to [cite paper] and performs each step used in the analysis.

## CAPE setup
Set up all the libraries and functions we will use


```r
needed.packages <- c("evd", "Matrix", "fdrtool", "shape", "corpcor", "RColorBrewer", 
"doParallel","foreach", "caTools", "stringr", "abind", "propagate", "knitr", "pander", "igraph",
"RColorBrewer", "gProfileR", "pheatmap", "here", "biomaRt")

for(i in 1:length(needed.packages)){library(needed.packages[i], character.only = TRUE)}
```

Set up directories.


```r
project.dir <- here()
code.dir <- here("code")
data.dir <- here("data")
results.dir <- here("results")
```

Source all the code used


```r
all.code.dir <- list.files(code.dir, full.names = TRUE)
for(i in 1:length(all.code.dir)){
	all.fun <- list.files(all.code.dir[i], pattern = ".R", full.names = TRUE)
	for(i in 1:length(all.fun)){source(all.fun[i])}
	}
```



```r
hum = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
```


The following code reads in raw data and builds the cape objects. It is not run here.

```r
build.ssc.cape.data()
```

Instead we will read in the built object carrying the phenotype data,


```r
SScPop <- readRDS(here("data", "project_data", "cross.RData"))
```

and the genotype data. We also converted the genotypes to a dominant coding, which 
gives us more patients for all marker pair genotypes.


```r
all.var <- ls()
geno.loaded <- as.logical(length(which(all.var == "geno")))
if(!geno.loaded){
	geno <- readRDS(here("data", "project_data", "geno.RData"))
	geno[which(geno >= 0.5)] <- 1
}
```

The figure below shows the first three principle components of the genotype matrix plotted
against each other in pairs. There is no evidence of population structure in these plots.
We thus did not correct for population structure.


We used sex as a covariate. 


```r
SScPop <- pheno2covar(SScPop, "sex")
```

## Look at trait distributions
We then selected the traits that we are interested in: 


```r
SScPop <- select.pheno(SScPop, c("centromere", "nucleolar", "pf1_fvcp", "pf1_dlcp"))
```

We looked at two patterns of autoantibody staining. The first is a centromeric pattern 
("centromere"). This pattern is associated with sclerodactyly and Raynaud's phenomenon 
\cite{centromereTitre} as well as other CREST symptomes \cite{centromereCREST}. The 
centromeric pattern of autoantibody is also associated with the limited form of SSc, 
and with patients at lower risk of internal organ involvement \cite{SScAutoAA}. The 
second autoantibody staining pattern is nucleolar, which is associated with interstitial 
lung disease and pulmonary hypertension \cite{nucleolarLung}. The distributions among the
patients were as follows:


```r
d.t <- table(as.data.frame(SScPop$pheno[,c("centromere", "nucleolar")]))
d.t
```

```
##           nucleolar
## centromere   0   1
##          0 379 189
##          1 252   5
```

```r
aa.v <- matrix(as.vector(d.t))
aa.mat <- aa.v/sum(aa.v)
aa.fig.width = g3.fig.min.width
aa.fig.height <- aa.fig.width/5.5
pdf(here("results", "Figures", "Fig_AA_proportions.pdf"), width = aa.fig.width, height = aa.fig.height)
par(mar = c(1.5,1,1.5,1), ps = 9, mgp = c(3, 0.3, 0))
barplot(aa.mat, col = c("#7fc97f", "#beaed4", "#fdc086", "#386cb0"), horiz = TRUE, axes = FALSE)
axis(1)
par(xpd = TRUE)
y.val = 1.7
text(y = y.val, x = 0.25, labels = paste0("None (", aa.v[1], ")"))
text(y = y.val, x = 0.6, labels = paste0("ACA (", aa.v[2], ")"))
text(y = y.val, x = 0.89, labels = paste0("ANA (", aa.v[3], ")"))
text(y = y.val, x = 0.99, labels = paste0("Both\n(", aa.v[4], ")"), srt = 45, adj = 0)
dev.off()
```

```
## quartz 
##      2
```

We also looked at two lung function tests: forced vital capacity percent predicted: 
"pf1\_fvcp", and diffusion lung capacity predicted: "pf1\_dlcp." Forced vital capacity
is the the volume of air that can be forcibly blown out after a full inspiration. The percent
predicted is the percent of the predicted value for a person with similar characteristics, such
as height, age, and sex. 

Diffusion lung capacity is a measure of the efficiency of gas transfer from the lungs to the
blood. As with FVCP, the percentage is based on the percentage of the predicted value for a
person with similar characteristics such as height, age, and sex.

Here we look at some properties of the traits.


```r
#assign variable names to all traits for easier use
centromere <- SScPop$pheno[,"centromere"]
centromere[which(centromere == 0)] <- "-"
centromere[which(centromere == "1")] <- "C"

nucleolar <- SScPop$pheno[,"nucleolar"]
nucleolar[which(nucleolar == 0)] <- "-"
nucleolar[which(nucleolar == "1")] <- "N"

sex <- SScPop$p.covar.table[,1]
sex[which(sex == 0)] <- "F"
sex[which(sex == "1")] <- "M"

FVCP <- SScPop$pheno[,"pf1_fvcp"]
DLCP <- SScPop$pheno[,"pf1_dlcp"]
```

Both traits were roughly normally distributed with maximum values ranging to above 100\%:


```r
par(mfrow = c(1,3))
hist(FVCP, main = "Distribution of FVCP", xlim = c(0, max(FVCP, na.rm = TRUE)), xlab = "FVCP")
hist(DLCP, main = "Distribution of DLCP", xlim = c(0, max(DLCP, na.rm = TRUE)), xlab = "DLCP")
correlation <- cor.test(FVCP, DLCP)
r <- signif(correlation$estimate, 2)
p <- signif(correlation$p.value, 2)
plot(FVCP, DLCP, pch = 16, main = paste0("r = ", r, ", p = ", p))
```

![](SSc_files/figure-html/trait_dist-1.png)<!-- -->


```
## quartz 
##      2
```


We normalized the traits for our analysis using a rank Z transform:


```r
SScPop <- norm.pheno(SScPop)

FVCP <- SScPop$pheno[,"pf1_fvcp"]
DLCP <- SScPop$pheno[,"pf1_dlcp"]

par(mfrow = c(1,3))
hist(FVCP, main = "Distribution of FVCP", xlab = "FVCP")
hist(DLCP, main = "Distribution of DLCP", xlab = "DLCP")
correlation <- cor.test(FVCP, DLCP)
r <- signif(correlation$estimate, 2)
p <- signif(correlation$p.value, 2)
plot(FVCP, DLCP, pch = 16, main = paste0("r = ", r, ", p = ", p))
```

![](SSc_files/figure-html/norm_dist-1.png)<!-- -->

We looked to see if there were any associations between the autoantibody patterns and the 
lung function tests.


```r
#function to make boxplot and test association
plot.association <- function(aa, lung_trait){
	aa.name <- deparse(substitute(aa))
	trait.name <- deparse(substitute(lung_trait))
	aa <- as.factor(aa)
	test <- t.test(lung_trait~aa)
	pval <- signif(test$p.value, 2)
	boxplot(lung_trait~aa, names = levels(aa), ylab = trait.name, 
	main = paste(trait.name, "by", aa.name, "\np =", pval), cex.lab = 1.5, cex.axis = 1.5)
	return(pval)
	}

par(mfrow = c(2,2))
fcp <- plot.association(centromere, FVCP)
dcp <- plot.association(centromere, DLCP)
fnp <- plot.association(nucleolar, FVCP)
dnp <- plot.association(nucleolar, DLCP)
```

![](SSc_files/figure-html/aa_lung-1.png)<!-- -->

The presence of centromeric autoantibody staining has previously been associated with 
lower risk of internal organ involvement in SSc \cite{centromereTitre}, while nucleolar
autoantibody staining has been associated with increased risk of interstitial lung disease
and pulmonary hypertension. In support of these findings, the presence of centromeric autoantibody 
staining in this population was associated with higher FVCP, while the presence of the nucleolar
autoantibody was associated with a slight decrease in FVCP. However, there was no association 
between either autoantibody and DLCP. 

We also looked at the effects of sex:


```r
par(mfrow = c(1,2))
fsp <- plot.association(sex, FVCP)
dsp <- plot.association(sex, DLCP)
```

![](SSc_files/figure-html/sex_lung-1.png)<!-- -->


```
## quartz 
##      2
```


Sex had a main effect on FVCP, but no DLCP. Males had a significantly lower FVCP than females.

In addition for main effects, we looked for interaction effects between the autoantibodies, and 
between the autoantibodies and sex.


```r
factor.interaction <- function(trait, factor1, factor2, plot.type = c("box", "line")){
	plot.type = plot.type[1]
	f1.name <- deparse(substitute(factor1))
	f2.name <- deparse(substitute(factor2))
	trait.name <- deparse(substitute(trait))
	if(plot.type == "box"){
		test <- anova(aov(trait~as.factor(factor1)*as.factor(factor2)))
		pval <- signif(test$"Pr(>F)"[3], 2)
		stats <- boxplot(trait~as.factor(factor1)*as.factor(factor2), ylab = trait.name, 
		main = paste(trait.name, "by", f1.name, "and", f2.name, "\np =", pval), cex.lab = 1.5, cex.axis = 1.5)
		par(xpd = TRUE)
		text(x = 1:4, y = min(trait, na.rm = TRUE)*1.45, labels = stats$n)
		par(xpd = FALSE)
		}else{
		test <- lm(trait~as.factor(factor1)*as.factor(factor2))
		results <- summary(test)$coefficients 
		pval <- signif(results[nrow(results),ncol(results)], 2)
		no.na <- which(!is.na(trait))
		interaction.plot(factor1[no.na], factor2[no.na], trait[no.na], xlab = f1.name, ylab = trait.name, trace.label = f2.name, main = paste("p =", pval))
		}
		invisible(pval)
	}
```

A plot of both autoantibodies together shows that there is no interaction between the effects
of the autoantibodies, although, patients with only the centromere autoantibody trend toward the
highest lung function.


```r
par(mfrow = c(1,2))
fcn.p <- factor.interaction(FVCP, centromere, nucleolar)
dcn.p <- factor.interaction(DLCP, centromere, nucleolar)
```

![](SSc_files/figure-html/aa_int-1.png)<!-- -->

There is an interaction between sex and the nucleolar autoantibody. While patients positive for
the nucleolar autoantibody tend to have a slightly lower FVCP on average, this effect is driven
primarily by the female patients. Male patients with the nucleolar autoantibody tend to have 
slightly higher FVCP than males without the nucleolar autoantibody.


```r
par(mfrow = c(2,2))
factor.interaction(FVCP, centromere, sex)
factor.interaction(FVCP, nucleolar, sex)
fcs.p <- factor.interaction(FVCP, centromere, sex, "line")
fns.p <- factor.interaction(FVCP, nucleolar, sex, "line")
```

![](SSc_files/figure-html/sex_int-1.png)<!-- -->

There was no interaction between sex and the autoantibodies affecting DLCP, although males
positive for each autoantibody trended toward higher DLCP. DLCP in females was largely unaffected
by the presence of the autoantibodies.


```r
par(mfrow = c(2,2))
factor.interaction(DLCP, sex, nucleolar)
factor.interaction(DLCP, sex, centromere)
dns.p <- factor.interaction(DLCP, nucleolar, sex, "line")
dcs.p <- factor.interaction(DLCP, centromere, sex, "line")
```

![](SSc_files/figure-html/sex_int_dlcp-1.png)<!-- -->


```r
int.lines <- function(y,x,split.by, lwd = 3, pval, errors = TRUE, yaxis = TRUE, 
names.cex = 1){
	split.levels <- levels(as.factor(split.by))
	x.levels <- levels(as.factor(x))
	se.mat <- mean.mat <- matrix(NA, nrow = length(split.levels), ncol = length(x.levels))
	rownames(mean.mat) <- rownames(se.mat) <- split.levels
	colnames(mean.mat) <- colnames(se.mat) <- x.levels
	for(i in 1:length(split.levels)){
		for(j in 1:length(x.levels)){
		idx <- intersect(which(x == x.levels[j]), which(split.by == split.levels[i]))
		mean.mat[i,j] <- mean(y[idx], na.rm = TRUE)
		se.mat[i,j] <- sd(y[idx], na.rm = TRUE)/sqrt(length(idx))
		}
	}
	if(errors){
		ymin <- min(mean.mat-se.mat);ymax <- max(mean.mat+se.mat)
	}else{
		ymin <- min(mean.mat);ymax <- max(mean.mat)
	}
	
	plot.height = ymax - ymin
	plot.new()
	plot.window(xlim = c(0.5, ncol(mean.mat)+0.5), ylim = c(ymin-(plot.height*0.05), ymax+(plot.height*0.05)))
	for(i in 1:nrow(mean.mat)){
		points(x = 1:ncol(mean.mat), y = mean.mat[i,], lty = i, type = "l", lwd = lwd)
		if(errors){
		segments(x0 = 1:ncol(mean.mat), y0 = c(mean.mat[i,]+se.mat[i,]), y1 = c(mean.mat[i,]-se.mat[i,]), lwd = lwd)
		}
	}
	par(xpd = TRUE)
	for(i in 1:ncol(mean.mat)){
		if(colnames(mean.mat)[i] == "-"){
			text(x = i, y = (ymin-(plot.height*0.1)), labels = expression(symbol("\306")), cex = names.cex)
		}else{
			text(x = i, y = (ymin-(plot.height*0.1)), labels = colnames(mean.mat)[i], cex = names.cex)
		}
	}
	if(yaxis){axis(2, cex.axis = cex.lab)}
	#text(x = 0, y = mean(c(ymin, ymax)), labels = deparse(substitute(y)), cex = cex.lab, srt = 90)
	text(x = 1.5, y = ymax, labels = paste("p =", signif(pval)))
}

make.legend <- function(split.by, cex = 2, lwd = 2, x = 0.2, y = 0.9){
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,1))
    split.levels <- levels(as.factor(split.by))
    split.levels[which(split.levels == "-")] <- expression(symbol("\306"))
    legend(x,y,lty = 1:length(split.levels), legend = split.levels, cex = cex, lwd = lwd, horiz = TRUE, bty = "n")
	# draw.rectangle(0,1,0,1)
}


pdf(here("results", "Figures", "Fig_Lung_Interactions_All.pdf"), width = 10, height = 6)
layout.mat <- matrix(c(0, 1:3, 0, 4:14), nrow = 4, byrow = TRUE)
layout(layout.mat, widths = c(0.2, rep(1, 3)), heights = c(0.3, 0.3, rep(1, 3)))

par(mar = c(1, 1, 1, 1))
plot.text("Centromere by Nucleolar", cex = cex.lab)
plot.text("Centromere by Sex", cex = cex.lab)
plot.text("Nucleolar by Sex", cex = cex.lab)

par(mar = c(0, 0, 0, 0))
make.legend(nucleolar, lwd = 1.5, cex = 1.5)
make.legend(sex, lwd = 1.5, cex = 1.5)
make.legend(sex, lwd = 1.5, cex = 1.5)

plot.text("FVCP", srt = 90, cex = cex.lab)
par(mar = c(2,3,0,0))
int.lines(y = FVCP, x = centromere, split.by = nucleolar, pval = fcn.p, errors = TRUE)
int.lines(FVCP, centromere, sex, pval = fcs.p, errors = TRUE, yaxis = FALSE)
int.lines(FVCP, nucleolar, sex, pval = fns.p, errors = TRUE, yaxis = FALSE)

par(mar = c(1, 1, 1, 1))
plot.text("DLCP", srt = 90, cex = cex.lab)
par(mar = c(2,3,0,0))
int.lines(DLCP, centromere, nucleolar, pval = dcn.p)
int.lines(DLCP, centromere, sex, pval = dcs.p)
int.lines(DLCP, nucleolar, sex, pval = dns.p)
dev.off()
```

```
## quartz 
##      2
```

```r
#also plot just the significant and almost significant interactions
fig.int.height <- g3.fig.min.width
```

![](SSc_files/figure-html/int_fig-1.png)<!-- -->

```r
fig.int.width <- fig.int.height*0.45
pdf(here("results", "Figures", "Fig_Lung_Trait_Interactions_Sig.pdf"), width = fig.int.height, height = fig.int.width)
# quartz(width = fig.int.height, height = fig.int.width)
layout.mat <- matrix(c(1:6), nrow = 2, byrow = TRUE)
layout(layout.mat, widths = c(1,1), heights = c(0.3, 1))
# layout.show(9)
legend.cex <- 1;legend.lwd = 1
par(mar = c(0, 0, 0, 0), ps = 9)
legend.x <- 0.25;legend.y <- 0.8
make.legend(sex, lwd = legend.lwd, cex = legend.cex, x = legend.x, y = legend.y)
make.legend(sex, lwd = legend.lwd, cex = legend.cex, x = legend.x, y = legend.y)
make.legend(sex, lwd = legend.lwd, cex = legend.cex, x = legend.x, y = legend.y)

par(mar = c(4,3,0,0));ylab.line = 1.8; xlab.line = 1.3
int.lines(FVCP, nucleolar, sex, pval = fns.p, errors = TRUE, yaxis = TRUE, lwd = 1, names.cex = 1.5)
mtext(side = 2, "FVCP", line = ylab.line);mtext(side = 1, "Nucleolar", line = xlab.line)
int.lines(DLCP, centromere, sex, pval = dcs.p, lwd = 1, names.cex = 1.5)
mtext(side = 2, "DLCP", line = ylab.line);mtext(side = 1, "Centromere", line = xlab.line)
int.lines(DLCP, nucleolar, sex, pval = dns.p, lwd = 1, names.cex = 1.5)
mtext(side = 2, "DLCP", line = ylab.line);mtext(side = 1, "Nucleolar", line = xlab.line)
dev.off()
```

```
## quartz 
##      2
```

## Single marker effects
We first tested the main effect of each SNP with a linear model. There were significant
main effects on 


## Filter SNPs
Because CAPE cannot test all pairs of markers exhaustively, we need to filter the markers 
to those that are most likely to have interactions in this data set. CAPE uses information 
from pleiotropic markers to identify epistasis. We thus selected markers based on both 
epistasis and pleiotropy criteria. We used MatrixEpistasis [@Zhu:dn], which is an ultra-fast 
matrix-based method that can test all the pairs in the data set exhaustively. We set the p 
value filter at 1e-6 and use sex as a covariate (specified above).

The next function runs MatrixEpistasis and asks for pairs with p values below a given cutoff. 
It saves the results for each phenotype. After running MatrixEpistasis, it goes through the
phenotype pairs and gets the unique SNPs from the most significant pairs across all traits. 
This function takes a long time to run, so I am not running it here.

The parameter to note is the p value. I used a cutoff of p = 1e-6, so MatrixEpistasis will 
only return pairs that have that p value or less.


```r
selected.SNPs <- filterSNPs(SScPop, geno, target.markers = 1500, ref.allele = "A", 
pval = 1e-6, chunk.size = 5000, n.cores = 4)
```

This function saves a file for each phenotype that includes the top SNP pairs for each 
phenotype based on the p value threshold. By running filterSNPs again, but skipping the 
MatrixEpistasis step, we can take a closer look at these pairs. For example, we can see which 
epistatic pairs are also pleiotropic, meaning they influence more than one phenotype.


```r
filterSNPs(SScPop, geno, run.matrixEpistasis = FALSE, target.markers = 1500)
```

There isn't much overlap among the traits. There are 52 pairs shared between the two lung 
phenotypes.
filterSNPs keeps all SNPs in these pairs that are both epistatic and pleiotropic. Then it 
moves on to the single-marker stats. Single markers can influence multiple traits when paired 
with other markers. Even it the marker pair is different, individual markers might epistatically 
influence multiple traits. flterSNPs finds the SNPs that are the most pleiotropic in the 
following way:

While we have not reached our target number of markers:
Start at N = the number of traits tested.

1) Find all SNPs that influence N traits. 
2) If we can add all of these SNPs to our list, add them. 
3) If this will make the list too big, sort the SNPs by their representation in the traits. 
For example, a SNP that influenced one trait in twelve pairs will be prioritized over a SNP 
that influenced one trait in five pairs. 
4) Add SNPs in order of their influence until we have reached the target number of markers.

We then run cape using this list of filtered SNPs and the following parameters:


```r
param.file <- here("results", "cape.parameters.txt")
param.table <- readParameters(param.file)
kable(param.table[which(!is.na(param.table[,1])),])
```

                          x                                      
------------------------  ---------------------------------------
traits                    centromere nucleolar pf1_fvcp pf1_dlcp 
covariates                sex                                    
traits.scaled             TRUE                                   
traits.normalized         TRUE                                   
scan.what                 Eigentraits                            
eig.which                 1 2                                    
pval.correction           fdr                                    
ref.allele                A                                      
singlescan.perm           0                                      
marker.selection.method   from.list                              
SNPfile                   filteredSNPs.txt                       
max.pair.cor              0.5                                    
pairscan.null.size        1500000                                

Before we run cape, we also want to check for population structure to decide whether we 
should do a correction for population structure. To test for structure, we performed single
marker regression at each marker to look for associations with each trait. Under the null
hypothesis, p values are uniformly distributed. If there is population structure influencing 
the results, we expect to see a substantial deviation from this distribution creating inflated
significance for a substantial number of markers. On the other hand, if the p value distributions
are largely uniform with a few markers deviating from this distribution, there is no population
structure to correct for.

The following figure shows the expected distribution of the -log p values compared to the 
theoretical distribution. There are a few deviations from the expected for both the nucleolar 
and centromere traits indicating that there are a few SNPs with significant main effects.
There is no deviation from the expected distribution at all for the two lung traits. Taken
together, these results indicate that there is no effect of population structure in this 
cohort. Therefore, we did not correct for population structure.



```r
#need to debug this


regression.file <- here("results", "single.regression.RData")
sub.geno <- get.geno(SScPop, geno)
geno.mat <- sub.geno[,2,]

if(!file.exists(regression.file)){
	get.lin.p <- function(y,x){
		num.y <- length(unique(y[which(!is.na(y))]))
		if(num.y == 2){fam = "binomial"}
		if(num.y > 2){fam = "gaussian"}
		y <- y - min(y, na.rm = TRUE)
		y <- y/max(y, na.rm = TRUE)
		model <- glm(y~x, family = fam)
		modelc <- coefficients(summary(model))
		p <- modelc[nrow(modelc), ncol(modelc)]
		return(p)
	}

	all.p <- matrix(NA, ncol = ncol(SScPop$pheno), nrow = ncol(geno.mat))
	colnames(all.p) <- colnames(SScPop$pheno)
	for(i in 1:ncol(SScPop$pheno)){
		cat(i, "\n")
		pvals <- unlist(lapply_pb(1:ncol(geno.mat), function(x) get.lin.p(SScPop$pheno[,i], geno.mat[,x])))
		#pvals <- apply(geno.mat, 2, function(x) get.lin.p(SScPop$pheno[,i], x))
		all.p[,i] <- pvals
	}

	save(all.p, file = regression.file)
}else{
	all.p <- readRDS(regression.file)
}

jpeg(here("results", "Figures", "Supplemental_Figure_Population_Structure.jpg"), 
width = 7, height= 7, units = "in", res = 300)
par(mfrow = c(2,2))
for(i in 1:ncol(all.p)){
	qqnorm(qnorm(all.p[,i]), main = colnames(SScPop$pheno)[i])
	qqline(qnorm(all.p[,i]))
}
dev.off()
```

Here we are running CAPE with two of four eigentraits. These eigentraits contrast the two 
lung function tests with each of the autoantibodies in turn.


```r
eig.pop <- get.eigentraits(SScPop)
```

```
## Removing 417 individuals with missing phenotypes.
```

```r
var.accounted <- plotSVD(eig.pop)
```

![](SSc_files/figure-html/eig.traits-1.png)<!-- -->


```
## quartz 
##      2
```


```r
system(paste("cp filteredSNPs.txt", results.dir))
final.pop <- run.cape(SScPop, geno)
```


```r
pop.data <- here("results", "cross.RData")
SScPop <- readRDS(pop.data)
```

These SNPs filtered for both pleiotropy and epistasis formed a CAPE network.


```r
plotNetworkDO(SScPop, show.alleles = FALSE)
```

![](SSc_files/figure-html/net_plot-1.png)<!-- -->

The table of significant SNP interactions is as follows:


```r
var.inf <- writeVariantInfluences(SScPop, include.main.effects = FALSE, write.file = FALSE)
knitr::kable(as.data.frame(var.inf[,c(1,2,4,5,10,14)]))
```



Source         Chr   Target         Chr   Effect               p.adjusted 
-------------  ----  -------------  ----  -------------------  -----------
rs1449292_B    3     rs792399_B     17    -1.00948265115757    0          
rs1449292_B    3     rs556874_B     3     -0.869869093605696   0          
rs2055566_B    2     rs3746070_B    19    0.818228099493505    0          
rs6574751_B    14    rs3746070_B    19    0.770632024760673    0          
rs7033108_B    9     rs6466555_B    7     1.15341633886639     0          
rs7033108_B    9     rs10245192_B   7     1.15341633886639     0          
rs2065936_B    10    rs2734510_B    11    1.29903849272547     0          
rs17694404_B   17    rs6105966_B    20    1.24791388253724     0          
sex            0     rs7854383_B    9     1.41845393260669     0          
rs7033108_B    9     rs3746070_B    19    0.806855601041902    0          
rs9607931_B    22    rs3746070_B    19    0.815874044242893    0          

We used [SNPNexus](http://www.snp-nexus.org) to learn more about the individual SNPs. 
In particular we identified the nearest (or overlapping) gene for each SNP. We assigned 
each SNP to its nearest gene. The nodes are labeled with the gene nearest to the SNP. 
Node color indicates whether the SNP overlaps the gene, is upstream, or downstream, and 
node shape indicates whether the gene is protein coding, or another type of gene. Edge color
indicates positive (enhancing) and negative (suppressing) interactions.



```r
snp.file <- here("data", "ncsnp_11708.txt")
snp.gene.table <- as.matrix(read.table(snp.file, sep = "\t", stringsAsFactors = FALSE, header = TRUE))
snp.net <- SNPNexus.to.gene.net(SScPop, snp.gene.table, vertex.size = 40, label.x.shift = 0.3,
label.y.shift = 0.1)
```

![](SSc_files/figure-html/snp_nexus-1.png)<!-- -->


```
## quartz 
##      2
```

The two most significant interactions are suppressing and create a three-node subnetwork between
three SNPs that overlap protein coding genes: WNT5A, RBMS3, and MSI2. WNT signaling is known 
to be involved in SSc, and for the remainder of this study, we focus on this top subnetwork. 


```r
top.ints <- var.inf[1:2, c("Source", "Target")]
top.int.genes <- top.ints #replace the SNP names with gene names
top.int.genes[,1] <- snp.gene.table[match(gsub("_B", "", top.ints[,1]), snp.gene.table[,2]), 5]
top.int.genes[,2] <- snp.gene.table[match(gsub("_B", "", top.ints[,2]), snp.gene.table[,2]), 5]
```

The SNPs interacting in this network are the following:


```r
top.snps <- unique(as.vector(top.ints))
cat(top.snps, sep = "\t")
```

```
## rs1449292_B	rs792399_B	rs556874_B
```

```r
top.genes <- unique(as.vector(top.int.genes))
cat(top.genes, sep = "\t")
```

```
## RBMS3	MSI2	WNT5A
```

## Gene Expression

The following data set measured gene expression in SSc and control lungs 
[GSE76808](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76808). 
Samples were obtained from open lung biopy in 28 consecutive patients with 
SSc-related interstitial lung disease and four controls. We looked for 
differential expression of the genes in our epistatic network.

We used the "Analyze with GEO2R" button on the GEO website to get the data downloaded and 
normalized properly. Then we added custom code to look for expression differences in 
individual genes. 

The following functions download the GEO data set, normalize it, and write out the expression 
data as well as a gene ID table in /Users/atyler/Documents/Projects/Scleroderma/paper/data

They then read in the expression data and plot expression differences in the genes identified
in the epistatic network.

We first looked at expression in the lung:


```r
dataset <- "GSE76808" #Association of Interferon- and Transforming Growth Factor-Regulated 
					  #Genes and Macrophage Activation With Systemic Sclerosis-Related 
					  #Progressive Lung Fibrosis
```


```r
geo.fun <- match.fun(paste0("Download_SSc_", dataset))
geo.fun(path = data.dir)
```




The expression values are normalized by patient, and there do not appear to be any 
batch effects in these data.


```r
expr <- read.table(here("data", paste0("Expression_Table_", dataset, ".txt")), stringsAsFactors = FALSE, sep = "\t", header = TRUE)
gene.id.table <- read.table(here("data", paste0("Top_table_", dataset, ".txt")), stringsAsFactors = FALSE, sep = "\t", header = TRUE)
pt.class <- sapply(strsplit(colnames(expr), "_"), function(x) x[2])
pt.cols <- brewer.pal(8, "Pastel1")
pt.col <- rep(pt.cols[1], length(pt.class))
pt.col[which(pt.class == "SSc")] <- pt.cols[2]
boxplot(expr, col = pt.col, las = 2)
```

![](SSc_files/figure-html/view_data-1.png)<!-- -->


```r
deg.expr <- DEG.in.GEO(top.genes, gene.id.table, expr, return.data = "all")
```

![](SSc_files/figure-html/plot_expression_lung-1.png)<!-- -->


```r
fig.width = g3.fig.min.width
fig.height = g3.fig.min.width/2
pdf(here("results", "Figures", "Fig_DEG.pdf"), width = fig.width, height = fig.height)
par(mfrow = c(1,2), mar = c(1,3,2,0))
for(i in 1:length(deg.expr)){
	if(!is.null(deg.expr[[i]])){
		gene.expr <- deg.expr[[i]][[1]]
		ymax <- max(unlist(gene.expr), na.rm = TRUE)
		ymin <- min(unlist(gene.expr), na.rm = TRUE)
		plot.height <- ymax - ymin
		boxplot(gene.expr, axes = FALSE, main = names(deg.expr)[i], cex = 0.7, cex.main = 0.9)
		axis(2, at = signif(segment.region(ymin, ymax, 4, "ends"), 2))
		par(xpd = TRUE)
		text(x = c(1,2), y = (ymin - (plot.height*0.1)), labels = names(gene.expr))
		if(i == 1){
			mtext("Expression (A.U.)", side = 2, line = 2)
		}
		par(xpd = FALSE)
		stripchart(gene.expr, vertical = TRUE, pch = 16, add = TRUE)
	}
}
dev.off()
```

```
## quartz 
##      2
```

Both WNT5A and RBMS3 were upregulated in SSc lung samples. MSI2 was not present in this data set.


## Gene coexpression
Identifying clusters of co-expressed genes can provide clues to the functions of genes with unknown
functions. We used this principle to identify possible functional roles of the two interactions 
in the top subnetwork. We used the two pairs of genes as query genes in the Search-Based Exploration of 
Expression Compendium (SEEK) \cite{seek}. We restricted our search to subsets of experiments with 
potentially informative subjects (Table XXX table showing terms used and number of data sets 
associated with each term). We found all genes that were significantly correlated (p < 0.05) with 
the query genes in each subset of data sets, and looked for functional enrichment among the genes 
using gprofiler \cite{gprofiler}.


```r
#set up objects to hold the gprofiler for both pairs
#name them by the interaction pair
pair.sig.enrichments <- vector(mode = "list", length = nrow(top.int.genes))
names(pair.sig.enrichments) <- apply(top.int.genes, 1, function(x) paste0(x, collapse = "_"))
for(gp in 1:nrow(top.int.genes)){
	pair.dir <- here("data", "SEEK", names(pair.sig.enrichments)[gp])
	pair.sig.enrichments[[gp]] <- SEEK.enrichment(top.int.genes[gp,], pair.dir, pval = 0.05)
	}
saveRDS(pair.sig.enrichments, here("results", "SEEK_enrichment.RDS"))
```

To identify terms that are SSc-related, we used the 
[comparative toxicogenomic database](http://ctdbase.org) has data sets predicting 
associations between diseases and GO terms 
[here](http://ctdbase.org/downloads/#phenotypediseases). We downloaded this data set on 
August 16, 2019.

Diseases are associated with GO terms. We hypothesized that related diseases have 
highly overlapping sets of GO terms. We further hypothesized that if the SEEK co-expressed
genes are related to Scleroderma, that the enriched terms will highly overlap with 
SSc-related GO terms.


```r
all.var <- ls()
bp.loaded <- as.logical(length(which(all.var == "bp.table")))
cc.loaded <- as.logical(length(which(all.var == "cc.table")))
mf.loaded <- as.logical(length(which(all.var == "mf.table")))

if(!bp.loaded){
    bp.table <- read.csv(
        here("data", "CTD", "CTD_Phenotype-Disease_biological_process_associations.csv"), 
        stringsAsFactors = FALSE, comment.char = "#", header = FALSE)
}

if(!cc.loaded){
    cc.table <- read.csv(
        here("data", "CTD", "CTD_Phenotype-Disease_cellular_component_associations.csv"), 
        stringsAsFactors = FALSE, comment.char = "#", header = FALSE)
}

if(!mf.loaded){
    mf.table <- read.csv(
        here("data", "CTD", "CTD_Phenotype-Disease_molecular_function_associations.csv"), 
        stringsAsFactors = FALSE, comment.char = "#", header = FALSE)
}
```

We performed the SEEK analysis separately for two gene-gene interactions:
WNT5A - RBMS3 and MSI2 - RBMS3. Here we looked at each of these results sets 
separately.


```r
seek.results <- readRDS(here("results", "SEEK_enrichment.RDS"))
```


```r
disease.id.table <- unique(rbind(cc.table[,c(3,4)], mf.table[,c(3,4)], bp.table[,c(3,4)]))
```

We created an incidence matrix relating diseases (in rows) to GO terms (in columns).


```r
#create a two-column matrix linking all diseases and GO terms
mat.bp <- cbind(bp.table[,3], bp.table[,2])
mat.cc <- cbind(cc.table[,3], cc.table[,2])
mat.mf <- cbind(mf.table[,3], mf.table[,2])

#include our gene pairs and the terms they are connected to.
pair1.terms <- lapply(seek.results[[1]], function(x) x[,"term.id"])
pair1.mat <- lapply(pair1.terms, function(x) cbind(rep(names(seek.results)[1], length(x)), x))

pair2.terms <- lapply(seek.results[[2]], function(x) x[,"term.id"])
pair2.mat <- lapply(pair2.terms, function(x) cbind(rep(names(seek.results)[2], length(x)), x))

#full.mat <- rbind(mat.bp, mat.cc, mat.mf, pair1.mat[[seek.set]], pair2.mat[[seek.set]])
full.mat <- rbind(mat.bp, mat.cc, mat.mf)

inc.mat.file <- here("results", "Incidence.matrix.RDS")

if(!file.exists(inc.mat.file)){
    inc.mat <- incidence.matrix(full.mat, verbose = interactive)
    saveRDS(inc.mat, inc.mat.file)
}else{
    inc.mat <- readRDS(inc.mat.file)
}
```

## Overlap between GO terms {.tabset .tabset-fade .tabset-pills}

We tested the significance of the overlap for each pair across all GEO scenarios using
a Fisher's Exact Test. The following plots all show the same data in different configurations.
Each cell represents the significance of the overlap between the enriched GO terms for the 
gene pair in the GEO context listed. The values shown in each cell are the p values. The
asterisks indicate which p values are significant at a 0.05/16 significance level.

The RBMS3-WNT5A gene pair GO sets were enriched for SSc-related terms more often than 
the RBMS3-MSI2 gene pair GO sets. The GO terms 


```r
ssc.locale <- grep("Scleroderma", rownames(inc.mat))
ssc.go.locale <- which(colSums(inc.mat[ssc.locale,]) > 0)
ssc.go.ids <- colnames(inc.mat)[ssc.go.locale]

overlap.p.table <- matrix(NA, ncol = length(seek.results), nrow = length(seek.results[[1]]))
colnames(overlap.p.table) <- names(seek.results)
rownames(overlap.p.table) <- names(seek.results[[1]])

all.contingencies <- vector(mode = "list", length = length(seek.results))
names(all.contingencies) <- names(seek.results)

for(i in 1:length(seek.results)){
    contingencies <- vector(mode = "list", length = length(seek.results[[1]]))
    names(contingencies) <- names(seek.results[[1]])
    for(j in 1:length(seek.results[[1]])){
        pair.go.ids <- seek.results[[i]][[j]][,"term.id"]
        test.results <- test.overlap(ssc.go.ids, pair.go.ids, colnames(inc.mat))
        contingencies[[j]] <- test.results[[1]]
        overlap.p.table[j,i] <- test.results[[2]]
    }
    all.contingencies[[i]] <- contingencies
}
```

### Clustered by -log10 p value


```r
bonferroni.p <- 0.05/length(overlap.p.table)
display.sig <- matrix(ifelse(overlap.p.table <= bonferroni.p, "*", ""), 
ncol = ncol(overlap.p.table))

alt.display <- display.sig
for(i in 1:nrow(overlap.p.table)){
    alt.display[i,] <- paste(signif(overlap.p.table[i,], 2), display.sig[i,])
}

pheatmap(-log10(overlap.p.table), display_numbers = alt.display,
fontsize_number = 12, fontsize_col = 12, fontsize_row = 12)
```

### Original order 

```r
#quartz(width = 4, height = 4)
pheatmap(-log10(overlap.p.table), display_numbers = alt.display, 
fontsize_number = 12, cluster_rows = FALSE, cluster_cols = FALSE, 
fontsize_col = 12, fontsize_row = 12)
```

## Proportion SSc-related 
The following figures show overlaps between SSc-related GO terms and the GO terms that 
were enriched when looking at genes co-expressed with the two gene pairs.
Each bar represents a different subset of GEO experiments, identified either by keyword
or tissue type. 


```r
prop.ssc <- lapply(all.contingencies, function(x) lapply(x, function(y) y[,2]))
prop.ssc.mats <- lapply(prop.ssc, function(x) Reduce("cbind", x))

bar.cols <- c("#f0f0f0", "#636363")

par(mfrow = c(1,2))
for(i in 1:length(prop.ssc.mats)){
    barplot(prop.ssc.mats[[i]][,8:1], names = rev(names(seek.results[[1]])), 
    main = names(seek.results)[i], col = bar.cols, horiz = TRUE, las = 2)
    legend("topright", fill = bar.cols, legend = c("Not Associated", "SSc Associated"), 
	cex = 0.5)
}
```

![](SSc_files/figure-html/GO_prop-1.png)<!-- -->

## Enrichment terms associated with SSc {.tabset .tabset-fade .tabset-pills}

The following plots show a more detailed version of the above information. Here the 
$-log_{10}(p)$ are shown for all individual terms that are shared between at least 
one SEEK search and SSc.


```r
term.table <- Reduce("rbind", lapply(seek.results, function(x) Reduce("rbind", lapply(x, function(y) y[,c("term.id", "term.name")]))))
pair.table <- vector(mode = "list", length = length(seek.results))
names(pair.table) <- names(seek.results)
ssc.locale <- which(rownames(inc.mat) == "Scleroderma, Systemic")
for(p in 1:length(seek.results)){
    inc.table <- matrix(0, ncol = length(seek.results[[p]]), nrow = nrow(term.table))
    colnames(inc.table) <- names(seek.results[[p]])
    rownames(inc.table) <- term.table[,1]
    pair.go <- lapply(seek.results[[p]], function(x) x[,"term.id"])
    pair.p <- lapply(seek.results[[p]], function(x) x[,"p.value"])
    
    for(i in 1:length(seek.results[[p]])){
        go.locale.col <- match(pair.go[[i]], colnames(inc.mat))
        ssc.go <- inc.mat[ssc.locale,go.locale.col]
        names(ssc.go) <- pair.go[[i]]
        ssc.go[which(is.na(ssc.go))] <- 0
        go.locale.row <- match(pair.go[[i]], term.table[,1])
        ssc.related <- which(ssc.go == 1)
        seek.locale <- match(names(ssc.go)[ssc.related], seek.results[[p]][[i]][,"term.id"])
        seek.p <- seek.results[[p]][[i]][seek.locale,"p.value"]
        ssc.go[ssc.related] <- seek.p
        inc.table[go.locale.row,i] <- ssc.go
    }
    pair.table[[p]] <- inc.table
}

combined.table <- Reduce("cbind", pair.table)
with.hits <- which(rowSums(combined.table) > 0)
pair.go.table <- combined.table[with.hits,]
    
#translate the IDs to names
term.locale <- match(rownames(pair.go.table), term.table[,1])
rownames(pair.go.table) <- term.table[term.locale,2]
go.ids <- term.table[term.locale,1]
go.names <- term.table[term.locale,2]

#get gene lists for each GO ID
go.gene.file <- here("results", "GO.genes.RData")
if(!file.exists(go.gene.file)){
    go.genes <- lapply(go.ids, function(x) 
    getBM(c("external_gene_name", "ensembl_gene_id"), "go", x, hum))
    names(go.genes) <- go.ids
    saveRDS(go.genes, go.gene.file)
}else{
    go.genes <- readRDS(go.gene.file)
}

#filter GO terms to those with fewer than a set number of genes
max.genes <- 300
n.genes <- sapply(go.genes, nrow)
small.sets <- which(n.genes <= max.genes)

log.p.table <- -log10(pair.go.table[small.sets,])
log.p.table[which(!is.finite(log.p.table))] <- 0
```

### Enriched terms associated with SSc by category

```r
pheatmap(sqrt(log.p.table), cluster_cols = FALSE)
```

### Terms clustered by keyword


```r
go.rds <- here("data", "Ontologies", "GO.RDS")
if(!file.exists(go.rds)){
    go.list <- read.obo(here("data", "Ontologies", "go.obo"))
    saveRDS(go.list, go.rds)
}else{
    go.list <- readRDS(go.rds)
}

parent.list <- lapply(go.ids, function(x) get.parents(x, go.list, "ID"))
names(parent.list) <- rownames(log.p.table)

keyword.list <- list(
    "ECM/\nCytoskeleton" = c("matrix", "cytoskeleton", "extracellular", "organization", "size", "filopodium"),
    "Collagen\nMetabolism" = c("collagen", "cytokine", "development"),
    "Rho protein\nsignal\ntransduction" = c("signal", "surface"),
    "Binding" = c("binding"),
    "Transcription" = c("transcription", "RNA"),
    "Cell Adhesion" = c("adhesion"),
    "Wounding" = c("wounding"),
    "Vascular" = c("vasoconstriction"),
    "Growth Factor" = c("growth factor"))    


grouped.GO <- group.GO.by.keyword(go.terms = rownames(log.p.table), keyword.list)
row.order <- match(unlist(grouped.GO), rownames(log.p.table))
pheatmap(sqrt(log.p.table[row.order,]), cluster_cols = FALSE, cluster_rows = FALSE)
```

```
## quartz 
##      2
```


```r
highlight.terms <- unique(unlist(keyword.list))
```
Supplemental\_Figure\_SEEK\_enrichment.pdf shows the $-log_{10}$ of the p value for each
term in each subset of data sets. Terms containing strings matching the strings above are
colored blue.


```r
pdf(here("results", "Figures", "Supplemental_Figure_SEEK_enrichment.pdf"), width = 11, height = 8)
for(gp in 1:length(pair.sig.enrichments)){
	all.enrich.sig <- pair.sig.enrichments[[gp]]
	for(i in 1:length(all.enrich.sig)){
		all.enrich.sig[[i]] <- plot.enrichment.vis(all.enrich.sig[[i]], 
		num.terms = nrow(all.enrich.sig[[i]]), plot.label = 
		paste("Enrichment of Genes Significantly Correlated with\n", names(pair.sig.enrichments)[gp],
		"in", names(all.enrich.sig)[i]),
		highlight.terms = highlight.terms)
		pair.sig.enrichments[[gp]] <- all.enrich.sig
		}
	}
dev.off()
```

```
## quartz 
##      2
```



## Clinical significance

The top two SNP interactions affect the traits in the following way:


```r
par(mar = c(6,7,6,3))
for(i in 1:nrow(top.ints)){
    plot.effects(SScPop, top.ints[i,1], top.ints[i,2], plot.type = "b", 
	error.bars = "se", marker1.label = top.int.genes[i,1], 
	marker2.label = top.int.genes[i,2], genotype.labels = c("O", "+"))
    }
```

![](SSc_files/figure-html/int_effects-1.png)<!-- -->![](SSc_files/figure-html/int_effects-2.png)<!-- -->


```r
#plot effects as figures as well
fig.width = g3.fig.max.width
fig.height <- fig.width/4
for(i in 1:nrow(top.ints)){
	pdf(here("results", "Figures", paste0("Fig_Effects", i, ".pdf")), 
	width = fig.width, height = fig.height)
	par(mar = c(3,5,2,3))
    plot.effects(SScPop, top.ints[i,1], top.ints[i,2], plot.type = "b", 
    error.bars = "se", marker1.label = top.int.genes[i,1], 
    marker2.label = top.int.genes[i,2], pheno.labels = c("ACA", "ANA", "FVCP", "DLCP"),
	genotype.labels = c("O", "+"), cex.axis = 1.2, cex.text = 1.2)
	dev.off()
    }
```

The interaction between RBMS3 and MSI2 has less than additive effects on both autoantibodies,
and no real effect on the lung traits. The RBMS3 and WNT5A interaction has a less than additive
effect on the presence of the nucleolar autoantibody as well as the two lung traits. Less than
additive effects suggest that the interacting genes operate in a single pathway, and that the
effects of the alternate alleles are redundant. In both interactions, both alternate alleles 
also affect all traits in the same direction, either negatively or positively.

## Positive interaction effects
We also looked at the top positive interaction effect. The most significant positive 
interaction was between rs2055566, which is located in an intron of FARP2, and rs3746070,
which is located 736 bp upstream of S1PR4. 

These two variants had additive effects on all traits except DLCP. Neither variant
had a significant main effect on DLCP, but individuals with both alternate alleles,
had improved DLCP relative to other individuals. 



```r
top.pos <- var.inf[3,c("Source", "Target"),drop=FALSE]
top.pos.genes <- top.pos #replace the SNP names with gene names
top.pos.genes[,1] <- snp.gene.table[match(gsub("_B", "", top.pos[,1]), snp.gene.table[,2]), 5]
top.pos.genes[,2] <- snp.gene.table[match(gsub("_B", "", top.pos[,2]), snp.gene.table[,2]), 8]

   for(i in 1:nrow(top.pos)){
	plot.effects(SScPop, top.pos[i,1], top.pos[i,2], plot.type = "b", 
	error.bars = "se", marker1.label = top.pos.genes[i,1], 
	marker2.label = top.pos.genes[i,2], genotype.labels = c("O", "+"))
}
```

![](SSc_files/figure-html/top_pos-1.png)<!-- -->

S1PR4 encodes the sphingosine-1-phosphate receptor 4. Sphingosine-1-phosphate is known
to affect SSc pathogenesis. It is elevated in the serum of SSc individuals and is capable
of producing many of the abnormalities seen in SSc [@pattanaik2010role,tokumura2009elevated]. 
Additionally, sphingosine-1-Phosphate Receptor 5 Modulates early-stage processes during 
fibrogenesis in a mouse model of SSc [@schmidt2017sphingosine], and Also, the selective 
S1P1 receptor modulator cenerimod attenuates murine sclerodermatous models 
[@kano2019attenuation].

FARP2 encodes FERM, RhoGEF and pleckstrin domain protein 2 and is involved in semaphorin
signaling. Semaphorin signaling was originally identified as being involved in axon guidance,
and has since been shown to be involved in processes highly related to SSc, such as 
angiogenesis, vasculogenesis, and immune regulation. Semaphorin signaling has also beein 
implicated in multiple human diseases including SSc [@besliu2011peripheral,rimar2015semaphorin].



## References














