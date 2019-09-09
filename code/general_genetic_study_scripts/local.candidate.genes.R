#This function finds genes that are significantly associated with
#a peak locus and with a trait of interest. These are potential
#causal candidate genes.

local.candidate.genes <- function(out.var, exp.var, mediators, covar = NULL,
sig.p = 1e-6, test.allele.effects = TRUE, plot.results = FALSE){
	#========================================================
	#adjust variables for covariates.
	#========================================================
    if(!is.null(covar)){
        cat("Adjusting variables for covariates...\n")
        adj.out <- adjust(out.var, covar)
        adj.med <- adjust(mediators, covar)
    }else{
        adj.out <- out.var
        adj.med <- mediators
    }

	#========================================================
	#do some quick checks for variance in the input matrices
	#========================================================
	
	test.var <- round(apply(adj.out, 2, function(x) var(x, na.rm = TRUE)), 2)
	if(any(test.var == 0)){
		stop("some outcome variables have 0 variance.")
		}
	
	test.var <- round(apply(adj.med, 2, function(x) var(x, na.rm = TRUE)), 2)
	if(any(test.var == 0)){
		zero.locale <- which(test.var == 0)
		adj.med <- adj.med[,-zero.locale]
		warning(paste("removing", length(zero.locale), "mediators with zero variance"))
		}
	
	test.var <- round(apply(adj.med, 1, function(x) var(x, na.rm = TRUE)), 4)
	na.locale <- which(is.na(test.var))
	zero.locale <- which(test.var == 0)
	if(length(na.locale) > 0 || length(zero.locale) > 0){
		remove.locale <- unique(c(na.locale, zero.locale))
		adj.med <- adj.med[-remove.locale,]
		adj.out <- adj.out[-remove.locale,]
		exp.var <- exp.var[-remove.locale]
		warning(paste("removing", length(remove.locale), "individuals with zero expression variance"))
		}


    model.p <- function(model){
        f <- summary(model)$fstatistic
	    p <- pf(f[1],f[2],f[3],lower.tail=F)
        return(p)
    }

    #calculate the relationship between the independent variable and the mediator(s)
    ind.to.med.model <- apply(adj.med, 2, function(x) lm(x~exp.var))
    ind.to.med.p <- sapply(ind.to.med.model, model.p)
    #qqunif.plot(int.to.med.p)

    #calculate the relationship between the mediators and the dependent variable
    candidate.genes <- vector(mode = "list", length = ncol(adj.out))
    names(candidate.genes) <- colnames(adj.out)

    for(i in 1:ncol(out.var)){
        med.to.dep.model <- apply(adj.med, 2, function(x) lm(adj.out[,i]~x))
        med.to.dep.p[,i] <- sapply(med.to.dep.model, model.p)
        #qqunif.plot(med.to.dep.p[,i])
        if(plot.results){
            plot(-log10(ind.to.med.p), -log10(med.to.dep.p[,i]), 
            xlab = "Genotype -> Mediator -log10(P)",
            ylab = "Mediator -> Trait -log10(P)")
            abline(v = -log10(sig.p), h = -log10(sig.p))
        }
        sig.med <- intersect(which(ind.to.med.p < sig.p), which(med.to.dep.p < sig.p))

        #if we are looking at allele effects, test the relationship between the
        #independent variable and the dependent variable
        #and compare the allele effects for each of the significant candidates
        if(test.allele.effects){
            ind.to.dep.model <- apply(adj.out, 2, function(x) lm(x~exp.var))
            allele.effects <- lapply(ind.to.dep.model, coef)
            allele.med.effects <- lapply(ind.to.med.model[sig.med], coef)
            allele.effects.cor.mat <- matrix(NA, ncol = ncol(adj.out), nrow = length(sig.med))
            colnames(allele.effects.cor.mat) <- colnames(adj.out)
            rownames(allele.effects.cor.mat) <- colnames(adj.med)[sig.med]
            for(j in 1:length(ind.to.dep.model)){
                for(k in 1:length(sig.med)){
                    allele.effects.cor.mat[k,j] <- cor(allele.effects[[i]], allele.med.effects[[k]])
                }
            }
        candidate.genes[[i]] <- allele.effects.cor.mat
        }else{
        candidate.genes[[i]] <- colnames(mediators)[sig.med]
        }

        
        }

    return(candidate.genes)


}
