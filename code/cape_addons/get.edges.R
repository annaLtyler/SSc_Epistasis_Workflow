#This function generates an edge list for get.network2 from
#a cape object at a given p value

get.edges <- function(data.obj, p.or.q = 0.05){
	
	all.net.data <- data.obj$var.to.var.p.val
	pheno.tables <- data.obj$max.var.to.pheno.influence

	sig.inter.locale <- which(all.net.data[,7] <= p.or.q)
	sig.inter <- all.net.data[sig.inter.locale,1:2,drop=FALSE]
	
	sig.main <- NULL
	for(ph in 1:length(pheno.tables)){
		sig.main.locale <- which(pheno.tables[[ph]][,7] <= p.or.q)
		sig.main.edges <- cbind(pheno.tables[[ph]][sig.main.locale,1], rep(names(pheno.tables)[ph], length(sig.main.locale)))
		sig.main <- rbind(sig.main, sig.main.edges)
		}
	
	full.edge.list <- rbind(sig.inter, sig.main)
	return(full.edge.list)

	}
