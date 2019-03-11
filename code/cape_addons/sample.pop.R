#This function subsets a data object for testing

sample.pop <- function(data.obj, num.ind = NULL, num.gene = NULL){
	
	if(is.null(num.ind) && is.null(num.gene)){
		stop("either num.ind or num.gene must be specified.")
		}
	
	if(!is.null(num.gene)){
		gene.sample <- sort(sample(1:length(data.obj$marker.names), num.gene))
		data.obj$marker.names <- data.obj$marker.names[gene.sample]
		data.obj$marker.num <- 1:num.gene
		data.obj$chromosome <- data.obj$chromosome[gene.sample]
		data.obj$marker.location <- data.obj$marker.location[gene.sample]
		}
	
	if(!is.null(num.ind)){
		ind.sample <- sort(sample(1:dim(data.obj$pheno)[1], num.ind))
		data.obj$pheno <- data.obj$pheno[ind.sample,]
		}
	
	return(data.obj)
	
}