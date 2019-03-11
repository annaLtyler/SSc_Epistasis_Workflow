#This function moves specified markers from the g.covar matrix and table to the genotype matrix

covar2marker <- function(data.obj, markers){

	covar.info <- data.obj$g.covar
	covar.data <- data.obj$g.covar.table
	
	if(is.character(markers)){
		marker.locale <- match(markers, covar.info[1,])
		}else{
		marker.locale <- markers	
		}
	
	#get the number of the marker	
	marker.num <- colnames(covar.info)[marker.locale]
	marker.name <- covar.info[1,marker.locale]
	marker.chr <- covar.info[2,marker.locale]
	marker.pos <- covar.info[3,marker.locale]
	
	for(m in 1:length(markers)){
		less.pos <- which(as.numeric(data.obj$marker.num) < as.numeric(marker.num[m]))
		greater.pos <- which(as.numeric(data.obj$marker.num) > as.numeric(marker.num[m]))
		
		if(!is.null(data.obj$geno)){
			data.obj$geno <- cbind(data.obj$geno[,less.pos], covar.data[,marker.locale[m]], data.obj$geno[,greater.pos])
			}
		data.obj$marker.names <- c(data.obj$marker.names[less.pos], marker.name[m], data.obj$marker.names[greater.pos])
		data.obj$marker.num <- c(data.obj$marker.num[less.pos], marker.num[m], data.obj$marker.num[greater.pos])
		data.obj$chromosome <- c(data.obj$chromosome[less.pos], marker.chr[m], data.obj$chromosome[greater.pos])
		data.obj$marker.location <- c(data.obj$marker.location[less.pos], marker.pos[m], data.obj$marker.location[greater.pos])
		}	
		
	data.obj$g.covar <- covar.info[,-marker.locale,drop=FALSE]
	data.obj$g.covar.table <- covar.data[,-marker.locale,drop=FALSE]
	
	if(dim(data.obj$g.covar)[2] == 0){
		data.obj$g.covar <- NULL
		data.obj$g.covar.table <- NULL
		}
		
	return(data.obj)
	
	
}