#This function compares the effects of markers from two different singlescan objects


compare.singlescan <- function(singlescan.obj1, singlescan.obj2, label1 = NULL, label2 = NULL){
	
	
	single.effects1 <- singlescan.obj1$singlescan.results
	single.effects2 <- singlescan.obj2$singlescan.results	
	
	num.pheno <- length(single.effects1)

	
	if(is.null(label1)){label1 = "singlescan 1"}
	if(is.null(label2)){label2 = "singlescan 2"}
	
	layout.mat <- get.layout.mat(num.pheno)
	layout(layout.mat)
	for(p in 1:num.pheno){
		std.effects1 <- single.effects1[[p]][,"t.stat"]
		std.effects2 <- single.effects2[[p]][,"t.stat"]
		plot(std.effects1, std.effects2, xlab = label1, ylab = label2, main = names(single.effects1)[p])
		abline(0,1)
		}

}