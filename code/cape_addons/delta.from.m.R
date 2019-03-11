#This function calculates delta values from the m12 and m21 values


delta.from.m <- function(data.obj, p.or.q = 0.05){
	
	int.effects <- data.obj$var.to.var.p.val
	sig.locale <- which(int.effects[,7] <= p.or.q)
	sig.markers <- int.effects[sig.locale,]

	ms <- data.obj$var.to.var.influences
	
	m.mat <- matrix(NA, nrow = dim(sig.markers)[1], ncol = 2)
	
	for(i in 1:dim(sig.markers)[1]){
		pair.locale <- intersect(c(which(ms[,1] == sig.markers[i,1]), which(ms[,2] == sig.markers[i,1])), c(which(ms[,1] == sig.markers[i,2]), which(ms[,2] == sig.markers[i,2])))
		m.mat[i,1] <- ms[pair.locale,"m12"]
		m.mat[i,2] <- ms[pair.locale,"m21"]
		}	

	calc.delta <- function(m12_m21){
		A = matrix(c(1, -1*m12_m21[1], -1*m12_m21[2], 1), byrow = TRUE, ncol = 2)
		b = m12_m21
		deltas <- solve(A, b)
		return(deltas)
		}
	
	all.deltas <- t(apply(m.mat, 1, calc.delta))
	colnames(all.deltas) <- c("delta1", "delta2")

	final.table <- 	cbind(sig.markers[,1:2], all.deltas)
	
	return(final.table)
}