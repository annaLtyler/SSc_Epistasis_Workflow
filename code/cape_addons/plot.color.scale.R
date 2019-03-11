#This function plots scales for color ramps
#It uses the same color matrix as imageWithText
#If two colors are used, the function pastes them
#together with the light colors in the middle

plot.color.scale <- function(min.x, max.x, col, mid.point = 0, steepness = 1, num.col = 256, main = ""){
		
		if(length(col) > 2){
			stop("This function can only handle two colors maximum")
			}
	
		# #color order: same as in arguments, light to dark 
		all.col.ref <- matrix(
		c("#e5f5f9", "#99d8c9", "#2ca25f",
		"#efedf5", "#bcbddc", "#756bb1",
		"#fee8c8", "#fdbb84", "#e34a33",
		"#deebf7", "#9ecae1", "#3182bd",
		"#f0f0f0", "#bdbdbd", "#636363"), nrow = 3, byrow = FALSE)
		

	#color order: same as in arguments, light to dark 
		#a slightly darker matrix (just purple and green)
		# all.col.ref <- matrix(
		# c("#ccece6", "#66c2a4", "#238b45",
		# "#dadaeb", "#9e9ac8", "#6a51a3",
		# "#fee8c8", "#fdbb84", "#e34a33",
		# "#deebf7", "#9ecae1", "#3182bd",
		# "#f0f0f0", "#bdbdbd", "#636363"), nrow = 3, byrow = FALSE)

		colnames(all.col.ref) <- c("green", "purple", "orange", "blue", "gray")

		# ColorLevels <- seq(min.x, max.x, length=256)
		ColorLevels <- exp.color.fun(min.x, max.x, steepness = 3, num.cols = 256)

			col1.locale <- which(colnames(all.col.ref) == col[1])
			mypal.neg <- colorRampPalette(all.col.ref[3:1,col1.locale]) #put the dark on bottom for the first color

			if(length(col == 2)){
				col2.locale <- which(colnames(all.col.ref) == col[2])
				mypal.pos <- colorRampPalette(all.col.ref[1:3,col2.locale]) #put the dark on top for the second
				}else{
				mypal.pos <- NULL
				}
			ColorRamp <- c(mypal.neg(length(which(ColorLevels < mid.point))), mypal.pos(length(which(ColorLevels >= mid.point))))
		
		
		par(mar = c(3,3,5,2))
		image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp, xlab="",ylab="",xaxt="n", cex.axis = 2, main = main)


}