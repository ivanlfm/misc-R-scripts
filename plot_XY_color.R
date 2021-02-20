#plots a XY graph with dots colored by value


#sets working directory, reads data
	setwd("c:/a/plot")
	data <- read.csv("supports_free_model.csv")

#defines limits for axis and for the color ramp
	xlim <- c(0,38)
	ylim <- c(0,1)
	#colors for the color ramp
		color1 <- "white"
		color2 <- "black"
		nbreaks <- 10

# defines n of plots
	nplots <- length(data[1,])

for(i in 2:nplots) {
	#defines which columns from the csv you want to plot
		columnx <- 1
		columny <- i
		
	#defines a color ramp with 10 breaks, ranging from ylim[1] (has color1) to ylim[2] (has color2)
		Pal <- colorRampPalette(c(color1,color2))
		data$Col <- Pal(nbreaks)[as.numeric(cut(c(ylim[1],ylim[2],data[,columny]),breaks = nbreaks)[3:102])]
		
	#plots and saves as pdf
		pdf(file=paste(names(data)[columnx],"_",names(data)[columny],".pdf",sep=""))
				plot(data[,columnx], data[,columny], xlim = xlim, ylim = ylim, xlab = names(data[columnx]), ylab= names(data[columny]), pch= 21, col = "black", bg=data$Col, cex=2)
				abline(v=15, col ="gray", lwd=2.5, lty=2)
				#legend, have to work this out fully
					legend_image <- as.raster(matrix(Pal(10), ncol=1))
					rasterImage(legend_image, max(xlim)-0.1, max(ylim), max(xlim)-2,min(ylim), interpolate = TRUE)
			dev.off()
	}