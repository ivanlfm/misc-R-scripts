#plots histogram and density

#set working directory
setwd("c:/a")

#read data
data <- read.csv("Galapagos_dispersal.csv")
names(data)

xlim=100 #limit for x axis
ylim=0.2 #limit for y axis
breaks=20 #number of bins for histogram
density_over_histogram=TRUE #plots density over histogram; FALSE plots in separate graph

#for each column in data, plots a kernel density plot and a histogram and saves them to PDFs
for(a in 1:length(names(data))){
	density <- density(data[,a][!is.na(data[,a])], bw=2)  
	pdf(file=paste(names(data)[a],"_histogram.pdf",sep=""))
		hist(data[,a][!is.na(data[,a])],xlab = names(data)[a],freq=FALSE,ylim=c(0,ylim), xlim=c(0,xlim),breaks=breaks)
		abline(v=15, col ="gray", lwd=2.5, lty=2)
		if (density_over_histogram) lines(density)
	dev.off()
	if (density_over_histogram==FALSE){
		pdf(file=paste(names(data)[a],"_density.pdf",sep=""))
			plot(density, xlab = names(data)[a], ylim=c(0,ylim), xlim=c(0,xlim))
		dev.off()
	}
}
