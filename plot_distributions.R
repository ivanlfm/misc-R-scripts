### plots individual maps for each species from a database
### colors them by attribute (e.g. source of data)

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(dplyr)

setwd("c:/a")
genusname <- "Micrathena"

#prepares map with desired countries
countrylist<-scan("america.txt",as.character())
lat<-map_data("world",region=countrylist)
idx<-which(lat$long<=-75 & lat$lat<=-20)
set.seed(42)

#reads records, generate list of unique species
records<-read.table("micrathena.txt",h=T,sep="\t",stringsAsFactors = F)
species<-unique(records$species)
plot.records <- data.frame(
  long = records$long,
  lat = records$lat,
  species = records$species,
  source = records$source,
  stringsAsFactors = FALSE
  )  

#plots map with all records
g1<-ggplot() + geom_polygon(data = lat[-idx,], aes(x=long, y = lat, group = group),fill='white',colour="black") +
			coord_fixed(1.5) +
			geom_point(data = plot.records, aes(x = long, y = lat), color = "black", size = 0.5) +
			geom_point(data = plot.records, aes(x = long, y = lat), color = "red", size = 0.4) 

filename <-  paste(genusname, ".pdf", sep = "")
pdf(filename)		
print(g1)
dev.off()


#prints individual maps for each species
for(i in 1:length(species)) {
	sp <- species[i]
	filename <-  paste(sp, ".jpg", sep = "")
	sp.record <- filter(plot.records, species == sp)
	#categorizes records by source -- this is awfully clumsy!!!
			literature <- filter(sp.record , source == "Levi1985" | source == "Bonaldo1990" | source == "Lise1995")
			revised <- filter(sp.record , source == "UFMG" | source == "UFPI" | source == "Magalhaes2011" | source == "Magalhaes2017")
			gbif <- filter(sp.record , source == "GBIF")
			inat <- filter(sp.record , source == "iNaturalist")
			ibsp <- filter(sp.record , source == "Cizaukas")
	g1<-ggplot() + geom_polygon(data = lat[-idx,], aes(x=long, y = lat, group = group),fill='white',colour="black") +
			coord_fixed(1.5) +
			geom_point(data = literature, aes(x = long, y = lat), color = "black", size = 7) +
			geom_point(data = revised, aes(x = long, y = lat), color = "red", size = 7) +
			geom_point(data = gbif, aes(x = long, y = lat), color = "blue", size = 2) +
			geom_point(data = inat, aes(x = long, y = lat), color = "yellow", size = 2) +
			geom_point(data = ibsp, aes(x = long, y = lat), color = "green", size = 2) +

	jpeg(filename, width = 1500, height = 1500, quality = 100)		
	print(g1)
	dev.off()
	}
