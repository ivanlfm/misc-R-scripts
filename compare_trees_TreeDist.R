### this script compares the numer of nodes of trees, collect nodal support, and estimate distances between trees

### loads packages
library(ape)
library(phylotate)
library(phangorn)
library(phytools)
library(TreeDist)
setwd("C:/micrathena_supports")

#defines outgroups
og <- c("Actinosoma_pentacanthum_166",
"Alpaida_veniliae_173",
"Araneus_diadematus_GenBank",
"Araneus_venatrix_GenBank",
"Argiope_argentata_GenBank",
"Aspidolasius_branicki_200",
"Cyclosa_conica_GenBank",
"Encyosaccus_sexmaculatus_202",
"Gasteracantha_cancriformis_GenBank",
"Hypognatha_cryptocephala_001",
"Hypognatha_cryptocephala_002",
"Hypognatha_sp_153",
"Pronous_peje_JCG",
"Pronous_tuberculifer_209",
"Spinepeira_schlingeri_201",
"Verrucosa_arenata_GenBank",
"Wagneriana_dimastophora_GenBank",
"Wagneriana_maseta_163",
"Xylethrus_superbus_IFM0083",
"Zygiella_atrica_GenBank")

### reads trees
coiH3 <- read_annotated("coih3.tre", format="nexus")
mag <- read_annotated("ma-g.tre", format="nexus")
mao <- read_annotated("ma-o.tre", format="nexus")
mug <- read_annotated("mu-g.tre", format="nexus")
muo <- read_annotated("mu-o.tre", format="nexus")

### prunes trees to include only terminals of the tree with less terminals
species <- coiH3$tip.label
p.mag <- keep.tip(mag, species)
p.mao <- keep.tip(mao, species)
p.mug <- keep.tip(mug, species)
p.muo <- keep.tip(muo, species)

#removes outgroups
coiH3 <- drop.tip(coiH3, og)
p.mag <- drop.tip(p.mag, og)
p.mao <- drop.tip(p.mao, og)
p.mug <- drop.tip(p.mug, og)
p.muo <- drop.tip(p.muo, og)

### count1s the number of nodes
Nnodes <- c(coiH3$Nnode, p.mag$Nnode, p.mao$Nnode, p.mug$Nnode, p.muo$Nnode)
names(Nnodes) <- c("coiH3", "mag", "mao", "mug", "muo")
Nnodes
write.table(Nnodes, file = "Nnodes.txt")

### checks distances between trees
distances <- matrix(NA, nrow = 4, ncol = 1)
rownames(distances) <- c("coiH3xMAG", "coiH3xMAO", "coiH3xMUG", "coiH3xMUO")
distances[1,] <- TreeDistance(coiH3, p.mag)
distances[2,] <- TreeDistance(coiH3, p.mao)
distances[3,] <- TreeDistance(coiH3, p.mug)
distances[4,] <- TreeDistance(coiH3, p.muo)
write.table(distances, "distances_TreeDist.txt")

### if trees has polytomies,
### randomly resolves polytomies X times
### later calculations will take a long time if number of replications is high!
replications <- 100
# creates empty lists of trees...
dichot_coiH3 <- vector(mode = "list", length = replications)
class(dichot_coiH3) <- "multiPhylo"
dichot_p.mag <- vector(mode = "list", length = replications)
class(dichot_p.mag) <- "multiPhylo"
dichot_p.mao <- vector(mode = "list", length = replications)
class(dichot_p.mao) <- "multiPhylo"
dichot_p.mug <- vector(mode = "list", length = replications)
class(dichot_p.mug) <- "multiPhylo"
dichot_p.muo <- vector(mode = "list", length = replications)
class(dichot_p.muo) <- "multiPhylo"
# ... and populates them with randomly resolved trees
for(i in 1:replications) {
	dichot_coiH3[[i]] <- multi2di(coiH3, random  = TRUE)
	dichot_p.mag[[i]] <- multi2di(p.mag, random  = TRUE)
	dichot_p.mao[[i]] <- multi2di(p.mao, random  = TRUE)
	dichot_p.mug[[i]] <- multi2di(p.mug, random  = TRUE)
	dichot_p.muo[[i]] <- multi2di(p.muo, random  = TRUE)
}
#saves trees just in case
write.nexus(dichot_coiH3, file = "dichot_coiH3.nex")
write.nexus(dichot_p.mag, file = "dichot_p.mag.nex")
write.nexus(dichot_p.mao, file = "dichot_p.mao.nex")
write.nexus(dichot_p.mug, file = "dichot_p.mug.nex")
write.nexus(dichot_p.muo, file = "dichot_p.muo.nex")

### then compares distances among randomly resolved trees
dist.coiH3.pmag <- matrix(NA, nrow = replications*replications, ncol = 1)
count2 <- 1
count1 <- 1
for(i in 1:replications) {
	for(i in 1:replications) {
	dist.coiH3.pmag[count1,] <- TreeDistance(dichot_coiH3[[count2]], dichot_p.mag[[i]])
	count1 <- count1+1 
	}
	count2 <- count2+1
}

dist.coiH3.pmao <- matrix(NA, nrow = replications*replications, ncol = 1)
count2 <- 1
count1 <- 1
for(i in 1:replications) {
	for(i in 1:replications) {
	dist.coiH3.pmao[count1,] <- TreeDistance(dichot_coiH3[[count2]], dichot_p.mao[[i]])
	count1 <- count1+1 
	}
	count2 <- count2+1
}

dist.coiH3.pmug <- matrix(NA, nrow = replications*replications, ncol = 1)
count2 <- 1
count1 <- 1
for(i in 1:replications) {
	for(i in 1:replications) {
	dist.coiH3.pmug[count1,] <- TreeDistance(dichot_coiH3[[count2]], dichot_p.mug[[i]])
	count1 <- count1+1 
	}
	count2 <- count2+1
}

dist.coiH3.pmuo <- matrix(NA, nrow = replications*replications, ncol = 1)
count2 <- 1
count1 <- 1
for(i in 1:replications) {
	for(i in 1:replications) {
	dist.coiH3.pmuo[count1,] <- TreeDistance(dichot_coiH3[[count2]], dichot_p.muo[[i]])
	count1 <- count1+1 
	}
	count2 <- count2+1
}

distDich <- list(dist.coiH3.pmag, dist.coiH3.pmao, dist.coiH3.pmug, dist.coiH3.pmuo)
names(distDich) <- c("mag", "mao", "mug", "muo")
write.table(distDich, "distDich.txt")

### collects and print supports from the trees (supporst include also de 100% associated to terminal tips,
### so I get only the second part of the vector
coiH3_supps <- coiH3$node.comment[(length(coiH3$tip.label)+1):length(coiH3$node.comment)]
mag_supps <- mag$node.comment[(length(mag$tip.label)+1):length(mag$node.comment)]
mao_supps <- mao$node.comment[(length(mao$tip.label)+1):length(mao$node.comment)]
mug_supps <- mug$node.comment[(length(mug$tip.label)+1):length(mug$node.comment)]
muo_supps <- muo$node.comment[(length(muo$tip.label)+1):length(muo$node.comment)]

write(coiH3_supps, "coiH3.txt")
write(mag_supps, "mag.txt")
write(mao_supps, "mao.txt")
write(mug_supps, "mug.txt")
write(muo_supps, "muo.txt")