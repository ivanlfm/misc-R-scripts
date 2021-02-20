### this script combines trees from the posterior distributions of different MCMC runs, then makes a random sample of 1000 trees
### resolves polytomies randomly to make them fully binary
### prunes unwanted taxa

require(ape)
setwd("C:/biogeobears")

#defines root; use this if your trees are unrooted
#root <- "Hypochilus_pococki_Genbank"

#list of unwanted tips to drop
unwanted_tips<-c("Chimerarachne_yingi",
		"Heptathela_kimurai_Genbank",
		"Amaurobius_similis_Genbank",
		"Antrodiaetus_unicolor_Genbank",
		"Archoleptoneta_sp_Genbank",
		"Ariadna_boesenbergi_Genbank",
		"Artema_atlanta_Genbank",
		"Calisoga_longitarsis_Genbank",
		"Drymusa_rengan_Genbank",
		"Eoplectreurys_gertschi",
		"Hickmania_troglodytes_Genbank",
		"Kibramoa_madrona_Genbank",
		"Montsecarachne_amicorum",
		"Otiothops_birabeni_Genbank",
		"Palaeothele_montceauensis",
		"Rosamygale_grauvogeli",
		"Seppo_koponeni",
		"Stegodyphus_mimosarum_Genbank",
		"Thaida_peculiaris_Genbank",
		"Uroctea_durandi_Genbank",
		"Zamilia_quattuormammillae")

### prepares the consensus tree
#consensus of TE-alltaxa-clock-noCOI3-notreeage-FBD1 and TE-alltaxa-clock-noCOI3-notreeage-FBD2
#randomly resolves polytomies and gives a negligible length to  zero-length branches
consensus <- read.tree("TE-alltaxa-clock-noCOI3-notreeage-FBD-1+2.nwk")
consensus <- ladderize(consensus, right = TRUE)
consensus <- multi2di(consensus, random  = TRUE)
consensus$edge.length[consensus$edge.length==0]<-1e-8
write.tree(consensus, "consensus.nwk")
pruned_consensus <- drop.tip(consensus, tip=unwanted_tips)
write.tree(pruned_consensus, "pruned_consensus.nwk")

write.tree(pruned_trees, "pruned_trees.nwk")

### prepares 1000 input trees from trees of the stationary phase of the MCMC

#assuming you have sample 10000 generations, 50% burn-in
n_max_trees = 10001
#burn-in
n_burned_trees = 5001

#combining trees from TE-alltaxa-clock-noCOI3-notreeage-FBD1 and TE-alltaxa-clock-noCOI3-notreeage-FBD2
# IMPORTANT: if trees are from the posterior distribution of MrBayes, ages must be converted prior by running the burntrees.pl script
tree1<-read.tree("1.nwk")
tree2<-read.tree("2.nwk")
tree3<-read.tree("3.nwk")
tree4<-read.tree("4.nwk")

#files have been already burned-in by Burntree, so these lines are commented
#tree1<-tree1[n_burned_trees+1:n_max_trees]
#tree2<-tree2[n_burned_trees+1:n_max_trees]
#tree3<-tree3[n_burned_trees+1:n_max_trees]
#tree4<-tree4[n_burned_trees+1:n_max_trees]

#combines the four files
combined_trees <- c(tree1, tree2, tree3,tree4)
class(combined_trees) <- "multiPhylo"

#samples 1000 trees from the posterior distribution - commented because files have been thinned by BurnTree
#combined_trees <- sample(combined_trees, size = 1000)
combined_trees<-lapply(combined_trees, ladderize, right=TRUE) #ladderizes trees
class(combined_trees) <- "multiPhylo"
write.tree(combined_trees, "combined_trees.nwk")
write.nexus(combined_trees, file = "combined_trees.nex")

#drops unwanted tips (outgroups, etc)
pruned_trees<-lapply(combined_trees, drop.tip, tip=unwanted_tips)
class(pruned_trees) <- "multiPhylo"

#resolves polytomies and zero-length branches
n_polytomous <- 0
for(i in 1:length(pruned_trees)) {
	tree <- pruned_trees[[i]]
		if(is.binary(tree)){} else{
			tree <- multi2di(tree, random = TRUE)
			tree$edge.length[tree$edge.length==0]<-1e-8
			pruned_trees[[i]] <- tree
			n_polytomous <- n_polytomous+1
		}
	}

n_polytomous #number of trees forcibly resolved	
write.tree(pruned_trees, "pruned_trees.nwk")
write.nexus(pruned_trees, file = "pruned_trees.nex")