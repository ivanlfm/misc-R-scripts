### this script combines trees from the posterior distributions of different MCMC runs, then makes a random sample of 1000 trees
### resolves polytomies randomly to make them fully binary
### prunes unwanted taxa

require(ape)
setwd("C:/a")
records <- read.table("all_records_pruned.txt")
species <- levels(records$V1)

### prepares the maximum clade credibility tree
#randomly resolves polytomies and gives a negligible length to  zero-length branches
consensus <- read.nexus("TE-red-BD-exp-FosSec-COI3JC-16SHKY+g-Zyg.1+2+3_MCCT.tre")
consensus <- ladderize(consensus, right = TRUE)
consensus <- multi2di(consensus, random  = TRUE)
consensus$edge.length[consensus$edge.length==0]<-1e-8
pruned_consensus <- keep.tip(consensus, tip=species)
write.tree(pruned_consensus, "TE-red-BD-exp-FosSec-COI3JC-16SHKY+g-Zyg.1+2+3_MCCT.tre.nwk")

posterior_trees <- read.nexus("TE-red-BD-exp-FosSec-COI3JC-16SHKY+g-Zyg.1+2+3.trees")
class(posterior_trees) <- "multiPhylo"

#samples 1000 trees from the posterior distribution
posterior_trees <- sample(posterior_trees, size = 1000)
posterior_trees<-lapply(posterior_trees, ladderize, right=TRUE) #ladderizes trees
class(posterior_trees) <- "multiPhylo"
write.nexus(posterior_trees, file = "1000_posterior_trees.nex")

#drops unwanted tips (outgroups, etc)
pruned_trees<-lapply(posterior_trees, keep.tip, tip=species)
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
write.nexus(pruned_trees, file = "1000_pruned_trees.nex")