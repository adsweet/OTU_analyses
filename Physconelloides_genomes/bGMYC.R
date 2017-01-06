###############################################################################################
#Use this script to sample random trees from a distribution of trees (i.e. a BEAST .trees file)
#and run bGMYC to estimate optimal OTU clustering
#Used in Sweet et al. in prep
################################################################################################

library(ape)
library(bGMYC)

#READ IN RANDOM TREES
trees <- read.nexus('MT_burnin.trees')
tree <- sample(trees, size=100)
write.tree(tree, "MT_random100.trees")

#ROOT ALL TREES
trees.root = NULL
for (i in 1:length(tree)) {
  trees.root[[i]] =  root(tree[[i]], c("Campanulotes_bidentatus_compar"))
}

#REMOVE OUTGROUPS
trees.noout = NULL
for (i in 1:length(trees.root)) {
  trees.noout[[i]] =  drop.tip(trees.root[[i]], c("Campanulotes_bidentatus_compar"))
}

#ASSESS MCMC AND BURNIN LENGTH
bgmyc.singlephy(trees.noout[[1]], mcmc=50000, burnin=1, thinning=10, t1=2, t2=35, start=c(1,1,25))->result.single

#RUN BGMYC
result.multi <- bgmyc.multiphylo(trees.noout, mcmc=20000, burnin=10000, thinning=10, t1=2, t2=35, start=c(1,1,25))
plot(result.multi)
result.spec <- bgmyc.spec(result.multi, file='Physcon_bgmyc_results.csv')
result.probmat <- spec.probmat(result.multi)
plot(result.probmat,trees.noout[[1]])
out <- bgmyc.point(result.probmat, 0.05)
out <- data.frame(out)
write(out, "Physcon_point_0.05.txt")
