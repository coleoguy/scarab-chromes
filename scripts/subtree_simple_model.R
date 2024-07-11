library(chromePlus)
library(diversitree)
library(phytools)
# subtrees
trees <- read.tree('../data/final_100trees')
chrom <- read.csv('../data/final_chrom.csv')
sca.tip <- c()
luc.tip <- c()
pas.tip <- c()
for (i in 1:length(chrom$Family)){
  if (chrom$Family[i] == 'Scarabaeidae'){
    sca.tip <- c(sca.tip,chrom$Species[i])
  }
  if (chrom$Family[i] == 'Lucanidae'){
    luc.tip <- c(luc.tip,chrom$Species[i])
  }
  if (chrom$Family[i] == 'Passalidae'){
    pas.tip <- c(pas.tip,chrom$Species[i])
  }
}
sca.tree <- keep.tip.multiPhylo(trees, sca.tip)
luc.tree <- keep.tip.multiPhylo(trees, luc.tip)
pas.tree <- keep.tip.multiPhylo(trees, pas.tip)
sca.chrom <- chrom[chrom$Species %in% sca.tip,]
luc.chrom <- chrom[chrom$Species %in% luc.tip,]
pas.chrom <- chrom[chrom$Species %in% pas.tip,]
# write.tree(sca.tree, file = '../data/subtrees_sca')
# write.tree(luc.tree, file = '../data/subtrees_luc')
# write.tree(pas.tree, file = '../data/subtrees_pas')

# subtree analyses
# define pars
iter <- 100
prior <- make.prior.exponential(r = 2)
# variables to hold results and depths

# read in data
sub <- c('sca','luc','pas')
for (i in 1:3){
  results <- vector(mode = "list", length = 100)
  tree.depth <- vector(mode = "numeric", length = 100)
  chrom <- eval(parse(text = paste(sub[i],'chrom',sep = '.')))
  trees <- eval(parse(text = paste(sub[i],'tree',sep = '.')))
  chrom$gen.prob <- 1
  for(j in 1:nrow(chrom)){
    if(chrom$SCS[j] %in% c('XY','NeoXY')){
      chrom$gen.prob[j] <- 0
    }
  }
  
  # get the range of chromosome number
  rng <- c(range(chrom$Chroms, na.rm = T)[1] - 1,
           range(chrom$Chroms, na.rm = T)[2] + 1)
  
  # make a probability matrix for chromosome number
  
  chrom.mat <- datatoMatrix(x = chrom[,c(2,3,5)],
                            range = rng,
                            hyper = F)
  for(k in 1:100){
    print(k)
    tree.depth[k] <- max(branching.times(trees[[k]]))
    trees[[k]]$edge.length <- trees[[k]]$edge.length / tree.depth[k]
    
    # make the likelihood function
    lik <- make.mkn(tree = trees[[k]],
                    states = chrom.mat,
                    k = ncol(chrom.mat),
                    strict = FALSE,
                    control = list(method = "ode"))
    argnames(lik)
    con.lik <- constrainMkn(data = chrom.mat,
                            lik = lik,
                            polyploidy = F,
                            hyper = F,
                            constrain = list(drop.demi = T,
                                             drop.poly = T))
    # safe
    argnames(con.lik)
    if(k == 1){
      temp <- c()
      temp <- mcmc(lik = con.lik,
                   x.init = runif(min=0, max=1,
                                  n=length(argnames(con.lik))),
                   prior = prior,
                   nsteps = 100,
                   w = 1)
      # get values for w
      w <- diff(sapply(temp[11:100, 2:(length(argnames(con.lik))+1)], quantile, c(.05, .95)))
    }
    # run MCMC
    results[[k]] <- mcmc(lik = con.lik,
                         x.init = runif(min=0, max=1,
                                        n=length(argnames(con.lik))),
                         nsteps = iter,
                         w = w,
                         prior = prior)
  }

  # convert rate to millions of years
  tree.depth <-c()
  for(l in 1:100){
    tree.depth[l] <- max(branching.times(trees[[l]]))
  }  
  
  for (m in 1:length(results)){
    results[[m]][,2:3] <- results[[m]][,2:3] / (tree.depth[m]*100)
  }
  #save data
  write.csv(do.call(rbind,results), file = paste(paste('../results/simple',sub[i], sep = '_'),'csv', sep = '.'), row.names = F)
}
