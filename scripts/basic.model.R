# Sean
# load required libraries
library(chromePlus)
library(diversitree)
library(phytools)

############################
### model with XO and XY ###
############################

# chromosome number evolution
# define pars
iter <- 100
prior <- make.prior.exponential(r = 2)

# variables to hold results and depths
results <- vector(mode = "list", length = 100)
tree.depth <- vector(mode = "numeric", length = 100)

# read in data
chrom <- read.csv('../data/final_chrom.csv')
trees <- read.tree("../data/final_100trees")
chrom$gen.prob <- 1
for(i in 1:nrow(chrom)){
  if(chrom$SCS[i] %in% c('XY','NeoXY')){
    chrom$gen.prob[i] <- 0
  }
}
# get the range of chromosome number
rng <- c(range(chrom$Chroms, na.rm = T)[1] - 1,
         range(chrom$Chroms, na.rm = T)[2] + 1)
# make a probability matrix for chromosome number
chrom.mat <- datatoMatrix(x = chrom[,c(2,3,5)],
                          range = rng,
                          hyper = F)
for(i in 1:100){
  print(i)
  tree.depth[i] <- max(branching.times(trees[[i]]))
  trees[[i]]$edge.length <- trees[[i]]$edge.length / tree.depth[i]
  
  # make the likelihood function
  lik <- make.mkn(tree = trees[[i]],
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
  if(i == 1){
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
  results[[i]] <- mcmc(lik = con.lik,
                       x.init = runif(min=0, max=1,
                                      n=length(argnames(con.lik))),
                       nsteps = iter,
                       w = w,
                       prior = prior)
}

# convert rate to millions of years
tree.depth <-c()
trees <- read.tree("../data/final_100trees")
for(i in 1:100){
  tree.depth[i] <- max(branching.times(trees[[i]]))
}
for (i in 1:length(results)){
  results[[i]][,2:3] <- results[[i]][,2:3] / (tree.depth[i]*100)
}
#write the results
saveRDS(results, file = '../results/chrom_number_model_result.rds')

####################################
### model with XO, XY, and NeoXY ###
####################################

# mcmc to estimate rate
trees <- read.tree('../data/final_100trees')
# chrom data 
chrom <- read.csv('../data/final_chrom.csv')
for (i in 1:nrow(chrom)){
  if(chrom$SCS[i] == 'XY'){
    chrom$gen.prob[i] <- 1
  }
  if(chrom$SCS[i] == 'NeoXY'){
    chrom$gen.prob[i] <- 2
  }
  if(chrom$SCS[i] == 'XO'){
    chrom$gen.prob[i] <- 3
  }
}
rownames(chrom) <- chrom$Species
chrom <- chrom[,-1]
scs <- setNames(chrom$gen.prob,row.names(chrom))
# prior
prior <- make.prior.exponential(2)
iter <- 100
# variables to hold results and depths
results <- vector(mode = "list", length = 100)
tree.depth <- vector(mode = "numeric", length = 100)
for ( i in 1:100){
  print(i)
  tree.depth[i] <- max(branching.times(trees[[i]]))
  trees[[i]]$edge.length <- trees[[i]]$edge.length / tree.depth[i]
  # make the likelihood function
  lik <- make.mkn(tree = trees[[i]], states = scs, k=3, strict =T, control = list(method='ode'))
  argnames(lik)
  # safe
  if (i == 1){
    temp <- c()
    temp <- mcmc(lik = lik,
                 x.init = runif(min=0, max=1,n=length(argnames(lik))),
                 prior = prior,
                 nsteps = 100,
                 w = 1)
    # get values for w
    w <- diff(sapply(temp[11:100, 2:(length(argnames(lik))+1)], quantile, c(.05, .95)))
  }
  
  #run mcmc
  results[[i]] <- mcmc(lik = lik,
                       x.init = runif(min=0, max=1,
                                      n=length(argnames(lik))),
                       nsteps = iter,
                       w = w,
                       prior = prior)
}

# convert rate to mya from hundred of mya
tree.depth <-c()
trees <- read.tree("../data/final_100trees")
for(i in 1:100){
  tree.depth[i] <- max(branching.times(trees[[i]]))
}
for (i in 1:length(results)){
  results[[i]][,2:7] <- results[[i]][,2:7] / (tree.depth[i]*100)
}
# write the results
# saveRDS(results, file = '../results/simple_model_scs.rds')

