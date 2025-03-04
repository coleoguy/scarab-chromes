# load required libraries
library(chromePlus)
library(diversitree)
library(phytools)
library(plyr)
############################
### model with XO and XY ###
############################

# chromosome number evolution

# read in data
allchrom <- read.csv('../data/SpeciesChromList.csv')
trees <- read.tree("../data/final100trees")
# define pars
iter <- 100
prior <- make.prior.exponential(r = 2)
# variables to hold results and depths
results <- vector(mode = "list", length = 100)
tree.depth <- vector(mode = "numeric", length = 100)

for(i in 1:100){
  print(i)
  ##################
  ### chrom data ###
  ##################
  # genus level might have different data 
  # so we need to random pick data everytime 
  # subset chrom data
  tip.names <- trees[[1]]$tip.label
  chrom <- data.frame()
  # generate chrom data 
  for (j in 1:length(tip.names)){
    if (tip.names[j] %in% allchrom$Name){
      if (length(which(allchrom$Name == tip.names[j])) == 1){
        hitdat <- allchrom[which(allchrom$Name == tip.names[j]),c(1,4,6,8)]
        colnames(hitdat) <- c("Family","Species","Chroms","SCS")
        chrom <-rbind(chrom, hitdat) 
      }
      if (length(which(allchrom$Name == tip.names[j])) > 1){
        pick <- sample(length(which(allchrom$Name == tip.names[j])),1)
        hitdat <- allchrom[which(allchrom$Name == tip.names[j])[pick],c(1,4,6,8)]
        colnames(hitdat) <- c("Family","Species","Chroms","SCS")
        chrom <-rbind(chrom, hitdat)
      }
    }
    if (tip.names[j] %in% allchrom$Genus){
      if (length(which(allchrom$Genus == tip.names[j])) == 1){
        hitdat <- allchrom[which(allchrom$Genus == tip.names[j]),c(1,2,6,8)]
        colnames(hitdat) <- c("Family","Species","Chroms","SCS")
        chrom <-rbind(chrom, hitdat)
      }
      if (length(which(allchrom$Genus == tip.names[j])) > 1){
        pick <- (sample(length(which(allchrom$Genus == tip.names[j])),1))
        hitdat <- allchrom[which(allchrom$Genus == tip.names[j])[pick],c(1,2,6,8)]
        colnames(hitdat) <- c("Family","Species","Chroms","SCS")
        chrom <-rbind(chrom, hitdat)
      }
    }
  }
  chrom$SCS <- sub("XXXXXO", "XO", chrom$SCS)
  chrom$gen.prob <- 1
  for(k in 1:nrow(chrom)){
    if(chrom$SCS[k] %in% c('XY','NeoXY')){
      chrom$gen.prob[k] <- 0
    }
  }
  # get the range of chromosome number
  rng <- c(range(chrom$Chroms, na.rm = T)[1] - 1,
           range(chrom$Chroms, na.rm = T)[2] + 1)
  # make a probability matrix for chromosome number
  chrom.mat <- datatoMatrix(x = chrom[,c(2,3,5)],
                            range = rng,
                            hyper = F)
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
trees <- read.tree("../data/final100trees")
for(i in 1:100){
  tree.depth[i] <- max(branching.times(trees[[i]]))
}
for (i in 1:length(results)){
  results[[i]][,2:3] <- results[[i]][,2:3] / (tree.depth[i]*100)
}
#write the results
# saveRDS(results, file = '../results/chrom_number_model_result.rds')

####################################
### model with XO, XY, and NeoXY ###
####################################
trees <- read.tree('../data/final100trees')
allchrom <- read.csv('../data/SpeciesChromList.csv')
# prior
prior <- make.prior.exponential(2)
iter <- 100
# variables to hold results and depths
results <- vector(mode = "list", length = 100)
tree.depth <- vector(mode = "numeric", length = 100)
for ( i in 1:100){
  print(i)
  #### chrom data ###
  tip.names <- trees[[1]]$tip.label
  chrom <- data.frame()
  for (j in 1:length(tip.names)){
    if (tip.names[j] %in% allchrom$Name){
      if (length(which(allchrom$Name == tip.names[j])) == 1){
        hitdat <- allchrom[which(allchrom$Name == tip.names[j]),c(1,4,6,8)]
        colnames(hitdat) <- c("Family","Species","Chroms","SCS")
        chrom <-rbind(chrom, hitdat) 
      }
      if (length(which(allchrom$Name == tip.names[j])) > 1){
        pick <- sample(length(which(allchrom$Name == tip.names[j])),1)
        hitdat <- allchrom[which(allchrom$Name == tip.names[j])[pick],c(1,4,6,8)]
        colnames(hitdat) <- c("Family","Species","Chroms","SCS")
        chrom <-rbind(chrom, hitdat)
      }
    }
    if (tip.names[j] %in% allchrom$Genus){
      if (length(which(allchrom$Genus == tip.names[j])) == 1){
        hitdat <- allchrom[which(allchrom$Genus == tip.names[j]),c(1,2,6,8)]
        colnames(hitdat) <- c("Family","Species","Chroms","SCS")
        chrom <-rbind(chrom, hitdat)
      }
      if (length(which(allchrom$Genus == tip.names[j])) > 1){
        pick <- (sample(length(which(allchrom$Genus == tip.names[j])),1))
        hitdat <- allchrom[which(allchrom$Genus == tip.names[j])[pick],c(1,2,6,8)]
        colnames(hitdat) <- c("Family","Species","Chroms","SCS")
        chrom <-rbind(chrom, hitdat)
      }
    }
  }
  chrom$SCS <- sub("XXXXXO", "XO", chrom$SCS)
  for (k in 1:nrow(chrom)){
    if(chrom$SCS[k] == 'XY'){
      chrom$gen.prob[k] <- 1
    }
    if(chrom$SCS[k] == 'NeoXY'){
      chrom$gen.prob[k] <- 2
    }
    if(chrom$SCS[k] == 'XO'){
      chrom$gen.prob[k] <- 3
    }
  }
  rownames(chrom) <- chrom$Species
  chrom <- chrom[,-1]
  scs <- setNames(chrom$gen.prob,row.names(chrom))
  ###################
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
