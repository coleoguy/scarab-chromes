# loading required library
library(chromePlus)
library(diversitree)
library(phytools)
# subtrees
trees <- read.tree('../data/final100trees')
allchrom <- read.csv('../data/SpeciesChromList.csv')
sca.tip <- c()
luc.tip <- c()
pas.tip <- c()
tip.names <- trees[[1]]$tip.label
for (i in 1:length(tip.names)){
  if (tip.names[i] %in% allchrom$Name){
    if (allchrom$Family[which(allchrom$Name == tip.names[i])[1]] == "Scarabaeidae"){
      sca.tip <- c(sca.tip, tip.names[i])
    }
    if (allchrom$Family[which(allchrom$Name == tip.names[i])[1]] == "Lucanidae"){
      luc.tip <- c(luc.tip, tip.names[i])
    }
    if (allchrom$Family[which(allchrom$Name == tip.names[i])[1]] == "Passalidae"){
      pas.tip <- c(pas.tip, tip.names[i])
    }
  }
  if (tip.names[i] %in% allchrom$Genus){
    if (allchrom$Family[which(allchrom$Genus == tip.names[i])[1]] == "Scarabaeidae"){
      sca.tip <- c(sca.tip, tip.names[i])
    }
    if (allchrom$Family[which(allchrom$Genus == tip.names[i])[1]] == "Lucanidae"){
      luc.tip <- c(luc.tip, tip.names[i])
    }
    if (allchrom$Family[which(allchrom$Genus == tip.names[i])[1]] == "Passalidae"){
      pas.tip <- c(pas.tip, tip.names[i])
    }
  }
}
# create sub trees
sca.tree <- keep.tip.multiPhylo(trees, sca.tip)
luc.tree <- keep.tip.multiPhylo(trees, luc.tip)
pas.tree <- keep.tip.multiPhylo(trees, pas.tip)

# subtree analyses
# define pars
iter <- 100
prior <- make.prior.exponential(r = 2)
sub <- c('sca','luc','pas')
for (i in 1:3){
  # subtree
  results <- vector(mode = "list", length = 100)
  tree.depth <- vector(mode = "numeric", length = 100)
  trees <- eval(parse(text = paste(sub[i],'tree',sep = '.')))
  for(k in 1:100){
    print(k)
    # chrom data
    # create subchrom
    chrom <- data.frame()
    tip.names <- eval(parse(text = paste(sub[i],'tip',sep = '.')))
    # generate chrom data 
    for (n in 1:length(tip.names)){
      if (tip.names[n] %in% allchrom$Name){
        if (length(which(allchrom$Name == tip.names[n])) == 1){
          hitdat <- allchrom[which(allchrom$Name == tip.names[n]),c(1,4,6,8)]
          colnames(hitdat) <- c("Family","Species","Chroms","SCS")
          chrom <-rbind(chrom, hitdat) 
        }
        if (length(which(allchrom$Name == tip.names[n])) > 1){
          pick <- sample(length(which(allchrom$Name == tip.names[n])),1)
          hitdat <- allchrom[which(allchrom$Name == tip.names[n])[pick],c(1,4,6,8)]
          colnames(hitdat) <- c("Family","Species","Chroms","SCS")
          chrom <-rbind(chrom, hitdat)
        }
      }
      if (tip.names[n] %in% allchrom$Genus){
        if (length(which(allchrom$Genus == tip.names[n])) == 1){
          hitdat <- allchrom[which(allchrom$Genus == tip.names[n]),c(1,2,6,8)]
          colnames(hitdat) <- c("Family","Species","Chroms","SCS")
          chrom <-rbind(chrom, hitdat)
        }
        if (length(which(allchrom$Genus == tip.names[n])) > 1){
          pick <- (sample(length(which(allchrom$Genus == tip.names[n])),1))
          hitdat <- allchrom[which(allchrom$Genus == tip.names[n])[pick],c(1,2,6,8)]
          colnames(hitdat) <- c("Family","Species","Chroms","SCS")
          chrom <-rbind(chrom, hitdat)
        }
      }
    }
    chrom$SCS <- sub("XXXXXO", "XO", chrom$SCS)
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
  # write.csv(do.call(rbind,results), file = paste(paste('../results/sub',sub[i], sep = '_'),'csv', sep = '.'), row.names = F)
}
