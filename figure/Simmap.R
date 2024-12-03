# simmap
source('../scripts/functions.R')
library(fields)
library(viridis)
library(ape)
library(chromePlus)

#####################
# chromosome number #
#####################

# randomly pick one tree 
tree <- read.tree('../data/final_100trees')[[37]]
# chrom data 
allchrom <- read.csv('../data/SpeciesChromList.csv')
# subset chrom 
tip.names <- tree$tip.label
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
#transform trees to mya from hundred of mya 
tree$edge.length <- tree$edge.length * 100
rng <- c(range(chrom$Chroms, na.rm = T)[1] - 1,
         range(chrom$Chroms, na.rm = T)[2] + 1)
chrom$gen.prob <- 1
for(i in 1:nrow(chrom)){
  if(chrom$SCS[i] %in% c('XY','NeoXY')){
    chrom$gen.prob[i] <- 0
  }
}
# make a Q matrix 
c <- rng[2]-rng[1]+1
Q <- matrix(data = 0, c,c)
# make row and column name vectors
cname<-c(rng[1]:rng[2])
colnames(Q)<-row.names(Q)<-cname
for (i in 1:c){
  for (j in 1:c){
    # define transition
    # number change
    if (i+1==j){
      Q[i,j] <- 1
    }
    if (i-1==j){
      Q[i,j] <- 2
    }
  }
}
#re order the tips
chrom$Species %in% tree$tip.label
chrom.s <- chrom
for (i in 1:length(tree$tip.label)){
  temp <- which(chrom$Species == tree$tip.label[i])
  chrom.s[i,] = chrom[temp,]
}
chrom.s$gen.prob <- 1
for(i in 1:nrow(chrom)){
  if(chrom.s$SCS[i] %in% c('XY','NeoXY')){
    chrom.s$gen.prob[i] <- 0
  }
}
chrom.mat_ <- datatoMatrix(x = chrom.s[,c(2,3,5)],
                           range = rng,
                           hyper = F)
results <- readRDS('../results/chrom_number_model_result.rds')
first_tree <- results[[37]]
b_result <- first_tree[51:100,]
model <- Q
Q[Q == 1] <- mean(b_result$asc1)
Q[Q == 2] <- mean(b_result$desc1)
diag(Q) <- -(rowSums(Q))
colnames(Q) <- rownames(Q) <- 1:19
colnames(model) <- rownames(model) <- 1:19
colnames(chrom.mat_) <- 1:19
test <- make.simmap2(tree = tree, x = chrom.mat_, model = model,Q = Q,nsim = 1,pi = "fitzjohn",rejmax = 1000000,rejint = 100000, monitor = T)
cols<-setNames(rev(viridis(n=19, option = 'H',begin = 0 )),
               c(1:19))
# fix simmap 
# only if there are rejections out of limits

dat.for.fix <- chrom.s[,c(2,3)]
dat.for.fix$Chroms <- dat.for.fix$Chroms -2
test.fixed <-fix.simmap(test,dat.for.fix,model)
test <- test.fixed[[1]]
#plot
plotSimmap(test, cols,fsize = .003, ftype = 'i',outline = F, lwd = 2, type = 'fan')
arc.cladelabels(node=439,text="Passalidae",offset=5,mark.node=FALSE)
arc.cladelabels(node=426,text="Lucanidae",offset=5,mark.node=FALSE)
arc.cladelabels(node=c(232,318),text="Scarabaeidae",offset=5,mark.node=FALSE)

# color bar
num_colors <- 19
values <- matrix(seq(3, 21, length.out = num_colors), ncol = 1)
plot(10, 10, type = "n", xlim = c(0, 2), ylim = c(0, 2), xlab = "", ylab = "")
image.plot(0, 1, values, col = colorRampPalette(rev(viridis(n=19, option = 'H',begin = 0 )))(num_colors), 
           axes = F, xlab = "", ylab = "", legend.only = T)

#########################
# sex chromosome system #
#########################

# checking divergence
results <- readRDS('../results/simple_model_scs.rds')
plot(results[[1]]$p, type = 'l', ylim = c(-130, -70))
for (i in 2:100){
  lines(results[[i]]$p)
}

# making simmap
tree <- read.tree('../data/final_100trees')[[37]]
results <- readRDS('../results/simple_model_scs.rds')
# chrom data 
allchrom <- read.csv('../data/SpeciesChromList.csv')
# subset chrom 
tip.names <- tree$tip.label
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
#transform trees to mya
tree$edge.length <- tree$edge.length * 100
rng <- c(range(chrom$Chroms, na.rm = T)[1] - 1,
         range(chrom$Chroms, na.rm = T)[2] + 1)
scs <- c()
for(i in 1:length(chrom$SCS)){
  hit <- which(chrom$Species == tree$tip.label[i])
  scs[i] <- chrom$SCS[hit]
}
names(scs) <- tree$tip.label
first_tree <- results[[37]]
b_result <- first_tree[51:100,]
Q<- matrix(c(0,3,5,1,0,6,2,4,0), 3)
colnames(Q) <- rownames(Q) <- c('NeoXY','XO','XY')
model <- Q
Q[1,2] <- mean(b_result$q23)
Q[1,3] <- mean(b_result$q21)
Q[2,1] <- mean(b_result$q32)
Q[2,3] <- mean(b_result$q31)
Q[3,1] <- mean(b_result$q12)
Q[3,2] <- mean(b_result$q13)
diag(Q) <- -rowSums(Q)
# make.simmap2 with the rate estimates from mcmc
test <- make.simmap2(tree=tree, x=scs ,Q=Q, pi="fitzjohn",nsim=1,model=model,monitor=T,rejmax=100000000)
cols <- setNames(viridis(3), c('NeoXY','XY','XO'))
plotSimmap(test, fsize =0.0003, type = 'fan', colors =cols)
arc.cladelabels(node=439,text="Passalidae",offset=5,mark.node=FALSE)
arc.cladelabels(node=426,text="Lucanidae",offset=5,mark.node=FALSE)
arc.cladelabels(node=c(232,318),text="Scarabaeidae",offset=5,mark.node=FALSE)
add.simmap.legend(leg =c('NeoXY','XY','XO'), colors = viridis(3), vertical = T, prompt = F, y =-105, x =130)
