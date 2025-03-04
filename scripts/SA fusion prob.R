# propotion of SA-fusion
library(evobiR)
library(ape)
library(phytools)
library(viridis)
source('functions.R')
library(coda)
# read in data
allchrom <- read.csv('../data/SpeciesChromList.csv')
trees <- read.tree("../data/final100trees")
tip.names <- trees[[1]]$tip.label
chrom <- data.frame()
# generate chrom data (random pick)
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
for(i in 1:length(trees)){
  trees[[i]]$edge.length <- trees[[i]]$edge.length * 100
}
# random pick one
# for plotting I chose 37 so keep it consistent
tree <- trees[[37]]
#making Q matrix
rng <- range(chrom$Chroms)
rng.len <- rng[2]-rng[1]+1
chrom.mat <- matrix(data = 0, nrow = (rng.len*3), ncol = (rng.len*3))
# fill in the Q matrix
# 12 parameters
# fill in +1 and -1 in same scs
for (i in 1:(rng.len*3)){
  for (j in 1:(rng.len*3)){
    # Neo to others
    if (i <= rng.len){
      # to Neo
      if (j <= rng.len){
        if (i +1 == j){
          chrom.mat[i,j] <- 1
        }
        if (i == j+1){
          chrom.mat[i,j] <- 2
        }
      }
      # to XO
      if (j > rng.len & j <= (rng.len*2)){
        if (i + rng.len == j){
          chrom.mat[i,j] <- 7
        }
      }
      # to XY
      if (j > (rng.len*2) & j <= (rng.len*3)){
        if (i + (rng.len*2) == j){
          chrom.mat[i,j] <- 8
        }
      }
    }
    # XO to other 
    if (i > rng.len & i <= (rng.len*2)){
      # to Neo
      if (j <= rng.len){
        if (i == j + rng.len +1){
          chrom.mat[i,j] <- 10
        }
      }
      # to XO
      if (j > rng.len & j <= (rng.len*2)){
        if ( i +1 == j){
          chrom.mat[i,j] <- 3
        }
        if ( i == j + 1){
          chrom.mat[i,j] <- 4
        }
      }
      # to XY 
      if (j > (rng.len*2)){
        if (i + (rng.len) == j +1){
          chrom.mat[i,j] <- 11
        }
      }
    }
    # XY to other
    if (i > (rng.len*2)){
      # to Neo
      if (j <= rng.len){
        if (i == j + (rng.len*2) + 1){
          chrom.mat[i,j] <- 12
        }
      }
      # to XO
      if (j > rng.len & j <= (rng.len*2)){
        if (j + rng.len == i) {
          chrom.mat[i,j] <- 9
        }
      }
      # to XY 
      if (j > (rng.len*2)){
        if (j == i+1 ){
          chrom.mat[i,j] <- 5
        }
        if (j +1 == i ){
          chrom.mat[i,j] <- 6
        }
      }
    }
  }
}
rownames(chrom.mat) <- colnames(chrom.mat) <- 1:51
# make sim state
chrom$sim.state <- NA
for (i in 1:length(chrom$SCS)){
  if (chrom$SCS[i] == 'NeoXY'){
    chrom$sim.state[i] <- chrom$Chroms[i] - rng[1] +1
  }
  if (chrom$SCS[i] == 'XO'){
    chrom$sim.state[i] <- chrom$Chroms[i] + (rng.len) - rng[1]+1
  }
  if (chrom$SCS[i] == 'XY'){
    chrom$sim.state[i] <- chrom$Chroms[i] + (rng.len*2) - rng[1]+1
  }
}
x <- matrix(0,  nrow=nrow(chrom), rng.len*3)
rownames(x) <- chrom$Species
for(i in 1:nrow(x)){
  x[i, chrom$sim.state[i]] <- 1
}
colnames(x) <- 1:(rng.len*3)
sim <- 100
test <- make.simmap2(tree, x = x, model = chrom.mat, pi = 'fitzjohn', nsim = sim,rejmax = 1000000,rejint = 100000, monitor = T )
# saveRDS(test, file = '../results/simmap.rds')
#
counts <- describe.simmap2(test)$count
## AA fusion (desc)
# column that are fusions
need.col <- c()
for (i in 2:rng.len){
  need.col <- c(need.col,paste(i,i-1, sep = ','))
}
for (i in 19:(rng.len*2)){
  need.col <- c(need.col,paste(i,i-1, sep = ','))
}
for (i in 36:(rng.len*3)){
  need.col <- c(need.col,paste(i,i-1, sep = ','))
}
AAfusion <- rowSums(counts[,which(colnames(counts) %in% need.col,T)])
## SAfusion
need.col <- c()
for (i in (rng.len+1+1):(rng.len*2)){
  need.col <- c(need.col,paste(i,i-18, sep = ','))
}
for (i in 36:51){
  need.col <- c(need.col,paste(i,i-35,sep = ','))
}
for (i in 19:34){
  need.col <- c(need.col,paste(i,i+16, sep = ','))
}
SAfusion <- rowSums(counts[,which(colnames(counts) %in% need.col,T)])
totalfusion <- rowSums(counts[,colnames(counts)])
obspropSA <- SAfusion / (SAfusion + AAfusion)
expSA <- c()
for(i in 1:100){
  print(i)
  times <- describe.simmap(test[[i]])$times[2, ]
  expSA[i] <- sum(Pfsa(Da =  6, scs = "XY") * times[c(1, 35)],
                  Pfsa(Da =  8, scs = "XY") * times[c(2, 36)],
                  Pfsa(Da = 10, scs = "XY") * times[c(3, 37)],
                  Pfsa(Da = 12, scs = "XY") * times[c(4, 38)],
                  Pfsa(Da = 14, scs = "XY") * times[c(5, 39)],
                  Pfsa(Da = 16, scs = "XY") * times[c(6, 40)],
                  Pfsa(Da = 18, scs = "XY") * times[c(7, 41)],
                  Pfsa(Da = 20, scs = "XY") * times[c(8, 42)],
                  Pfsa(Da = 22, scs = "XY") * times[c(9, 43)],
                  Pfsa(Da = 24, scs = "XY") * times[c(10, 44)],
                  Pfsa(Da = 26, scs = "XY") * times[c(11, 45)],
                  Pfsa(Da = 28, scs = "XY") * times[c(12, 46)],
                  Pfsa(Da = 30, scs = "XY") * times[c(13, 47)],
                  Pfsa(Da = 32, scs = "XY") * times[c(14, 48)],
                  Pfsa(Da = 36, scs = "XY") * times[c(15, 50)],
                  Pfsa(Da = 38, scs = "XY") * times[c(16, 51)],
                  #XO
                  Pfsa(Da =  7, scs = "XO") * times[18],
                  Pfsa(Da =  9, scs = "XO") * times[19],
                  Pfsa(Da = 11, scs = "XO") * times[20],
                  Pfsa(Da = 13, scs = "XO") * times[21],
                  Pfsa(Da = 15, scs = "XO") * times[22],
                  Pfsa(Da = 17, scs = "XO") * times[23],
                  Pfsa(Da = 19, scs = "XO") * times[24],
                  Pfsa(Da = 21, scs = "XO") * times[25],
                  Pfsa(Da = 23, scs = "XO") * times[26],
                  Pfsa(Da = 25, scs = "XO") * times[27],
                  Pfsa(Da = 27, scs = "XO") * times[28],
                  Pfsa(Da = 29, scs = "XO") * times[29],
                  Pfsa(Da = 31, scs = "XO") * times[30],
                  Pfsa(Da = 33, scs = "XO") * times[31],
                  Pfsa(Da = 35, scs = "XO") * times[32],
                  Pfsa(Da = 37, scs = "XO") * times[33],
                  Pfsa(Da = 39, scs = "XO") * times[34]
                  )
}
mean(expSA)
mean(obspropSA)
# plot
cols <- viridis(2, begin = 0.5, alpha = 0.65)
plot(density(expSA, bw = .01),
     xlim = c(.15, 0.65), main = "",
     xlab = "Proportion of Sex-Autosome Fusion")
polygon(density(expSA, bw = .01),
        col = cols[1])
lines(density(obspropSA))
polygon(density(obspropSA),
        col = cols[2])
x=0.58
points(x=x,y=40, pch =16, col = cols[1])
text(x=x,y=40, pos = 4, labels = "Expected", cex = 0.9)
hpd <- HPDinterval(as.mcmc(expSA))
lines(y=c(-0.5,-0.5), x=hpd[1:2], lwd=2,col=cols[1])
points(x=x,y=38, pch =16, col = cols[2])
text(x=x,y=38, pos = 4, labels = "Inferred", cex = 0.9)
hpd <- HPDinterval(as.mcmc(obspropSA))
lines(y=c(-0.5,-0.5), x=hpd[1:2], lwd=2,col=cols[2])
# save as PDF 6x6


## simmap plot
library(fields)
# NeoXY, XO, XY
cols <- setNames(c(colorRampPalette(c("firebrick1", "firebrick4"))(17),
                   colorRampPalette(c("green", "green4"))(17),
                   colorRampPalette(c("dodgerblue", "dodgerblue4"))(17)), 
                 c(1:51))
plotSimmap(test[[1]], type = 'fan', fsize = 0.03, colors = cols)
sexmode <- as.factor(setNames(chrom$SCS, chrom$Species))
sexmode <- to.matrix(sexmode, levels(sexmode))
sexmode <- sexmode[test[[1]]$tip.label,]
tiplabels(pie = sexmode, piecol = palette()[c(2,3,4)], cex = 0.15)
# Create a color palette
color_pal <- colorRampPalette(c("firebrick1", "firebrick4"))(17)
data_matrix <- matrix(1:17, nrow = 17, ncol = 1)
image.plot(z = data_matrix, col = color_pal, axes = FALSE, legend.only = TRUE,
           legend.width = 1,
           legend.shrink = 1, legend.mar = 3)

color_pal <- colorRampPalette(c("green", "green4"))(17)
data_matrix <- matrix(1:17, nrow = 17, ncol = 1)
image.plot(z = data_matrix, col = color_pal, axes = FALSE, legend.only = TRUE,
           legend.width = 1,
           legend.shrink = 1, add = TRUE, legend.mar = 6,
           axis.args = list(labels = FALSE, tick = FALSE))

color_pal <- colorRampPalette(c("dodgerblue", "dodgerblue4"))(17)
data_matrix <- matrix(1:17, nrow = 17, ncol = 1)
image.plot(z = data_matrix, col = color_pal, axes = FALSE, legend.only = TRUE,
           legend.width = 1,
           legend.shrink = 1, add = TRUE, legend.mar = 9,
           axis.args = list(labels = FALSE, tick = FALSE))