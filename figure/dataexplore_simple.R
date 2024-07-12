library(coda)
library(viridis)


library(chromePlus)
library(ape)
library(beeswarm)
library(phytools)
library (diversitree)
source('functions.R')
library(fields)
library(plotrix)

## checking convergence
all_luc <- read.csv('../results/sub_luc.csv')
all_pass <- read.csv('../results/sub_pas.csv')
all_sca <- read.csv('../results/sub_sca.csv')

plot(all_luc$p[1:100], type = 'l', ylim = c(-60,-35), main = 'Lucanidae', ylab = '')
for (i in 1:(length(all_luc$p)/100-1)){
  index <- seq(101, length(all_luc$p),by = 100)
  start <- index[i]
  end <- index[i]+100-1
  lines(all_luc$p[start:end])
}

plot(all_pass$p[1:100], type = 'l', ylim = c(-120, -60), main = 'Passalidae', ylab = '')
for (i in 1:(length(all_pass$p)/100-1)){
  index <- seq(101, length(all_pass$p),by = 100)
  start <- index[i]
  end <- index[i]+100-1
  lines(all_pass$p[start:end])
}

plot(all_sca$p[1:100], type = 'l', ylim = c(-210,-165), main = 'Scarabeidae', ylab = '')
for (i in 1:(length(all_sca$p)/100-1)){
  index <- seq(101, length(all_sca$p),by = 100)
  start <- index[i]
  end <- index[i]+100-1
  lines(all_sca$p[start:end])
}

## subset data 
## post burnin

sub_luc <- all_luc[all_luc$i == c(51:100),]
sub_pass <- all_pass[all_pass$i == c(51:100),]
sub_sca <- all_sca[all_sca$i == c(51:100),]

# asc (fission)
# scarab
cols <- viridis(3, option = 'D',alpha = 0.7, begin = 0.45)
plot(density((sub_sca$asc1)),main ='',xlab='Fission (/MY)',
     ylim =c(-20,700), xlim=c(0,0.07),)
polygon(density(sub_sca$asc1),col=cols[1])
hpd <- HPDinterval(as.mcmc(sub_sca$asc))
y <- -7
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=cols[1])
# lucanidae
lines(density(sub_luc$asc1))
polygon(density(sub_luc$asc1),col=cols[2])
hpd <- HPDinterval(as.mcmc(sub_luc$asc))
y <- -15
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=cols[2])
# passalidae
lines(density((sub_pass$asc1)))
polygon(density(sub_pass$asc1),col=cols[3])
hpd <- HPDinterval(as.mcmc(sub_pass$asc))
y<- -22
cols <- viridis(3, option = 'D',alpha = 1, begin = 0.45)
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=cols[3])

points(x=0.053,y=700, col = cols[1], pch= 16,cex=1)
text(x=0.053,y=700,pos= 4, "Scarabaeidae", cex=1)
points(x=0.053,y=660, col = cols[2], pch= 16,cex =1)
text(x=0.053,y=660,pos= 4, "Lucanidae", cex=1)
points(x=0.053,y=620, col = cols[3], pch= 16,cex=1)
text(x=0.053,y=620,pos= 4, "Passalidae", cex=1)
# save PDF 6x6

#desc
cols <- viridis(3, option = 'D',alpha = 0.7, begin = 0.45)
plot(density(sub_sca$desc1),main ='',xlab='Fusion (/MY)',
     xlim= c(0,0.09), ylim =c(-5,380))
polygon(density(sub_sca$desc1),col=cols[1])
hpd <- HPDinterval(as.mcmc(sub_sca$desc))
y <- -5
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=cols[1])
lines(density((sub_luc$desc1)))
polygon(density(sub_luc$desc1),col=cols[2])
hpd <- HPDinterval(as.mcmc(sub_luc$desc))
y <- -10
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=cols[2])
lines(density((sub_pass$desc1)))
polygon(density(sub_pass$desc1),col=cols[3])
hpd <- HPDinterval(as.mcmc(sub_pass$desc))
y <-  -15
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=cols[3])
cols <- viridis(3, option = 'D',alpha = 1, begin = 0.45)
points(x=0.07,y=370, col = cols[1], pch= 16)
text(x=0.07,y=370,pos= 4, "Scarabaeidae")
points(x=0.07,y=350, col = cols[2], pch= 16)
text(x=0.07,y=350,pos= 4, "Lucanidae")
points(x=0.07,y=330, col = cols[3], pch= 16)
text(x=0.07,y=330,pos= 4, "Passalidae")
# save as PDF 6x6 


# compare chromosome number and chromosome type
type_num_dat <- read.csv('../data/chrom.data/SpeciesChromList.csv')
df <- data.frame()
for (i in 1:length(type_num_dat$Family)){
  if (type_num_dat$Family[i] %in% c('Scarabaeidae', 'Passalidae','Lucanidae')){
    if (!type_num_dat$SCS[i] == ''){
      df <- rbind(df,type_num_dat[i,])
    }
  }
}
# XXXXXO
df <- df[-21,]
cols <- viridis(3, option = 'D',alpha = 0.7, begin = 0.45)
famcolor <- c('Scarabaeidae' = cols[1],"Lucanidae" = cols[2],"Passalidae" = cols[3])
beeswarm( df$autosome.haploid~df$SCS,
          method = c('center'),
         cex =2, pch = 16, spacing = 0.6,
         xlab = 'Sex Chromosome System', ylab = 'Haploid Autosome Number', pwcol = famcolor[df$Family],
         corral = c("random")
         )
legend("topleft", legend = c( 'Scarabaeidae',"Lucanidae","Passalidae"),
       col = cols, pch = 19, cex = 1,bty = "n")
## plot the genus that have XY and Neo-XY and XY XO
## compare to their chromosome number change 
# chrom v sex chrom plot
dat <- read.csv("../data/chrom.data/SpeciesChromList.csv")
dat <- dat[!is.na(dat$autosome.haploid),]
gen <- unique(dat$Genus)
XY <- NeoXY <- rep(NA,length(gen))
pdat <- data.frame(gen,XY,NeoXY)
for(i in 1:nrow(pdat)){
  for(j in 2:3){
    pdat[i,j] <- mean(dat$autosome.haploid[dat$Genus == pdat$gen[i] & 
                                             dat$SCS == colnames(pdat)[j]])
  }
}
pdat <- pdat[complete.cases(pdat),]
pdat <- pdat[order(pdat$XY),]
xs1 <- c(1,1,1, seq(from=.95, by=.03,length.out=4))

plot(y=pdat$XY, x=xs1,
     xlim=c(.75,2.25),
     ylim=c(4.5,9.5),
     xaxt="n",
     xlab="Sex Chromosome System", ylab = 'Haploid Autosome Number ( x\u0305 )',
     cex =1.5
     )
cols <- viridis(7, option = 'D',alpha = 0.7, begin = 0)
text(x=1.98,y=9.6, "Phanaeus", adj=c(0,0.5))
text(x=1.98,y=9.45, "Dorcus", adj=c(0,0.5))
text(x=1.98,y=9.3, "Oryctes", adj=c(0,0.5))
text(x=1.98,y=9.15, "Deltochilum", adj=c(0,0.5))
text(x=1.98,y=9, "Haplidia", adj=c(0,0.5))
text(x=1.98,y=8.85, "Phileurus", adj=c(0,0.5))
text(x=1.98,y=8.7, "Phyllognathus",adj=c(0,0.5))

lines(x=c(1.87,1.95), y=rep(9.6,2), col = cols[1],lwd = 3)
lines(x=c(1.87,1.95), y=rep(9.45,2), col = cols[2],lwd = 3, lty =3)
lines(x=c(1.87,1.95), y=rep(9.3,2), col = cols[3],lwd = 3)
lines(x=c(1.87,1.95), y=rep(9.15,2), col = cols[4],lwd = 3)
lines(x=c(1.87,1.95), y=rep(9,2), col = cols[5],lwd = 3)
lines(x=c(1.87,1.95), y=rep(8.85,2), col = cols[6],lwd = 3)
lines(x=c(1.87,1.95), y=rep(8.7,2), col = cols[7],lwd = 3)
axis(side=1,at=c(1,2), c("XY","NeoXY"))
xs2 <- c(2, 1+xs1[4], 1+xs1[5],2,1+xs1[6],2,1+xs1[7])
points(y=pdat$NeoXY, x=xs2, cex=1.5)

# cols[2] <- "red"
for(i in 1:7){
  if (i != 2){
    lines(x=c(xs1[i],xs2[i]), y=c(pdat$XY[i],pdat$NeoXY[i]),
          col=cols[i], lwd = 3)
  }
  if (i == 2){
    lines(x=c(xs1[i],xs2[i]), y=c(pdat$XY[i],pdat$NeoXY[i]),
          col=cols[i], lwd = 3, lty =3)
  }
}

###chrom number###
#make simmap (only chrom number)
tree <- read.tree('../data/final_100trees')[[37]]
chrom <- read.csv('../data/final_chrom.csv')
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
results <- readRDS('../results/simple_model_scs.rds')
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
# fix simmap
dat.for.fix <- chrom.s[,c(1,2)]
dat.for.fix$Chroms <- dat.for.fix$Chroms -2
test.fixed <-fix.simmap(test,dat.for.fix,model)
# write.simmap(test.fixed[[1]], file = '../results/simmap_chrom_num', map.order = 'right-to-left')
test <- read.simmap(file = '../results/simmap_chrom_num',format = 'phylip')
cols<-setNames(rev(viridis(n=19, option = 'H',begin = 0 )),
               c(1:19))
plotSimmap(test,cols,fsize = 0.003, ftype = 'i',outline = F, lwd = 2, type = 'fan')

plotSimmap(test,cols, fsize = .2, ftype = 'i',outline = F, lwd = 1)
ape::nodelabels(cex = 0.1)
# drop_test$tip.label
# ape::nodelabels(cex = 0.1)
# need to drop 
test$tip.label
whole_dat <- read.csv('../data/chrom.data/SpeciesChromList.csv')
fam.col <- c()
for (i in 1:length(test$tip.label)){
  if (test$tip.label[i] %in% whole_dat$Name){
    fam.col <- c(fam.col,whole_dat$Family[which(test$tip.label[i] == whole_dat$Name)[1]])
  }
  if (test$tip.label[i] %in% whole_dat$Genus){
    fam.col <- c(fam.col,whole_dat$Family[which(test$tip.label[i] == whole_dat$Genus)[1]])
  }
  
}
test2 <- test
test2$tip.label <- fam.col
plotSimmap(test2, cols,fsize = .5, ftype = 'i',outline = F, lwd = 2, type = 'fan')
plotSimmap(test, cols,fsize = .5, ftype = 'i',outline = F, lwd = 2, type = 'fan')
ape::nodelabels(cex = 0.1)
# drop those tips 
# Diplotaxis obscura
# Autoserica
# Maladera
# Lichnanthe rathvoni
# Amphimallon majale
# Phyllophaga delata
# Phyllophaga fusca
# Phyllophaga pleei
# Cotalpa lanigera
# Glaresis
# Homaloplia
# Pelidnota punctata
# Ectinohoplia rufipes

#plot
plotSimmap(test, cols,fsize = .003, ftype = 'i',outline = F, lwd = 2, type = 'fan')
arc.cladelabels(node=422,text="Passalidae",offset=5,mark.node=FALSE)
arc.cladelabels(node=434,text="Lucanidae",offset=5,mark.node=FALSE)
arc.cladelabels(node=c(238,325),text="Scarabaeidae",offset=5,mark.node=FALSE)




# color bar
num_colors <- 19
values <- matrix(seq(3, 21, length.out = num_colors), ncol = 1)
plot(10, 10, type = "n", xlim = c(0, 2), ylim = c(0, 2), xlab = "", ylab = "")
image.plot(0, 1, values, col = colorRampPalette(rev(viridis(n=19, option = 'H',begin = 0 )))(num_colors), 
           axes = F, xlab = "", ylab = "", legend.only = T)

###########
### SCS ###
###########

# checking divergence
results <- readRDS('../results/simple_model_scs.rds')
plot(results[[1]]$p, type = 'l', ylim = c(-110, -70))
for (i in 2:100){
  lines(results[[i]]$p)
}
chrom <- read.csv('../data/final_chrom.csv')
trees <- read.tree("../data/final_100trees")
#transform trees to mya from hundred of mya 
for(i in 1:length(trees)){
  trees[[i]]$edge.length <- trees[[i]]$edge.length * 100
}
scs <- c()
for(i in 1:length(chrom$SCS)){
  hit <- which(chrom$Species == trees[[37]]$tip.label[i])
  scs[i] <- chrom$SCS[hit]
}
names(scs) <- trees[[37]]$tip.label
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

# transition between 3 different sex chromosome systems are differnt

# makeing simmap with ARD model
test <- make.simmap(trees[[37]], x=scs, model = 'ARD' , pi="fitzjohn")
# makeing simmap with the rate estimates from mcmc
test <- make.simmap(tree=trees[[37]], x=scs ,Q=Q, pi="fitzjohn",nsim=100,model=model)
# make.simmap2 with the rate estimates from mcmc
test <- make.simmap2(tree=trees[[37]], x=scs ,Q=Q, pi="fitzjohn",nsim=1,model=model,monitor=T,rejmax=100000000)

# the issue might have be the low rate and transition between short branch length
# so I did plot the tree and sex chromosome systems and drop some suspecious tips 
plotSimmap(test,fsize = 0.003, ftype = 'i',outline = F, lwd = 2, type = 'fan')

cols <- setNames(viridis(3), c('NeoXY','XY','XO'))
plotSimmap(test, fsize =0.0003, type = 'fan', colors =cols)
target_node <- 496  # Replace with the ID of your target node
descendant <- getDescendants(test, target_node)
drop_tip <- test$tip.label[descendant]
target_node <- 499  # Replace with the ID of your target node
descendant <- getDescendants(test, target_node)
drop_tip <- c(drop_tip,test$tip.label[descendant])
drop_test <- drop.tip.simmap(test, drop_tip)
target_node <- 463  # Replace with the ID of your target node
descendant <- getDescendants(drop_test, target_node)
drop_tip <- drop_test$tip.label[descendant]
target_node <- 459  # Replace with the ID of your target node
descendant <- getDescendants(drop_test, target_node)
drop_tip <- c(drop_tip,test$tip.label[descendant])
drop_test <- drop.tip.simmap(drop_test, drop_tip)
#plot
plotSimmap(drop_test, cols,fsize = .003, ftype = 'i',outline = F, lwd = 2, type = 'fan')
arc.cladelabels(node=449,text="Passalidae",offset=5,mark.node=FALSE)
arc.cladelabels(node=434,text="Lucanidae",offset=5,mark.node=FALSE)
arc.cladelabels(node=c(238,325),text="Scarabaeidae",offset=5,mark.node=FALSE)
add.simmap.legend(leg =c('NeoXY','XY','XO'), colors = viridis(3), vertical = T, prompt = F, y =-1.1, x =1.3)

