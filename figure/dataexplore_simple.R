library(coda)
library(viridis)
library(beeswarm)
library(ape)
library(chromePlus)
library(phytools)
library (diversitree)
library(plotrix)
# chromosome number evo all dataset
chrom.number.results <- do.call(rbind,readRDS('../results/chrom_number_model_result.rds'))
plot(chrom.number.results$p[1:100], type = 'l',ylim = c(-410,-310), main = '', ylab = '')
# checking convergence time 
for (i in 1:(length(chrom.number.results$p)/100-1)){
  index <- seq(101, length(chrom.number.results$p),by = 100)
  start <- index[i]
  end <- index[i]+100-1
  lines(chrom.number.results$p[start:end])
}
# post burn-in 
# asc
post.chrom.number.results <- chrom.number.results[chrom.number.results$i == c(51:100),]
# average 
mean(post.chrom.number.results$asc1)
mean(post.chrom.number.results$desc1)
cols <- viridis(2, option = 'D',alpha = 0.7, begin = 0.45)
plot(density((post.chrom.number.results$asc1)),main ='',xlab='Rate(/MY)',
     ylim =c(-10,255), xlim=c(0.002,0.017),)
polygon(density(post.chrom.number.results$asc1),col=cols[1])
hpd <- HPDinterval(as.mcmc(post.chrom.number.results$asc))
y <- -5
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=cols[1])
# desc
lines(density(post.chrom.number.results$desc1))
polygon(density(post.chrom.number.results$desc1),col=cols[2])
hpd <- HPDinterval(as.mcmc(post.chrom.number.results$desc))
y <- -10
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=cols[2])
points(0.015, 255, pch = 16, col = cols[1])
text(0.015, 255, pos = 4, labels = 'Fission', cex = 0.9)
points(0.015, 240, pch = 16, col = cols[2])
text(0.015, 240, pos = 4, labels = 'Fussion', cex = 0.9)
# PDF 6x6



# sub-tree
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

###############
# subset data #
###############

# post burnin
sub_luc <- all_luc[all_luc$i == c(51:100),]
sub_pass <- all_pass[all_pass$i == c(51:100),]
sub_sca <- all_sca[all_sca$i == c(51:100),]
### asc (fission) ###
# scarab
cols <- viridis(3, option = 'D',alpha = 0.7, begin = 0.45)
plot(density((sub_sca$asc1)),main ='',xlab='Fission (/MY)',
     ylim =c(-20,520), xlim=c(0,0.07),)
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
points(x=0.053,y=520, col = cols[1], pch= 16,cex=1)
text(x=0.053,y=520,pos= 4, "Scarabaeidae", cex=1)
points(x=0.053,y=490, col = cols[2], pch= 16,cex =1)
text(x=0.053,y=490,pos= 4, "Lucanidae", cex=1)
points(x=0.053,y=460, col = cols[3], pch= 16,cex=1)
text(x=0.053,y=460,pos= 4, "Passalidae", cex=1)
# save PDF 6x6

### desc ###
cols <- viridis(3, option = 'D',alpha = 0.7, begin = 0.45)
plot(density(sub_sca$desc1),main ='',xlab='Fusion (/MY)',
     xlim= c(0,0.09), ylim =c(-5,340))
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
points(x=0.07,y=340, col = cols[1], pch= 16)
text(x=0.07,y=340,pos= 4, "Scarabaeidae")
points(x=0.07,y=320, col = cols[2], pch= 16)
text(x=0.07,y=320,pos= 4, "Lucanidae")
points(x=0.07,y=300, col = cols[3], pch= 16)
text(x=0.07,y=300,pos= 4, "Passalidae")
# save as PDF 6x6 

# scs and number distribution
type_num_dat <- read.csv('../data/SpeciesChromList.csv')
df <- data.frame()
for (i in 1:length(type_num_dat$Family)){
  if (type_num_dat$Family[i] %in% c('Scarabaeidae', 'Passalidae','Lucanidae')){
    if (!type_num_dat$SCS[i] == ''){
      df <- rbind(df,type_num_dat[i,])
    }
  }
}
# remove XXXXXO
df <- df[-21,]
cols <- viridis(3, option = 'D',alpha = 0.7, begin = 0.45)
famcolor <- c('Scarabaeidae' = cols[1],"Lucanidae" = cols[2],"Passalidae" = cols[3])
beeswarm( df$autosome.haploid~df$SCS,
          method = c('center'),
         cex =2, pch = 16, spacing = 0.6,
         xlab = 'Sex Chromosome System', ylab = 'Haploid Autosome Number', pwcol = famcolor[df$Family],
         corral = c("random"),
         cex.lab = 1.2
         )
points(0.5,21, col = cols[1], pch = 16)
text(0.5,21, pos = 4, labels = 'Scarabaeidae')
points(0.5,20, col = cols[2], pch = 16)
text(0.5,20, pos = 4, labels = 'Lucanidae')
points(0.5,19, col = cols[3], pch = 16)
text(0.5,19, pos = 4, labels = 'Passalidae')
# save PDF 6x6

## plot the genus that have XY and Neo-XY and XY XO
## compare to their chromosome number change 
# chrom v sex chrom plot
dat <- read.csv("../data/SpeciesChromList.csv")
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
cols <- viridis(7, option = 'D',alpha = 0.7, begin = 0)
plot(y=pdat$XY, x=xs1,
     xlim=c(.75,2.25),
     ylim=c(4.5,9.5),
     xaxt="n",
     xlab="Sex Chromosome System", ylab = 'Haploid Autosome Number ( x\u0305 )',
     cex =1.5, pch = 16,
     cex.lab = 1.2,
     col = cols
     )
cols <- viridis(7, option = 'D',alpha = 0.7, begin = 0)
x <- 1.98
y <- 9.6
y.dis <- 0.17
cex <- 0.9
text(x=x,y=y, "Phanaeus", adj=c(0,0.5), cex = cex)
text(x=x,y=y-y.dis, "Dorcus", adj=c(0,0.5), cex = cex)
text(x=x,y=y-2*y.dis, "Oryctes", adj=c(0,0.5), cex = cex)
text(x=x,y=y-3*y.dis, "Deltochilum", adj=c(0,0.5), cex = cex)
text(x=x,y=y-4*y.dis, "Haplidia", adj=c(0,0.5), cex = cex)
text(x=x,y=y-5*y.dis, "Phileurus", adj=c(0,0.5), cex = cex)
text(x=x,y=y-6*y.dis, "Phyllognathus",adj=c(0,0.5), cex = cex)
x <- 1.9
x2 <- 1.97
lwd <- 2
lines(x=c(x,x2), y=rep(y,2), col = cols[1],lwd = lwd)
lines(x=c(x,x2), y=rep(y-y.dis,2), col = cols[2],lwd = lwd, lty =3)
lines(x=c(x,x2), y=rep(y-2*y.dis,2), col = cols[3],lwd = lwd)
lines(x=c(x,x2), y=rep(y-3*y.dis,2), col = cols[4],lwd = lwd)
lines(x=c(x,x2), y=rep(y-4*y.dis,2), col = cols[5],lwd = lwd)
lines(x=c(x,x2), y=rep(y-5*y.dis,2), col = cols[6],lwd = lwd)
lines(x=c(x,x2), y=rep(y-6*y.dis,2), col = cols[7],lwd = lwd)
axis(side=1,at=c(1,2), c("XY","NeoXY"))
xs2 <- c(2, 1+xs1[4], 1+xs1[5],2,1+xs1[6],2,1+xs1[7])
points(y=pdat$NeoXY, x=xs2, cex=1.5, col = cols, pch =16)
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
# save PDF 6x6

