## plot the genus that have XY and Neo-XY and XY XO
## compare to their chromosome number change 
# chrom v sex chrom plot
library(viridis)
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
title(main = "(B)", adj = 0, line = 0.5)
# save PDF 6x6