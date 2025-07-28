library(viridis)
library(coda)
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

plot(all_pass$p[1:100], type = 'l', ylim = c(-120, -50), main = 'Passalidae', ylab = '')
for (i in 1:(length(all_pass$p)/100-1)){
  index <- seq(101, length(all_pass$p),by = 100)
  start <- index[i]
  end <- index[i]+100-1
  lines(all_pass$p[start:end])
}

plot(all_sca$p[1:100], type = 'l', ylim = c(-270,-150), main = 'Scarabeidae', ylab = '')
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

### desc ###
cols <- viridis(3, option = 'D',alpha = 0.7, begin = 0.45)
hpdcols <- viridis(3, option = 'D', begin = 0.45)
plot(density(sub_sca$desc1),main ='',xlab='Fusion (MY)',
     xlim= c(0,0.09), ylim =c(-10,260))
title(main = "(A)", adj = 0, line = 0.5)
polygon(density(sub_sca$desc1),col=cols[1])
hpd <- HPDinterval(as.mcmc(sub_sca$desc))
y <- -5
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=hpdcols[1])
lines(density((sub_luc$desc1)))
polygon(density(sub_luc$desc1),col=cols[2])
hpd <- HPDinterval(as.mcmc(sub_luc$desc))
y <- -10
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=hpdcols[2])
lines(density((sub_pass$desc1)))
polygon(density(sub_pass$desc1),col=cols[3])
hpd <- HPDinterval(as.mcmc(sub_pass$desc))
y <-  -15
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=hpdcols[3])
x = 0.07
cex = 1.2 
points(x=x,y=260, col = hpdcols[1], pch= 16, cex = cex )
text(x=x,y=260,pos= 4, "Scarabaeidae")
points(x=x,y=245, col = hpdcols[2], pch= 16, cex = cex)
text(x=x,y=245,pos= 4, "Lucanidae")
points(x=x,y=230, col = hpdcols[3], pch= 16, cex = cex)
text(x=x,y=230,pos= 4, "Passalidae")
# save as PDF 6x6 

### asc (fission) ###
# scarab
cols <- viridis(3, option = 'D',alpha = 0.7, begin = 0.45)
hpdcols <- viridis(3, option = 'D', begin = 0.45)
plot(density((sub_sca$asc1)),main ='',xlab='Fission (MY)',
     ylim =c(-10,260), xlim=c(0,0.09))
title(main = "(B)", adj = 0, line = 0.5)
polygon(density(sub_sca$asc1),col=cols[1])
hpd <- HPDinterval(as.mcmc(sub_sca$asc))
y <- -5
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=hpdcols[1])
# lucanidae
lines(density(sub_luc$asc1))
polygon(density(sub_luc$asc1),col=cols[2])
hpd <- HPDinterval(as.mcmc(sub_luc$asc))
y <- -10
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=hpdcols[2])
# passalidae
lines(density((sub_pass$asc1)))
polygon(density(sub_pass$asc1),col=cols[3])
hpd <- HPDinterval(as.mcmc(sub_pass$asc))
y<- -15
lines(y=c(y,y), x=hpd[1:2], lwd=2,col=hpdcols[3])
x = 0.07
cex = 1.2
points(x=x,y=260, col = hpdcols[1], pch= 16,cex = cex )
text(x=x,y=260,pos= 4, "Scarabaeidae", cex=1)
points(x=x,y=245, col = hpdcols[2], pch= 16,cex = cex )
text(x=x,y=245,pos= 4, "Lucanidae", cex=1)
points(x=x,y=230, col = hpdcols[3], pch= 16,cex = cex )
text(x=x,y=230,pos= 4, "Passalidae", cex=1)
# save PDF 6x6

