# fusion and fision 
# chromosome number evo all dataset
library(coda)
chrom.number.results <- do.call(rbind,readRDS('../results/chrom_number_model_result.rds'))
plot(chrom.number.results$p[1:100], type = 'l',ylim = c(-410,-280), main = '', ylab = '')
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
     ylim =c(-10,300), xlim=c(0.002,0.015),)
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
x = 0.0133
points(x, 300, pch = 16, col = cols[1])
text(x, 300, pos = 4, labels = 'Fission', cex = 0.9)
points(x, 285, pch = 16, col = cols[2])
text(x, 285, pos = 4, labels = 'Fussion', cex = 0.9)
# PDF 6x6
