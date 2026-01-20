library(viridis)
library(beeswarm)
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
cex = 2
beeswarm( df$autosome.haploid~df$SCS,
          method = c('center'),
          cex = cex , pch = 16, spacing = 0.6,
          xlab = 'Sex Chromosome System', ylab = 'Haploid Autosome Number', pwcol = famcolor[df$Family],
          corral = c("random"),
          cex.lab = 1.2
)
cex = 1.5
cols <- viridis(3, option = 'D', begin = 0.45)
points(0.5,21, col = cols[1], pch = 16, cex = cex )
text(0.5,21, pos = 4, labels = 'Scarabaeidae')
points(0.5,20, col = cols[2], pch = 16, cex = cex )
text(0.5,20, pos = 4, labels = 'Lucanidae')
points(0.5,19, col = cols[3], pch = 16, cex = cex )
text(0.5,19, pos = 4, labels = 'Passalidae')
title(main = "(A)", adj = 0, line = 0.5)
# save PDF 6x6