# getting subset chromosome data
library(dplyr)
dat <- read_csv('SpeciesChromList.csv')
length(dat$Name)
tree <- read.tree('final_100trees')
tips <- tree[[1]]$tip.label
final <- c()
# record more than 1
count <- 0
for (i in 1:length(tips)){
  if (tips[i] %in% dat$Name){
    hitdat <- dat[which(dat$Name == tips[i]),c(1,4,6,8)]
    # de dup data
    colnames(hitdat) <- c("Family","Species","Chroms","SCS")
    hitdat <- hitdat %>%
      distinct( Chroms, SCS, .keep_all = TRUE)
    if (length(hitdat$Family) == 1){
      count = count + 1
    }
    final <-rbind(final, hitdat)
  }
  if (tips[i] %in% dat$Genus){
    hitdat <- dat[which(dat$Genus == tips[i]),c(1,4,6,8)]
    colnames(hitdat) <- c("Family","Species","Chroms","SCS")
    hitdat <- hitdat %>%
      distinct( Chroms, SCS, .keep_all = TRUE)
    if (length(hitdat$Family) == 1){
      count = count + 1
    }
    final <-rbind(final, hitdat)
  }
}
write_csv(final, file = "../data/sub_chrom.csv")
