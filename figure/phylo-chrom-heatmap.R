# genus tree
library(ape)
library(phytools)
library(viridis)
library(RColorBrewer)
trees <- read.tree('../data/final100trees')
chrom <- read.csv('../data/SpeciesChromList.csv')
# random pick one tree
set.seed(6)
tree <- trees[[sample(100,1)]]
genus <- c()
for (i in 1:length(tree$tip.label)){
  # full name
  if (tree$tip.label[i] %in% chrom$Name){
    genus <- c(genus, chrom$Genus[which(tree$tip.label[i] == chrom$Name)][1])
  }
  # genus name
  if (tree$tip.label[i] %in% chrom$Genus){
    genus <- c(genus, chrom$Genus[which(tree$tip.label[i] == chrom$Genus)][1])
  }
}
genus <- unique(genus)
# random pick tips 
tips <- c()
for (i in 1:length((genus))){
  # randome sample from all hit
  tips <- c(tips, tree$tip.label[grepl(genus[i], tree$tip.label)][sample(sum(grepl(genus[i], tree$tip.label)),1)])
}
# trim tree
trim.genus.tree <- keep.tip(tree, tip = tips)
# change tip label only genus level
for (i in 1:length(trim.genus.tree$tip.label)){
  trim.genus.tree$tip.label[i] <- sub("_.*", "", trim.genus.tree$tip.label[i])
    
}
# creating data table
table <- matrix(data = 0, nrow = length(tips), ncol = (range(chrom$male.haploid)[2]-range(chrom$male.haploid)[1]+1))
table <- data.frame(table)
row.names(table) <- genus
colnames(table) <- range(chrom$male.haploid)[1]:range(chrom$male.haploid)[2]
for (i in 1:length(genus)){
  for (j in 1:length(chrom$male.haploid[which(chrom$Genus == genus[i])])){
    num <- 0
    num <- chrom$male.haploid[which(chrom$Genus == genus[i])][j]
    table[i,num - range(chrom$male.haploid)[1] +1] <- table[i,num - range(chrom$male.haploid)[1] +1] +1
  }
}
# remove the last two row because of missing data
table <- table[,-c(19,18)]
# reorder tree 
phylo.heatmap(trim.genus.tree, table, fsiz = c(0.0025,0.5,0.4),
              split=c(0.7,0.3),
              legend = T,labels = T,
              pts = F, lwd = 0.5,
              #colors = heat.colors(n= 300)[300:1])
              colors = viridis(300, direction = -1, option = "F"))
