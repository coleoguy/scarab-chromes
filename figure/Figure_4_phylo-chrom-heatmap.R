# genus tree
library(ape)
library(phytools)
library(viridis)
trees <- read.tree('../data/final100trees')
chrom <- read.csv('../data/SpeciesChromList.csv')
# random pick one tree
set.seed(4)
tree <- trees[[sample(100,1)]]
# trimming genus level
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

fam <- c()
for (i in 1:length(trim.genus.tree$tip.label)){
  if (trim.genus.tree$tip.label[i] %in% chrom$Genus){
    fam <- c(fam, chrom$Family[which(trim.genus.tree$tip.label[i] == chrom$Genus)][1])
  }
}

# reorder tip 
tip_order <- trim.genus.tree$tip.label[order(fam)]
tree_ordered <- rotateConstr(trim.genus.tree, tip_order)


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
cols <- viridis(40, direction = -1, option = "F")
cols[1] <- 'white'
phylo.heatmap(tree_ordered, table, fsiz = c(0.0022,0.5,0.5),
              split=c(0.7,0.3),
              legend = T,labels = T,
              legend.pos = c(0.02, 0.01),
              pts = F, lwd = 0.5,
              colors = cols)
famcols <- viridis(3, option = 'D', begin = 0.45)
fam <- c(rep("Lucanidae",8),rep("Passalidae",12),rep("Scarabaeidae",93))
fam <- gsub("Scarabaeidae", famcols[1], fam)
fam <- gsub("Passalidae", famcols[2], fam)
fam <- gsub("Lucanidae", famcols[3],fam)
# Add colored tip labels manually
tiplabels(pch = 19, col = fam, cex = 0.5, adj = 0.5)
x = -0.3
y = -0.015
s = 0.035
points(x,y, pch = 19, col = famcols[1])
text(x,y, label = "Scarabaeidae", pos =4)
points(x,y-s,pch = 19, col = famcols[2])
text(x,y-s, label = "Passalidae", pos =4)
points(x,y-2*s,pch = 19, col = famcols[3])
text(x,y-2*s, label = "Lucanidae", pos =4)


# draw grid line
lwd = 0.05
grid.col <- 'gray87'
y0=-0.0045
y1=1.0045
x0 = 1.1424
x1 = 0.501
# find the limit first
# segments(x0 = x1, y0 = y0, x1 = x1, y1 = y1, lwd = lwd, col = grid.col)
# segments(x0 = x0, y0 = y0, x1 = x0, y1 = y1, lwd = lwd, col = grid.col)
# segments(x0 = 0.5, y0 = y0, x1 = 1.143, y1 = y0, lwd = lwd, col = grid.col)
# segments(x0 = 0.5, y0 = y1, x1 = 1.143, y1 = y1, lwd = lwd, col = grid.col)

for (i in 1:(length(trim.genus.tree_ordered$tip.label)+1)){
  z = seq(from = y0, to = y1, length.out = length(trim.genus.tree_ordered$tip.label)+1)[i]
  segments(x0 = x1, y0 = z, x1 = x0, y1 = z, lwd = lwd, col = grid.col)
}

for (i in 1:(ncol(table)+1)){
  z = seq(from = x0, to = x1, length.out = ncol(table)+1)[i]
  segments(x0 = z, y0 = y0, x1 = z, y1 = y1, lwd = lwd,col = grid.col)
}
# export 600x600

