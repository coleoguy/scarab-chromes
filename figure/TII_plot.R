# Sean Chien
# TII plot
library(seqinr)
###TII by Mesquite
## bs500 
datMesquite<- read.csv("../data/TII_Mesquite.csv")
# plot
plot(x = (1:length(datMesquite$ID)),y = datMesquite$TII[order(datMesquite$TII)],
     xlab = "Taxa", ylab = "Taxonomic Instability Index")
th <- 200000
abline(h = th, col = 'red')
rogue_list <- datMesquite$ID[datMesquite$TII > th]
length(rogue_list)/length(datMesquite$ID)
# keep the outgroup
rogue_list <- rogue_list[!rogue_list == 'Enochrus falcarius']
