################################################################################
# Genes in cluters of biological processes
################################################################################
# State graphical variables
arco <- rainbow(9)
palette(arco)
par(pch = 16)

#genes de cada pathway
n <- NULL
for (i in 1:length(PDS$genesinpathway)){
  ng <- length(PDS$genesinpathway[[i]])
  n <- c(n, ng)
}
maxn <- max(n) # Maximo número de genes en pathways

#Matrix of genes per pathway
gip <- NULL
for (i in 1:length(PDS$genesinpathway)){
  z <- PDS$genesinpathway[[i]]
  length(z) <- maxn
  gip <- rbind(gip, z)
}
row.names(gip) <- row.names(PDSmatrix)

## Frequencia de genes en grupos de vías
phi <- 1.61803398
#BPG1 genes
npaths <- nBP[myrw == 7] # Numeros de pathways pertenecientes al grupo
gcluster <- sort(gip[npaths,])
gcluster <- gsub(" ", "", gcluster, fixed = T)
tableBPG1 <- sort(table(gcluster), decreasing = T)
tableBPG1 <- tableBPG1/length(npaths)
jpeg("Bio_Process_G1_genesfreq.jpg",
     width = phi*600,
     height = 600,
     units = "px"
     
)
barplot(tableBPG1, col = 1, main = "Frecuencia de genes en cluster", ylab = "Frecuencia", 
        xlab = "Genes" )
dev.off()


#BPG2 genes
npaths <- nBP[myrw == 9] # Numeros de pathways pertenecientes al grupo
length(npaths)
gcluster <- sort(gip[npaths,])
gcluster <- gsub(" ", "", gcluster, fixed = T)
tableBPG2 <- sort(table(gcluster), decreasing = T)
tableBPG2 <- tableBPG2/length(npaths)
jpeg("Bio_Process_G2_genesfreq.jpg",
     width = phi*600,
     height = 600,
     units = "px"
     
)
barplot(tableBPG2, col = 2, main = "Frecuencia de genes en cluster", ylab = "Frecuencia", 
        xlab = "Genes")
dev.off()

#BPG3 genes
npaths <- nBP[myrw == 1] # Numeros de pathways pertenecientes al grupo
length(npaths)
gcluster <- sort(gip[npaths,])
gcluster <- gsub(" ", "", gcluster, fixed = T)
tableBPG3 <- sort(table(gcluster), decreasing = T)
tableBPG3 <- tableBPG3/length(npaths)
jpeg("Bio_Process_G3_genesfreq.jpg",
     width = phi*600,
     height = 600,
     units = "px"
     
)
barplot(tableBPG3, col = 3, main = "Frecuencia de genes en cluster", ylab = "Frecuencia", 
        xlab = "Genes")
dev.off()

#BPG4 genes
npaths <- nBP[myrw == 3] # Numeros de pathways pertenecientes al grupo
length(npaths)
gcluster <- sort(gip[npaths,])
gcluster <- gsub(" ", "", gcluster, fixed = T)
tableBPG4 <- sort(table(gcluster), decreasing = T)
tableBPG4 <- tableBPG4/length(npaths)
jpeg("Bio_Process_G4_genesfreq.jpg",
     width = phi*600,
     height = 600,
     units = "px"
     
)
barplot(tableBPG4, col = 4, main = "Frecuencia de genes en cluster", ylab = "Frecuencia", 
        xlab = "Genes")
dev.off()

#BPG5 genes
npaths <- nBP[myrw == 5] # Numeros de pathways pertenecientes al grupo
length(npaths)
gcluster <- sort(gip[npaths,])
gcluster <- gsub(" ", "", gcluster, fixed = T)
tableBPG5 <- sort(table(gcluster), decreasing = T)
tableBPG5 <- tableBPG5/length(npaths)
jpeg("Bio_Process_G5_genesfreq.jpg",
     width = phi*600,
     height = 600,
     units = "px"
     
)
barplot(tableBPG5, col = 5, main = "Frecuencia de genes en cluster", ylab = "Frecuencia", 
        xlab = "Genes" )
dev.off()

#BPG6 genes
npaths <- nBP[myrw == 2] # Numeros de pathways pertenecientes al grupo
length(npaths)
gcluster <- sort(gip[npaths,])
gcluster <- gsub(" ", "", gcluster, fixed = T)
tableBPG6 <- sort(table(gcluster), decreasing = T)
tableBPG6 <- tableBPG6/length(npaths)
jpeg("Bio_Process_G6_genesfreq.jpg",
     width = phi*600,
     height = 600,
     units = "px"
     
)
barplot(tableBPG6, col = 6, main = "Frecuencia de genes en cluster", ylab = "Frecuencia", 
        xlab = "Genes")
dev.off()

#BPG7 genes
npaths <- nBP[myrw == 8] # Numeros de pathways pertenecientes al grupo
length(npaths)
gcluster <- sort(gip[npaths,])
gcluster <- gsub(" ", "", gcluster, fixed = T)
tableBPG7 <- sort(table(gcluster), decreasing = T)
tableBPG7 <- tableBPG7/length(npaths)
jpeg("Bio_Process_G7_genesfreq.jpg",
     width = phi*600,
     height = 600,
     units = "px"
     
)
barplot(tableBPG7, col = 7, main = "Frecuencia de genes en cluster", ylab = "Frecuencia", 
        xlab = "Genes")
dev.off()


#BPG8 genes
npaths <- nBP[myrw == 4] # Numeros de pathways pertenecientes al grupo
length(npaths)
gcluster <- sort(gip[npaths,])
gcluster <- gsub(" ", "", gcluster, fixed = T)
tableBPG8 <- sort(table(gcluster), decreasing = T)
tableBPG8 <- tableBPG8/length(npaths)
jpeg("Bio_Process_G8_genesfreq.jpg",
     width = phi*600,
     height = 600,
     units = "px"
     
)
barplot(tableBPG8, col = 8, main = "Frecuencia de genes en cluster", ylab = "Frecuencia", 
        xlab = "Genes")
dev.off()

#BPG9 genes
npaths <- nBP[myrw == 6] # Numeros de pathways pertenecientes al grupo
length(npaths)
gcluster <- sort(gip[npaths,])
gcluster <- gsub(" ", "", gcluster, fixed = T)
tableBPG9 <- sort(table(gcluster), decreasing = T)
tableBPG9 <- tableBPG9/length(npaths)
jpeg("Bio_Process_G9_genesfreq.jpg",
     width = phi*600,
     height = 600,
     units = "px"
     
)
barplot(tableBPG9, col = 9, main = "Frecuencia de genes en cluster", ylab = "Frecuencia", 
        xlab = "Genes")
dev.off()

# Guardar tablas con 10 genes más repetidos de cada cluster

write.table(x = head(tableBPG1, 10), file = "Bio_Process_G1_top10genes.txt",
            sep = "\t", quote = F, col.names = F )

write.table(x = head(tableBPG2, 10), file = "Bio_Process_G2_top10genes.txt",
            sep = "\t", quote = F, col.names = F )

write.table(x = head(tableBPG3, 10), file = "Bio_Process_G3_top10genes.txt",
            sep = "\t", quote = F, col.names = F )

write.table(x = head(tableBPG4, 10), file = "Bio_Process_G4_top10genes.txt",
            sep = "\t", quote = F, col.names = F )

write.table(x = head(tableBPG5, 10), file = "Bio_Process_G5_top10genes.txt",
            sep = "\t", quote = F, col.names = F )

write.table(x = head(tableBPG6, 10), file = "Bio_Process_G6_top10genes.txt",
            sep = "\t", quote = F, col.names = F )

write.table(x = head(tableBPG7, 10), file = "Bio_Process_G7_top10genes.txt",
            sep = "\t", quote = F, col.names = F )

write.table(x = head(tableBPG8, 10), file = "Bio_Process_G8_top10genes.txt",
            sep = "\t", quote = F, col.names = F )

write.table(x = head(tableBPG9, 10), file = "Bio_Process_G9_top10genes.txt",
            sep = "\t", quote = F, col.names = F )
