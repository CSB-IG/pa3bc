################################################################################
# 8.5 First neighbors filtering
# Author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Loading apoptosis and autophagy genes
genes <- as.matrix(read.delim(file = "genesGO_apop_autop.txt", header = F, 
                              sep = "\t"))
genes <- sort(as.character(genes))

# FOR LOOP: Load networks in SIF format, extract interactions with at least one
# autophagy or apoptosis node, print process results, and write new SIF file
for (i in 0:4){
  print(paste("Red", i))
  net <- read.delim(file = paste("SQnetwork_", i, "_0.999_DPI_0.1.sif", sep=""),
                    header = T, sep = "\t")
  net <- as.matrix(net)
  print(paste("Original number of interactions", length(net[,1])))
  inds1 <- which(net[,1] %in% genes)
  netx <- net[inds1,]
  inds2 <- which(net[,3] %in% genes)
  nety <- net[inds2,]
  netw <- rbind(netx,nety)
  dupl <- duplicated.matrix(netw) # Finding duplicated interactions
  netw <- netw[!dupl, ] # Removing duplicated interactions
  print(paste("Interacciones finales", length(netw[,1]),
              "=", signif(length(netw[,1])/length(net[,1]),2)))
  write.table(x = netw, file = paste("Final_network_", i , "_0.999_DPI_0.1.sif",
                                     sep = ""),
              sep = "\t", quote = FALSE, row.names = F)
  print("SIF file saved!")
}