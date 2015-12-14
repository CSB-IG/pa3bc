################################################################################
# Converting from square matrix to SIF file (for use in Cytoscape)
# Author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Install/load packages
source("http://bioconductor.org/biocLite.R")
if (!require("igraph")) {
  biocLite("igraph", ask = FALSE)
  library(igraph)
}

#Loop ADJ -> GRAPH -> SIF (from network 0 to 4)
for (i in 0:4) {
  nw <- as.matrix(read.table(paste("SQnetwork_", i, "_0.999_DPI_0.1.txt", 
                                   sep= ""), header=T, sep="\t", row.names = 1))
  nw <- nw[sort(rownames(nw)), sort(colnames(nw))]
  colnames(nw) <- rownames(nw)
  sifNET <- graph.adjacency(adjmatrix= nw, mode='undirected', diag=F,weighted=T)
  xsif <- get.data.frame(sifNET)
  xsif <- xsif[, c(1,3,2)]
  length(xsif[,1])
  write.table(xsif, file = paste("SQnetwork_",i,"_0.999_DPI_0.1.sif"), 
              col.names= T, row.names= F, quote = F, sep = "\t")
  print(paste("Network", i, "done."))
}