################################################################################
# Applying Data Processing Inequality from MINET pakage to filtered network
# Based on wonderful: Diana Drago https://github.com/dianadrago
# Author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Load/Install packages
if (!require("minet")) {
  biocLite("minet", ask =FALSE)
  library(minet)
}

# DPI loop for 5 networks from 0 to 4
for (i in 0:4) {
  nw <- as.matrix(read.table(paste("SQnetwork_", i, "_0.999.txt", sep= ""),
                             header=T, sep="\t", row.names = 1))
  nw <- nw[sort(rownames(nw)), sort(colnames(nw))]
  colnames(nw) <- rownames(nw)
  print(isSymmetric(nw))
  # Apply DPI
  DPInw <- aracne(mim = nw, eps = 0.1)
  # Print resulting number of edges
  edges <- sum(nw > 0)
  edgesPRIME <- sum(DPInw > 0)
  print(paste("Resultados red",i))
  print(paste("Numero de aristas en red original =", edges))
  print(paste ("Numero de aristas en red filtrada =", edgesPRIME))
  print(paste ("Porcentaje de aristas conservados", edgesPRIME*100/edges))
  # Write new ADJ
  write.table(DPInw, file = paste("SQnetwork_",i,"_0.999_DPI_0.1.txt", sep =""), 
              sep = "\t", col.names= NA, quote = F)
}