################################################################################
# Filtering of square matrices by same quantile value (0.999)
# Author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Filtering network FOR loop
for (i in 0:4){
  # Loading square matrix file
  net <- as.matrix(read.delim(file = paste("SQnetwork_", i , ".txt", sep= ""),
                header = T, sep = "\t", row.names = 1))
  edges <- sum(net > 0) # Number of edges in original network
  pct <- 0.999 # Cutting percentile
  mmi <- quantile(net, pct) # Minimum Mutual Information to keep value
  # Filtering MI values lower than "mmi"
  netPRIME <- net
  netPRIME[netPRIME <= mmi] <- 0
  edgesPRIME <- sum(netPRIME > 0)
  # Print process results
  print(paste("Number of edges in original network =", edges))
  print(paste ("Number of edges in filtered matrix =", edgesPRIME))
  print(paste ("Percentage of conserved edges =", edgesPRIME*100/edges))
  # Write filtered network to TXT file
  write.table(netPRIME, file = paste("SQnetwork_", i, "_", pct, ".txt", sep = ""),
              sep = "\t", col.names= NA, row.names= T, quote = F )
}