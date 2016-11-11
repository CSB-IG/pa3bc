################################################################################
# Filtering of square matrices by Apoptosis and autophagy genes
# Author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Load bioprocess genes
genes <- as.character(as.matrix(read.delim(file = "genes_all.txt", header = F)))

# Filtering network FOR loop
for (i in 0:4){
  # Load square network
  net <- as.matrix(read.delim(file = paste("SQnetwork_", i , ".txt", sep= ""),
                            header = T, sep = "\t", row.names = 1))
  colnames(net) <- rownames(net)
  # Select bioprocess genes present in network
  igenes <- intersect(genes, rownames(net))
  # Select interactions in network
  selnet <- net[igenes,igenes]
  # Write filtered network to TXT file
  write.table(selnet, file = paste("SQnetwork_", i,"_Apo-Aut", ".txt", sep = ""),
            sep = "\t", col.names= NA, row.names= T, quote = F )
  print(paste("Network", i, "done!" )) # Progress message
}
