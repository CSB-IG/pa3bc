################################################################################
# Analisis de subprocesos en los subgrupos
## Author: Angel Garcia-Campos https://github.com/AngelCampos
################################################################################

################################################################################
# Loading and labeling Pathifier results
################################################################################
phi <- 1.61803398 # phi constant

load ("PDS.RData")
PDSmatrix <- mapply(FUN = c, PDS$scores)
PDSmatrix <- t(PDSmatrix)

## Generating and assigning labels for pathways used in the analysis
# Read input pathways from .txt
pathways <- read.delim("pathways.txt", header = F)
rownames(pathways) <- pathways[,2] # Naming rownames

# Extract NUMBER.ID of pathways that Pathifier succesfully useds
pathwaysInPDS <- rownames(PDSmatrix)  

# Selecting the NAMES of pathways
paths <- as.vector(pathways[,1])      # Extracting pathway names in vector
names(paths) <- pathways[,2]          # Nombrar vector con pathways
labels <- paths[pathwaysInPDS]        # Selects names of pathways as labels
rownames(PDSmatrix) <- labels         # assign labels as row names in PDSmatrix

# Loading subgroup labeling info
subg <- as.matrix(read.delim("Samples_by_subgroup.txt",header=F, row.names = 1))
subg[870,] <- 0 #Removal of tumoral sample that groups with healthy samples
# Loading cluster label info
clus <- as.matrix(read.delim("BioProcess_by_cluster.txt",header=F,row.names= 1))
# Loading sample names
snames <- as.character(as.matrix(read.delim("sample_names.txt", header = F)))
colnames(PDSmatrix) <- snames

# Reorder PDSmatrix

# Orden de bio procesos arriba->abajo c(7,9,1,3,5,2,8,4,6)
# Orden de subgrupos izquierda->derecha c(1,4,5,3,2)
rord <- c(rownames(PDSmatrix[clus == 7,]),
          unname(rownames(PDSmatrix)[which(clus == 9)]),
          rownames(PDSmatrix[clus == 1,]), rownames(PDSmatrix[clus == 3,]),
          rownames(PDSmatrix[clus == 5,]), rownames(PDSmatrix[clus == 2,]),
          rownames(PDSmatrix[clus == 8,]), rownames(PDSmatrix[clus == 4,]),
          rownames(PDSmatrix[clus == 6,]))
cord <- c(colnames(PDSmatrix[,subg == 1]), colnames(PDSmatrix[,subg == 4]), 
          colnames(PDSmatrix[,subg == 5]), colnames(PDSmatrix[,subg == 3]), 
          colnames(PDSmatrix[,subg == 2]))
PDSord <- PDSmatrix[rord, cord] # PDS matrix reordered 

## Mean and median PDS for each process in TUMORAL samples

meanPDS <- sort(apply(PDSmatrix[,-c(1:61)], 1, mean), decreasing = T)
medianPDS <- sort(apply(PDSmatrix[,-c(1:61)], 1, median), decreasing = T)


write.table(x = meanPDS, file = "meanPDS_tumoral.txt", sep = "\t", quote = F,
            col.names = F)
write.table(x = medianPDS, file = "medianPDS_tumoral.txt", sep = "\t", quote = F,
            col.names = F)


