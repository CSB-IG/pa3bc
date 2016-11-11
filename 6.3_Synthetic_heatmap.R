################################################################################
# Heatmap Sintético Grupos tumorales y grupos de procesos
### Author: Angel García-Campos https://github.com/AngelCampos
################################################################################

################################################################################
## Load Pathifier results and turn into a matrix
################################################################################

load ("PDS.RData")
PDSmatrix <- mapply(FUN = c, PDS$scores)
PDSmatrix <- t(PDSmatrix)

################################################################################
## Generating and assigning labels for pathways used in the analysis
################################################################################

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

###############################################################################
### Installing and/or loading required packages "gplots" "RColorBrewer"
###############################################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

################################################################################
## Creating Palette
################################################################################

# creates a own color palette passing from blue, green yellow to red
my_palette <- colorRampPalette(c("blue", "cyan", "chartreuse1", "yellow", 
                                 "red", "firebrick4"))(n = 1000)

################################################################################
## Clustering Methods
################################################################################

row.distance = dist(PDSmatrix, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(PDSmatrix), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

###############################################################################
# Subgroup labeling
###############################################################################

mycl <- cutree(col.cluster, k= 6)
mycl[mycl==2] <- 1; mycl[mycl==3] <- 2; mycl[mycl==4] <- 3; mycl[mycl==5] <- 4
mycl[mycl==6] <- 5
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
mycolhc <- mycolhc[as.vector(mycl)]
#Escribir resultados en archivo
write.table(x = as.matrix(mycl), file = "Samples_by_subgroup.txt", sep = "\t",
            quote = F, col.names = F)

myrw <- cutree(row.cluster, k = 9)
myrowhc <- rainbow(length(unique(myrw)), start=0.1, end=0.9)
myrowhc <- myrowhc[as.vector(myrw)]
# Write res to file
write.table(x = as.matrix(myrw), file = "BioProcess_by_cluster.txt", sep = "\t",
            quote = F, col.names = F)

################################################################################
## Dividing the processes in groups 
################################################################################

png(paste("Synthetic_Heatmap_9.png", sep = ""), # Name of png file       
    width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # font size

heatmap.2(PDSmatrix,
          main = "Subgrupos tumorales / Grupos de procesos",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,21),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 0.8,           # size of color key
          ## Color labeling columns (Opt. RowSideColors for rows)
          ColSideColors= mycolhc,
          RowSideColors= myrowhc    # Apotosis = orange; autophagy =  yellow
)

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c(as.character(1:5)), # category labels
       col = unique(mycolhc),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)

dev.off()               # close the PNG device

################################################################################
# Separación y promedio de datos de subgrupos
################################################################################
# Subgrupo 1
x1 <- NULL
for(i in 1:9){
  m <- mean(PDSmatrix[myrw == i, mycl == 1])
  x1 <- rbind(x1,m)
  print(x1)
}
# Subgrupo 2
x2 <- NULL
for(i in 1:9){
  m <- mean(PDSmatrix[myrw == i, mycl == 2])
  x2 <- rbind(x2,m)
  print(x2)
}
# Subgrupo 3
x3 <- NULL
for(i in 1:9){
  m <- mean(PDSmatrix[myrw == i, mycl == 3])
  x3 <- rbind(x3,m)
  print(x3)
}
# Subgrupo 4
x4 <- NULL
for(i in 1:9){
  m <- mean(PDSmatrix[myrw == i, mycl == 4])
  x4 <- rbind(x4,m)
  print(x4)
}
# Subgrupo 5
x5 <- NULL
for(i in 1:9){
  m <- mean(PDSmatrix[myrw == i, mycl == 5])
  x5 <- rbind(x5,m)
  print(x5)
}

#Unir por columnas
sumPDS <- cbind(x1,x2,x3,x4,x5)
sumPDSPRIME <- sumPDS[c(7,9,1,3,5,2,8,4,6), c(1,4,5,3,2)]
rownames(sumPDSPRIME) <- c("BPG1","BPG2","BPG3","BPG4","BPG5","BPG6","BPG7",
                           "BPG8","BPG9")
colnames(sumPDSPRIME) <- c("Healthy","basal+her2","lumB","lumA","lumB+her2")

sumPDSsignif <- signif(sumPDSPRIME, 2) # Round results to 2 decimals

################################################################################
## Plotting the SYNTHETIC Heatmap!! (where all colorful things happen... again)
################################################################################
# creates a own color palette passing from blue, green yellow to red
new_palette <- colorRampPalette(c("blue", "cyan", "chartreuse1", "yellow", 
                                 "red", "firebrick4"))(n = 100)

png(paste("Synthetic_Heatmap.png", sep = ""), # Name of png file       
    width = 7.5 * 400,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 17)        # font size

heatmap.2(sumPDSPRIME,
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          key = FALSE,
          margins= c(9,5),     # widens margins around plot
          col=new_palette,        # use on color palette defined earlier 
          Rowv = "none", # apply selected clustering method
          Colv = "none", # apply selected clustering method
          # lmat = c(10,10),
          lhei = c(10,100),
          lwid = c(10,100)
)

dev.off()               # close the PNG device

### SQUARED cells Synthetic Heatmap

png(paste("SQ_Synthetic_Heatmap.png", sep = ""), # Name of png file       
    width = 5.5 * 400,      # Easier re-scaling X*500 = Y pixels
    height = 9.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 17)        # font size

heatmap.2(sumPDSPRIME,
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          key = FALSE,
          margins= c(9,5),     # widens margins around plot
          col=new_palette,        # use on color palette defined earlier 
          Rowv = "none", # apply selected clustering method
          Colv = "none", # apply selected clustering method
          # lmat = c(10,10),
          lhei = c(10,100),
          lwid = c(10,100)
)

dev.off()               # close the PNG device

### Version with quantitative data notes

png(paste("SQ_Synthetic_Heatmap_cuantitative.png", sep = ""), # Name of png file
    width = 5.5 * 400,      # Easier re-scaling X*500 = Y pixels
    height = 9.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 17)        # font size

heatmap.2(sumPDSPRIME,
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          key = FALSE,
          margins= c(9,5),     # widens margins around plot
          col=new_palette,        # use on color palette defined earlier 
          Rowv = "none", # apply selected clustering method
          Colv = "none", # apply selected clustering method
          # lmat = c(10,10),
          lhei = c(10,100),
          lwid = c(10,100),
          cellnote = sumPDSsignif,
          notecol = "black"
)

dev.off()               # close the PNG device

################################################################################
# Informacion de entrecruzamiento por cada subproceso
################################################################################

# Grado de interseccion para cada subproceso
nBP <- c(1:96)
xsyn <- NULL
for (i in 1:9) {
  BPGg <- NULL
  BPGg <- unlist(PDS$genesinpathway[nBP[myrw == i]])
  names(BPGg) <- NULL
  BPGg <- unique(gsub(" ", "", BPGg, fixed = TRUE))
  p <- intersect(BPGg, genes1)
  q <- intersect(BPGg, genes2)
  if (length(p)>length(q)) {
    xt <- length(q)/length(p)
    xsyn <- c(xsyn, xt)
  }
  if (length(q)>length(p)) {
    xsyn <- length(p)/length(q)
    xt <- length(q)/length(p)
    xsyn <- c(xsyn, xt)
  }
}
xsyn <- xsyn[c(7,9,1,3,5,2,8,4,6)] #Re order
xsyn <- matrix(data = xsyn, 9, 1)
xsyn <- cbind(xsyn, seq(0, 1, length.out = 9))
rownames(xsyn) <- rownames(sumPDSPRIME)
colnames(xsyn) <- c("X-t", "g-s")
round_xsyn <- round(xsyn, 2)

# creates a own color palette passing from blue, green yellow to red
gray_scale <- colorRampPalette(c("white", "black"))(n = 50)

png(paste("SQ_Synthetic_Xtalk.png", sep = ""), # Name of png file       
    width = 3 * 400,      # Easier re-scaling X*500 = Y pixels
    height = 9 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 20)        # font size

heatmap.2(xsyn,
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          key = FALSE,
          # margins= c(5,20),     # widens margins around plot
          col= gray_scale,        # use on color palette defined earlier 
          Rowv = "none", # apply selected clustering method
          Colv = "none", # apply selected clustering method
          # lmat = c(10,10),
          lhei = c(10,100),
          lwid = c(10,100),
          cellnote = round_xsyn,
          notecol = "black"
          )

dev.off()


################################################################################
# Procesos biologicos y genes de cada grupo
################################################################################
# Reordering subgroups bioprocess number

# Escribir un archivo con los bio procesos de cada grupo
#1
bpn <- rownames(PDSdataframe[myrw == 7,])
write(x = bpn, file = paste("Bio_Process_G", 1, ".txt", sep = ""),sep = "\t")
#2
bpn <- rownames(PDSdataframe[myrw == 9,])
write(x = bpn, file = paste("Bio_Process_G", 2, ".txt", sep = ""),sep = "\t")
#3
bpn <- rownames(PDSdataframe[myrw == 1,])
write(x = bpn, file = paste("Bio_Process_G", 3, ".txt", sep = ""),sep = "\t")
#4
bpn <- rownames(PDSdataframe[myrw == 3,])
write(x = bpn, file = paste("Bio_Process_G", 4, ".txt", sep = ""),sep = "\t")
#5
bpn <- rownames(PDSdataframe[myrw == 5,])
write(x = bpn, file = paste("Bio_Process_G", 5, ".txt", sep = ""),sep = "\t")
#6
bpn <- rownames(PDSdataframe[myrw == 2,])
write(x = bpn, file = paste("Bio_Process_G", 6, ".txt", sep = ""),sep = "\t")
#7
bpn <- rownames(PDSdataframe[myrw == 8,])
write(x = bpn, file = paste("Bio_Process_G", 7, ".txt", sep = ""),sep = "\t")
#8
bpn <- rownames(PDSdataframe[myrw == 4,])
write(x = bpn, file = paste("Bio_Process_G", 8, ".txt", sep = ""),sep = "\t")
#9
bpn <- rownames(PDSdataframe[myrw == 6,])
write(x = bpn, file = paste("Bio_Process_G", 9, ".txt", sep = ""),sep = "\t")

c(7,9,1,3,5,2,8,4,6)
# Obtener la lista de genes de cada grupo
#1
BPGg <- unlist(PDS$genesinpathway[nBP[myrw == 7]])
write(x = BPGg, file = paste("Bio_Process_G", 1 , "_genes.txt", sep = ""))
#2
BPGg <- unlist(PDS$genesinpathway[nBP[myrw == 9]])
write(x = BPGg, file = paste("Bio_Process_G", 2 , "_genes.txt", sep = ""))
#3
BPGg <- unlist(PDS$genesinpathway[nBP[myrw == 1]])
write(x = BPGg, file = paste("Bio_Process_G", 3 , "_genes.txt", sep = ""))
#4
BPGg <- unlist(PDS$genesinpathway[nBP[myrw == 3]])
write(x = BPGg, file = paste("Bio_Process_G", 4 , "_genes.txt", sep = ""))
#5
BPGg <- unlist(PDS$genesinpathway[nBP[myrw == 5]])
write(x = BPGg, file = paste("Bio_Process_G", 5 , "_genes.txt", sep = ""))
#6
BPGg <- unlist(PDS$genesinpathway[nBP[myrw == 2]])
write(x = BPGg, file = paste("Bio_Process_G", 6 , "_genes.txt", sep = ""))
#7
BPGg <- unlist(PDS$genesinpathway[nBP[myrw == 8]])
write(x = BPGg, file = paste("Bio_Process_G", 7 , "_genes.txt", sep = ""))
#8
BPGg <- unlist(PDS$genesinpathway[nBP[myrw == 4]])
write(x = BPGg, file = paste("Bio_Process_G", 8 , "_genes.txt", sep = ""))
#9
BPGg <- unlist(PDS$genesinpathway[nBP[myrw == 6]])
write(x = BPGg, file = paste("Bio_Process_G", 9 , "_genes.txt", sep = ""))


# List with all genes used in PATHIFIER analysis
genes0 <- union(genes5, genes6)
write(x = genes0, file = paste("Bio_Process_G", 0 , "_genes.txt", sep = ""))