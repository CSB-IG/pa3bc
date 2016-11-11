################################################################################
# Analisis de subprocesos en los subgrupos
## Author: Angel Garcia-Campos https://github.com/AngelCampos
################################################################################

################################################################################
# Installing and/or loading required packages "gplots" "RColorBrewer"
################################################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("Biobase")) {
  install.packages("Biobase", dependencies = TRUE)
  library(Biobase)
}

################################################################################
## Creating Palette
################################################################################

# creates a own color palette passing from blue, green yellow to red
my_palette <- colorRampPalette(c("blue", "cyan", "chartreuse1", "yellow", 
                                 "red", "firebrick4"))(n = 1000)

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

# Sub matrices of respective subgroups (from heatmap)
PDS_0 <- PDSmatrix[,subg==1]
PDS_1 <- PDSmatrix[,subg==4]
PDS_2 <- PDSmatrix[,subg==5]
PDS_3 <- PDSmatrix[,subg==3]
PDS_4 <- PDSmatrix[,subg==2]

rM_0 <- rowMeans(PDS_0); rM_1 <- rowMeans(PDS_1)
rM_2 <- rowMeans(PDS_2); rM_3 <- rowMeans(PDS_3)
rM_4 <- rowMeans(PDS_4)

# Bind row means of subgroups
meanPDS <- as.matrix(data.frame(SG0 = rM_0, SG1 = rM_1, SG2 = rM_2, SG3 = rM_3,
                                SG4 =rM_4))
rownames(meanPDS) <- rownames(PDSmatrix)

# Median PDS
rmed_0 <- rowMedians(PDS_0); rmed_1 <- rowMedians(PDS_1)
rmed_2 <- rowMedians(PDS_2); rmed_3 <- rowMedians(PDS_3)
rmed_4 <- rowMedians(PDS_4)
medianPDS <- as.matrix(data.frame(SG0 = rmed_0, SG1 = rmed_1, SG2 = rmed_2, 
                                  SG3 = rmed_3, SG4 =rmed_4))
rownames(medianPDS) <- rownames(PDSmatrix)

################################################################################
## Clustering Methods
################################################################################

row.distance = dist(meanPDS, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(meanPDS), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")
myrw <- cutree(row.cluster, k = 10)
myrowhc <- rainbow(length(unique(myrw)), start=0.1, end=0.9)
myrowhc <- myrowhc[as.vector(myrw)]

################################################################################
# Heatmap Mean PDS
################################################################################
png(paste("MeanPDS_Heatmap.png", sep = ""), # Name of png file       
    width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # font size
heatmap.2(meanPDS,
          main = "PDS Promedio para cada subgrupo",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,40),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 1,           # size of color key
          RowSideColors= c(           # Grouping col-samples into two different
            rep("slateblue", 81),    # categories, Samples 1-61: blue
            rep("magenta", 15))    # Samples 62-880
)
# Legend for Color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Apoptosis", "Autofagia"), # category labels
       col = c("slateblue","magenta"),  # color key
       lty= 1,          # line style
       lwd = 8, cex = 3     # line width
)
dev.off()               # close the PNG device

# Plot PDS boxplot
png(filename = "PDS_by_group.png", width = 600*phi, height = 600, pointsize = 18)
boxplot(meanPDS, col =brewer.pal(5, "Set1"), ylim = c(0,1), 
        main = "Mean PDS by subgroup" )
dev.off()

# PDS por subgrupo
SG0_PDS <- sort(meanPDS[,"SG0"], decreasing = T)
SG1_PDS <- sort(meanPDS[,"SG1"], decreasing = T)
SG2_PDS <- sort(meanPDS[,"SG2"], decreasing = T)
SG3_PDS <- sort(meanPDS[,"SG3"], decreasing = T)
SG4_PDS <- sort(meanPDS[,"SG4"], decreasing = T)

write.table(SG0_PDS, file = "SG0_PDS.txt",quote = F, sep = "\t", col.names = F )
write.table(SG1_PDS, file = "SG1_PDS.txt",quote = F, sep = "\t", col.names = F )
write.table(SG2_PDS, file = "SG2_PDS.txt",quote = F, sep = "\t", col.names = F )
write.table(SG3_PDS, file = "SG3_PDS.txt",quote = F, sep = "\t", col.names = F )
write.table(SG4_PDS, file = "SG4_PDS.txt",quote = F, sep = "\t", col.names = F )

# Order matrix in hierarchical clusterin order

meanPDS <- meanPDS[rev(row.cluster$order),]
write.table(meanPDS, file = "ALL_PDS.txt",quote = F, sep = "\t", col.names = F )

# Numero de genes por pathway analizado
gnumber <- NULL
for (i in 1:96){
  x <- length(unlist(PDS$genesinpathway[i]))
  gnumber <- c(gnumber, x)
}

# Genes por pathway analizado
gnames <- NULL
for (i in 1:96){
  x <- unlist(PDS$genesinpathway[i])
  length(x) <- max(gnumber)
  gnames <- rbind(gnames, x)
}
colnames(gnames) <- NULL
rownames(gnames) <- rownames(PDSmatrix)

write.table(x = gnames, file = "Genes_in_PDS.txt", sep = "\t",
            quote = F, col.names = F, na = "")

## Mediana de los PDS promedio
medPDS <- NULL
for (i in 1:nrow(meanPDS)){
  x <- median(meanPDS[i, 2:5])
  medPDS <- c(medPDS, x)
}
names(medPDS) <- rownames(meanPDS); medPDS <- sort(medPDS, decreasing = T)
med_PDS <- meanPDS[names(medPDS),]
write.table(med_PDS, file = "Median_PDS_bySG.txt", quote = F, sep = "\t", 
            col.names = F)

## Desviacion estándar de los procesos biológicos
sdPDS <- NULL
for (i in 1:nrow(meanPDS)){
  x <- sd(meanPDS[i, 2:5])
  sdPDS <- c(sdPDS, x)
}
names(sdPDS) <- rownames(meanPDS); sdPDS <- sort(sdPDS, decreasing = T)
sd_PDS <- meanPDS[names(sdPDS),]
write.table(sd_PDS, file = "SDordered_PDS_bySG.txt", quote = F, sep = "\t", 
            col.names = F)

################################################################################
## Heatmap para desviación estandar (MÁS cambios)

row.distance = dist(sd_PDS[1:10,], method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(sd_PDS[1:10,]), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")
#BioProcess labels
bplabs <- c("slateblue","magenta","slateblue","slateblue","slateblue",
            "magenta","slateblue","slateblue","slateblue","slateblue")
png(filename = "SD_Heatmap_more.png", width = 1000, height = 800, pointsize = 16)
heatmap.2(sd_PDS[1:10,], 
          main = "10 BioProcesses - Most SD",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,40),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 1,           # size of color key
          RowSideColors= bplabs
)
# Legend for Color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Apoptosis", "Autofagia"), # category labels
       col = c("slateblue","magenta"),  # color key
       lty= 1,          # line style
       lwd = 8, unit    # line width
)
dev.off()

## Heatmap para desviación estandar (MENOS cambios)

row.distance = dist(sd_PDS[87:96,], method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(sd_PDS[87:96,]), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

#BioProcess labels
bplabs <- c("slateblue","slateblue","slateblue","slateblue","slateblue",
            "slateblue","slateblue","slateblue","slateblue","slateblue")

#Heatmap in PNG
png(filename = "SD_Heatmap_less.png", width = 1000, height = 800, pointsize = 16)
heatmap.2(sd_PDS[87:96,],
          main = "10 BioProcesses - Least SD",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,40),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 1,           # size of color key
          RowSideColors= bplabs
)
# Legend for Color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Apoptosis", "Autofagia"), # category labels
       col = c("slateblue","magenta"),  # color key
       lty= 1,          # line style
       lwd = 8, unit    # line width
)
dev.off()

## Heatmap para Valores por mediana (MAS desregulación)

row.distance = dist(med_PDS[1:10,], method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(med_PDS[1:10,]), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

#BioProcess labels
bplabs <- c("slateblue","magenta","slateblue","slateblue","magenta",
            "magenta","magenta","slateblue","slateblue","slateblue")

png(filename = "MED_Heatmap_more.png", width = 1000, height = 800, pointsize = 16)
heatmap.2(med_PDS[1:10,],
          main = "10 most deregulated BioProc",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,40),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 1,           # size of color key
          RowSideColors= bplabs
)
# Legend for Color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Apoptosis", "Autofagia"), # category labels
       col = c("slateblue","magenta"),  # color key
       lty= 1,          # line style
       lwd = 8, unit    # line width
)
dev.off()

## Heatmap para Valores por mediana (MENOS desregulación)

row.distance = dist(med_PDS[87:96,], method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(med_PDS[87:96,]), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")
#BioProcess labels
bplabs <- c("slateblue","slateblue","slateblue","slateblue","slateblue",
            "magenta","slateblue","slateblue","slateblue","slateblue")

png(filename = "MED_Heatmap_less.png", width = 1000, height = 800, pointsize = 16)
heatmap.2(med_PDS[87:96,],
          main = "10 less deregulated BioProc",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,40),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 1,           # size of color key
          RowSideColors= bplabs
)
# Legend for Color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Apoptosis", "Autofagia"), # category labels
       col = c("slateblue","magenta"),  # color key
       lty= 1,          # line style
       lwd = 8, unit    # line width
)
dev.off()

### Regulación positiva y negativa

pn <- c(grep("positive", rownames(meanPDS)), grep("negative", rownames(meanPDS)))
pn <- sort(rownames(meanPDS)[pn])
pn <- pn[-c(13,14,17,18,19,20,21,24,25,26)]
# Color labels
clabs <- pn
clabs[grep("positive", clabs)] <- "firebrick1"
clabs[grep("negative", clabs)] <- "dodgerblue1"

rsc <- c(grep("positive", pn), -grep("negative", pn))
rsc <- 
  # Heatmap
  pn_PDS <- meanPDS[pn,]

# Hierar Clust
row.distance = dist(pn_PDS, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")
col.distance = dist(t(pn_PDS), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

png(filename = "Positive-Negative_PDS.png", width = 1000, height = 800, pointsize = 16)
heatmap.2(pn_PDS,
          main = "Regulación Positiva-Negativa",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,40),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 1,           # size of color key
          RowSideColors = clabs
)
# Legend for Color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Regulación Positiva", "Regulación Negativa"), # category labels
       col = c("firebrick1","dodgerblue1"),  # color key
       lty= 1,          # line style
       lwd = 8, unit    # line width
)
dev.off()

#### Comparacion PDS regulación positiva vs. negativa

pn_PDSmatrix <- PDSmatrix[pn,-c(1:61)]
pn_PDSmed <- apply(pn_PDSmatrix, 1, median)
pn_PDSmean <- apply(pn_PDSmatrix, 1, mean)

write.table(x = pn_PDSmean, file = "PN_regulation_PDS_MEAN.txt", quote = F, sep = "\t",
            col.names = F)
write.table(x = pn_PDSmed, file = "PN_regulation_PDS_MEDIAN.txt", quote = F, sep = "\t",
            col.names = F)

# Binomial test on negative regulation being overrepresented
binom.test(12, 14, p = 0.5, alternative = "g", conf.level = 0.95) # Mean values
binom.test(11, 14, p = 0.5, alternative = "g", conf.level =0.95) # Median values


