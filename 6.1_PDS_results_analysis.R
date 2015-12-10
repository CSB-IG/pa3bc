###############################################################################
## Analisis de resultados
### Author: Angel García-Campos https://github.com/AngelCampos
###############################################################################

# State graphical variables
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch = 16)

# load/install pakages
if (!require("VennDiagram")) {
  install.packages("VennDiagram", dependencies = TRUE)
  library(VennDiagram)
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

###############################################################################
# Loading and labeling results data
###############################################################################

# Load Pathifier results and turn into a matrix

load ("PDS.RData")
PDSmatrix <- mapply(FUN = c, PDS$scores)
PDSmatrix <- t(PDSmatrix)

# Generating and assigning labels for pathways used in the analysis

## Genesets/Pathways
apop_autop <- as.matrix(read.delim(file = "PATHIFIER_genesets_apop_autop.gmt",
                                  header = F, sep = "\t", as.is = T))
apop <- as.matrix(read.delim(file = "PATHIFIER_genesets_apop.gmt",
                                  header = F, sep = "\t", as.is = T))
autop <- as.matrix(read.delim(file = "PATHIFIER_genesets_autop.gmt",
                                  header = F, sep = "\t", as.is = T))
array <- as.vector(as.matrix(read.delim(file = "probesets.txt", header = F)))

## Read input pathways from .txt
pathways <- read.delim("pathways.txt", header = F)
rownames(pathways) <- pathways[,2] # Naming rownames

## Extract NUMBER.ID of pathways that Pathifier succesfully useds
pathwaysInPDS <- rownames(PDSmatrix)  

## Selecting the NAMES of pathways
paths <- as.vector(pathways[,1])      # Extracting pathway names in vector
names(paths) <- pathways[,2]          # Nombrar vector con pathways
labels <- paths[pathwaysInPDS]        # Selects names of pathways as labels
rownames(PDSmatrix) <- labels         # assign labels as row names in PDSmatrix

write.table(x = labels, file = "pathways_in_PDS.txt", sep = "\t", col.names = F)

################################################################################
# Numero de pathways con 3 genes o mas
################################################################################

# All biological processes
x <- 0
for (i in 1:nrow(apop_autop)) {
  a <- apop_autop[i,]
  a[a == ""] <- NA
  x <- x + sum(length(na.omit(a))-1 >= 3)
  print(x)
}
# Autophagy
x <- 0
for (i in 124:nrow(apop_autop)) {
  a <- apop_autop[i,]
  a[a == ""] <- NA
  x <- x + sum(length(na.omit(a))-1 >= 3)
}
# Apoptosis
x <- 0
for (i in 1:123) {
  a <- apop_autop[i,]
  a[a == ""] <- NA
  x <- x + sum(length(na.omit(a))-1 >= 3)
}

################################################################################
# Venn diagrams 
################################################################################

# Genes in Gene Ontology
genes1 <- head(unique(as.vector(rbind(apop[,3:ncol(apop)]))), -1)
genes2 <- head(unique(as.vector(rbind(autop[,3:ncol(autop)]))), -1)

venn.diagram(x = list(Apoptosis = genes1, Autofagia = genes2), 
             filename = "Venn_GeneOntology_genes.tiff", col = "transparent",
             fill = c("darkorange1", "yellow1"), fontfamily = "serif",  
             fontface = "bold", rotation.degree = 45, margin = 0.1, cex =1.5,
             cat.cex = 2, cat.dist = 0.08)

genes0 <- sort(union(genes1, genes2))
genes0 <- genes0[2:length(genes0)]
write(x = genes0, file = "genesGO_apop_autop.txt", sep = " \t")

# Genes in microarray
genes3 <- intersect(genes1, array)
genes4 <- intersect(genes2, array)

venn.diagram(x = list(Apoptosis = genes3, Autofagia = genes4), 
             filename = "Venn_Array_genes.tiff", col = "transparent",
             fill = c("darkorange1", "yellow1"), fontfamily = "serif",  
             fontface = "bold", rotation.degree = 45, margin = 0.1, cex =1.5,
             cat.cex = 2, cat.dist = 0.08)

# Genes in PATHIFIER analysis

allg <- NULL
for (i in 1:96) {
  pg <- unlist(PDS$genesinpathway[i])
  names(pg) <- NULL
  allg <- c(allg, pg)
  allg <- unique(gsub(" ", "", allg, fixed = TRUE))
}

genes5 <- intersect(genes1, allg)
genes6 <- intersect(genes2, allg)

venn.diagram(x = list(Apoptosis = genes5, Autofagia = genes6), 
             filename = "Venn_Pathifier_genes.tiff", col = "transparent",
             fill = c("darkorange1", "yellow1"), fontfamily = "serif",  
             fontface = "bold", rotation.degree = 45, margin = 0.1, cex =1.5,
             cat.cex = 2, cat.dist = 0.08)

################################################################################
# PDS distribution
################################################################################

# Todas las muestras
phi <- 1.61803398 
jpeg("Distribución PDS.jpg",
     width = phi*600*0.7,
     height = 600*0.7,
     units = "px"
)
hist(PDSmatrix, col=4, main= "Distribución PDS", ylab= "Frecuencia", 
     xlab= "PDS", breaks = 25, ylim = c(0, 9000))
dev.off()

# Sanos
jpeg("Distribución PDS-Sanos.jpg",
     width = phi*600*0.7,
     height = 600*0.7,
     units = "px"
) # Name of png file
hist(PDSmatrix[,1:61], col=2, main= "Distribución PDS-Sanos", 
     ylab= "Frecuencia", xlab= "PDS", breaks =25, ylim = c(0, 9000))
dev.off()

# Enfermos
jpeg("Distribución PDS-Enfermos.jpg",
     width = phi*600*0.7,
     height = 600*0.7,
     units = "px"
)
hist(PDSmatrix[,62:880], col=1, main= "Distribución PDS-Enfermos", 
     ylab= "Frecuencia", xlab= "PDS", breaks =25, ylim = c(0, 9000))
dev.off()

################################################################################
# Jaccard index between processes
################################################################################

# Grado de interseccion para cada subproceso
BPg <- NULL
xtalk <- NULL
for (i in 1:96) {
  pg <- unlist(PDS$genesinpathway[i])
  names(pg) <- NULL
  BPg <- c(BPg, pg)
  BPg <- unique(gsub(" ", "", BPg, fixed = TRUE))
  p <- intersect(BPg, genes1)
  q <- intersect(BPg, genes2)
  if (length(p)>=length(q)) {
    xt <- length(q)/length(p)
    xtalk <- c(xtalk, xt)
  }
  if (length(q)>length(p)) {
    xtalk <- length(p)/length(q)
    xt <- length(q)/length(p)
    xtalk <- c(xtalk, xt)
  }
}

xtalk <- matrix(data = xtalk, nrow = length(xtalk), ncol = 1 )
xtalk <- cbind(xtalk, seq(0, 1, length.out = 96))
rownames(xtalk) <- rownames(PDSmatrix)

# Gray scale color pallette for Jaccard index
gray_scale <- colorRampPalette(c("white", "black"))(n = 50)

png("xtalk.png", # Name of png file       
    width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # font size

heatmap.2(xtalk,
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,21),     # widens margins around plot
          col=gray_scale,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          keysize= 0.8           # size of color key
          ## Color labeling columns (Opt. RowSideColors for rows)
)    
dev.off()               # close the PNG device

################################################################################
# Lists of genes pertaining ONLY to one category. Apop, Autop, or BOTH processes
################################################################################

genesautop <- setdiff(genes2, genes1)
genesapop <- setdiff(genes1, genes2)
genesinter <- intersect(genes1, genes2)

write(x = genesautop, file = "genes_autop.txt", sep = " \t")
write(x = genesapop, file = "genes_apop.txt", sep = " \t")
write(x = genesinter, file = "genes_inter.txt", sep = " \t")