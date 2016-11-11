################################################################################
# Intrinsic and extrinsic apoptosis
## Author: Angel Garcia-Campos https://github.com/AngelCampos
################################################################################

### Apoptosis Extrinseca e intrínseca

ei <- c(grep("intrin", rownames(meanPDS)), grep("extrin", rownames(meanPDS)))
ei <- sort(rownames(meanPDS)[ei]); names(ei) <- NULL
ei <- ei[c(7,9,26,32,27,33,28,34)]


# Color labels
clabs <- ei
clabs[grep("intrin", clabs)] <- "firebrick1"
clabs[grep("extrin", clabs)] <- "dodgerblue1"

### Heatmap
ei_PDS <- meanPDS[ei,]
# Hierar Clust
row.distance = dist(ei_PDS, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")
col.distance = dist(t(ei_PDS), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

png(filename = "Intrinsic_Extrinsic_PDS.png", width = 1000, height = 800)
heatmap.2(ei_PDS,
          main = "Apoptosis Intríseca-Extrínseca",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,40),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 0.8,           # size of color key
          RowSideColors = clabs
)
dev.off()

