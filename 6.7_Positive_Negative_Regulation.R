################################################################################
# PDS Regulación positiva y negativa 
## Author: Angel Garcia-Campos https://github.com/AngelCampos
################################################################################

################################################################################
### Installing and/or loading required packages "gplots" "RColorBrewer"
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
# creates a own color palette passing from blue, green yellow to red
my_palette <- colorRampPalette(c("blue", "cyan", "chartreuse1", "yellow", 
                                 "red", "firebrick4"))(n = 1000)
################################################################################
### Regulación positiva y negativa
################################################################################
pn <- c(grep("positive", rownames(meanPDS)), grep("negative",rownames(meanPDS)))
pn <- sort(rownames(meanPDS)[pn])
pn <- pn[-c(13,14,17,18,19,20,21,24,25,26)]
# Color labels
clabs <- pn
clabs[grep("positive", clabs)] <- "firebrick1"
clabs[grep("negative", clabs)] <- "dodgerblue1"
### Heatmap
pn_PDS <- meanPDS[pn,]
# Hierar Clust
row.distance = dist(pn_PDS, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")
col.distance = dist(t(pn_PDS), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")
# Write image
png(filename = "Positive-Negative_PDS.png", width = 1000, height = 800)
heatmap.2(pn_PDS,
          main = "Regulación Positiva-Negativa",  # heat map title
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
# Write image (NO Hier Clust)
png(filename = "Positive-Negative_PDS_noHC.png", width = 1000, height = 800,
    pointsize = 15)
heatmap.2(pn_PDS,
          # main = "Regulación Positiva-Negativa",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,40),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = "none", # apply selected clustering method
          Colv = "none", # apply selected clustering method
          keysize= 0.2           # size of color key
)
dev.off()

################################################################################
#### Comparacion PDS regulación positiva vs. negativa

pn_PDSmatrix <- PDSmatrix[pn,-c(1:61)]
pn_PDSmed <- apply(pn_PDSmatrix, 1, median)
pn_PDSmean <- apply(pn_PDSmatrix, 1, mean)

write.table(x = pn_PDSmean, file = "PN_regulation_PDS_MEAN.txt", quote = F, 
            sep = "\t", col.names = F)
write.table(x = pn_PDSmed, file = "PN_regulation_PDS_MEDIAN.txt", quote = F, 
            sep = "\t", col.names = F)

# Binomial test on negative regulation being overrepresented
binom.test(12, 14, p = 0.5, alternative = "g", conf.level = 0.95) # Mean values
binom.test(11, 14, p = 0.5, alternative = "g", conf.level =0.95) # Median values

################################################################################
##### Checar empalme de regulacion positiva y negativa
pn_genes <- gnames[pn, ]
pn_genes <- gsub(" ", "", pn_genes)
int <- seq(1,28, 2)
J <- NULL
for (i in int){
  x <- length(na.omit(intersect(pn_genes[i,], pn_genes[i+1,]))) / 
    length(na.omit(union(pn_genes[i,], pn_genes[i+1,])))
  J <- c(J, x)
}
names(J) <- substr(pn[int], start = 21, stop = 200)
J
View(as.matrix(J))
save.image() # Save session image
