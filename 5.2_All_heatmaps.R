################################################################################
## Heatmaps for Pathifier results in R
### Author: Angel García-Campos https://github.com/AngelCampos
#### Base by wonderful: Sebastian Raschka https://github.com/rasbt
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

################################################################################
# Loading and labeling Pathifier results
################################################################################

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

################################################################################
## Hierarchical clustering for non-square matrices
################################################################################

row.distance = dist(PDSmatrix, method = "euclidean") 
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(PDSmatrix), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

################################################################################
# Batch clusters
################################################################################
batch.number <- read.table(file = "batch_number.txt", header = T, row.names = 1)

myBATCHhc <- as.integer(batch.number)
myBATCHhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
myBATCHhc <- mycolhc[as.vector(mycl)] 

################################################################################
# PAM50 subtype information
################################################################################

subtype <- read.table(file = "PAM50subt_numeric.txt", header = T, row.names = 1)
myPAM50 <- as.integer(subtype)
myPAM50hc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
myPAM50hc <- mycolhc[as.vector(mycl)] 

################################################################################
# Plotting the Heatmaps!! (where all colorful things happen...)
################################################################################
# creates a own color palette passing from blue, green yellow to red
my_palette <- colorRampPalette(c("blue", "cyan", "chartreuse1", "yellow", 
                                 "red", "firebrick4"))(n = 1000)

## Healthy vs. tumoral heatmap #################################################
png("sanos_enfermos.png", # Name of png file       
    width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # font size
heatmap.2(PDSmatrix,
          main = "Sanos_Enfermos_ 880BC",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,21),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 0.8,           # size of color key
          ColSideColors=  c(           # Grouping row-samples into two different
            rep("dodgerblue", 61),    # categories, Apoptosis pathways
            rep("firebrick1", 819)),    # Autophagy pathways
          ## Color labeling rows
          RowSideColors= c(           # Grouping col-samples into two different
            rep("darkorange1", 81),    # categories, Samples 1-61: blue
            rep("yellow1", 15))    # Samples 62-880
)    
# Legend for ColumnSide color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Sanos", "Enfermos"), # category labels
       col = c("dodgerblue", "firebrick1"),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()               # close the PNG device

## Healthy vs. tumoral heatmap NO-hierarchical clustering ######################
png("sanos_enfermos_noHC.png", # Name of png file       
    width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # font size
heatmap.2(PDSmatrix,
          main = "",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,21),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = "None", # apply selected clustering method
          Colv = "None", # apply selected clustering method
          key = FALSE,
          keysize= 0.8,           # size of color key
          ColSideColors=  c(           # Grouping row-samples into two different
            rep("dodgerblue", 61),    # categories, Apoptosis pathways
            rep("firebrick1", 819)),    # Autophagy pathways
          RowSideColors= c(           # Grouping col-samples into two different
            rep("darkorange1", 81),    # categories, Samples 1-61: blue
            rep("yellow1", 15))    # Samples 62-880
)    
# Legend for ColumnSide color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Sanos", "Enfermos"), # category labels
       col = c("dodgerblue", "firebrick1"),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()               # close the PNG device

## Batch number ################################################################

png("lotes.png", # Name of png file       
    width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # font size
heatmap.2(PDSmatrix,
          main = "Lotes",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          key = TRUE,
          margins= c(10,21),     # widens margins around plot
          col= my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 0.8,           # size of color key
          ColSideColors= myBATCHhc,
          RowSideColors= c(           # Grouping col-samples into two different
            rep("darkorange1", 81),    # categories, Samples 1-61: blue
            rep("yellow1", 15)),    # Samples 62-880
)    
## Legend for ColumnSide color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("1","2","3","4","5","6","7","8","9","10"), # category labels
       col = rainbow(length(unique(mycl)), start=0.1, end=0.9),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()     

## Batch number NO-hierarchical clustering #####################################

png("lotes_noHC.png", # Name of png file       
    width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # font size
heatmap.2(PDSmatrix,
          main = "",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          key = FALSE,
          keysize = 0.8,
          margins= c(10,21),     # widens margins around plot
          col= my_palette,        # use on color palette defined earlier 
          Rowv = "none", # apply selected clustering method
          Colv = "none", # apply selected clustering method
          ## Color labeling columns (Opt. RowSideColors for rows)
          ColSideColors= myBATCHhc,
          RowSideColors= c(           # Grouping col-samples into two different
            rep("darkorange1", 81),    # categories, Samples 1-61: blue
            rep("yellow1", 15))    # Samples 62-880
)    
## Legend for ColumnSide color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("1","2","3","4","5","6","7","8","9","10"), # category labels
       col = rainbow(length(unique(mycl)), start=0.1, end=0.9),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()               # close the PNG device

## PAM50 subtype ###############################################################

png("subtipos_PAM50.png", # Name of png file       
    width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # font size
heatmap.2(PDSmatrix,
          main = "Subtipos tumorales PAM50",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          key = TRUE,
          margins= c(10,21),     # widens margins around plot
          col= my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 0.8,           # size of color key
          ColSideColors= myPAM50hc,
          RowSideColors= c(           # Grouping col-samples into two different
            rep("darkorange1", 81),    # categories, Samples 1-61: blue
            rep("yellow1", 15))    # Samples 62-880
)    
## Legend for ColumnSide color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Normals", "Basal", "LumA", "LumB", "Her2"), # category labels
       col = rainbow(length(unique(mycl)), start=0.1, end=0.9),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()               # close the PNG device

## PAM50 subtype NO-hierarchical clustering ####################################

png("subtipos_PAM50.png", # Name of png file       
    width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
    height = 7.5 * 400,     # 6 x 400 = 2400 px
    units = "px",         # px (Pixels = default), in (inches), cm or mm
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # font size
heatmap.2(PDSmatrix,
          main = "Subtipos tumorales PAM50",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          key = TRUE,
          margins= c(10,21),     # widens margins around plot
          col= my_palette,        # use on color palette defined earlier 
          Rowv = "none", # apply selected clustering method
          Colv = "none", # apply selected clustering method
          keysize= 0.8,           # size of color key
          ColSideColors= myPAM50hc,
          RowSideColors= c(           # Grouping col-samples into two different
            rep("darkorange1", 81),    # categories, Samples 1-61: blue
            rep("yellow1", 15))    # Samples 62-880
)    
## Legend for ColumnSide color labeling 
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Normals", "Basal", "LumA", "LumB", "Her2"), # category labels
       col = rainbow(length(unique(mycl)), start=0.1, end=0.9),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()               # close the PNG device

# 5 Subgroups ##################################################################

# Subgroup labeling
mycl <- cutree(col.cluster, k= 6)
mycl <-sub(2,1L,mycl); mycl <- sub(3,2L,mycl); mycl <- sub(4,3L,mycl)
mycl <- sub(5,4L,mycl); mycl <- sub(6,5L,mycl); mycl <- as.integer(mycl)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
mycolhc <- mycolhc[as.vector(mycl)]

# Heatmap
png(paste("Apop+Autop_Subgrupos.png", sep = ""), # Name of png file       
      width = 7.5 * 500,      # Easier re-scaling X*500 = Y pixels
      height = 7.5 * 400,     # 6 x 400 = 2400 px
      units = "px",         # px (Pixels = default), in (inches), cm or mm
      res = 300,            # 300 pixels per inch
      pointsize = 6)        # font size
heatmap.2(PDSmatrix,
          main = "Subgrupos: Apoptosis + Autofagia",  # heat map title
          density.info= "none",  # turns off density plot inside color legend
          trace= "none",         # turns off trace lines inside the heat map
          margins= c(10,21),     # widens margins around plot
          col=my_palette,        # use on color palette defined earlier 
          Rowv = as.dendrogram(row.cluster), # apply selected clustering method
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize= 0.8,           # size of color key
          ColSideColors= mycolhc,
          RowSideColors= c(           # Grouping col-samples into two different
            rep("darkorange1", 81),    # categories, Samples 1-61: blue
            rep("yellow1", 15))    # Samples 62-880
)
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c(as.character(1:5)), # category labels
       col = unique(mycolhc),  # color key
       lty= 1,          # line style
       lwd = 5, unit    # line width
)
dev.off()               # close the PNG device