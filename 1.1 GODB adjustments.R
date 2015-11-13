################################################################################
# Adjustments to Gene Ontology DB downloaded from Bader et al.,
# http://baderlab.org/GeneSets
# Merico, Daniele, et al. "Enrichment map: a network-based method for gene-set 
# enrichment visualization and interpretation." PloS one 5.11 (2010)
# Script author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Load Gene Ontology DB, remove genesets with no GO identifier
con <- "http://download.baderlab.org/EM_Genesets/March_24_2015/Human/symbol/Human_GOBP_AllPathways_no_GO_iea_March_24_2015_symbol.gmt"
no_col <- max(count.fields(con, sep="\t"))-1
GODB <- read.delim(con, col.names = 1:no_col, header = F, fill = T)
GODB <- GODB[2962:15008,]

# Modify labels to remove white spaces
labels <- as.character(GODB[,1])
x <- gsub(" ", "_", labels)
y <- gsub("%", "_", x)
GODB[,1] <- y

# Function: Extract GO ids in format "GO:XXXXXXXX"
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Extract labels and assign to rownames
etiquetas <- substrRight(labels, 10)
rownames(GODB) <- etiquetas

# Write adjusted GODB to TXT file
write.table(x = GODB, file = "GO_biologicalprocess_adjusted.txt", sep = "\t",
            row.names = TRUE, col.names = FALSE )