################################################################################
# Script Genes de cada proceso biologico Apoptosis// Atofagia//Ambos
# Author: Miguel Angel García Campos - Github: angel.campos
################################################################################

# Read data
GSapo <- as.matrix(read.delim(file = "PATHIFIER_genesets_apop.gmt", header = F))
rownames(GSapo) <- GSapo[,1]; GSapo <- GSapo[,-2:-1]
GSaut <- as.matrix(read.delim(file = "PATHIFIER_genesets_autop.gmt",header = F))
rownames(GSaut) <- GSaut[,1]; GSaut <- GSaut[,-2:-1]

# All genes
gapo <- na.omit(unique(as.character(GSapo))); gapo <- gapo[gapo != ""]
gaut <- na.omit(unique(as.character(GSaut))); gaut <- gaut[gaut != ""]
gint <- intersect(gapo,gaut)

# Write gene lists
write(x = gapo, file = "genes_apoptosis.txt",sep = "\t")
write(x = gaut, file = "genes_autofagia.txt",sep = "\t")
write(x = gint, file = "genes_interseccion.txt",sep = "\t")
