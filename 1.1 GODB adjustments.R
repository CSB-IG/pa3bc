################################################################################
# Ajustes a base de datos de GeneOntology descargada de Bader et al.,
# http://baderlab.org/GeneSets
# Citar: Enrichment map: a network-based method for gene-set enrichment
# visualization and interpretation.
# Autor: Miguel Angel Garcia Campos
################################################################################

# Cargar matrix de bases de datos, remover genesets que no son GeneOntology
con <- "http://download.baderlab.org/EM_Genesets/March_24_2015/Human/symbol/Human_GOBP_AllPathways_no_GO_iea_March_24_2015_symbol.gmt"
no_col <- max(count.fields(con, sep="\t"))-1
GODB <- read.delim(con, col.names = 1:no_col, header = F, fill = T)
GODB <- GODB[2962:15008,]

# Modificar las etiquetas para quitar espacios en blanco
labels <- as.character(GODB[,1])
x <- gsub(" ", "_", labels)
y <- gsub("%", "_", x)
GODB[,1] <- y

#Funcion para extraer las categorias GO en formato "GO:XXXXXXXX"
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Extraer las etiquetas y asignar valores numericos a genesets no GO
etiquetas <- substrRight(labels, 10)
rownames(GODB) <- etiquetas

# Escribir a archivo TXT
write.table(x = GODB, file = "GO_biologicalprocess_adjusted.txt", sep = "\t",
            row.names = TRUE, col.names = FALSE )