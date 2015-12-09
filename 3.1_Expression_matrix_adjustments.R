################################################################################
# Cutoff expression matrix to annotated gene symbols and give format for use
# in ARACNE, PATHIFIER and PAM50 analyses
# Author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Read gene expression data
exp.matrix <- read.delim(file ="EXPMATRIX_SANOS_ENFERMOS_COLAPSED.txt", 
                         as.is = T, row.names = 1)
exp.matrix <- as.matrix(exp.matrix)

# Cutoff expression matrix to annotated GeneSymbols

all_genes <- read.delim(file= "probesets.txt", header = F, sep = "\t")
all_genes <- as.matrix(all_genes)
ALLG_emat <- exp.matrix[all_genes,]

# Write expression matrix to TXT file
write.table(x = ALLG_emat, file = "ARACNE_matriz_de_expresion_reducida.txt",
            sep = "\t", row.names = TRUE, col.names = NA)

## Cutoff expression matrix to genes in GO genesets#############################
GO_genes <- read.delim(file= "genesets_apop_autop.txt", header = F, sep = "\t", 
           row.names = 1)
genes <- as.matrix(GO_genes[,3:ncol(GO_genes)])
genes <- unique(as.character(genes))
GOgenes_exp_mat <- as.data.frame(exp.matrix)[genes,]
GOgenes_exp_mat <- na.omit(GOgenes_exp_mat)

# Give GOgenes_exp_mat PATHIFIER gene expression data format
# 61 healthy individual samples and 819 tumoral samples
status <- matrix(data = c(rep(1,61), rep(0, 819)), nrow = 1)
status <- as.character(status)
GOg_exp_mat <- rbind(status, GOgenes_exp_mat)
rownames(GOg_exp_mat)[1] <- c("NORMALS")
samples <- colnames(exp.matrix)
GOg_exp_mat <- rbind(samples,GOg_exp_mat)
rownames(GOg_exp_mat)[1] <- c("NAME")

# Write expression matrix to TXT file - for use in PATHIFIER
write.table(x= GOg_exp_mat, file = "PATHIFIER_matriz_de_expresion_reducida.txt",
            sep = "\t", row.names = TRUE, col.names = FALSE)

## Cutoff expression matrix to genes in PAM50###################################
pam50_genes <- as.matrix(read.delim("pam50_genes.txt", header = F))
pam50_genes <- unique(as.character(pam50_genes))
pam50_expm <- as.data.frame(exp.matrix)[pam50_genes,]
pam50_expm <- na.omit(pam50_expm)

#Write to TXT file
write.table(x = pam50_expm, file = "PAM50_matriz_de_expresion_reducida.txt", 
            sep = "\t", col.names = NA, quote = FALSE)