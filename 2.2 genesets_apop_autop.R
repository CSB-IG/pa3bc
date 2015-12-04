################################################################################
# Get genes from selected GO ids and construct geneset list in GMT file format
# Author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Read adjusted GODB
GODB <- read.delim("GO_biologicalprocess_adjusted.txt", header = F, 
                   as.is = TRUE, row.names = 1)

# Read selected apoptosis and autophagy GO ids
apop_autop <- read.delim(file = "GO.apoptotic_process+autophagy.txt", as.is = T,
                         header = FALSE)

# Extract genesets from Bader et al., GODB
GODB_apop_autop <- GODB[apop_autop[,1],]

# Wipeout blank columns (non-found: no genes annotated), add NaNs in 2nd column
m <- na.omit(GODB_apop_autop)
m[,2] <- "NaN"
m <- as.matrix(m)

# Write genesets to TXT file 
write.table(x = m, file = "genesets_apop_autop.txt", row.names = TRUE, 
            col.names = F, sep = "\t")

# Write genesets to GMT file (for use in GSEA and PATHIFIER)
pathifier_labels <- paste(rownames(m), m[,1],sep = "_")
m[,1] <- pathifier_labels
write.table(x = m, file = "PATHIFIER_genesets_apop_autop.gmt", row.names = F, 
            col.names = F, sep = "\t", quote = FALSE)

################################################################################
# Only Apoptosis genesets
################################################################################
# Read selected apoptosis and autophagy GO ids
apop_autop <- read.delim(file = "GO.apoptotic_process_WITHOUT_ctsapop.txt", 
                         as.is = T, header = FALSE)

# Extract genesets from Bader et al., GODB
GODB_apop_autop <- GODB[apop_autop[,1],]

# Wipeout blanck columns (non-found: no genes annotated), add NaNs in 2nd column
m <- na.omit(GODB_apop_autop)
m[,2] <- "NaN"
m <- as.matrix(m)

# Write genesets to GMT file (for use in GSEA and PATHIFIER)
pathifier_labels <- paste(rownames(m), m[,1],sep = "_")
m[,1] <- pathifier_labels
write.table(x = m, file = "PATHIFIER_genesets_apop.gmt", row.names = F, 
            col.names = F, sep = "\t", quote = FALSE)

################################################################################
# Only autophagy genesets
################################################################################
# Read selected apoptosis and autophagy GO ids
apop_autop <- read.delim(file = "GO.autophagy.txt", as.is = T,
                         header = FALSE)

# Extract genesets from Bader et al., GODB
GODB_apop_autop <- GODB[apop_autop[,1],]

# Wipeout blanck columns (non-found: no genes annotated), add NaNs in 2nd column
m <- na.omit(GODB_apop_autop)
m[,2] <- "NaN"
m <- as.matrix(m)

# Write genesets to GMT file (for use in GSEA and PATHIFIER)
pathifier_labels <- paste(rownames(m), m[,1],sep = "_")
m[,1] <- pathifier_labels
write.table(x = m, file = "PATHIFIER_genesets_autop.gmt", row.names = F, 
            col.names = F, sep = "\t", quote = FALSE)
