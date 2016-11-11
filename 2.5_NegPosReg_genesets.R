################################################################################
# Los genes de regulacion negativa y positiva se empalman en los genesets?
# Author: Miguel Angel García Campos - Github: angel.campos
################################################################################

gs <- as.matrix(read.delim(file = "PATHIFIER_genesets_apop_autop.gmt", header = F, 
                 row.names = 1)[,-1])

pn_gs <- gs[c(grep("POSITIVE", rownames(gs)), grep("NEGATIVE", rownames(gs))),]
pn_gs <- pn_gs[sort(rownames(pn_gs)),]
pn_gs <- pn_gs[-c(5,10,15,16,17,24,27,39,40,42,44,58),]


(length(intersect(as.character(pn_gs[1,]), as.character(pn_gs[2,])))-2) / 
  (length(union(as.character(pn_gs[1,]), as.character(pn_gs[2,]))) -2)

(length(intersect(as.character(pn_gs[3,]), as.character(pn_gs[4,])))-2) / 
  (length(union(as.character(pn_gs[3,]), as.character(pn_gs[4,]))) -2)

(length(intersect(as.character(pn_gs[5,]), as.character(pn_gs[6,])))-2) / 
  (length(union(as.character(pn_gs[5,]), as.character(pn_gs[6,]))) -2)

(length(intersect(as.character(pn_gs[7,]), as.character(pn_gs[8,])))-2) / 
  (length(union(as.character(pn_gs[7,]), as.character(pn_gs[8,]))) -2)

(length(intersect(as.character(pn_gs[37,]), as.character(pn_gs[38,])))-2) / 
  (length(union(as.character(pn_gs[37,]), as.character(pn_gs[38,]))) -2)

(length(intersect(as.character(pn_gs[39,]), as.character(pn_gs[40,])))-2) / 
  (length(union(as.character(pn_gs[39,]), as.character(pn_gs[40,]))) -2)

