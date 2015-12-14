################################################################################
# Expression matrix subgrouping
### Author: Angel Garcia-Campos https://github.com/AngelCampos
###############################################################################

#Loading Pathifier DATA
load ("PDS.RData")
PDSmatrix <- mapply(FUN = c, PDS$scores)
PDSmatrix <- t(PDSmatrix)

# Clustering Methods
col.distance = dist(t(PDSmatrix), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")

# Subgroup labeling of samples
mycl <- cutree(col.cluster, k= 6)

# Loading complete expression matrix
exp_mat <- read.delim(file = "ARACNE_matriz_de_expresion_reducida.txt",
	row.names = 1)

# Write expression matrices files in format "exp_X.txt" from 0-6
for(i in 1:6) {
  s <- mycl == i
  em <- as.matrix(exp_mat[,s])
  write.table(em, paste("exp_",i,".txt", sep = ""), sep = "\t", col.names = NA)
}

# Merge healthy individuals expression matrix
exp_mat1 <- read.delim("exp_1.txt", row.names = 1)
exp_mat1 <- exp_mat1[,-19] # RM Tumoral sample that clustets with normal samples
exp_mat2 <- read.delim("exp_2.txt", row.names = 1)
healthy <- cbind(exp_mat1, exp_mat2)

write.table(healthy, "exp_0.txt", sep = "\t", col.names = NA)