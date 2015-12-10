################################################################################
# Hypergeometric tests for testing subgroups against molecular subtypes
### Author: Angel Garcia-Campos https://github.com/AngelCampos
###############################################################################

## Table Rows = Subgroups, Colums = Molecular subtypes
# Ordering data
subtype_ch <- read.delim(file = "subtypes_pam50.txt")
subgroup <- as.data.frame(t(mycl))

# Total de muestras por subtipo
r1 <- cbind(sum(subtype_ch == "basal"), sum(subtype_ch == "lumA"), 
            sum(subtype_ch == "lumB"), sum(subtype_ch == "her2"))

# Obtener numero de muestras por subtipo en cada subgrupo
y1 <- NULL
for (i in c("basal", "lumA", "lumB", "her2", "normal")){
  g <- as.logical(subtype_ch == i) & as.logical(subgroup == 2)
  y <- sum(g)
  y1 <- cbind(y1, y) 
}
y2 <- NULL
for (i in c("basal", "lumA", "lumB", "her2", "normal")){
  g <- as.logical(subtype_ch == i) & as.logical(subgroup == 3)
  y <- sum(g)
  y2 <- cbind(y2, y) 
}
y3 <- NULL
for (i in c("basal", "lumA", "lumB", "her2", "normal")){
  g <- as.logical(subtype_ch == i) & as.logical(subgroup == 4)
  y <- sum(g)
  y3 <- cbind(y3, y) 
}
y4 <- NULL
for (i in c("basal", "lumA", "lumB", "her2", "normal")){
  g <- as.logical(subtype_ch == i) & as.logical(subgroup == 5)
  y <- sum(g)
  y4 <- cbind(y4, y) 
}

sub_table <- rbind(y1, y2, y3, y4)
rownames(sub_table) <- c("SG1","SG2","SG3","SG4")
colnames(sub_table) <- c("Basal","lumA", "lumB", "her2", "normal")
write.table(x = sub_table, file = "Tabla_subtipos_VS_subgrupos.txt", 
            sep = " \t", col.names = NA, quote = F)


### Hypergeometric tests

library(stats)
# m: number of succeses in population
# n: number of failures in population
# k: sample size
# x: number of succeses in sample
pop <- sum(sub_table) # Total population

# basal
pv1 <- NULL
for (i in 1:4){
  m <- sum(sub_table[,1])
  n <- pop - m
  k <- sum(sub_table[i,])
  x <- sub_table[i,1]
  pv <- phyper(x, m, n, k, lower.tail = F)
  pv1 <- c(pv1, pv)
}
# lumA
pv2 <- NULL
for (i in 1:4){
  m <- sum(sub_table[,2])
  n <- pop - m
  k <- sum(sub_table[i,])
  x <- sub_table[i,2]
  pv <- phyper(x, m, n, k, lower.tail = F)
  pv2 <- signif(c(pv2, pv), 5)
}
# lumB
pv3 <- NULL
for (i in 1:4){
  m <- sum(sub_table[,3])
  n <- pop - m
  k <- sum(sub_table[i,])
  x <- sub_table[i,3]
  pv <- phyper(x, m, n, k, lower.tail = F)
  pv3 <- signif(c(pv3, pv), 5)
}
# her2
pv4 <- NULL
for (i in 1:4){
  m <- sum(sub_table[,4])
  n <- pop - m
  k <- sum(sub_table[1,])
  x <- sub_table[i,4]
  pv <- phyper(x, m, n, k, lower.tail = F)
  pv4 <- signif(c(pv4, pv), 5)
}

pv_table <- cbind(pv1, pv2, pv3, pv4)
rownames(pv_table) <- c("SG1","SG2","SG3","SG4")
colnames(pv_table) <- c("basal", "lumA", "lumB", "her2")
write.table(x = pv_table, file = "P-value_table.txt", sep = "\t", 
            col.names = NA, quote = F)
