################################################################################
## Differential expression analysis MA: HGU133-A
### Author: Angel Garcia-Campos https://github.com/AngelCampos
################################################################################

################################################################################
### 1) Install and/or load required packages
################################################################################

source("http://bioconductor.org/biocLite.R")

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("limma")) {
  install.packages("limma", dependencies = TRUE)
  library(limma)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

################################################################################
### 2) Load expression data
################################################################################
# Load expression data
e.data <- read.delim(file = "EXPMATRIX_SANOS_ENFERMOS_COLAPSED.txt", 
                     row.names = 1)

# Load Subgroup data
subg <- read.delim(file = "Muestras_por_Subgrupos.txt", header = F)
subg <- as.matrix(t(subg))
subg <- as.vector(subg)
names(subg) <- colnames(e.data)

################################################################################
### 3) Differential Expression Analysis (LIMMA) & Volcano plots
################################################################################

design = matrix(rep(0,880*5), nrow=880) #Design matrix
colnames(design) = c("Healthy", "Subgroup_1", "Subgroup_2", "Subgroup_3"
                     , "Subgroup_4")
for (i in 1:5){
  design[subg == i,i] <- 1 # Genera la matriz de diseño a partir de la lista de subgrupos
}
# Muestra 870 esta en el subgrupo de sanos!!!!!!!!
# Remover muestra 870
e.data <- e.data[,-870]
subg <- subg[-870]
design <- design[-870,]

# Matriz de contrastes
cont.matrix = makeContrasts(Healthy - Subgroup_1, Healthy - Subgroup_2, 
                            Healthy - Subgroup_3, Healthy - Subgroup_4,
                            levels = design)
fit = lmFit(e.data,design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

#Get the results for every contrast
for (i in 1:4){
  name <- paste("tt", i, sep = "")
  assign(name, topTable(fit2, coef = i, adjust = "fdr", number = 14862))
}

contrasts = c("Healthy - Subgroup_1", "Healthy - Subgroup_2", 
              "Healthy - Subgroup_3", "Healthy - Subgroup_4")

coef = c("Healthy - Subgroup_1", "Healthy - Subgroup_2", 
         "Healthy - Subgroup_3", "Healthy - Subgroup_4")

# Volcano plots PDF report
pdf(file="volcano_plots.pdf", height=7, width=10)

# Parametros gráficos
M = 1; B = 5; y0 = -10; y1 = 800; x0 = -5; x1 = 5

for (i in 1:4){
  volcanoplot(fit2, col = "blue", ylim = c(y0,y1), xlim = c(x0,x1),
              coef = i, main = contrasts[i], cex.lab = 1.3)
  par(new = T)
  abline(v = -M, col = "darkgray")
  par(new = T)
  abline(v = M, col="darkgray")
  par(new = T)
  abline(h = B, col = "black")
  par(new = T)
  ind1 = abs(fit2$coef[,i])>M
  ind2 = fit2$lods[,i] >B
  ind3 = (fit2$coef[,i]>M & fit2$lods[,i]>B)
  ind4 = (fit2$coef[,i]< -M & fit2$lods[,i]>B)
  x = as.matrix(fit2$coef[ind1,i])
  y = as.matrix(fit2$lods[ind1,i])
  plot(x, y, col="magenta", ylim = c(y0,y1), xlim = c(x0,x1),
       main = "", pch = "*", xlab = "", ylab = "", cex.lab = 1.2)
  x = as.matrix(fit2$coef[ind2,i])
  y = as.matrix(fit2$lods[ind2,i])
  par(new=T)
  plot(x, y, col = "orange", ylim = c(y0,y1), xlim = c(x0,x1),
       main = "", pch = 19, xlab = "", ylab = "", cex.lab = 1.2)
  x = as.matrix(fit2$coef[ind3,i])
  y = as.matrix(fit2$lods[ind3,i])
  par(new = T)
  plot(x, y, col = "red",  ylim =c (y0,y1), xlim = c(x0,x1),
       main = "", pch = 19, xlab = "", ylab = "",cex.lab = 1.2)
  x = as.matrix(fit2$coef[ind4,i])
  y = as.matrix(fit2$lods[ind4,i])
  par(new = T)
  plot(x, y, col = "green3",  ylim = c(y0,y1), xlim = c(x0,x1),
       main = "", pch = 19, xlab = "", ylab = "", cex.lab = 1.2)
}
dev.off()

## Write results
write.table(tt1, file = "Dif_exp.Healthy-SG1.txt", quote = F, col.names = NA,
            row.names = T, sep = "\t" )
write.table(tt2, file = "Dif_exp.Healthy-SG2.txt", quote = F, col.names = NA,
            row.names = T, sep = "\t" )
write.table(tt3, file = "Dif_exp.Healthy-SG3.txt", quote = F, col.names = NA,
            row.names = T, sep = "\t" )
write.table(tt4, file = "Dif_exp.Healthy-SG4.txt", quote = F, col.names = NA,
            row.names = T, sep = "\t" )