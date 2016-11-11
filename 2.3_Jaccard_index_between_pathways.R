################################################################################
# Scripts para calcular el Jaccard index de los pathways de los procesos de 
# apoptosis y autofagia
# Autor: Miguel Angel García Campos Github:/AngelCampos
# With help of wonderful: Guillermo de Anda Jauregui GH:/guillermodeandajauregui
################################################################################

# Load apoptosis and autophagy pathways in gmt format
apop <- read.delim(file = "PATHIFIER_genesets_apop.gmt", header = F)
autop <- read.delim(file = "PATHIFIER_genesets_autop.gmt", header = F)

# Generate functions Jaccar & Jaccard-Matrix (by Guillermo de Anda)
jaccard<-function(a,b)
{
  x<-intersect(a,b)
  y<-union(a,b)
  nx<-length(x)
  ny<-length(y)
  J<-as.numeric(nx/ny)
  return(J)
}
jaccard_matrix <- function(x,y){
  matriz <- matrix(nrow = length(x), ncol = length(y))  
  colnames(matriz) <- names(y)
  rownames(matriz) <- names(x)
  for (i in seq_along(x)){
    for (j in seq_along(y)){
      alfa <- x[[i]]
      beta <- y[[j]]
      matriz[i,j] <- jaccard(alfa,beta)
    }
  }
  return(matriz)  
}

# Prepare pathways
labels_apo <- substring(apop[,1], 1,10); rownames(apop) <- labels_apo
labels_aut <- substring(autop[,1], 1,10); rownames(autop) <- labels_aut

apoptosis <- list()
for (i in rownames(apop)){
  x <- as.character(as.matrix(apop[i, c(-1,-2)]))
  x <- x[!is.na(x)]
  x <- x[x != ""]
  if (length(x) >= 3) {apoptosis[[i]] <- x}
}

autofagia <- list()
for (i in rownames(autop)){
  x <- as.character(as.matrix(autop[i, c(-1,-2)]))
  x <- x[!is.na(x)]
  x <- x[x != ""]
  if (length(x) >= 3) {autofagia[[i]] <- x}
  rm(x)
}

# Calculate Jaccard matrices (Apo-Aut, apo-apo, aut-aut) & draw in PNG format
# Apoptosis vs Autophagy
apo_aut <- jaccard_matrix(apoptosis, autofagia)
png(filename = "J_apo-aut.png", width = 2000*1.618, height = 2000, pointsize =35)
heatmap(apo_aut, scale = "none"); dev.off()
#Apoptosis vs Apoptosis
apo_apo <- jaccard_matrix(apoptosis, apoptosis)
png(filename = "J_apo-apo.png", width = 2000*1.618, height = 2000, pointsize =35)
heatmap(apo_apo, scale = "none", Rowv = NA, Colv = NA); dev.off()
# Autophagy vs autophagy
aut_aut <- jaccard_matrix(autofagia, autofagia)
png(filename = "J_aut-aut.png", width = 2000*1.618, height = 2000, pointsize =35)
heatmap(aut_aut, scale = "none"); dev.off()

# Genes Apoptosis y autofagia

autG <- unique(unlist(autofagia))
apoG <- unique(unlist(apoptosis))
