################################################################################
# ADJs to Square matrix
# Based ond wonderful: Diana Drago https://github.com/dianadrago
################################################################################

# Install/load packages
if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}
if (!require("minet")) {
  biocLite("minet", ask =FALSE)
  library(minet)
}
if (!require("igraph")) {
  install.packages("igraph", dependencies = TRUE)
  library(igraph)
}

# Getting file names
adjs <- sort(list.files(pattern="*.adj"))

# Number of ADJ must be equal to number ofgenes
print(paste("Number of ADJ to be coerced =", length(adjs)))

# Reads the last line from file take values as characters
# converts to matrix and assign first col as rownames 

adjtosqmtx <- function(g){
  a <- scan(file = g, skip = 17, what = "")
  m <- matrix(unlist(a[-(1)]), ncol = 2, byrow = T)
  rownames(m) <- c(m[,1])
  m <- t(m[,2])
  m
}

# Merges ADJS in square matrix and fills empty spaces with NA
M <- lapply(adjs, adjtosqmtx)
x <- t(rbind.fill.matrix(M))

# Replace all NA with 0
x[is.na(x)] <- 0

# Turn character matrix to numeric
class(x) <- "numeric"

# Set colnames as the gene name from ADJs name
names <- gsub('.{6}$','', adjs)
colnames(x) <- names

# Select rownames as colnames and orders them alphabetically
x <- x[sort(rownames(x)), sort(colnames(x)) ]

# Control: class must be 'matrix'
# Matrix must be symmetric
print(paste("Is object matrix? =", class(x))
print(paste("Is matrix square? =", isSymmetric(x))

# save squared matrix
if isSymmetric(x) {
write.table(x, file = "SQnetwork.txt", sep = "\t", col.names= NA, row.names= T,
            quote = F )
}