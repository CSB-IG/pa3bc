################################################################################
# Input variables for the subtype prediction script (PAM50)
# Original author: JS Parker et al., J Clin Onc 27(8):1160-7 (2009)
# Original code: https://genome.unc.edu/pubsup/breastGEO/PAM50.zip
#
################################################################################

# Install/load packages
source("http://bioconductor.org/biocLite.R")
if (!require("ctc")) {
  biocLite("ctc", ask = FALSE)
  library(ctc)
}
if (!require("heatmap.plus")) {
  install.packages("heatmap.plus", dependencies = TRUE)
  library(heatmap.plus)
}
if (!require("impute")) {
  biocLite("impute", ask = FALSE)
  library(impute)
}
library(utils)

# Download and unzip PAM50 files
temp <- tempfile()
download.file("https://genome.unc.edu/pubsup/breastGEO/PAM50.zip", temp)
unzip(temp)

# the location of unchanging files such as the function library and main program
paramDir<- paste(getwd(), "/bioclassifier_R", sep= "/")
# the location of the data matrix, and where output will be located
inputDir<- getwd() 

# the input data matrix as a tab delimited text file
inputFile<- "PAM50_matriz_de_expresion_reducida.txt" 
# short name that will be used for output files
short<-"880_BC_SAMPLES" 

calibrationParameters<- NA 	# the column of the "mediansPerDataset.txt" file to 
							# use for calibration; NA will force centering 
							# within the test set & -1 will not do any
							# adjustment (when adjustment performed by used)

hasClinical<-FALSE 	# may include tumor size as second row, with 'T' as the gene
					# name, and encoded as binary (0 for size <= 2cm or 1 for 
					# size > 2cm) set this variable to FALSE if tumor size is 
					# not available

collapseMethod<-"mean" 	# can be mean or iqr (probe with max iqr is selected)
						# typically, mean is preferred for long oligo and
						# iqr is preferred for short oligo platforms


################################################################################
# run the assignment algorithm
################################################################################

source(paste(paramDir,"subtypePrediction_functions.R",sep="/"))
source(paste(paramDir,"subtypePrediction_distributed.R",sep="/"))

################################################################################
# Convert PAM50 results to numeric vector (better handling for heatmap)
# Normal = 1, Basal = 2, LumA = 3, LumB = 4, Her2 = 5
x <- read.table(file = "subtypes_pam50.txt", header = T, sep = "\t",
				row.names = 1)
p <- as.matrix(x)

subt <- c("1","2","3","4","5")
subt <- c(1,2,3,4,5)

q <- replace(x = p ,list = "Normal" == p, values = subt[1])
q <- replace(x = q ,list = "Basal" == p, values = subt[2])
q <- replace(x = q ,list = "LumA" == p, values = subt[3])
q <- replace(x = q ,list = "LumB" == p, values = subt[4])
q <- replace(x = q ,list = "Her2" == p, values = subt[5])

q <- t(as.matrix((as.numeric(q))))
colnames(q) <- colnames(p)
rownames(q) <- rownames(p)

# Write as tab separated TXT file
write.table(x = q, file = "PAM50subt_numeric.txt", sep = "\t")

#Remplazo de subtipos
# Normal = 1, Basal = 2, LumA = 3, LumB = 4, Her2 = 5
x <- read.table(file = "subtypes_pam50.txt", header = T, sep = "\t", row.names = 1)
p <- as.matrix(x)

subt <- c("1","2","3","4","5")
subt <- c(1,2,3,4,5)

q <- replace(x = p ,list = "Normal" == p, values = subt[1])
q <- replace(x = q ,list = "Basal" == p, values = subt[2])
q <- replace(x = q ,list = "LumA" == p, values = subt[3])
q <- replace(x = q ,list = "LumB" == p, values = subt[4])
q <- replace(x = q ,list = "Her2" == p, values = subt[5])

q <- t(as.matrix((as.numeric(q))))
colnames(q) <- colnames(p)
rownames(q) <- rownames(p)

# Write as tab separated TXT file
write.table(x = q, file = "PAM50subt_numeric.txt", sep = "\t")
