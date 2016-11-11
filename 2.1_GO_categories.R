################################################################################
# Merge GO ids of apoptosis and autophagy biological process, and remove cell-
# type specific apoptosis processes
# Author: Miguel Angel Garcia-Campos - https://github.com/AngelCampos
################################################################################

# Load/install packages
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

# Load GO ids (apoptosis, autophagy, cell-type specific apoptosis)
apop <- read.table(file = "GO.apoptotic_process.txt", header = F, sep = "\t", 
           row.names = 1, as.is = T)
autop <- read.table(file = "GO.autophagy.txt", header = F, sep = "\t", 
                   row.names = 1, as.is = T)
cts.apop <- read.table(file = "GO.cell-type specific apoptotic process.txt",
                       header = F, sep = "\t", row.names = 1, as.is = T)

# Remove CTSA GO ids
just_apop <- apop[!rownames(apop) %in% rownames(cts.apop), ]

# Write pathways of apoptosis WITHOUT CTS apoptosis TXT file
write.table(x = just_apop, file = "GO.apoptotic_process_WITHOUT_ctsapop.txt",
            sep = "\t", row.names = T, col.names = F)

# Merge resulting matrices
apop_autop <- bind_rows(just_apop, autop)

# Name rows with GO ids
x <- rownames(just_apop)
y <- rownames(autop)
z <- c(x,y)
rownames(apop_autop) <- z

# Write final set of GO ids to a TXT file
write.table(x = apop_autop, file = "GO.apoptotic_process+autophagy.txt",
            sep = "\t", row.names = T, col.names = F)