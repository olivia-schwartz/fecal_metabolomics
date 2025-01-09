
library(tidyverse)

quant_table <- read.csv("C:/Users/Olivia.Schwartz/OneDrive - University of Denver/DU/Aron Lab/Experiments/20250106_NAU_HILIC/mzmine/20250107_iimn_gnps_quant.csv", sep = ",") #mzmine output
metadata <- read.csv('C:/Users/Olivia.Schwartz/OneDrive - University of Denver/DU/Aron Lab/Experiments/20250106_NAU_HILIC/mzmine/20250106_NAU_HILIC_md_18O.csv', header = T, check.names = F, sep = ",")

colnames(quant_table)<-gsub(".Peak.area","",colnames(quant_table)) #Remove ".Peak.area." 

names <- c(colnames(quant_table)[1:13],metadata$filename)

new_quant <- quant_table[,c(colnames(quant_table) %in% names)]

setdiff(metadata$filename, colnames(new_quant))
setdiff(colnames(new_quant), metadata$filename)

colnames(new_quant)[14:ncol(new_quant)] <- paste(colnames(new_quant)[14:ncol(new_quant)], " Peak area", sep = "")

write.csv(quant_table, file="C:/Users/Olivia.Schwartz/OneDrive - University of Denver/DU/Aron Lab/Experiments/20250106_NAU_HILIC/mzmine/20250107_iimn_gnps_quant_18O.csv")