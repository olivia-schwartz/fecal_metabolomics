options(install.packages.compile.from.source="never")

if (!require("pacman")) install.packages("pacman") #Installing pacman if not present
pacman::p_load("tidyverse", "KODAMA", "devtools") # install and load the necessary packages

# this 'bisoreg' package installations works well with RStudio instead of Jupyter Notebook
devtools::install_github("cran/bisoreg") 

Directory <- normalizePath(readline("Enter the path of the folder with input files on the output box: "),"/",mustWork=FALSE)
setwd(Directory)

# reading the data from GitHub 
# ft_md_url <- 'https://raw.githubusercontent.com/Functional-Metabolomics-Lab/FBMN-STATS/main/R/outputs_R_Notebook/others/2023-09-07_Ft_md_merged.csv'

# Read the merged feature table with metadata
# Each row represents a sample with columns being metadata attributes and features
ft_merged <-  read.csv("ft_md_merge_batch_correct.csv", header = TRUE, check.names = FALSE) 
corrected_data_title <- "20241105_fecal_correct"

#Putting X at the beginning of each feature
colnames(ft_merged)[15:ncol(ft_merged)] <- paste("X",colnames(ft_merged)[15:ncol(ft_merged)])
colnames(ft_merged)<-gsub(" ","",colnames(ft_merged)) #Remove space
colnames(ft_merged)

head(ft_merged, n=2) #returns the first 2 rows of md
dim(ft_merged) #returns the number of rows and columns

ft_merged <- ft_merged[,-1] #excluding the 1st column

filename_col <- 'Row.names'
batch_info_col <- 'ATTRIBUTE_Batch'
qc_info_col <- 'ATTRIBUTE_Type'
injection_col <- 'ATTRIBUTE_Injection_order'

colnames(ft_merged)[1:20]

# Renaming columns based on user input
ft_merged <- ft_merged %>%
  rename(
    "filename" = !!sym(filename_col),
    "ATTRIBUTE_Batch" = !!sym(batch_info_col),
    "ATTRIBUTE_QCInfo" = !!sym(qc_info_col),
    "ATTRIBUTE_Injection_order" = !!sym(injection_col)
  )

table(ft_merged$ATTRIBUTE_QCInfo)

# Extract relevant columns: filename, batch information, and feature intensity columns
ft_merged2 <- ft_merged %>% 
  select(`filename`, `ATTRIBUTE_Batch`, starts_with("X")) 
head(ft_merged2, n=2)
colnames(ft_merged2)[1:20]

# Step 1: Calculate the overall mean for each feature
# This excludes the first two columns (Row.names, ATTRIBUTE_Batch)
feature_means <- colMeans(ft_merged2[, -(1:2)])

# Step 2: Calculate batchwise feature means
batch_means <- ft_merged2[, -1] %>%  # excluding Row.names column
  group_by(`ATTRIBUTE_Batch`) %>%
  summarise_all(mean) %>%
  column_to_rownames('ATTRIBUTE_Batch') %>%
  as.data.frame()

# Step 3: Correct for inter-batch effects
batch_groups <- ft_merged2 %>%
  group_split(`ATTRIBUTE_Batch`) %>%
  lapply(function(x) {
    column_to_rownames(x, 'filename')
  })

# Check the dimensions of each batch group
sapply(batch_groups, dim)

# List to store corrected data
corrected_batches <- list()

# Apply correction for each batch group
for (i in 1:length(batch_groups)){
  corrected_batches[[i]] <- sweep(batch_groups[[i]][,-1], 2, as.numeric(batch_means[i,] + 1), "/")
  corrected_batches[[i]] <- sweep(corrected_batches[[i]], 2, as.numeric(feature_means + 1), "*")
}

# Combine all corrected batches into a single data frame
corrected_data <- bind_rows(corrected_batches)

# Now, the dataframe "corrected_data" is corrected for inter-batch effects, 
# making it easier to correct for intra-batch effects in the subsequent steps.
 
colnames(corrected_data)<-gsub("X","",colnames(corrected_data)) #Remove space
write.csv(corrected_data, paste0(corrected_data_title,".csv"))

data_trans <- corrected_data %>%
  t() %>%
  as.data.frame()
library(tibble)
data_trans <- tibble::rownames_to_column(data_trans, "ID")

data_col <- select(data_trans,c(1))
data_col <- data_col %>%
  separate_wider_delim(ID, delim = "_", names = c("row ID", "row m/z","row retention time"), cols_remove = FALSE)

merge_trans <- data_col %>% 
  inner_join(.,data_trans)

merge_trans <- select(merge_trans,-c(4))

write.csv(merge_trans, paste0(corrected_data_title,"_trans.csv"))

