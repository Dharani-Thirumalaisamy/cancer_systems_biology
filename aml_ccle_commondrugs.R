library(MASS)
library(dplyr)
library(tidyr)
library(skimr)
#library(tidyverse)
library(data.table)
library(readr)
#library(synapser)
library(purrr)

# scaling the train data 
data_scaling <- function(filtered_data){
  
  scaled_filtered_data <- scale(filtered_data[2:ncol(filtered_data)], center = TRUE, scale = TRUE)
  
  return(scaled_filtered_data)
}

print(paste("Start time", Sys.time()))

### Extracting the data from synapse

### CCLE DATASET 

ccle_e <- read_csv("Spring'23-Cancer Systems Biology/Project/ccle_aml_alldrugs_e.csv", col_names = TRUE)
ccle_dr <- read_csv("Spring'23-Cancer Systems Biology/Project/ccle_pancancer_commondrugs_ic50.csv", col_names = TRUE)

ccle_e_complete <- ccle_e %>% select_if(~ !any(is.na(.))) 
print("Expression data extracted")

ccle_dr_complete <- ccle_dr[complete.cases(ccle_dr), ]
print("Drug data extracted")

ccle_e_complete$rn <- toupper(gsub("[.]","",ccle_e_complete$rn))
ccle_e_complete$rn <- toupper(gsub("-","",ccle_e_complete$rn))
ccle_e_complete$rn <- toupper(gsub(" ","",ccle_e_complete$rn))
ccle_e_complete$rn <- toupper(gsub("[()]","",ccle_e_complete$rn))

names(ccle_e_complete) <- toupper(gsub("[.]","",names(ccle_e_complete)))
names(ccle_e_complete) <- toupper(gsub("-","",names(ccle_e_complete)))
names(ccle_e_complete) <- toupper(gsub(" ","",names(ccle_e_complete)))
names(ccle_e_complete) <- toupper(gsub("[()]","",names(ccle_e_complete)))

ccle_dr_complete$RN <- toupper(gsub("[.]","",ccle_dr_complete$RN))
ccle_dr_complete$RN <- toupper(gsub("-","",ccle_dr_complete$RN))
ccle_dr_complete$RN <- toupper(gsub(" ","",ccle_dr_complete$RN))
ccle_dr_complete$RN <- toupper(gsub("[()]","",ccle_dr_complete$RN))

common_rows <- as.data.frame(Reduce(intersect, list(ccle_e_complete$RN, ccle_dr_complete$RN)))
colnames(common_rows) <- "rn"

colnames(ccle_e_complete)[1] <- "rn"
colnames(ccle_dr_complete)[1] <- "rn"

common_filtered_exp <- left_join(common_rows, ccle_e_complete, by="rn")
#list(common_rows, ccle_e_complete) %>% reduce(left_join, by="rn") 
common_filtered_drug <- left_join(common_rows, ccle_dr_complete, by="rn")
#list(common_rows, ccle_dr_complete) %>% reduce(left_join, by="rn")


# scaling the data 
ccle_features_drug <- data_scaling(as.data.frame(common_filtered_drug))
rownames(ccle_features_drug) <- as.matrix(common_rows)
ccle_drug_complete <- ccle_features_drug[complete.cases(ccle_features_drug), ]


ccle_features_exp <- data_scaling(as.data.frame(common_filtered_exp))
rownames(ccle_features_exp) <- as.matrix(common_rows)
ccle_exp_complete <- ccle_features_exp[, colSums(is.na(ccle_features_exp))==0]


if(!exists("GBGFAexperiment")) {
  source("gbgfa/GBGFA.R")
}

opts <- getDefaultOpts()
DATANAME <- "e__ic50_ccle_24_aml_commondrugs"
K <- 24
OUTPUTDIR <- "cancer_systems/aml/commondrugs"

if (sum(is.na(ccle_drug_complete)) == 0 && 
    sum(is.na(ccle_exp_complete)) == 0 && sum(is.nan(ccle_exp_complete)) == 0) 
{
  print("Model Learning...")
  model <- GBGFAexperiment(list(ccle_exp_complete, ccle_drug_complete), 24, opts)
}

concatNames <- paste(DATANAME,K,sep="_")
filename = paste(OUTPUTDIR, concatNames,".Rdata", sep="")

record <- NULL
record$params <- concatNames
record$model <- model

save(record, file = filename)
