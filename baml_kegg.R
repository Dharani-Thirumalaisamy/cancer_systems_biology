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

baml_e <- read_csv("exacloud_folder/baml_e_kegg.csv", col_names = TRUE)
baml_dr <- read_csv("exacloud_folder/baml_commondrugs_ic50.csv", col_names = TRUE)

baml_e_complete <- baml_e %>% select_if(~!any(is.na(.)))
baml_dr_complete <- baml_dr[complete.cases(baml_dr), ]

baml_e_complete$rn <- toupper(gsub("[.]","",baml_e_complete$rn))
baml_e_complete$rn <- toupper(gsub("-","",baml_e_complete$rn))
baml_e_complete$rn <- toupper(gsub(" ","",baml_e_complete$rn))
baml_e_complete$rn <- toupper(gsub("[()]","",baml_e_complete$rn))

names(baml_e_complete) <- toupper(gsub("[.]","",names(baml_e_complete)))
names(baml_e_complete) <- toupper(gsub("-","",names(baml_e_complete)))
names(baml_e_complete) <- toupper(gsub(" ","",names(baml_e_complete)))
names(baml_e_complete) <- toupper(gsub("[()]","",names(baml_e_complete)))

baml_dr_complete$RN <- toupper(gsub("[.]","",baml_dr_complete$RN))
baml_dr_complete$RN <- toupper(gsub("-","",baml_dr_complete$RN))
baml_dr_complete$RN <- toupper(gsub(" ","",baml_dr_complete$RN))
baml_dr_complete$RN <- toupper(gsub("[()]","",baml_dr_complete$RN))

common_rows <- as.data.frame(Reduce(intersect, list(baml_e_complete$RN, baml_dr_complete$RN)))
colnames(common_rows) <- "rn"

colnames(baml_e_complete)[1] <- "rn"
colnames(baml_dr_complete)[1] <- "rn"

baml_e_complete <- baml_e_complete[, !duplicated(colnames(baml_e_complete))]

common_filtered_exp <- left_join(common_rows, baml_e_complete, by="rn")
#list(common_rows, baml_e_complete) %>% reduce(left_join, by="rn")

common_filtered_drug <- left_join(common_rows, baml_dr_complete, by="rn")
#list(common_rows, baml_dr_complete) %>% reduce(left_join, by="rn")

# scaling the data

baml_features_exp <- data_scaling(as.data.frame(common_filtered_exp))
rownames(baml_features_exp) <- as.matrix(common_rows)
baml_exp_complete <- baml_features_exp[, colSums(is.na(baml_features_exp))==0]

baml_features_drug <- data_scaling(as.data.frame(common_filtered_drug))
rownames(baml_features_drug) <- as.matrix(common_rows)
baml_drug_complete <- baml_features_drug[complete.cases(baml_features_drug), ]


if(!exists("GBGFAexperiment")) {
  source("gbgfa/GBGFA.R")
}

opts <- getDefaultOpts()
DATANAME <- "e__ic50_baml_24_kegg"
K <- 24
OUTPUTDIR <- "trained_models/aml/"

if (sum(is.na(baml_drug_complete)) == 0 &&
    sum(is.na(baml_exp_complete)) == 0 && sum(is.nan(baml_exp_complete)) == 0)
{
  print("Model Learning...")
  model <- GBGFAexperiment(list(baml_exp_complete, baml_drug_complete), 24, opts)
}

concatNames <- paste(DATANAME,K,sep="_")
filename = paste(OUTPUTDIR, concatNames,".Rdata", sep="")

record <- NULL
record$params <- concatNames
record$model <- model

save(record, file = filename)