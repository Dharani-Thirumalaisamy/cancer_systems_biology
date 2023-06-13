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
### GDSC DATASET 
gdsc_e <- read_csv("GDSC-CCLE/molecular_data_allcommons/sanger_exp_modified.csv", col_names = TRUE)
gdsc_dr <- read_csv("Spring'23-Cancer Systems Biology/Project/gdsc_pancancer_commondrugs_ic50.csv", col_names = TRUE)

## GDSC_all_drugs_pan-cancer_learning

## GDSC
gdsc_e_complete <- gdsc_e %>% select_if(~!any(is.na(.))) 
gdsc_dr_complete <- gdsc_dr[complete.cases(gdsc_dr), ]

gdsc_e_complete$rn <- toupper(gsub("[.]","",gdsc_e_complete$rn))
gdsc_e_complete$rn <- toupper(gsub("-","",gdsc_e_complete$rn))
gdsc_e_complete$rn <- toupper(gsub(" ","",gdsc_e_complete$rn))
gdsc_e_complete$rn <- toupper(gsub("[()]","",gdsc_e_complete$rn))

names(gdsc_e_complete) <- toupper(gsub("[.]","",names(gdsc_e_complete)))
names(gdsc_e_complete) <- toupper(gsub("-","",names(gdsc_e_complete)))
names(gdsc_e_complete) <- toupper(gsub(" ","",names(gdsc_e_complete)))
names(gdsc_e_complete) <- toupper(gsub("[()]","",names(gdsc_e_complete)))

gdsc_dr_complete$RN <- toupper(gsub("[.]","",gdsc_dr_complete$RN))
gdsc_dr_complete$RN <- toupper(gsub("-","",gdsc_dr_complete$RN))
gdsc_dr_complete$RN <- toupper(gsub(" ","",gdsc_dr_complete$RN))
gdsc_dr_complete$RN <- toupper(gsub("[()]","",gdsc_dr_complete$RN))

common_rows <- as.data.frame(Reduce(intersect, list(gdsc_e_complete$RN, gdsc_dr_complete$RN)))
colnames(common_rows) <- "rn"

colnames(gdsc_e_complete)[1] <- "rn"
colnames(gdsc_dr_complete)[1] <- "rn"

common_filtered_exp <- left_join(common_rows, gdsc_e_complete, by="rn")
#list(common_rows, gdsc_e_complete) %>% reduce(left_join, by="rn") 
common_filtered_drug <- left_join(common_rows, gdsc_dr_complete, by="rn")
#list(common_rows, gdsc_dr_complete) %>% reduce(left_join, by="rn")


# scaling the data 

gdsc_features_exp <- data_scaling(as.data.frame(common_filtered_exp))
rownames(gdsc_features_exp) <- as.matrix(common_rows)
gdsc_exp_complete <- gdsc_features_exp[, colSums(is.na(gdsc_features_exp))==0]


gdsc_features_drug <- data_scaling(as.data.frame(common_filtered_drug))
rownames(gdsc_features_drug) <- as.matrix(common_rows)
gdsc_drug_complete <- gdsc_features_drug[complete.cases(gdsc_features_drug), ]


if(!exists("GBGFAexperiment")) {
  source("gbgfa/GBGFA.R")
}

opts <- getDefaultOpts()
DATANAME <- "e__ic50_gdsc_24_pancancer_commondrugs"
K <- 24
OUTPUTDIR <- "cancer_systems/pancancer/commondrugs"

if (sum(is.na(gdsc_drug_complete)) == 0 && 
    sum(is.na(gdsc_exp_complete)) == 0 && sum(is.nan(gdsc_exp_complete)) == 0) 
{
  print("Model Learning...")
  model <- GBGFAexperiment(list(gdsc_exp_complete, gdsc_drug_complete), 24, opts)
}

concatNames <- paste(DATANAME,K,sep="_")
filename = paste(OUTPUTDIR, concatNames,".Rdata", sep="")

record <- NULL
record$params <- concatNames
record$model <- model

save(record, file = filename)
