# CCLE 
# GDSC
# BetaAML
# LeeAML
# Tavor AML
# swedish and FIMM - not available 

# libraries 
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(tidyverse)
library(synapser)
library(data.table)
library(PharmacoGx)
library(readr)
library(purrr)

library(MASS)
library(dplyr)
library(tidyr)
library(skimr)
#library(tidyverse)
library(data.table)
library(readr)
#library(synapser)
library(purrr)
library(janitor)


# EXPRESSION DATASET

# CCLE dataset 
#ccle_e <- downloadPSet("CCLE_2015", timeout = 24000)

# GDSC dataset 
#gdsc_e <- downloadPSet("GDSC_2020(v2-8.2)", timeout = 24000)

# BAML dataset
baml_pset_e <- downloadPSet("BeatAML_2018", timeout = 24000)

# Tavor dataset
#tavor_e <- downloadPSet("Tavor_2020", timeout = 24000)
#tavor_e_df <- as.data.frame(tavor_e@molecularProfiles$rnaseq@assays@data@listData$expression)
#tavor_e_df_t <- as.data.frame(t(tavor_e_df))

# LeeAML dataset 
#leeaml_e <- read_excel("LeeAML.xlsx", sheet=1)
#leeaml_e <- leeaml_e[2:18811, ]
#names(leeaml_e) <- leeaml_e[1, ]
#leeaml_e <- leeaml_e[-1, ]
#leeaml_e <- leeaml_e[, -2]

#leeaml_e_t <- as.data.frame(t(leeaml_e))
#names(leeaml_e_t) <- leeaml_e_t[1, ]
#leeaml_e_t <- leeaml_e_t[-1, ]

# DRUG SENSITIVITY DATASET 

# CCLE dataset 
#ccle_e <- downloadPSet("CCLE_2015", timeout = 24000)
#x <- read_csv("Expression_Public_22Q4.csv")

# GDSC dataset 
#gdsc_e <- downloadPSet("GDSC_2020(v2-8.2)", timeout = 24000)

# BAML dataset
#baml_e <- downloadPSet("BeatAML_2018", timeout = 24000)
baml_e <- readRDS("Spring'23-Cancer Systems Biology/Project/PSet_BeatAML.rds")

# Tavor dataset
#tavor_e <- downloadPSet("Tavor_2020", timeout = 24000)
#tavor_dr_response <- tavor_e@treatmentResponse$profiles$ic50_recomputed
#tavor_df_compound <- 


# LeeAML dataset 
#leeaml_dr <- read_excel("LeeAML.xlsx", sheet=3)
#leeaml_dr <- leeaml_dr[, -32]

#leeaml_dr_t <- as.data.frame(t(leeaml_dr))
#names(leeaml_dr_t) <- leeaml_dr_t[1, ]
#leeaml_dr_t <- leeaml_dr_t[-1, ]

##### LeeAML has some NA values in the data set. Delete that. Also, format the column names. 
##### Tavor data -> e has 43 cell lines but the drug data has 2668 cell lines. 

#################################################################################

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

ccle_e <- read_csv("gene_expression_data.csv", col_names = TRUE)
ccle_dr <- read_csv("drug_response_data.csv", col_names = TRUE)
ccle_mut <- read_csv("mutation_data.csv", col_names = TRUE)
ccle_cp <- read_csv("Spring'23-CancerSystemsBiology/Project/ccle_pan_cp.csv", col_names = TRUE) #readRDS(synGet("syn3525238")$path)

ccle_mapping <- read_csv("Spring'23-CancerSystemsBiology/Project/Expression_Public_22Q4.csv", col_names = TRUE)
# map with lineage 2 


### GDSC DATASET 
gdsc_e <- read_csv("GDSC-CCLE/molecular_data_allcommons/sanger_exp_modified.csv", col_names = TRUE)
gdsc_dr <- read_excel("GDSC-CCLE/NIHMS81104-supplement-DataS1.xlsx", sheet = 1)
gdsc_mut <- read_csv("GDSC-CCLE/molecular_data_allcommons/mutation_sanger_modified.csv", col_names = TRUE)
gdsc_cp <- read_csv("Spring'23-CancerSystemsBiology/Project/gdsc_pan_cp.csv", col_names = TRUE)  #readRDS(synGet("syn2343774")$path)

gdsc_mapping <- read_csv("Spring'23-Cancer Systems Biology/Project/SANGER_Cell_Lines_Annot.csv", col_names = TRUE)


### Beta AML 

## EXPRESSION DATASET
x <- as.data.frame(baml_pset_e@molecularProfiles$rnaseq@assays@data@listData$exprs)
ensemble_ids <- as.data.frame(baml_pset_e@molecularProfiles$rnaseq@elementMetadata$gene_id)
gene_symbol <- as.data.frame(baml_pset_e@molecularProfiles$rnaseq@elementMetadata$Symbol)

gene_id_symbol_mapping <- cbind(ensemble_ids, gene_symbol)
colnames(gene_id_symbol_mapping) <- c("rn", "gene")

setDT(x, keep.rownames = TRUE)

x_mapped <- left_join(x, gene_id_symbol_mapping, by="rn")
x_mapped <- x_mapped[, -1]
x_mapped <- x_mapped %>% relocate(gene)
x_mapped <- x_mapped %>% distinct(gene, .keep_all = TRUE)

bml_e <- as.data.frame(t(x_mapped))
bml_e <- bml_e %>% row_to_names(row_number = 1)

setDT(bml_e, keep.rownames = TRUE)
write_csv(bml_e, "Spring'23-Cancer Systems Biology/Project/baml_expression.csv", quote ="none")
baml_e <- read_csv("Spring'23-CancerSystemsBiology/Project/baml_expression.csv", col_names = TRUE)

## DRUG DATASET 
bml_dr_id <- as.data.frame(baml_pset_e@treatmentResponse$info$sampleid) 
bml_dr_name <- as.data.frame(baml_pset_e@treatmentResponse$info$treatmentid) 
bml_dr_ic50 <- as.data.frame(baml_pset_e@treatmentResponse$profiles$ic50_published)

bml_dr_mapped <- cbind(bml_dr_id, bml_dr_name, bml_dr_ic50)
colnames(bml_dr_mapped) <- c("rn", "drug", "ic50")
bml_dr_mapped <- unique(bml_dr_mapped)
bml_dr_mapped$drug <- toupper(gsub("[.]","",bml_dr_mapped$drug))
bml_dr_mapped$drug <- toupper(gsub("-","",bml_dr_mapped$drug))
bml_dr_mapped$drug <- toupper(gsub("[()]","",bml_dr_mapped$drug))
bml_dr_mapped$drug <- toupper(gsub(" ","",bml_dr_mapped$drug))

#bml_dr_mapped$ic50 <- as.numeric(bml_dr_mapped$ic50)

#temp_df <- head(bml_dr_mapped, 20)

bml_dr_wider <- pivot_wider(bml_dr_mapped, names_from = "drug", values_from = "ic50", values_fill = 0)

#bml_dr_ic50_complete <- bml_dr_wider[complete.cases(bml_dr_wider), ]
#bml_dr_ic50_complete <- bml_dr_wider %>% select_if(~ !any(is.na(.))) 

write_csv(bml_dr_wider, "Spring'23-Cancer Systems Biology/Project/baml_ic50_drug_data.csv", quote ="none")

baml_dr_ic50 <- read_csv("Spring'23-Cancer Systems Biology/Project/baml_ic50_drug_data.csv", col_names = TRUE)

## MUTATION DATASET 

baml_mut_rownames <- as.data.frame(baml_e@molecularProfiles$mutationchosen@colData@rownames)
baml_mut_genes <- as.data.frame(baml_e@molecularProfiles$mutationall@rowRanges@elementMetadata@listData$symbol)

################################################################################

### Pan-cancer dataset - ALL DRUGS

## CCLE 
ccle_e_complete <- ccle_e %>% select_if(~ !any(is.na(.))) 
print("Expression data extracted")

ccle_dr_complete <- ccle_dr[complete.cases(ccle_dr), ]
print("Drug data extracted")

ccle_mut_complete <- ccle_mut %>% select_if(~ !any(is.na(.))) 
print("Mutation data extracted")

ccle_cp_complete <- ccle_cp %>% select_if(~ !any(is.na(.))) 
print("Copy number alteration data extracted")

# copy number alteration
setDT(ccle_cp_complete, keep.rownames = TRUE)

names(ccle_cp_complete) <- toupper(gsub("[.]", "", names(ccle_cp_complete)))
names(ccle_cp_complete) <- toupper(gsub("-", "", names(ccle_cp_complete)))
names(ccle_cp_complete) <- toupper(gsub("[()]", "", names(ccle_cp_complete)))
names(ccle_cp_complete) <- toupper(gsub(" ", "", names(ccle_cp_complete)))

ccle_cp_complete$RN <- toupper(gsub("[.]", "", ccle_cp_complete$RN))
ccle_cp_complete$RN <- toupper(gsub("-", "", ccle_cp_complete$RN))
ccle_cp_complete$RN <- toupper(gsub("[()]", "", ccle_cp_complete$RN))
ccle_cp_complete$RN <- toupper(gsub(" ", "", ccle_cp_complete$RN))

ccle_cp_complete <- distinct(ccle_cp_complete)
ccle_cp_complete <- ccle_cp_complete %>% distinct(RN, .keep_all = TRUE)
colnames(ccle_cp_complete)[1] <- "rn"

write_csv(ccle_cp_complete, "Spring'23-CancerSystemsBiology/Project/ccle_pan_cp.csv", quote ="none")
################################################################################

## GDSC_all_drugs_pan-cancer_learning

## GDSC
gdsc_e_complete <- gdsc_e %>% select_if(~!any(is.na(.))) 

gdsc_df <- gdsc_dr[, c(1, 2, 8)]
gdsc_df$unified.name <- toupper(gdsc_df$unified.name)
colnames(gdsc_df) <- c("name", "compounds", "IC50")
gdsc_df$IC50 <- as.numeric(gdsc_df$IC50)
gdsc_df <- gdsc_df[complete.cases(gdsc_df), ]

gdsc_ic50 <- gdsc_df %>% pivot_wider(names_from = "compounds", values_from = "IC50")

gdsc_ic50$name[1:9] <- paste('X', gdsc_ic50$name[1:9], sep="")

gdsc_ic50_rownames <- vector()
temp_list <- as.list(gdsc_ic50$name)

for (i in temp_list){
  gdsc_ic50_rownames[i] <- strsplit(i, "_")[[1]][1]
} 

gdsc_ic50 <- gdsc_ic50 %>% mutate("rn" = c(gdsc_ic50_rownames))
gdsc_ic50 <- gdsc_ic50[, -1]
gdsc_ic50 <- gdsc_ic50 %>% relocate(rn)

gdsc_ic50[277, 1] <- "X2313287" #2313287
gdsc_ic50[278, 1] <- "X647V" #647V

names(gdsc_ic50) <- toupper(gsub("[.]", "", names(gdsc_ic50)))
names(gdsc_ic50) <- toupper(gsub("-", "", names(gdsc_ic50)))
names(gdsc_ic50) <- toupper(gsub("NUTLIN3A","NUTLIN3",names(gdsc_ic50)))
names(gdsc_ic50) <- toupper(gsub("17AAG","X17AAG",names(gdsc_ic50)))

gdsc_ic50_complete <- gdsc_ic50[complete.cases(gdsc_ic50), ]
write_csv(gdsc_ic50_complete, "Spring'23-Cancer Systems Biology/Project/sanger_drug_data_ic50.csv", quote ="none")

gdsc_ic50_complete <- read_csv("Spring'23-Cancer Systems Biology/Project/sanger_drug_data_ic50.csv", col_names = TRUE) 

#### copy number 
gdsc_cp_complete <- gdsc_cp %>% select_if(~ !any(is.na(.))) 
print("Copy number alteration data extracted")

# copy number alteration
setDT(gdsc_cp_complete, keep.rownames = TRUE)

names(gdsc_cp_complete) <- toupper(gsub("[.]", "", names(gdsc_cp_complete)))
names(gdsc_cp_complete) <- toupper(gsub("-", "", names(gdsc_cp_complete)))
names(gdsc_cp_complete) <- toupper(gsub("[()]", "", names(gdsc_cp_complete)))
names(gdsc_cp_complete) <- toupper(gsub(" ", "", names(gdsc_cp_complete)))

gdsc_cp_complete$RN <- toupper(gsub("[.]", "", gdsc_cp_complete$RN))
gdsc_cp_complete$RN <- toupper(gsub("-", "", gdsc_cp_complete$RN))
gdsc_cp_complete$RN <- toupper(gsub("[()]", "", gdsc_cp_complete$RN))
gdsc_cp_complete$RN <- toupper(gsub(" ", "", gdsc_cp_complete$RN))

gdsc_cp_complete <- distinct(gdsc_cp_complete)
gdsc_cp_complete <- gdsc_cp_complete %>% distinct(RN, .keep_all = TRUE)
colnames(gdsc_cp_complete)[1] <- "rn"

write_csv(gdsc_cp_complete, "Spring'23-CancerSystemsBiology/Project/gdsc_pan_cp.csv", quote ="none")


################################################################################
################################################################################

## common_drugs_pan-cancer_learning
names(ccle_dr_complete) <- toupper(gsub("[.]","",names(ccle_dr_complete)))

names(gdsc_ic50_complete) <- toupper(gsub("[.]","",names(gdsc_ic50_complete)))

names(baml_dr_ic50) <- toupper(gsub("[.]","",names(baml_dr_ic50)))
colnames(baml_dr_ic50)[89] <- "X17AAG"
colnames(baml_dr_ic50)[78] <- "AZD6244"
colnames(baml_dr_ic50)[102] <- "AZD0530"
colnames(baml_dr_ic50)[116] <- "PD0332991"
colnames(baml_dr_ic50)[60] <- "PF2341066"

#common_cols_2 <- as.data.frame(Reduce(intersect, list(names(ccle_dr_complete), 
#                                                    names(gdsc_ic50_complete))))

common_cols <- as.data.frame(Reduce(intersect, list(names(ccle_dr_complete), 
                                                    names(gdsc_ic50_complete), 
                                                    names(baml_dr_ic50))))
colnames(common_cols) <- "drug"

## CCLE_common_drugs_pan-cancer_learning
ccle_dr_new = subset(ccle_dr_complete, select = c(common_cols$drug))
write_csv(ccle_dr_new, "Spring'23-Cancer Systems Biology/Project/ccle_pancancer_commondrugs_ic50.csv", quote ="none")

## GDSC_common_drugs_pan-cancer_learning
gdsc_dr_new = subset(gdsc_ic50_complete, select = c(common_cols$drug))
write_csv(gdsc_dr_new, "Spring'23-Cancer Systems Biology/Project/gdsc_pancancer_commondrugs_ic50.csv", quote ="none")

## baml_common_drugs_pan-cancer_learning
baml_dr_new = subset(baml_dr_ic50, select = c(common_cols$drug))
write_csv(baml_dr_new, "Spring'23-Cancer Systems Biology/Project/baml_commondrugs_ic50.csv", quote ="none")

################################################################################

# ONLY AML - ALL DRUGS

### CCLE 
ccle_mapping_df <- ccle_mapping[, c(2, 4)]
ccle_mapping_df$cell_line_display_name[1:23] <- paste('X', ccle_mapping_df$cell_line_display_name[1:23], sep="")
ccle_mapping_df <- ccle_mapping_df %>% filter(lineage_2 == "Acute Myeloid Leukemia")

# exp data
common_lines <- as.data.frame(intersect(ccle_e_complete$rn, ccle_mapping_df$cell_line_display_name))
colnames(common_lines) <- "rn"
ccle_e_new <- left_join(common_lines, ccle_e_complete, by="rn") 
write_csv(ccle_e_new, "Spring'23-Cancer Systems Biology/Project/ccle_aml_alldrugs_e.csv", quote ="none")
#list(common_lines, as.data.frame(ccle_e_complete)) %>% reduce(left_join, by="rn") 

# mut data
common_lines <- as.data.frame(intersect(ccle_mut_complete$rn, ccle_mapping_df$cell_line_display_name))
colnames(common_lines) <- "rn"
ccle_mut_new <- left_join(common_lines, ccle_mut_complete, by="rn") 
write_csv(ccle_mut_new, "Spring'23-Cancer Systems Biology/Project/ccle_aml_h.csv", quote ="none")

# copy number data
common_lines <- as.data.frame(intersect(ccle_cp_complete$RN, ccle_mapping_df$cell_line_display_name))
colnames(common_lines) <- "rn"
colnames(ccle_cp_complete)[1] <- "rn"
ccle_cp_complete <- distinct(ccle_cp_complete)
ccle_cp_complete <- ccle_cp_complete %>% distinct(rn, .keep_all = TRUE)
ccle_cp_new <- left_join(common_lines, ccle_cp_complete, by="rn") 
write_csv(ccle_cp_new, "Spring'23-CancerSystemsBiology/Project/ccle_aml_cp.csv", quote ="none")

### GDSC

# exp data
common_lines <- as.data.frame(intersect(gdsc_e_complete$rn, ccle_mapping_df$cell_line_display_name))
colnames(common_lines) <- "rn"
gdsc_e_new <- left_join(common_lines, gdsc_e_complete, by="rn") 
write_csv(gdsc_e_new, "Spring'23-Cancer Systems Biology/Project/gdsc_aml_alldrugs_e.csv", quote ="none")

# mut data
gdsc_mut_complete <- gdsc_mut %>% select_if(~!any(is.na(.))) 

common_lines <- as.data.frame(intersect(gdsc_mut_complete$rn, ccle_mapping_df$cell_line_display_name))
colnames(common_lines) <- "rn"
gdsc_mut_new <- left_join(common_lines, gdsc_mut_complete, by="rn") 
write_csv(gdsc_mut_new, "Spring'23-Cancer Systems Biology/Project/gdsc_aml_mut.csv", quote ="none")

# copy number data
common_lines <- as.data.frame(intersect(gdsc_cp_complete$RN, ccle_mapping_df$cell_line_display_name))
colnames(common_lines) <- "rn"
colnames(gdsc_cp_complete)[1] <- "rn"
gdsc_cp_complete <- distinct(gdsc_cp_complete)
gdsc_cp_complete <- gdsc_cp_complete %>% distinct(rn, .keep_all = TRUE)
gdsc_cp_new <- left_join(common_lines, gdsc_cp_complete, by="rn") 
write_csv(gdsc_cp_new, "Spring'23-CancerSystemsBiology/Project/gdsc_aml_cp.csv", quote ="none")

################################################################################

# Gene set curated - pathways curated 

#### ALL PATHWAYS
curated_gene <- read_csv("Spring'23-CancerSystemsBiology/Project/MSigDB_curated.csv", col_names = TRUE)
#curated_gene <- curated_gene[, -1]
list_curated_gene <- as.list(curated_gene)
list_curated_gene <- as.vector(unlist(list_curated_gene))
unique_curated_gene <- as.data.frame(unique(list_curated_gene))
#unique_curated_gene <- list_curated_gene[!duplicated(list_curated_gene)]

unique_curated_gene$`unique(list_curated_gene)` <- toupper(gsub("[.]","",unique_curated_gene$`unique(list_curated_gene)`))
unique_curated_gene$`unique(list_curated_gene)` <- toupper(gsub("[()]","",unique_curated_gene$`unique(list_curated_gene)`))
unique_curated_gene$`unique(list_curated_gene)` <- toupper(gsub("-","",unique_curated_gene$`unique(list_curated_gene)`))
unique_curated_gene$`unique(list_curated_gene)` <- toupper(gsub(" ","",unique_curated_gene$`unique(list_curated_gene)`))

#### KEGG CURATED 
curated_gene_kegg <- read_csv("Spring'23-CancerSystemsBiology/Project/MSigDB_kegg.csv", col_names = TRUE)
#curated_gene <- curated_gene[, -1]
list_kegg_gene <- as.list(curated_gene_kegg)
list_kegg_gene <- as.vector(unlist(list_kegg_gene))
unique_kegg_gene <- as.data.frame(unique(list_kegg_gene))
#unique_curated_gene <- list_curated_gene[!duplicated(list_curated_gene)]

unique_kegg_gene$`unique(list_kegg_gene)` <- toupper(gsub("[.]","",unique_kegg_gene$`unique(list_kegg_gene)`))
unique_kegg_gene$`unique(list_kegg_gene)` <- toupper(gsub("[()]","",unique_kegg_gene$`unique(list_kegg_gene)`))
unique_kegg_gene$`unique(list_kegg_gene)` <- toupper(gsub("-","",unique_kegg_gene$`unique(list_kegg_gene)`))
unique_kegg_gene$`unique(list_kegg_gene)` <- toupper(gsub(" ","",unique_kegg_gene$`unique(list_kegg_gene)`))

#### REACTOME CURATED 
curated_gene_reactome <- read_csv("Spring'23-CancerSystemsBiology/Project/MSigDB_reactome.csv", col_names = TRUE)
#curated_gene <- curated_gene[, -1]
list_reactome_gene <- as.list(curated_gene_reactome)
list_reactome_gene <- as.vector(unlist(list_reactome_gene))
unique_reactome_gene <- as.data.frame(unique(list_reactome_gene))
#unique_curated_gene <- list_curated_gene[!duplicated(list_curated_gene)]

unique_reactome_gene$`unique(list_reactome_gene)` <- toupper(gsub("[.]","",unique_reactome_gene$`unique(list_reactome_gene)`))
unique_reactome_gene$`unique(list_reactome_gene)` <- toupper(gsub("[()]","",unique_reactome_gene$`unique(list_reactome_gene)`))
unique_reactome_gene$`unique(list_reactome_gene)` <- toupper(gsub("-","",unique_reactome_gene$`unique(list_reactome_gene)`))
unique_reactome_gene$`unique(list_reactome_gene)` <- toupper(gsub(" ","",unique_reactome_gene$`unique(list_reactome_gene)`))

#### BAML 
### exp curation
names(baml_e) <- toupper(gsub("-", "", names(baml_e)))
names(baml_e) <- toupper(gsub("[.]", "", names(baml_e)))
names(baml_e) <- toupper(gsub("[()]", "", names(baml_e)))
names(baml_e) <- toupper(gsub(" ", "", names(baml_e)))

baml_e <- as.data.frame(baml_e)
colnames(baml_e)[1] <- "rn"

rownames(baml_e) <- baml_e[, 1]
baml_e <- baml_e[, -1]

baml_e_cols <- as.data.frame(names(baml_e))
#bml_e_cols <- bml_e_cols[-1, ]

bml_e_cols_all <- as.data.frame(intersect(baml_e_cols$`names(baml_e)`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(bml_e_cols_all) <- "genes"
bml_e_curated = subset(baml_e, select = bml_e_cols_all$genes)
setDT(bml_e_curated, keep.rownames = TRUE)
write_csv(bml_e_curated, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/baml_e_mgidb.csv", quote = "none")

bml_e_cols_kegg <- as.data.frame(intersect(baml_e_cols$`names(baml_e)`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(bml_e_cols_kegg ) <- "genes"
bml_e_curated_kegg = subset(baml_e, select = bml_e_cols_kegg$genes)
setDT(bml_e_curated_kegg, keep.rownames = TRUE)
write_csv(bml_e_curated_kegg, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/baml_e_kegg.csv", quote = "none")

bml_e_cols_reactome <- as.data.frame(intersect(baml_e_cols$`names(baml_e)`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(bml_e_cols_reactome ) <- "genes"
bml_e_curated_reactome = subset(baml_e, select = bml_e_cols_reactome$genes)
setDT(bml_e_curated_reactome, keep.rownames = TRUE)
write_csv(bml_e_curated_reactome, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/baml_e_reactome.csv", quote = "none")

#### CCLE 

### exp curation
names(ccle_e_complete) <- toupper(gsub("-", "", names(ccle_e_complete)))
names(ccle_e_complete) <- toupper(gsub("[.]", "", names(ccle_e_complete)))
names(ccle_e_complete) <- toupper(gsub("[()]", "", names(ccle_e_complete)))
names(ccle_e_complete) <- toupper(gsub(" ", "", names(ccle_e_complete)))

ccle_e_complete <- as.data.frame(ccle_e_complete)
colnames(ccle_e_complete)[1] <- "rn"
#df <- ccle_e_complete
rownames(ccle_e_complete) <- ccle_e_complete[, 1]
ccle_e_complete <- ccle_e_complete[, -1]

ccle_e_cols <- as.data.frame(names(ccle_e_complete))
#ccle_e_cols <- ccle_e_cols[-1, ]

common_e_cols <- as.data.frame(intersect(ccle_e_cols$`names(ccle_e_complete)`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_e_cols) <- "genes"
ccle_e_curated = subset(ccle_e_complete, select = common_e_cols$genes)
setDT(ccle_e_curated, keep.rownames = TRUE)
write_csv(ccle_e_curated, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_exp_mgidb.csv", quote = "none")

common_e_cols_kegg <- as.data.frame(intersect(ccle_e_cols$`names(ccle_e_complete)`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_e_cols_kegg ) <- "genes"
ccle_e_curated_kegg = subset(ccle_e_complete, select = common_e_cols_kegg$genes)
setDT(ccle_e_curated_kegg, keep.rownames = TRUE)
write_csv(ccle_e_curated_kegg, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_exp_kegg.csv", quote = "none")

common_e_cols_reactome <- as.data.frame(intersect(ccle_e_cols$`names(ccle_e_complete)`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_e_cols_reactome ) <- "genes"
ccle_e_curated_reactome = subset(ccle_e_complete, select = common_e_cols_reactome$genes)
setDT(ccle_e_curated_reactome, keep.rownames = TRUE)
write_csv(ccle_e_curated_reactome, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_exp_reactome.csv", quote = "none")

#### exp aml curation
ccle_e_aml <- read_csv("Spring'23-CancerSystemsBiology/Project/ccle_aml_alldrugs_e.csv", col_names = TRUE)

names(ccle_e_aml) <- toupper(gsub("-", "", names(ccle_e_aml)))
names(ccle_e_aml) <- toupper(gsub("[.]", "", names(ccle_e_aml)))
names(ccle_e_aml) <- toupper(gsub("[()]", "", names(ccle_e_aml)))
names(ccle_e_aml) <- toupper(gsub(" ", "", names(ccle_e_aml)))

ccle_e_aml <- as.data.frame(ccle_e_aml)
colnames(ccle_e_aml)[1] <- "rn"

rownames(ccle_e_aml) <- ccle_e_aml[, 1]
ccle_e_aml <- ccle_e_aml[, -1]

ccle_e_cols_aml <- as.data.frame(names(ccle_e_aml))
#ccle_e_cols_aml <- as.data.frame(ccle_e_cols_aml[-1, ])

common_e_cols_aml <- as.data.frame(intersect(ccle_e_cols_aml$`names(ccle_e_aml)`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_e_cols_aml) <- "genes"
ccle_e_curated_aml = subset(ccle_e_aml, select = common_e_cols_aml$genes)
setDT(ccle_e_curated_aml, keep.rownames = TRUE)
write_csv(ccle_e_curated_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_exp_aml_mgidb.csv", quote = "none")

common_e_cols_kegg_aml <- as.data.frame(intersect(ccle_e_cols_aml$`names(ccle_e_aml)`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_e_cols_kegg_aml) <- "genes"
ccle_e_curated_kegg_aml = subset(ccle_e_aml, select = common_e_cols_kegg_aml$genes)
setDT(ccle_e_curated_kegg_aml, keep.rownames = TRUE)
write_csv(ccle_e_curated_kegg_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_exp_aml_kegg.csv", quote = "none")

common_e_cols_reactome_aml <- as.data.frame(intersect(ccle_e_cols_aml$`names(ccle_e_aml)`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_e_cols_reactome_aml) <- "genes"
ccle_e_curated_reactome_aml = subset(ccle_e_aml, select = common_e_cols_reactome_aml$genes)
setDT(ccle_e_curated_reactome_aml, keep.rownames = TRUE)
write_csv(ccle_e_curated_reactome_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_exp_aml_reactome.csv", quote = "none")


### mut curation
names(ccle_mut_complete) <- toupper(gsub("-", "", names(ccle_mut_complete)))
names(ccle_mut_complete) <- toupper(gsub("[.]", "", names(ccle_mut_complete)))
names(ccle_mut_complete) <- toupper(gsub("[()]", "", names(ccle_mut_complete)))
names(ccle_mut_complete) <- toupper(gsub(" ", "", names(ccle_mut_complete)))

ccle_mut_complete <- as.data.frame(ccle_mut_complete)
colnames(ccle_mut_complete)[1] <- "rn"

rownames(ccle_mut_complete) <- ccle_mut_complete[, 1]
ccle_mut_complete <- ccle_mut_complete[, -1]

ccle_mut_cols <- as.data.frame(names(ccle_mut_complete))
#ccle_mut_cols <- as.data.frame(ccle_mut_cols[-1, ])

common_mut_cols <- as.data.frame(intersect(ccle_mut_cols$`names(ccle_mut_complete)`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_mut_cols) <- "genes"
ccle_mut_curated = subset(ccle_mut_complete, select = common_mut_cols$genes)
setDT(ccle_mut_curated, keep.rownames = TRUE)
write_csv(ccle_mut_curated, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_mut_mgidb.csv", quote = "none")

common_mut_cols_kegg <- as.data.frame(intersect(ccle_mut_cols$`names(ccle_mut_complete)`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_mut_cols_kegg ) <- "genes"
ccle_mut_curated_kegg = subset(ccle_mut_complete, select = common_mut_cols_kegg$genes)
setDT(ccle_mut_curated_kegg, keep.rownames = TRUE)
write_csv(ccle_mut_curated_kegg, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_mut_kegg.csv", quote = "none")

common_mut_cols_reactome <- as.data.frame(intersect(ccle_mut_cols$`names(ccle_mut_complete)`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_mut_cols_reactome ) <- "genes"
ccle_mut_curated_reactome = subset(ccle_mut_complete, select = common_mut_cols_reactome$genes)
setDT(ccle_mut_curated_reactome, keep.rownames = TRUE)
write_csv(ccle_mut_curated_reactome, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_mut_reactome.csv", quote = "none")


#### mut aml curation
ccle_mut_aml <- read_csv("Spring'23-CancerSystemsBiology/Project/ccle_aml_mut.csv", col_names = TRUE)

names(ccle_mut_aml) <- toupper(gsub("-", "", names(ccle_mut_aml)))
names(ccle_mut_aml) <- toupper(gsub("[.]", "", names(ccle_mut_aml)))
names(ccle_mut_aml) <- toupper(gsub("[()]", "", names(ccle_mut_aml)))
names(ccle_mut_aml) <- toupper(gsub(" ", "", names(ccle_mut_aml)))

ccle_mut_aml <- as.data.frame(ccle_mut_aml)
colnames(ccle_mut_aml)[1] <- "rn"

rownames(ccle_mut_aml) <- ccle_mut_aml[, 1]
ccle_mut_aml <- ccle_mut_aml[, -1]

ccle_mut_cols_aml <- as.data.frame(names(ccle_mut_aml))
#ccle_mut_cols_aml <- as.data.frame(ccle_mut_cols_aml[-1, ])

common_mut_cols_aml <- as.data.frame(intersect(ccle_mut_cols_aml$`names(ccle_mut_aml)`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_mut_cols_aml) <- "genes"
ccle_mut_curated_aml = subset(ccle_mut_aml, select = common_mut_cols_aml$genes)
setDT(ccle_mut_curated_aml, keep.rownames = TRUE)
write_csv(ccle_mut_curated_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_mut_aml_mgidb.csv", quote = "none")

common_mut_cols_kegg_aml <- as.data.frame(intersect(ccle_mut_cols_aml$`names(ccle_mut_aml)`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_mut_cols_kegg_aml) <- "genes"
ccle_mut_curated_kegg_aml = subset(ccle_mut_aml, select = common_mut_cols_kegg_aml$genes)
setDT(ccle_mut_curated_kegg_aml, keep.rownames = TRUE)
write_csv(ccle_mut_curated_kegg_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_mut_aml_kegg.csv", quote = "none")

common_mut_cols_reactome_aml <- as.data.frame(intersect(ccle_mut_cols_aml$`names(ccle_mut_aml)`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_mut_cols_reactome_aml) <- "genes"
ccle_mut_curated_reactome_aml = subset(ccle_mut_aml, select = common_mut_cols_reactome_aml$genes)
setDT(ccle_mut_curated_reactome_aml, keep.rownames = TRUE)
write_csv(ccle_mut_curated_reactome_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_mut_aml_reactome.csv", quote = "none")


### copy number curation
names(ccle_cp_complete) <- toupper(gsub("-", "", names(ccle_cp_complete)))
names(ccle_cp_complete) <- toupper(gsub("[.]", "", names(ccle_cp_complete)))
names(ccle_cp_complete) <- toupper(gsub("[()]", "", names(ccle_cp_complete)))
names(ccle_cp_complete) <- toupper(gsub(" ", "", names(ccle_cp_complete)))

colnames(ccle_cp_complete)[1] <- "rn"
ccle_cp_cols <- as.data.frame(names(ccle_cp_complete))
ccle_cp_cols <- as.data.frame(ccle_cp_cols[-1, ])

common_cp_cols <- as.data.frame(intersect(ccle_cp_cols$`ccle_cp_cols[-1, ]`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_cp_cols) <- "genes"
ccle_cp_curated = subset(ccle_cp_complete, select = common_cp_cols$genes)
write_csv(ccle_cp_curated, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_cp_mgidb.csv", quote = "none")

common_cp_cols_kegg <- as.data.frame(intersect(ccle_cp_cols$`ccle_cp_cols[-1, ]`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_cp_cols_kegg ) <- "genes"
ccle_cp_curated_kegg = subset(ccle_cp_complete, select = common_cp_cols_kegg$genes)
write_csv(ccle_cp_curated_kegg, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_cp_kegg.csv", quote = "none")

common_cp_cols_reactome <- as.data.frame(intersect(ccle_cp_cols$`ccle_cp_cols[-1, ]`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_cp_cols_reactome ) <- "genes"
ccle_cp_curated_reactome = subset(ccle_cp_complete, select = common_cp_cols_reactome$genes)
write_csv(ccle_cp_curated_reactome, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_cp_reactome.csv", quote = "none")


#### copy number aml curation
ccle_cp_aml <- read_csv("Spring'23-CancerSystemsBiology/Project/ccle_aml_cp.csv", col_names = TRUE)

names(ccle_cp_aml) <- toupper(gsub("-", "", names(ccle_cp_aml)))
names(ccle_cp_aml) <- toupper(gsub("[.]", "", names(ccle_cp_aml)))
names(ccle_cp_aml) <- toupper(gsub("[()]", "", names(ccle_cp_aml)))
names(ccle_cp_aml) <- toupper(gsub(" ", "", names(ccle_cp_aml)))

colnames(ccle_cp_aml)[1] <- "rn"
ccle_cp_cols_aml <- as.data.frame(names(ccle_cp_aml))
ccle_cp_cols_aml <- as.data.frame(ccle_cp_cols_aml[-1, ])

common_cp_cols_aml <- as.data.frame(intersect(ccle_cp_cols_aml$`ccle_cp_cols_aml[-1, ]`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_cp_cols_aml) <- "genes"
ccle_cp_curated_aml = subset(ccle_cp_aml, select = common_cp_cols_aml$genes)
write_csv(ccle_cp_curated_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_cp_aml_mgidb.csv", quote = "none")

common_cp_cols_kegg_aml <- as.data.frame(intersect(ccle_cp_cols_aml$`ccle_cp_cols_aml[-1, ]`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_cp_cols_kegg_aml) <- "genes"
ccle_cp_curated_kegg_aml = subset(ccle_cp_aml, select = common_cp_cols_kegg_aml$genes)
write_csv(ccle_cp_curated_kegg_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_cp_aml_kegg.csv", quote = "none")

common_cp_cols_reactome_aml <- as.data.frame(intersect(ccle_cp_cols_aml$`ccle_cp_cols_aml[-1, ]`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_cp_cols_reactome_aml) <- "genes"
ccle_cp_curated_reactome_aml = subset(ccle_cp_aml, select = common_cp_cols_reactome_aml$genes)
write_csv(ccle_cp_curated_reactome_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_cp_aml_reactome.csv", quote = "none")


#### GDSC 
### exp curation
names(gdsc_e_complete) <- toupper(gsub("-", "", names(gdsc_e_complete)))
names(gdsc_e_complete) <- toupper(gsub("[.]", "", names(gdsc_e_complete)))
names(gdsc_e_complete) <- toupper(gsub("[()]", "", names(gdsc_e_complete)))
names(gdsc_e_complete) <- toupper(gsub(" ", "", names(gdsc_e_complete)))

gdsc_e_complete <- as.data.frame(gdsc_e_complete)
colnames(gdsc_e_complete)[1] <- "rn"

rownames(gdsc_e_complete) <- gdsc_e_complete[, 1]
gdsc_e_complete <- gdsc_e_complete[, -1]

gdsc_e_cols <- as.data.frame(names(gdsc_e_complete))
#gdsc_e_cols <- as.data.frame(gdsc_e_cols[-1, ])

common_e_cols_gdsc <- as.data.frame(intersect(gdsc_e_cols$`names(gdsc_e_complete)` , unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_e_cols_gdsc) <- "genes"
gdsc_e_curated = subset(gdsc_e_complete, select = common_e_cols_gdsc$genes)
setDT(gdsc_e_curated, keep.rownames = TRUE)
write_csv(gdsc_e_curated, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_exp_mgidb.csv", quote = "none")

common_e_cols_kegg_gdsc <- as.data.frame(intersect(gdsc_e_cols$`names(gdsc_e_complete)`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_e_cols_kegg_gdsc ) <- "genes"
gdsc_e_curated_kegg = subset(gdsc_e_complete, select = common_e_cols_kegg_gdsc$genes)
setDT(gdsc_e_curated_kegg, keep.rownames = TRUE)
write_csv(gdsc_e_curated_kegg, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_e_kegg.csv", quote = "none")

common_e_cols_reactome_gdsc <- as.data.frame(intersect(gdsc_e_cols$`names(gdsc_e_complete)`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_e_cols_reactome_gdsc ) <- "genes"
gdsc_e_curated_reactome = subset(gdsc_e_complete, select = common_e_cols_reactome_gdsc$genes)
setDT(gdsc_e_curated_reactome, keep.rownames = TRUE)
write_csv(gdsc_e_curated_reactome, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_e_reactome.csv", quote = "none")

#### exp aml curation
gdsc_e_aml <- read_csv("Spring'23-CancerSystemsBiology/Project/gdsc_aml_alldrugs_e.csv", col_names = TRUE)

names(gdsc_e_aml) <- toupper(gsub("-", "", names(gdsc_e_aml)))
names(gdsc_e_aml) <- toupper(gsub("[.]", "", names(gdsc_e_aml)))
names(gdsc_e_aml) <- toupper(gsub("[()]", "", names(gdsc_e_aml)))
names(gdsc_e_aml) <- toupper(gsub(" ", "", names(gdsc_e_aml)))

gdsc_e_aml <- as.data.frame(gdsc_e_aml)
colnames(gdsc_e_aml)[1] <- "rn"

rownames(gdsc_e_aml) <- gdsc_e_aml[, 1]
gdsc_e_aml <- gdsc_e_aml[, -1]

gdsc_e_cols_aml <- as.data.frame(names(gdsc_e_aml))
#gdsc_e_cols_aml <- as.data.frame(gdsc_e_cols_aml[-1, ])

common_e_cols_aml_gdsc <- as.data.frame(intersect(gdsc_e_cols_aml$`names(gdsc_e_aml)` , unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_e_cols_aml_gdsc) <- "genes"
gdsc_e_curated_aml = subset(gdsc_e_aml, select = common_e_cols_aml_gdsc$genes)
setDT(gdsc_e_curated_aml, keep.rownames = TRUE)
write_csv(gdsc_e_curated_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_exp_aml_mgidb.csv", quote = "none")

common_e_cols_kegg_aml_gdsc <- as.data.frame(intersect(gdsc_e_cols_aml$`names(gdsc_e_aml)`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_e_cols_kegg_aml_gdsc) <- "genes"
gdsc_e_curated_kegg_aml = subset(gdsc_e_aml, select = common_e_cols_kegg_aml_gdsc$genes)
setDT(gdsc_e_curated_kegg_aml, keep.rownames = TRUE)
write_csv(gdsc_e_curated_kegg_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_e_aml_kegg.csv", quote = "none")

common_e_cols_reactome_aml_gdsc <- as.data.frame(intersect(gdsc_e_cols_aml$`names(gdsc_e_aml)`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_e_cols_reactome_aml_gdsc) <- "genes"
gdsc_e_curated_reactome_aml = subset(gdsc_e_aml, select = common_e_cols_reactome_aml_gdsc$genes)
setDT(gdsc_e_curated_reactome_aml, keep.rownames = TRUE)
write_csv(gdsc_e_curated_reactome_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_e_aml_reactome.csv", quote = "none")

### mut curation

names(gdsc_mut_complete) <- toupper(gsub("-", "", names(gdsc_mut_complete)))
names(gdsc_mut_complete) <- toupper(gsub("[.]", "", names(gdsc_mut_complete)))
names(gdsc_mut_complete) <- toupper(gsub("[()]", "", names(gdsc_mut_complete)))
names(gdsc_mut_complete) <- toupper(gsub(" ", "", names(gdsc_mut_complete)))

gdsc_mut_complete <- as.data.frame(gdsc_mut_complete)
colnames(gdsc_mut_complete)[1] <- "rn"

rownames(gdsc_mut_complete) <- gdsc_mut_complete[, 1]
gdsc_mut_complete <- gdsc_mut_complete[, -1]

gdsc_mut_cols <- as.data.frame(names(gdsc_mut_complete))
#gdsc_mut_cols <- as.data.frame(gdsc_mut_cols[-1, ])

common_mut_cols_gdsc <- as.data.frame(intersect(gdsc_mut_cols$`names(gdsc_mut_complete)`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_mut_cols_gdsc) <- "genes"
gdsc_mut_curated_gdsc = subset(gdsc_mut_complete, select = common_mut_cols_gdsc$genes)
setDT(gdsc_mut_curated_gdsc, keep.rownames = TRUE)
write_csv(gdsc_mut_curated_gdsc, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_mut_mgidb.csv", quote = "none")

common_mut_cols_kegg_gdsc <- as.data.frame(intersect(gdsc_mut_cols$`names(gdsc_mut_complete)`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_mut_cols_kegg_gdsc ) <- "genes"
gdsc_mut_curated_kegg = subset(gdsc_mut_complete, select = common_mut_cols_kegg_gdsc$genes)
setDT(gdsc_mut_curated_kegg, keep.rownames = TRUE)
write_csv(gdsc_mut_curated_kegg, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_mut_kegg.csv", quote = "none")

common_mut_cols_reactome_gdsc <- as.data.frame(intersect(gdsc_mut_cols$`names(gdsc_mut_complete)`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_mut_cols_reactome_gdsc ) <- "genes"
gdsc_mut_curated_reactome = subset(gdsc_mut_complete, select = common_mut_cols_reactome_gdsc$genes)
setDT(gdsc_mut_curated_reactome, keep.rownames = TRUE)
write_csv(gdsc_mut_curated_reactome, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_mut_reactome.csv", quote = "none")

#### mut aml curation
gdsc_mut_aml <- read_csv("Spring'23-CancerSystemsBiology/Project/gdsc_aml_mut.csv", col_names = TRUE)

names(gdsc_mut_aml) <- toupper(gsub("-", "", names(gdsc_mut_aml)))
names(gdsc_mut_aml) <- toupper(gsub("[.]", "", names(gdsc_mut_aml)))
names(gdsc_mut_aml) <- toupper(gsub("[()]", "", names(gdsc_mut_aml)))
names(gdsc_mut_aml) <- toupper(gsub(" ", "", names(gdsc_mut_aml)))

gdsc_mut_aml <- as.data.frame(gdsc_mut_aml)
colnames(gdsc_mut_aml)[1] <- "rn"

rownames(gdsc_mut_aml) <- gdsc_mut_aml[, 1]
gdsc_mut_aml <- gdsc_mut_aml[, -1]

gdsc_mut_cols_aml <- as.data.frame(names(gdsc_mut_aml))
#gdsc_mut_cols_aml <- as.data.frame(gdsc_mut_cols_aml[-1, ])

common_mut_cols_aml_gdsc <- as.data.frame(intersect(gdsc_mut_cols_aml$`names(gdsc_mut_aml)`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_mut_cols_aml_gdsc) <- "genes"
gdsc_mut_curated_aml = subset(gdsc_mut_aml, select = common_mut_cols_aml_gdsc$genes)
setDT(gdsc_mut_curated_aml, keep.rownames = TRUE)
write_csv(gdsc_mut_curated_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_mut_aml_mgidb.csv", quote = "none")

common_mut_cols_kegg_aml_gdsc <- as.data.frame(intersect(gdsc_mut_cols_aml$`names(gdsc_mut_aml)`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_mut_cols_kegg_aml_gdsc) <- "genes"
gdsc_mut_curated_kegg_aml = subset(gdsc_mut_aml, select = common_mut_cols_kegg_aml_gdsc$genes)
setDT(gdsc_mut_curated_kegg_aml, keep.rownames = TRUE)
write_csv(gdsc_mut_curated_kegg_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_mut_aml_kegg.csv", quote = "none")

common_mut_cols_reactome_aml_gdsc <- as.data.frame(intersect(gdsc_mut_cols_aml$`names(gdsc_mut_aml)`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_mut_cols_reactome_aml_gdsc) <- "genes"
gdsc_mut_curated_reactome_aml = subset(gdsc_mut_aml, select = common_mut_cols_reactome_aml_gdsc$genes)
setDT(gdsc_mut_curated_reactome_aml, keep.rownames = TRUE)
write_csv(gdsc_mut_curated_reactome_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_mut_aml_reactome.csv", quote = "none")

### copy number curation
names(gdsc_cp_complete) <- toupper(gsub("-", "", names(gdsc_cp_complete)))
names(gdsc_cp_complete) <- toupper(gsub("[.]", "", names(gdsc_cp_complete)))
names(gdsc_cp_complete) <- toupper(gsub("[()]", "", names(gdsc_cp_complete)))
names(gdsc_cp_complete) <- toupper(gsub(" ", "", names(gdsc_cp_complete)))

colnames(gdsc_cp_complete)[1] <- "rn"
gdsc_cp_cols <- as.data.frame(names(gdsc_cp_complete))
gdsc_cp_cols <- as.data.frame(gdsc_cp_cols[-1, ])

common_cp_cols_gdsc <- as.data.frame(intersect(gdsc_cp_cols$`gdsc_cp_cols[-1, ]`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_cp_cols_gdsc) <- "genes"
gdsc_cp_curated_gdsc = subset(gdsc_cp_complete, select = common_cp_cols_gdsc$genes)
write_csv(gdsc_cp_curated_gdsc, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_cp_mgidb.csv", quote = "none")

common_cp_cols_kegg_gdsc <- as.data.frame(intersect(gdsc_cp_cols$`gdsc_cp_cols[-1, ]`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_cp_cols_kegg_gdsc ) <- "genes"
gdsc_cp_curated_kegg = subset(gdsc_cp_complete, select = common_cp_cols_kegg_gdsc$genes)
write_csv(gdsc_cp_curated_kegg, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_cp_kegg.csv", quote = "none")

common_cp_cols_reactome_gdsc <- as.data.frame(intersect(gdsc_cp_cols$`gdsc_cp_cols[-1, ]`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_cp_cols_reactome_gdsc ) <- "genes"
gdsc_cp_curated_reactome = subset(gdsc_cp_complete, select = common_cp_cols_reactome_gdsc$genes)
write_csv(gdsc_cp_curated_reactome, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_cp_reactome.csv", quote = "none")

#### copy number aml curation
gdsc_cp_aml <- read_csv("Spring'23-CancerSystemsBiology/Project/gdsc_aml_cp.csv", col_names = TRUE)

names(gdsc_cp_aml) <- toupper(gsub("-", "", names(gdsc_cp_aml)))
names(gdsc_cp_aml) <- toupper(gsub("[.]", "", names(gdsc_cp_aml)))
names(gdsc_cp_aml) <- toupper(gsub("[()]", "", names(gdsc_cp_aml)))
names(gdsc_cp_aml) <- toupper(gsub(" ", "", names(gdsc_cp_aml)))

colnames(gdsc_cp_aml)[1] <- "rn"
gdsc_cp_cols_aml <- as.data.frame(names(gdsc_cp_aml))
gdsc_cp_cols_aml <- as.data.frame(gdsc_cp_cols_aml[-1, ])

common_cp_cols_aml_gdsc <- as.data.frame(intersect(gdsc_cp_cols_aml$`gdsc_cp_cols_aml[-1, ]`, unique_curated_gene$`unique(list_curated_gene)`))
colnames(common_cp_cols_aml_gdsc) <- "genes"
gdsc_cp_curated_aml = subset(gdsc_cp_aml, select = common_cp_cols_aml_gdsc$genes)
write_csv(gdsc_cp_curated_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_cp_aml_mgidb.csv", quote = "none")

common_cp_cols_kegg_aml_gdsc <- as.data.frame(intersect(gdsc_cp_cols_aml$`gdsc_cp_cols_aml[-1, ]`, unique_kegg_gene$`unique(list_kegg_gene)`))
colnames(common_cp_cols_kegg_aml_gdsc) <- "genes"
gdsc_cp_curated_kegg_aml = subset(gdsc_cp_aml, select = common_cp_cols_kegg_aml_gdsc$genes)
write_csv(gdsc_cp_curated_kegg_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_cp_aml_kegg.csv", quote = "none")

common_cp_cols_reactome_aml_gdsc <- as.data.frame(intersect(gdsc_cp_cols_aml$`gdsc_cp_cols_aml[-1, ]`, unique_reactome_gene$`unique(list_reactome_gene)`))
colnames(common_cp_cols_reactome_aml_gdsc) <- "genes"
gdsc_cp_curated_reactome_aml = subset(gdsc_cp_aml, select = common_cp_cols_reactome_aml_gdsc$genes)
write_csv(gdsc_cp_curated_reactome_aml, "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_cp_aml_reactome.csv", quote = "none")

################################################################################

# CSV TO RDS 

ccle_aml_e <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_aml_alldrugs_e.csv", col_names = TRUE))
toRDS_ccle_aml_e <- ccle_aml_e[, -1]
rownames(toRDS_ccle_aml_e) <- ccle_aml_e[, 1]
saveRDS(as.data.frame(toRDS_ccle_aml_e), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_aml_e.RDS")

ccle_aml_h <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_aml_mut.csv", col_names = TRUE))
toRDS_ccle_aml_h <- ccle_aml_h[, -1]
rownames(toRDS_ccle_aml_h) <- ccle_aml_h[, 1]
saveRDS(as.data.frame(toRDS_ccle_aml_h), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_aml_h.RDS")

ccle_aml_cp <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_aml_cp.csv", col_names = TRUE))
toRDS_ccle_aml_cp <- ccle_aml_cp[, -1]
rownames(toRDS_ccle_aml_cp) <- ccle_aml_cp[, 1]
saveRDS(as.data.frame(toRDS_ccle_aml_cp), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_aml_cp.RDS")

ccle_pan_cp <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_pan_cp.csv", col_names = TRUE))
toRDS_ccle_pan_cp <- ccle_pan_cp[, -1]
rownames(toRDS_ccle_pan_cp) <- ccle_pan_cp[, 1]
saveRDS(as.data.frame(toRDS_ccle_pan_cp), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_pan_cp.RDS")

ccle_aml_dr <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_pancancer_commondrugs_ic50.csv", col_names = TRUE))
toRDS_ccle_aml_dr <- ccle_aml_dr[, -1]
rownames(toRDS_ccle_aml_dr) <- ccle_aml_dr[, 1]
saveRDS(as.data.frame(toRDS_ccle_aml_dr), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/ccle_aml_dr.RDS")


gdsc_aml_e <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_aml_alldrugs_e.csv", col_names = TRUE))
toRDS_gdsc_aml_e <- gdsc_aml_e[, -1]
rownames(toRDS_gdsc_aml_e) <- gdsc_aml_e[, 1]
saveRDS(as.data.frame(toRDS_gdsc_aml_e), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_aml_e.RDS")

gdsc_aml_h <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_aml_mut.csv", col_names = TRUE))
toRDS_gdsc_aml_h <- gdsc_aml_h[, -1]
rownames(toRDS_gdsc_aml_h) <- gdsc_aml_h[, 1]
saveRDS(as.data.frame(toRDS_gdsc_aml_h), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_aml_h.RDS")

gdsc_aml_cp <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_aml_cp.csv", col_names = TRUE))
toRDS_gdsc_aml_cp <- gdsc_aml_cp[, -1]
rownames(toRDS_gdsc_aml_cp) <- gdsc_aml_cp[, 1]
saveRDS(as.data.frame(toRDS_gdsc_aml_cp), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_aml_cp.RDS")

gdsc_pan_cp <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_pan_cp.csv", col_names = TRUE))
toRDS_gdsc_pan_cp <- gdsc_pan_cp[, -1]
rownames(toRDS_gdsc_pan_cp) <- gdsc_pan_cp[, 1]
saveRDS(as.data.frame(toRDS_gdsc_pan_cp), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_pan_cp.RDS")

gdsc_aml_dr <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_pancancer_commondrugs_ic50.csv", col_names = TRUE))
toRDS_gdsc_aml_dr <- gdsc_aml_dr[, -1]
rownames(toRDS_gdsc_aml_dr) <- gdsc_aml_dr[, 1]
saveRDS(as.data.frame(toRDS_gdsc_aml_dr), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/gdsc_aml_dr.RDS")


baml_e <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/baml_expression.csv", col_names = TRUE))
toRDS_baml_e <- baml_e[, -1]
rownames(toRDS_baml_e) <- baml_e[, 1]
saveRDS(as.data.frame(toRDS_baml_e), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/baml_e.RDS")

baml_dr <- as.data.frame(read_csv("Spring'23-CancerSystemsBiology/Project/exacloud_folder/baml_commondrugs_ic50.csv", col_names = TRUE))
toRDS_baml_dr <- baml_dr[, -1]
rownames(toRDS_baml_dr) <- baml_dr[, 1]
saveRDS(as.data.frame(toRDS_baml_dr), "Spring'23-CancerSystemsBiology/Project/exacloud_folder/baml_dr.RDS")


