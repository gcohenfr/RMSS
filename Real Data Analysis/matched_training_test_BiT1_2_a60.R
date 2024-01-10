rm(list=ls())
library(readxl)
library(dplyr)
library(stringr)
library(readr)

#setwd "ValidationData"
pull <- function(x,y) {x[,if(is.name(substitute(y))) deparse(substitute(y)) else y, drop = FALSE][[1]]}

### Read BiT 1 and 2 patient info
BiT1_BiT2_raw <- read_excel("BiT1 & BiT2 CAV MRM_2013.02.01_review by Mustafa_2015.01.06.xlsx")
BiT1_BiT2_named <- BiT1_BiT2_raw %>% select('Cohort', 'BIOBANK_ID', 'STUDY_ID', 'TIMEPOINT', 'DAYS_SINCE_TX','Max %DS LAD_clean', 'Max %DS All', 'Max %DS LAD, RCA, LCX')
BiT1_BiT2_numbered_1 <- BiT1_BiT2_raw[, 1:6] # Columns A-F in Excel
BiT1_BiT2_numbered_2 <- BiT1_BiT2_raw[, 21:25] # Columns U-Y in Excel
BiT1_BiT2 <- bind_rows(BiT1_BiT2_named, BiT1_BiT2_numbered_1, BiT1_BiT2_numbered_2)
rm(BiT1_BiT2_raw, BiT1_BiT2_named, BiT1_BiT2_numbered_1, BiT1_BiT2_numbered_2)

### Get response info
patients_w_ID <- BiT1_BiT2 %>% select(BIOBANK_ID) %>% na.omit %>% pull(BIOBANK_ID)
# Keep patients with id and response
BiT1_BiT2_id_resp <- BiT1_BiT2 %>% filter(BIOBANK_ID %in% patients_w_ID) %>% select(BIOBANK_ID, STUDY_ID,TIMEPOINT,DAYS_SINCE_TX,"Max %DS LAD_clean") %>% na.omit
# > length(unique(BiT1_BiT2_id_resp$STUDY_ID))
# [1] 82 #82 patients, 100 samples
#Some patients got multiple angios at different timepoints (i.e., multiple responses per patient).
# If a blood sample is taken at the time of the angio, they also have multiple BIOBANK_IDs (one per blood sample). 

# The variable TIMEPOINT (in response file and second part part of BIOBANK_ID in peptide files) is useful to identify 
# samples used in the discovery study.

### Discovery set 
load("../parsed_training_data.rda")
discovery.samples<-data.frame(sample=rownames(x))
discovery.samples$STUDY_ID <- as.character(str_extract(discovery.samples$sample, "[^W]*")) 
discovery.samples$weeks <- as.character(str_extract(discovery.samples$sample, "[W][0-9]{2,3}")) 
discovery.samples$weeks <- as.numeric(str_extract(discovery.samples$weeks, "[0-9]{2,3}"))
discovery.samples$days <- discovery.samples$weeks*7

# Use week in discovery info to match samples 
# weeks converted to days, match to the closest sample in time using "DAYS_SINCE_TX"
BiT1_BiT2_id_resp<-merge(BiT1_BiT2_id_resp,discovery.samples,by="STUDY_ID",all.x = T)
BiT1_BiT2_id_resp$days_diff<-abs(BiT1_BiT2_id_resp$days-BiT1_BiT2_id_resp$DAYS_SINCE_TX)
train_biobank<-BiT1_BiT2_id_resp %>% group_by(STUDY_ID) %>% filter(days_diff==min(days_diff)) %>% 
  select(sample,BIOBANK_ID,STUDY_ID)#all have a response!

### Read peptides for BiT 1 and 2
data_path<-""
peptides_plate1 <- read_excel(paste0(data_path, "BiT1 - BiT3 CAV/2013/20130321_PROOF 0019_CAV plate#1_Accuracy_RR_NAT Resp.xlsx"), sheet = "PRF019_CAV_plt#1_RR")
peptides_plate2 <- read_excel(paste0(data_path, "BiT1 - BiT3 CAV/2013/20130321_PROOF 0019_CAV plate#2_Accuracy_RR_NAT Resp.xlsx"), sheet = "PRF019_CAV_plt#2_RR")
peptides_plate3 <- read_excel(paste0(data_path, "BiT1 - BiT3 CAV/2013/20130321_PROOF 0019_CAV plate#3_Accuracy_RR_NAT Resp.xlsx"), sheet = "PRF019_CAV_plt#3_RR")

num_cols_peptides_file <- ncol(peptides_plate1)
num_obs_peptides_file <- nrow(peptides_plate1)
plates <- 1:3
for (plate in plates){
  initial_name <- paste0('peptides_plate', plate) 
  final_name <- paste0(initial_name, "_final") 
  names <- eval(parse(text=initial_name))[-1, 3] # This is Biobank ID_TIMEPOINT
  names <- as.character(sapply(names, function(x){ str_extract(x, "[^_]*")})) # Keep only first part of the identifier
  peptides <- as.matrix(apply(eval(parse(text=initial_name))[-1, 8:num_cols_peptides_file], 2, as.numeric))# These are the peptides
  colnames(peptides) <- sapply(colnames(peptides), function(x){ str_trim(str_extract(x, "[^\\.]*"))}) # Keep only aminoacid sequence in the names
  assign(final_name, data.frame(BIOBANK_ID = names, peptides, plate = rep(plate, num_obs_peptides_file - 1), stringsAsFactors = FALSE)) # add plate ID
}

BiT1_BiT2_peptides <- bind_rows(peptides_plate1_final, peptides_plate2_final, peptides_plate3_final) #merge everything into one df

rm(peptides_plate1, peptides_plate2, peptides_plate3, peptides_plate1_final,
   peptides_plate2_final, peptides_plate3_final, names, final_name, initial_name, num_cols_peptides_file, num_obs_peptides_file,
   plate, plates)

#--------------------------------------------------------------
#Validation of proteins discovered using weighted PENSE(M)
###### Make feature matrices for BiT1 and 2; filter
#CAV_QC <- read_csv("BiT1 - BiT3 CAV/2013/CAV QC.csv")
#useful_peptides <- CAV_QC$Peptide[!is.na(CAV_QC$MM_peptides)] #peptides associated with proteins selected in the first analysis
#--------------------------------------------------------------

#--------------------------------------------------------------
#Validation of proteins discovered using unweighted PENSE(M)- in paper
###### Make feature matrices for BiT1 and 2; filter
CAV_QC <- read_csv("BiT1 - BiT3 CAV/2013/CAV QC.csv")
useful_peptides <- CAV_QC$Peptide[!is.na(CAV_QC$MM_uw_a60)] #peptides associated with proteins selected in the first analysis


# keep only peptides of BIOBANK IDS with response, and select peptides from 1st study
peptides_w_response <- BiT1_BiT2_peptides %>% filter(BIOBANK_ID %in% BiT1_BiT2_id_resp$BIOBANK_ID) %>%
  select(BIOBANK_ID, useful_peptides, plate)

# My 4-NAs rule is not filtering peptide "ELSHLPSLYDYSAYR" that has very bad quality (0 in almost all samples),
# and other peptides with very low signal. Manually inspected: 10 left, 2 more can be cosidered later (QC=2)
# There is one NA left in this set for peptide LIDQYGTHYLQSGSLGGEYR, which can be imputed with the minimum value (or 0.002349774,value of the other peptide below LOD)) 

#Based on correlations between peptides of the same protein peptide LGEVNTYAGDLQK from PGC23 
# and TNQVNSGGVLLR of PGC103 are left out (QC=2)

good_qc_pept <- CAV_QC %>% filter(QC_60 ==1) %>% select(Peptide, PGC) #check this!
peptides_w_response <- peptides_w_response %>% select(BIOBANK_ID, good_qc_pept$Peptide, plate) #10 peptides left

# See which peptides correspond to which proteins
match_pept_prot <- CAV_QC %>% select(Peptide, PGC) %>% na.omit %>% filter(Peptide %in% colnames(peptides_w_response)) #check this!
PGC <- match_pept_prot %>% pull(PGC) %>% unique %>% sort
proteins <- c()
for (pgc in PGC){
  which_pept <- match_pept_prot %>% filter(PGC==pgc) %>% pull(Peptide)
  aux_pept <- peptides_w_response %>% select(which_pept) %>% apply(., 1, mean, na.rm=TRUE) #each protein is the average of its peptides
  proteins <- cbind(proteins, aux_pept)
}
colnames(proteins) <- PGC
rownames(proteins) <- peptides_w_response %>% pull(BIOBANK_ID) #100 samples, 7 proteins!

# Replace NaN by column minimum
index_nans <- apply(proteins, 2, function(x) { nans <- which(is.nan(x)) })
col_mins <- apply(proteins, 2, min, na.rm=TRUE)
for (j in 1:length(index_nans)){
  if (length(index_nans[j]) > 0){
    bad_row <- index_nans[[j]]
  }
  proteins[bad_row, j] <- col_mins[j]
}

# Add response to the feature matrix 
response_all<-BiT1_BiT2_id_resp %>% select(sample,BIOBANK_ID,`Max %DS LAD_clean`)
proteins<-data.frame(BIOBANK_ID=row.names(proteins),proteins)
dat_all<-merge(proteins,response_all,by="BIOBANK_ID")

train<-dat_all %>% filter(BIOBANK_ID %in% train_biobank$BIOBANK_ID)

#Test: keep samples from patients not included in the training set
test_biobank<- BiT1_BiT2_id_resp %>% filter(!STUDY_ID %in% train_biobank$STUDY_ID) %>% 
  select(BIOBANK_ID)

test<-dat_all %>% filter(BIOBANK_ID %in% test_biobank$BIOBANK_ID)

save(file='trainBiT1_testBiT2_a60.Rdata', list=c('train', 'test','train_biobank'))