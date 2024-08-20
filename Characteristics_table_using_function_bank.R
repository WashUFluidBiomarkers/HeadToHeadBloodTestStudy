#### Characteristics Table
#### Coded: Benjamin Saef
#### Date as of commenting and documentation - 05/28/2024

library(tidyr)
library(data.table)
library(RVAideMemoire)
library(pROC)
library(RColorBrewer)
library(MuMIn)
source("./Functions_for_Analysis.R")


##Set up project folders and 
out_project <- "/FNIH_Project_2/Paper1" #Folder where you want all output written to, a results and plot folder will be written here by default
out_results <- paste0(out_project,"/results/")
out_plots <- paste0(out_project,"/plots/")

dir.create(out_project,showWarnings = F)
dir.create(out_results,showWarnings = F)
dir.create(out_plots,showWarnings = F)
##Read files 
data_file <- "/FNIH_Project_2/DATA_STUDY_2/cross_sectional_data_2024_04_23.csv" 

last_crossectional <- read.csv(data_file)

#Fixed up race variable coding
last_crossectional <- as.data.frame(last_crossectional)
last_crossectional$RACE <- ifelse(last_crossectional$RACE==2,3,last_crossectional$RACE)
last_crossectional$RACE <- ifelse(last_crossectional$RACE==0,2,last_crossectional$RACE)

#Create Interval Variables

last_crossectional$PLASMA_CDR_INT <- as.numeric(abs(difftime(last_crossectional$EXAMDATE,last_crossectional$VISDATE_CDR  ,units = "days"))/365.25)
last_crossectional$PLASMA_TAU_INT <- as.numeric(abs(difftime(last_crossectional$EXAMDATE,last_crossectional$SCANDATE_TAU,units = "days"))/365.25)
last_crossectional$PLASMA_AMY_INT <- as.numeric(abs(difftime(last_crossectional$EXAMDATE,last_crossectional$SCANDATE_AMY  ,units = "days"))/365.25)
last_crossectional$PLASMA_MRI_INT <- as.numeric(abs(difftime(last_crossectional$EXAMDATE,last_crossectional$SCANDATE_ATROPHY  ,units = "days"))/365.25)
last_crossectional$PLASMA_CSF_INT <- as.numeric(abs(difftime(last_crossectional$EXAMDATE,last_crossectional$EXAMDATE_CSF  ,units = "days"))/365.25)
last_crossectional[which(is.na(last_crossectional$PTAU_over_ABETA42)),"PLASMA_CSF_INT"] <- NA


#list of variable IDs to create characteristic table entries for 
varb_1 <- c("AGE","PTGENDER","APOE_genotype","APOE4_positive","PTEDUCAT","CDR_SOB","PLASMA_CDR_INT","MMSCORE","RACE",
            "ABETA42","TAU","PTAU","PTAU_over_ABETA42",#CSF
            "PLASMA_CSF_INT",
            "C2N_plasma_Abeta42","C2N_plasma_Abeta40","C2N_plasma_Abeta42_Abeta40","C2N_plasma_ptau217","C2N_plasma_ptau217_ratio","C2N_plasma_nptau217",#C2N
            "Roche_plasma_Ab42","Roche_plasma_Ab40","Roche_plasma_Ab42_Ab40","Roche_plasma_ptau181","Roche_plasma_GFAP","Roche_plasma_NfL",#Elecsys
            "Fuji_plasma_Ab42","Fuji_plasma_Ab40","Fuji_plasma_Ab42_Ab40","Fuji_plasma_ptau217",#Lumipulse
            "AlzPath_plasma_ptau217" ,"Janssen_plasma_ptau217",#Simoa
            "QX_plasma_Ab42","QX_plasma_Ab40","QX_plasma_Ab42_Ab40","QX_plasma_ptau181","QX_plasma_GFAP","QX_plasma_NfL", #Quanterix
            "PLASMA_AMY_INT",
            "CENTILOIDS","CENTILOIDS_10_sens",
            "PLASMA_TAU_INT",
            "MesialTemporal","TAU_MesialTemporal_10",
            "TemporoParietal","TAU_TemporoParietal_10",
            "atrophy","atrophy_10","PLASMA_MRI_INT") #Imaging Variables

#list of variable types to go in that table - please see Functions_for_Analysis.R to get more detailed breakdown of types
type_1 <- c("A","P","4","P","A","A","A","A","R",
            "A","A","A","A",#CSF
            "A",
            "A","A","A","A","A","A",
            "A","A","A","A","A","A",         
            "A","A","A","A",
            "A","A",
            "A","A","A","A","A","A",
            "A",
            "A","P",
            "A",
            "A","P",
            "A","P",
            "A","P","A")

#list of labels for the variables to go into the table
label_1 <- c("Age (Years)","Sex (% Female)","APOE4 Genotype (22/23/24/33/34/44)","APOE4 (% Carrier)","Years of Education",
             "CDR Sum of Boxes","Plasma Collection to CDR Evaluation Time Interval","Mini-Mental State Exam Score","Race (Black/White/Other)",
              "CSF Ab42","CSF Total Tau","CSF p-Tau 181","CSF p-Tau 181 / Ab42",#CSF
             "Time between Plasma Collection and CSF Collection",
             "Plasma Ab42","Plasma Ab40","Plasma Ab42/Ab40","Plasma p-tau217","Plasma p-tau217 Ratio","np-tau217 (pg/ml)", #C2N
             "Plasma Ab42","Plasma Ab40","Plasma Ab42/Ab40","Plasma p-Tau 181","Plasma GFAP","Plasma Neurofilament Light",#Elecsys
             "Plasma Ab42","Plasma Ab40","Plasma Ab42/Ab40","Plasma p-tau217 Concentration",#Lumipulse
             "Plasma p-tau217 Concentration","Plasma p-tau217 Concentration",#Simoa
             "Plasma Ab42","Plasma Ab40","Plasma Ab42/Ab40","Plasma p-Tau 181","Plasma GFAP","Plasma Neurofilament Light", #Quanterix
             "Time between Plasma Collection and Amyloid PET imaging",
             "Amyloid PET Centiloid","Amyloid PET Positivity (Centiloid > 37)",
             "Time between Plasma Collection and Tau PET imaging",
             "Tau PET Early","Tau PET Early Positivity",
             "Tau PET Late","Tau PET Late Positivity",
             "Cortical signature volume","Brain Atrophy Positivity","Time between Plasma Collection and MRI")



####Function call - Function is in Functions_for_Analysis.R - Create full Data characteristics tables
char_table_full_bin(last_crossectional,
                    "CENTILOIDS_10",
                    "Amyloid PET",
                    varb_1,
                    type_1,
                    label_1,
                    paste0(out_results,"Characteristics_Table_by_Centiloid20_positivity_05212024.xlsx"))


##Didn't make a function specifically for this, so decided to make this row separately

last_crossectional$CDR_CAT <- NA
last_crossectional[which(last_crossectional$CDR==0),"CDR_CAT"] <- "0"
last_crossectional[which(last_crossectional$CDR==0.5),"CDR_CAT"] <- "0.5"
last_crossectional[which(last_crossectional$CDR>0.5),"CDR_CAT"] <- "1+"

chisq.test(last_crossectional[,"CENTILOIDS_10"], 
           as.factor(last_crossectional[,"CDR_CAT"]), 
           correct=FALSE)


####
## CSF complete cases - characteristics tables
####
data_file_csf_only <- "/FNIH_Project_2/DATA_STUDY_2/cross_sectional_data_2024_04_30b_CSF_ONLY.csv" 
CSF_only_data <- read.csv(data_file_csf_only)

#Creation of interval variables
CSF_only_data$PLASMA_CDR_INT <- as.numeric(abs(difftime(CSF_only_data$EXAMDATE,CSF_only_data$VISDATE_CDR  ,units = "days"))/365.25)
CSF_only_data$PLASMA_TAU_INT <- as.numeric(abs(difftime(CSF_only_data$EXAMDATE,CSF_only_data$SCANDATE_TAU,units = "days"))/365.25)
CSF_only_data$PLASMA_AMY_INT <- as.numeric(abs(difftime(CSF_only_data$EXAMDATE,CSF_only_data$SCANDATE_AMY  ,units = "days"))/365.25)
CSF_only_data$PLASMA_MRI_INT <- as.numeric(abs(difftime(CSF_only_data$EXAMDATE,CSF_only_data$SCANDATE_ATROPHY  ,units = "days"))/365.25)
CSF_only_data$PLASMA_CSF_INT <- as.numeric(abs(difftime(CSF_only_data$EXAMDATE,CSF_only_data$EXAMDATE_CSF  ,units = "days"))/365.25)

##didn't make separate function for CDR Category, so needed to make this custom to copy into table
table(CSF_only_data$CDR_CAT)
table(CSF_only_data[which(CSF_only_data$CENTILOIDS_10==0),]$CDR_CAT)
table(CSF_only_data[which(CSF_only_data$CENTILOIDS_10==1),]$CDR_CAT)
chisq.test(CSF_only_data[,"CENTILOIDS_10"], 
           as.factor(CSF_only_data[,"CDR_CAT"]), 
           correct=FALSE)


####Function call - Function is in Functions_for_Analysis.R
char_table_full_bin(CSF_only_data,
                    "CENTILOIDS_10",
                    "Amyloid PET",
                    varb_1,
                    type_1,
                    label_1,
                    paste0(out_results,"Characteristics_Table_by_CENT20_05162024_CSFONLY.xlsx"))

####
## Tau PET complete cases - characteristics tables
####

TAU_only_data <- subset(last_crossectional,!is.na(MesialTemporal))

table(TAU_only_data$CDR_CAT)
table(TAU_only_data[which(TAU_only_data$CENTILOIDS_10==0),]$CDR_CAT)
table(TAU_only_data[which(TAU_only_data$CENTILOIDS_10==1),]$CDR_CAT)

chisq.test(TAU_only_data[,"CENTILOIDS_10"], 
           as.factor(TAU_only_data[,"CDR_CAT"]), 
           correct=FALSE)

char_table_full_bin(TAU_only_data,
                    "CENTILOIDS_10",
                    "Amyloid PET",
                    varb_1,
                    type_1,
                    label_1,
                    paste0(out_results,"Characteristics_Table_by_CENT20_05212024_TAUONLY.xlsx"))


####
## Atrophy complete cases - characteristics tables
####
ATROPHY_only_data <- subset(last_crossectional,!is.na(atrophy))
ATROPHY_only_data$LP_ATROPHY_INT <- as.numeric(abs(difftime(ATROPHY_only_data$EXAMDATE,ATROPHY_only_data$SCANDATE_ATROPHY,units = "days"))/365.25)

table(ATROPHY_only_data$CDR_CAT)
table(ATROPHY_only_data[which(ATROPHY_only_data$CENTILOIDS_10==0),]$CDR_CAT)
table(ATROPHY_only_data[which(ATROPHY_only_data$CENTILOIDS_10==1),]$CDR_CAT)


chisq.test(ATROPHY_only_data[,"CENTILOIDS_10"], 
           as.factor(ATROPHY_only_data[,"CDR_CAT"]), 
           correct=FALSE)


char_table_full_bin(ATROPHY_only_data,
                    "CENTILOIDS_10",
                    "Amyloid PET",
                    varb_1,
                    type_1,
                    label_1,
                    paste0(out_results,"Characteristics_Table_by_CENT20_05212024_ATROPHYONLY.xlsx"))