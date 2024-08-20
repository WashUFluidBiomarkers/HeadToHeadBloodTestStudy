#####Spearman Analyses
#### Coded: Benjamin Saef
#### Date as of commenting and documentation - 05/28/2024


library(tidyr)
library(data.table)
library(RVAideMemoire)
library(pROC)
library(RColorBrewer)
library(MuMIn)
source("./Functions_for_Analysis.R")



#NO YEARS OF EDUCATION
#Make partial comparisons to top comparing partial
#Do algorithmic check on top performer before setting reference - new function?
#Remove all Quanterix only individuals

####Set seed - bootstrapping is involved so make sure
set.seed(12345)




#keep this script in the same folder as the "Functions_for_Analysis.R"
out_project <- "/FNIH_Project_2/Paper1/Spearman_Analysis" #Folder where you want all output written to, a results and plot folder will be written here by default
out_results <- paste0(out_project,"/results/")
out_plots <- paste0(out_project,"/plots/")

dir.create(out_project,showWarnings = F)
dir.create(out_results,showWarnings = F)
dir.create(out_plots,showWarnings = F)

data_file <- "/FNIH_Project_2/DATA_STUDY_2/cross_sectional_data_2024_04_24b.csv"  # cross sectional dataset

###read in dataset
data.in <- read.csv(data_file)

####This was done so I could verify the analysis by running analysis code in a different programming language
####write.csv(data.in,"/FNIH_Project_2/DATA_STUDY_2/cross_sectional_data_2024_04_24b_FOR_SAS.csv",na = ".",row.names = F)
sex_var <- "PTGENDER" #variable with sex
age_var <- "AGE" #variable with age
APOE_variable <- "APOE_genotype" #Variable with APOE genotypes

#analysis vars
#C2N
pTau217_c2n <- "C2N_plasma_ptau217_ratio" #Variable representing p-tau217 ratio - c2n
pTau217conc_c2n <- "C2N_plasma_ptau217" #Variable representing p-tau217 concentration - c2n
abeta42_c2n <- "C2N_plasma_Abeta42" #Variable representing Ab42 - c2n
abeta40_c2n <- "C2N_plasma_Abeta42_Abeta40" #Variable representing Ab40 - c2n
abeta4240_c2n <- "C2N_plasma_Abeta42_Abeta40" #Variable representing Ab42/Ab40 - c2n

#Elecsys
abeta42_esys <- "Roche_plasma_Ab42" #Variable representing Ab42 - elecsys
abeta40_esys <- "Roche_plasma_Ab40" #Variable representing Ab40 - elecsys
abeta4240_esys <- "Roche_plasma_Ab42_Ab40" #Variable representing Ab42/Ab40 - elecsys
ptau_181_esys <- "Roche_plasma_ptau181" #Variable representing ptau181 - elecsys
GFAP_esys <- "Roche_plasma_GFAP" #Variable representing GFAP - elecsys
NfL_esys <- "Roche_plasma_NfL" #Variable representing NfL - elecsys


#Lumipulse

abeta42_lumi <- "Fuji_plasma_Ab42" #Variable representing Ab42 - lumipulse
abeta40_lumi <- "Fuji_plasma_Ab40" #Variable representing Ab40 - lumipulse
abeta4240_lumi <- "Fuji_plasma_Ab42_Ab40" #Variable representing Ab42/Ab40 - lumipulse
pTau217_lumi <- "Fuji_plasma_ptau217" #Variable representing p-tau217 ratio - lumipulse

#Simoa 
alzpath_ptau217 <- "AlzPath_plasma_ptau217"  #Variable representing p-tau217 - AlzPath
janssen_ptau217 <- "Janssen_plasma_ptau217"  #Variable representing p-tau217 - Janssen

#Quanterix - Simoa
abeta42_quan <- "QX_plasma_Ab42" #Variable representing Ab42 - elecsys
abeta40_quan <- "QX_plasma_Ab40" #Variable representing Ab40 - elecsys
abeta4240_quan <- "QX_plasma_Ab42_Ab40" #Variable representing Ab42/Ab40 - elecsys
ptau_181_quan <- "QX_plasma_ptau181" #Variable representing ptau181 - elecsys
GFAP_quan <- "QX_plasma_GFAP" #Variable representing GFAP - elecsys
NfL_quan <- "QX_plasma_NfL" #Variable representing NfL - elecsys

#Variabels for PET metrics
centiloid_Variable <- "CENTILOIDS" #Variable representing Centiloid
early_tau <- "MesialTemporal" #Variable representing Tau PET Centaur
late_tau <- "TemporoParietal"
atrophy <- "atrophy"
cdr_sum_of_boxes <- "CDR_SOB"
MMSE <- "MMSCORE"

#Making sure that the allele recordings fit neatly into the APOE levels we have selected
data.in[which(data.in[,APOE_variable]=="43"),APOE_variable] <- "34"
data.in[which(data.in[,APOE_variable]=="42"),APOE_variable] <- "24"

###Make sure APOE 33 is the reference level
data.in$APOE_Cat <- factor(data.in[,APOE_variable],
                              levels=c("33","34","24","22","23","44"))


#Covariates
c_var <- c(age_var, #Age
           sex_var, #Sex
           "APOE_Cat") #APOE category variable

#####Analyses groups per platform
#Groupings if necessary - these are repeated for the purpose of comparisons within group keeping measures from same company together
l_var <- c( "c2n", #Ab42/40
            "c2n", #ptau217
            "c2n",
            "roche", #Ab42/40
            "roche", #ptau181
            "roche", #GFAP
            "roche", #NFL
            "fuji",#Ab42/40
            "fuji",#ptau217
            "alz",#ptau217
            "jan",#ptau217
            "quan", #Ab42/40
            "quan", #ptau181
            "quan", #GFAP
            "quan"#NFL
)

##list out the different analytes - Labels
analytes <- c(paste0("A","\U03B2","42/A","\U03B2","40"), #Ab42/40
              "p-tau217 ratio (%)", #ptau217
              "p-tau217 (pg/ml)", #ptau217
              paste0("A","\U03B2","42/A","\U03B2","40"), #Ab42/40
              "p-tau181 (pg/ml)", #ptau181
              "GFAP (ng/ml)", #GFAP
              "NfL (pg/mL)", #NFL
              paste0("A","\U03B2","42/A","\U03B2","40"), #Ab42/40
              "p-tau217 (pg/ml)", #ptau217
              "p-tau217 (pg/ml)", #ptau217
              "p-tau217 (pg/ml)", #ptau217
              paste0("A","\U03B2","42/A","\U03B2","40"), #Ab42/40
              "p-tau181 (pg/ml)", #ptau181
              "GFAP (pg/ml)", #GFAP
              "NfL (pg/mL)" #NFL
)
#list out the variables representing those analytes.
i_var <- c(abeta4240_c2n, #Ab42/40 - c2n
           pTau217_c2n, #ptau217 - c2n
           pTau217conc_c2n, #ptau217 - c2n
           abeta4240_esys,#Ab42/40 - elecsys
           ptau_181_esys,#ptau181 - elecsys
           GFAP_esys,#GFAP - elecsys
           NfL_esys, #NFL - elecsys
           abeta4240_lumi,#Ab42/40 - lumipulse
           pTau217_lumi, #ptau217 - lumipulse
           alzpath_ptau217,#ptau217 - alzpath
           janssen_ptau217,#ptau217 - lumipulse
           abeta4240_quan,#Ab42/40 - quanterix
           ptau_181_quan,#ptau181 - quanterix
           GFAP_quan,#GFAP - quanterix
           NfL_quan #NFL - quanterix
)

###output variable vector creation
o_vars <- c(centiloid_Variable, # Centiloid
            early_tau,#Early Tau
            late_tau, #Late Tau
            atrophy, #Atrophy
            cdr_sum_of_boxes, #CDR sum of boxes
            MMSE) #MMSE (not actually used)


spdiff_Partial_rho <- function(data, indices) { ###sort of a double check on the function creation for the bootstratpping comparisons
  d <- data[indices,] # allows boot to select sample
  c_cov <- c(age_var,sex_var,"APOE_Cat")   #change to reflect data
  cor_xz <- RVAideMemoire::pcor.test(d$x, d$z,d[,c_cov], method = "spearman",nrep = 1)
  cor_yz <- RVAideMemoire::pcor.test(d$y, d$z, d[,c_cov],method = "spearman",nrep = 1)
  diff <- abs(cor_xz$estimate) - abs(cor_yz$estimate)
  
  return(diff)
}
###the above function also exists in Functions_for_Analysis.R


#####Centiloid

###Full Cohort

###Determine the comparison standard before running
comparison_standard <- Eval_comp_standard(data.in,o_vars[1],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(data.in,o_vars[1],c_var,i_var)
ptm <- proc.time()
table_out_centiloid_full <-  spearman_plasma_runs(data.in, # Data inputted 
                        o_vars[1],  # variable A - continuous variable being correlated
                        i_var,# vector of strings - continuous variables being correlated with variable A
                        comparison_standard,
                        comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                        c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                        analytes,
                        l_var)

#####Centiloid 20 - positive
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==1),o_vars[1],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==1),o_vars[1],c_var,i_var)
table_out_centiloid_amy20pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==1), # Data inputted 
                                                      o_vars[1],  # variable A - continuous variable being correlated
                                                      i_var,# vector of strings - continuous variables being correlated with variable A
                                                      comparison_standard,
                                                      comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                      c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                      analytes,
                                                      l_var)


####Centiloid 20 - negative
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==0),o_vars[1],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==0),o_vars[1],c_var,i_var)
ptm <- proc.time()
table_out_centiloid_amy20neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==0), # Data inputted 
                                                      o_vars[1],  # variable A - continuous variable being correlated
                                                      i_var,# vector of strings - continuous variables being correlated with variable A
                                                      comparison_standard,
                                                      comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                      c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                      analytes,
                                                      l_var)
proc.time() - ptm

####Centiloid 37 - positive
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==1),o_vars[1],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==1),o_vars[1],c_var,i_var)
table_out_centiloid_amy37pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==1), # Data inputted 
                                                      o_vars[1],  # variable A - continuous variable being correlated
                                                      i_var,# vector of strings - continuous variables being correlated with variable A
                                                      comparison_standard,
                                                      comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                      c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                      analytes,
                                                      l_var)


####Centiloid 37 negative
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==0),o_vars[1],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==0),o_vars[1],c_var,i_var)
table_out_centiloid_amy37neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==0), # Data inputted 
                                                      o_vars[1],  # variable A - continuous variable being correlated
                                                      i_var,# vector of strings - continuous variables being correlated with variable A
                                                      comparison_standard,
                                                      comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                      c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                      analytes,
                                                      l_var)

####CDR - positive
comparison_standard <- Eval_comp_standard(subset(data.in,CDR_10==1),o_vars[1],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CDR_10==1),o_vars[1],c_var,i_var)
table_out_centiloid_cdrpos <-  spearman_plasma_runs(subset(data.in,CDR_10==1), # Data inputted 
                                                   o_vars[1],  # variable A - continuous variable being correlated
                                                   i_var,# vector of strings - continuous variables being correlated with variable A
                                                   comparison_standard,
                                                   comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                   c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                   analytes,
                                                   l_var)

####CDR - negative
comparison_standard <- Eval_comp_standard(subset(data.in,CDR_10==0),o_vars[1],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CDR_10==0),o_vars[1],c_var,i_var)
table_out_centiloid_cdrneg <-  spearman_plasma_runs(subset(data.in,CDR_10==0), # Data inputted 
                                                    o_vars[1],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                    comparison_standard,
                                                    comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                    c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)


##### This function is only in this code, it loses use once the above 'analytes' vector is corrected.
##### putting it here for documentation's purpose
relabel_spearman_res <- function(DF){
  


  DF$Biomarker <-recode(DF$Biomarker, "p-Tau 217 Ratio" = "p-tau217 ratio (%)", 
                        "p-Tau 217 concentration" = "p-tau217 (pg/ml)",
                        "Neurofilament Light"="NfL (pg/mL)",
                        "GFAP"="GFAP (pg/ml)",
                        "p-Tau 181"="p-tau181 (pg/ml)")

  DF[which(DF$group=="roche" & DF$Biomarker=="GFAP (pg/ml)"),"Biomarker"] <- "GFAP (ng/ml)"
  DF[which(DF$Biomarker==paste0("Plasma A","\U03B2","42/A","\U03B2","40")),"Biomarker"] <- paste0("A","\U03B2","42/A","\U03B2","40")
  
  DF$company <- factor(recode(DF$group,"quan"="Quanterix","alz"="ALZpath","c2n"="C2N",
                              "fuji"="Fujirebio", "jan"="Janssen","roche" = "Roche"))
  levels(DF$company) <- c("ALZpath","C2N", "Fujirebio",  "Janssen", "Quanterix","Roche")
  DF$company <-factor(DF$company, levels = c("C2N", "Fujirebio", "ALZpath", "Janssen", "Roche", "Quanterix"))
  
return(DF)

}


#####compile all the results data into a single list to save out and use in further output
list_of_cent_correlation_dataset <- list("Full Data"=table_out_centiloid_full,
                                         "Cent20 Positive Data"=table_out_centiloid_amy20pos,
                                         "Cent20 Negative Data"=table_out_centiloid_amy20neg,
                                         "Cent37 Positive Data"=table_out_centiloid_amy37pos,
                                         "Cent37 Negative Data"=table_out_centiloid_amy37neg,
                                         "CDR Positive Data"=table_out_centiloid_cdrpos,
                                         "CDR Negative Data"=table_out_centiloid_cdrneg)
####Due to length of table creation - save output data upon finishing 
save(list_of_cent_correlation_dataset,
     file=paste0(out_results,"cent_correlations.rdata"))



########## deprecated - no longer necessary
#
#for(i in 1:length(list_of_cent_correlation_dataset)){
#  
#  
#  list_cent_new[[i]] <-  relabel_spearman_res(list_of_cent_correlation_dataset[[i]])
#  
#}
#names(list_cent_new) <- names(list_of_cent_correlation_dataset)
##########
###Create Clean Plots
list_cent_new <- list()
list_of_cent_clean <- list()
for(i_clean in 1: length(list_of_cent_correlation_dataset)){
  DF <- list_of_cent_correlation_dataset[[i_clean]]
  
  DF$company <- factor(recode(DF$group,"quan"="Quanterix","alz"="ALZpath","c2n"="C2N",
                              "fuji"="Fujirebio", "jan"="Janssen","roche" = "Roche")) ### rename group to company, set as factor
  levels(DF$company) <- c("ALZpath","C2N", "Fujirebio",  "Janssen", "Quanterix","Roche") #set factor levels properly
  DF$company <-factor(DF$company, levels = c("C2N", "Fujirebio", "ALZpath", "Janssen", "Roche", "Quanterix")) # order factor properly
  
  list_cent_new[[i_clean]] <- DF ###add dataset with proper variables to go into plot creation/cleaning to new list
  list_of_cent_clean[[i_clean]]  <- clean_spearman_table(DF)
}

names(list_of_cent_clean) <- names(list_of_cent_correlation_dataset)

###Write out publication format tables
write.xlsx(list_of_cent_clean,
           file=paste0(out_results,"Centiloid_Correlations_05132024.xlsx"))

#####Tau Early

###Full cohort

###determine comparison standard for cohort, then run
comparison_standard <- Eval_comp_standard(data.in,o_vars[2],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(data.in,o_vars[2],c_var,i_var)
ptm <- proc.time()
table_out_tauearly_full <-  spearman_plasma_runs(data.in, # Data inputted 
                                                  o_vars[2],  # variable A - continuous variable being correlated
                                                  i_var,# vector of strings - continuous variables being correlated with variable A
                                                  comparison_standard,
                                                 comparison_standard_partial,#variable in vars_to_cor vector that is used as standard for comparisons
                                                  c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                  analytes,
                                                  l_var)
proc.time() - ptm

#### Centiloid 20 positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==1),o_vars[2],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==1),o_vars[2],c_var,i_var)
table_out_tauearly_amy20pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==1), # Data inputted 
                                                      o_vars[2],  # variable A - continuous variable being correlated
                                                      i_var,# vector of strings - continuous variables being correlated with variable A
                                                     comparison_standard,
                                                     comparison_standard_partial,#variable in vars_to_cor vector that is used as standard for comparisons
                                                     c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                      analytes,
                                                      l_var)

proc.time() - ptm

#### Centiloid 20 negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==0),o_vars[2],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==0),o_vars[2],c_var,i_var)
table_out_tauearly_amy20neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==0), # Data inputted 
                                                      o_vars[2],  # variable A - continuous variable being correlated
                                                      i_var,# vector of strings - continuous variables being correlated with variable A
                                                     comparison_standard,
                                                     comparison_standard_partial,#variable in vars_to_cor vector that is used as standard for comparisons
                                                     c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                      analytes,
                                                      l_var)

proc.time() - ptm

#### Centiloid 37 positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==1),o_vars[2],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==1),o_vars[2],c_var,i_var)
table_out_tauearly_amy37pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==1), # Data inputted 
                                                      o_vars[2],  # variable A - continuous variable being correlated
                                                      i_var,# vector of strings - continuous variables being correlated with variable A
                                                     comparison_standard,
                                                     comparison_standard_partial,#variable in vars_to_cor vector that is used as standard for comparisons
                                                     c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                      analytes,
                                                      l_var)
proc.time() - ptm

#### Centiloid 37 negative

ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==0),o_vars[2],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==0),o_vars[2],c_var,i_var)
table_out_tauearly_amy37neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==0), # Data inputted 
                                                      o_vars[2],  # variable A - continuous variable being correlated
                                                      i_var,# vector of strings - continuous variables being correlated with variable A
                                                     comparison_standard,
                                                     comparison_standard_partial,#variable in vars_to_cor vector that is used as standard for comparisons
                                                     c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                      analytes,
                                                      l_var)

proc.time() - ptm

#### CDR positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CDR_10==1),o_vars[2],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CDR_10==1),o_vars[2],c_var,i_var)
table_out_tauearly_cdrpos <-  spearman_plasma_runs(subset(data.in,CDR_10==1), # Data inputted 
                                                    o_vars[2],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                   comparison_standard,
                                                   comparison_standard_partial,#variable in vars_to_cor vector that is used as standard for comparisons
                                                   c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)
proc.time() - ptm

#### CDR negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CDR_10==0),o_vars[2],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CDR_10==0),o_vars[2],c_var,i_var)
table_out_tauearly_cdrneg <-  spearman_plasma_runs(subset(data.in,CDR_10==0), # Data inputted 
                                                    o_vars[2],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                   comparison_standard,
                                                   comparison_standard_partial,#variable in vars_to_cor vector that is used as standard for comparisons
                                                   c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)
proc.time() - ptm

####Add all results tables to list
list_of_tauearly_correlation_dataset <- list("Full Data"=table_out_tauearly_full,
                                         "Cent20 Positive Data"=table_out_tauearly_amy20pos,
                                         "Cent20 Negative Data"=table_out_tauearly_amy20neg,
                                         "Cent37 Positive Data"=table_out_tauearly_amy37pos,
                                         "Cent37 Negative Data"=table_out_tauearly_amy37neg,
                                         "CDR Positive Data"=table_out_tauearly_cdrpos,
                                         "CDR Negative Data"=table_out_tauearly_cdrneg)
####Due to length of table creation - save output data upon finishing 
save(list_of_tauearly_correlation_dataset,
     file=paste0(out_results,"tauearly_correlations.rdata"))


####Deprecated - no longer useful after fixing analytes vector
#list_tauearly_new <- list()
#for(i in 1:length(list_of_tauearly_correlation_dataset)){
#  
#  
#  list_tauearly_new[[i]] <-  relabel_spearman_res(list_of_tauearly_correlation_dataset[[i]])
#  
#}
#
#names(list_tauearly_new) <- names(list_of_tauearly_correlation_dataset)

write.xlsx(list_of_tauearly_correlation_dataset,
           file=paste0(out_results,"TauEarly_Correlations.xlsx"))

###Create Clean Plots
list_tauearly_new <- list()
list_of_tauearly_clean <- list()
for(i_clean in 1: length(list_of_tauearly_correlation_dataset)){
  DF <- list_of_tauearly_correlation_dataset[[i_clean]]
  
  DF$company <- factor(recode(DF$group,"quan"="Quanterix","alz"="ALZpath","c2n"="C2N",
                              "fuji"="Fujirebio", "jan"="Janssen","roche" = "Roche")) ### rename group to company, set as factor
  levels(DF$company) <- c("ALZpath","C2N", "Fujirebio",  "Janssen", "Quanterix","Roche") #set factor levels properly
  DF$company <-factor(DF$company, levels = c("C2N", "Fujirebio", "ALZpath", "Janssen", "Roche", "Quanterix")) # order factor properly
  
  list_tauearly_new[[i_clean]] <- DF ###add dataset with proper variables to go into plot creation/cleaning to new list
  list_of_tauearly_clean[[i_clean]]  <- clean_spearman_table(DF)
}

names(list_of_tauearly_clean) <- names(list_of_tauearly_correlation_dataset)

#write out publication format tables
write.xlsx(list_of_tauearly_clean,
           file=paste0(out_results,"TauEarly_Correlations_05112024.xlsx"))

#####


####TAU LATE

####Full Data
#### Determine comparison standards for each cohort then run
comparison_standard <- Eval_comp_standard(data.in,o_vars[3],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(data.in,o_vars[3],c_var,i_var)

ptm <- proc.time()
table_out_TauLate_full <-  spearman_plasma_runs(data.in, # Data inputted 
                                                 o_vars[3],  # variable A - continuous variable being correlated
                                                 i_var,# vector of strings - continuous variables being correlated with variable A
                                                 comparison_standard,
                                                comparison_standard_partial,#variable in vars_to_cor vector that is used as standard for comparisons
                                                 c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                 analytes,
                                                 l_var)
proc.time() - ptm


##### Amyloid 20 positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==1),o_vars[3],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==1),o_vars[3],c_var,i_var)

table_out_TauLate_amy20pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==1), # Data inputted 
                                                     o_vars[3],  # variable A - continuous variable being correlated
                                                     i_var,# vector of strings - continuous variables being correlated with variable A
                                                    comparison_standard,
                                                    comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                     c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                     analytes,
                                                     l_var)

proc.time() - ptm


##### Amyloid 20 negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==0),o_vars[3],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==0),o_vars[3],c_var,i_var)

table_out_TauLate_amy20neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==0), # Data inputted 
                                                     o_vars[3],  # variable A - continuous variable being correlated
                                                     i_var,# vector of strings - continuous variables being correlated with variable A
                                                    comparison_standard,
                                                    var_comp_p = comparison_standard_partial,
                                                     c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                     analytes,
                                                     l_var)
proc.time() - ptm

##### Amyloid 37 positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==1),o_vars[3],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==1),o_vars[3],c_var,i_var)

table_out_TauLate_amy37pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==1), # Data inputted 
                                                     o_vars[3],  # variable A - continuous variable being correlated
                                                     i_var,# vector of strings - continuous variables being correlated with variable A
                                                    comparison_standard,
                                                    comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                     c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                     analytes,
                                                     l_var)
proc.time() - ptm

##### Amyloid 37 negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==0),o_vars[3],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==0),o_vars[3],c_var,i_var)

table_out_TauLate_amy37neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==0), # Data inputted 
                                                     o_vars[3],  # variable A - continuous variable being correlated
                                                     i_var,# vector of strings - continuous variables being correlated with variable A
                                                    comparison_standard,
                                                    comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                     c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                     analytes,
                                                     l_var)
proc.time() - ptm

##### CDR positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CDR_10==1),o_vars[3],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CDR_10==1),o_vars[3],c_var,i_var)

table_out_TauLate_cdrpos <-  spearman_plasma_runs(subset(data.in,CDR_10==1), # Data inputted 
                                                   o_vars[3],  # variable A - continuous variable being correlated
                                                   i_var,# vector of strings - continuous variables being correlated with variable A
                                                  comparison_standard,
                                                  comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                   c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                   analytes,
                                                   l_var)
proc.time() - ptm

##### CDR negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CDR_10==0),o_vars[3],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CDR_10==0),o_vars[3],c_var,i_var)

table_out_TauLate_cdrneg <-  spearman_plasma_runs(subset(data.in,CDR_10==0), # Data inputted 
                                                   o_vars[3],  # variable A - continuous variable being correlated
                                                   i_var,# vector of strings - continuous variables being correlated with variable A
                                                   comparison_standard,
                                                  comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                   c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                   analytes,
                                                   l_var)

proc.time() - ptm


#####Deprecated as analytes vector fixed
#list_taulate_new <- list()
#for(i in 1:length(list_of_TauLate_correlation_dataset)){
#  
#  
#  list_taulate_new[[i]] <-  relabel_spearman_res(list_of_TauLate_correlation_dataset[[i]])
#  
#}
#
#names(list_taulate_new) <- names(list_of_TauLate_correlation_dataset)

list_of_TauLate_correlation_dataset <- list("Full Data"=table_out_TauLate_full,
                                             "Cent20 Positive Data"=table_out_TauLate_amy20pos,
                                             "Cent20 Negative Data"=table_out_TauLate_amy20neg,
                                             "Cent37 Positive Data"=table_out_TauLate_amy37pos,
                                             "Cent37 Negative Data"=table_out_TauLate_amy37neg,
                                             "CDR Positive Data"=table_out_TauLate_cdrpos,
                                             "CDR Negative Data"=table_out_TauLate_cdrneg)
####Due to length of table creation - save output data upon finishing 
save(list_of_TauLate_correlation_dataset,
     file=paste0(out_results,"taulate_correlations.rdata"))
write.xlsx(list_of_TauLate_correlation_dataset,
           file=paste0(out_results,"TauLate_Correlations.xlsx"))




###Create Clean Plots
list_taulate_new <- list()
list_of_taulate_clean <- list()
for(i_clean in 1: length(list_of_TauLate_correlation_dataset)){
  DF <- list_of_TauLate_correlation_dataset[[i_clean]]
  
  DF$company <- factor(recode(DF$group,"quan"="Quanterix","alz"="ALZpath","c2n"="C2N",
                              "fuji"="Fujirebio", "jan"="Janssen","roche" = "Roche"))
  levels(DF$company) <- c("ALZpath","C2N", "Fujirebio",  "Janssen", "Quanterix","Roche")
  DF$company <-factor(DF$company, levels = c("C2N", "Fujirebio", "ALZpath", "Janssen", "Roche", "Quanterix"))
  
  list_taulate_new[[i_clean]] <- DF
  list_of_taulate_clean[[i_clean]]  <- clean_spearman_table(DF)
}

names(list_of_taulate_clean) <- names(list_of_TauLate_correlation_dataset)
write.xlsx(list_of_taulate_clean,
           file=paste0(out_results,"TauLate_Correlations_05112024.xlsx"))





#####Full data
####Select comparison standards for each cohort - then run
comparison_standard <- Eval_comp_standard(data.in,o_vars[4],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(data.in,o_vars[4],c_var,i_var)
ptm <- proc.time()
table_out_Atrophy_full <-  spearman_plasma_runs(data.in, # Data inputted 
                                                o_vars[4],  # variable A - continuous variable being correlated
                                                i_var,# vector of strings - continuous variables being correlated with variable A
                                                comparison_standard,
                                                comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                analytes,
                                                l_var)
proc.time() - ptm

#####Centiloid 20 positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==1),o_vars[4],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==1),o_vars[4],c_var,i_var)

table_out_Atrophy_amy20pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==1), # Data inputted 
                                                    o_vars[4],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                    comparison_standard,
                                                    comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                    c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)
proc.time() - ptm

#### Centiloid 20 negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==0),o_vars[4],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==0),o_vars[4],c_var,i_var)
table_out_Atrophy_amy20neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==0), # Data inputted 
                                                    o_vars[4],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                    comparison_standard,
                                                    comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                    c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)

proc.time() - ptm

##### Centiloid 37 Positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==1),o_vars[4],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==1),o_vars[4],c_var,i_var)

table_out_Atrophy_amy37pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==1), # Data inputted 
                                                    o_vars[4],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                    comparison_standard,
                                                    comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                    c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)
proc.time() - ptm

##### Centiloid 37 negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==0),o_vars[4],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==0),o_vars[4],c_var,i_var)

table_out_Atrophy_amy37neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==0), # Data inputted 
                                                    o_vars[4],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                    comparison_standard,
                                                    comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                    c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)
proc.time() - ptm

##### CDR positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CDR_10==1),o_vars[4],i_var)

comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CDR_10==1),o_vars[4],c_var,i_var)

table_out_Atrophy_cdrpos <-  spearman_plasma_runs(subset(data.in,CDR_10==1), # Data inputted 
                                                  o_vars[4],  # variable A - continuous variable being correlated
                                                  i_var,# vector of strings - continuous variables being correlated with variable A
                                                  comparison_standard,
                                                  comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                  c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                  analytes,
                                                  l_var)
proc.time() - ptm

##### CDR negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CDR_10==0),o_vars[4],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CDR_10==0),o_vars[4],c_var,i_var)

table_out_Atrophy_cdrneg <-  spearman_plasma_runs(subset(data.in,CDR_10==0), # Data inputted 
                                                  o_vars[4],  # variable A - continuous variable being correlated
                                                  i_var,# vector of strings - continuous variables being correlated with variable A
                                                  comparison_standard,
                                                  comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                  c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                  analytes,
                                                  l_var)

#### Assemble outputs into list
list_of_Atrophy_correlation_dataset <- list("Full Data"=table_out_Atrophy_full,
                                            "Cent20 Positive Data"=table_out_Atrophy_amy20pos,
                                            "Cent20 Negative Data"=table_out_Atrophy_amy20neg,
                                            "Cent37 Positive Data"=table_out_Atrophy_amy37pos,
                                            "Cent37 Negative Data"=table_out_Atrophy_amy37neg,
                                            "CDR Positive Data"=table_out_Atrophy_cdrpos,
                                            "CDR Negative Data"=table_out_Atrophy_cdrneg)



####Due to length of table creation - save output data upon finishing 
save(list_of_Atrophy_correlation_dataset,
     file=paste0(out_results,"atrophy_correlations.rdata"))

#### Clean up only necessary before analytes vector correction
#list_atrophy_new <- list()
#for(i in 1:length(list_of_Atrophy_correlation_dataset)){
#  
#  
#  list_atrophy_new[[i]] <-  relabel_spearman_res(list_of_Atrophy_correlation_dataset[[i]])
#  
#}
#
#names(list_atrophy_new) <- names(list_of_Atrophy_correlation_dataset)
#####

write.xlsx(list_of_Atrophy_correlation_dataset,
           file=paste0(out_results,"Atrophy_Correlations.xlsx"))


###Create Clean Plots
list_atrophy_new <- list()
list_of_atrophy_clean <- list()
for(i_clean in 1: length(list_of_Atrophy_correlation_dataset)){
  DF <- list_of_Atrophy_correlation_dataset[[i_clean]]
  
  DF$company <- factor(recode(DF$group,"quan"="Quanterix","alz"="ALZpath","c2n"="C2N",
                              "fuji"="Fujirebio", "jan"="Janssen","roche" = "Roche"))
  levels(DF$company) <- c("ALZpath","C2N", "Fujirebio",  "Janssen", "Quanterix","Roche")
  DF$company <-factor(DF$company, levels = c("C2N", "Fujirebio", "ALZpath", "Janssen", "Roche", "Quanterix"))
  list_atrophy_new[[i_clean]] <- DF
  
  list_of_atrophy_clean[[i_clean]]  <- clean_spearman_table(DF)
}

names(list_of_atrophy_clean) <- names(list_of_Atrophy_correlation_dataset)
write.xlsx(list_of_atrophy_clean,
           file=paste0(out_results,"Atrophy_Correlations_05132024.xlsx"))

#### Ended up not being used for final paper, here for documationation purposesd


#Atrophyplots <- list()
#for(i_plot in 1:length(list_atrophy_new)){
#  Atrophyplots[[i_plot]] <- spearman_gg_formanuscript(list_atrophy_new[[i_plot]],breaks1 = plot_breaks[[4]],indicatorline = -0.25,rev_axis = TRUE)
#}

#for(i_plot in 1:length(list_atrophy_new)){
#  ggsave(paste0(out_plots,names(list_atrophy_new)[i_plot],"_Atrophy_Spearman_Forest.pdf"),
#         plot = Atrophyplots[[i_plot]],
#         device=cairo_pdf)#writes out spearman plot
#  
#}




###CDR Sumbox

####Create comparison standards for unadjusted and adjusted spearmans
comparison_standard <- Eval_comp_standard(data.in,o_vars[5],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(data.in,o_vars[5],c_var,i_var)

ptm <- proc.time()

####Run spearman table

### Full Data
table_out_sumbox_full <-  spearman_plasma_runs(data.in, # Data inputted 
                                                o_vars[5],  # variable A - continuous variable being correlated
                                                i_var,# vector of strings - continuous variables being correlated with variable A
                                               comparison_standard,
                                               comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                analytes,
                                                l_var)

proc.time() - ptm

###### Centiloid 20 positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==1),o_vars[5],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==1),o_vars[5],c_var,i_var)

table_out_sumbox_amy20pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==1), # Data inputted 
                                                    o_vars[5],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                   comparison_standard,
                                                   comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                    c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)
proc.time() - ptm

####Centiloid 20 Negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10==0),o_vars[5],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10==0),o_vars[5],c_var,i_var)
table_out_sumbox_amy20neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10==0), # Data inputted 
                                                    o_vars[5],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                   comparison_standard,
                                                   comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                    c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)
proc.time() - ptm

#### Centiloid 37 Positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==1),o_vars[5],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==1),o_vars[5],c_var,i_var)
table_out_sumbox_amy37pos <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==1), # Data inputted 
                                                    o_vars[5],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                   comparison_standard,
                                                   comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                    c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)
proc.time() - ptm

####Centiloid 37 Negative
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CENTILOIDS_10_sens==0),o_vars[5],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CENTILOIDS_10_sens==0),o_vars[5],c_var,i_var)
table_out_sumbox_amy37neg <-  spearman_plasma_runs(subset(data.in,CENTILOIDS_10_sens==0), # Data inputted 
                                                    o_vars[5],  # variable A - continuous variable being correlated
                                                    i_var,# vector of strings - continuous variables being correlated with variable A
                                                   comparison_standard,
                                                   comparison_standard_partial, #variable in vars_to_cor vector that is used as standard for comparisons
                                                    c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                    analytes,
                                                    l_var)
proc.time() - ptm

#####CDR positive
ptm <- proc.time()
comparison_standard <- Eval_comp_standard(subset(data.in,CDR_10==1),o_vars[5],i_var)
comparison_standard_partial <- Eval_comp_standard_partial(subset(data.in,CDR_10==1),o_vars[5],c_var,i_var)

table_out_sumbox_cdrpos <-  spearman_plasma_runs(subset(data.in,CDR_10==1), # Data inputted 
                                                  o_vars[5],  # variable A - continuous variable being correlated
                                                  i_var,# vector of strings - continuous variables being correlated with variable A
                                                  comparison_standard,
                                                 comparison_standard_partial,#variable in vars_to_cor vector that is used as standard for comparisons
                                                  c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                                  analytes,
                                                  l_var)
proc.time() - ptm
#####CDR negative - unnecessary for sum of boxes correlations
#table_out_sumbox_cdrneg <-  spearman_plasma_runs(subset(data.in,CDR_10==0), # Data inputted 
#                                                  o_vars[5],  # variable A - continuous variable being correlated
#                                                  i_var,# vector of strings - continuous variables being correlated with variable A
#                                                  comparison_standard, #variable in vars_to_cor vector that is used as standard for comparisons
#                                                 c_var,    # covariates - variables used as covariates (age, sex, APOE normally)
#                                                  analytes,
#                                                  l_var)


#####Create list of results datasets
list_of_sumbox_correlation_dataset <- list("Full Data"=table_out_sumbox_full,
                                            "Cent20 Positive Data"=table_out_sumbox_amy20pos,
                                            "Cent20 Negative Data"=table_out_sumbox_amy20neg,
                                            "Cent37 Positive Data"=table_out_sumbox_amy37pos,
                                            "Cent37 Negative Data"=table_out_sumbox_amy37neg,
                                            "CDR Positive Data"=table_out_sumbox_cdrpos)

write.xlsx(list_of_sumbox_correlation_dataset,
           file=paste0(out_results,"sumbox_Correlations.xlsx"))

####Due to length of table creation - save output data upon finishing 
save(list_of_sumbox_correlation_dataset,
     file=paste0(out_results,"cdrsumbox_correlations.rdata"))

###Create Clean tables - These will create publication quality tables
list_of_sumbox_clean <- list()
list_cdr_new <- list()
for(i_clean in 1: length(list_of_sumbox_correlation_dataset)){
  DF <- list_of_sumbox_correlation_dataset[[i_clean]] ###Asign Dataset to work on
  
  DF$company <- factor(recode(DF$group,"quan"="Quanterix","alz"="ALZpath","c2n"="C2N",
                              "fuji"="Fujirebio", "jan"="Janssen","roche" = "Roche")) ###Recode group to company
  levels(DF$company) <- c("ALZpath","C2N", "Fujirebio",  "Janssen", "Quanterix","Roche")
  DF$company <-factor(DF$company, levels = c("C2N", "Fujirebio", "ALZpath", "Janssen", "Roche", "Quanterix")) #Order company factor
  list_cdr_new[[i_clean]] <- DF # this changed after the analytes vector was corrected
  
  list_of_sumbox_clean[[i_clean]]  <- clean_spearman_table(DF) ###Call the cleaning function
}

names(list_of_sumbox_clean) <- names(list_of_sumbox_correlation_dataset) #make sure all the dataset labels are correct

###Write out clean - publication formatted datasets
write.xlsx(list_of_sumbox_clean,
           file=paste0(out_results,"CDR_Sumbox_Correlations_05132024.xlsx"))

####

cdr_index <- index_biomarker_forest(table_out_cdr, #dataset to form basis for index
                       which(colnames(table_out_cdr)=="rho"), #column number for the metric being used - basis
                       which(colnames(table_out_cdr)=="group"),#column number for group variable - basis
                       "C2N", #way to distinguish groups within groups (basically, given this is nested, this tells which to highlight)
                       which(colnames(table_out_cdr)=="Biomarker"), #Biomarker identification variable where the function should search for a way to distinguish
                      l_var[-which(i_var==comparison_standard)])


cdr_plot <- Forest_plot_create(table_out_cdr, #creates forest plots for spearmans
                               cdr_index,
                               "CDR Sum of Boxes Correlations",
                               "Spearman Rho",
                               reverse_scales=FALSE)

##########These were more necessary when the 
#list_cdr_new <- list()
#for(i in 1:length(list_of_sumbox_correlation_dataset)){
#  
#  
#  list_cdr_new[[i]] <-  relabel_spearman_res(list_of_sumbox_correlation_dataset[[i]])
  
#}

#names(list_cdr_new) <- names(list_of_sumbox_correlation_dataset)

save.image("/FNIH_Project_2/Current_Analysis_data.rdata")
###########################
#### Plots ##### All Functions can be found in Functions_for_Analysis.R
###########################

plot_breaks <- list() ##create list of breaks for each plot
plot_breaks[[1]] <- c(0.20, 0.40, 0.60, 0.80)
plot_breaks[[2]] <- c(0, 0.2, 0.4, 0.6, 0.8)
plot_breaks[[3]] <- c(0, 0.2, 0.4, 0.6, 0.8)
plot_breaks[[4]] <- c(0.2, 0,-0.2, -0.40)
plot_breaks[[5]] <- c(0.1, 0.3, 0.5)
names(plot_breaks) <- c("Centiloid","TauEarly","TauLate","Atrophy","CDR_SB")

refline_curr <- c(0.6,0.4,-0.2,0.3) ##create vector of positions for reference line (line is arbitrary for visulization)
## we did try several solutions to make a line that would have some meaningful purpose within the dataset but it ultimately didn't
## work out with our range of values

####Full cohort plots
cent_index <- index_biomarker_forest(list_cent_new[[1]]) ### input results dataset

##### Centiloid
centiloid_plot_all <- spearman_gg_formanuscript(list_cent_new[[1]], ### Results Date
                                                breaks1 = plot_breaks$Centiloid, #### Breaks
                                                df_plot_index = cent_index, #### Index from Biomarker Forest index function
                                                indicatorline = refline_curr[1], ###Reference line for plot
                                                labels_all = TRUE, ####Indicates whether or not the labels will show up on this plot
                                                title_plot = "Amyloid PET", #### Title of plot
                                                amytag=TRUE) #### Specific paramenter to set axis
centiloid_plot_all <- centiloid_plot_all+ggtitle("Amyloid PET")


#### Tau Early
tauearly_plot_all <- spearman_gg_formanuscript(list_tauearly_new[[1]],breaks1 = plot_breaks$TauEarly,
                                               df_plot_index = cent_index,indicatorline = refline_curr[2],labels_all = FALSE,title_plot = "Early tau PET")
tauearly_plot_all <- tauearly_plot_all+ggtitle("Early tau PET")


#### Atrophy
atrophy_plot_all <- spearman_gg_formanuscript(list_atrophy_new[[1]],breaks1 = plot_breaks$Atrophy,
                                              df_plot_index = cent_index,indicatorline = refline_curr[3],rev_axis = TRUE,labels_all = FALSE,
                                              title_plot = "Cortical thickness")
atrophy_plot_all <- atrophy_plot_all+ggtitle("Cortical thickness")


##### CDR Sum of Boxes
cdr_plot_all <- spearman_gg_formanuscript(list_cdr_new[[1]],breaks1 = plot_breaks$CDR_SB,
                                          df_plot_index = cent_index,indicatorline = refline_curr[4],labels_all = FALSE,title_plot = "Dementia severity")

cdr_plot_all <- cdr_plot_all+ggtitle("Dementia severity")
plot_out <- centiloid_plot_all | tauearly_plot_all | atrophy_plot_all | cdr_plot_all

ggsave(paste0(out_plots,"Full_cohort_ATN.jpg"),
       plot=plot_out,
       device=jpeg,width = 8,height=4.02)
####Amyloid Positive cohort plots
cent_index <- index_biomarker_forest(list_cent_new[[2]])

##### Centiloid
centiloid_plot_amypos <- spearman_gg_formanuscript(list_cent_new[[2]],breaks1 = c(plot_breaks$Centiloid,0,-0.2),
                                                   df_plot_index = cent_index,indicatorline = 0.6,labels_all = TRUE,title_plot = "Amyloid PET")
centiloid_plot_amypos <- centiloid_plot_amypos+ggtitle("Amyloid PET")

##### Tau Early
tauearly_plot_amypos <- spearman_gg_formanuscript(list_tauearly_new[[2]],breaks1 = c(plot_breaks$TauEarly,-0.2),
                                                  df_plot_index = cent_index,indicatorline = 0.4,labels_all = FALSE,title_plot = "Early tau PET")
tauearly_plot_amypos <- tauearly_plot_amypos+ggtitle("Early tau PET")

##### Atrophy
atrophy_plot_amypos <- spearman_gg_formanuscript(list_atrophy_new[[2]],breaks1 = c(plot_breaks$Atrophy,-0.6),
                                                 df_plot_index = cent_index,indicatorline = -0.2,rev_axis = TRUE,labels_all = FALSE,
                                                 title_plot = "Cortical thickness")
atrophy_plot_amypos <- atrophy_plot_amypos+ggtitle("Cortical thickness")


##### CDR Sum of boxes
cdr_plot_amypos <- spearman_gg_formanuscript(list_cdr_new[[2]],breaks1 = c(plot_breaks$CDR_SB,-0.1),
                                             df_plot_index = cent_index,indicatorline = 0.3,labels_all = FALSE,title_plot = "Dementia severity")

cdr_plot_amypos <- cdr_plot_amypos+ggtitle("Dementia severity")
plot_out <- centiloid_plot_amypos | tauearly_plot_amypos | atrophy_plot_amypos | cdr_plot_amypos

ggsave(paste0(out_plots,"Amyloid_PET_Positive_cohort_ATN.jpg"),
       plot=plot_out,
       device=jpeg,width = 8,height=4.02)


####Amyloid Negative cohort plots
cent_index <- index_biomarker_forest(list_cent_new[[3]])

###### Centiloid
centiloid_plot_amyneg <- spearman_gg_formanuscript(list_cent_new[[3]],breaks1 = c(plot_breaks$Centiloid,0,-0.2),
                                                   df_plot_index = cent_index,indicatorline = 0.6,labels_all = TRUE,title_plot = "Amyloid PET")
centiloid_plot_amyneg <- centiloid_plot_amyneg+ggtitle("Amyloid PET")

##### Tau Early
tauearly_plot_amyneg <- spearman_gg_formanuscript(list_tauearly_new[[3]],breaks1 = c(plot_breaks$TauEarly,-0.2),
                                                  df_plot_index = cent_index,indicatorline = 0.4,labels_all = FALSE,title_plot = "Early tau PET")
tauearly_plot_amyneg <- tauearly_plot_amyneg+ggtitle("Early tau PET")

##### Atrophy
atrophy_plot_amyneg <- spearman_gg_formanuscript(list_atrophy_new[[3]],breaks1 = c(plot_breaks$Atrophy,-0.6),
                                                 df_plot_index = cent_index,indicatorline = -0.2,rev_axis = TRUE,labels_all = FALSE,
                                                 title_plot = "Cortical thickness")
atrophy_plot_amyneg <- atrophy_plot_amyneg+ggtitle("Cortical thickness")

##### CDR Sum of Boxes
cdr_plot_amyneg <- spearman_gg_formanuscript(list_cdr_new[[3]],breaks1 = c(plot_breaks$CDR_SB,-0.1),
                                             df_plot_index = cent_index,indicatorline = 0.3,labels_all = FALSE,title_plot = "Dementia severity")

cdr_plot_amyneg <- cdr_plot_amyneg+ggtitle("Dementia severity")
plot_out <- centiloid_plot_amyneg | tauearly_plot_amyneg | atrophy_plot_amyneg | cdr_plot_amyneg

ggsave(paste0(out_plots,"Amyloid_PET_negative_cohort_ATN.jpg"),
       plot=plot_out,
       device=jpeg,width = 8,height=4.02)



######Cognitively Impaired Cohort Plots
cent_index <- index_biomarker_forest(list_cent_new[[6]])

#####Centiloid
centiloid_plot_cdrpos <- spearman_gg_formanuscript(list_cent_new[[6]],breaks1 = plot_breaks$Centiloid,
                                                   df_plot_index = cent_index,indicatorline = 0.6,labels_all = TRUE,title_plot = "Amyloid PET",
                                                   amytag=TRUE)
centiloid_plot_cdrpos <- centiloid_plot_cdrpos+ggtitle("Amyloid PET")

##### Tau Early
tauearly_plot_cdrpos <- spearman_gg_formanuscript(list_tauearly_new[[6]],breaks1 = plot_breaks$TauEarly,
                                                  df_plot_index = cent_index,indicatorline = 0.4,labels_all = FALSE,title_plot = "Early tau PET")
tauearly_plot_cdrpos <- tauearly_plot_cdrpos+ggtitle("Early tau PET")


###### Atrophy
atrophy_plot_cdrpos <- spearman_gg_formanuscript(list_atrophy_new[[6]],breaks1 = plot_breaks$Atrophy,
                                                 df_plot_index = cent_index,indicatorline = -0.2,rev_axis = TRUE,labels_all = FALSE,
                                                 title_plot = "Cortical thickness")
atrophy_plot_cdrpos <- atrophy_plot_cdrpos+ggtitle("Cortical thickness")

######CDR Sum of Boxes
cdr_plot_cdrpos <- spearman_gg_formanuscript(list_cdr_new[[6]],breaks1 = plot_breaks$CDR_SB,
                                             df_plot_index = cent_index,indicatorline = 0.3,labels_all = FALSE,title_plot = "Dementia severity")

cdr_plot_cdrpos <- cdr_plot_cdrpos+ggtitle("Dementia severity")
plot_out <- centiloid_plot_cdrpos | tauearly_plot_cdrpos | atrophy_plot_cdrpos | cdr_plot_cdrpos

ggsave(paste0(out_plots,"Cognitively_impaired_cohort_ATN.jpg"),
       plot=plot_out,
       device=jpeg,width = 8,height=4.02)




######Cognitively Unimpaired Cohort Plots
cent_index <- index_biomarker_forest(list_cent_new[[7]])

##### Centiloid
centiloid_plot_cdrneg <- spearman_gg_formanuscript(list_cent_new[[7]],breaks1 = c(0,plot_breaks$Centiloid),
                                                   df_plot_index = cent_index,indicatorline = 0.6,labels_all = TRUE,title_plot = "Amyloid PET",
                                                   amytag=TRUE)
centiloid_plot_cdrneg <- centiloid_plot_cdrneg+ggtitle("Amyloid PET")


#####Tau Early
tauearly_plot_cdrneg <- spearman_gg_formanuscript(list_tauearly_new[[7]],breaks1 = plot_breaks$TauEarly,
                                                  df_plot_index = cent_index,indicatorline = 0.4,labels_all = FALSE,title_plot = "Early tau PET")
tauearly_plot_cdrneg <- tauearly_plot_cdrneg+ggtitle("Early tau PET")

#### Atrophy
atrophy_plot_cdrneg <- spearman_gg_formanuscript(list_atrophy_new[[7]],breaks1 = c(0.2,0.0,-0.2,-0.4,-0.6),
                                                 df_plot_index = cent_index,indicatorline = -0.2,rev_axis = TRUE,labels_all = FALSE,
                                                 title_plot = "Cortical thickness")
atrophy_plot_cdrneg <- atrophy_plot_cdrneg+ggtitle("Cortical thickness")

plot_out <- centiloid_plot_cdrneg | tauearly_plot_cdrneg | atrophy_plot_cdrneg

ggsave(paste0(out_plots,"Cognitively_unimpaired_cohort_ATN.jpg"),
       plot=plot_out,
       device=jpeg,width = 8,height=4.02)











####Centiloid All PLots
cent_index <- index_biomarker_forest(list_cent_new[[6]])

#### Amyloid Positive
centiloid_plot_amypos <- spearman_gg_formanuscript(list_cent_new[[2]],breaks1 = c(plot_breaks$Centiloid,0,-0.2),
                                                   df_plot_index = cent_index,indicatorline = 0.6,labels_all = FALSE,title_plot = "Amyloid PET\n positive")
#### Amyloid Negative
centiloid_plot_amyneg <- spearman_gg_formanuscript(list_cent_new[[3]],breaks1 = c(plot_breaks$Centiloid,0,-0.2),
                                                   df_plot_index = cent_index,indicatorline = 0.6,labels_all = FALSE,title_plot = "Amyloid PET\n negative")


#### CDR Positive
centiloid_plot_cdrpos <- spearman_gg_formanuscript(list_cent_new[[6]],breaks1 = c(plot_breaks$Centiloid,0,-0.2),
                                                   df_plot_index = cent_index,indicatorline = 0.6,labels_all = TRUE,title_plot = "Cognitively\n impaired",
                                                   amytag=TRUE)
#### CDR Negative
centiloid_plot_cdrneg <- spearman_gg_formanuscript(list_cent_new[[7]],breaks1 = c(plot_breaks$Centiloid,0,-0.2),
                                                   df_plot_index = cent_index,indicatorline = 0.6,labels_all = FALSE,title_plot = "Cognitively\n unimpaired",
                                                   amytag=TRUE)

#### Combine plots
plot_out <-  centiloid_plot_cdrpos | centiloid_plot_cdrneg | centiloid_plot_amypos | centiloid_plot_amyneg 


#### write out plots
ggsave(paste0(out_plots,"Centiloid_all_ATN.png"),
       plot=plot_out,
       device=png,dpi = 1200,width = 8,height=4.02)


