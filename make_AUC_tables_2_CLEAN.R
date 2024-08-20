# FNIH Project 2
# make AUC tables
#### Coded: Kellen Petersen
#### Date as of commenting and documentation - 05/28/2024

# Load libraries
library(tidyverse)
library(data.table)
library(RVAideMemoire)
library(pROC)
library(RColorBrewer)
library(MuMIn)
library(here)
library(cowplot)
library(ggeasy)

#here::i_am("make_AUC_tables_2_CLEAN.R")

# Load data
data_file <- here("adni_data","cross_sectional_data_CLEAN.csv")
data.in00 <- read.csv(data_file)

# Define output folder
out_folder <- "results_clean"

# Choose outcome
info1 <- "outcomeCENT20"
#info1 <- "outcomeCENT37"
#info1 <- "outcomeEarlyTau"
#info1 <- "outcomeLateTau"
#info1 <- "outcomeAtrophy"
#info1 <- "outcomeImpairment"

# Choose cohort
info2 <- "cohortAll"
#info2 <- "cohortCU"
#info2 <- "cohortCI"
#info2 <- "cohortAB20n"
#info2 <- "cohortAB20p"
#info2 <- "cohortAB37n"
#info2 <- "cohortAB37p"

# Make output file names
outfile_singlecutoff <- paste0(info1,"_",info2,"_singleCutoff",".csv")
outfile_AUCs <- paste0(info1,"_",info2,"_AUCs",".csv")
outfile_twoCuts <- paste0(info1,"_",info2,"_twoCuts",".csv")

# Make cohort dataframe
data.in0 <- data.in00
#data.in0 <- data.in00 %>% filter(CDR == 0)
#data.in0 <- data.in00 %>% filter(CDR > 0)
#data.in0 <- data.in00 %>% filter(CENTILOIDS_10 == 0)
#data.in0 <- data.in00 %>% filter(CENTILOIDS_10 == 1)
#data.in0 <- data.in00 %>% filter(CENTILOIDS_10_sens == 0)
#data.in0 <- data.in00 %>% filter(CENTILOIDS_10_sens == 1)

# Define binary variable for cognitive impairment
data.in0$impaired_10 <- ifelse(data.in0$CDR > 0, 1, 0)

# Rename dataframe
data.in.cohort <- data.in0

# Parameters (covariates included?)
covars_TF <- FALSE

# Outcome variables
vars_outcome <- c("CENTILOIDS_10")
#vars_outcome <- c("CENTILOIDS_10_sens")
#vars_outcome <- c("TAU_MesialTemporal_10")
#vars_outcome <- c("TAU_TemporoParietal_10")
#vars_outcome <- c("atrophy_10")
#vars_outcome <- c("impaired_10")

#Covariates
age_var <- "AGE" #variable with age
sex_var <- "PTGENDER" #variable with sex
educ_var <- "PTEDUCAT" #variable with education
APOE_variable <- "APOE_genotype" #Variable with APOE genotypes
data.in.cohort$APOE_Cat <- factor(data.in.cohort[[APOE_variable]], 
                                  levels = c("33", "34", "24", "22", "23", "44"))
list_covariates <- list(c(age_var, sex_var, "APOE_Cat"))
if (covars_TF) {
  vars_covariates <- list_covariates
} else {
  vars_covariates <- NULL
}

# Assays Lists
assays_c2n       <- c("C2N_plasma_Abeta42_Abeta40", "C2N_plasma_ptau217", "C2N_plasma_ptau217_ratio")
assays_fuji      <- c("Fuji_plasma_Ab42_Ab40","Fuji_plasma_ptau217")
assays_alzpath   <- c("AlzPath_plasma_ptau217")
assays_janssen   <- c("Janssen_plasma_ptau217")
assays_roche     <- c("Roche_plasma_Ab42_Ab40","Roche_plasma_ptau181","Roche_plasma_NfL","Roche_plasma_GFAP")
assays_quanterix <- c("QX_plasma_Ab42_Ab40","QX_plasma_ptau181","QX_plasma_NfL","QX_plasma_GFAP")
assays_csfratio  <- c("PTAU_over_ABETA42")

# Create list of model indices
combos_c2n       <- list(assays_c2n[c(3,1)],assays_c2n[c(2,1)],
                         assays_c2n[c(3)],assays_c2n[c(2)],assays_c2n[c(1)])              # 5 models
combos_fuji      <- list(assays_fuji[c(2,1)],assays_fuji[c(2)],assays_fuji[c(1)])         # 3 models
combos_alzath    <- list(assays_alzpath)                                                  # 1 model
combos_janssen   <- list(assays_janssen)                                                  # 1 model
combos_roche     <- list(assays_roche[c(2,1,4,3)],assays_roche[c(2,1,3)],
                         assays_roche[c(2,1)],
                         assays_roche[c(2)],assays_roche[c(1)],
                         assays_roche[c(4)],assays_roche[c(3)])                           # 7 models
combos_quanterix <- list(assays_quanterix[c(2,1,4,3)],assays_quanterix[c(2,1,3)],
                         assays_quanterix[c(2,1)],
                         assays_quanterix[c(2)],assays_quanterix[c(1)],
                         assays_quanterix[c(4)],assays_quanterix[c(3)])                   # 7 models
combos_csfratio   <- list(assays_csfratio)                                                # 1 model

combos_all <- list()
combos_all[[1]] <- combos_c2n[1]
combos_all[[2]] <- combos_c2n[2]
combos_all[[3]] <- combos_c2n[3]
combos_all[[4]] <- combos_c2n[4]
combos_all[[5]] <- combos_c2n[5]

combos_all[[6]] <- combos_fuji[1]
combos_all[[7]] <- combos_fuji[2]
combos_all[[8]] <- combos_fuji[3]

combos_all[[9]] <- combos_alzath[1]

combos_all[[10]] <- combos_janssen[1]

combos_all[[11]] <- combos_roche[1]
combos_all[[12]] <- combos_roche[2]
combos_all[[13]] <- combos_roche[3]
combos_all[[14]] <- combos_roche[4]
combos_all[[15]] <- combos_roche[5]
combos_all[[16]] <- combos_roche[6]
combos_all[[17]] <- combos_roche[7]

combos_all[[18]] <- combos_quanterix[1]
combos_all[[19]] <- combos_quanterix[2]
combos_all[[20]] <- combos_quanterix[3]
combos_all[[21]] <- combos_quanterix[4]
combos_all[[22]] <- combos_quanterix[5]
combos_all[[23]] <- combos_quanterix[6]
combos_all[[24]] <- combos_quanterix[7]

Combos <- combos_all

# Company label
C_c2n <- rep("C2N", 5)
C_fuji <- rep("Fujirebio", 3)
C_alzpath <- rep("AlzPath", 1)
C_janssen <- rep("Janssen", 1)
C_roche <- rep("Roche", 7)
C_5csfratio <- rep("CSF", 1)
C_quanterix <- rep("Quanterix", 7)
C <- c(C_c2n, C_fuji, C_alzpath, C_janssen, C_roche, C_quanterix)
CC <- C

# Analyte label
A_c2n <- c("p-tau217 ratio + Ab42/Ab40","p-tau217 + Ab42/Ab40",
           "p-tau217 ratio","p-tau217","Ab42/Ab40")
A_fuji <- c("p-tau217 + Ab42/Ab40", "p-tau217", "Ab42/Ab40")
A_alzpath <- c("p-tau217")
A_janssen <- c("p-tau217")
A_roche <- c("p-tau181 + Ab42/Ab40 + GFAP + NfL",
             "p-tau181 + Ab42/Ab40 + NfL",
             "p-tau181 + Ab42/Ab40",
             "p-tau181","Ab42/Ab40",
             "GFAP","NfL")
A_quanterix <- c("p-tau181 + Ab42/Ab40 + GFAP + NfL",
                 "p-tau181 + Ab42/Ab40 + NfL",
                 "p-tau181 + Ab42/Ab40",
                 "p-tau181","Ab42/Ab40",
                 "GFAP","NfL")
A <- c(A_c2n, A_fuji, A_alzpath, A_janssen, A_roche, A_quanterix)
AA <- A

# List of single analytes
single_analytes <- c("p-tau217 ratio",
                     "Ab42/Ab40",
                     "p-tau217",
                     "p-tau181",
                     "GFAP", 
                     "NfL",
                     "CSF PTAU/AB42")

OUT_df <- list()
MODELS_df <- list()
k <- 0
data.in <- data.in.cohort

# Loop through the combinations of predictors
for (i_outcome in 1:length(vars_outcome)) {
  for (i_assay_combos in 1:length(Combos)) {
    
    outcome <- vars_outcome[i_outcome]
    pred1 <- Combos[[i_assay_combos]]
    pred2 <- vars_covariates
    preds <- c(unlist(pred1), unlist(pred2))
    
    formula_str <- paste(outcome, " ~ ", paste(preds, collapse = " + "))
    model_formula <- as.formula(formula_str)
    
    temp_data <- data.in[complete.cases(data.in[,preds]), c(outcome,preds)] %>% 
      na.omit()
    
    mdl_logistic <- glm(model_formula, data = temp_data, family = "binomial")
    mdl_pred <- predict(mdl_logistic, data = temp_data, type = "response")
    
    roc_tmp <- pROC::roc(temp_data[[outcome]],mdl_pred,plot=FALSE,print.auc=FALSE) 
    youden_index <- coords(roc_tmp, "best", ret="all")
    optimal_cutoff <- youden_index$threshold

    T <- data.frame(roc_tmp$thresholds, roc_tmp$sensitivities, roc_tmp$specificities)
    auc_tmp <- pROC::auc(roc_tmp)
    ci_AUC  <- ci.auc(roc_tmp) #get confidence interval
    
    sensitivity <- youden_index$sensitivity
    specificity <- youden_index$specificity
    
    # Sensitivity = 90%
    II_sensX <- which(roc_tmp$sensitivities >= .90)
    I_sensX <- max(II_sensX)
    sensX_threshold <- roc_tmp$thresholds[I_sensX]
    sensX_specificity <- roc_tmp$specificities[I_sensX]
    
    # Specificity = 90%
    II_specX <- which(roc_tmp$specificities >= .90)
    I_specX <- min(II_specX)
    specX_threshold <- roc_tmp$thresholds[I_specX]
    specX_sensitivities <- roc_tmp$sensitivities[I_specX]
    
    num_intermediate <- sum( (mdl_pred > sensX_threshold) & (mdl_pred < specX_threshold) )
    num_intermediate
    
    out_v <- c()
    out_v[1] <- CC[i_assay_combos] # company label
    out_v[2] <- round(roc_tmp$auc,3) # AUC rounded to 3 decimals
    out_v[3] <- round(ci_AUC[1],3) # AUC lower bound
    out_v[4] <- round(ci_AUC[3],3) # AUC upper bound
    out_v[5] <- paste0(formatC(sensitivity, digits=3, format="fg"))
    out_v[6] <- paste0(formatC(specificity, digits=3, format="fg"))
    out_v[7] <- paste0(formatC(youden_index$ppv, digits=3, format="fg"))
    out_v[8] <- paste0(formatC(youden_index$npv, digits=3, format="fg"))
    out_v[9] <- paste0(formatC(youden_index$accuracy, digits=3, format="fg"))
    
    tmp <- formatC((log(optimal_cutoff / (1 - optimal_cutoff)) - mdl_logistic$coefficients[1]) / mdl_logistic$coefficients[2], digits=4, format="fg")
    if (AA[i_assay_combos] %in% single_analytes) {out_v[10] <- tmp} 
    else {out_v[10] <- -99}
    
    out_v[11] <- AA[i_assay_combos] # analyte label
    
    out_v[12] <- formatC(sensX_threshold, digits=3, format="fg")
    out_v[13] <- formatC(sensX_specificity, digits=3, format="fg")
    
    out_v[14] <- formatC(specX_threshold, digits=3, format="fg")
    out_v[15] <- formatC(specX_sensitivities, digits=3, format="fg")
    
    out_v[16] <- num_intermediate
    
    out_v[17] <- formula_str
    out_v[18] <- outcome
    out_v[19] <- paste(pred1, collapse = " + ")
    out_v[20] <- covars_TF
    
    dd <- dim(temp_data)
    out_v[21] <- dd[1]
    
    brierScore <- mean((mdl_pred-temp_data[[outcome]])^2)
    out_v[22] <- brierScore
    
    model <- list()
    model[[1]] <- mdl_logistic
    model[[2]] <- roc_tmp
    
    k <- i_assay_combos + length(combos_all)*(i_outcome-1)
    OUT_df[[k]] <- out_v
    MODELS_df[[k]] <- model
  }
}



################################################################################


df_out <- t(data.frame(OUT_df))
df_out <- data.frame(df_out)
rownames(df_out) <- NULL
colnames(df_out) <- c("company","AUC","AUC_Lower","AUC_Upper",
                      "sensitivity","specificity","ppv","npv","accuracy",
                      "cutoff_value","analyte",
                      "90% Sens Threshold","90% Sens - Specificity",
                      "90% Spec Threshold","90% Spec - Sensitivity",
                      "Num_Intermediate","formula","outcome","predictors",
                      "covars_TF","n_cases","Brier_Score")

df_out$AUC <- as.numeric(df_out$AUC)
df_out$AUC_Lower <- as.numeric(df_out$AUC_Lower)
df_out$AUC_Upper <- as.numeric(df_out$AUC_Upper)
df_out$Num_Intermediate <- as.numeric(df_out$Num_Intermediate)
df_out$n_cases <- as.numeric(df_out$n_cases)
df_out <- df_out %>% 
  mutate(intermediate_results = paste(Num_Intermediate, " (", round(Num_Intermediate/n_cases,2)*100, "%)", sep = ""))
DF <- df_out

DFB <- DF %>% 
  mutate(index_0 = row_number()) %>% 
  filter(outcome == vars_outcome[1] ) %>% 
  arrange(company,AUC) %>% 
  group_by(company) %>%
  mutate(index = row_number())
DFB$company <- factor(DFB$company, levels = c("C2N", "Fujirebio", "AlzPath","Janssen", "Roche", "Quanterix","CSF"))

DFB <- DFB %>% 
  mutate(companyNums = as.numeric(company)) %>% 
  arrange(companyNums,desc(AUC)) %>% 
  ungroup() %>% 
  mutate(plot_index = desc(row_number()))

DFB$analyte <- factor(DFB$analyte,
                      levels = c("p-tau217 ratio + Ab42/Ab40",
                                 "p-tau217 ratio",
                                 "Ab42/Ab40",
                                 "p-tau217 + Ab42/Ab40",
                                 "p-tau217",
                                 "p-tau181 + Ab42/Ab40 + GFAP + NfL",
                                 "p-tau181 + Ab42/Ab40 + NfL",
                                 "p-tau181 + Ab42/Ab40",
                                 "p-tau181",
                                 "GFAP", "NfL",
                                 "CSF PTAU/AB42"))

DFB2 <- DFB %>%
  mutate(AUC_CI = paste0(format(AUC,nsmall=3)," (",format(AUC_Lower,nsmall=3), "-",format(AUC_Upper,nsmall=3),")")) %>%
  select("company","analyte",
         "AUC_CI",
         "sensitivity","specificity",
         "90% Sens Threshold","90% Sens - Specificity",
         "90% Spec Threshold","90% Spec - Sensitivity",
         "intermediate_results")



################################################################################
#  Single analyte table
################################################################################
T1 <- DFB %>% 
  mutate(AUC_CI = paste0(formatC(AUC, digits=3, format="fg")," (",
                         formatC(AUC_Lower, digits=3, format="fg"), "-",
                         format(AUC_Upper, digits=3, format="fg"),")")) %>% 
  select("company","analyte","AUC",
         "AUC_CI","cutoff_value",
         "sensitivity","specificity","accuracy","Brier_Score","ppv","npv") %>% 
  filter(analyte %in% single_analytes) %>% 
  rename("Platform" = "company",
         "Analyte" = "analyte",
         "AUCraw" = "AUC",
         "AUC" = "AUC_CI",
         "Cutoff Value" = "cutoff_value",
         "Sensitivity" = "sensitivity",
         "Specificity" = "specificity",
         "Accuracy" = "accuracy",
         "Brier Score" = "Brier_Score",
         "PPV" = "ppv",
         "NPV" = "npv",) %>% 
  filter(Platform != "CSF") 

T2 <- T1 %>%
  group_by(Platform) %>%
  mutate(max_AUCraw = max(AUCraw)) %>%
  ungroup() %>% 
  arrange(desc(max_AUCraw), Platform, desc(AUCraw)) %>% 
  select(-max_AUCraw,-AUCraw)

write.csv(T2, file = here(out_folder,outfile_singlecutoff), row.names = FALSE)





################################################################################
################## Make AUC Table: with & w/o covariates  ######################
################################################################################

AUCunadj <- DFB %>% 
  mutate(AUC_CI = paste0(format(AUC,nsmall=3)," (",format(AUC_Lower,nsmall=3), "-",format(AUC_Upper,nsmall=3),")")) %>%
  select("company","analyte","AUC_CI","AUC") %>% 
  group_by(company) %>%
  mutate(max_AUC = max(AUC)) %>%
  ungroup() %>%
  arrange(desc(max_AUC),company,desc(AUC))%>%
  select(-max_AUC,-AUC) %>% 
  rename("AUC (unadjusted)" = "AUC_CI")

covars_TF <- TRUE
if (covars_TF) {
  vars_covariates <- list_covariates
} else {
  vars_covariates <- NULL
}

OUT_5adj <- list()
MODELS_5adj <- list()
k <- 0
data.in <- data.in.cohort
# Loop through the combinations of predictors
for (i_outcome in 1:length(vars_outcome)) {
  for (i_assay_combos in 1:length(Combos)) {
    
    outcome <- vars_outcome[i_outcome]
    pred1 <- Combos[[i_assay_combos]]
    pred2 <- vars_covariates
    preds <- c(unlist(pred1), unlist(pred2))
    
    formula_str <- paste(outcome, " ~ ", paste(preds, collapse = " + "))
    model_formula <- as.formula(formula_str)
    
    temp_data <- data.in[complete.cases(data.in[,preds]), c(outcome,preds)] %>% 
      na.omit()
    
    mdl_logistic <- glm(model_formula, data = temp_data, family = "binomial")
    mdl_pred <- predict(mdl_logistic, data = temp_data, type = "response")
    
    roc_tmp <- pROC::roc(temp_data[[outcome]],mdl_pred,plot=FALSE,print.auc=FALSE)
    youden_index <- coords(roc_tmp, "best", ret="all")
    optimal_cutoff <- youden_index$threshold
    
    T <- data.frame(roc_tmp$thresholds, roc_tmp$sensitivities, roc_tmp$specificities)
    auc_tmp <- pROC::auc(roc_tmp)
    ci_AUC  <- ci.auc(roc_tmp) #get confidence interval
    
    sensitivity <- youden_index$sensitivity
    specificity <- youden_index$specificity
    
    # Sensitivity = 90%
    II_sensX <- which(roc_tmp$sensitivities >= .90)
    I_sensX <- max(II_sensX)
    sensX_threshold <- roc_tmp$thresholds[I_sensX]
    sensX_specificity <- roc_tmp$specificities[I_sensX]
    
    # Specificity = 90%
    II_specX <- which(roc_tmp$specificities >= .90)
    I_specX <- min(II_specX)
    specX_threshold <- roc_tmp$thresholds[I_specX]
    specX_sensitivities <- roc_tmp$sensitivities[I_specX]
    
    num_intermediate <- sum( (mdl_pred > sensX_threshold) & (mdl_pred < specX_threshold) )
    
    out_v <- c()
    out_v[1] <- CC[i_assay_combos] # company label
    out_v[2] <- formatC(roc_tmp$auc, digits=3, format="fg") # AUC rounded to 2 decimals
    out_v[3] <- formatC(ci_AUC[1], digits=3, format="fg") # AUC lower bound
    out_v[4] <- formatC(ci_AUC[3], digits=3, format="fg") # AUC upper bound
    out_v[5] <- formatC(sensitivity, digits=3, format="fg")
    out_v[6] <- formatC(specificity, digits=3, format="fg")
    out_v[7] <- formatC(youden_index$ppv, digits=3, format="fg")
    out_v[8] <- formatC(youden_index$npv, digits=3, format="fg")
    out_v[9] <- formatC(youden_index$accuracy, digits=3, format="fg")
    
    tmp <- formatC((log(optimal_cutoff / (1 - optimal_cutoff)) - mdl_logistic$coefficients[1]) / mdl_logistic$coefficients[2],digits=4,format="fg")
    if (AA[i_assay_combos] %in% single_analytes) {out_v[10] <- tmp} 
    else {out_v[10] <- -99}
    out_v[11] <- AA[i_assay_combos] # analyte label
    
    out_v[12] <- formatC(sensX_threshold, digits=3, format="fg")
    out_v[13] <- formatC(sensX_specificity, digits=3, format="fg")
    
    out_v[14] <- formatC(specX_threshold, digits=3, format="fg")
    out_v[15] <- formatC(specX_sensitivities, digits=3, format="fg")
    
    out_v[16] <- num_intermediate
    
    out_v[17] <- formula_str
    out_v[18] <- outcome
    out_v[19] <- paste(pred1, collapse = " + ")
    out_v[20] <- covars_TF
    dd <- dim(temp_data)
    out_v[21] <- dd[1]
    
    brierScore <- mean((mdl_pred-temp_data[[outcome]])^2)
    out_v[22] <- brierScore
    
    model <- list()
    model[[1]] <- mdl_logistic
    model[[2]] <- roc_tmp
    
    k <- i_assay_combos + length(combos_all)*(i_outcome-1)
    OUT_5adj[[k]] <- out_v
    MODELS_5adj[[k]] <- model
  }
}



# dataframce for models adjusted (a) for covariates
df_outa <- t(data.frame(OUT_5adj))
df_outa <- data.frame(df_outa)
rownames(df_outa) <- NULL
colnames(df_outa) <- c("company","AUC","AUC_Lower","AUC_Upper",
                       "sensitivity","specificity","ppv","npv","accuracy",
                       "cutoff_value","analyte",
                       "90% Sens Threshold","90% Sens - Specificity",
                       "90% Spec Threshold","90% Spec - Sensitivity",
                       "Num_Intermediate","formula","outcome","predictors",
                       "covars_TF","n_cases","Brier_Score")

df_outa$AUC <- as.numeric(df_outa$AUC)
df_outa$AUC_Lower <- as.numeric(df_outa$AUC_Lower)
df_outa$AUC_Upper <- as.numeric(df_outa$AUC_Upper)
df_outa$Num_Intermediate <- as.numeric(df_outa$Num_Intermediate)
df_outa$n_cases <- as.numeric(df_outa$n_cases)
df_outa <- df_outa %>% 
  mutate(intermediate_results = paste(Num_Intermediate, " (", round(Num_Intermediate/n_cases,2)*100, "%)", sep = ""))

DFa <- df_outa

DFBa <- DFa %>% 
  mutate(index_0 = row_number()) %>% 
  filter(outcome == vars_outcome[1] ) %>% 
  arrange(company,AUC) %>% 
  group_by(company) %>%
  mutate(index = row_number())
DFBa$company <- factor(DFBa$company, 
                       levels = c("C2N", "Fujirebio", "AlzPath","Janssen", "Roche", "Quanterix","CSF"))

DFBa <- DFBa %>% 
  mutate(companyNums = as.numeric(company)) %>% 
  arrange(companyNums,desc(AUC)) %>% 
  ungroup() %>% 
  mutate(plot_index = desc(row_number()))
DFBa$analyte <- factor(DFBa$analyte,
                       levels = c("p-tau217 ratio + Ab42/Ab40",
                                  "p-tau217 ratio",
                                  "Ab42/Ab40",
                                  "p-tau217 + Ab42/Ab40",
                                  "p-tau217",
                                  "p-tau181 + Ab42/Ab40 + GFAP + NfL",
                                  "p-tau181 + Ab42/Ab40 + NfL",
                                  "p-tau181 + Ab42/Ab40",
                                  "p-tau181",
                                  "GFAP", "NfL",
                                  "CSF PTAU/AB42"))



DFBa2 <- DFBa %>% 
  mutate(AUC_CI = paste0(AUC," (",AUC_Lower, "-",AUC_Upper,")")) %>% 
  select("company","analyte",
         "AUC_CI",
         "sensitivity","specificity",
         "90% Sens Threshold","90% Sens - Specificity",
         "90% Spec Threshold","90% Spec - Sensitivity",
         "intermediate_results")
AUCadj <- DFBa %>% 
  mutate(AUC_CI = paste0(format(AUC,nsmall=3)," (",format(AUC_Lower,nsmall=3), "-",format(AUC_Upper,nsmall=3),")")) %>%
  select("company","analyte","AUC_CI","AUC") %>% 
  group_by(company) %>%
  mutate(max_AUC = max(AUC)) %>%
  ungroup() %>%
  arrange(desc(max_AUC),company,desc(AUC))%>%
  select(-max_AUC,-AUC) %>% 
  rename("AUC (adjusted)" = "AUC_CI")

AUCboth <- full_join(AUCunadj,AUCadj, by = c("company","analyte"))
AUCboth <- AUCboth %>% 
  filter(company != "CSF") %>%
  rename(Platform = company,
         Analytes = analyte)

AUCboth$p_unadj <- NA
AUCboth$p_adj <- NA

REF <- data.frame(Platform = C, Analytes = A) %>% 
  mutate(index = row_number()) %>% 
  filter(Platform != "CSF")
Platforms <- unique(AUCboth$Platform)
Platforms <- Platforms[Platforms != c("AlzPath","Janssen")]

# Loop through each company and perform Delong's test
for (platform in Platforms) {
  # Filter the data for the current company
  platform_data <- AUCboth[AUCboth$Platform == platform, ]
  
  # Extract the AUC values
  ref_platform <- platform_data$Platform[1]
  ref_analyte <- platform_data$Analytes[1]
  
  i_ref <- which(REF$Platform == ref_platform & REF$Analytes == ref_analyte)
  roc0 <- MODELS_df[[i_ref]][[2]]
  
  a_list <- platform_data$Analytes
  a_list <- a_list[-1]
  ps0 <- rep(NA, length(a_list))
  ps00 <- rep(NA, length(a_list))
  ps <- rep(NA, length(a_list))
  psAdj <- rep(NA, length(a_list))
  
  j <- 1
  for (a in a_list) {
    i_a <- which(REF$Platform == ref_platform & REF$Analytes == a)

    roc1 <- MODELS_df[[i_a]][[2]]
    t <- roc.test(roc0, roc1, method = "delong")
    ps0[j] <- t$p.value
    j <- j + 1
  }
  p_adjusted <- p.adjust(ps0, method="BH")
  pRaw <- c(NA, ps0)
  pAdj <- c(NA, p_adjusted)
  
  a_list <- platform_data$Analytes
  k0 <- which(AUCboth$Platform == platform & AUCboth$Analytes == a_list[1])
  kN <- which(AUCboth$Platform == platform & AUCboth$Analytes == a_list[length(a_list)])
  AUCboth$p_unadj[k0:kN] <- pRaw
  AUCboth$p_adj[k0:kN] <- pAdj
}

AUCbest <- AUCboth %>%
  mutate(original_order = row_number()) %>% 
  group_by(Platform) %>%
  slice(1) %>% 
  ungroup() %>% 
  arrange(original_order) %>%
  select(Platform,Analytes,-original_order)
AUCbest$p_best <- NA

for (i in 2:nrow(AUCbest)) {
  ref_platform <- AUCbest$Platform[1]
  ref_analyte <- AUCbest$Analytes[1]
  i_ref <- which(REF$Platform == ref_platform & REF$Analytes == ref_analyte)

  roc0 <- MODELS_df[[i_ref]][[2]]
  
  comp_platform <- AUCbest$Platform[i]
  comp_analyte <- AUCbest$Analytes[i]
  i_comp <- which(REF$Platform == comp_platform & REF$Analytes == comp_analyte)

  roc1 <- MODELS_df[[i_comp]][[2]]
  
  t <- roc.test(roc0, roc1, method = "delong")
  p <- t$p.value
  AUCbest$p_best[i] <- p
}

pbest <- AUCbest$p_best
pbest <- pbest[-1]

pbest_adjusted <- p.adjust(pbest, method="BH")
AUCbest$p_best_adj <- c(NA, pbest_adjusted)

AUCboth2 <- AUCboth %>% 
  left_join(AUCbest, by = c("Platform","Analytes"))

T_auc_noCov <- AUCboth2 %>% 
  mutate(p_unadj = ifelse(p_unadj<0.0001,"<0.0001",
                          ifelse(p_unadj>0.1,round(p_unadj,2),
                                 ifelse(is.na(p_unadj),NA,format(round(p_unadj, 4), scientific = FALSE)))),
         p_adj = ifelse(p_adj<0.0001,"<0.0001",
                        ifelse(p_adj>0.1,round(p_adj,2),
                               ifelse(is.na(p_adj),NA,format(round(p_adj, 4), scientific = FALSE)))),
         p_best = ifelse(p_best<0.0001,"<0.0001",
                         ifelse(p_best>0.1,round(p_best,2),
                                ifelse(is.na(p_best),NA,format(round(p_best, 4), scientific = FALSE)))),
         p_best_adj = ifelse(p_best_adj<0.0001,"<0.0001",
                             ifelse(p_best_adj>0.1,round(p_best_adj,2),
                                    ifelse(is.na(p_best_adj),NA,format(round(p_best_adj, 4), scientific = FALSE))))) %>% 
  rename("No Cov - unadjusted p=" = p_unadj,
         "No Cov - adjusted p=" = p_adj,
         "No Cov - Best unadjusted  p=" = p_best,
         "No Cov - Best adjusted p=" = p_best_adj)

################################################################################
################################################################################

AUCadj <- AUCadj %>% 
  rename(Platform = company, Analytes = analyte) %>% 
  filter(Platform != "CSF")

AUCboth <- AUCadj %>% #### Replace AUCboth with AUCadj
  select(Platform,Analytes)
AUCboth$p_unadjCov <- NA
AUCboth$p_adjCov <- NA

REF <- data.frame(Platform = C, Analytes = A) %>% 
  mutate(index = row_number()) %>% 
  filter(Platform != "CSF")
Platforms <- unique(AUCadj$Platform)
Platforms <- Platforms[Platforms != c("AlzPath","Janssen")]

# Loop through each company and perform Delong's test
for (platform in Platforms) {
  # Filter the data for the current company
  platform_data <- AUCboth[AUCboth$Platform == platform, ]
  
  # Extract the AUC values
  ref_platform <- platform_data$Platform[1]
  ref_analyte <- platform_data$Analytes[1]
  
  i_ref <- which(REF$Platform == ref_platform & REF$Analytes == ref_analyte)
  roc0 <- MODELS_5adj[[i_ref]][[2]]
  
  a_list <- platform_data$Analytes
  a_list <- a_list[-1]
  ps0 <- rep(NA, length(a_list))
  ps00 <- rep(NA, length(a_list))
  ps <- rep(NA, length(a_list))
  psAdj <- rep(NA, length(a_list))
  
  j <- 1
  for (a in a_list) {
    i_a <- which(REF$Platform == ref_platform & REF$Analytes == a)

    roc1 <- MODELS_5adj[[i_a]][[2]]
    t <- roc.test(roc0, roc1, method = "delong")
    ps0[j] <- t$p.value
    j <- j + 1
  }
  p_adjusted <- p.adjust(ps0, method="BH")
  pRaw <- c(NA, ps0)
  pAdj <- c(NA, p_adjusted)
  
  a_list <- platform_data$Analytes
  k0 <- which(AUCboth$Platform == platform & AUCboth$Analytes == a_list[1])
  kN <- which(AUCboth$Platform == platform & AUCboth$Analytes == a_list[length(a_list)])
  AUCboth$p_unadjCov[k0:kN] <- pRaw
  AUCboth$p_adjCov[k0:kN] <- pAdj
}

AUCbest <- AUCboth %>%
  mutate(original_order = row_number()) %>% 
  group_by(Platform) %>%
  slice(1) %>% 
  ungroup() %>% 
  arrange(original_order) %>%
  select(Platform,Analytes,-original_order)
AUCbest$p_bestCov <- NA

for (i in 2:nrow(AUCbest)) {
  ref_platform <- AUCbest$Platform[1]
  ref_analyte <- AUCbest$Analytes[1]
  i_ref <- which(REF$Platform == ref_platform & REF$Analytes == ref_analyte)
  roc0 <- MODELS_5adj[[i_ref]][[2]]
  
  comp_platform <- AUCbest$Platform[i]
  comp_analyte <- AUCbest$Analytes[i]
  i_comp <- which(REF$Platform == comp_platform & REF$Analytes == comp_analyte)

  roc1 <- MODELS_5adj[[i_comp]][[2]]
  
  t <- roc.test(roc0, roc1, method = "delong")
  p <- t$p.value
  AUCbest$p_bestCov[i] <- p
}

pbest <- AUCbest$p_bestCov
pbest <- pbest[-1]

pbest_adjusted <- p.adjust(pbest, method="BH")
AUCbest$p_best_adjCov <- c(NA, pbest_adjusted)

AUCboth2 <- AUCboth %>% 
  left_join(AUCbest, by = c("Platform","Analytes"))

T_auc_Cov <- AUCboth2 %>% 
  mutate(p_unadjCov = ifelse(p_unadjCov<0.0001,"<0.0001",
                             ifelse(p_unadjCov>0.1,round(p_unadjCov,2),
                                    ifelse(is.na(p_unadjCov),NA,format(round(p_unadjCov, 4), scientific = FALSE)))),
         p_adjCov = ifelse(p_adjCov<0.0001,"<0.0001",
                           ifelse(p_adjCov>0.1,round(p_adjCov,2),
                                  ifelse(is.na(p_adjCov),NA,format(round(p_adjCov, 4), scientific = FALSE)))),
         p_bestCov = ifelse(p_bestCov<0.0001,"<0.0001",
                            ifelse(p_bestCov>0.1,round(p_bestCov,2),
                                   ifelse(is.na(p_bestCov),NA,format(round(p_bestCov, 4), scientific = FALSE)))),
         p_best_adjCov = ifelse(p_best_adjCov<0.0001,"<0.0001",
                                ifelse(p_best_adjCov>0.1,round(p_best_adjCov,2),
                                       ifelse(is.na(p_best_adjCov),NA,format(round(p_best_adjCov, 4), scientific = FALSE))))) %>% 
  rename("Cov - unadjusted p=" = p_unadjCov,
         "Cov - adjusted p=" = p_adjCov,
         "Cov - Best unadjusted  p=" = p_bestCov,
         "Cov - Best adjusted p=" = p_best_adjCov)

T_auc <- full_join(T_auc_noCov, T_auc_Cov, by = c("Platform","Analytes"))

T_auc2 <- T_auc %>% 
  select(Platform, Analytes, 
         "AUC (unadjusted)", "No Cov - adjusted p=", "No Cov - Best adjusted p=",
         "AUC (adjusted)","Cov - adjusted p=","Cov - Best adjusted p=")

write.csv(T_auc2, file = here(out_folder,outfile_AUCs), row.names = FALSE)

