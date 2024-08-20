# FNIH Project 2
# get AUCs for covariates-only models
#### Coded: Kellen Petersen
#### Date as of commenting and documentation - 05/28/2024

# Load libraries
library(tidyverse)
library(pROC)
library(here)

#here::i_am("get_AUC_covariatesOnly.R")

data_file <- here("adni_data","cross_sectional_data_CLEAN.csv")
data.in00 <- read.csv(data_file)

# Define binary variable for cognitive impairment
data.in00$impaired_10 <- ifelse(data.in00$CDR > 0, 1, 0)

# Define cohort
data.in0 <- data.in00
# data.in0 <- data.in00 %>% filter(CDR == 0)
# data.in0 <- data.in00 %>% filter(CDR > 0)
# data.in0 <- data.in00 %>% filter(CENTILOIDS_10 == 0)
# data.in0 <- data.in00 %>% filter(CENTILOIDS_10 == 1)
# data.in0 <- data.in00 %>% filter(CENTILOIDS_10_sens == 0)
# data.in0 <- data.in00 %>% filter(CENTILOIDS_10_sens == 1)

# Parameters for covariates (included or not?)
covars_TF <- TRUE

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
data.in0$APOE_Cat <- factor(data.in0[[APOE_variable]], levels = c("33", "34", "24", "22", "23", "44"))
list_covariates <- list(c(age_var, sex_var, "APOE_Cat"))
if (covars_TF) {
  vars_covariates <- list_covariates
} else {
  vars_covariates <- NULL
}

# Dataset that goes into the model
data.in <- data.in0

# Make model formula
outcome <- vars_outcome
pred1 <- NULL
pred2 <- vars_covariates
preds <- c(unlist(pred1), unlist(pred2))
formula_str <- paste(outcome, " ~ ", paste(preds, collapse = " + "))
model_formula <- as.formula(formula_str)

# Make relevant dataset
temp_data <- data.in[complete.cases(data.in[,preds]), c(outcome,preds)] %>% 
  na.omit()

# Fit logistic regression model
mdl_logistic <- glm(model_formula, data = temp_data, family = "binomial")
mdl_pred <- predict(mdl_logistic, data = temp_data, type = "response")

# ROC and AUC
roc_tmp <- pROC::roc(temp_data[[outcome]],mdl_pred,plot=FALSE,print.auc=FALSE) #this will get the ROC output
auc_tmp <- pROC::auc(roc_tmp)
ci_AUC  <- ci.auc(roc_tmp)

# Output
out_v <- c()
out_v[1] <- round(roc_tmp$auc,3) # AUC rounded to 3 decimals
out_v[2] <- round(ci_AUC[1],3) # AUC lower bound
out_v[3] <- round(ci_AUC[3],3) # AUC upper bound

outp <- paste0("AUC (CI): ", format(out_v[1],nsmall=3)," (",format(out_v[2],nsmall=3), "-",format(out_v[3],nsmall=3),")")
print(outp)
