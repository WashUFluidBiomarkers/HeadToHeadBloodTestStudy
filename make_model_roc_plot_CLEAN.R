# FNIH Project 2
# make models for ROC plots
#### Coded: Kellen Petersen
#### Date as of commenting and documentation - 05/28/2024

make_model <- function(data.in0,covars_TF,vars_outcome){
  data.in <- data.in0
  
  #Covariates
  age_var <- "AGE" #variable with age
  sex_var <- "PTGENDER" #variable with sex
  educ_var <- "PTEDUCAT" #variable with education
  
  APOE_variable <- "APOE_genotype" #Variable with APOE genotypes
  data.in$APOE_Cat <- factor(data.in[[APOE_variable]], levels = c("33", "34", "24", "22", "23", "44"))
  list_covariates <- list(c(age_var, sex_var, "APOE_Cat"))
  
  if (covars_TF) {
    vars_covariates <- list_covariates
  } else {
    vars_covariates <- NULL
  }

  # Assays Lists for each platform
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
  A_5csfratio <- c("CSF PTAU/AB42")
  A_quanterix <- c("p-tau181 + Ab42/Ab40 + GFAP + NfL",
                   "p-tau181 + Ab42/Ab40 + NfL",
                   "p-tau181 + Ab42/Ab40",
                   "p-tau181","Ab42/Ab40",
                   "GFAP","NfL")
  A_qcsfratio <- c("CSF PTAU/AB42")
  
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
  
  
  # Run the models
  OUT_df <- list()
  MODELS_df <- list()
  k <- 0
  data.in <- data.in
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
      brierScore <- mean((mdl_pred-temp_data[[outcome]])^2)/dim(temp_data)[1]
      
      roc_tmp <- pROC::roc(temp_data[[outcome]],mdl_pred,plot=FALSE,print.auc=FALSE) #this will get the ROC output
      
      # YOUDEN INDEX
      youden_index <- pROC::coords(roc_tmp, "best", ret="all")
      optimal_cutoff <- youden_index$threshold

      T <- data.frame(roc_tmp$thresholds, roc_tmp$sensitivities, roc_tmp$specificities)
      auc_tmp <- pROC::auc(roc_tmp)
      ci_AUC  <- ci.auc(roc_tmp) #get confidence interview 
      
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
      out_v[2] <- round(roc_tmp$auc,3) # AUC rounded to 2 decimals
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
      out_v[22] <- brierScore
      
      model <- list()
      model[[1]] <- mdl_logistic
      model[[2]] <- roc_tmp
      
      k <- i_assay_combos + length(combos_all)*(i_outcome-1)

      OUT_df[[k]] <- out_v
      MODELS_df[[k]] <- model
    }
    
    # Run a single time for the covariates only
    outcome <- vars_outcome[i_outcome]
    pred2 <- vars_covariates
    preds <- c(unlist(list_covariates))
    
    formula_str <- paste(outcome, " ~ ", paste(preds, collapse = " + "))
    model_formula <- as.formula(formula_str)
    
    temp_data <- data.in[complete.cases(data.in[,preds]), c(outcome,preds)] %>% 
      na.omit()
    
    mdl_logistic <- glm(model_formula, data = temp_data, family = "binomial")
    mdl_pred <- predict(mdl_logistic, data = temp_data, type = "response")
    brierScore <- mean((mdl_pred-temp_data[[outcome]])^2)
    
    roc_tmp <- pROC::roc(temp_data[[outcome]],mdl_pred,plot=FALSE,print.auc=FALSE) #this will get the ROC output
    
    youden_index <- pROC::coords(roc_tmp, "best", ret="all")
    optimal_cutoff <- youden_index$threshold
    
    T <- data.frame(roc_tmp$thresholds, roc_tmp$sensitivities, roc_tmp$specificities)
    auc_tmp <- pROC::auc(roc_tmp)
    ci_AUC  <- ci.auc(roc_tmp) #get confidence intervals
    
    sensitivity <- youden_index$sensitivity
    specificity <- youden_index$specificity
    
    # Sensitivity = 95%
    II_sensX <- which(roc_tmp$sensitivities >= .90)
    I_sensX <- max(II_sensX)
    sensX_threshold <- roc_tmp$thresholds[I_sensX]
    sensX_specificity <- roc_tmp$specificities[I_sensX]
    
    # Specificity = 95%
    II_specX <- which(roc_tmp$specificities >= .90)
    I_specX <- min(II_specX)
    specX_threshold <- roc_tmp$thresholds[I_specX]
    specX_sensitivities <- roc_tmp$sensitivities[I_specX]
    
    num_intermediate <- sum( (mdl_pred > sensX_threshold) & (mdl_pred < specX_threshold) )
    num_intermediate
    
    out_v <- c()
    out_v[1] <- "covariates" # company label
    out_v[2] <- round(roc_tmp$auc,3) # AUC rounded to 2 decimals
    out_v[3] <- round(ci_AUC[1],3) # AUC lower bound
    out_v[4] <- round(ci_AUC[3],3) # AUC upper bound
    out_v[5] <- paste0(formatC(sensitivity, digits=3, format="fg"))
    out_v[6] <- paste0(formatC(specificity, digits=3, format="fg"))
    out_v[7] <- paste0(formatC(youden_index$ppv, digits=3, format="fg"))
    out_v[8] <- paste0(formatC(youden_index$npv, digits=3, format="fg"))
    out_v[9] <- paste0(formatC(youden_index$accuracy, digits=3, format="fg"))
    
    out_v[10] <- -99
    
    out_v[11] <- "covariates" # analyte label

    out_v[12] <- formatC(sensX_threshold, digits=3, format="fg")
    out_v[13] <- formatC(sensX_specificity, digits=3, format="fg")
    
    out_v[14] <- formatC(specX_threshold, digits=3, format="fg")
    out_v[15] <- formatC(specX_sensitivities, digits=3, format="fg")
    
    out_v[16] <- num_intermediate
    
    out_v[17] <- formula_str
    out_v[18] <- outcome
    out_v[19] <- " "
    out_v[20] <- covars_TF
    out_v[21] <- dd[1]
    out_v[22] <- brierScore
    
    model <- list()
    model[[1]] <- mdl_logistic
    model[[2]] <- roc_tmp
    
    k <- 1 + length(combos_all) #+ length(combos_all)*(i_outcome-1)
    
    OUT_df[[k]] <- out_v
    MODELS_df[[k]] <- model
  }
  
  # Make output into dataframe
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
  DFB$company <- factor(DFB$company, levels = c("C2N", "Fujirebio", "AlzPath","Janssen", "Roche", "Quanterix","CSF","covariates"))
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
                                   "CSF PTAU/AB42",
                                   "covariates"))
  DFB2 <- DFB %>%
    mutate(AUC_CI = paste0(format(AUC,nsmall=3)," (",format(AUC_Lower,nsmall=3), "-",format(AUC_Upper,nsmall=3),")")) %>%
    select("company","analyte",
           "AUC_CI",
           "sensitivity","specificity",
           "90% Sens Threshold","90% Sens - Specificity",
           "90% Spec Threshold","90% Spec - Sensitivity",
           "intermediate_results")
  
  # Final output dataframe T
  T <- DFB %>% 
    select(company,AUC,AUC_Lower,AUC_Upper,analyte,Brier_Score) %>%
    mutate(original_order = row_number()) %>% 
    group_by(company) %>%
    mutate(max_AUC = max(AUC)) %>%
    ungroup() %>%
    arrange(desc(max_AUC), company, desc(AUC)) %>% 
    select(-max_AUC) %>% 
    filter(company != "CSF")
  
  REF <- data.frame(company = C, analyte = A) %>% 
    mutate(index = row_number()) %>% 
    filter(company != "CSF")
  
  return(list(T = T, MODELS_df = MODELS_df, REF = REF))
}
