# FNIH Project 2
# Function to make AUC Forest plots
#### Coded: Kellen Petersen
#### Date as of commenting and documentation - 05/31/2024

getAUCs_withLabels <- function(T,breaks1,indicatorline,rev_axis=FALSE,y1,y2){
  df <- T %>% 
    select(company,AUC,AUC_Lower,AUC_Upper,analyte)
  
  df$AUC <- as.numeric(df$AUC)
  df$AUC_Lower <- as.numeric(df$AUC_Lower)
  df$AUC_Upper <- as.numeric(df$AUC_Upper)
  DF <- df %>%
    mutate(index_0 = row_number()) %>%
    group_by(company) %>%
    mutate(index = row_number()) %>%
    ungroup() %>%
    mutate(plot_index = desc(row_number()))
  DF$analyte_c <- DF$analyte
  DF$analyte_c <- factor(DF$analyte_c,levels=c("p-tau217 ratio + Ab42/Ab40",
                                               "p-tau217 ratio",
                                               "p-tau217 + Ab42/Ab40",
                                               "p-tau217",
                                               "Ab42/Ab40",
                                               "p-tau181 + Ab42/Ab40 + GFAP + NfL",
                                               "p-tau181 + Ab42/Ab40 + NfL",
                                               "p-tau181 + Ab42/Ab40",
                                               "p-tau181",
                                               "GFAP",
                                               "NfL"))
  AB <- paste0("A","\U03B2","42/A","\U03B2","40")
  DF$analyte_l <- DF$analyte
  DF$analyte_l <- as.character(DF$analyte_l)
  DF[which(DF$analyte=="p-tau217 + Ab42/Ab40"),"analyte_l"] <- paste0("p-tau217 + ",AB)
  DF[which(DF$analyte=="p-tau217 ratio + Ab42/Ab40"),"analyte_l"] <- paste0("%p-tau217 + ",AB)
  DF[which(DF$analyte=="p-tau217 ratio"),"analyte_l"] <- paste0("%p-tau217")
  DF[which(DF$analyte=="Ab42/Ab40"),"analyte_l"] <- AB
  DF[which(DF$analyte=="p-tau181 + Ab42/Ab40"),"analyte_l"] <- paste0("p-tau181 + ",AB)
  DF[which(DF$analyte=="p-tau181 + Ab42/Ab40 + NfL"),"analyte_l"] <- paste0("p-tau181 + ",AB," + NfL")
  DF[which(DF$analyte=="p-tau181 + Ab42/Ab40 + GFAP + NfL"),"analyte_l"] <- paste0("p-tau181 + ",AB," + GFAP + NfL")
  
  analyte_colors <- c("p-tau217 ratio + Ab42/Ab40" = "black",
                      "p-tau217 + Ab42/Ab40" = "darkgray",
                      
                      "p-tau217 ratio" = "#145A32",
                      "p-tau217" = "#008000",
                      "p-tau181" = "#00FF00",
                      
                      "Ab42/Ab40" = "blue",
                      "NfL" = "maroon",
                      "GFAP" = "red", 
                      
                      "p-tau181 + Ab42/Ab40 + GFAP + NfL" = "purple",
                      "p-tau181 + Ab42/Ab40 + NfL" = "orange",
                      "p-tau181 + Ab42/Ab40" = "cyan"
  )
  
  # min_y <- round(min(DF$AUC_Lower),digits = 2) 
  min_y <- y1
  print(min_y)
  min_y_text <- min_y-0.1
  text_start <- min_y_text
  round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
  max_y <-  round_any(max(DF$AUC_Upper) ,accuracy = 0.2 , f = ceiling) 
  max_y <- max(c(max(breaks1),max_y))
  if(rev_axis==TRUE){
    max_y <- max_y+0.1
    text_start <- max_y
  }
  
  company_count <- table(DF$company)
  complab <- DF %>% 
    group_by(company) %>%
    slice(1) %>%
    ungroup() %>% 
    arrange(desc(plot_index)) 
  
  C2N_Index_Point <- as.numeric(complab[complab$company=="C2N","plot_index"])
  print(C2N_Index_Point)
  Fujirebio_Index_Point <- as.numeric(complab[complab$company=="Fujirebio","plot_index"])
  ALZPath_Index_Point <- as.numeric(complab[complab$company=="AlzPath","plot_index"])
  Janssen_Index_Point <- as.numeric(complab[complab$company=="Janssen","plot_index"])
  Roche_Index_Point <- as.numeric(complab[complab$company=="Roche","plot_index"])
  Quanterix_Index_Point <- as.numeric(complab[complab$company=="Quanterix","plot_index"])
  
  Indices <- complab$plot_index
  Endpoint <- as.numeric(Quanterix_Index_Point - company_count["Quanterix"])
  title_plot <- ""
  
  test_plot <-  ggplot(data=DF,
                                aes(x = plot_index, #uses index as the x-axis (flipped to y later)
                                    y = AUC, # y-axis is spearman rho, flipped to x-axis later
                                    ymin = min_y, ymax = max_y ))+ #limited of y-axis (flipped to x later)
    geom_rect(aes(xmin = Indices[1]+0.5, xmax = Indices[2]+0.5, ymin = min_y_text, ymax = max_y),
              fill = "#E5E4E2", alpha = 0.04) +
    geom_rect(aes(xmin = Indices[3]+0.5, xmax = Indices[4]+0.5, ymin = min_y_text, ymax = max_y),
              fill = "#E5E4E2", alpha = 0.04) +
    geom_rect(aes(xmin = Indices[5]+0.5, xmax = Indices[6]+0.5, ymin = min_y_text, ymax = max_y),
              fill = "#E5E4E2", alpha = 0.04) +
    geom_hline(aes(fill="black"),yintercept =indicatorline, linetype=2)+ #creates horizontal (flipped to vertical) line at 0
    geom_point(aes(col=analyte_c))+ #sets the coloring based on whether its the summary measure or not
    geom_errorbar(aes(ymin=(AUC_Lower), ymax=(AUC_Upper),col=analyte_c),width = 0, cex = 1,size=4)+ #formatting of forest lines
    geom_text(aes(x = C2N_Index_Point, y = text_start, label = "C2N", hjust = 0),size=3,family="Calibri") +
    geom_text(aes(x = Fujirebio_Index_Point, y = text_start, label = "Fujirebio", hjust = 0),size=3,family="Calibri") +
    geom_text(aes(x = ALZPath_Index_Point, y = text_start, label = "ALZpath", hjust = 0),size=3,family="Calibri") +
    geom_text(aes(x = Janssen_Index_Point, y = text_start, label = "Janssen", hjust = 0),size=3,family="Calibri") +
    geom_text(aes(x = Roche_Index_Point, y = text_start, label = "Roche", hjust = 0),size=3,family="Calibri") +
    geom_text(aes(x = Quanterix_Index_Point, y = text_start, label = "Quanterix", hjust = 0),size=3,family="Calibri")+
    labs(x = " ", 
         y = "AUC", #y axis will have no label, x axis will be spearman rho
         title = title_plot)+
    theme_classic()+ #gives the title based on the function call variable
    theme(plot.title=element_text(size=16,face="bold"), #these are all text formatting
          axis.ticks.y=element_line(size = 1, color = "black"),
          axis.text.x=element_text(size=7,face="bold"),
          axis.text.y=element_text(size=8),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_blank(),
          legend.position =  "none")+ #no legend
    scale_x_continuous(breaks=DF$plot_index,labels=DF$analyte_l)+
    scale_color_manual(values=analyte_colors)+
    scale_fill_manual(values = analyte_colors)+ #makes sure dots and lines are colored based on color code
    coord_flip() 
  
  test_plot <- test_plot +
    scale_y_continuous(limits = c(y1-.1,y2),
                       breaks = breaks1,
                       labels = scales::number_format(accuracy = 0.001))
  return(test_plot)
}








getAUCs_withoutLabels <- function(T,breaks1,indicatorline,rev_axis=FALSE,y1,y2){
  df <- T %>% 
    select(company,AUC,AUC_Lower,AUC_Upper,analyte)
  
  df$AUC <- as.numeric(df$AUC)
  df$AUC_Lower <- as.numeric(df$AUC_Lower)
  df$AUC_Upper <- as.numeric(df$AUC_Upper)
  DF <- df %>%
    mutate(index_0 = row_number()) %>%
    group_by(company) %>%
    mutate(index = row_number()) %>%
    ungroup() %>%
    mutate(plot_index = desc(row_number()))
  DF$analyte_c <- DF$analyte
  DF$analyte_c <- factor(DF$analyte_c,levels=c("p-tau217 ratio + Ab42/Ab40",
                                               "p-tau217 ratio",
                                               "p-tau217 + Ab42/Ab40",
                                               "p-tau217",
                                               "Ab42/Ab40",
                                               "p-tau181 + Ab42/Ab40 + GFAP + NfL",
                                               "p-tau181 + Ab42/Ab40 + NfL",
                                               "p-tau181 + Ab42/Ab40",
                                               "p-tau181",
                                               "GFAP",
                                               "NfL"))
  AB <- paste0("A","\U03B2","42/A","\U03B2","40")
  DF$analyte_l <- DF$analyte
  DF$analyte_l <- as.character(DF$analyte_l)
  DF[which(DF$analyte=="p-tau217 + Ab42/Ab40"),"analyte_l"] <- paste0("p-tau217 + ",AB)
  DF[which(DF$analyte=="p-tau217 ratio + Ab42/Ab40"),"analyte_l"] <- paste0("p-tau217 ratio + ",AB)
  DF[which(DF$analyte=="Ab42/Ab40"),"analyte_l"] <- AB
  DF[which(DF$analyte=="p-tau181 + Ab42/Ab40"),"analyte_l"] <- paste0("p-tau181 + ",AB)
  DF[which(DF$analyte=="p-tau181 + Ab42/Ab40 + NfL"),"analyte_l"] <- paste0("p-tau181 + ",AB," + NfL")
  DF[which(DF$analyte=="p-tau181 + Ab42/Ab40 + GFAP + NfL"),"analyte_l"] <- paste0("p-tau181 + ",AB," + GFAP + NfL")
  
  analyte_colors <- c("p-tau217 ratio + Ab42/Ab40" = "black",
                      "p-tau217 + Ab42/Ab40" = "darkgray",
                      
                      "p-tau217 ratio" = "#145A32",
                      "p-tau217" = "#008000",
                      "p-tau181" = "#00FF00",
                      
                      "Ab42/Ab40" = "blue",
                      "NfL" = "maroon",
                      "GFAP" = "red", 
                      
                      "p-tau181 + Ab42/Ab40 + GFAP + NfL" = "purple",
                      "p-tau181 + Ab42/Ab40 + NfL" = "orange",
                      "p-tau181 + Ab42/Ab40" = "cyan"
  )
  
  
  min_y <- round(min(DF$AUC_Lower),digits = 2) 
  min_y <- y1
  min_y_text <- min_y-0
  text_start <- min_y_text
  round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
  max_y <-  round_any(max(DF$AUC_Upper) ,accuracy = 0.05 , f = ceiling) 
  print(max_y)
  max_y <- 0
  max_y <- max(c(max(breaks1),max_y))
  print(max_y)
  if(rev_axis==TRUE){
    max_y <- max_y+0.1
    text_start <- max_y
  }
  
  company_count <- table(DF$company)
  complab <- DF %>% 
    group_by(company) %>%
    slice(1) %>%
    ungroup() %>% 
    arrange(desc(plot_index)) 
  
  C2N_Index_Point <- as.numeric(complab[complab$company=="C2N","plot_index"])
  Fujirebio_Index_Point <- as.numeric(complab[complab$company=="Fujirebio","plot_index"])
  ALZPath_Index_Point <- as.numeric(complab[complab$company=="AlzPath","plot_index"])
  Janssen_Index_Point <- as.numeric(complab[complab$company=="Janssen","plot_index"])
  Roche_Index_Point <- as.numeric(complab[complab$company=="Roche","plot_index"])
  Quanterix_Index_Point <- as.numeric(complab[complab$company=="Quanterix","plot_index"])
  
  Indices <- complab$plot_index
  Endpoint <- as.numeric(Quanterix_Index_Point - company_count["Quanterix"])
  title_plot <- ""
  
  test_plot <-  ggplot(data=DF,
                                aes(x = plot_index, #uses index as the x-axis (flipped to y later)
                                    y = AUC, # y-axis is spearman rho, flipped to x-axis later
                                    ymin = min_y, ymax = max_y ))+ #limited of y-axis (flipped to x later)
    geom_rect(aes(xmin = Indices[1]+0.5, xmax = Indices[2]+0.5, ymin = min_y_text, ymax = max_y),
              fill = "#E5E4E2", alpha = 0.04) +
    geom_rect(aes(xmin = Indices[3]+0.5, xmax = Indices[4]+0.5, ymin = min_y_text, ymax = max_y),
              fill = "#E5E4E2", alpha = 0.04) +
    geom_rect(aes(xmin = Indices[5]+0.5, xmax = Indices[6]+0.5, ymin = min_y_text, ymax = max_y),
              fill = "#E5E4E2", alpha = 0.04) +
    geom_hline(aes(fill="black"),yintercept =indicatorline, linetype=2)+ #creates horizontal (flipped to vertical) line at 0
    geom_point(aes(col=analyte_c))+ #sets the coloring based on whether its the summary measure or not
    geom_errorbar(aes(ymin=(AUC_Lower), ymax=(AUC_Upper),col=analyte_c),width = 0, cex = 1,size=4)+ #formatting of forest lines
    labs(x = " ", 
         y = "AUC", #y axis will have no label, x axis will be spearman rho
         title = title_plot)+
    theme_classic()+ #gives the title based on the function call variable
    theme(plot.title=element_text(size=16,face="bold"), #these are all text formatting
          axis.ticks.y=element_line(size = 1, color = "black"),
          axis.text.x=element_text(size=7,face="bold"),
          axis.text.y=element_text(size=8),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_blank(),
          legend.position =  "none")+ #no legend
    scale_x_continuous(breaks=DF$plot_index,labels=DF$analyte_l)+
    scale_color_manual(values=analyte_colors)+
    scale_fill_manual(values = analyte_colors)+ #makes sure dots and lines are colored based on color code
    theme(axis.text.y=element_blank(), # Remove x axis labels
          axis.ticks.y=element_blank()) + # Remove x axis ticks
    coord_flip() #+#+ 
  
  test_plot <- test_plot +
    scale_y_continuous(limits = c(y1,y2),
                       breaks = breaks1,
                       labels = scales::number_format(accuracy = 0.001))
  
  return(test_plot)
}






make_model <- function(data.in0,covars_TF,vars_outcome){
  data.in.5 <- data.in0
  #Covariates
  age_var <- "AGE" #variable with age
  sex_var <- "PTGENDER" #variable with sex
  educ_var <- "PTEDUCAT" #variable with education
  
  APOE_variable <- "APOE_genotype" #Variable with APOE genotypes
  data.in.5$APOE_Cat <- factor(data.in.5[[APOE_variable]], levels = c("33", "34", "24", "22", "23", "44"))
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
  A_qcsfratio <- c("CSF PTAU/AB42")
  A <- c(A_c2n, A_fuji, A_alzpath, A_janssen, A_roche, A_quanterix)
  AA <- A

  single_analytes <- c("p-tau217 ratio",
                       "Ab42/Ab40",
                       "p-tau217",
                       "p-tau181",
                       "GFAP", 
                       "NfL",
                       "CSF PTAU/AB42")

  OUT_5 <- list()
  MODELS_5 <- list()
  k <- 0
  data.in <- data.in.5
  # Loop through the combinations of predictors
  for (i_outcome in 1:length(vars_outcome)) {
    for (i_assay_combos in 1:length(Combos)) {
      
      outcome <- vars_outcome[i_outcome]
      pred1 <- Combos[[i_assay_combos]]
      pred2 <- vars_covariates
      preds <- c(unlist(pred1), unlist(pred2))
      
      print(i_assay_combos)
      print(preds)
      
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
      out_v[5] <- AA[i_assay_combos] # analyte label
      
      model <- list()
      model[[1]] <- mdl_logistic
      model[[2]] <- roc_tmp
      
      k <- i_assay_combos + length(combos_all)*(i_outcome-1)
      print(i_assay_combos)
      print(i_outcome)
      print(k)
      print(optimal_cutoff)
      OUT_5[[k]] <- out_v
      MODELS_5[[k]] <- model
    }
  }
  
  df_out <- t(data.frame(OUT_5))
  df_out <- data.frame(df_out)
  rownames(df_out) <- NULL
  colnames(df_out) <- c("company","AUC","AUC_Lower","AUC_Upper","analyte")
  
  df_out$AUC <- as.numeric(df_out$AUC)
  df_out$AUC_Lower <- as.numeric(df_out$AUC_Lower)
  df_out$AUC_Upper <- as.numeric(df_out$AUC_Upper)
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
  
  T <- DFB %>% 
    select(company,AUC,AUC_Lower,AUC_Upper,analyte) %>%
    group_by(company) %>%
    mutate(max_AUC = max(AUC)) %>%
    ungroup() %>%
    arrange(desc(max_AUC), company, desc(AUC)) %>% 
    select(-max_AUC) %>% 
    filter(company != "CSF")
  
  return(T)
}