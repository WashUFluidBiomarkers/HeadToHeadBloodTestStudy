#Function Bank and Variable Definitions for FNIH Project 2
#Crossectional Part of Study
#### Coded: Benjamin Saef
#### Date as of commenting and documentation - 05/28/2024
library(tidyr)
library(data.table)
library(RVAideMemoire)
library(pROC)
library(RColorBrewer)
library(MuMIn)

####### Variable creation - centralized location, just so that it can be made sure of across all scripts

#analysis vars
#Covariates
sex_var <- "PTGENDER" #variable with sex
age_var <- "AGE" #variable with age
APOE_variable <- "APOE_genotype" #Variable with APOE genotypes

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






#Function if needed to bring p-values in line with how they are outputted in SAS

pub_pval <- function(pval){
  
  if(pval < 0.0001){ pubpval <- "<0.0001" 
  } else if(pval > 0.001){
    pubpval <- paste(signif(pval,2))
  }else{ pubpval <- paste(signif(pval,4))}
  return(pubpval)
}
#Function to Compare Rho
spdiff_rho <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  cor_xz <- cor(d$x, d$z, method = "spearman",use = "complete.obs")
  cor_yz <- cor(d$y, d$z, method = "spearman",use = "complete.obs")
  diff <- abs(cor_xz) - abs(cor_yz)
  
  return(diff)
}
#Function to compare partial rhos
#Compare Rhos - needs fixing dependent on age, sex and APOE variables
#the covariate variables are currently hard coded because out how the bootstrap function works


spdiff_Partial_rho <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  c_cov <- c("AGE","PTGENDER","APOE_Cat")   #change to reflect data
  cor_xz <- RVAideMemoire::pcor.test(d$x, d$z,d[,c_cov], method = "spearman",nrep = 1)
  cor_yz <- RVAideMemoire::pcor.test(d$y, d$z, d[,c_cov],method = "spearman",nrep = 1)
  diff <- abs(cor_xz$estimate) - abs(cor_yz$estimate)
  
  return(diff)
}
#Functions to create Characteristics Tables

###Creates entry focused around how many participants are positive for a specific variable (or just which one has higher numerical code that would go 
#second in a table)
positivity_dist_entry <- function(testvar,
                                  binvar,
                                  data_pos,
                                  data_neg,
                                  dataset,
                                  entryname){
  
  full_dist <- paste0(table(unlist(dataset[,testvar]))[2],
                      ", ",
                      round(100*table(unlist(dataset[,testvar]))[2]/length(unlist(dataset[,testvar])),2),
                      "%")
  pos_dist <- paste0(table(unlist(data_pos[,testvar]))[2],
                     ", ",
                     round(100*table(unlist(data_pos[,testvar]))[2]/length(unlist(data_pos[,testvar])),2),
                     "%")
  neg_dist <- paste0(table(unlist(data_neg[,testvar]))[2],
                     ", ",
                     round(100*table(unlist(data_neg[,testvar]))[2]/length(unlist(data_neg[,testvar])),2),
                     "%")
  
  chisq_test_Res <- chisq.test(unlist(dataset[,binvar]), 
                               unlist(dataset[,testvar]), 
                               correct=FALSE)
  pval_pre <- chisq_test_Res$p.value
  if(pval_pre < 0.0001){ pval <- "<0.0001" 
  } else{ pval <- paste(signif(pval_pre,1))}
  
  
  posit_dist <-  c(entryname,
                   sum(!is.na(dataset[,testvar])),
                   full_dist,
                   sum(!is.na(data_neg[,testvar])),
                   neg_dist,
                   sum(!is.na(data_pos[,testvar])),
                   pos_dist,
                   pval)
  
  return(posit_dist)
  
}




####creates an entry focused around the median and first and third quantiles for continuous data (usually non-normal)

create_analyte_dist_entry <-  function(Varid,binvar,g1_dataset,g2_dataset,full_dataset,labelid){
  g1_dataset <- g1_dataset
  g2_dataset <- g2_dataset
  full_dataset <- full_dataset
  pos_dat <- g1_dataset[which(!is.na(g1_dataset[,Varid])),]
  neg_dat <- g2_dataset[which(!is.na(g2_dataset[,Varid])),]
  all_dat <- full_dataset[which(!is.na(full_dataset[,Varid])),]
  all_sum <- signif(quantile(all_dat[,Varid],type=3), digits = 3)
  pos_sum <- signif(quantile(pos_dat[,Varid],type=3), digits = 3)
  neg_sum <- signif(quantile(neg_dat[,Varid],type=3), digits = 3)
  formula_of_model <- paste(Varid,binvar,sep="~")  
  test_data <-  kruskal.test(as.formula(formula_of_model),data=all_dat)
  pval_pre <- test_data$p.value
  if(pval_pre < 0.0001){ pval <- "<0.0001" 
  } else{ pval <- paste(signif(pval_pre,4))}
  entry_dist <- c(labelid,
                  sum(!is.na(all_dat[,Varid])),
                  paste0(all_sum["50%"]," (",
                         all_sum["25%"],
                         "-",all_sum["75%"],")"),
                  sum(!is.na(neg_dat[,Varid])),
                  paste0(neg_sum["50%"]," (",
                         neg_sum["25%"],
                         "-",neg_sum["75%"],")"),
                  sum(!is.na(pos_dat[,Varid])),
                  paste0(pos_sum["50%"]," (",
                         pos_sum["25%"],
                         "-",pos_sum["75%"],")"),
                  pval)
  
  return(entry_dist)
}

###Creates an entry focused around the means and standard deviation for the continuous data

create_mean_entry <-  function(Varid,binvar,g1_dataset,g2_dataset,full_dataset,labelid){
  
  pos_dat <- g1_dataset[which(!is.na(g1_dataset[,Varid])),]
  neg_dat <- g2_dataset[which(!is.na(g2_dataset[,Varid])),]
  all_dat <-full_dataset[which(!is.na(full_dataset[,Varid])),]
  all_sum <- signif(mean(all_dat[,Varid],na.rm=T), digits = 3)
  pos_sum <- signif(mean(pos_dat[,Varid],na.rm=T), digits = 3)
  neg_sum <- signif(mean(neg_dat[,Varid],na.rm=T), digits = 3)
  all_sd <- signif(sd(all_dat[,Varid],na.rm=T), digits = 3)
  pos_sd <- signif(sd(pos_dat[,Varid],na.rm=T), digits = 3)
  neg_sd <- signif(sd(neg_dat[,Varid],na.rm=T), digits = 3)
  
  formula_of_model <- paste(Varid,binvar,sep="~")  
  test_data <-  kruskal.test(as.formula(formula_of_model),data=all_dat)
  pval_pre <- test_data$p.value
  if(pval_pre < 0.0001){ pval <- "<0.0001" 
  } else{ pval <- paste(signif(pval_pre,4))}
  entry_dist <- c(labelid,
                  sum(!is.na(all_dat[,Varid])),
                  paste0(all_sum," (",
                         all_sd,")"),
                  sum(!is.na(neg_dat[,Varid])),
                  paste0(neg_sum," (",
                         neg_sd,")"),
                  sum(!is.na(pos_dat[,Varid])),
                  paste0(pos_sum," (",
                         pos_sd,")"),
                  pval)
  
  return(entry_dist)
}

###The following two are circumstantial and only to be used for their exact purpose. They may require re-coding if things are differently 
#coded in the target dataset

#Race (Black/White/Other) - need to be coded as 2/1/3 respectively
create_race_entry <- function(Varid,binvar,g1_dataset,g2_dataset,full_dataset,labelid){

  pos_dat <- g1_dataset[which(!is.na(g1_dataset[,Varid])),]
  neg_dat <- g2_dataset[which(!is.na(g2_dataset[,Varid])),]
  all_dat <-full_dataset[which(!is.na(full_dataset[,Varid])),]
  
  
full_dist <- paste0(table(all_dat[,Varid])["2"],"/",
                         table(all_dat[,Varid])["1"],"/",
                         table(all_dat[,Varid])["3"])
neg_dist <- paste0(table(neg_dat[,Varid])["2"],"/",
                           table(neg_dat[,Varid])["1"],"/",
                           table(neg_dat[,Varid])["3"])
pos_dist <- paste0(table(pos_dat[,Varid])["2"],"/",
                           table(pos_dat[,Varid])["1"],"/",
                           table(pos_dat[,Varid])["3"])

chisq_race <- chisq.test(all_dat[,binvar], 
                         as.factor(all_dat[,Varid]), 
                         correct=FALSE)

racial_characteristcs <- c("Race (Black/White/Other)",
                           sum(!is.na(all_dat[,Varid])),
                           full_dist,
                           sum(!is.na(neg_dat[,Varid])),
                           neg_dist,
                           sum(!is.na(pos_dat[,Varid])),
                           pos_dist,
                           signif(chisq_race$p.value,2))

return(racial_characteristcs)

}


###Creating this just in case we decide to do APOE in a different way this time
create_apoe_entry <- function(Varid,binvar,g1_dataset,g2_dataset,full_dataset,labelid){
  
  pos_dat <- g1_dataset[which(!is.na(g1_dataset[,Varid])),]
  neg_dat <- g2_dataset[which(!is.na(g2_dataset[,Varid])),]
  all_dat <-full_dataset[which(!is.na(full_dataset[,Varid])),]
  
  
  full_dist <- paste0(table(all_dat[,Varid])["22"],"/",
                      table(all_dat[,Varid])["23"],"/",
                      table(all_dat[,Varid])["24"],"/",
                      table(all_dat[,Varid])["33"],"/",
                      table(all_dat[,Varid])["34"],"/",
                      table(all_dat[,Varid])["44"])
  neg_dist <- paste0(table(neg_dat[,Varid])["22"],"/",
                     table(neg_dat[,Varid])["23"],"/",
                     table(neg_dat[,Varid])["24"],"/",
                     table(neg_dat[,Varid])["33"],"/",
                     table(neg_dat[,Varid])["34"],"/",
                     table(neg_dat[,Varid])["44"])
  pos_dist <- paste0(table(pos_dat[,Varid])["22"],"/",
                     table(pos_dat[,Varid])["23"],"/",
                     table(pos_dat[,Varid])["24"],"/",
                     table(pos_dat[,Varid])["33"],"/",
                     table(pos_dat[,Varid])["34"],"/",
                     table(pos_dat[,Varid])["44"])
  
  chisq_apoe <- chisq.test(all_dat[,binvar], 
                           as.factor(all_dat[,Varid]), 
                           correct=FALSE)
  
  apoe_characteristcs <- c("APOE Genotype (22/23/24/33/34/44)",
                             sum(!is.na(all_dat[,Varid])),
                             full_dist,
                             sum(!is.na(neg_dat[,Varid])),
                             neg_dist,
                             sum(!is.na(pos_dat[,Varid])),
                             pos_dist,
                             signif(chisq_apoe$p.value,2))
  
  return(apoe_characteristcs)
  
}




#### Pulling all functions together to construct tables

#dataset - input data - needs to contain all variables in binvar and varlist
#binvar and binvarlab - the binary variable which will determine the columns after the "all" columns and the label will determine
#how the columns are named. both of these need to be strings
#varlist - list of variables, simple as it sounds, it needs to be a vector of strings
#typelist - here's where things get a little complicated. this needs to be the exact same length as the varlist, as it's what
#determines how each entry is handled
#Here are the type listings
#A - this is the code for any analyte/non-normal variable. It will output an entry with the median and interquartile range
#T - this code is very similar to the previous one, it's only use specifically for Tracers when dealing with different measurement standards for
#SUVR, the namelist entry, where the labels are, will be used to determine which tracer this entry will represent 
#P - this is the code that is used to create entries based on binary variables. Specifically the "positive" status that is measured
#is whatever R would place in the second slot of the variable when placed into the table function
#This is usually used for things like disease or biomarker status, but can be used for other things as well
#namelist - list of labels that go along with varlist and typelist, needs to be in order corresponding to those two
#filout - output destination for resulting excel table
char_table_full_bin <- function(dataset,binvar,
                                binvarlab,varlist,
                                typelist,namelist,
                                fileout){
  require(openxlsx)
  require(tidyverse)
  looplength <- length(varlist)
  dnames <- c("All",paste0(binvarlab," Negative"),
              paste0(binvarlab," Positive"))
  
  d1 = dataset[which(dataset[,binvar]==1),]
  d2 = dataset[which(dataset[,binvar]==0),]
  characteristics_table <- data.frame(matrix(nrow=looplength,ncol=(length(dnames)*2)+2))
  colnames(characteristics_table) <- c("Characteristic",
                                       paste(rep(dnames,each=2),
                                             c("n","Values"),sep=" "),"p=")
  
  for(i in 1:looplength){
    vartype <- typelist[i]
    varid <- varlist[i]
    label_name <- namelist[i]
    if(vartype=="A"){
      newentry  <- create_analyte_dist_entry(varid,
                                             binvar,
                                             d1,
                                             d2,
                                             dataset,
                                             label_name)
    }else if(vartype=="P"){
      newentry  <- positivity_dist_entry(varid,
                                         binvar,
                                         d1,
                                         d2,
                                         dataset,
                                         label_name)
    } 
    
    
    else if(vartype=="T"){
      newentry <- create_mean_entry(varid,
                                    binvar,
                                    d1[which(d1$Tracer==label_name),],
                                    d2[which(d2$Tracer==label_name),],
                                    dataset[which(dataset$Tracer==label_name),],
                                    paste0(label_name," SUVR"))
    } else if(vartype=="M"){
      newentry  <- create_mean_entry(varid,
                                     binvar,
                                     d1,
                                     d2,
                                     dataset,
                                     label_name)
    } else if(vartype=="R"){
      newentry  <- create_race_entry(varid,
                        binvar,
                        d1,
                        d2,
                        dataset,
                        label_name)
    } else if(vartype=="4"){
      newentry  <-create_apoe_entry(varid,
                        binvar,
                        d1,
                        d2,
                        dataset,
                        label_name)
    }
    characteristics_table[i,] <- newentry   
  }
  
  write.xlsx(characteristics_table,
             fileout,
             rowNames=FALSE,
             overwrite = TRUE)
  return(characteristics_table)
  
}







###Functions to run spearmans and partial spearmans 


#Function to evaluate the comparison standard for the Spearman comparisons
#Just outputs the top analyte variable name
Eval_comp_standard <- function(data_input, #Dataset
                               var_a, #Phenotype - outcome variable
                               vars_list){ #list of analytes
  
  length_analyte <- length(vars_list)
  out_res <- as.data.frame(matrix(nrow=length_analyte,ncol=2))
  colnames(out_res) <- c("varname","rho")
  
  for(i_analyte in 1:length_analyte){
    
    mdl_rho=cor.test(data_input[,vars_list[i_analyte]],
                     data_input[,var_a], method = "spearman")
    out_res$varname[i_analyte]=vars_list[i_analyte]
    out_res$rho[i_analyte]=mdl_rho$estimate
    
  }
  
  out_res <- out_res[order(-abs(out_res$rho)),]
  TopVar <- out_res[1,"varname"]
  return(TopVar)
}

####
#Function to evaluate the comparison standard for the Partial Spearman comparisons
#Just outputs the top analyte variable name
Eval_comp_standard_partial <- function(data_input,#Dataset
                               var_a, #Phenotype - outcome variable
                               cov_var, #Covariates - Age, Gender, APOE genotype
                               vars_list){ #list of analytes
  df <- data_input
  length_analyte <- length(vars_list)
  out_res <- as.data.frame(matrix(nrow=length_analyte,ncol=2))
  colnames(out_res) <- c("varname","rho")
  
  for(i_analyte in 1:length_analyte){
    df2 <- df[which(complete.cases(df[,c(var_a,vars_list[i_analyte],cov_var)])),]
    mdl_rho=RVAideMemoire::pcor.test(df2[,var_a],
                                     df2[,vars_list[i_analyte]], 
                                     df2[,cov_var],method = "spearman",
                                     conf.level = 0.95, nrep = 2)
    out_res$varname[i_analyte]=vars_list[i_analyte]
    out_res$rho[i_analyte]=mdl_rho$estimate
    
  }
  
  out_res <- out_res[order(-abs(out_res$rho)),]
  TopVar <- out_res[1,"varname"]
  return(TopVar)
}



###
spearman_plasma_runs <- function(data_input, # Data inputted 
                                 var_a,  # variable A - continuous variable being correlated
                                 vars_to_cor,# vector of strings - continuous variables being correlated with variable A
                                 var_comp, #variable in vars_to_cor vector that is used as standard for comparisons (unadjusted)
                                 var_comp_p, #variable in vars_to_cor vector that is used as standard for comparisons (adjusted)
                                 cov_var,    # covariates - variables used as covariates (age, sex, APOE normally)
                                 analytes_labels, #labels for analytes in vars to cor
                                 group_lab){ #groups for analytes in vars to cor
  require(RVAideMemoire)
  require(boot)
 df  <- data_input
 print_headers <- c("Biomarker","group","Rho","Rho_L","Rho_H","p",
                    paste("pComp",var_comp,sep="_"),"Partial_Rho","PRho_L","PRho_H","p_partial",paste("Par_pComp",var_comp_p,sep="_"))
 results_print=data.frame(matrix(nrow = length(vars_to_cor),
                                 ncol = length(print_headers)))
 colnames(results_print) <- print_headers
 
 
 
 
 for(i_plasma in 1:length(vars_to_cor)){
   results_print$Biomarker[i_plasma]=analytes_labels[i_plasma]
   results_print$group[i_plasma]=group_lab[i_plasma]
   ###Spearman (unadjusted) Correlation Statistics below
   mdl_rho=cor.test(df[,vars_to_cor[i_plasma]],
                    df[,var_a], method = "spearman")
   b=RVAideMemoire::spearman.ci(df[,vars_to_cor[i_plasma]], 
                                df[,var_a],
                                nrep = 1000,conf.level = 0.95)
 
   df2=data_input[,c(var_a,var_comp,vars_to_cor[i_plasma])]
   colnames(df2)=c("z","x","y")
   
   boot_results = boot(data = df2, statistic = spdiff_rho, R = 12500)
   results_under_H0 <- boot_results$t - mean(boot_results$t)
   boot_pvalue <- mean(abs(results_under_H0) >= abs(boot_results$t0))
   results_print[i_plasma,paste("pComp",var_comp,sep="_")]=pub_pval(boot_pvalue)
   

   results_print$Rho[i_plasma]=mdl_rho$estimate
   results_print$Rho_L[i_plasma]=b$conf.int[1]
   results_print$Rho_H[i_plasma]=b$conf.int[2]  
   results_print$p[i_plasma]=paste(pub_pval(mdl_rho$p.value))

   ###Partial Correlation Statistics below
   df_for_cor=df[complete.cases(df[,c(var_a,vars_to_cor[i_plasma],cov_var)]),]
 mdl_rho=RVAideMemoire::pcor.test(df_for_cor[,var_a],
                                  df_for_cor[,vars_to_cor[i_plasma]], 
                                  df_for_cor[,cov_var],method = "spearman",
                                  conf.level = 0.95, nrep = 1000)

 df2=df_for_cor[,c(var_a,var_comp_p,vars_to_cor[i_plasma],cov_var)]
 colnames(df2)=c("z","x","y",cov_var)
 
 set.seed(12345)
 boot_results = boot(data = df2, statistic = spdiff_Partial_rho, R = 12500)
 boot_ci_norm <- boot.ci(boot_results, type = "norm")
 results_under_H0 <- boot_results$t - mean(boot_results$t)
 boot_pvalue <- mean(abs(results_under_H0) >= abs(boot_results$t0))
 

 results_print$Partial_Rho[i_plasma]=mdl_rho$estimate
 results_print$PRho_L[i_plasma]=mdl_rho$conf.int[1]
 results_print$PRho_H[i_plasma]=mdl_rho$conf.int[2]
 results_print$p_partial[i_plasma]=paste(pub_pval(mdl_rho$p.value))
 results_print[i_plasma,paste("Par_pComp",var_comp_p,sep="_")]=pub_pval(boot_pvalue)
 }
 results_print <- results_print[order(results_print$group,-abs(results_print$Rho)),]
 return(results_print)
 
}

#The resulting dataset will include the following columns:
#Biomarker - Biomarker Name (from the analytes_labels vector)
#Rho - Spearman correlation statistic
#Rho_L - spearman rho - low confidence interval bound
#Rho_H - spearman rho - high confidence interval bound
#p - spearman correlation p-value
#pComp_'Analyte' - comparison p-value for the spearman correlations to the lead 'Analyte' chosen at the top of the function
#Partial_Rho - Partial Spearman's Rho
#PRho_L - Partial Rho - Low confidence interval bound
#PRho_H - Partial Rho - High confidence interval bound
#p_partial - partial correlation p-value
#Par_pComp_'Analyte' - comparison p-value for the partial correlations to the lead 'Analyte' chosen at the top of the function


###Functions to create forest plots

#creating function to create index. this will be output
#the main output will be the Analyte variable, which will be used as a way to merge with the dataset
#that's about to be plotted
#The summary indication variable "Analyte.group"
#The group_factor variable, which will be used to order the groups themselves (intergroup, not intragroup)
#The index variable, which will be used to determine the ordering within each group
## The functionality of these has the potential to fall apart if the Analyte variables between
## Datasets are not formatted the same.
#I'm also super hoping data availability between analytes doesn't become a super large issue
#Don't worry about the variable being labeled as AUC, you'll be able to label the axis in the final product
index_biomarker_forest <- function(DF1){
  
  require(dplyr)
  df <- data.frame(DF1[,c("company","group","Rho","Biomarker")])
  rownames(df) <- NULL
  df$Rho <- as.numeric(df$Rho)
  df$company <- factor(df$company)
  df$company <-factor(df$company, levels = c(unique(df[order(-abs(df$Rho)),"company"])))
  DF <- df %>% group_by(group) %>%
    arrange(company,desc(abs(Rho))) %>% 
    ungroup() %>% 
    mutate(plot_index = desc(row_number()))

  return(DF[,c("Biomarker","company","plot_index")])
}


####Create Output forest plot from spearman function
##As of right now, is company agnostic, but not analyte agnostic due to hard coded analyte labels

spearman_gg_formanuscript <- function(DF1, #output from spearman_plasma_runs function
                                      breaks1, # a vector representing the axis ticks
                                      indicatorline, # position dashed line that goes on the horizontal axis (coordinate flip - so technically y)
                                      df_plot_index, # output from index_biomarker_forest function, used to order the forest plot
                                      title_plot="", # Put the title of the forest plot here
                                      rev_axis=FALSE, #if you'd like to reverse the horizontal Axis
                                      labels_all=FALSE, # whether or not you want the axis ticks for analytes and the company labels
                                      amytag=FALSE){ #only flag this if wanting to hard set max of the range
  
  require(dplyr)
  df <- data.frame(DF1[,c("company","Rho","Rho_L","Rho_H","Biomarker")])
  df <- merge(df,df_plot_index,by=c("Biomarker","company"))
  rownames(df) <- NULL
  colnames(df) <- c("analyte","company","AUC","AUC_Lower","AUC_Upper","plot_index")
  
  
  company_count <- table(df_plot_index$company)
  Index1 <- -1
  Company1 <- names(company_count)[1]
  Index2 <- Index1 - company_count[1]
  Company2 <- names(company_count)[2]
  Index3 <- Index2 - company_count[2]
  Company3 <- names(company_count)[3]
  Index4 <- Index3 - company_count[3]
  Company4<- names(company_count)[4]
  Index5 <- Index4 - company_count[4]
  Company5 <- names(company_count)[5]
  Index6 <- Index5 - company_count[5]
  Company6 <- names(company_count)[6]
  Index7 <- Index6 - company_count[6]
  

  df$AUC <- as.numeric(df$AUC)
  df$AUC_Lower <- as.numeric(df$AUC_Lower)
  df$AUC_Upper <- as.numeric(df$AUC_Upper)
  df$company <- factor(df$company)
  
  df$company <-factor(df$company, levels = c(unique(df[order(-abs(df$AUC)),"company"])))
  DF <- df
  
  DF$analyte_c <- NA
  DF[which(DF$analyte=="p-tau217 (pg/ml)"),"analyte_c"] <- "p-Tau 217"
  DF[grep(paste0("A","\U03B2","42/A","\U03B2","40"),DF$analyte),"analyte_c"] <- "Ab42/Ab40"
  DF[which(DF$analyte=="p-tau217 ratio (%)"),"analyte_c"] <- "p-Tau 217 ratio"
  DF[grep("GFAP",DF$analyte),"analyte_c"] <- "GFAP"
  DF[grep("NfL",DF$analyte),"analyte_c"] <- "NfL"
  DF[grep("p-tau181",DF$analyte),"analyte_c"] <- "p-Tau 181"
  
  DF$analyte_c <- factor(DF$analyte_c,levels=c("p-Tau 217 ratio",
                                               "p-Tau 217",
                                               "Ab42/Ab40",
                                               "p-Tau 181",
                                               "GFAP","NfL"))
  
  
  analyte_colors <- c("Ab42/Ab40" = "blue",
                      "p-Tau 217" = "#008000", 
                      "p-Tau 181" = "#00FF00", 
                      "NfL" = "maroon", 
                      "GFAP" = "red",  
                      "p-Tau 217 ratio" = "#145A32")
  
  
  
DF[which(DF$analyte_c=="Ab42/Ab40"),"AUC"] <- DF[which(DF$analyte_c=="Ab42/Ab40"),"AUC"]*-1
DF[which(DF$analyte_c=="Ab42/Ab40"),"AUC_Lower"] <- DF[which(DF$analyte_c=="Ab42/Ab40"),"AUC_Lower"]*-1
DF[which(DF$analyte_c=="Ab42/Ab40"),"AUC_Upper"] <- DF[which(DF$analyte_c=="Ab42/Ab40"),"AUC_Upper"]*-1


round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

min_break_val <- min(breaks1)
min_y <-  round_any(min(c(DF$AUC_Lower,DF$AUC_Upper,min_break_val)) ,accuracy = 0.01 , f = floor) 
min_y_text <- min_y-0.1
text_start <- min_y_text



max_y <-  round_any(max(DF$AUC_Upper) ,accuracy = 0.2 , f = ceiling) 
max_y <- max(c(max(breaks1),max_y))
if(amytag==TRUE){
  max_y <- 0.9 
}

if(rev_axis==TRUE){
  max_y <- max_y+0.1
  text_start <- max_y
  
}




  Spearman_test_plot <-  ggplot(data=DF,
                                aes(x = plot_index, #uses index as the x-axis (flipped to y later)
                                    y = AUC, # y-axis is spearman rho, flipped to x-axis later
                                    ymin = min_y, ymax = max_y ))+ #limited of y-axis (flipped to x later)
    geom_rect(aes(xmin = Index1+0.5, xmax = Index2+0.5, ymin = min_y_text, ymax = max_y),
              fill = "#E5E4E2", alpha = 0.04) +
    geom_rect(aes(xmin = Index3+0.5, xmax = Index4+0.5, ymin = min_y_text, ymax = max_y),
              fill = "#E5E4E2", alpha = 0.04) +
    geom_rect(aes(xmin = Index5+0.5, xmax = Index6+0.5, ymin = min_y_text, ymax = max_y),
              fill = "#E5E4E2", alpha = 0.04) +
    geom_hline(aes(fill="black"),yintercept =indicatorline, linetype=2)+ #creates horizontal (flipped to vertical) line at 0
    geom_point(aes(col=analyte_c))+ #sets the coloring based on whether its the summary measure or not
    geom_errorbar(aes(ymin=(AUC_Lower), ymax=(AUC_Upper),col=analyte_c),width = 0, cex = 1,size=4)+ #formatting of forest lines
      labs(x = "", y = expression(paste("Spearman "   ,rho)))+theme_classic()+ #gives the title based on the function call variable
    theme(plot.title = element_text(hjust = 0.5,size = 10), #these are all text formatting
          axis.ticks.y=element_line(size = 1, color = "black"),
          axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6),
          axis.title=element_text(size=7,face="bold"),
          strip.text.y = element_blank(),
          legend.position =  "none")+ #no legend
    scale_color_manual(values=analyte_colors)+
    scale_fill_manual(values = analyte_colors)+ #makes sure dots and lines are colored based on color code
    coord_flip()+ggtitle(title_plot)
  
  if(rev_axis==TRUE){
    Spearman_test_plot <- Spearman_test_plot+scale_y_reverse(breaks=breaks1)
    
  }else {
    Spearman_test_plot <- Spearman_test_plot+scale_y_continuous(breaks=breaks1)
    
  }
  
  if(labels_all==TRUE){
Spearman_test_plot <- Spearman_test_plot+
      geom_text(aes(x = Index1, y = text_start, label = Company1, hjust = 0),size=3,family="Calibri") +
      geom_text(aes(x = Index2, y = text_start, label = Company2, hjust = 0),size=3,family="Calibri") +
      geom_text(aes(x = Index3, y = text_start, label = Company3, hjust = 0),size=3,family="Calibri") +
      geom_text(aes(x = Index4, y = text_start, label = Company4, hjust = 0),size=3,family="Calibri") +
      geom_text(aes(x = Index5, y = text_start, label = Company5, hjust = 0),size=3,family="Calibri") +
      geom_text(aes(x = Index6, y = text_start, label = Company6, hjust = 0),size=3,family="Calibri")+
     scale_x_continuous(breaks=DF$plot_index,labels=DF$analyte)
  }
  if(labels_all==FALSE){
    Spearman_test_plot <- Spearman_test_plot + theme(axis.text.y=element_blank(), 
                                                     axis.ticks.y=element_blank())
    
  }
  
  return(Spearman_test_plot)
  
}


####This function's sole purpose is to clean up output from the spearman_plasma_runs functions specifically tuned to our dataset and variables
## Will output a dataset with the following columnts and structure:
##company - self explanatory - company who made the assay	
##measure	- which analyte is being represented
##Rho	- Spearman's Rho with the particular analyte and the outcome, will also include 95% confidence interval - Rho (Rho Lower Bound - Rho Upper Bound)
##pcomp	- comparison between the top correlator Rho and the analyte - will read REFERENCE if analyte is the reference
##pRho	- Partial Spearman's Rho with the particular analyte and the outcome, will also include 95% confidence interval - Rho (Rho Lower Bound - Rho Upper Bound)
##parcomp 	- comparison between the top correlator partial Rho and the analyte - will read REFERENCE if analyte is the reference



clean_spearman_table <- function(DF1){ # output from dataset
  

  pparnormal_comp <- grep("Par_pComp_",colnames(DF1))
  pnormal_comp  <- grep("pComp_",colnames(DF1))[1]
  Reference_bio  <- DF1[which(abs(DF1$Rho)==max(abs(DF1$Rho))),"Biomarker"]
  Reference_comp  <- as.character(DF1[which(abs(DF1$Rho)==max(abs(DF1$Rho))),"company"])
  Pref_bio   <- DF1[which(abs(DF1$Partial_Rho)==max(abs(DF1$Partial_Rho))),"Biomarker"]
  Pref_com   <- as.character(DF1[which(abs(DF1$Partial_Rho)==max(abs(DF1$Partial_Rho))),"company"])
  
 spearman_out  <- data.frame(matrix(nrow=nrow(DF1),ncol=6))

 colnames(spearman_out) <- c("company","measure","Rho","pcomp","pRho","parcomp")
 
 spearman_out$company <- as.factor(DF1$company)
 spearman_out$company <- factor(spearman_out$company,levels = c("C2N", "Fujirebio", "ALZpath", "Janssen", "Roche", "Quanterix"))
 
 
 spearman_out$measure <- DF1$Biomarker
 spearman_out$Rho <- paste0(format(round(DF1$Rho,3),nsmall=3)," (",
                            format(round(DF1$Rho_L,3),nsmall=3)," to ",
                            format(round(DF1$Rho_H,3), nsmall=3),")")
 
 
 
 spearman_out[which(spearman_out$company==Reference_comp & spearman_out$measure==Reference_bio),"pcomp"] <- "REFERENCE"

DF_comp <-  DF1[which(!(spearman_out$company==Reference_comp & spearman_out$measure==Reference_bio)),]
DF_comp <- as.data.frame(DF_comp)
DF_comp$P_prep <- DF_comp[,pnormal_comp]
DF_comp[which(DF_comp[,pnormal_comp]=="<0.0001"),"P_prep"]  <- 0.0001
DF_comp$pcomp  <- format(round(p.adjust(as.numeric(DF_comp[,"P_prep"]),method = "BH"),4),nsmall=4)
DF_comp[which(DF_comp[,pnormal_comp]=="<0.0001"),"pcomp"] <- "<0.0001"
 spearman_out[which(spearman_out$company==Pref_com & spearman_out$measure==Pref_bio),"parcomp"] <- "REFERENCE"
for(i_comparisons in 1:nrow(DF_comp)){
  
topind <-  which(spearman_out$company==DF_comp$company[i_comparisons] & spearman_out$measure==DF_comp$Biomarker[i_comparisons])
spearman_out[topind,"pcomp"] <-  DF_comp$pcomp[i_comparisons]
  
}


spearman_out$pRho <- paste0(format(round(DF1$Partial_Rho,3),nsmall=3)," (",
                           format(round(DF1$PRho_L,3),nsmall=3)," to ",
                           format(round(DF1$PRho_H,3), nsmall=3),")")



DF_comp <-  DF1[which(!(spearman_out$company==Pref_com & spearman_out$measure==Pref_bio)),]
DF_comp <- as.data.frame(DF_comp)
DF_comp$Partial_prep <- DF_comp[,ppartial_comp]
DF_comp[which(DF_comp[,ppartial_comp]=="<0.0001"),"Partial_prep"]  <- 0.0001
DF_comp$ppart  <- format(round(p.adjust(as.numeric(DF_comp[,"Partial_prep"]),method = "BH"),4),nsmall=4)
DF_comp[which(DF_comp[,ppartial_comp]=="<0.0001"),"ppart"] <- "<0.0001"

for(i_comparisons in 1:nrow(DF_comp)){
  
  topind <-  which(spearman_out$company==DF_comp$company[i_comparisons] & spearman_out$measure==DF_comp$Biomarker[i_comparisons])
  spearman_out[topind,"parcomp"] <-  DF_comp$ppart[i_comparisons]
  
}
spearman_out[which(spearman_out$company==Pref_com & spearman_out$measure==Pref_bio),"parcomp"] <- "REFERENCE"
spearman_out <- spearman_out[order(spearman_out$company),]

return(spearman_out)
}


##creates heatmap from data and a "ref_data" table
##the ref data must include the following:
# analyte code - code that represents which analyte the measure represents
# var_id - the variable name for the analyte
# fullname - the label for the variable that you want to show up on the heatmap
## data.in is just the dataset
## order_cor_analyte is what analyte's correlations do you want to determine the order of the heatmap


create_cor_heatmap_for_project <- function(data.in,ref_data,order_cor_analyte="C2N_plasma_ptau217_ratio"){
  require(tidyverse)
  require(Hmisc)
  require(reshape2)
  
  mydata.cor = cor(data.in[,ref_data$var_id], method = c("spearman"),use = "pairwise.complete.obs")
  
  
  
  ref_data$ptau_correlation <- abs(data.frame(mydata.cor)[,order_cor_analyte])
  ref_reorderd <-  ref_data %>% group_by(analyte_code) %>% mutate(mx = max(ptau_correlation)) %>% 
    arrange(mx,ptau_correlation)
  ref_reorderd <- data.frame(ref_reorderd)
  
  
  mydata.cor <- mydata.cor[ref_reorderd$var_id,ref_reorderd$var_id]
  melted_cormat <- melt(mydata.cor)
  
  
  
  melted_cormat$Var1_Lab <- NA
  melted_cormat$Var2_Lab <- NA
  for(i_cycle in 1:nrow(ref_reorderd)){
    melted_cormat[which(melted_cormat$Var1==ref_reorderd$var_id[i_cycle]),"Var1_Lab"] <- ref_reorderd$fullname[i_cycle]
    melted_cormat[which(melted_cormat$Var2==ref_reorderd$var_id[i_cycle]),"Var2_Lab"] <- ref_reorderd$fullname[i_cycle]
  }
  
  melted_cormat$Var1_Lab <- factor(melted_cormat$Var1_Lab,levels=rev(ref_reorderd$fullname))
  melted_cormat$Var2_Lab <- factor(melted_cormat$Var2_Lab,levels=ref_reorderd$fullname)
  
  heatcolor_pallette <- c("#0087B3","white","red","black")
  color_scale <- c(0,0.5,0.999, 1)
  
  heatmap_of_cor <-  ggplot(melted_cormat, aes(Var1_Lab, Var2_Lab)) +
    geom_tile(aes(fill = abs(value))) +
    geom_text(aes(label = format(round(value, 2), nsmall = 2)),size=2) +
    scale_fill_gradientn(colours=heatcolor_pallette,values=color_scale)+
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 6, hjust = 1),axis.text.y=element_text(size=6),
          axis.title = element_blank())
  return(heatmap_of_cor)
}


##creates heatmap from data and a "ref_data" table
##the ref data must include the following:
# analyte code - code that represents which analyte the measure represents
# var_id - the variable name for the analyte
# fullname - the label for the variable that you want to show up on the heatmap
## data.in is just the dataset
## ref order is how you want to order the analyte groups
## order_cor_analyte is what analyte's correlations do you want to determine the order of the heatmap within the analyte groups

create_cor_heatmap_for_project_custom <- function(data.in,ref_data,reforder,order_cor_analyte="C2N_plasma_ptau217_ratio"){
  require(tidyverse)
  mydata.cor = cor(data.in[,ref_data$var_id], method = c("spearman"),use = "pairwise.complete.obs")
  
  
  
  ref_data$ptau_correlation <- abs(data.frame(mydata.cor)[,order_cor_analyte])
  ref_data$analyte_code <- factor(ref_data$analyte_code,levels=reforder)
  ref_reorderd <-  ref_data[order(ref_data$analyte_code,abs(ref_data$ptau_correlation)),]
  ref_reorderd <- data.frame(ref_reorderd)
  
  
  mydata.cor <- mydata.cor[ref_reorderd$var_id,ref_reorderd$var_id]
  melted_cormat <- melt(mydata.cor)
  
  
  
  melted_cormat$Var1_Lab <- NA
  melted_cormat$Var2_Lab <- NA
  for(i_cycle in 1:nrow(ref_reorderd)){
    melted_cormat[which(melted_cormat$Var1==ref_reorderd$var_id[i_cycle]),"Var1_Lab"] <- ref_reorderd$fullname[i_cycle]
    melted_cormat[which(melted_cormat$Var2==ref_reorderd$var_id[i_cycle]),"Var2_Lab"] <- ref_reorderd$fullname[i_cycle]
  }
  
  melted_cormat$Var1_Lab <- factor(melted_cormat$Var1_Lab,levels=rev(ref_reorderd$fullname))
  melted_cormat$Var2_Lab <- factor(melted_cormat$Var2_Lab,levels=ref_reorderd$fullname)
  
  heatcolor_pallette <- c("#0087B3","white","red","black")
  color_scale <- c(0,0.5,0.999, 1)
  
  heatmap_of_cor <-  ggplot(melted_cormat, aes(Var1_Lab, Var2_Lab)) +
    geom_tile(aes(fill = abs(value))) +
    geom_text(aes(label = format(round(value, 2), nsmall = 2)),size=2) +
    scale_fill_gradientn(colours=heatcolor_pallette,values=color_scale)+
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 6, hjust = 1),axis.text.y=element_text(size=6),
          axis.title = element_blank())
  return(heatmap_of_cor)
}

###

make_scatter_figure <- function(data.in, #input data
                                ref_data, #similar to ref_data for heatmaps, but also need group variable with the codes for the company
                                metric,  # outcome variable ID
                                results_dataset, #results dataset from spearman_plasma_runs
                                phenotype_name, #outcome name
                                rev_axis=F){ #do you want to reverse the x-axis
  
  require(tidyr)
  require(cowplot)
  
  
  graph_data <- data.in[which(!is.na(data.in[,metric])),]
  pheno_name <- phenotype_name
  top_performers <- results_dataset %>% group_by(group) %>% filter(abs(Rho)==max(abs(Rho)))
  top_performers <- as.data.frame(top_performers)
  
  
  graph_data$participant_condition <- NA
  graph_data[which(graph_data$CENTILOIDS_10==0 & graph_data$CDR_10==0), "participant_condition"] <- "A" #Blue: amyloid PET negative, CDR 0
  graph_data[which(graph_data$CENTILOIDS_10==0 & graph_data$CDR_10==1), "participant_condition"] <- "B" #Yellowish green: amyloid PET negative, CDR >0
  graph_data[which(graph_data$CENTILOIDS_10==1 & graph_data$CDR_10==0), "participant_condition"] <- "C" #Orange: amyloid PET positive, CDR 0
  graph_data[which(graph_data$CENTILOIDS_10==1 & graph_data$CDR_10==1), "participant_condition"] <- "D" #Red: amyloid PET positive, CDR >0
  
  
  color_palette_scatter <- c("A" = "blue",
                             "B" = "#9acd32",
                             "C" = "orange",
                             "D" = "red")
  
  
  
  i_group <- "c2n"
  i_lab <- top_performers[which(top_performers$group==i_group),"Biomarker"]
  best_company <- ref_data[which(ref_data$group=="c2n" & ref_data$label==i_lab),"var_id"]
  
  scatter_c2n <- ggplot(data=graph_data %>%
                          arrange(participant_condition),
                        aes_string(x=metric,y=best_company))+
    geom_point(aes(colour=factor(participant_condition)),size=2,alpha=0.5)+
    scale_color_manual(values = color_palette_scatter)+
    theme_classic()+xlab(pheno_name)+ylab(paste("C2N \n",i_lab)) + 
    theme(legend.position="none",axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
          axis.title=element_text(size=12))
  
  i_group <- "fuji"
  i_lab <- top_performers[which(top_performers$group==i_group),"Biomarker"]
  best_company <- ref_data[which(ref_data$group=="fuji" & ref_data$label==i_lab),"var_id"]
  
  scatter_fuji <- ggplot(data=graph_data %>%
                           arrange(participant_condition),aes_string(x=metric,y=best_company))+
    geom_point(aes(colour=factor(participant_condition)),size=2,alpha=0.5)+
    scale_color_manual(values = color_palette_scatter)+
    theme_classic()+xlab(pheno_name)+ylab(paste("Fujirebio \n",i_lab)) + 
    theme(legend.position="none",axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
          axis.title=element_text(size=12))
  
  
  i_group <- "jan"
  i_lab <- top_performers[which(top_performers$group==i_group),"Biomarker"]
  best_company <- ref_data[which(ref_data$group=="jan" & ref_data$label==i_lab),"var_id"]
  
  scatter_janssen <- ggplot(data=graph_data %>%
                              arrange(participant_condition),aes_string(x=metric,y=best_company))+
    geom_point(aes(colour=factor(participant_condition)),size=2,alpha=0.5)+
    scale_color_manual(values = color_palette_scatter)+
    theme_classic()+xlab(pheno_name)+ylab(paste("Janssen \n",i_lab)) + 
    theme(legend.position="none",axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
          axis.title=element_text(size=12))
  
  i_group <- "roche"
  i_lab <- top_performers[which(top_performers$group==i_group),"Biomarker"]
  best_company <- ref_data[which(ref_data$group=="roche" & ref_data$label==i_lab),"var_id"]
  
  scatter_roche <- ggplot(data=graph_data %>%
                            arrange(participant_condition),aes_string(x=metric,y=best_company))+
    geom_point(aes(colour=factor(participant_condition)),size=2,alpha=0.5)+
    scale_color_manual(values = color_palette_scatter)+
    theme_classic()+xlab(pheno_name)+ylab(paste("Roche \n",i_lab)) + 
    theme(legend.position="none",axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
          axis.title=element_text(size=12))
  
  i_group <- "alz"
  i_lab <- top_performers[which(top_performers$group==i_group),"Biomarker"]
  best_company <- ref_data[which(ref_data$group=="alz" & ref_data$label==i_lab),"var_id"]
  
  scatter_alzpath  <- ggplot(data=graph_data %>%
                               arrange(participant_condition),aes_string(x=metric,y=best_company))+
    geom_point(aes(colour=factor(participant_condition)),size=2,alpha=0.5)+
    scale_color_manual(values = color_palette_scatter)+
    theme_classic()+xlab(pheno_name)+ylab(paste("ALZpath \n",i_lab)) + 
    theme(legend.position="none",axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
          axis.title=element_text(size=12))
  
  i_group <- "quan"
  i_lab <- top_performers[which(top_performers$group==i_group),"Biomarker"]
  best_company <- ref_data[which(ref_data$group=="quan" & ref_data$label==i_lab),"var_id"]
  
  scatter_quanterix  <- ggplot(data=graph_data %>%
                                 arrange(participant_condition),
                               aes_string(x=metric,y=best_company))+
    geom_point(aes(colour=factor(participant_condition)),size=2,alpha=0.5)+
    scale_color_manual(values = color_palette_scatter)+
    theme_classic()+xlab(pheno_name)+ylab(paste("Quanterix \n",i_lab)) + 
    theme(legend.position="none",axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
          axis.title=element_text(size=12))
  
  
  if(rev_axis==T){
    scatter_c2n <- scatter_c2n + scale_x_reverse()
    scatter_fuji <- scatter_fuji + scale_x_reverse()
    scatter_alzpath <- scatter_alzpath + scale_x_reverse()
    scatter_janssen <- scatter_janssen + scale_x_reverse()
    scatter_roche <- scatter_roche + scale_x_reverse()
    scatter_quanterix <- scatter_quanterix + scale_x_reverse()
    
  }
  
  order_of_plots <- top_performers[order(abs(top_performers$Rho),decreasing = T),"group"]
  plot_set_up <- list("alz"=scatter_alzpath,
                      "c2n"=scatter_c2n,
                      "fuji"=scatter_fuji,
                      "jan"=scatter_janssen,
                      "quan"=scatter_quanterix,
                      "roche"=scatter_roche)
  
  
  plot_ordered <- plot_set_up[order_of_plots]
  
  plot_out <- plot_grid(plot_ordered[[1]],plot_ordered[[2]],
                        plot_ordered[[3]],plot_ordered[[4]],
                        plot_ordered[[5]],plot_ordered[[6]],nrow = 3,ncol = 2,labels = "AUTO")
  
  return(plot_out)
}


