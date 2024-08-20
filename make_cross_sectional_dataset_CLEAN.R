# FNIH Project 2
# make cross-sectional dataset
#### Coded: Kellen Petersen
#### Date as of commenting and documentation - 05/28/2024


# Load libraries
library(tidyverse)
library(here)

#here::i_am("make_cross_sectional_dataset_CLEAN.R")

out_project <- here() #Folder where you want all output written to, a results and plot folder will be written here by default
out_results <- paste0(out_project,"/results/")
out_plots <- paste0(out_project,"/plots/")

data_file <- here("adni_data","merged_data_CLEAN.csv") # file name + path here
data.in00 <- read.csv(data_file)

# Change -1 in cdr file to NA
data.in000 <- data.in00 %>% 
  mutate(CDR = ifelse(CDR == -1, NA, CDR)) %>% 
  mutate(CDR_SOB = ifelse(CDR == 0, 0, CDR_SOB))

# GROUP
data.in0 <- data.in000

################################################################################
c2n <- read.csv(here("adni_data","Study_2_unblinded_datasets","C2N_precivityresults_formatted_unblind_simple.csv"))
roche<- read.csv(here("adni_data","Study_2_unblinded_datasets","U_Gothenburg_Elecsys results_formatted_unblind_simple.csv"))
fuji <- read.csv(here("adni_data","Study_2_unblinded_datasets","Indiana_U_Lumipulse Results_formatted_unblind_data_simple.csv"))
alzpath<- read.csv(here("adni_data","Study_2_unblinded_datasets","Quanterix_ALZpath_pTau217_results_formatted_unblind_simple.csv"))
janssen<- read.csv(here("adni_data","Study_2_unblinded_datasets","Quanterix_Janssen_pTau217_results_formatted_unblind_simple.csv"))
quanterix<- read.csv(here("adni_data","Study_2_unblinded_datasets","Quanterix_Simoa_pTau181_results__formatted_unblind_simple.csv"))
quanterix2<- read.csv(here("adni_data","Study_2_unblinded_datasets","Quanterix_Simoa_N4PE_results_formatted_unblind_simple.csv"))

IDs_other5 <- unique(c2n$RID)
IDs_quant <- unique(quanterix$RID)

ID_both <- intersect(IDs_other5, IDs_quant)
ID_only_other5 <- setdiff(IDs_other5, IDs_quant)
ID_only_quant <- setdiff(IDs_quant, IDs_other5)
ID_all <- union(ID_both, union(ID_only_other5, ID_only_quant))

ID_other5 <- union(ID_only_other5, ID_both)
ID_quant <- union(ID_only_quant, ID_both)

ID_2tp <- data.in0 %>% group_by(RID) %>% summarise(n_time = n()) %>% filter(n_time == 2) %>% pull(RID)
ID_3tp <- data.in0 %>% group_by(RID) %>% summarise(n_time = n()) %>% filter(n_time == 3) %>% pull(RID)
ID_4tp <- data.in0 %>% group_by(RID) %>% summarise(n_time = n()) %>% filter(n_time == 4) %>% pull(RID)
ID_5tp <- data.in0 %>% group_by(RID) %>% summarise(n_time = n()) %>% filter(n_time == 5) %>% pull(RID)

vars_ID <- c("RID","EXAMDATE","RID_EXAMDATE")
vars_c2n <- c("C2N_plasma_Abeta40", "C2N_plasma_Abeta42",
              "C2N_plasma_Abeta42_Abeta40", "C2N_plasma_ptau217",
              "C2N_plasma_nptau217", "C2N_plasma_ptau217_ratio")
vars_fuji <- c("Fuji_plasma_ptau217", "Fuji_plasma_Ab40", "Fuji_plasma_Ab42",
               "Fuji_plasma_Ab42_Ab40") 
vars_alzpath <- c("AlzPath_plasma_ptau217") 
vars_janssen <- c("Janssen_plasma_ptau217")
vars_roche <- c("Roche_plasma_Ab40", "Roche_plasma_Ab42", "Roche_plasma_GFAP",
                "Roche_plasma_NfL", "Roche_plasma_ptau181", "Roche_plasma_Ab42_Ab40")
vars_quanterix <- c("QX_plasma_ptau181", "QX_plasma_Ab40", "QX_plasma_Ab42",
                    "QX_plasma_GFAP", "QX_plasma_NfL", "QX_plasma_Ab42_Ab40")
vars_plasma5 <- c(vars_c2n, vars_fuji, vars_alzpath, vars_janssen, vars_roche)
vars_plasma <- c(vars_c2n, vars_fuji, vars_alzpath, vars_janssen, vars_roche, vars_quanterix)

vars_cdr <- c("VISDATE_CDR", "CDR", "CDR_SOB", "CDR_10")
vars_dx <- c("EXAMDATE_DIAGNOSIS", "DX_123")
vars_nptests <- c("VISDATE_MMSE", "MMSCORE")
vars_CSF <- c("EXAMDATE_CSF", "PTAU_over_ABETA42")
vars_AB <- c("SCANDATE_AMY", "CENTILOIDS", "CENTILOIDS_10", "CENTILOIDS_10_sens")
vars_TAU <- c("SCANDATE_TAU", "MesialTemporal", "TemporoParietal", 
              "TAU_MesialTemporal_10", "TAU_TemporoParietal_10")
vars_MRI <- c("SCANDATE_ATROPHY", "atrophy", "atrophy_10")
vars_demo <- c("PTDOB", "PTGENDER", "PTEDUCAT", "PTRACCAT", "PTETHCAT", "RACE", 
               "APOE4_count", "APOE4_positive", "APOE_genotype", "AGE")

# How many for each visit?
data.in0_ind <- data.in0 %>% 
  group_by(RID) %>% 
  mutate(index_RID_timepoint = row_number()) %>% 
  ungroup()

################################################################################
#IDs <- unique(data.in0_ind$RID)
IDs <- ID_only_other5
vars_complete <- c(vars_ID, vars_plasma5, vars_cdr, vars_AB,vars_demo)
# IDs <- unique(data.in0$RID)
tmp1 <- data.in0_ind %>% 
  filter(RID %in% IDs) %>% 
  select(vars_complete) %>% 
  filter(complete.cases(.))
cross1 <- data.in0_ind %>%
  filter(RID_EXAMDATE %in% tmp1$RID_EXAMDATE) %>%
  group_by(RID) %>%
  filter(EXAMDATE == max(EXAMDATE))

################################################################################
# ID_both <- intersect(IDs_other5, IDs_quant)
IDs <- ID_both
vars_complete <- c(vars_ID, vars_plasma5, vars_quanterix, vars_cdr, vars_AB,vars_demo)
#IDs <- unique(data.in0$RID)
tmp2 <- data.in0_ind %>% 
  filter(RID %in% IDs) %>% 
  select(vars_complete) %>% 
  filter(complete.cases(.))
cross2 <- data.in0_ind %>%
  filter(RID_EXAMDATE %in% tmp2$RID_EXAMDATE) %>%
  group_by(RID) %>%
  filter(EXAMDATE == max(EXAMDATE))

################################################################################
# # ID_only_quant <- setdiff(IDs_quant, IDs_other5)
# IDs <- ID_only_quant
# vars_complete <- c(vars_ID, vars_quanterix, vars_cdr, vars_AB,vars_demo)
# #IDs <- unique(data.in0$RID)
# tmp3 <- data.in0_ind %>% 
#   filter(RID %in% IDs) %>% 
#   select(vars_complete) %>% 
#   filter(complete.cases(.))
# cross3 <- data.in0_ind %>%
#   filter(RID_EXAMDATE %in% tmp3$RID_EXAMDATE) %>%
#   group_by(RID) %>%
#   filter(EXAMDATE == max(EXAMDATE))

cross <- rbind(cross1,cross2) %>% 
  arrange(RID, EXAMDATE)

II <- which(cross$RID %in% union(ID_both,ID_only_other5))
cross$other5YN <- 0
cross$other5YN[II] <- 1

II <- which(cross$RID %in% union(ID_both,ID_only_quant))
cross$quantYN <- 0
cross$quantYN[II] <- 1

cross2 <- cross %>% 
  mutate(RACE = ifelse(RACE==1,1,
                       ifelse(RACE==2,3,
                              ifelse(RACE==0,2,-99))))

# Write to file
write.csv(cross2, here("adni_data","cross_sectional_data_CLEAN.csv"), row.names = FALSE)
