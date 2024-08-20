# FNIH Project 2
# combine all plasma data
#### Coded: Kellen Petersen
#### Date as of commenting and documentation - 05/28/2024

# Merge unblinded data
library(tidyverse)
library(lubridate)
library(here)

#here::i_am("make_plasma_data_CLEAN.R")

# Load files
c2n <- read.csv(here("adni_data","Study_2_unblinded_datasets","C2N_precivityresults_formatted_unblind_simple.csv"))
roche<- read.csv(here("adni_data","Study_2_unblinded_datasets","U_Gothenburg_Elecsys results_formatted_unblind_simple.csv"))
fuji <- read.csv(here("adni_data","Study_2_unblinded_datasets","Indiana_U_Lumipulse Results_formatted_unblind_data_simple.csv"))
alzpath<- read.csv(here("adni_data","Study_2_unblinded_datasets","Quanterix_ALZpath_pTau217_results_formatted_unblind_simple.csv"))
janssen<- read.csv(here("adni_data","Study_2_unblinded_datasets","Quanterix_Janssen_pTau217_results_formatted_unblind_simple.csv"))
quanterix<- read.csv(here("adni_data","Study_2_unblinded_datasets","Quanterix_Simoa_pTau181_results__formatted_unblind_simple.csv"))
quanterix2<- read.csv(here("adni_data","Study_2_unblinded_datasets","Quanterix_Simoa_N4PE_results_formatted_unblind_simple.csv"))


#############################################################################
# Make datasets "wide"
# Calculate AB ratio variables as needed

c2n_wide <- c2n %>% 
  select(RID,EXAMDATE,VISCODE2,TESTNAME,TESTVALUE) %>% 
  pivot_wider(names_from = "TESTNAME", values_from = "TESTVALUE")
names(c2n_wide) <- c("RID","EXAMDATE","VISCODE2",
                     "C2N_plasma_Abeta40", "C2N_plasma_Abeta42",
                     "C2N_plasma_Abeta42_Abeta40",
                     "C2N_plasma_ptau217","C2N_plasma_nptau217",
                     "C2N_plasma_ptau217_ratio")


roche_wide <- roche %>% 
  select(RID,EXAMDATE,VISCODE2,TESTNAME,TESTVALUE) %>% 
  pivot_wider(names_from = "TESTNAME", values_from = "TESTVALUE")
names(roche_wide) <- c("RID","EXAMDATE","VISCODE2",
                       "Roche_plasma_Ab40","Roche_plasma_Ab42",
                       "Roche_plasma_GFAP","Roche_plasma_NfL",
                       "Roche_plasma_ptau181")
roche_wide$Roche_plasma_Ab40 <- as.numeric(roche_wide$Roche_plasma_Ab40)
roche_wide$Roche_plasma_Ab42 <- as.numeric(roche_wide$Roche_plasma_Ab42)
roche_wide$Roche_plasma_GFAP <- as.numeric(roche_wide$Roche_plasma_GFAP)
roche_wide$Roche_plasma_NfL <- as.numeric(roche_wide$Roche_plasma_NfL)
roche_wide$Roche_plasma_ptau181 <- as.numeric(roche_wide$Roche_plasma_ptau181)
roche_wide <- roche_wide %>% 
  mutate(Roche_plasma_Ab42_Ab40 = Roche_plasma_Ab42/Roche_plasma_Ab40 * 0.001)


fuji_wide <- fuji %>% 
  select(RID,EXAMDATE,VISCODE2,TESTNAME,TESTVALUE) %>% 
  pivot_wider(names_from = "TESTNAME", values_from = "TESTVALUE")
names(fuji_wide) <- c("RID","EXAMDATE","VISCODE2",
                      "Fuji_plasma_ptau217",
                      "Fuji_plasma_Ab40","Fuji_plasma_Ab42")
fuji_wide <- fuji_wide %>% 
  mutate(Fuji_plasma_Ab42_Ab40 = Fuji_plasma_Ab42/Fuji_plasma_Ab40)


alzpath_wide <- alzpath %>% 
  select(RID,EXAMDATE,VISCODE2,TESTNAME,TESTVALUE) %>% 
  pivot_wider(names_from = "TESTNAME", values_from = "TESTVALUE")
names(alzpath_wide) <- c("RID","EXAMDATE","VISCODE2",
                         "AlzPath_plasma_ptau217")
alzpath_wide$AlzPath_plasma_ptau217 <- as.numeric(alzpath_wide$AlzPath_plasma_ptau217)


janssen_wide <- janssen %>% 
  select(RID,EXAMDATE,VISCODE2,TESTNAME,TESTVALUE) %>% 
  pivot_wider(names_from = "TESTNAME", values_from = "TESTVALUE")
names(janssen_wide) <- c("RID","EXAMDATE","VISCODE2",
                         "Janssen_plasma_ptau217")


quanterix_wide <- quanterix %>% 
  select(RID,EXAMDATE,VISCODE2,TESTNAME,TESTVALUE) %>% 
  pivot_wider(names_from = "TESTNAME", values_from = "TESTVALUE")
names(quanterix_wide) <- c("RID","EXAMDATE","VISCODE2",
                           "QX_plasma_ptau181")
quanterix_wide$QX_plasma_ptau181 <- as.numeric(quanterix_wide$QX_plasma_ptau181)


quanterix2_wide <- quanterix2 %>% 
  select(RID,EXAMDATE,VISCODE2,TESTNAME,TESTVALUE) %>% 
  pivot_wider(names_from = "TESTNAME", values_from = "TESTVALUE")
names(quanterix2_wide) <- c("RID","EXAMDATE","VISCODE2",
                            "QX_plasma_Ab40","QX_plasma_Ab42",
                            "QX_plasma_GFAP","QX_plasma_NfL")
quanterix2_wide$QX_plasma_Ab40 <- as.numeric(quanterix2_wide$QX_plasma_Ab40)
quanterix2_wide$QX_plasma_Ab42 <- as.numeric(quanterix2_wide$QX_plasma_Ab42) # 1378 missing Ab42
quanterix2_wide$QX_plasma_GFAP <- as.numeric(quanterix2_wide$QX_plasma_GFAP)
quanterix2_wide$QX_plasma_NfL <- as.numeric(quanterix2_wide$QX_plasma_NfL)
quanterix2_wide$RID <- as.factor(quanterix2_wide$RID)
quanterix2_wide <- quanterix2_wide %>% 
  mutate(QX_plasma_Ab42_Ab40 = QX_plasma_Ab42 / QX_plasma_Ab40)


# Combine Quanterix datasets
QX <- merge(quanterix_wide, quanterix2_wide, by=c("RID","EXAMDATE"), all=TRUE)
QX <- QX %>% arrange(RID,EXAMDATE)


# Remove VISCODEs
c2n_wide_b <- c2n_wide %>% select(-VISCODE2)
fuji_wide_b <- fuji_wide %>% select(-VISCODE2)
alzpath_wide_b <- alzpath_wide %>% select(-VISCODE2)
janssen_wide_b <- janssen_wide %>% select(-VISCODE2)
roche_wide_b <- roche_wide %>% select(-VISCODE2)
QX_b <- QX


# Merge
df  <- merge(c2n_wide_b, fuji_wide_b, by=c("RID","EXAMDATE"), all = TRUE)
df2 <- merge(df, alzpath_wide_b, by=c("RID","EXAMDATE"), all = TRUE)
df3 <- merge(df2, janssen_wide_b, by=c("RID","EXAMDATE"), all = TRUE)
df4 <- merge(df3, roche_wide_b, by=c("RID","EXAMDATE"), all = TRUE)
df5 <- merge(df4, QX_b, by=c("RID","EXAMDATE"), all = TRUE)


# Add RID_EXAMDATE and arrange
df5$EXAMDATE <- as.Date(df5$EXAMDATE, format = "%m/%d/%y")
df6 <- df5 %>% arrange(RID,EXAMDATE) %>% 
  mutate(RID_EXAMDATE = paste(RID,EXAMDATE,sep = "_")) %>% 
  relocate(RID_EXAMDATE, .after = 2)


# Write to csv
write.csv(df6, here("adni_data","Study_2_unblinded_datasets","plasma_master_CLEAN.csv"), row.names=FALSE)