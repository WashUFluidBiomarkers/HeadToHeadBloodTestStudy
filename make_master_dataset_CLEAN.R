# FNIH Project 2
# make longitudinal dataset
#### Coded: Kellen Petersen
#### Date as of commenting and documentation - 05/28/2024

# load libraries
library(tidyverse)
library(here)

# here::i_am("make_master_dataset_CLEAN.R")

# Parameters
year <- 365.25
cut_years <- 1
cut_centiloids <- 20
cut_centiloids_sens <- 37
cut_tau_mesialtemporal <- 1.32826246868779 # mixEM
cut_tau_temporoparietal <- 1.26937256044181 # mixEM
cut_atrophy <- 2.57231252744212 # mixEM


################################################################################
# DEMOGRAPHICS
demo0 <- read.csv(here("adni_data", "PTDEMOG_25Mar2024.csv"))
demo <- demo0 %>% 
  mutate(AGE = difftime(as.Date(VISDATE,"%m/%d/%Y"),as.Date(PTDOB,"%m/%d/%Y"), units = "days")) %>% 
  mutate(VISDATE_DEMO = VISDATE) %>% 
  select(RID,VISDATE_DEMO,PTDOB,AGE,PTGENDER,PTEDUCAT,PTRACCAT,PTETHCAT)
demo$AGE <- round(as.numeric(demo$AGE)/year,1)

demo <- demo %>% 
  mutate(RACE = ifelse(PTRACCAT == 4, 0, 
                       ifelse(PTRACCAT == 5, 1,
                              ifelse(PTRACCAT == -4, -4, 2 ))))
demo$VISDATE_DEMO <- format(as.Date(demo$VISDATE_DEMO, "%m/%d/%Y"), "%Y-%m-%d")
demo$PTDOB <- format(as.Date(demo$PTDOB, "%m/%d/%Y"), "%Y-%m-%d")


################################################################################
# APOE4: 3 different variables
apo0 <- read.csv(here("adni_data", "APOERES_02Apr2024.csv"))
apo0$APOE4_count <- rowSums(apo0[, c("APGEN1","APGEN2")] == 4)
apo0 <- apo0 %>% mutate(APOE4_positive = ifelse(APOE4_count > 0, 1, 0))
apo0$APOE_genotype <- paste(apo0$APGEN1, apo0$APGEN2, sep = "")
#apoe <- subset(apo, VISCODE == 'sc', c("RID", "APOE4_count", "APOE4_positive", "APOE_genotype"))
#apoe2 <- apoe[!duplicated(apoe$RID, fromLast = TRUE),]
apo <- apo0 %>% 
  select(RID,APOE4_count, APOE4_positive, APOE_genotype)


################################################################################
# DIAGNOSIS
dx0 <- read.csv(here("adni_data", "DXSUM_PDXCONV_02Apr2024.csv"))
dx <- dx0 %>% 
  mutate(EXAMDATE_DIAGNOSIS = EXAMDATE) %>% 
  mutate(DX_123 = DIAGNOSIS) %>% 
  select(RID, EXAMDATE_DIAGNOSIS, DX_123)
dx$EXAMDATE_DIAGNOSIS <- as.Date(dx$EXAMDATE_DIAGNOSIS, format = "%Y-%m-%d")


################################################################################
# CDR
cdr0 <- read.csv(here("adni_data", "CDR_25Mar2024.csv"))
cdr <- cdr0 %>% 
  mutate(VISDATE_CDR = VISDATE) %>%
  mutate(CDR = CDGLOBAL) %>%
  mutate(CDR_10 = ifelse(CDR == -1, -1, 
                         ifelse(CDR == 0, 0, 1))) %>% 
  mutate(CDR_SOB = CDMEMORY + CDORIENT + CDJUDGE + CDCOMMUN + CDHOME + CDCARE) %>% 
  select(RID, VISDATE_CDR, CDR, CDR_SOB, CDR_10)
cdr$VISDATE_CDR <- as.Date(cdr$VISDATE_CDR, format = "%Y-%m-%d")


################################################################################
# MMSE
mmse0 <- read.csv(here("adni_data", "MMSE_02Apr2024.csv"))

mmse0$VISCODE2 <- gsub("f", "bl", mmse0$VISCODE2)
mmse0$VISCODE2 <- gsub("sc", "bl", mmse0$VISCODE2)
mmse <- mmse0 %>% 
  mutate(VISDATE_MMSE = VISDATE) %>%
  select(RID, VISCODE2, VISDATE_MMSE, MMSCORE) %>% 
  arrange(RID, VISDATE_MMSE)
mmse$VISDATE_MMSE <- as.Date(mmse$VISDATE_MMSE, format = "%Y-%m-%d")


################################################################################
# CSF
# No PTID in the CSF file
csf0 <- read.csv(here("adni_data", "UPENNBIOMK_ROCHE_ELECSYS_02Apr2024.csv"))
csf <- csf0 %>% 
  mutate(EXAMDATE_CSF = EXAMDATE) %>%
  mutate(RUNDATE_CSF = RUNDATE) %>% 
  mutate(ABETA42_over_40  = ABETA42/ABETA40) %>%
  mutate(PTAU_over_ABETA42 = PTAU/ABETA42) %>% 
  select(RID, EXAMDATE_CSF, ABETA42, PTAU, TAU, ABETA42_over_40, PTAU_over_ABETA42)
csf$EXAMDATE_CSF <- as.Date(csf$EXAMDATE_CSF, format = "%Y-%m-%d")


################################################################################
# AMYLOID PET
amy0 <- read.csv(here("adni_data", "UCBERKELEY_AMY_6MM_25Mar2024.csv"))
amy <- amy0 %>% 
  mutate(SCANDATE_AMY = SCANDATE) %>%
  mutate(PROCESSDATE_AMY = PROCESSDATE) %>% 
  select(RID, SCANDATE_AMY,	CENTILOIDS) %>% 
  mutate(CENTILOIDS_10 = ifelse(CENTILOIDS > cut_centiloids, 1, 0)) %>% 
  mutate(CENTILOIDS_10_sens = ifelse(CENTILOIDS > cut_centiloids_sens, 1, 0))
amy$SCANDATE_AMY <- as.Date(amy$SCANDATE_AMY, format = "%Y-%m-%d")


################################################################################
# TAU PET
tau0 <- read.csv(here("adni_data", "ADNI_Tau_ROIs.csv"))
tau <- tau0 %>% 
  mutate(SCANDATE_TAU = SCANDATE) %>% 
  select(RID,SCANDATE_TAU,MesialTemporal,TemporoParietal) %>% 
  mutate(TAU_MesialTemporal_10 = ifelse(MesialTemporal > cut_tau_mesialtemporal, 1, 0)) %>% 
  mutate(TAU_TemporoParietal_10 = ifelse(TemporoParietal > cut_tau_temporoparietal, 1, 0))
tau$SCANDATE_TAU <- as.Date(tau$SCANDATE_TAU, format = "%Y-%m-%d")


################################################################################
# ATROPHY
atrophy0 <- read.csv(here("adni_data", "ADNI_metaROI_thickness_harmonized_.csv"))
atrophy <- atrophy0 %>% 
  mutate(SCANDATE_ATROPHY = EXAMDATE) %>%
  mutate(atrophy = metaROI.AgeAdj) %>%
  select(RID,SCANDATE_ATROPHY,atrophy) %>% 
  mutate(atrophy_10 = ifelse(atrophy < cut_atrophy, 1, 0))
atrophy$SCANDATE_ATROPHY <- as.Date(atrophy$SCANDATE_ATROPHY, format = "%Y-%m-%d")


################################################################################
# REMOVE 0 dataframes
rm(list = c("demo0", "apo0","dx0", "cdr0", "mmse0", "csf0", "amy0", "tau0", "atrophy0"))


################################################################################
# SORT
demo <- demo %>% arrange(RID,VISDATE_DEMO)
apo <- apo %>% arrange(RID)
cdr <- cdr %>% arrange(RID, VISDATE_CDR)
dx <- dx %>% arrange(RID, EXAMDATE_DIAGNOSIS)
mmse <- mmse %>% arrange(RID, VISDATE_MMSE)
csf <- csf %>% arrange(RID, EXAMDATE_CSF)
amy <- amy %>% arrange(RID, SCANDATE_AMY)
tau <- tau %>% arrange(RID, SCANDATE_TAU)
atrophy <- atrophy %>% arrange(RID, SCANDATE_ATROPHY)


################################################################################
# Load MASTER PLASMA file
plasma <- read.csv(here("adni_data", "Study_2_unblinded_datasets", "plasma_master_CLEAN.csv"))
plasma$EXAMDATE <- as.Date(plasma$EXAMDATE, format = "%Y-%m-%d")
IDs_plasma <- unique(plasma$RID)


################################################################################
# MERGE
m_1 <- merge(plasma, cdr, by = c("RID"), all = TRUE) %>%
  mutate(difftime_years = as.numeric(abs(difftime(EXAMDATE, VISDATE_CDR, units = "days")/year))) %>%
  group_by(RID_EXAMDATE) %>%
  slice(which.min(difftime_years)) %>%
  mutate(VISDATE_CDR = ifelse(all(difftime_years > cut_years), NA, VISDATE_CDR)) %>%
  mutate(CDR = ifelse(all(difftime_years > cut_years), NA, CDR)) %>%
  mutate(CDR_SOB = ifelse(all(difftime_years > cut_years), NA, CDR_SOB)) %>%
  mutate(CDR_10 = ifelse(all(difftime_years > cut_years), NA, CDR_10)) %>% 
  mutate(difftime_years = ifelse(all(difftime_years > cut_years), NA, difftime_years)) %>% 
  arrange(RID, EXAMDATE) %>% 
  select(-difftime_years)
m_1$VISDATE_CDR <- as.Date(m_1$VISDATE_CDR, format = "%Y-%m-%d")


m_2 <- merge(m_1, dx, by = c("RID"), all = TRUE) %>%
  mutate(difftime_years = as.numeric(abs(difftime(EXAMDATE, EXAMDATE_DIAGNOSIS, units = "days")/year))) %>%
  group_by(RID_EXAMDATE) %>%
  slice(which.min(difftime_years)) %>%
  mutate(EXAMDATE_DIAGNOSIS = ifelse(all(difftime_years > cut_years), NA, EXAMDATE_DIAGNOSIS)) %>%
  mutate(DX_123 = ifelse(all(difftime_years > cut_years), NA, DX_123))%>% 
  mutate(difftime_years = ifelse(all(difftime_years > cut_years), NA, difftime_years)) %>% 
  arrange(RID, EXAMDATE) %>% 
  select(-difftime_years)
m_2$EXAMDATE_DIAGNOSIS <- as.Date(m_2$EXAMDATE_DIAGNOSIS, format = "%Y-%m-%d")


m_4 <- merge(m_2, mmse, by = c("RID"), all = TRUE) %>%
  mutate(difftime_years = as.numeric(abs(difftime(EXAMDATE, VISDATE_MMSE, units = "days")/year))) %>%
  group_by(RID_EXAMDATE) %>%
  slice(which.min(difftime_years)) %>%
  mutate(VISDATE_MMSE = ifelse(all(difftime_years > cut_years), NA, VISDATE_MMSE)) %>%
  mutate(MMSCORE = ifelse(all(difftime_years > cut_years), NA, MMSCORE))%>% 
  mutate(difftime_years = ifelse(all(difftime_years > cut_years), NA, difftime_years)) %>% 
  arrange(RID, EXAMDATE) %>% 
  select(-difftime_years)
m_4$VISDATE_MMSE <- as.Date(m_4$VISDATE_MMSE, format = "%Y-%m-%d")


m_5 <- merge(m_4, csf, by = c("RID"), all = TRUE) %>%
  mutate(difftime_years = as.numeric(abs(difftime(EXAMDATE, EXAMDATE_CSF, units = "days")/year))) %>%
  mutate(difftime_years = ifelse(is.na(difftime_years),1000,difftime_years)) %>% 
  group_by(RID_EXAMDATE) %>%
  slice(which.min(difftime_years)) %>%
  mutate(EXAMDATE_CSF = ifelse(all(difftime_years > cut_years), NA, EXAMDATE_CSF)) %>%
  mutate(PTAU_over_ABETA42 = ifelse(all(difftime_years > cut_years), NA, PTAU_over_ABETA42))%>% 
  mutate(difftime_years = ifelse(all(difftime_years > cut_years), NA, difftime_years)) %>% 
  arrange(RID, EXAMDATE) %>% 
  select(-difftime_years) %>% 
  filter(RID %in% IDs_plasma)
m_5$EXAMDATE_CSF <- as.Date(m_5$EXAMDATE_CSF, format = "%Y-%m-%d")


m_6 <- merge(m_5, amy, by = c("RID"), all = TRUE) %>%
  mutate(difftime_years = as.numeric(abs(difftime(EXAMDATE, SCANDATE_AMY, units = "days")/year))) %>%
  mutate(difftime_years = ifelse(is.na(difftime_years),1000,difftime_years)) %>%
  group_by(RID_EXAMDATE) %>%
  slice(which.min(difftime_years)) %>%
  mutate(SCANDATE_AMY = ifelse(all(difftime_years > cut_years), NA, SCANDATE_AMY)) %>%
  mutate(CENTILOIDS = ifelse(all(difftime_years > cut_years), NA, CENTILOIDS))%>% 
  mutate(CENTILOIDS_10 = ifelse(all(difftime_years > cut_years), NA, CENTILOIDS_10))%>% 
  mutate(CENTILOIDS_10_sens = ifelse(all(difftime_years > cut_years), NA, CENTILOIDS_10_sens))%>% 
  mutate(difftime_years = ifelse(all(difftime_years > cut_years), NA, difftime_years)) %>% 
  arrange(RID, EXAMDATE) %>% 
  select(-difftime_years) %>% 
  filter(RID %in% IDs_plasma)
m_6$SCANDATE_AMY <- as.Date(m_6$SCANDATE_AMY, format = "%Y-%m-%d")


m_7 <- merge(m_6, tau, by = c("RID"), all = TRUE) %>%
  mutate(difftime_years = as.numeric(abs(difftime(EXAMDATE, SCANDATE_TAU, units = "days")/year))) %>%
  mutate(difftime_years = ifelse(is.na(difftime_years),1000,difftime_years)) %>%
  group_by(RID_EXAMDATE) %>%
  slice(which.min(difftime_years)) %>%
  mutate(SCANDATE_TAU = ifelse(all(difftime_years > cut_years), NA, SCANDATE_TAU)) %>%
  mutate(MesialTemporal = ifelse(all(difftime_years > cut_years), NA, MesialTemporal))%>% 
  mutate(TemporoParietal = ifelse(all(difftime_years > cut_years), NA, TemporoParietal))%>% 
  mutate(TAU_MesialTemporal_10 = ifelse(all(difftime_years > cut_years), NA, TAU_MesialTemporal_10))%>% 
  mutate(TAU_TemporoParietal_10 = ifelse(all(difftime_years > cut_years), NA, TAU_TemporoParietal_10))%>% 
  mutate(difftime_years = ifelse(all(difftime_years > cut_years), NA, difftime_years)) %>% 
  arrange(RID, EXAMDATE) %>% 
  select(-difftime_years)%>% 
  filter(RID %in% IDs_plasma)
m_7$SCANDATE_TAU <- as.Date(m_7$SCANDATE_TAU, format = "%Y-%m-%d")


m_8 <- merge(m_7, atrophy, by = c("RID"), all = TRUE) %>%
  mutate(difftime_years = as.numeric(abs(difftime(EXAMDATE, SCANDATE_ATROPHY, units = "days")/year))) %>%
  mutate(difftime_years = ifelse(is.na(difftime_years),1000,difftime_years)) %>%
  group_by(RID_EXAMDATE) %>%
  slice(which.min(difftime_years)) %>%
  mutate(SCANDATE_ATROPHY = ifelse(all(difftime_years > cut_years), NA, SCANDATE_ATROPHY)) %>%
  mutate(atrophy = ifelse(all(difftime_years > cut_years), NA, atrophy))%>% 
  mutate(atrophy_10 = ifelse(all(difftime_years > cut_years), NA, atrophy_10))%>% 
  mutate(difftime_years = ifelse(all(difftime_years > cut_years), NA, difftime_years)) %>% 
  arrange(RID, EXAMDATE) %>% 
  select(-difftime_years)%>% 
  filter(RID %in% IDs_plasma)
m_8$SCANDATE_ATROPHY <- as.Date(m_8$SCANDATE_ATROPHY, format = "%Y-%m-%d")


# ADD demo and apo at end
demo_noage <- demo %>% select(-AGE) %>% 
  group_by(RID) %>%
  slice_head(n = 1)
demo_noage_apo <- merge(demo_noage, apo, by = c("RID"), all = TRUE) %>% 
  arrange(RID, VISDATE_DEMO) %>% 
  select(-VISDATE_DEMO)

m_9 <- merge(m_8, demo_noage_apo, by = c("RID"), all=FALSE) %>% 
  arrange(RID, EXAMDATE)
m_9$PTDOB <- as.Date(m_9$PTDOB, format = "%Y-%m-%d")
m_9 <- m_9 %>% 
  mutate(AGE = as.numeric(abs(difftime(EXAMDATE, PTDOB, units = "days")/year))) %>%
  mutate(AGE = round(AGE,1)) %>% 
  arrange(RID, EXAMDATE)

# Write to csv
write.csv(m_9, here("adni_data","merged_data_CLEAN.csv"), row.names = FALSE)
