# FNIH Project 2
# make AUC Forest plots
#### Coded: Kellen Petersen
#### Date as of commenting and documentation - 05/31/2024

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
library(scales)
library(patchwork)

# Load functions
source("functions_forest.R")

# Load data
here::i_am("make_forest_using_function.R")
here()
# Load data
data_file <- here("cross_sectional_data_CLEAN.csv")
data.in00 <- read.csv(data_file)

# Choose cohort and coavariates
data.in0 <- data.in00 # Entire cohort
covars_TF <- FALSE # without covariates

# Make models for outcome 1
vars_outcome <- c("CENTILOIDS_10")
T1 <- make_model(data.in0,covars_TF,vars_outcome)
breaks1 = c(0, 0.50, 0.733, 1)
AUC_indicatorline = 0.733
Y <- c(.4, 1)
p1 <- getAUCs_withLabels(T1,breaks1,AUC_indicatorline,rev_axis=FALSE,Y[1],Y[2])
p1
T1$key <- paste(T1$company, T1$analyte)

# Make models for outcome 2
vars_outcome <- c("TAU_MesialTemporal_10")
T2 <- make_model(data.in0,covars_TF,vars_outcome)
T2$key <- paste(T2$company, T2$analyte)
order_index2 <- match(T1$key, T2$key)
T2_ordered <- T2[order_index2, ]
T2_ordered$key <- NULL
breaks1 =c(0, 0.50, 0.767, 1)
AUC_indicatorline = 0.767
Y <- c(.35, 1)
p2 <- getAUCs_withoutLabels(T2_ordered,breaks1,AUC_indicatorline,rev_axis=FALSE,Y[1],Y[2])

# Make models for outcome 3
vars_outcome <- c("atrophy_10")
T3 <- make_model(data.in0,covars_TF,vars_outcome)
T3$key <- paste(T3$company, T3$analyte)
order_index3 <- match(T1$key, T3$key)
T3_ordered <- T3[order_index3, ]
T3_ordered$key <- NULL
breaks1 = c(0.50, 0.758, 1)
AUC_indicatorline = 0.758
Y <- c(0.35, 1)
p3 <- getAUCs_withoutLabels(T3_ordered,breaks1,AUC_indicatorline,rev_axis=FALSE,Y[1],Y[2])

# Make models for outcome 4
data.in0$impaired_10 <- ifelse(data.in0$CDR > 0, 1, 0)
data.in0$impaired_10 <- as.factor(data.in0$impaired_10)
vars_outcome <- c("impaired_10")
T4 <- make_model(data.in0,covars_TF,vars_outcome)
T4$key <- paste(T4$company, T4$analyte)
order_index4 <- match(T1$key, T4$key)
T4_ordered <- T4[order_index4, ]
T4_ordered$key <- NULL
breaks1 = c(0, 0.45, 0.559, .75)
AUC_indicatorline = 0.559
Y <- c(.35, .85)
p4 <- getAUCs_withoutLabels(T4_ordered,breaks1,AUC_indicatorline,rev_axis=FALSE,Y[1],Y[2])


# Add titles
p1b <- p1 +
  ggtitle("Amyloid PET") +
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        plot.subtitle = element_text(hjust = 0.5))
p2b <- p2 +
  ggtitle("Early tau PET") +
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        plot.subtitle = element_text(hjust = 0.5))
p3b <- p3 + 
  ggtitle("Cortical thickness") +
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        plot.subtitle = element_text(hjust = 0.5))
p4b <- p4 +
  ggtitle("Cognitive impairment") +
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        plot.subtitle = element_text(hjust = 0.5))


# Combine plots
pp <- p1b | p2b | p3b | p4b
pp

# Save plot
# ggsave("Forest_AUC_CLEAN.png", pp, device = "png", dpi = 500, width = 11, height = 6, units = "in")