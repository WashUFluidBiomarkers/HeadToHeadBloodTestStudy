# FNIH Project 2
# make ROC plots
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
library(showtext)
library(ggtext)  # Load the ggtext package
library(ggpubr)  # Load the ggpubr package

#here::i_am("make_roc_plots_CLEAN.R")

# Load make_model file
source(here("make_model_roc_plot_CLEAN.R"))

# Load data
data_file <- here("adni_data","cross_sectional_data_2024_04_30.csv")
data.in00 <- read.csv(data_file)

# Choose group (here: entire cohort)
df <- data.in00

# Parameters: with covariates (True) or not (False)
covars_TF <- FALSE

# Outcome variables
vars_outcome <- c("CENTILOIDS_10")
#vars_outcome <- c("TAU_MesialTemporal_10")
#vars_outcome <- c("atrophy_10")
#vars_outcome <- c("impaired_10")

# Run Make Models
result1 <- make_model(df,covars_TF,vars_outcome)
T1 <- result1$T
M1 <- result1$MODELS_df
REF <- result1$REF

# Re-arrange data frame and choose best performing models
Tbest <- T1 %>%
  group_by(company) %>%
  slice(1) %>% 
  ungroup() %>% 
  arrange(desc(AUC))

# Make labels
Tbest$analyte2 <- Tbest$analyte
Tbest$analyte2 <- gsub("b", "β", Tbest$analyte2)
Tbest$company <- gsub("AlzPath", "ALZpath", Tbest$company)
Tbest$labels <- paste0(Tbest$company, ": ", Tbest$analyte2)
order <- Tbest$original_order
order <- c(order,26)
labels <- Tbest$labels
labels[which(labels == "covariates: covariates")] <- "Covariates (age, sex, APOE genotype)"

# Make ROC plot list
roc_list <- setNames(
  object = list(
    M1[[order[1]]][[2]],
    M1[[order[2]]][[2]],
    M1[[order[3]]][[2]],
    M1[[order[4]]][[2]],
    M1[[order[5]]][[2]],
    M1[[order[6]]][[2]],
    M1[[25]][[2]]
  ),
  nm = labels
)

# Colors for the legend
company_colors <- c("blue", "darkgreen", "red", "purple", "darkorange", "brown", "black")
company_names <- c("C2N", "Fujirebio", "ALZpath", "Janssen", "Roche", "Quanterix", "covariates")

# Create a new column 'colors' in your data frame
Tbest$colors <- company_colors[match(Tbest$company, company_names)]
COLS <- Tbest$colors

# Plotting with ggroc and adding a legend
title0 <- "Amyloid PET"
roc_plot <- ggroc(roc_list, linetype = 1, size = 1) +
  geom_abline(intercept = 1, slope = 1, linetype = "dotted", color = "black") +
  scale_color_manual(values = COLS, name = "Best performing models") +
  theme_cowplot(12) +
  coord_fixed() +
  labs(x = "NPA", y = "PPA", title = title0) + 
  theme(plot.title = element_text(hjust = 0.5))
roc_plot <- ggpar(roc_plot,
                  legend.text = ggtext::element_markdown(), 
                  legend.title = ggtext::element_markdown())

# Save plot
ggsave(here("results_clean","roc_Amyloid_CLEAN.jpg"), roc_plot, device = "jpg", dpi = 500, width = 11, height = 6, units = "in")