###Spearman Scatter Plots
#### Coded: Benjamin Saef
#### Date as of commenting and documentation - 05/28/2024


#Scatter plots for best performing analytes from each platform (Benjamin) versus amyloid PET (x-axis):
#  Yellowish green: amyloid PET negative, CDR >0
#Blue: amyloid PET negative, CDR 0
#Orange: amyloid PET positive, CDR 0
#Red: amyloid PET positive, CDR >0
#a.	Amyloid PET Centiloid
#b.	Early Tau PET
#c.	Brain atrophy
#d.	Cognitive impairment
library(tidyr)
library(data.table)
library(RVAideMemoire)
library(pROC)
library(RColorBrewer)
library(MuMIn)


####Set project directory
out_project <- "/FNIH_Project_2/Paper1/Spearman_Analysis" #Folder where you want all output written to, a results and plot folder will be written here by default
out_plots <- paste0(out_project,"/plots/scatter/")


dir.create(out_plots,showWarnings = F)

data_file <- "/FNIH_Project_2/DATA_STUDY_2/cross_sectional_data_2024_04_24b.csv"  # input cross-sectional data

data.in <- read.csv(data_file)

#analysis vars
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

 
l_var <- c( "ab4240", #Ab42/40
            "pt217", #ptau217
            "APS", #APS
            "ab4240", #Ab42/40
            "pt181", #ptau181
            "gfap", #GFAP
            "nfl", #NFL
            "ab4240",#Ab42/40
            "pt217",#ptau217
            "pt217",#ptau217
            "pt217",#ptau217
            "ab4240", #Ab42/40
            "pt181", #ptau181
            "gfap", #GFAP
            "nfl"#NFL
)

##list out the different analytes - Labels
analytes <- c(paste0("A","\U03B2","42/A","\U03B2","40"), #Ab42/40
              "p-tau217 ratio (%)", #ptau217
              "p-tau217 (pg/ml)", #ptau217
              paste0("A","\U03B2","42/A","\U03B2","40"), #Ab42/40
              "p-tau181 (pg/ml)", #ptau181
              "GFAP (ng/ml)", #GFAP
              "NfL (pg/mL)", #NFL
              paste0("A","\U03B2","42/A","\U03B2","40"), #Ab42/40
              "p-tau217 (pg/ml)", #ptau217
              "p-tau217 (pg/ml)", #ptau217
              "p-tau217 (pg/ml)", #ptau217
              paste0("A","\U03B2","42/A","\U03B2","40"), #Ab42/40
              "p-tau181 (pg/ml)", #ptau181
              "GFAP (pg/ml)", #GFAP
              "NfL (pg/mL)" #NFL
)
#list out the variables representing those analytes.
i_var <- c(abeta4240_c2n, #Ab42/40 - c2n
           pTau217_c2n, #ptau217 - c2n
           APS_c2n,#aps - c2n
           abeta4240_esys,#Ab42/40 - elecsys
           ptau_181_esys,#ptau181 - elecsys
           GFAP_esys,#GFAP - elecsys
           NfL_esys, #NFL - elecsys
           abeta4240_lumi,#Ab42/40 - lumipulse
           pTau217_lumi, #ptau217 - lumipulse
           alzpath_ptau217,#ptau217 - alzpath
           janssen_ptau217,#ptau217 - lumipulse
           abeta4240_quan,#Ab42/40 - quanterix
           ptau_181_quan,#ptau181 - quanterix
           GFAP_quan,#GFAP - quanterix
           NfL_quan #NFL - quanterix
)


ref_data <- data.frame(group=l_var,var_id=i_var,label=analytes)

###Create Scatter plot figure for Centiloid results
load("/FNIH_Project_2/Paper1/Spearman_Analysis/results/cent_correlations.rdata") #load in results data from Spearman_Analysis.R code
plot_out <- make_scatter_figure(data.in,
                    ref_data,
                    o_vars[1],
                    list_of_cent_correlation_dataset[[1]],
                    "Amyloid PET Centiloid")
ggsave(paste0(out_plots,"Centiloid_Scatter_ATN.png"),
       plot=plot_out,
       device=png,dpi = 1200,width = 8,height=12)

###Create Scatter plot figure for Early Tau results
load("/FNIH_Project_2/Paper1/Spearman_Analysis/results/tauearly_correlations.rdata") #load in results data from Spearman_Analysis.R code
plot_out <- make_scatter_figure(data.in,
                                ref_data,
                                o_vars[2],
                                list_of_tauearly_correlation_dataset[[1]],
                                "Early Tau PET")
ggsave(paste0(out_plots,"EarlyTau_Scatter_ATN.png"),
       plot=plot_out,
       device=png,dpi = 1200,width = 8,height=12)

load("/FNIH_Project_2/Paper1/Spearman_Analysis/results/taulate_correlations.rdata") #load in results data from Spearman_Analysis.R code
###Create Scatter plot figure for Late Tau results
plot_out <- make_scatter_figure(data.in,
                                ref_data,
                                o_vars[3],
                                list_of_TauLate_correlation_dataset[[1]],
                                "Late Tau PET")
ggsave(paste0(out_plots,"LateTau_Scatter_ATN.png"),
       plot=plot_out,
       device=png,dpi = 1200,width = 8,height=12)


load("/FNIH_Project_2/Paper1/Spearman_Analysis/results/atrophy_correlations.rdata") #load in results data from Spearman_Analysis.R code
###Create Scatter plot figure for Atrophy/Cortical Thickness
plot_out <- make_scatter_figure(data.in,
                                ref_data,
                                o_vars[4],
                                list_of_Atrophy_correlation_dataset[[1]],
                                "Cortical thickness",rev_axis = T)
ggsave(paste0(out_plots,"Atrophy_Scatter_ATN.png"),
       plot=plot_out,
       device=png,dpi = 1200,width = 8,height=12)

load("/FNIH_Project_2/Paper1/Spearman_Analysis/results/cdrsumbox_correlations.rdata") #load in results data from Spearman_Analysis.R code
###Create Scatter plot figure for CDR Sum of boxes
plot_out <- make_scatter_figure(data.in,
                                ref_data,
                                o_vars[5],
                                list_of_sumbox_correlation_dataset[[1]],
                                "Dementia severity [CDR Sum of Boxes]")
ggsave(paste0(out_plots,"CDRSumbox_Scatter_ATN.png"),
       plot=plot_out,
       device=png,dpi = 1200,width = 8,height=12)


