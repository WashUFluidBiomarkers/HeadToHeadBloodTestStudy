####Heatmap creation
#### Coded: Benjamin Saef
#### Date as of commenting and documentation - 05/28/2024


load("/FNIH_Project_2/Current_Analysis_data.rdata")
library(Hmisc)
library(reshape2)
out_project <- "/FNIH_Project_2/Paper1/Spearman_Analysis" #Folder where you want all output written to, a results and plot folder will be written here by default
out_plots <- paste0(out_project,"/plots/heatmaps/")
dir.create(out_plots,showWarnings = F)

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

ptau217_ab4240_lumi <- "" #Variable representing AB42/AB40 + p-tau217 ratio

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


a_var <- c( "ab4240", #Ab42/40
            "pt217", #ptau217
            "pt217", #pt217
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

l_var <- c( "c2n", #Ab42/40
            "c2n", #ptau217
            "c2n",
            "roche", #Ab42/40
            "roche", #ptau181
            "roche", #GFAP
            "roche", #NFL
            "fuji",#Ab42/40
            "fuji",#ptau217
            "alz",#ptau217
            "jan",#ptau217
            "quan", #Ab42/40
            "quan", #ptau181
            "quan", #GFAP
            "quan"#NFL
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
           pTau217conc_c2n, #ptau217 - c2n
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



ref_data <- data.frame(group=l_var,var_id=i_var,label=analytes,analyte_code=a_var)
ref_data$company <- recode(ref_data$group,"quan"="Quanterix","alz"="ALZpath","c2n"="C2N",
                           "fuji"="Fujirebio", "jan"="Janssen","roche" = "Roche")
ref_data$fullname <- paste(ref_data$company,ref_data$label)



#  a.	Correlation matrix A (matrix with colors and numbers of unadjusted Spearman correlations)
#i.	Entire cohort
#ii.	Color is based on absolute value (red higher, blue lower)
#iii.	All plasma measures versus one another


order_of_analytes <- rev(c("pt217","pt181","gfap","nfl","ab4240"))

####Creation of heatmap - analyte groups ordered by max absolute rho value
full_heatmap <- create_cor_heatmap_for_project(data.in,ref_data = ref_data)
####Creation of heatmap - analyte groups ordered by custom preference
full_heatmap_reorder <- create_cor_heatmap_for_project_custom(data.in,ref_data = ref_data,reforder = order_of_analytes)

ggsave(paste0(out_plots,"heatmap_full_reorder.png"),
       plot=full_heatmap_reorder,
       device=png,dpi = 1200,width = 6,height=6)



#b.	Correlation matrix B (matrix with colors and numbers of unadjusted Spearman correlations)
#i.	Amyloid PET positive cohort
#ii.	Color is based on absolute value (red higher, blue lower)
#iii.	All plasma measures versus one another
####Creation of heatmap - analyte groups ordered by max absolute rho value
amyloid_positive_heatmap <- create_cor_heatmap_for_project(subset(data.in,CENTILOIDS_10==1),ref_data = ref_data)
ggsave(paste0(out_plots,"heatmap_amyloid_pos.jpg"),
       plot=amyloid_positive_heatmap,
       device=jpeg,width = 6,height=6)

####Creation of heatmap - analyte groups ordered by custom preference
amyloid_positive_heatmap_reorder <- create_cor_heatmap_for_project_custom(subset(data.in,CENTILOIDS_10==1),ref_data = ref_data,reforder = order_of_analytes)

ggsave(paste0(out_plots,"heatmap_amypos_reorder.png"),
       plot=amyloid_positive_heatmap_reorder,
       device=png,dpi = 1200,width = 6,height=6)

#c.	Correlation matrix C (matrix with colors and numbers of unadjusted Spearman correlations)
#i.	Cognitively impaired cohort
#ii.	Color is based on absolute value (red higher, blue lower)
#iii.	All plasma measures versus one another
####Creation of heatmap - analyte groups ordered by max absolute rho value
cdr_positive_heatmap <- create_cor_heatmap_for_project(subset(data.in,CDR_10==1),ref_data = ref_data)
ggsave(paste0(out_plots,"heatmap_cdr_pos.jpg"),
       plot=cdr_positive_heatmap,
       device=jpeg,width = 6,height=6)

####Creation of heatmap - analyte groups ordered by custom preference
cdr_positive_heatmap_reorder <- create_cor_heatmap_for_project_custom(subset(data.in,CDR_10==1),ref_data = ref_data,reforder = order_of_analytes)

ggsave(paste0(out_plots,"heatmap_cdr_pos_reorder.png"),
       plot=cdr_positive_heatmap_reorder,
       device=png,dpi = 1200,width = 6,height=6)


#d.	Correlation matrix D (matrix with colors and numbers of unadjusted Spearman correlations)
#i.	Cognitively unimpaired cohort
#ii.	Color is based on absolute value (red higher, blue lower)
#iii.	All plasma measures versus one another
####Creation of heatmap - analyte groups ordered by max absolute rho value
cdr_negative_heatmap <- create_cor_heatmap_for_project(subset(data.in,CDR_10==0),ref_data = ref_data)
ggsave(paste0(out_plots,"heatmap_cdr_neg.jpg"),
       plot=cdr_negative_heatmap,
       device=jpeg,width = 6,height=6)

####Creation of heatmap - analyte groups ordered by custom preference
cdr_negative_heatmap_reorder <- create_cor_heatmap_for_project_custom(subset(data.in,CDR_10==0),ref_data = ref_data,reforder = order_of_analytes)

ggsave(paste0(out_plots,"heatmap_cdr_neg_reorder.png"),
       plot=cdr_negative_heatmap_reorder,
       device=png,dpi = 1200,width = 6,height=6)


