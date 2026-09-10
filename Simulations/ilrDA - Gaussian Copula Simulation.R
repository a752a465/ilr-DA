#################################################################################################################
### Title: SSDA - Cell Proportion Simulation (Gaussian Copula) V1 ----
### Current Analyst: Alexander Alsup ----
### Last Updated: (07/10/2026)  ----
### Notes:  ----
#################################################################################################################



# PACKAGES -----
packages <- c("writexl","readxl","ggplot2","MASS","tidyverse","dplyr")
lapply(packages, require, character.only = TRUE)
library(DirichletReg)
# library(logratiolasso)
library(Hotelling)
# library(harmonicmeanp)
library(mvtnorm)
###########################################################################
# CLEAR LIST ----
rm(list=ls())
###########################################################################
# GLOBAL ARGUMENTS ----

## Number of Iterations 
Nsims=1000 #; Nsims = 100

## Sample size PER GROUP
sample.size=c(20,40,60,80,100) 
# sample.size=c(10,20,30,40,50,60,70,80,90,100)  # A more dense sample size set
# sample.size=c(20,40,80) # For Bug-fixing

## A scaling factor for the Poisson means
scale=50 

## Models to Run
run_models = c("Poisson"
               ,"Negative_Bino"
               ,"ilrDA"
               ,"Dirch_Reg"
)

# Data Generating function
# Possible arguments: "poisson" ; "neg_bino" ; "gaus_cop"
data_function="gaus_cop"

## Number of Permutations to use for ilr-DA Permutation testing
n_perm=200

## Iteration Message for tracking
it_msg = seq(0,Nsims,by=100)

## Content Saving Argument
# Boolean, save results to xlsx?
save_results = TRUE

###########################################################################
# FILEPATHS ----
Figures_path <- "4. Results/Figures/"
Results_path <- "4. Results/xlsx results/"
# clipr::write_clip(your_data_frame)
###########################################################################
# FUNCTIONS ----
source("3. Code/ilrDA - Functions v0.5.R")
###########################################################################
# SIMPLE SIM -----
names(simparams)
Time_0 <- Sys.time()
Model_Results_Complete <- data.frame()

## 5 Cells - 1 DA - Small Effect ----
set.seed(8895)
sim_config <- list(
  Nsims=Nsims,sample.size=sample.size,scale=scale,
  run_models = run_models,it_msg = it_msg,
  cell_mus = simparams$cell_5_mu_1DA_SmallEffect,
  n_perm = 200,data_function=data_function
)

# Sim Function
Model_Results <- Small_Sim_v1(config=sim_config)
Model_Results_Complete <- bind_rows(Model_Results_Complete,
                                    data.frame(Model_Results) %>% mutate(Effect="Small", CellTypes=5, DA=1))

## 5 Cells - 2 DA - Small Effect ----
set.seed(8895)
sim_config <- list(
  Nsims=Nsims,sample.size=sample.size,scale=scale,
  run_models = run_models,it_msg = it_msg,
  cell_mus = simparams$cell_5_mu_2DA_SmallEffect,
  n_perm = 200,data_function=data_function
)

# Sim Function
Model_Results <- Small_Sim_v1(config=sim_config)
Model_Results_Complete <- bind_rows(Model_Results_Complete,
                                    data.frame(Model_Results) %>% mutate(Effect="Small", CellTypes=5, DA=2))

## 12 Cells - 1 DA - Small Effect ----
set.seed(8895)
sim_config <- list(
  Nsims=Nsims,sample.size=sample.size,scale=scale,
  run_models = run_models,it_msg = it_msg,
  cell_mus = simparams$cell_12_mu_1DA_SmallEffect,
  n_perm = 200,data_function=data_function
)

# Sim Function
Model_Results <- Small_Sim_v1(config=sim_config)
Model_Results_Complete <- bind_rows(Model_Results_Complete,
                                    data.frame(Model_Results) %>% mutate(Effect="Small", CellTypes=12, DA=1))

## 12 Cells - 2 DA - Small Effect ----
set.seed(8895)
sim_config <- list(
  Nsims=Nsims,sample.size=sample.size,scale=scale,
  run_models = run_models,it_msg = it_msg,
  cell_mus = simparams$cell_12_mu_2DA_SmallEffect,
  n_perm = 200,data_function=data_function
)

# Sim Function
Model_Results <- Small_Sim_v1(config=sim_config)
Model_Results_Complete <- bind_rows(Model_Results_Complete,
                                    data.frame(Model_Results) %>% mutate(Effect="Small", CellTypes=12, DA=2))

## Simulation Complete ----

Time_1 <- Sys.time() ; Time_Elapsed = round(difftime(Time_1,Time_0,units="mins"),2)
cat(paste0("\n","Simulation Complete. ",Time_Elapsed," minutes","\n"))

###########################################################################
# Data Saving -----
if(save_results==TRUE){
  write_xlsx(Model_Results_Complete,path=paste0(Results_path,"GaussianCop_",gsub("-","_",Sys.Date()),".xlsx"))
}

# Post-Processing ----
celltype_12 <- names(simparams$cell_12_mu_1DA_SmallEffect$cell_mu0)
celltype_5 <- names(simparams$cell_5_mu_1DA_SmallEffect$cell_mu0)
unique(Model_Results_Complete$Model)

# Renaming
Model_Results_Complete1 <- Model_Results_Complete %>%
  dplyr::filter(Model %in% c("ilrDA","Poisson","Dirch_Reg"))%>%
  mutate(Model=case_when(Model=="Poisson"~"Poisson GLM",
                         Model=="ilrDA"~"ilrDA",
                         Model=="SSDA"~"SSDA",
                         Model=="LR_Lasso"~"Logratio Lasso",
                         Model=="Dirch_Reg"~"Dirichlet Regression",
                         TRUE~NA))

# 12 Cell Type 1 DA
cell12da1 <- Model_Results_Complete1%>%
  dplyr::select(CellTypes,DA,Effect,N,mean_FDR,Model,all_of(celltype_12))%>%
  arrange(DA,Model,Effect,Model,N)%>%
  dplyr::filter(DA==1 & CellTypes==12)

# 12 Cell Type 2 DA
cell12da2 <- Model_Results_Complete1%>%
  dplyr::select(CellTypes,DA,Effect,N,mean_FDR,Model,all_of(celltype_12))%>%
  arrange(DA,Model,Effect,Model,N)%>%
  dplyr::filter(DA==2 & CellTypes==12)

# 5 Cell Type 1 DA
cell5da1 <- Model_Results_Complete1 %>%
  dplyr::select(CellTypes,DA,Effect,N,mean_FDR,Model,all_of(celltype_5))%>%
  arrange(DA,Model,Effect,Model,N)%>%
  dplyr::filter(DA==1 & CellTypes==5)

# 5 Cell Type 2 DA
cell5da2 <- Model_Results_Complete1 %>%
  dplyr::select(CellTypes,DA,Effect,N,mean_FDR,Model,all_of(celltype_5))%>%
  arrange(DA,Model,Effect,Model,N)%>%
  dplyr::filter(DA==2 & CellTypes==5)

###########################################################################