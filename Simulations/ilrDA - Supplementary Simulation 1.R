#################################################################################################################
### Title: ilr-DA: Cell Proportion Simulation (Negative Binomial) V1 ----
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
Nsims=100 #; Nsims = 100

## Sample size PER GROUP
sample.size=c(40,80,120,160,200)
# sample.size=c(10,20,30,40,50,60,70,80,90,100)  # A more dense sample size set
# sample.size=c(20,40,80) # For Bug-fixing

## A scaling factor for the Poisson means
scale=50 

## Models to Run
run_models = c("Negative_Bino"
               ,"ilrDA"
               #,"Poisson"
               #,"Dirch_Reg"
)

# Data Generating function
  # Possible arguments: "poisson" ; "neg_bino" ; "gaus_cop"
data_function="neg_bino"

## Number of Permutations to use for ilr-DA Permutation testing
n_perm=200

## Iteration Message for tracking
it_msg = seq((Nsims/10),Nsims,length=10)

## Content Saving Argument
# Boolean, save results to xlsx?
save_results = TRUE
save_plots = FALSE

## SIMULATION PARAMS -----
suppsim1.simparams <- list(
  #### 5 Cell Types : 1 DA : Small Effect Size (Neut)  ----
  cell_5_mu_1DA_SmallEffect = list(
    cell_mu0 = c("Neutrophils"=2.8,"Lymphocytes"=2.4,"Monocytes"=0.44,"Eosinophils"=0.30,"Basophils"=0.04),
    cell_mu1 = c("Neutrophils"=-0.03,"Lymphocytes"=0,"Monocytes"=0,"Eosinophils"=0,"Basophils"=0)
  ),
  #### 5 Cell Types : 2 DA : Small Effect Size (Neut & Mono)  ----
  cell_5_mu_2DA_SmallEffect = list(
    cell_mu0 = c("Neutrophils"=2.8,"Lymphocytes"=2.4,"Monocytes"=0.44,"Eosinophils"=0.30,"Basophils"=0.04),
    cell_mu1 = c("Neutrophils"=-0.03,"Lymphocytes"=0,"Monocytes"=+0.005,"Eosinophils"=0,"Basophils"=0)
  )
)


###########################################################################
# FILEPATHS ----
Figures_path <- "4. Results/Figures/"
Results_path <- "4. Results/xlsx results/"
# clipr::write_clip(your_data_frame)
###########################################################################
# FUNCTIONS ----

## Main Functions ----
source("3. Code/ilrDA - Functions v0.5.R")

## Supplementary Functions ----
### Supplementary Data Generation with Continuous X ----
SuppSim1_data_generation = function(cell_mu0,cell_mu1,X,N,type=data_function,size_nbinom=25){
  # Start of Function
  cell_mu0_matrix = matrix(rep(cell_mu0, each=N), nrow=N, ncol=length(cell_mu0))
  cell_mu1_matrix = matrix(rep(cell_mu1, each=N), nrow=N, ncol=length(cell_mu1))
  
  cell_mu_matrix = cell_mu0_matrix + outer(X, cell_mu1)
  ncell = length(cell_mu0)
  # Cell Counts
  counts <- matrix(NA, nrow=N, ncol=length(cell_mu0))
  if(type=="poisson"){
    # Poisson Distribution
    for(i in 1:N){
      counts[i,] <- rpois(length(cell_mu0), lambda=cell_mu_matrix[i,]) + 0.1
    }
  } else if(type=="neg_bino"){
    # Negative Binomial Distribution
    for(i in 1:N){
      counts[i,] <- rnbinom(n=length(cell_mu0), mu=cell_mu_matrix[i,], size=size_nbinom) + 0.1
    }
  }
  # Cell Data Bound
  data.cells <- counts; colnames(data.cells) <- names(cell_mu0)
  # Cell Proportions
  data.proportions <- data.cells/rowSums(data.cells)
  return(list("Counts"=data.cells,"Proportions"=data.proportions))
}

### Supplementary Sim Wrapper  -----
Small_SupSim_v1 <- function(config=config){
  with(config,{
    
    ## Simulation Parameters ----
    cell_mu0 = cell_mus$cell_mu0*scale ;  cell_mu1 = cell_mus$cell_mu1*scale
    power_cell = which(cell_mu1 != 0)
    fdr_cell = which(cell_mu1 == 0)
    ncell=length(cell_mu0)
    celltype=names(cell_mu0)
    
    ## Compute FDR Function ----
    # After each model's result_i is computed, before appending:
    compute_FDR_i <- function(result_i, power_cell, ncell){
      declared <- as.integer(result_i[1:ncell])
      R <- sum(declared,na.rm=TRUE)
      if(R == 0) return(0)
      V <- sum(declared,na.rm=TRUE) - sum(declared[power_cell],na.rm=TRUE)  # false discoveries
      V / R
    }
    
    ## Compute F-score Function ----
    compute_fscore <- function(result_i, power_cell, celltype){
      # declared: binary vector (0/1 or TRUE/FALSE) of length K
      # true_DA_idx: indices of truly DA cell types
      declared=result_i
      true_DA <- as.integer(seq_along(celltype) %in% power_cell)
      
      TP <- sum(declared == 1 & true_DA == 1)  # correctly declared DA
      FP <- sum(declared == 1 & true_DA == 0)  # falsely declared DA
      FN <- sum(declared == 0 & true_DA == 1)  # missed DA
      
      precision <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
      sensitivity <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
      
      fscore <- ifelse(is.na(precision) | is.na(sensitivity) | 
                         (precision + sensitivity) == 0, NA,
                       2 * (precision * sensitivity) / (precision + sensitivity))
      
      return(fscore)
    }
    
    ## Preparing Empty Results Frames ----
    Search_Candidates = c()
    Model_Results = c()
    Time_0 <- Sys.time()
    
    ## Simulation ----
    for(k in 1:length(sample.size)){
      ### Outer Loop Specification ----
      
      # K Level Results
      Model_Results_k = c()
      Search_Candidates_k = c()
      
      # K-Level Parameters
      N_k=sample.size[k]
      cat(paste0("**Outer Loop** Sample Size: ",N_k,"\n"))
      cat(paste0("Iteration: "))
      ids = paste0("P_",1:N_k)
      
      ### Inner Loop ----
      for(i in 1:Nsims){
        
        # Continuous Predictor
        #X=floor(runif(n=N_k,18,55))
        X <- round(rnorm(N_k, mean=40, sd=10))
        X <- pmax(pmin(X, 65), 18)  # truncate to [18,65]
        
        # Data Generation
        toy_data <- SuppSim1_data_generation(cell_mu0,cell_mu1,X,N_k,type=data_function,size_nbinom=25)
        props <- toy_data$Proportions
        counts <- toy_data$Counts
        
        ### Poisson Counts Check  ----
        if("Poisson" %in% sim_config$run_models){
          result_i <- Poisson_Check(counts,X,p_correction="BH")
          Corrected_Results=result_i$Corrected_Results
          Uncorrected_Results=result_i$Uncorrected_Results
          # Corrected Results
          Result_search_i <- c(Corrected_Results,
                               compute_FDR_i(Corrected_Results, power_cell, ncell),
                               compute_fscore(Corrected_Results, power_cell, celltype),
                               "Poisson_Corrected",
                               0)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
          # Uncorrected Results
          Result_search_i <- c(Uncorrected_Results,
                               compute_FDR_i(Uncorrected_Results, power_cell, ncell),
                               compute_fscore(Uncorrected_Results, power_cell, celltype),
                               "Poisson",
                               0)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
        }
        
        ### ilrDA: ILR-basis permutation search ----
        if("ilrDA" %in% sim_config$run_models){
          result_i <- ilr_DA_General(props, celltype, X, covars=NULL,perms=500,test_function=sim_config$test_function,omnibus=FALSE,p_correction="BH", perm_test=TRUE)
          Corrected_Results=result_i$Corrected_Results
          Uncorrected_Results=result_i$Uncorrected_Results
          # Corrected Results
          Result_search_i <- c(Corrected_Results,
                               compute_FDR_i(Corrected_Results, power_cell, ncell),
                               compute_fscore(Corrected_Results, power_cell, celltype),
                               "ilrDA_Corrected",
                               0)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
          # Uncorrected Results
          Result_search_i <- c(Uncorrected_Results,
                               compute_FDR_i(Uncorrected_Results, power_cell, ncell),
                               compute_fscore(Uncorrected_Results, power_cell, celltype),
                               "ilrDA",
                               0)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
        }
        
        ### Dirichlet Regression: Common Model ----
        if("Dirch_Reg" %in% sim_config$run_models){
          result_i <- Dirichlet_DA(props=props, X=X, celltype=celltype,p_correction="BH")
          Corrected_Results=result_i$Corrected_Results
          Uncorrected_Results=result_i$Uncorrected_Results
          # Corrected Results
          Result_search_i <- c(Corrected_Results,
                               compute_FDR_i(Corrected_Results, power_cell, ncell),
                               compute_fscore(Corrected_Results, power_cell, celltype),
                               "Dirch_Reg_Corrected",
                               result_i$error)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
          # Uncorrected Results
          Result_search_i <- c(Uncorrected_Results,
                               compute_FDR_i(Uncorrected_Results, power_cell, ncell),
                               compute_fscore(Uncorrected_Results, power_cell, celltype),
                               "Dirch_Reg",
                               result_i$error)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
        }
        
        ### Negative Binomial Counts Check ----
        if("Negative_Bino" %in% sim_config$run_models){
          result_i <- NB_Check(X, counts, celltype,p_correction="BH")
          Corrected_Results=result_i$Corrected_Results
          Uncorrected_Results=result_i$Uncorrected_Results
          # Corrected Results
          Result_search_i <- c(Corrected_Results,
                               compute_FDR_i(Corrected_Results, power_cell, ncell),
                               compute_fscore(Corrected_Results, power_cell, celltype),
                               "Negative_Binom_Corrected",
                               0)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
          # Uncorrected Results
          Result_search_i <- c(Uncorrected_Results,
                               compute_FDR_i(Uncorrected_Results, power_cell, ncell),
                               compute_fscore(Uncorrected_Results, power_cell, celltype),
                               "Negative_Binom",
                               0)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
        }
        
        ### ilrDA: Boostrap ----
        if("ilrDA_Boot" %in% sim_config$run_models){
          # Full SBP W Search
          Search <- ilr_DA_Boot(props=props, celltype=celltype,B_num=200, X=X,search_type="HMP")
          Result_search_i <- as.integer(celltype %in% Search$DA_cells)
          FDR_i <- compute_FDR_i(Result_search_i, power_cell, ncell)
          fscore_i <- compute_fscore(Result_search_i, power_cell, celltype)
          Result_search_i = c(Result_search_i,FDR_i,fscore_i,"ilrDA_boot",0)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
        }
        
        ### ilrDA GLM: ILR-basis permutation search using GLM parameterization ----
        if("ilrDA_GLM" %in% sim_config$run_models){
          # Full SBP W Search
          Search <- ilr_DA_GLM(props,celltype,X,covars=NULL,general_formula="~X",perm_test=TRUE,search_type="HMP",perms=200)
          Result_search_i <- as.integer(celltype %in% Search$DA_cells)
          FDR_i <- compute_FDR_i(Result_search_i, power_cell, ncell)
          fscore_i <- compute_fscore(Result_search_i, power_cell, celltype)
          Result_search_i = c(Result_search_i,FDR_i,fscore_i,"ilrDA_GLM",0)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
        }
        
        ### Logratio Lasso -----
        if("LR Lasso" %in% sim_config$run_models){
          # Full SBP W Search
          Search <- logratiolasso_DA(props, X, celltype)
          Result_search_i <- as.integer(celltype %in% Search)
          FDR_i <- compute_FDR_i(Result_search_i, power_cell, ncell)
          fscore_i <- compute_fscore(Result_search_i, power_cell, celltype)
          Result_search_i = c(Result_search_i,FDR_i,fscore_i,"LR_Lasso",0)
          Model_Results_k = rbind(Model_Results_k, Result_search_i)
        }
        
        ### SSDA ---- 
        if("SSDA" %in% sim_config$run_models){
          SSDA_result_i <- Fisher_Test_Full(ncell=ncell,data.proportions=props,group=X,N=N,cell.type=celltype,
                                            power_cell=power_cell,fdr_cell=fdr_cell,FDR_control="BH",FDR_control_type="Sub Hypothesis")$Model.Test
          FDR_i <- compute_FDR_i(SSDA_result_i, power_cell, ncell)
          fscore_i <- compute_fscore(SSDA_result_i, power_cell, celltype)
          SSDA_result_i = c(SSDA_result_i,FDR_i,fscore_i,"SSDA",0)
          Model_Results_k <- rbind(Model_Results_k, SSDA_result_i)
        }
        
        ### Dirichlet Regression: Alternate Model ----
        if("Dirch_Reg_Alt" %in% sim_config$run_models){
          result_i <- Dirichlet_DA_Alt(X, props, celltype,model_type="Alt")
          FDR_i <- compute_FDR_i(result_i$result, power_cell, ncell)
          fscore_i <- compute_fscore(result_i$result, power_cell, celltype)
          result_i = c(result_i$result,FDR_i,fscore_i,"Dirch_Reg_Alt",result_i$error)
          Model_Results_k <- rbind(Model_Results_k, result_i)
        }
        
        ### ANCOM-BC2 -----
        if("ANCOM_BC2" %in% sim_config$run_models){
          result_i <- ANCOMBC2_test(props,X)
          FDR_i <- compute_FDR_i(result_i, power_cell, ncell)
          fscore_i <- compute_fscore(result_i, power_cell, celltype)
          result_i = c(result_i,FDR_i,fscore_i,"ANCOM_BC2",0)
          Model_Results_k <- rbind(Model_Results_k, result_i)
        }
        
        ### Message  ---- 
        if(i %in% it_msg){
          cat(paste0(i,"|"))
        }
        
      } # Inner Loop End
      
      cat(paste0("\n"))
      
      ### Appending New Results  ---- 
      Model_Results_k = Model_Results_k[!is.na(Model_Results_k[,1]),]
      Model_Results_k = cbind(Model_Results_k,N_k)
      Model_Results = rbind(Model_Results,Model_Results_k)
      
    } # Outer Loop End
    
    ## Results Prep ----
    Time_1 <- Sys.time() ; Time_Elapsed = round(difftime(Time_1,Time_0,units="mins"),2)
    cat(paste0("\n","Simulation Complete. ",Time_Elapsed," minutes","\n"))
    # Cleaning and Finalizing Results
    colnames(Model_Results)=c(celltype,"FDR","F_score","Model","Error_Message","N")
    
    # Calculating Mean FDR
    Mean_FDR <- data.frame(Model_Results) %>% 
      mutate(FDR=as.numeric(FDR),
             FDR=case_when(is.na(FDR)~0,TRUE~FDR),
             Error_Message=as.numeric(Error_Message),
             N=as.numeric(N),
             F_score=case_when(is.na(F_score)~"0",TRUE~F_score),
             F_score=as.numeric(F_score))%>%
      group_by(Model,N)%>%
      dplyr::summarise(mean_FDR = mean(FDR, na.rm=TRUE),
                       mean_error = mean(Error_Message, na.rm=TRUE),
                       F_score = mean(F_score, na.rm=TRUE),
                       .groups = "drop")
    # Cleaning Results
    Model_Results_final <- data.frame(Model_Results) %>%
      mutate(across(all_of(celltype),as.numeric),
             FDR = as.numeric(FDR))%>%
      pivot_longer(cols=all_of(celltype),names_to="celltype",values_to="detection")%>%
      group_by(N,celltype,Model)%>%
      dplyr::summarise(detection=mean(detection,na.rm=TRUE),.groups = "drop_last")%>%
      ungroup()%>%
      pivot_wider(names_from="celltype",values_from = "detection")%>%
      mutate(N=as.numeric(N))%>%
      arrange(N,Model)
    # Joining to Mean FDR
    Model_Results_final <- Model_Results_final %>%
      left_join(Mean_FDR,by=c("Model","N"))
    ## Return ----
    return("Model_Results_Final"=Model_Results_final)
  })
}

### ilr-DA Test Function ----
ilr_DA_test_fn = function(lr_vals, X, covars=NULL){
  # Define Analysis Dataframe
  data <- data.frame(lr=lr_vals, X=X)
  if(!is.null(covars)) data <- cbind(data, covars)
  # Fit Model
  lm <- summary(glm(lr ~ ., data=data, family=gaussian()))
  coef_row <- coef(lm)["X",]
  list(t_stat=coef_row[3], p_value=coef_row[4])
}
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
  cell_mus = suppsim1.simparams$cell_5_mu_1DA_SmallEffect,
  n_perm = 200,data_function=data_function,
  test_function=ilr_DA_test_fn
)

# Sim Function
Model_Results <- Small_SupSim_v1(config=sim_config)
Model_Results_Complete <- bind_rows(Model_Results_Complete,
                                    data.frame(Model_Results) %>% mutate(Effect="Small", CellTypes=5, DA=1))

## 5 Cells - 2 DA - Small Effect ----
set.seed(8895)
sim_config <- list(
  Nsims=Nsims,sample.size=sample.size,scale=scale,
  run_models = run_models,it_msg = it_msg,
  cell_mus = suppsim1.simparams$cell_5_mu_2DA_SmallEffect,
  n_perm = 200,data_function=data_function,
  test_function=ilr_DA_test_fn
)

# Sim Function
Model_Results <- Small_SupSim_v1(config=sim_config)
Model_Results_Complete <- bind_rows(Model_Results_Complete,
                                    data.frame(Model_Results) %>% mutate(Effect="Small", CellTypes=5, DA=2))
## Simulation Complete ----

Time_1 <- Sys.time() ; Time_Elapsed = round(difftime(Time_1,Time_0,units="mins"),2)
cat(paste0("\n","Simulation Complete. ",Time_Elapsed," minutes","\n"))

###########################################################################
# Data Saving -----
if(save_results==TRUE){
  write_xlsx(Model_Results_Complete,path=paste0(Results_path,"SuppSim1_",gsub("-","_",Sys.Date()),".xlsx"))
}

# Post-Processing ----
celltype_5 <- names(simparams$cell_5_mu_1DA_SmallEffect$cell_mu0)
unique(Model_Results_Complete$Model)

# Renaming
Model_Results_Complete1 <- Model_Results_Complete %>%
  dplyr::filter(Model %in% c("Negative_Binom","Negative_Binom_Corrected","ilrDA","ilrDA_Corrected"))%>%
  mutate(Model=case_when(Model=="Negative_Binom"~"Negative Binomial GLM",
                         Model=="ilrDA"~"ilrDA",
                         Model=="Negative_Binom_Corrected"~"Negative Binomial GLM with Correction",
                         Model=="ilrDA_Corrected"~"ilr-DA with Correction",
                         TRUE~NA))

# 5 Cell Type 1 DA
cell5da1 <- Model_Results_Complete1 %>%
  dplyr::select(CellTypes,DA,Effect,N,mean_FDR,F_score,Model,all_of(celltype_5))%>%
  arrange(DA,Model,Effect,Model,N)%>%
  dplyr::filter(DA==1 & CellTypes==5)

# 5 Cell Type 2 DA
cell5da2 <- Model_Results_Complete1 %>%
  dplyr::select(CellTypes,DA,Effect,N,mean_FDR,F_score,Model,all_of(celltype_5))%>%
  arrange(DA,Model,Effect,Model,N)%>%
  dplyr::filter(DA==2 & CellTypes==5)

###########################################################################
# VISUALIZATIONS -----
## Overall prep ----

# Figure Paths and Subtitles
models_selected = unique(Model_Results_Complete1$Model)
sim_type <- c("Negative Binomial"="Negative_Bino","Poisson"="Poisson","Gaussian Copula"="Gauss")
Biplot_title = paste0("Supplementary Simulation 1 Results: ")
subtitle_figure = paste0(sim_type," Data Generating Function")
#
Figures_path_modified = paste0(Figures_path,sim_type,"/")

# Breaks for Scales
N_breaks = unique(Model_Results_Complete1$N)
y_breaks=seq(from=0,to=1,by=0.25)

# Cell Names
rename_vec <- c("Neutrophils" = "Neutrophils",
                "Lymphocytes" = "Lymphocytes",
                "Eosinophils" = "Eosinophils",
                "Basophils" = "Basophils",
                "CD4_Naive" = "CD4 Naive",
                "CD4_Mem" = "CD4 Memory",
                "CD8_Naive" = "CD8 Naive",
                "CD8_Mem" = "CD8 Memory",
                "B_Naive" =  "B Naive",
                "B_Mem" = "B Memory",
                "Treg" = "T Regulatory",
                "NK" = "Natural Killer")

# Discrete Aesthetics
scale_shape_manual_vals = c("ilrDA"=17,"ilr-DA with Correction"=17,
                            "Poisson GLM"=16,"Poisson GLM (MC)"=16,
                            "Dirichlet Regression - Common"=15,"Dirichlet Regression - Common (MC)"=15,
                            "Negative Binomial GLM"=12,"Negative Binomial GLM with Correction"=12)
# Color of the Points and Lines
scale_color_manual_vals = c("ilrDA"="#CD9600","ilr-DA with Correction"="#CD9600",
                            "Poisson GLM"="#4075d6","Poisson GLM (MC)"="#4075d6",
                            "Dirichlet Regression - Common"="darkorchid4","Dirichlet Regression - Common (MC)"="darkorchid4",
                            "Negative Binomial GLM"="#8B0000","Negative Binomial GLM with Correction"="#8B0000")

# Linetype
scale_linetype_manual_vals = c("ilrDA"="solid","ilr-DA with Correction"="dashed",
                               "Poisson GLM"="solid","Poisson GLM (MC)"="dashed",
                               "Dirichlet Regression - Common"="solid","Dirichlet Regression - Common (MC)"="dashed",
                               "Negative Binomial GLM"="solid","Negative Binomial GLM with Correction"="dashed")

# Linewidth
scale_linewidth_manual_vals = c("ilrDA"=1.2,"ilr-DA with Correction"=1.2,
                                "Poisson GLM"=0.5,"Poisson GLM (MC)"=0.5,
                                "Dirichlet Regression - Common"=0.5,"Dirichlet Regression - Common (MC)"=0.5,
                                "Negative Binomial GLM"=0.5,"Negative Binomial GLM with Correction"=0.5)

# Point Size
scale_size_manual_vals = c("ilrDA"=2,"ilr-DA with Correction"=2,
                           "Poisson GLM"=1,"Poisson GLM (MC)"=1,
                           "Dirichlet Regression - Common"=1,"Dirichlet Regression - Common (MC)"=1,
                           "Negative Binomial GLM"=1,"Negative Binomial GLM with Correction"=1)

### Bi-Plot Function ----
bi_plot <- function(data=list(cell12da1,cell12da2)){
  p <- bind_rows(data)%>%
    dplyr::select(CellTypes,DA,N,Model,F_score,mean_FDR)%>%
    dplyr::filter(Model %in% models_selected) %>%
    mutate(Model=factor(Model,levels=models_selected))%>%
    pivot_longer(cols=c("F_score","mean_FDR"),names_to="Metric",values_to="Value")%>%
    mutate(Metric=case_when(Metric=="F_score"~"F1 Score",Metric=="mean_FDR"~"FDR",TRUE~NA),
           Metric=factor(Metric,levels=c("F1 Score","FDR")),
           `DA Elements`=DA)%>%
    ggplot(aes(x=N,y=Value,linetype=Model,shape=Model,color=Model))+
    geom_line(linewidth=1,alpha=0.7)+
    geom_point(size=2.5,alpha=1)+
    # Labels
    labs(title=Biplot_title,
         subtitle = title_suffix,
         x="Total Sample Size",
         y=element_blank())+
    # Color & Shape of Lines by Model
    scale_shape_manual(values=scale_shape_manual_vals)+
    scale_color_manual(values=scale_color_manual_vals)+
    scale_linetype_manual(values=scale_linetype_manual_vals)+
    scale_linewidth_manual(values=scale_linewidth_manual_vals,guide="none")+
    scale_size_manual(values=scale_size_manual_vals,guide="none")+
    # Scale of Axes
    scale_x_continuous(breaks=N_breaks)+
    scale_y_continuous(limits=c(0,1),breaks=y_breaks)+
    # Wrapping
    facet_grid(`DA Elements`~Metric,scales="fixed",labeller=label_both)+
    # Theme
    theme_bw()+
    theme(legend.position="none",
          strip.text=element_text(size=9),
          #axis.text.x = element_blank(),
          panel.grid.minor = element_blank()
    )
  #
  if(save_plots==TRUE){
    file_name <-paste0(Figures_path_modified,results_type,figname_prefix,".png")
    ggsave(file_name,plot=p, width = 8, height = 6, units="in",dpi=400)
  }
  p
}
################################
## Biplot ----
# Title Suffix
title_suffix <- "Negative Binomial DGF, 100 Simulations, 5 Cell Types"
figname_prefix <- "_Suppsim1_Cell5_BiPlot"
bi_plot(data=list(cell5da1,cell5da2))