#################################################################################################################
### Title: SSDA Update - Functions ----
### Current Analyst: Alexander Alsup ----
### Last Updated: (06/12/2026)  ----
### Notes:  ----
#################################################################################################################

# Packages - (ANCOM) -----

if("ANCOM_BC2" %in% run_models){
  library(ANCOMBC)
  library(phyloseq)
  library(microbiome)
}

# Harmonic Mean P-value
library(harmonicmeanp)

#################################################################################################################

# Paths & Parameters -----


#################################################################################################################

# Cell Type Mus ----
###################################################
## TOY DATA -----

toyparams <- list(
#### 5 Cell Types: 3 DA ----
cell_5_mu_3DA = list(
  cell_mu0 = c("Neutrophils"=600,"Lymphocytes"=300,"Monocytes"=50,"Eosinophils"=30,"Basophils"=20),
  cell_mu1 = c("Neutrophils"=400,"Lymphocytes"=350,"Monocytes"=50,"Eosinophils"=30,"Basophils"=10)
),

#### 5 Cell Types: 2 DA ----
cell_5_mu_2DA = list(
  cell_mu0 = c("Neutrophils"=600,"Lymphocytes"=300,"Monocytes"=50,"Eosinophils"=30,"Basophils"=20),
  cell_mu1 = c("Neutrophils"=400,"Lymphocytes"=300,"Monocytes"=50,"Eosinophils"=30,"Basophils"=10)
),


#### 5 Cell Types: 1 DA ----
cell_5_mu_1DA = list(
  cell_mu0 = c("Neutrophils"=600,"Lymphocytes"=300,"Monocytes"=50,"Eosinophils"=30,"Basophils"=20),
  cell_mu1 = c("Neutrophils"=400,"Lymphocytes"=300,"Monocytes"=50,"Eosinophils"=30,"Basophils"=20)
),


#### 5 Cell Types: 0 DA ----
cell_5_mu_0DA = list(
  cell_mu0 = c("Neutrophils"=600,"Lymphocytes"=300,"Monocytes"=50,"Eosinophils"=30,"Basophils"=20),
  cell_mu1 = c("Neutrophils"=600,"Lymphocytes"=300,"Monocytes"=50,"Eosinophils"=30,"Basophils"=20)
),

#### 12 Cell Types (2 DA) ----
cell_12_mu_2DA = list(
  cell_mu0 = c("Neutrophils"=600, "Lymphocytes"=300, "Monocytes"=50, 
               "Eosinophils"=30, "Basophils"=20, "NK_cells"=40,
               "CD4T"=200, "CD8T"=100, "Bcells"=80, 
               "Tregs"=15, "pDC"=10, "mDC"=12),
  cell_mu1 = c("Neutrophils"=400, "Lymphocytes"=400, "Monocytes"=50, 
               "Eosinophils"=30, "Basophils"=20, "NK_cells"=40,
               "CD4T"=200, "CD8T"=100, "Bcells"=80, 
               "Tregs"=15, "pDC"=10, "mDC"=12)
),

#### 20 Cell Type Mu (1 DA) ----
cell_20_mu_1DA = list(
  cell_mu0 = c("Neutrophils"=600, "Lymphocytes"=300, "Monocytes"=50,
               "Eosinophils"=30, "Basophils"=20, "NK_cells"=40,
               "CD4T"=200, "CD8T"=100, "Bcells"=80,
               "Tregs"=15, "pDC"=10, "mDC"=12,
               "Classical_Mono"=35, "NonClassical_Mono"=15,
               "NKT"=12, "Naive_CD4"=80, "Memory_CD4"=120,
               "Naive_CD8"=40, "Memory_CD8"=60, "Plasma"=8)
  ,
  cell_mu1 = c("Neutrophils"=400, "Lymphocytes"=300, "Monocytes"=50,
               "Eosinophils"=30, "Basophils"=20, "NK_cells"=40,
               "CD4T"=200, "CD8T"=100, "Bcells"=80,
               "Tregs"=15, "pDC"=10, "mDC"=12,
               "Classical_Mono"=35, "NonClassical_Mono"=15,
               "NKT"=12, "Naive_CD4"=80, "Memory_CD4"=120,
               "Naive_CD8"=40, "Memory_CD8"=60, "Plasma"=8)
),

#### 20 Cell Type Mu (0 DA) ----
cell_20_mu_0DA = list(
  cell_mu0 = c("Neutrophils"=600, "Lymphocytes"=300, "Monocytes"=50,
               "Eosinophils"=30, "Basophils"=20, "NK_cells"=40,
               "CD4T"=200, "CD8T"=100, "Bcells"=80,
               "Tregs"=15, "pDC"=10, "mDC"=12,
               "Classical_Mono"=35, "NonClassical_Mono"=15,
               "NKT"=12, "Naive_CD4"=80, "Memory_CD4"=120,
               "Naive_CD8"=40, "Memory_CD8"=60, "Plasma"=8)
  ,
  cell_mu1 = c("Neutrophils"=600, "Lymphocytes"=300, "Monocytes"=50,
               "Eosinophils"=30, "Basophils"=20, "NK_cells"=40,
               "CD4T"=200, "CD8T"=100, "Bcells"=80,
               "Tregs"=15, "pDC"=10, "mDC"=12,
               "Classical_Mono"=35, "NonClassical_Mono"=15,
               "NKT"=12, "Naive_CD4"=80, "Memory_CD4"=120,
               "Naive_CD8"=40, "Memory_CD8"=60, "Plasma"=8)
)
) # End toy params list
###################################################
## SIM DATA ----
simparams <- list(
  ### 5 Cell Types ----
  #### 5 Cell Types : 0 DA : Null Effect Size  ----
  cell_5_mu_0DA_NullEffect = list(
    cell_mu0 = c("Neutrophils"=2.8,"Lymphocytes"=2.4,"Monocytes"=0.44,"Eosinophils"=0.30,"Basophils"=0.04),
    cell_mu1 = c("Neutrophils"=2.8,"Lymphocytes"=2.4,"Monocytes"=0.44,"Eosinophils"=0.30,"Basophils"=0.04)
  ),
  #### 5 Cell Types : 1 DA : Small Effect Size (Neut)  ----
  cell_5_mu_1DA_SmallEffect = list(
    cell_mu0 = c("Neutrophils"=2.8,"Lymphocytes"=2.4,"Monocytes"=0.44,"Eosinophils"=0.30,"Basophils"=0.04),
    cell_mu1 = c("Neutrophils"=2.3,"Lymphocytes"=2.4,"Monocytes"=0.44,"Eosinophils"=0.30,"Basophils"=0.04)
  ),
  #### 5 Cell Types : 2 DA : Small Effect Size (Neut & Mono)  ----
  cell_5_mu_2DA_SmallEffect = list(
    cell_mu0 = c("Neutrophils"=2.8,"Lymphocytes"=2.4,"Monocytes"=0.44,"Eosinophils"=0.30,"Basophils"=0.04),
    cell_mu1 = c("Neutrophils"=1.3,"Lymphocytes"=2.4,"Monocytes"=0.64,"Eosinophils"=0.30,"Basophils"=0.04)
  ),
  
  #### 5 Cell Types: 2 DA : Small Effect Size (Neut & Baso) ----
  cell_5_mu_2DA_SmallEffect_BasoNeut = list(
    cell_mu0 = c("Neutrophils"=2.8,"Lymphocytes"=2.4,"Monocytes"=0.44,"Eosinophils"=0.30,"Basophils"=0.04),
    cell_mu1 = c("Neutrophils"=2.5,"Lymphocytes"=2.4,"Monocytes"=0.44,"Eosinophils"=0.30,"Basophils"=0.065)
  ),
  
  ### 12 Cell Types ----
  #### 12 Cell Types : 1 DA : Small Effect Size (Neut) ----
  cell_12_mu_1DA_SmallEffect = list(
    cell_mu0 = c("Neutrophils"=2.74,"Monocytes"=0.27,"Eosinophils"=0.16,"Basophils"=0.055,"CD4_Naive"=0.65,"CD4_Mem"=0.38,"CD8_Naive"=0.16,"CD8_Mem"=0.22,"B_Naive"=0.22,"B_Mem"=0.16,"Treg"=0.11,"NK"=0.33),
    cell_mu1 = c("Neutrophils"=2.29,"Monocytes"=0.27,"Eosinophils"=0.16,"Basophils"=0.055,"CD4_Naive"=0.65,"CD4_Mem"=0.38,"CD8_Naive"=0.16,"CD8_Mem"=0.22,"B_Naive"=0.22,"B_Mem"=0.16,"Treg"=0.11,"NK"=0.33)
  ),
  #### 12 Cell Types : 2 DA : Small Effect Size (Neut & Mono) ----
  cell_12_mu_2DA_SmallEffect = list(
    cell_mu0 = c("Neutrophils"=2.74,"Monocytes"=0.27,"Eosinophils"=0.16,"Basophils"=0.055,"CD4_Naive"=0.65,"CD4_Mem"=0.38,"CD8_Naive"=0.16,"CD8_Mem"=0.22,"B_Naive"=0.22,"B_Mem"=0.16,"Treg"=0.11,"NK"=0.33),
    cell_mu1 = c("Neutrophils"=2.29,"Monocytes"=0.37,"Eosinophils"=0.16,"Basophils"=0.055,"CD4_Naive"=0.65,"CD4_Mem"=0.38,"CD8_Naive"=0.16,"CD8_Mem"=0.22,"B_Naive"=0.22,"B_Mem"=0.16,"Treg"=0.11,"NK"=0.33)
  ),
  
  #### 12 Cell Types: 2 DA : Small Effect Size (Neut & Baso) ----
  cell_12_mu_2DA_SmallEffect_BasoNeut = list(
    cell_mu0 = c("Neutrophils"=2.74,"Monocytes"=0.27,"Eosinophils"=0.16,"Basophils"=0.055,"CD4_Naive"=0.65,"CD4_Mem"=0.38,"CD8_Naive"=0.16,"CD8_Mem"=0.22,"B_Naive"=0.22,"B_Mem"=0.16,"Treg"=0.11,"NK"=0.33),
    cell_mu1 = c("Neutrophils"=2.54,"Monocytes"=0.27,"Eosinophils"=0.16,"Basophils"=0.080,"CD4_Naive"=0.65,"CD4_Mem"=0.38,"CD8_Naive"=0.16,"CD8_Mem"=0.22,"B_Naive"=0.22,"B_Mem"=0.16,"Treg"=0.11,"NK"=0.33)
  ),
  
  #### 12 Cell Types : 3 DA : Small Effect Size ----
  cell_12_mu_3DA_SmallEffect = list(
    cell_mu0 = c("Neutrophils"=2.74,"Monocytes"=0.27,"Eosinophils"=0.16,"Basophils"=0.055,"CD4_Naive"=0.65,"CD4_Mem"=0.38,"CD8_Naive"=0.16,"CD8_Mem"=0.22,"B_Naive"=0.22,"B_Mem"=0.16,"Treg"=0.11,"NK"=0.33),
    cell_mu1 = c("Neutrophils"=2.29,"Monocytes"=0.37,"Eosinophils"=0.16,"Basophils"=0.055,"CD4_Naive"=0.5,"CD4_Mem"=0.38,"CD8_Naive"=0.16,"CD8_Mem"=0.22,"B_Naive"=0.22,"B_Mem"=0.16,"Treg"=0.11,"NK"=0.33)
  )
  
  ### 20 Cell Types ----
  
  ) # End Sim data List



############################################################
# Functions ----

## Data Generation ----
data.generation <- function(ncell=ncell,cell_mu0=cell_mu0,cell_mu1=cell_mu1,N1=N1,N0=N0,type=data_function,size_nbinom=25){
  ncell = length(cell_mu0)
  # Cell Counts
  if(type=="poisson"){
    count0 <- sapply(1:ncell, function(x) rpois(n=N0, lambda=cell_mu0[x]))+0.1
    count1 <- sapply(1:ncell, function(x) rpois(n=N1, lambda=cell_mu1[x]))+0.1 
  } else if(type=="neg_bino"){
    count0 <- sapply(1:ncell, function(x) rnbinom(n=N0, mu=cell_mu0[x], size=size_nbinom))+0.1
    count1 <- sapply(1:ncell, function(x) rnbinom(n=N1, mu=cell_mu1[x], size=size_nbinom))+0.1  
  } else if(type=="gaus_cop"){
    rho=0.4
    Sigma <- matrix(rho, ncell, ncell)
    diag(Sigma) <- 1
    # Group 1: Correlated Normal, Transform to Poisson
    Z <- rmvnorm(N0, sigma=Sigma)
    count0 <- matrix(qpois(pnorm(Z), lambda=rep(cell_mu0, each=N0)), nrow=N0, ncol=ncell)+0.1
    # Group 1: Correlated Normal, Transform to Poisson
    Z <- rmvnorm(N1, sigma=Sigma)
    count1 <- matrix(qpois(pnorm(Z), lambda=rep(cell_mu1, each=N1)), nrow=N1, ncol=ncell)+0.1
  }
  # Cell Data Bound
  data.cells <- rbind(count0,count1); colnames(data.cells) <- names(cell_mu0)
  # Cell Proportions
  data.proportions <- data.cells/rowSums(data.cells)
  return(list("Counts"=data.cells,"Proportions"=data.proportions))
}

## ALR Transformation ----
alr <- function(data_alr=data_alr){
  ALR_data = data.frame("id"=NA,"ALR_Numerator"=NA,"ALR_Denominator"=NA,"value"=NA)
  log_data <- log(data_alr)
  for(i in 1:length(cell_mu0)){
    ALR_i = data.frame(apply(log_data, 2, function(x) log_data[,i] - x))%>%
      mutate(id=ids)%>%
      pivot_longer(cols=names(cell_mu0),names_to="ALR_Denominator")%>%
      dplyr::filter(ALR_Denominator != names(cell_mu0)[i])%>%
      mutate(ALR_Numerator=names(cell_mu0)[i])
    ALR_data = bind_rows(ALR_data,ALR_i)
  }
  ALR_data <- ALR_data %>% dplyr::filter(!is.na(id))
  return(ALR_data)
}

## ILR Transformation ----
ilr_contrast <- function(props,numerator_set,denominator_set,scale=TRUE){
  if(scale==TRUE){
    r <- length(numerator_set)
    s <- length(denominator_set)
    gm_A <- apply(props[,numerator_set,drop=FALSE], 1, function(x) exp(mean(log(x))))
    gm_B <- apply(props[,denominator_set,drop=FALSE], 1, function(x) exp(mean(log(x))))
    sqrt((r*s)/(r+s))*log(gm_A/ gm_B)
  } else {
    gm_A <- apply(props[,numerator_set,drop=FALSE], 1, function(x) exp(mean(log(x))))
    gm_B <- apply(props[,denominator_set,drop=FALSE], 1, function(x) exp(mean(log(x))))
    log(gm_A/ gm_B)
  }

}

## SLR Transformation ----
slr_contrast <- function(props,numerator_set,denominator_set){
  gm_A <- apply(props[,numerator_set,drop=FALSE], 1, function(x) sum(x,na.rm=TRUE))
  gm_B <- apply(props[,denominator_set,drop=FALSE], 1, function(x) sum(x,na.rm=TRUE))
  log(gm_A/ gm_B)
}

## Test a given contrast ----
test_contrasts <- function(props, group, contrast_list) {
  bind_rows(mapply(function(nums, denoms) {
    num_ordered <- nums[order(nums)]
    denom_ordered <- denoms[order(denoms)]
    ilr_vals <- ilr_contrast(props, num_ordered, denom_ordered)
    t_test <- t.test(ilr_vals[group==1], ilr_vals[group==0])
    data.frame(
      numerator = paste(num_ordered, collapse="+"),
      denominator = paste(denom_ordered, collapse="+"),
      n_denom = length(denom_ordered),
      mean_diff = t_test$estimate[1] - t_test$estimate[2],
      t_stat = t_test$statistic,
      p_value = t_test$p.value,
      row.names = NULL
    )
  }, contrast_list$numerator_set, contrast_list$denom_set, SIMPLIFY=FALSE))
}

## Simulation Function v1 ----
Small_Sim_v1 <- function(config=config){
  with(config,{
    
    ## Simulation Parameters ----
    cell_mu0 = cell_mus$cell_mu0*scale ;  cell_mu1 = cell_mus$cell_mu1*scale
    power_cell = which(cell_mu1 != cell_mu0)
    fdr_cell = which(cell_mu1 == cell_mu0)
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
      N0=sample.size[k] ; N1=sample.size[k] ; N=N0+N1
      cat(paste0("**Outer Loop** Sample Size: ",N,"\n"))
      cat(paste0("Iteration: "))
      X = c(rep(0,times=N0),rep(1,times=N1))
      ids = paste0("P_",1:N)
      base_data <- data.frame("id"=ids,"X"=X)
      
      ### Inner Loop ----
      for(i in 1:Nsims){
        
        # Data Generation
        toy_data <- data.generation(ncell=ncell,cell_mu0=cell_mu0,cell_mu1=cell_mu1,N1=N1,N0=N0,type=data_function)
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
          result_i <- ilr_DA(props, celltype, X,perm_test=TRUE,perms=500,p_correction="BH")
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
      Model_Results_k = cbind(Model_Results_k,N)
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


############################################################
# Methods ----

## ilr-DA General Form Model v1 ----
### Define a Test Function
# Adjusted test function
test_fn_adjusted <- function(lr_vals, X, covars=NULL){
  # Define Analysis Dataframe
  data <- data.frame(lr=lr_vals, X=X)
  if(!is.null(covars)) data <- cbind(data, covars)
  # Fit Model
  lm <- summary(glm(lr ~ ., data=data, family=gaussian()))
  coef_row <- coef(lm)["X",]
  list(t_stat=coef_row[3], p_value=coef_row[4])
}

### ilr-DA general
ilr_DA_General <- function(props, celltype, X, covars=NULL,perms=500,test_function=NULL,omnibus=FALSE,p_correction="BH", perm_test=TRUE){
  
  ### Define test Function ----
  # Default test function - binary two-group t-test
  default_test_fn <- function(lr_vals, X, data=NULL){
    t_result <- t.test(lr_vals[X==1], lr_vals[X==0])
    list(t_stat=as.numeric(t_result$statistic), 
         p_value=t_result$p.value)
  }
  # Then at the top of ilr_DA_General, set the test function:
  if(is.null(test_function)){
    test_fn <- default_test_fn
  } else {
    test_fn <- test_function
  }
  
  ### Sub-functions -----
  # Prepare Test Contrast Function
  test_contrasts <- function(props, group, contrast_list) {
    bind_rows(mapply(function(nums, denoms) {
      num_ordered <- nums[order(nums)]
      denom_ordered <- denoms[order(denoms)]
      ilr_vals <- ilr_contrast(props, num_ordered, denom_ordered)
      result <- test_fn(lr_vals=ilr_vals, X=X, covars=covars)
      data.frame(
        numerator = paste(num_ordered, collapse="+"),
        denominator = paste(denom_ordered, collapse="+"),
        n_denom = length(denom_ordered),
        t_stat = result$t_stat,
        p_value = result$p_value,
        row.names = NULL
      )
    }, contrast_list$numerator_set, contrast_list$denom_set, SIMPLIFY=FALSE))
  }
  # Prepare SBP_Search_W Function
  SBP_Search_subfunction <- function(props, celltype, X, W_order){
    remaining <- celltype[W_order]  # ordered by W descending
    peeled <- c()
    all_results <- data.frame()
    level=1
    
    repeat{
      # Test all contrasts in remaining set
      nums <- remaining
      denoms <- lapply(remaining, function(x) setdiff(remaining, x))
      contrast_list <- list("numerator_set"=as.list(nums), "denom_set"=denoms)
      results_i <- test_contrasts(props, X, contrast_list)
      results_i$level <- level
      results_i$peeled_so_far <- paste(peeled, collapse="+")
      all_results <- bind_rows(all_results, results_i)
      
      # Check stopping rule
      flag_significance <- all(results_i$p_value >= 0.05)
      
      if(flag_significance | length(remaining) <= 2){
        # Clean terminal found
        return(list(
          "nonDA_cells" = remaining,
          "DA_cells" = peeled,
          "Results_Frame"=all_results
        ))
      }
      # Peel off highest W cell type still in remaining
      next_peel <- remaining[1]  # already W-ordered
      peeled <- c(peeled, next_peel)
      remaining <- remaining[-1]
      level <- level + 1
    }
  }
  # Prepare Permutation Test Function 
  Permutation.test.v2 <- function(props, X, DA_cells, nonDA_cells,all_celltypes, n_perm=2000, max_its=10,test_perm=FALSE){
    DA_candidates <- DA_cells
    nDA_candidates <- nonDA_cells
    Permutation_Results=data.frame()
    ## DA Candidates: Re-testing and Pruning ----
    if(length(nDA_candidates)>=1 & length(DA_candidates)>=1){
      p_value <- c()
      for(i in 1:length(DA_candidates)){
        ilr <- ilr_contrast(props, 
                            numerator_set=DA_candidates[i], 
                            denominator_set=nDA_candidates)
        result <- test_fn(lr_vals=ilr, X=X, covars=covars)
        t1_obs <- result$t_stat
        t1_p <- result$p_value
        # Permutation Testing
        if(test_perm==TRUE){
          t1_perm <- replicate(n_perm, {
            group_perm <- sample(X)
            test_fn(lr_vals=ilr, X=group_perm, covars=covars)$t_stat
          })
          perm_p <- mean(abs(t1_perm) >= abs(t1_obs))
          p_value[i] <- perm_p
        } else{
          p_value[i] <- t1_p
          perm_p <- NA
        }
        
        # Results
        results_i = c("Type"="Retest","Numerator"=DA_candidates[i],"Denominator"=paste0(nDA_candidates,collapse=":"),"T_obs"=t1_obs,"Raw_p"=t1_p,"Perm_p"=perm_p)
        Permutation_Results=bind_rows(Permutation_Results,results_i)
      }
      names(p_value) <- DA_candidates
      
      ## Pruning loop  ----
      its <- 1
      while(any(p_value > 0.05, na.rm=TRUE) & its < max_its & length(DA_candidates) > 0){
        
        # Move failing candidates to nonDA reference
        failing <- DA_candidates[which(p_value > 0.05)]
        nDA_candidates <- c(nDA_candidates, failing)
        DA_candidates <- setdiff(all_celltypes, nDA_candidates)
        
        if(length(DA_candidates) == 0) break
        
        p_value <- c()
        for(i in 1:length(DA_candidates)){
          ilr <- ilr_contrast(props, 
                              numerator_set=DA_candidates[i], 
                              denominator_set=nDA_candidates)
          
          # Observed t/p values
          result <- test_fn(lr_vals=ilr, X=X, covars=covars)
          t1_obs <- result$t_stat
          t1_p <- result$p_value
          # Permutation Testing
          if(test_perm==TRUE){
            t1_perm <- replicate(n_perm, {
              group_perm <- sample(X)
              test_fn(lr_vals=ilr, X=group_perm, covars=covars)$t_stat
            })
            perm_p <- mean(abs(t1_perm) >= abs(t1_obs))
            p_value[i] <- perm_p
          } else{
            p_value[i] <- t1_p
            perm_p <- NA
          }
          
          # Results
          results_i = c("Type"="Pruning","Numerator"=DA_candidates[i],"Denominator"=paste0(nDA_candidates,collapse=":"),"T_obs"=t1_obs,"Raw_p"=t1_p,"Perm_p"=perm_p)
          Permutation_Results=bind_rows(Permutation_Results,results_i)
        }
        names(p_value) <- DA_candidates
        its <- its + 1
      }
    }
    
    if(!exists("p_value")){
      p_value = c()
    }
    ## Return Results  ----
    return(list(
      "nonDA_cells" = nDA_candidates,
      "DA_cells" = DA_candidates,
      "P_values" = p_value,
      "Permutation_Results"=Permutation_Results
    ))
  }
  
  ### ILR-DA ----
  #### ALR Testing ----
  pairs <- which(upper.tri(matrix(0, length(celltype), length(celltype))), arr.ind=TRUE) # Generate unique ALRs only
  pval_matrix <- matrix(1, length(celltype), length(celltype)) # Compute p-values for each unique pair
  diag(pval_matrix) <- NA
  for(r in 1:nrow(pairs)){
    i <- pairs[r, 1]
    j <- pairs[r, 2]
    alr_ij <- log(props[, i]) - log(props[, j])
    p <- test_fn(lr_vals=alr_ij, X=X, covars=covars)$p_value
    pval_matrix[i, j] <- p
    pval_matrix[j, i] <- p  # mirror
  }
  dimnames(pval_matrix) <- list(celltype,celltype)
  W <- rowSums(pval_matrix <= 0.05, na.rm=TRUE) ;  names(W) <- celltype # W statistic per cell type
  Test_Matrix <- pval_matrix <= 0.05 ;  colnames(Test_Matrix) <- celltype ; rownames(Test_Matrix) <- celltype 
  search_order <- order(W, decreasing=TRUE) ; celltype_ordered <- celltype[search_order] # Search order
  
  #### Search Subfunction ----
  SBP_Search = SBP_Search_subfunction(props, celltype, X, W_order=search_order)
  nonDA_cells = SBP_Search$nonDA_cells ; DA_cells = SBP_Search$DA_cells 
  Results_Frame = SBP_Search$Results_Frame %>% dplyr::select(-denominator)
  
  #### Permutation Testing ----
  perm_test_result <- Permutation.test.v2(props, X, DA_cells=DA_cells, nonDA_cells=nonDA_cells,all_celltypes=celltype, n_perm=perms, max_its=3,test_perm=perm_test)
  nonDA_cells = perm_test_result$nonDA_cells ; DA_cells = perm_test_result$DA_cells ; Permutation_Results=perm_test_result$Permutation_Results
  
  
  ### P-value Correction ----
  # Check for whether any DA cells exist
  if(!is.null(DA_cells)){
    final_round <- as.data.frame(Permutation_Results) %>% group_by(Numerator) %>% slice_tail(n=1) %>% ungroup()
    model_pvals <- as.numeric(final_round$Perm_p) ; names(model_pvals) <- final_round$Numerator
    if(!is.null(p_correction)){
      model_pvals_corrected <- p.adjust(model_pvals, method=p_correction)
    }
    # Results
    Result_i <- as.integer(celltype %in% DA_cells) ;  names(Result_i) <- colnames(props)
    DA_corrected <- DA_cells[model_pvals_corrected <= 0.05]
    Result_i_corrected <- as.integer(celltype %in% DA_corrected) ;  names(Result_i_corrected) <- colnames(props)
  }else{
    Result_i <- rep(0,times=ncol(props)) ; names(Result_i) <- colnames(props)
    Result_i_corrected <- rep(0,times=ncol(props)) ; names(Result_i_corrected) <- colnames(props)
  }
  ### Return Results -----
  # Edge case for full Null composition
  if(all( celltype  %in% nonDA_cells )){DA_cells="Null"}
  # Results List
  return(list("nonDA_cells"=nonDA_cells,"DA_cells"=DA_cells,"Results_Frame"=Results_Frame,
              "W"=W,"Test_Matrix"=pval_matrix,"Permutation_Results"=Permutation_Results,
              "Corrected_Results"=Result_i_corrected,"Uncorrected_Results"=Result_i))
}


## ILR-based DA Model v1 ----
ilr_DA <- function(props, celltype, X, perm_test=TRUE,perms=500,omnibus=FALSE,p_correction="BH"){
  
  ### Sub-functions -----
  # Prepare Test Contrast Function
  test_contrasts <- function(props, group, contrast_list) {
    bind_rows(mapply(function(nums, denoms) {
      num_ordered <- nums[order(nums)]
      denom_ordered <- denoms[order(denoms)]
      ilr_vals <- ilr_contrast(props, num_ordered, denom_ordered)
      t_test <- t.test(ilr_vals[group==1], ilr_vals[group==0])
      data.frame(
        numerator = paste(num_ordered, collapse="+"),
        denominator = paste(denom_ordered, collapse="+"),
        n_denom = length(denom_ordered),
        mean_diff = t_test$estimate[1] - t_test$estimate[2],
        t_stat = t_test$statistic,
        p_value = t_test$p.value,
        row.names = NULL
      )
    }, contrast_list$numerator_set, contrast_list$denom_set, SIMPLIFY=FALSE))
  }
  # Prepare SBP_Search_W Function
  SBP_Search_subfunction <- function(props, celltype, X, W_order){
    remaining <- celltype[W_order]  # ordered by W descending
    peeled <- c()
    all_results <- data.frame()
    level=1
    
    repeat{
      # Test all contrasts in remaining set
      nums <- remaining
      denoms <- lapply(remaining, function(x) setdiff(remaining, x))
      contrast_list <- list("numerator_set"=as.list(nums), "denom_set"=denoms)
      results_i <- test_contrasts(props, X, contrast_list)
      results_i$level <- level
      results_i$peeled_so_far <- paste(peeled, collapse="+")
      all_results <- bind_rows(all_results, results_i)
      
      # Check stopping rule
      flag_significance <- all(results_i$p_value >= 0.05)
      
      if(flag_significance | length(remaining) <= 2){
        # Clean terminal found
        return(list(
          "nonDA_cells" = remaining,
          "DA_cells" = peeled,
          "Results_Frame"=all_results
        ))
      }
      # Peel off highest W cell type still in remaining
      next_peel <- remaining[1]  # already W-ordered
      peeled <- c(peeled, next_peel)
      remaining <- remaining[-1]
      level <- level + 1
    }
  }
  # Prepare Permutation Test Function 
  Permutation.test.v2 <- function(props, X, DA_cells, nonDA_cells,all_celltypes, n_perm=2000, max_its=10,test_perm=FALSE){
    DA_candidates <- DA_cells
    nDA_candidates <- nonDA_cells
    Permutation_Results=data.frame()
    ## DA Candidates: Re-testing and Pruning ----
    if(length(nDA_candidates)>=1 & length(DA_candidates)>=1){
      p_value <- c()
      for(i in 1:length(DA_candidates)){
        ilr <- ilr_contrast(props, 
                            numerator_set=DA_candidates[i], 
                            denominator_set=nDA_candidates)
        t1 <- t.test(ilr[X==1], ilr[X==0])
        t1_obs <- t1$statistic
        t1_p <- t1$p.value
        
        # Permutation Testing
        if(test_perm==TRUE){
          t1_perm <- replicate(n_perm, {
            group_perm <- sample(X)
            t.test(ilr[group_perm==1], ilr[group_perm==0])$statistic
          })
          perm_p <- mean(abs(t1_perm) >= abs(t1_obs))
          p_value[i] <- perm_p
        } else{
          p_value[i] <- t1_p
          perm_p <- NA
        }
        
        # Results
        results_i = c("Type"="Retest","Numerator"=DA_candidates[i],"Denominator"=paste0(nDA_candidates,collapse=":"),"T_obs"=t1_obs,"Raw_p"=t1_p,"Perm_p"=perm_p)
        Permutation_Results=bind_rows(Permutation_Results,results_i)
      }
      names(p_value) <- DA_candidates
      
      ## Pruning loop  ----
      its <- 1
      while(any(p_value > 0.05, na.rm=TRUE) & its < max_its & length(DA_candidates) > 0){
        
        # Move failing candidates to nonDA reference
        failing <- DA_candidates[which(p_value > 0.05)]
        nDA_candidates <- c(nDA_candidates, failing)
        DA_candidates <- setdiff(all_celltypes, nDA_candidates)
        
        if(length(DA_candidates) == 0) break
        
        p_value <- c()
        for(i in 1:length(DA_candidates)){
          ilr <- ilr_contrast(props, 
                              numerator_set=DA_candidates[i], 
                              denominator_set=nDA_candidates)
          
          # Observed t/p values
          t1 <- t.test(ilr[X==1], ilr[X==0])
          t1_obs <- t1$statistic
          t1_p <- t1$p.value
          
          # Permutation Testing
          if(test_perm==TRUE){
            t1_perm <- replicate(n_perm, {
              group_perm <- sample(X)
              t.test(ilr[group_perm==1], ilr[group_perm==0])$statistic
            })
            perm_p <- mean(abs(t1_perm) >= abs(t1_obs))
            p_value[i] <- perm_p
          } else{
            p_value[i] <- t1_p
            perm_p <- NA
          }
          
          # Results
          results_i = c("Type"="Pruning","Numerator"=DA_candidates[i],"Denominator"=paste0(nDA_candidates,collapse=":"),"T_obs"=t1_obs,"Raw_p"=t1_p,"Perm_p"=perm_p)
          Permutation_Results=bind_rows(Permutation_Results,results_i)
        }
        names(p_value) <- DA_candidates
        its <- its + 1
      }
    }
    
    if(!exists("p_value")){
      p_value = c()
    }
    ## Return Results  ----
    return(list(
      "nonDA_cells" = nDA_candidates,
      "DA_cells" = DA_candidates,
      "P_values" = p_value,
      "Permutation_Results"=Permutation_Results
    ))
  }

  ### ILR-DA ----
    #### ALR Testing ----
    pairs <- which(upper.tri(matrix(0, length(celltype), length(celltype))), arr.ind=TRUE) # Generate unique ALRs only
    pval_matrix <- matrix(1, length(celltype), length(celltype)) # Compute p-values for each unique pair
    diag(pval_matrix) <- NA
    for(r in 1:nrow(pairs)){
      i <- pairs[r, 1]
      j <- pairs[r, 2]
      alr_ij <- log(props[, i]) - log(props[, j])
      p <- t.test(alr_ij[X==1], alr_ij[X==0])$p.value
      pval_matrix[i, j] <- p
      pval_matrix[j, i] <- p  # mirror
    }
    dimnames(pval_matrix) <- list(celltype,celltype)
    W <- rowSums(pval_matrix <= 0.05, na.rm=TRUE) ;  names(W) <- celltype # W statistic per cell type
    Test_Matrix <- pval_matrix <= 0.05 ;  colnames(Test_Matrix) <- celltype ; rownames(Test_Matrix) <- celltype 
    search_order <- order(W, decreasing=TRUE) ; celltype_ordered <- celltype[search_order] # Search order
    
    #### Search Subfunction ----
    SBP_Search = SBP_Search_subfunction(props, celltype, X, W_order=search_order)
    nonDA_cells = SBP_Search$nonDA_cells ; DA_cells = SBP_Search$DA_cells 
    Results_Frame = SBP_Search$Results_Frame %>% dplyr::select(-denominator)
    
    #### Permutation Testing ----
    perm_test_result <- Permutation.test.v2(props, X, DA_cells=DA_cells, nonDA_cells=nonDA_cells,all_celltypes=celltype, n_perm=perms, max_its=3,test_perm=perm_test)
    nonDA_cells = perm_test_result$nonDA_cells ; DA_cells = perm_test_result$DA_cells ; Permutation_Results=perm_test_result$Permutation_Results
    
    
    ### P-value Correction ----
    # Check for whether any DA cells exist
    if(!is.null(DA_cells)){
      final_round <- as.data.frame(Permutation_Results) %>% group_by(Numerator) %>% slice_tail(n=1) %>% ungroup()
      model_pvals <- as.numeric(final_round$Perm_p) ; names(model_pvals) <- final_round$Numerator
      if(!is.null(p_correction)){
        model_pvals_corrected <- p.adjust(model_pvals, method=p_correction)
      }
      # Results
      Result_i <- as.integer(celltype %in% DA_cells) ;  names(Result_i) <- colnames(props)
      DA_corrected <- DA_cells[model_pvals_corrected <= 0.05]
      Result_i_corrected <- as.integer(celltype %in% DA_corrected) ;  names(Result_i_corrected) <- colnames(props)
    }else{
      Result_i <- rep(0,times=ncol(props)) ; names(Result_i) <- colnames(props)
      Result_i_corrected <- rep(0,times=ncol(props)) ; names(Result_i_corrected) <- colnames(props)
    }
    ### Return Results -----
    # Edge case for full Null composition
    if(all( celltype  %in% nonDA_cells )){DA_cells="Null"}
    # Results List
    return(list("nonDA_cells"=nonDA_cells,"DA_cells"=DA_cells,"Results_Frame"=Results_Frame,
                "W"=W,"Test_Matrix"=pval_matrix,"Permutation_Results"=Permutation_Results,
                "Corrected_Results"=Result_i_corrected,"Uncorrected_Results"=Result_i))
}



## Poisson Check ----
Poisson_Check <- function(counts,X,p_correction="BH"){
  counts_2 <- floor(counts)
  model_pvals <- NULL
  ncell <- ncol(counts)
  for(j in 1:ncell){
    lm <- summary(glm(formula = counts_2[,j] ~ as.factor(X), data = data.frame(X = X), family = poisson(link = "log")))
    model_pvals[j] <- coef(lm)[2,4]
  }
  # P-value Correction 
  if(!is.null(p_correction)){
    model_pvals_corrected <- p.adjust(model_pvals, method=p_correction)
  }
  # Results
  Result_i = as.integer(model_pvals<=0.05) ;  names(Result_i) <- colnames(counts)
  Result_i_corrected = as.integer(model_pvals_corrected<=0.05) ;  names(Result_i_corrected) <- colnames(counts)
  return(list("Corrected_Results"=Result_i_corrected,"Uncorrected_Results"=Result_i))
}

## Negative Binomial Check ----
NB_Check <- function(X, counts, celltype,p_correction="BH"){
  counts_2 <- floor(counts)
  model_pvals <- tryCatch({
    # Fit Model
    suppressWarnings({
      model_pvals <- sapply(1:length(celltype), function(k){
        lm <- summary(glm.nb(counts_2[,k] ~ X))
        coef(lm)[2,4]
      })
    })
    names(model_pvals) <- colnames(counts)
    model_pvals
  }, error = function(e) rep(1, length(celltype)))
  # P-value Correction 
  if(!is.null(p_correction)){
    model_pvals_corrected <- p.adjust(model_pvals, method=p_correction)
  }
  # Results
  Result_i = as.integer(model_pvals<=0.05) ;  names(Result_i) <- colnames(counts)
  Result_i_corrected = as.integer(model_pvals_corrected<=0.05) ;  names(Result_i_corrected) <- colnames(counts)
  return(list("Corrected_Results"=Result_i_corrected,"Uncorrected_Results"=Result_i))
}

## Dirichlet Regression - Common Model - with trycatch ----
## Dirichlet DA with trycatch
Dirichlet_DA <- function(X, props, celltype,p_correction="BH"){
  # Prepare Data
  ncell <- length(celltype)
  x_dir <- as.factor(X)
  data <- data.frame(x_dir=x_dir,props)
  data$y <- DR_data(data[,2:(ncell+1)])
  # Fit Models
  result <- tryCatch({
    lm <- DirichReg(y ~ x_dir, data)
    lm <- summary(lm)
    # extract p-values
    coef_matrix <- lm$coef.mat
    model_pvals <- coef_matrix[seq(2,nrow(coef_matrix),2),4]
    list("result"=model_pvals,"error"=0)
  }, error = function(e){
    # Return NAs on convergence failure
    result =rep(NA, ncell) ; names(result) <- celltype
    list("result"=result,"error"=1)
  })
  # P-value Correction
  model_pvals <- result$result
  if(!is.null(p_correction)){
    model_pvals_corrected <- p.adjust(model_pvals, method=p_correction)
  }
  # Results
  Result_i = as.integer(model_pvals<=0.05) ;  names(Result_i) <- colnames(props)
  Result_i_corrected = as.integer(model_pvals_corrected<=0.05) ;  names(Result_i_corrected) <- colnames(props)
  return(list("Corrected_Results"=Result_i_corrected,"Uncorrected_Results"=Result_i,
              "pvals"=model_pvals,"pvals_corrected"=model_pvals_corrected,
              "error"=result$error))
}
## SSDA ----
# Fisher Statistic Calculation 
fisher_statistic_calc <- function(p_values){
  if(length(na.omit(p_values))==0){
    fisher_statistic=NA
  } else{
    p_values=p_values[!is.na(p_values)]
    fisher_statistic=-2*sum(log(p_values)) 
  }
  return(fisher_statistic)
}
# Fisher Test Full
Fisher_Test_Full <- function(ncell=ncell,data.proportions=data.proportions,group=group,N=N,cell.type=cell.type,power_cell=power_cell,fdr_cell=fdr_cell,FDR_control=FDR_control,FDR_control_type=FDR_control_type){
  # Output Objects
  model_test <- NULL
  model_pvals <- NULL
  model_pvals_corrected <- NULL
  # SUB HYPOTHESIS P-VALUES
  alr_pvals<- data.frame(matrix(NA,nrow=ncol(data.proportions),ncol=ncol(data.proportions)))
  dimnames(alr_pvals) <- list(cell.type,cell.type)
  for(i in 1:ncell){
    y1 <- log(data.proportions)
    y = apply(y1, 2, function(x) x - y1[,i])  
    x <- as.factor(group) 
    alr_lm <- summary(lm(formula = y ~ x, data = data.frame(X = x)))
    alr_pvals[i,] <- c(unlist(lapply(coef(alr_lm), function(x)x[2,4])))
  }
  diag(alr_pvals) <- NA
  alr_pvals_test <- alr_pvals # Save ALR P-value Object for test
  if(FDR_control_type == "Sub Hypothesis" & FDR_control!="None"){
    alr_pvals_test <- apply(alr_pvals, 1,function(x) p.adjust(x,method=FDR_control))
  }
  ## Subcomposition Test ----
  sequence_matrix <- matrix(NA,nrow=ncell,ncol=ncell)
  dimnames(sequence_matrix) <- list(paste("Step",seq(1,ncell,1)),cell.type)
  test_pvals <- c(0,0,0);i <- 1
  # Loop
  while(min(test_pvals,na.rm=TRUE)<=0.05){
    fisher_statistic_row <- apply(alr_pvals_test, 1,function(x) fisher_statistic_calc(x))
    df <- 2*length(na.omit(fisher_statistic_row))
    test_pvals <- pchisq(fisher_statistic_row,df,lower.tail=FALSE)
    if(FDR_control_type == "Fisher Test" & FDR_control!="None"){
      test_pvals <- p.adjust(test_pvals,method=FDR_control)
    }
    sequence_matrix[i,] <- test_pvals
    index_omit <- which(fisher_statistic_row==max(fisher_statistic_row,na.rm=TRUE))
    alr_pvals_test[index_omit,] <- NA;alr_pvals_test[,index_omit] <- NA
    i <- i+1
  } 
  #sequence_matrix <- as.data.frame(sequence_matrix[rowSums(sequence_matrix,na.rm=TRUE)!=0,])
  model_test <- apply(sequence_matrix,2,function(x) max(x,na.rm=TRUE))
  model_test <- ifelse(model_test<= 0.05,1,0)
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="Fisher","Power"=model_power,"FDR"=model_fdr)
  return(list("Result"=result,"Model.Test"=model_test)) 
}


## ANCOM-BC2 -----
ANCOMBC2_test <- function(props,X){
  # Phyloseq object
  X_factor <- ifelse(X==0,"Case","Control")
  X_factor <- factor(X_factor,levels=c("Control","Case"))
  data.proportions.ancombc2 <- props
  rownames(data.proportions.ancombc2) <- paste0("patient_",1:nrow(data.proportions.ancombc2))
  sample_data <- data.frame(sample=rownames(data.proportions.ancombc2),group=X_factor)
  rownames(sample_data) <- rownames(data.proportions.ancombc2)
  otu_table <- data.proportions.ancombc2
  ps <- phyloseq(otu_table(otu_table,taxa_are_rows=FALSE),sample_data(sample_data))
  # ANCOM BC2
  invisible(capture.output(
    suppressMessages(
      suppressWarnings(
        result <- tryCatch({
          ancombc2(data = ps, assay_name = "counts",
                   fix_formula = "group",
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   prv_cut = 0.10, s0_perc = 0.05,
                   group = "group",
                   alpha = 0.05, verbose = FALSE,
                   iter_control = list(tol = 1e-5, max_iter = 10,verbose = FALSE),
                   em_control = list(tol = 1e-5, max_iter = 10),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 50))
          
          
        }, error = function(e) rep(NA, length(celltype))
        )
      )
    )
  ))
  # Results
  res_prim = result$res
  Results_ANCOM <- as.integer(res_prim$diff_groupCase)
  names(Results_ANCOM) <- colnames(props)
  return(Results_ANCOM)
}

## Logratio Lasso ----
# Log ratio lasso - Approx Forward Stepwise
logratio_SW_DA <- function(props, X, celltype){
  w <- log(props)
  y <- as.numeric(X)
  N <- nrow(props)
  K <- length(celltype)
  
  # Adaptive k_max and n_folds based on composition size and sample size
  k_max <- min(3, K-1, floor(N/10))  # conservative k_max
  n_folds <- min(5, floor(N/5))       # ensure reasonable fold size
  
  # Guard against degenerate cases
  if(k_max < 1 | n_folds < 2) return(character(0))
  
  suppressWarnings({
    afs_model <- approximate_fs(w, y, k_max=k_max)
    afs_cv <- cv_approximate_fs(w, y, k_max=k_max, n_folds=n_folds)
  })
  
  k_best <- afs_cv$k_min
  beta <- afs_model$beta[, 1]
  
  return(beta)
}

# Log-ratio Lasso: Lasso CV
logratio_Lasso_DA <- function(props, X, celltype){
  w <- log(props)
  y <- as.numeric(X)
  centered_w <- scale(w, center=TRUE, scale=FALSE)
  centered_y <- X - mean(X)
  N <- nrow(props)
  K <- length(celltype)
  
  # Adaptive k_max and n_folds based on composition size and sample size
  k_max <- min(3, K-1, floor(N/10))  # conservative k_max
  n_folds <- min(5, floor(N/5))       # ensure reasonable fold size
  
  # Guard against degenerate cases
  if(k_max < 1 | n_folds < 2) return(character(0))
  
  # Two-stage procedure
  invisible(capture.output(
    ts_model <- cv_two_stage(centered_w, centered_y, k_max=k_max, n_folds=n_folds)
  ))
  beta <- ts_model$beta_min

  return(beta)
}



## Dirichlet Regression - Alt Model - with trycatch ----
Dirichlet_DA_Alt <- function(X, props, celltype,model_type="Alt"){
    ## Dataset Prep
    ncell <- length(celltype)
    x_dir <- as.factor(X)
    data <- data.frame(x_dir=x_dir,props)
    data$y <- DR_data(data[,2:(ncell+1)])
    tol1_model = 1e-5
    tol2_model = 1e-10
    iterlim_model = 1000
    
    ## ALR testing for reference cell
    if(model_type=="Alt"){
      pairs <- which(upper.tri(matrix(0, length(celltype), length(celltype))), arr.ind=TRUE) # Generate unique ALRs only
      pval_matrix <- matrix(1, length(celltype), length(celltype)) # Compute p-values for each unique pair
      diag(pval_matrix) <- NA
      for(r in 1:nrow(pairs)){
        i <- pairs[r, 1]
        j <- pairs[r, 2]
        alr_ij <- log(props[, i]) - log(props[, j])
        p <- t.test(alr_ij[X==1], alr_ij[X==0])$p.value
        pval_matrix[i, j] <- p
        pval_matrix[j, i] <- p  # mirror
      }
      W <- rowSums(pval_matrix <= 0.05, na.rm=TRUE) ;  names(W) <- celltype # W statistic per cell type
      base_cell_alt = order(W, decreasing=TRUE)[ncell] ; base_cell_name = celltype[base_cell_alt]
    }
    ## Model Fitting
    result <- tryCatch({
      if(model_type=="Alt"){
        # Alternate Model
        lm <- DirichReg(y ~ x_dir, data, base=base_cell_alt,
                        model="alternative",
                        control=list(iterlim=iterlim_model,tol1=tol1_model,tol2=tol2_model))
        lm <- summary(lm)
        # extract p-values
        coef_matrix <- lm$coef.mat
        model_pvals <- coef_matrix[seq(2,nrow(coef_matrix),2),4] ; names(model_pvals) <- setdiff(celltype,base_cell_name)
        model_pvals <- c(model_pvals,setNames(NA, base_cell_name))
        model_pvals <- model_pvals[celltype]
        result = as.integer(model_pvals <= 0.05) ; names(result) <- celltype
      } else if(model_type=="Common"){
        # Common Model
        lm <- DirichReg(y ~ x_dir, data,
                        control=list(iterlim=iterlim_model,tol1=tol1_model,tol2=tol2_model))
        lm <- summary(lm)
        # extract p-values
        coef_matrix <- lm$coef.mat
        model_pvals <- coef_matrix[seq(2,nrow(coef_matrix),2),4]
        result = as.integer(model_pvals <= 0.05)
        names(result) <- celltype
      }
      result = list("result"=result,"error"=0)
    },
    error = function(e){
      # Return NAs on convergence failure
      result =rep(NA, ncell) ; names(result) <- celltype
      result = list("result"=result,"error"=1)
    }
    # , warning = function(w){
    #   # Optionally catch warnings too
    #   rep(0, ncell)
    # }
    )
    return(result)
  }

## (Defunct)  ILR-based DA Model v2 - GLM Basis ----
ilr_DA_GLM <- function(props, celltype, X,covars=NULL, general_formula = "~X",
                       search_type="HMP", perm_test=TRUE,perms=500,omnibus=FALSE,max_its=10){
  ## Parameters -----
  alr_formula <- as.formula(paste0("alr_ij",general_formula))
  contrasts_formula <- as.formula(paste0("ilr_vals",general_formula))
  celltype <- colnames(props)
  
  ## Prepare Test Contrast Function -----
  test_contrasts <- function(props,X,covars, contrast_list,contrasts_formula) {
    test_contrast_df <- data.frame()
    contrasts_formula <- as.formula(contrasts_formula)
    for(i in 1:length(contrast_list$numerator_set)){
      # Ordering Numerators and Denominators
      nums_ordered <- contrast_list$numerator_set[[i]][order(contrast_list$numerator_set[[i]])]
      denom_ordered <- contrast_list$denom_set[[i]][order(contrast_list$denom_set[[i]])]
      # Calculating ILR Value
      ilr_vals <- ilr_contrast(props, nums_ordered, denom_ordered)
      #
      if(!is.null(covars)){
        data_test = data.frame(ilr_vals,X,covars)
      } else{
        data_test = data.frame(ilr_vals,X)
      }
      # Testing Via GLM
      lm <- summary(glm(contrasts_formula,data=data_test))$coefficients
      t <- lm[rownames(lm)=="X",3] ; p = lm[rownames(lm)=="X",4]
      # Results
      res_contrast <- c(numerator = paste(nums_ordered, collapse="+"),
                        denominator = paste(denom_ordered, collapse="+"),
                        n_denom = length(denom_ordered),
                        mean_diff = 0,
                        t_stat = t,
                        p_value = p)
      test_contrast_df <- bind_rows(test_contrast_df,res_contrast)
    }
    rownames(test_contrast_df) <- NULL
    return(test_contrast_df)
  }
  
  ## Omnibus Test ----
  if(omnibus==TRUE){
    # Omnibus Test
    ilr_coords <- compositions::ilr(props)
    H2.test <- hotelling.test(ilr_coords[X==1,], ilr_coords[X==0,])$pval
  } else{
    H2.test = 0.0001
  }
  
  ## SBP Search W and Testing ----
  if(H2.test > 0.05){
    #cat("Failed to Reject Omnibus Test. No DA Cell types detected.")
    return(list("nonDA_cells" = celltype,"DA_cells" = "Null", "Results_Frame"= "Null"))
  } else{
    ### Search Order ----
    # ALR Testing
    pairs <- which(upper.tri(matrix(0, length(celltype), length(celltype))), arr.ind=TRUE) # Generate unique ALRs only
    pval_matrix <- matrix(1, length(celltype), length(celltype)) # Compute p-values for each unique pair
    diag(pval_matrix) <- NA
    for(r in 1:nrow(pairs)){
      i <- pairs[r, 1]
      j <- pairs[r, 2]
      alr_ij <- log(props[, i]) - log(props[, j])
      #
      if(!is.null(covars)){
        data_test = data.frame(alr_ij,X,covars)
      } else{
        data_test = data.frame(alr_ij,X)
      }
      #
      lm <- summary(glm(alr_formula,data=data_test))$coefficients
      p <- lm[rownames(lm)=="X",4]
      pval_matrix[i, j] <- p
      pval_matrix[j, i] <- p  # mirror
    }
    dimnames(pval_matrix) <- list(celltype,celltype)
    
    # Search Order
    if(search_type=="W"){
      W <- rowSums(pval_matrix <= 0.05, na.rm=TRUE) ;  names(W) <- celltype # W statistic per cell type
      Test_Matrix <- pval_matrix <= 0.05 ;  colnames(Test_Matrix) <- celltype ; rownames(Test_Matrix) <- celltype 
      search_order <- order(W, decreasing=TRUE) ; celltype_ordered <- celltype[search_order] # Search order
    }else if(search_type=="HMP"){
      # Harmonic mean p-values
      W <- apply(pval_matrix, 1, function(row){
        p.hmp(row[!is.na(row)], L=sum(!is.na(row)))
      })
      names(W) <- celltype
      search_order <- order(W, decreasing=FALSE) ; celltype_ordered <- celltype[search_order] # Search order 
    }
    
    ### SBP Search -----
    W_order=search_order
    remaining <- celltype[W_order]  # ordered by W descending
    peeled <- c()
    all_results <- data.frame()
    level=1
    flag_significance=FALSE
    
    while(flag_significance==FALSE){
      # Test all contrasts in remaining set
      nums <- remaining
      denoms <- lapply(remaining, function(x) setdiff(remaining, x))
      contrast_list <- list("numerator_set"=as.list(nums), "denom_set"=denoms)
      results_i <- test_contrasts(props,X,covars, contrast_list=contrast_list, contrasts_formula=contrasts_formula)
      results_i$level <- level
      results_i$peeled_so_far <- paste(peeled, collapse="+")
      all_results <- bind_rows(all_results, results_i)
      # Check stopping rule
      flag_significance <- all(results_i$p_value >= 0.05)
      # Results 
      nonDA_cells=remaining
      DA_cells=peeled
      Results_Frame=all_results
      # Peel off highest W cell type still in remaining
      next_peel <- remaining[1]  # already W-ordered
      peeled <- c(peeled, next_peel)
      remaining <- remaining[-1]
      level <- level + 1
    }
    Results_Frame = Results_Frame %>% dplyr::select(-denominator)
    
    ### Permutation Testing ----
    DA_candidates <- DA_cells
    nDA_candidates <- nonDA_cells
    Permutation_Results=data.frame()
    
    #### DA Candidates: Re-testing and Pruning ----
    if(length(nDA_candidates)>=1 & length(DA_candidates)>=1){
      p_value <- c()
      for(i in 1:length(DA_candidates)){
        ilr_vals <- ilr_contrast(props, 
                                 numerator_set=DA_candidates[i], 
                                 denominator_set=nDA_candidates)
        # Testing Via GLM
        lm <- summary(glm(contrasts_formula,data = data.frame(ilr_vals = ilr_vals, X = X)))$coefficients
        t1_obs <- lm[rownames(lm)=="X",3] ; t1_p = lm[rownames(lm)=="X",4]
        # Permutation Testing
        if(perm_test==TRUE){
          t1_perm <- replicate(perms, {
            # Testing Via GLM
            group_perm <- sample(X)
            lm <- summary(glm(contrasts_formula,data = data.frame(ilr_vals = ilr_vals, X = group_perm)))$coefficients
            lm[rownames(lm)=="X",3]
          })
          perm_p <- mean(abs(t1_perm) >= abs(t1_obs))
          p_value[i] <- perm_p
        } else{
          p_value[i] <- t1_p
          perm_p <- NA
        }
        
        # Results
        results_i = c("Type"="Retest","Numerator"=DA_candidates[i],"Denominator"=paste0(nDA_candidates,collapse=":"),"T_obs"=t1_obs,"Raw_p"=t1_p,"Perm_p"=perm_p)
        Permutation_Results=bind_rows(Permutation_Results,results_i)
      }
      names(p_value) <- DA_candidates
      
      ## Pruning loop  ----
      its <- 1
      while(any(p_value > 0.05, na.rm=TRUE) & its < max_its & length(DA_candidates) > 0){
        
        # Move failing candidates to nonDA reference
        failing <- DA_candidates[which(p_value > 0.05)]
        nDA_candidates <- c(nDA_candidates, failing)
        DA_candidates <- setdiff(all_celltypes, nDA_candidates)
        
        if(length(DA_candidates) == 0) break
        
        p_value <- c()
        for(i in 1:length(DA_candidates)){
          ilr_vals <- ilr_contrast(props, 
                                   numerator_set=DA_candidates[i], 
                                   denominator_set=nDA_candidates)
          # Testing Via GLM
          lm <- summary(glm(contrasts_formula,data = data.frame(ilr_vals = ilr_vals, X = X)))$coefficients
          t1_obs <- lm[rownames(lm)=="X",3] ; t1_p = lm[rownames(lm)=="X",4]
          
          # Permutation Testing
          if(test_perm==TRUE){
            t1_perm <- replicate(n_perm, {
              # Testing Via GLM
              group_perm <- sample(X)
              lm <- summary(glm(contrasts_formula,data = data.frame(ilr_vals = ilr_vals, X = group_perm)))$coefficients
              lm[rownames(lm)=="X",3]
            })
            perm_p <- mean(abs(t1_perm) >= abs(t1_obs))
            p_value[i] <- perm_p
          } else{
            p_value[i] <- t1_p
            perm_p <- NA
          }
          
          # Results
          results_i = c("Type"="Pruning","Numerator"=DA_candidates[i],"Denominator"=paste0(nDA_candidates,collapse=":"),"T_obs"=t1_obs,"Raw_p"=t1_p,"Perm_p"=perm_p)
          Permutation_Results=bind_rows(Permutation_Results,results_i)
        }
        names(p_value) <- DA_candidates
        its <- its + 1
      }
      #
      if(!exists("p_value")){
        p_value = c()
      }
      P_values=p_value
      nonDA_cells = nDA_candidates
      DA_cells = DA_candidates
      
    }
    ## Returning Results ----
    # Edge case for full Null composition
    if(all( celltype  %in% nonDA_cells )){DA_cells="Null"}
    # Return Final Results
    return(list("nonDA_cells"=nonDA_cells,"DA_cells"=DA_cells,"Results_Frame"=Results_Frame,"W"=W,"Test_Matrix"=pval_matrix,"Permutation_Results"=Permutation_Results))
    
  }
}


## (Defunct) ILR-DA Bootstrapping -----
ilr_DA_Boot <- function(props, celltype, X, covars=NULL, general_formula = "~X",
                        search_type="HMP", alpha=0.05 , B_num=500 , Pruning=TRUE, max_its=10){
  
  ###  Creating Parameters ----
  lr_formula <- as.formula(paste0("lr",general_formula))
  celltype <- colnames(props)
  # Analysis Data
  if(!is.null(covars)){
    data_test = data.frame(X,covars)
  } else{
    data_test = data.frame(X)
  }
  
  ### Boostrap ILR Function ----
  bootstrap_ilr <- function(formula=lr_formula, data, B=500, alpha=0.05){
    alpha_boot = alpha
    formula=as.formula(formula)
    B_boot=B
    # Observed ILR T-value
    lm <- summary(glm(formula,data=data),family=gaussian())$coefficients
    t_obs <- lm[rownames(lm)=="X",3]
    # Bootstrap distribution
    t_boot <- replicate(B_boot, {
      # Resample WITH REPLACEMENT, preserving group structure
      idx0 <- sample(which(data$X==0), replace=TRUE)
      idx1 <- sample(which(data$X==1), replace=TRUE)
      idx <- c(idx0, idx1)
      data_boot <- data[idx, ]
      lm <- summary(glm(formula,data=data_boot))$coefficients
      lm[rownames(lm)=="X",3]
    })
    # Bootstrap confidence interval
    list("CI"=quantile(as.numeric(t_boot), c(alpha_boot/2, 1-alpha_boot/2)),"t_obs"=t_obs)
  }
  
  ### Search Order ----
  # ALR Testing
  pairs <- which(upper.tri(matrix(0, length(celltype), length(celltype))), arr.ind=TRUE) # Generate unique ALRs only
  pval_matrix <- matrix(1, length(celltype), length(celltype)) # Compute p-values for each unique pair
  diag(pval_matrix) <- NA
  for(r in 1:nrow(pairs)){
    i <- pairs[r, 1]
    j <- pairs[r, 2]
    alr_ij <- log(props[, i]) - log(props[, j])
    # Data
    data_test$lr <- alr_ij
    #
    lm <- summary(glm(lr_formula,data=data_test,family=gaussian()))$coefficients
    p <- lm[rownames(lm)=="X",4]
    pval_matrix[i, j] <- p
    pval_matrix[j, i] <- p  # mirror
  }
  dimnames(pval_matrix) <- list(celltype,celltype)
  
  # Search Order
  if(search_type=="W"){
    W <- rowSums(pval_matrix <= 0.05, na.rm=TRUE) ;  names(W) <- celltype # W statistic per cell type
    Test_Matrix <- pval_matrix <= 0.05 ;  colnames(Test_Matrix) <- celltype ; rownames(Test_Matrix) <- celltype 
    search_order <- order(W, decreasing=TRUE) ; celltype_ordered <- celltype[search_order] # Search order
  }else if(search_type=="HMP"){
    # Harmonic mean p-values
    W <- apply(pval_matrix, 1, function(row){
      p.hmp(row[!is.na(row)], L=sum(!is.na(row)))
    })
    names(W) <- celltype
    search_order <- order(W, decreasing=FALSE) ; celltype_ordered <- celltype[search_order] # Search order 
  }
  
  ### Build Helmert Contrast -----
  K <- length(celltype_ordered)
  contrast_list <- list(numerator_set=list(), denom_set=list())
  for(i in 1:(K-1)){
    contrast_list$numerator_set[[i]] <- celltype_ordered[i]
    contrast_list$denom_set[[i]] <- celltype_ordered[(i+1):K]
  }
  
  ### Coordinate Search -----
  W_order=search_order
  remaining <- celltype[W_order]  # ordered by W descending
  peeled <- c()
  all_results <- data.frame()
  level=1
  flag_significance=FALSE
  
  # Search!
  while(flag_significance==FALSE & length(remaining) >= 2){
    # Calculating ILR Value
    nums_ordered <- contrast_list$numerator_set[[level]][order(contrast_list$numerator_set[[level]])]
    denom_ordered <- contrast_list$denom_set[[level]][order(contrast_list$denom_set[[level]])]
    ilr_vals <- ilr_contrast(props, nums_ordered, denom_ordered)
    # Data & Bootstrap
    data_test$lr <- ilr_vals
    results_i <- bootstrap_ilr(formula=lr_formula,data=data_test, B=B_num)
    # Check stopping rule
    flag_significance <- !(results_i$CI[1] > 0 | results_i$CI[2] < 0)
    # Additional information
    #results_i <- unlist(results_i)
    results_i <- data.frame(numerator = paste(nums_ordered, collapse=":"),denominator = paste(denom_ordered, collapse=":"),n_denom = length(denom_ordered),
                            t_observed = results_i$t_obs,t_lower = results_i$CI[1],t_upper = results_i$CI[2],
                            alpha=alpha,level=level,peeled_so_far=paste(peeled, collapse=":"))
    rownames(results_i) <- NULL
    all_results <- bind_rows(all_results, results_i)
    # Results 
    nonDA_cells=remaining
    DA_cells=peeled
    # Peel off highest W cell type still in remaining
    next_peel <- remaining[1]  # already W-ordered
    peeled <- c(peeled, next_peel)
    remaining <- remaining[-1]
    level <- level + 1
  }
  Results_Frame = all_results %>% dplyr::select(-denominator)
  
  ### Confirmatory Testing ----
  DA_candidates <- DA_cells
  nDA_candidates <- nonDA_cells
  Confirmatory_Results=data.frame()
  
  #### DA Candidates: Re-testing and Pruning ----
  if(length(nDA_candidates)>=1 & length(DA_candidates)>=1){
    
    ##### Confirmatory Testing of DA Candidates -----
    p_value <- c()
    for(i in 1:length(DA_candidates)){
      # Calculating ILR Value
      ilr_vals <- ilr_contrast(props, 
                               numerator_set=DA_candidates[i], 
                               denominator_set=nDA_candidates)
      # Data & Bootstrap
      data_test$lr <- ilr_vals
      results_i <- bootstrap_ilr(formula=lr_formula,data=data_test, B=B_num)
      # Check stopping rule
      p_value[i] <- as.integer(results_i$CI[1] > 0 | results_i$CI[2] < 0)
      # Results
      #results_i <- unlist(results_i)
      results_i <- data.frame("Type"="Retest","Numerator"=DA_candidates[i],"Denominator"=paste0(nDA_candidates,collapse=":"),
                              t_observed = results_i$t_obs,t_lower = results_i$CI[1],t_upper = results_i$CI[2],
                              alpha=alpha)
      rownames(results_i) <- NULL
      Confirmatory_Results=bind_rows(Confirmatory_Results,results_i)
    }
    names(p_value) <- DA_candidates
    
    
    ##### Pruning loop  ----
    if(Pruning==TRUE){
      its <- 1
      
      while(any(p_value==0, na.rm=TRUE) & its < max_its & length(DA_candidates) > 0){
        
        # Move failing candidates to nonDA reference
        failing <- DA_candidates[which(p_value == 0)]
        nDA_candidates <- c(nDA_candidates, failing)
        DA_candidates <- setdiff(celltype, nDA_candidates)
        
        if(length(DA_candidates) == 0) break
        
        # LOOP
        p_value <- c()
        for(i in 1:length(DA_candidates)){
          # Calculating ILR Value
          ilr_vals <- ilr_contrast(props,numerator_set=DA_candidates[i],denominator_set=nDA_candidates)
          # Data & Bootstrap
          data_test$lr <- ilr_vals
          results_i <- bootstrap_ilr(formula=lr_formula,data=data_test, B=B_num)
          # Check stopping rule
          p_value[i] <- as.integer(results_i$CI[1] > 0 | results_i$CI[2] < 0)
          # Results
          #results_i <- unlist(results_i)
          results_i <- data.frame("Type"="Pruning","Numerator"=DA_candidates[i],"Denominator"=paste0(nDA_candidates,collapse=":"),
                                  t_observed = results_i$t_obs,t_lower = results_i$CI[1],t_upper = results_i$CI[2],
                                  alpha=alpha)
          rownames(results_i) <- NULL
          Confirmatory_Results=bind_rows(Confirmatory_Results,results_i)
        }
        names(p_value) <- DA_candidates
        its <- its + 1
      }
    }
    
    ##### Returning Results  ---- 
    if(!exists("p_value")){
      p_value = c()
    }
    P_values=p_value
    nonDA_cells = nDA_candidates
    DA_cells = DA_candidates
  }
  
  ### Returning Results ----
  # Edge case for full Null composition
  if(all( celltype  %in% nonDA_cells )){DA_cells="Null"}
  # Return Final Results
  return(list("nonDA_cells"=nonDA_cells,"DA_cells"=DA_cells,"Test_Matrix"=pval_matrix,"W"=W,"Search_Results"=Results_Frame,"Confirmatory_Results"=Confirmatory_Results))
}
