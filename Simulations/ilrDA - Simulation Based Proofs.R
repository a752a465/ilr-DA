#################################################################################################################
### Title: ilr-DA: Simulation Based Supports V1 ----
### Current Analyst: Alexander Alsup ----
### Last Updated: (07/10/2026)  ----
### Notes:  ----
#################################################################################################################

# PACKAGES -----
packages <- c("writexl","readxl","ggplot2","MASS","tidyverse","dplyr")
lapply(packages, require, character.only = TRUE)
library(DirichletReg)
#library(logratiolasso)
#library(Hotelling)
#library(harmonicmeanp)
###########################################################################
# CLEAR LIST ----
rm(list=ls())
###########################################################################
# GLOBAL ARGUMENTS ----
Nsims=1000
# sample.size=c(10,20,30,40,50,60,70,80) # This is sample size PER GROUP
sample.size=c(10,40,80,100,150) # Double Check - Delete later in favor of the sample sizes above
scale=50 # A scaling factor for the Poisson means
run_models = c("Poisson"
               ,"ilrDA"
               #,"Dirch_Reg"
)
n_perm=500
it_msg = seq(0,Nsims,by=100)
# Boolean, save plots?
save_plots = FALSE
# Boolean, save results to xlsx?
save_results = FALSE

###########################################################################
# FILEPATHS ----
Figures_path <- "4. Results/Results/Figures/"
Results_path <- "4. Results/Results/xlsx results/"

###########################################################################
# FUNCTIONS ----
source(".../ilrDA - Functions v0.5.R")

celltype_12 <- names(simparams$cell_12_mu_1DA_SmallEffect$cell_mu0)
celltype_5 <- names(simparams$cell_5_mu_1DA_SmallEffect$cell_mu0)
###########################################################################
# Support 1: The W-statistic is different in DA and null Cells -----

### Parameters -----
set.seed(2413)
#cell_mus = toyparams$cell_12_mu_2DA ; scale=1
cell_mus = simparams$cell_12_mu_2DA_SmallEffect
cell_mu0 = cell_mus$cell_mu0*scale ;  cell_mu1 = cell_mus$cell_mu1*scale
N0=50 ; N1=50 ; N=N0+N1
X = c(rep(0,times=N0),rep(1,times=N1))
#
power_cell = which(cell_mu1 != cell_mu0) ; fdr_cell = which(cell_mu1 == cell_mu0)
ncell=length(cell_mu0) ; celltype=names(cell_mu0)

### Simulation ----- 
W_results = data.frame()
SearchOrder_results = data.frame()
for(i in 1:Nsims){
  # Data Generation
  toy_data <- data.generation(ncell=ncell,cell_mu0=cell_mu0,cell_mu1=cell_mu1,N1=N1,N0=N0,type="neg_bino",size_nbinom=25)
  props <- toy_data$Proportions ; counts <- toy_data$Counts
  # ALR Testing
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
  search_order <- order(W, decreasing=TRUE) ; names(search_order) <- celltype
  #
  W_results = bind_rows(W_results,W)
  SearchOrder_results = bind_rows(SearchOrder_results,search_order)
}

### Visualization -----
SearchOrder_results %>%
  dplyr::select(Monocytes,Eosinophils)%>%
  pivot_longer(cols=c("Monocytes","Eosinophils"),names_to="Element",values_to="rank")%>%
  mutate(Element=case_when(Element=="Monocytes"~"DA Element",
                           Element=="Eosinophils"~"Null Element",
                           TRUE~NA))%>%
  group_by(Element)%>%
  count(rank)%>%
  mutate(prop=n/sum(n))%>%
  ungroup()%>%
  ggplot(aes(x=factor(rank),y=prop,fill=Element))+
  geom_bar(stat = "identity" ,color="black") +
  scale_y_continuous(labels = scales::percent) +
  # Labels
  labs(title="Support 1: W-Statistic Ranking Differentiates DA and Null Elements",
       subtitle = paste0("N = 500 simulations (1 = Highest W)"),
       y="Ranking Frequency",
       x="Rankings of Element by W-Statistic",
       )+
  # Facet Wrap
  facet_wrap(~Element,scales="fixed",labeller=label_both)+
  # Theme
  theme(legend.position="top",
        strip.text=element_text(size=9),
        #axis.text.x = element_blank(),
        panel.grid.minor = element_blank()
  )


###########################################################################
# Support 2: Coordinate Permutations at a given level need to ALL reject to define a null set -----

### Parameters -----
set.seed(2413)
#cell_mus = toyparams$cell_12_mu_2DA ; scale=1
cell_mus = simparams$cell_12_mu_2DA_SmallEffect
cell_mu0 = cell_mus$cell_mu0*scale ;  cell_mu1 = cell_mus$cell_mu1*scale
N0=100 ; N1=100 ; N=N0+N1
X = c(rep(0,times=N0),rep(1,times=N1))
#
power_cell = which(cell_mu1 != cell_mu0) ; fdr_cell = which(cell_mu1 == cell_mu0)
ncell=length(cell_mu0) ; celltype=names(cell_mu0)

### Simulation ----- 
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

# Simulate
Coordinate_results = data.frame()
for(i in 1:1000){
  # Data Generation
  toy_data <- data.generation(ncell=ncell,cell_mu0=cell_mu0,cell_mu1=cell_mu1,N1=N1,N0=N0,type="neg_bino",size_nbinom=25)
  props <- toy_data$Proportions ; counts <- toy_data$Counts
  # Coordinate permutations #1 - Contains DA Element
  remaining <- setdiff(celltype,c("Monocytes","Basophils"))
  nums <- remaining
  denoms <- lapply(remaining, function(x) setdiff(remaining, x))
  contrast_list <- list("numerator_set"=as.list(nums), "denom_set"=denoms)
  results_i <- test_contrasts(props, X, contrast_list)
  results_i <- results_i %>% dplyr::select(numerator,p_value) %>% mutate(Set="non-Null Set",iteration=i)
  Coordinate_results = bind_rows(Coordinate_results,results_i)
  # Coordinate permutations #2 - Does not contain DA Element
  remaining <- setdiff(celltype,c("Monocytes","Neutrophils"))
  nums <- remaining
  denoms <- lapply(remaining, function(x) setdiff(remaining, x))
  contrast_list <- list("numerator_set"=as.list(nums), "denom_set"=denoms)
  results_i <- test_contrasts(props, X, contrast_list)
  results_i <- results_i %>% dplyr::select(numerator,p_value) %>% mutate(Set="Null Set",iteration=i)
  Coordinate_results = bind_rows(Coordinate_results,results_i)
  if(i %in% seq(100,1000,100)){
    cat(paste0(i,"|"))
  }
}
head(Coordinate_results)

### Summarize ----
sum_df <- Coordinate_results %>%
  mutate(p_value_binary=ifelse(p_value<=0.05,1,0))

# False Null set when Neutrophils are not included in the coordinate permutations
typei_df <- sum_df %>%
  dplyr::filter(Set=="non-Null Set" & numerator != "Neutrophils")%>%
  group_by(iteration)%>%
  dplyr::summarise(sum_detection=sum(p_value_binary,na.rm=TRUE))%>%
  mutate(sum_detection=ifelse(sum_detection>=1,1,0))
paste0("False null set rate (No Neutrophils): ", 1-mean(typei_df$sum_detection,na.rm=TRUE))

# False Null set when Neutrophils are included in the coordinate permutations
typei_df <- sum_df %>%
  dplyr::filter(Set=="non-Null Set")%>%
  group_by(iteration)%>%
  dplyr::summarise(sum_detection=sum(p_value_binary,na.rm=TRUE))%>%
  mutate(sum_detection=ifelse(sum_detection>=1,1,0))
paste0("False null set rate (Neutrophils): ", 1-mean(typei_df$sum_detection,na.rm=TRUE))

# True Null set true detection rate
typei_df <- sum_df %>%
  dplyr::filter(Set=="Null Set")%>%
  group_by(iteration)%>%
  dplyr::summarise(sum_detection=sum(p_value_binary,na.rm=TRUE))%>%
  mutate(sum_detection=ifelse(sum_detection>=1,1,0))
paste0("True null set rate: ", mean(typei_df$sum_detection,na.rm=TRUE))

### Visualization: Non-Color distribution of P-values from coordinate permutations -----
Coordinate_results %>%
  dplyr::filter(numerator != "Neutrophils")%>%
  mutate(Numerator=case_when(numerator=="Neutrophils"))%>%
  ggplot(aes(x=p_value,fill=numerator))+
  geom_histogram(stat="bin",breaks=seq(0, 1, by=0.05),color="black")+
  # Labels
  labs(title="Distribution of P-values from testing ILR coordinates",
       x="Frequency across simulations",
       y="P-value of Coordinate",
       subtitle = paste0("1000 simulations, N=200 "))+
  # Facet Wrap
  facet_wrap(~Set,scales="fixed",labeller=label_both)+
  # Theme
  theme(legend.position="right",
        strip.text=element_text(size=9),
        #axis.text.x = element_blank(),
        panel.grid.minor = element_blank()
  )

### Visualization: Stacked distribution of P-values from coordinate permutations -----
Coordinate_results %>%
  #dplyr::filter(numerator != "Neutrophils")%>%
  mutate(Numerator=numerator)%>%
  ggplot(aes(x=p_value,fill=Numerator))+
  geom_histogram(stat="bin",position="stack",breaks=seq(0, 1, by=0.05),color="black")+
  coord_cartesian(xlim=c(0.01, 0.99))+
  # Labels
  labs(title="Distribution of P-values from testing ILR coordinates",
       x="Frequency across simulations",
       y="P-value of Coordinate",
       subtitle = paste0("1000 simulations, N=200 "))+
  # Facet Wrap
  facet_wrap(~Set,scales="fixed",labeller=label_both)+
  # Theme
  theme(legend.position="right",
        strip.text=element_text(size=9),
        #axis.text.x = element_blank(),
        panel.grid.minor = element_blank()
  )

###########################################################################
# Support 3: T-statistic distributions against a null set are different for true DA and non-DA  ----
### Notes -----
### Parameters -----
set.seed(2413)
#cell_mus = toyparams$cell_12_mu_2DA ; scale=1
cell_mus = simparams$cell_12_mu_2DA_SmallEffect
celltype = names(cell_mus$cell_mu0)
cell_mu0 = cell_mus$cell_mu0*scale ;  cell_mu1 = cell_mus$cell_mu1*scale
N=200 ; N0=N/2 ; N1=N/2
X = c(rep(0,times=N0),rep(1,times=N1))
n_perm=500
#
power_cell = which(cell_mu1 != cell_mu0) ; fdr_cell = which(cell_mu1 == cell_mu0)
ncell=length(cell_mu0) ; celltype=names(cell_mu0)

### Simulation: Distribution of t-statistics from Random Sampling ----
t1_obs_results = data.frame()
for(i in 1:1000){
  # Data Generation
  toy_data <- data.generation(ncell=ncell,cell_mu0=cell_mu0,cell_mu1=cell_mu1,N1=N1,N0=N0,type="neg_bino",size_nbinom=25)
  props <- toy_data$Proportions ; counts <- toy_data$Counts
  # Distribution of t_obs for True DA candidate w/a clean null set
  DA_candidate <- "Neutrophils"
  nDA_candidates <- setdiff(celltype,c("Eosinophils","Neutrophils","Monocytes"))
  ilr <- ilr_contrast(props, 
                      numerator_set=DA_candidate, 
                      denominator_set=nDA_candidates)
  t1_obs <- t.test(ilr[X==1], ilr[X==0])$statistic
  results_i <- c("Cell Type"="True DA Candidate","T_stat"=t1_obs)
  t1_obs_results <- bind_rows(t1_obs_results,results_i)
  # Distribution of t_obs for False DA candidate w/a clean null set
  DA_candidate <- "Eosinophils"
  nDA_candidates <- setdiff(celltype,c("Eosinophils","Neutrophils","Monocytes"))
  ilr <- ilr_contrast(props, 
                      numerator_set=DA_candidate, 
                      denominator_set=nDA_candidates)
  t1_obs <- t.test(ilr[X==1], ilr[X==0])$statistic
  results_i <- c("Cell Type"="False DA Candidate","T_stat"=t1_obs)
  t1_obs_results <- bind_rows(t1_obs_results,results_i)
  #
  if(i %in% seq(100,1000,100)){
    cat(paste0(i,"|"))
  }
}

### Visualization -----
t1_obs_results %>%
  mutate(T_stat.t=as.numeric(T_stat.t))%>%
  ggplot(aes(x=T_stat.t,fill=`Cell Type`,group=`Cell Type`))+
  geom_histogram(aes(y = after_stat(density)), 
                 stat = "bin", bins = 50, 
                 color = "black", alpha = 0.7,
                 position = "identity")+
  # Overlay theoretical t-distribution
    #stat_function(fun = dt, args = list(df = df),color = "black", linewidth = 1,linetype = "dashed") +
  # Labels
  labs(title="Simulated T-statistic distribution of a False DA candidate and True DA candidate",
       x="T-statistic value",
       y="Frequency",
       subtitle = paste0("1000 simulations, N=",N1*2))+
  # Facet Wrap
  #facet_wrap(~Set,scales="fixed",labeller=label_both)+
  # Theme
  theme(legend.position="right",
        strip.text=element_text(size=9),
        #axis.text.x = element_blank(),
        panel.grid.minor = element_blank()
  )


###########################################################################
# Support 4: Permutation Testing Differentiates true DA candidates and false DA candidates  ----

### Notes -----

### Parameters -----
seed.it <- floor(runif(1,100,1000))
set.seed(981)
#cell_mus = toyparams$cell_12_mu_2DA ; scale=1
cell_mus = simparams$cell_5_mu_1DA_SmallEffect
celltype = names(cell_mus$cell_mu0)
cell_mu0 = cell_mus$cell_mu0*scale ;  cell_mu1 = cell_mus$cell_mu1*scale
N=40 ; N0=N/2 ; N1=N/2
X = c(rep(0,times=N0),rep(1,times=N1))
n_perm=1000
#
power_cell = which(cell_mu1 != cell_mu0) ; fdr_cell = which(cell_mu1 == cell_mu0)
ncell=length(cell_mu0) ; celltype=names(cell_mu0)

### Simulation: Permutation Testing for a single data sample -----
pvals_simulation = data.frame()
Nsims=100
false_DA_cell = "Lymphocytes"

for(i in 1:Nsims){
  #### Data Generation ----
  toy_data <- data.generation(ncell=ncell,cell_mu0=cell_mu0,cell_mu1=cell_mu1,N1=N1,N0=N0,type="poisson",size_nbinom=25)
  props <- toy_data$Proportions ; counts <- toy_data$Counts
  
  #### Permutation Testing values for False DA candidate w/a clean null set ----
  DA_candidate <- false_DA_cell
  nDA_candidates <- setdiff(celltype,c(false_DA_cell,names(power_cell)))
  ilr <- ilr_contrast(props, 
                      numerator_set=DA_candidate, 
                      denominator_set=nDA_candidates)
  t1_obs_false <- t.test(ilr[X==1], ilr[X==0])$statistic
  p1_obs_false <- t.test(ilr[X==1], ilr[X==0])$p.value
  t1_perm_false <- replicate(n_perm, {
    group_perm <- sample(X)
    t.test(ilr[group_perm==1], ilr[group_perm==0])$statistic
  })
  perm_pval = data.frame("Cell"="False DA","Type"="Permutation","Pval"=mean(abs(t1_perm_false) >= abs(t1_obs_false)))
  obs_pval = data.frame("Cell"="False DA","Type"="Parametric","Pval"=p1_obs_false)
  pvals_simulation = bind_rows(pvals_simulation,perm_pval,obs_pval)
  
  #### Permutation Testing values for True DA candidate w/a clean null set ----
  DA_candidate <- names(power_cell)[1]
  nDA_candidates <- setdiff(celltype,c(false_DA_cell,names(power_cell)))
  ilr <- ilr_contrast(props, 
                      numerator_set=DA_candidate, 
                      denominator_set=nDA_candidates)
  t1_obs_true <- t.test(ilr[X==1], ilr[X==0])$statistic
  p1_obs_true <- t.test(ilr[X==1], ilr[X==0])$p.value
  t1_perm_true <- replicate(n_perm, {
    group_perm <- sample(X)
    t.test(ilr[group_perm==1], ilr[group_perm==0])$statistic
  })
  perm_pval = data.frame("Cell"="True DA","Type"="Permutation","Pval"=mean(abs(t1_perm_true) >= abs(t1_obs_true)))
  obs_pval = data.frame("Cell"="True DA","Type"="Parametric","Pval"=p1_obs_true)
  pvals_simulation = bind_rows(pvals_simulation,perm_pval,obs_pval)
  #### Message ----
  if(i %in% seq((Nsims/10),Nsims,length=10)){
    cat(paste0(i,"|"))
  }
}
pvals_simulation = pvals_simulation %>%
  arrange(Cell,Type)

### Summarization -----
pvals_simulation %>% 
  mutate(Pvals_binary=ifelse(Pval<=0.05,1,0))%>%
  group_by(Cell,Type)%>%
  summarize(detection_rate=mean(Pvals_binary,na.rm=TRUE))

### Visualization -----
# Combined Figure
pvals_simulation %>%
  ggplot(aes(x=Pval, fill=Type)) +
  geom_histogram(bins=20, color="black", alpha=0.7, position="identity") +
  # geom_vline(xintercept=t1_obs_true, color="darkblue", linewidth=1.2, linetype="dashed") +
  # geom_vline(xintercept=t1_obs_false, color="darkred", linewidth=1.2, linetype="dashed") +
  #geom_vline(data=obs_df, aes(xintercept=t_obs),color="black", linewidth=1.2, linetype="dashed")+
  facet_grid(Cell~Type, scales="free_y") +
  labs(title="Permutation null distribution vs observed t-statistic",
       x="T-statistic", y="Frequency",
       subtitle="500 permutations, Single representative dataset") +
  theme(legend.position="none",
        panel.grid.minor=element_blank())


###########################################################################