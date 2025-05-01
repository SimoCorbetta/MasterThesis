# Combination therapy single dose in simultaneous.
#We have 4 bacteria species,WT cells, phage resistant cells, antibiotic resistant cells and double resistant cells
#We model the latency period with 10 (n) transit compartments and transit rate equal to n/latency period
#we model collateral sensitivity of phage resistance as decrease of MIC of this strain
# we assume logistic growth curve and same maximum growth rate for all 4 strains(no fitness costs of resistance)
# phage predation modelled as saturated process, antibiotic having a constant (low) decay because we are modelling a in vitro scenario
# we modeled resistance with a fixed rate and we have an higher mutation rate for phage resistance
# than for antibiotic resistance.
rm(list=ls())
require(patchwork)
require(grid)
require(ggplot2) 
require(deSolve) 
require(dplyr)
require(doParallel)
require(tidyr)


output_dir <- "/home/corbettas/IntershipScriptsCTSingleDose_DoubleRes_CE_noFitnessCost_SimpleModel_modThresh_1302"
if (!file.exists(output_dir)) {dir.create (output_dir)}
combination <- function(time,X,parameters)
{
  mu_max_u  = parameters[["mu_max_u"]] # Maximum growth rate of WT cells
  mu_max_rA = parameters[["mu_max_rA"]] # Maximum growth rate for antibiotic resistant strain
  mu_max_rP = parameters[["mu_max_rP"]] # Maximum growth rate for phage resistant strain
  mu_max_rAP = parameters[["mu_max_rAP"]] # Maximum growth rate for double resistant strain
  qA = parameters[["qA"]] # rate of acquisition of antibiotic resistance
  qP = parameters[["qP"]] # rate of acquisition of phage resistance
  B_max   = parameters[["B_max"]] # Maximum carrying capacity
  kmax = parameters[["kmax"]] #Maximum antibiotic killing rate
  MIC = parameters[["MIC"]] # MIC of WT strain
  MIC_rP = parameters[["MIC_rP"]]  # MIC of phage resistant strain
  MIC_rA = parameters[["MIC_rA"]] # MIC of antibiotic and double resistant strain
  H = parameters[["H"]] # Hill coefficient
  y= parameters[["y"]] # antibiotic decay rate
  k = parameters[["k"]] # transit rate (/h)
  beta    = parameters[["beta"]] #linear phage adsorption rate
  P50     = parameters[["P50"]] # [phage] at which phage infection rate is half saturated
  b       = parameters[["b"]] # burst size
  tau     = parameters[["tau"]] # latency period
  Dec_P   = parameters[["Dec_P"]] # decay rate of phages
  
  Bu = X[[1]] # uninfected bacteria
  Bi = X[2:11] # infected bacteria
  Br_A =  X[[12]] # antibiotic resistant bacteria
  Br_i = X[13:22] # antibiotic resistant bacteria infected by phages
  Br_P = X[[23]] # phage resistant bacteria
  Br_AP = X[[24]] # double resistant bacteria
  L = X[[length(X)-2]] # phages
  P = X[[length(X)-1]] # lysed cells
  A = X[[length(X)]] # Antibiotic
  
  Bi_tot = sum(Bi) # total concentration of infected cells sensitive to the antibiotic
  Br_i_tot = sum(Br_i) # total concentration of infected cells resistant to  antibiotic
  
  P_eff = beta/(1+P/P50) # phage predation (Peff)
  B_growth_u = mu_max_u * (1-(Bu+Bi_tot+Br_A+Br_i_tot+Br_P+Br_AP)/B_max) # growth of uninfected WT cells
  B_growth_rA = mu_max_rA * (1-(Bu+Bi_tot+Br_A+Br_i_tot+Br_P+Br_AP)/B_max) # growth of antibiotic resistant cells
  B_growth_rP = mu_max_rP * (1-(Bu+Bi_tot+Br_A+Br_i_tot+Br_P+Br_AP)/B_max) # growth of phage resistant cells
  B_growth_rAP = mu_max_rAP * (1-(Bu+Bi_tot+Br_A+Br_i_tot+Br_P+Br_AP)/B_max) # growth of double resistant
  Keff_u = kmax* ((A/MIC)^H) / (((A/MIC)^H)+ ((kmax-mu_max_u)/mu_max_u)) # antibiotic killing of WT cells
  Keff_rP = kmax* ((A/MIC_rP)^H) / (((A/MIC_rP)^H)+ ((kmax-mu_max_rP)/mu_max_rP)) # antibiotic killing of WT cells
  Keff_rA = kmax* ((A/MIC_rA)^H) / (((A/MIC_rA)^H)+ ((kmax-mu_max_rA)/mu_max_rA)) # antibiotic killing of the antibiotic resistant strain
  Keff_rAP = kmax* ((A/MIC_rA)^H) / (((A/MIC_rA)^H)+ ((kmax-mu_max_rAP)/mu_max_rAP)) # antibiotic killing of the double resistant strain
  dBu = B_growth_u * Bu - Keff_u * Bu - P*P_eff*Bu - Bu*qA - Bu*qP
  dBr_A = B_growth_rA*Br_A + Bu*qA - P_eff*P*Br_A - Br_A*qP - Br_A * Keff_rA
  dBi = c(P_eff * Bu * P - k * Bi[1] - Keff_u* Bi[1]) # equation for newly infected bacteria
  if (length(Bi)>=2)
    for (i in 2:(length(Bi)))
    {
      dBi[i] = k * Bi[i-1] - k * Bi[i] - Keff_u*Bi[i] # equation for the other stages of infection
    }
  dBr_i = c(P_eff * Br_A * P - k * Br_i[1] -  Keff_rA *Br_i[1]  ) # equation for newly infected bacteria
  if (length(Br_i)>=2)
    for (i in 2:(length(Br_i)))
    {
      dBr_i[i] = k * Br_i[i-1] - k * Br_i[i] -Keff_rA* Br_i[i]  # equation for the other stages of infection
    }
  dBr_P = B_growth_rP*Br_P + Bu*qP - Keff_rP*Br_P - Br_P *qA # phage resistant strain
  dBr_AP = B_growth_rAP*Br_AP + qA*Br_P + qP*Br_A - Keff_rAP*Br_AP # double resistant strain
  dP = b * k * (Bi[length(Bi)]+Br_i[length(Br_i)]) - P_eff * (Bu+Br_A+ Bi_tot+Br_i_tot) * P - Dec_P * P # equation for phages
  dL = k * (Bi[length(Bi)]+Br_i[length(Br_i)]) # just to keep track of the number of lysed bacteria in case of "toxicity assays"
  dA = -y*A
  return(list(c(dBu, dBi,dBr_A,dBr_i,dBr_P,dBr_AP,dL,dP, dA)))
}

km<-c(3) # kmax
MIC<-c(1) #MIC of WT cells
MIC_r <- c(10) # MIC of antibiotic and double resistant strain
MIC_rP <- c(0.1,0.5,1) # MIC of phage resistant strain 
h<-c(2) # Hill coeff
d<-c(0.03) # antibiotic decay rate
d_p <- 0.07
m<- 0.7
latent <- c(0.4)
p50<-10^7.5
bmax <- 10^10
ab <- c(10^(-7.40)) # linear phage adsorption rate [1],[9]
burst <- c(150)# burst size tested values
mt_rate_A<- c(10^-9)
mt_rate_P<- c(10^-7)
n<- 10
transit <- n / latent

comb_temp <- expand.grid(kmax=km,MIC=MIC,MIC_rA=MIC_r,MIC_rP=MIC_rP,H=h,
                         mu_max_u  = m,
                         mu_max_rA = m,
                         mu_max_rP = m,
                         mu_max_rAP = m,
                         B_max   = bmax,
                         tau=latent,beta = ab,
                         b = burst,Dec_P = d_p,
                         P50 = p50,
                         k=transit,
                         qA=mt_rate_A,
                         qP=mt_rate_P,
                         y   = d)


params_list <- apply(comb_temp , 1, as.list)
param_list <- lapply(params_list, function(sublist) {
  values <- unlist(sublist)
  numeric_values <- as.numeric(values)
  names(numeric_values) <- names(values)
  return(numeric_values)
})
Bu_wt <- 10^7
rest_of_it <- 0
MOI<- c(0,10^(-4),10^(-3),0.01,0.1,0.5,1,5,10)
ant <- c(0,0.25,0.5,1,1.5,2,3,4,5)
x<- length(MOI)*length(ant)
initial_dosage <- expand.grid(P=MOI*10^7,A=ant)
initial_densities <- c(Bu=Bu_wt,Bi=numeric(n),Br_A=rest_of_it,Br_i=numeric(n),Br_P=rest_of_it,Br_AP=rest_of_it,L=rest_of_it)
initial_densities_mat <- as.data.frame(matrix(rep(initial_densities, x), nrow = x, byrow = TRUE))
colnames(initial_densities_mat) <- names(initial_densities)
initial_values <- as.data.frame(cbind(initial_densities_mat,initial_dosage))
initial_values_list <- apply(initial_values,1,as.list)
initial_list <- lapply(initial_values_list, function(sublist) {
  values <- unlist(sublist)
  numeric_values <- as.numeric(values)
  names(numeric_values) <- names(values)
  return(numeric_values)
})

# solving function
ode_result <-function(initial_values, params){
  times <- seq(0, 72, by = 0.1) # time step of 12 minutes
  
  out_temp <- deSolve::ode(y = initial_values, times = times,
                           func = combination, p = params,method="lsoda")
  out <- cbind(out_temp, 
               as.data.frame(t(params)), 
               as.data.frame(t(initial_values)))
  
  return(out) 
}
# run it on parallel to speed up the computations
#setting core for simulation 
num_cores <- 8
#Creates a set of copies of R running in parallel and communicating over sockets
cl <- makeCluster(num_cores)
#The registerDoParallel function is used to register the parallel backend with the foreach package.
registerDoParallel(cl)
start_exec = Sys.time()
#do parallel simulation

results <- foreach(
  initial_values = initial_list,
  .combine = "rbind"
)%dopar% {
  result <- lapply(param_list,
                   function(params) {
                     ode_result(initial_values, params)
                   })
  
  return(result)
}
stopCluster(cl)
end_exec = Sys.time()
runtime = end_exec - start_exec
print(runtime)
# we combine the results in a dataframe
start_exec_1 = Sys.time()
final_results <- do.call(rbind, lapply(results, function(res) {
  data.frame(res)
}))
end_exec_1 = Sys.time()
run_t_1 = end_exec_1-start_exec_1
print(run_t_1)
final_results$Bi_tot <- rowSums(final_results[, 3:12])
final_results$Br_i_tot <- rowSums(final_results[, 14:23])
final_results$A <- ifelse(final_results$A>10^-4,final_results$A,0)
# function for post processing Br_P
# with ode, we never reach concentration of 0 CFU, so we need to define a cutoff strategies to process the resistant strains
# having  a function only based on the value would not be enough, because at the start of the
# simulation the density of the resistant strain is 0
# we look across the time length of the simulation, and we slide across it with a sliding window.
# if the function is increasing, we go to the next window, if instead it is decreasing(meaning that the cells have been killed)
# and the density is below the threshold, we set it to 0 till the end of simulation
post_process_Br_P <- function(df,w,t){
  df <- df %>%
    filter(time>0)
  n=nrow(df)
  # print(n)
  for (i in seq(1, n, by = w)){
    first_min <- min(df$Br_P[i:w])
    s <-i+w/2 -1
    e <- s + w 
    if (s > n || e > n){
      break
    }
    foll_min <- min(df$Br_P[s  : e ])
    if (first_min< foll_min){
      next
    }
    else {
      if (foll_min< t){
        s0 <- which(df$Br_P[s :e] < t)[1]
        s0 <- s+ s0
        df$Br_P[s0 :n]<- 0
        break}
      else{
        next
      }
    }
  }
  return(df)
}
# post processing function for Br_A,same logic above
post_process_Br_A <- function(df,w,t){
  df <- df %>%
    filter(time>0)
  n=nrow(df)
  # print(n)
  for (i in seq(1, n, by = w)){
    first_min <- min(df$Br_A[i:w])
    s <-i+w/2 -1
    e <- s + w 
    if (s > n || e > n){
      break
    }
    foll_min <- min(df$Br_A[s  : e ])
    if (first_min< foll_min){
      next
    }
    else {
      if (foll_min< t){
        s0 <- which(df$Br_A[s :e] < t)[1]
        s0 <- s+ s0
        df$Br_A[s0 :n]<- 0
        break}
      else{
        next
      }
    }
  }
  return(df)
}
# post processing function for Br_i_tot
post_process_Br_i_tot <- function(df,w,t){
  df <- df %>%
    filter(time>0)
  n=nrow(df)
  # print(n)
  for (i in seq(1, n, by = w)){
    first_min <- min(df$Br_i_tot[i:w])
    s <-i+w/2 -1
    e <- s + w 
    if (s > n || e > n){
      break
    }
    foll_min <- min(df$Br_i_tot[s  : e ])
    if (first_min< foll_min){
      next
    }
    else {
      if (foll_min< t){
        s0 <- which(df$Br_i_tot[s :e] < t)[1]
        s0 <- s+ s0
        df$Br_i_tot[s0 :n]<- 0
        break}
      else{
        next
      }
    }
  }
  return(df)
}
# post processing function for double resistant strain Br_AP
post_process_Br_AP <- function(df,w,t){
  df <- df %>%
    filter(time>0)
  n=nrow(df)
  # print(n)
  for (i in seq(1, n, by = w)){
    first_min <- min(df$Br_AP[i:w])
    s <-i+w/2 -1
    e <- s + w 
    if (s > n || e > n){
      break
    }
    foll_min <- min(df$Br_AP[s  : e ])
    if (first_min< foll_min){
      next
    }
    else {
      if (foll_min< t){
        s0 <- which(df$Br_AP[s :e] < t)[1]
        s0 <- s+ s0
        df$Br_AP[s0 :n]<- 0
        break}
      else{
        next
      }
    }
  }
  return(df)
}
#Plotting functions
# this function creates a dataframe where we can store different quantities such as
# minimum CFU reached, density of single and double resistant at the end of the simulation...
# for each unique combination of phage and antibiotic
create_df <- function(df,mic_rp){
  d<- df %>%
    filter(MIC_rP==mic_rp)%>%
    group_by(P.1, A.1) %>%
    summarise(
      min_CFU = min(CFU, na.rm = TRUE),  
      CFU_at_24 = CFU[time == 24],
      CFU_at_72 = CFU[time == 72],
      Br_A_at_72 = Br_A[time==72],
      Br_P_at_72 = Br_P[time==72],
      Br_AP_at_72 = Br_AP[time==72],
      time_to_min_CFU = time[which.min(CFU)],  
      .groups = 'drop')
  d$TOut <- ifelse(d$CFU_at_72==0|d$min_CFU==0,1,0)
  return(d)}
forced_df <- function(df){
  forced_df<- df
  forced_df$CFU_at_72 <- ifelse(forced_df$min_CFU==0,0,forced_df$CFU_at_72)
  return(forced_df)
}
# plotting function to create heatmap of bacteria eradication for each combination
# of phage and antibiotic dose
TO <- function(df,title){
  heatmap_mild <- df %>%
    dplyr::select(P.1, A.1, TOut) %>%
    pivot_wider(names_from = A.1, values_from = TOut, values_fill = list(TOut = 0))  # fill with 0 if 
  df_long_mild <- heatmap_mild %>%
    pivot_longer(cols = -P.1, names_to = "A.1", values_to = "TOut")
  # Create the heatmap plot
  heatmap_plot_mild <- ggplot(df_long_mild, aes(x = A.1, y = factor(P.1), fill = factor(TOut))) +
    geom_tile(color = "black", size = 0.2) +
    scale_fill_manual(values = c("pink", "darkblue")) +  # Adjust colors (white for 0, darkblue for 1)
    theme_minimal() +
    theme_minimal() +
    labs(x = "Antibiotic Dose (mg)", y = "Phage Dose (PFU/g)", fill = "Bacteria eradication") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1, size = 20),  # Increase size of x-axis labels
      axis.title.x = element_text(size = 20),  # Increase size of x-axis title
      axis.text.y = element_text(size = 14),  # Reduce size of y-axis labels (adjust to your preference)
      axis.title.y = element_text(size = 24),  # Increase size of y-axis title (adjust to your preference)
      panel.grid.major = element_line(color = "black", size = 0.5),  # Major grid lines
      panel.grid.minor = element_line(color = "gray", size = 0.5)  # Minor grid lines
    )
  jpeg(filename = paste0(output_dir,title), 
       width = 1500, height = 1300, res = 150)  # Increased resolution for better quality
  print(heatmap_plot_mild)
  dev.off()
  return(heatmap_plot_mild)
}
# plotting function to create heatmap of log density of total bacteria population for each combination
# of phage and antibiotic dose
Log_CFU <- function(df,forced_df,title){
  df_long_mild_doubleres <- df %>%
    mutate(log_CFU_at_72 = ifelse(forced_df$CFU_at_72==0,-1,log10(CFU_at_72 + 1)))
  heatmap_mild_res <- df_long_mild_doubleres %>%
    dplyr::select(P.1, A.1,log_CFU_at_72) %>%
    pivot_wider(names_from = A.1, values_from = log_CFU_at_72)  # fill with 0 if 
  df_long_mild <- heatmap_mild_res %>%
    pivot_longer(cols = -P.1, names_to = "A.1", values_to = "log_CFU_at_72")
  df_long_mild$A.1 <- as.numeric(df_long_mild$A.1)
  df_long_mild <- df_long_mild %>%
    mutate(log_density = case_when(
      log_CFU_at_72 == -1 ~ "bacteria eradication",
      log_CFU_at_72 >= 0 & log_CFU_at_72 < 1 ~ "1",
      log_CFU_at_72 >= 1 & log_CFU_at_72 < 2 ~ "2",
      log_CFU_at_72 >= 2 & log_CFU_at_72 < 3 ~ "3",
      log_CFU_at_72 >= 3 & log_CFU_at_72 < 4 ~ "4",
      log_CFU_at_72 >= 4 & log_CFU_at_72 < 5 ~ "5",
      log_CFU_at_72 >= 5 & log_CFU_at_72 < 6 ~ "6",
      log_CFU_at_72 >= 6 & log_CFU_at_72 < 7 ~ "7",
      log_CFU_at_72 >= 7 & log_CFU_at_72 < 8 ~ "8",
      log_CFU_at_72 >= 8 & log_CFU_at_72 < 9 ~ "9",
      log_CFU_at_72 >= 9 & log_CFU_at_72  ~ "10",
      TRUE ~ "unknown"
    ))
  df_long_mild$log_density <- factor(df_long_mild$log_density, 
                                     levels = c(
                                       "bacteria eradication","1" ,
                                       "2", "3", "4", "5", "6", "7","8","9", "10"
                                     ))
  heatmap_plot_mild <- ggplot(df_long_mild, aes(x = factor(A.1), y = factor(P.1), fill = log_density)) +
    geom_tile(color = "black", size = 0.2) +
    scale_fill_manual(values = c("bacteria eradication" = "darkblue","1"="darkolivegreen1", "2" = "seagreen4","3"="lightyellow","4"="yellow","5"="orange","6"="orange3","7"="sienna","8"="sienna4","9"="brown", "10"="red")) +
    theme_minimal() +
    labs(x = "Antibiotic Dose (mg)", y = "Phage Dose (PFU/g)", fill = "Log density of total bacteria population") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1, size = 20),
      axis.title.x = element_text(size = 20),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 24),
      #legend.text = element_text(size = 14),
      #legend.title = element_text(size = 16),
      panel.grid.major = element_line(color = "black", size = 0.5),
      panel.grid.minor = element_line(color = "gray", size = 0.5)
    )
  # Save the heatmap as a JPEG
  jpeg(filename = paste0(output_dir,title), 
       width = 1800, height = 1600, res = 150)  # Adjusted resolution for better quality
  print(heatmap_plot_mild)
  dev.off()
  return(heatmap_plot_mild)
}
# here i apply my processing functions to the unprocessed dataframe
final_results_processed <- final_results %>%
  group_by(MIC_rP, A.1, P.1) %>%
  # Apply the post_process_Br function to each group
  group_modify(~ {
    df <- .x  # The dataframe for each group
    # Apply the post_process_Br function to the dataframe
    df <- post_process_Br_P(df,30,1)
    df <- post_process_Br_A(df,30,1)
    df <- post_process_Br_i_tot(df,30,1)
    df <- post_process_Br_AP(df,30,1)
    df <- df %>%
      mutate(Bu=case_when(
        Bu<1 ~ 0, 
        TRUE ~ Bu))
    for (i in 1:nrow(df)) {
      # Check if Bu, Br_A, and Br_P are 0, and Br_AP is below the threshold
      if (df$Bu[i] == 0 && df$Br_A[i] == 0 && df$Br_P[i] == 0 && df$Br_AP[i] < 1) {
        # Set Br_AP to 0 from this row till the end of the dataframe
        df$Br_AP[i:nrow(df)] <- 0
        break  # Exit the loop once the condition is met
      }
    }# this is needed for plotting purposes to force Br_AP to 0 in case of bacteria eradication
    # Add the CFU column based on the given condition
    df <- df %>%
      mutate(CFU = case_when(
        Bu + Bi_tot + Br_A + Br_i_tot+Br_P+Br_AP <= 1 ~ 0, 
        TRUE ~ Bu + Bi_tot + Br_A +Br_i_tot+Br_P+Br_AP       
      ))
    
    
    return(df)
  }) %>%
  ungroup()