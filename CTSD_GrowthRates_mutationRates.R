# in this case we want to evaluate the impact of fitness costs of resistance and mutation rates,
# in the setting of simultaneous administration of phage-antibiotic
rm(list=ls())
require(patchwork)
require(grid)
require(ggplot2) #plot library
require(deSolve) #ode solvers
require(dplyr)
require(doParallel)
require(tidyr)

#Implementation of fitness costs
## associated with resistance there s a fitness cost, modelled as reduction of the max growth rate
# we ll have 10 %,25 % and 40 % fold reduction
# we ll also have a scenario where double resistance  is more penalized than single resistance

#4 scenarios concerning the mutation rates :
# 1) qA = 10^-10/h & qP = 10^-8/h
# 2)  qA = qP = 10^-9 /h
#3) qA = 10^-8/h & qP = 10^-6/h
#4) qA = qP = 10^-7/h
# we have qP >= qA because in phage therapy is almost certain the emergence of phage resistance
# reference value of qP = 10^-7/h  & qA=10^-9
# we will evaluate these 2 parameters separately, meaning that when we vary the growth rates
# of the different strain we will use the reference values for the mutation rates,
# and when evaluate the mutation rates, we will not have any  fitness cost of resistance
output_dir <- "/home/corbettas/CTSD_CE_fitCosts_mutRates"
if (!file.exists(output_dir)) {dir.create (output_dir)}
combination <- function(time,X,parameters)
{
  mu_max_u  = parameters[["mu_max_u"]] # Maximum growth rate
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
  dP = b * k * (Bi[length(Bi)]+Br_i[length(Br_i)]) - P_eff * (Bu+Br_A+ Bi_tot+Br_i_tot) * P - Dec_P * P # equation for 
  dL = k * (Bi[length(Bi)]+Br_i[length(Br_i)]) # just to keep track of the number of lysed bacteria in case of "toxicity assays" for instance
  dA = -y*A
  return(list(c(dBu, dBi,dBr_A,dBr_i,dBr_P,dBr_AP,dL,dP, dA)))
}
km<-c(3) # kmax
MIC<-c(1) #MIC of WT cells
MIC_r <- c(10) # MIC of antibiotic and double resistant strain
MIC_rP <- c(0.1,0.5,1) # MIC of phage resistant strain 
h<-c(2) # Hill coeff
d<-c(0.03) # antibiotic decay rate
latent <- c(0.4)
ab <- c(10^(-7.4)) 
burst <- 150
m <- 0.7
bmax <- 10^10
p50 <-10^7.5 
dec_p <- 0.07
mx_r <- c(0.63,0.525,0.42)
qa <- c(10^-10,10^-9,10^-8,10^-7)
qp <- c(10^-9,10^-8,10^-7,10^-6)
n<- 10
transit <- n / latent

comb_temp <- expand.grid(kmax=km,MIC=MIC,MIC_rA=MIC_r,MIC_rP=MIC_rP,H=h,
                         mu_max_u  = m,
                         mu_max_rA = mx_r ,
                         mu_max_rP = mx_r,
                         mu_max_rAP = mx_r,
                         B_max   = bmax,
                         tau=latent,beta = ab,
                         b = burst,Dec_P = dec_p,
                         P50 = p50 ,
                         k=transit,
                         qA=qa,
                         qP=qp,
                         y   = d)

comb_temp <- comb_temp %>%
  filter(
    (mu_max_rA == 0.63 & mu_max_rP ==0.63 &  mu_max_rAP ==0.63) |
      (mu_max_rA == 0.525 & mu_max_rP ==0.525 &  mu_max_rAP ==0.525 ) |
      (mu_max_rA == 0.42 & mu_max_rP ==0.42 &  mu_max_rAP ==0.42 ) |
      (mu_max_rA == 0.63 & mu_max_rP ==0.63 &  mu_max_rAP ==0.525 )
  )
comb_temp <- comb_temp %>% filter(
  (qA == 10^-10 & qP ==10^-8) |
    (qA == 10^-9 & qP ==10^-9) |
    (qA == 10^-8 & qP ==10^-6)|
    (qA == 10^-7 & qP ==10^-7)
)
params_list <- apply(comb_temp , 1, as.list)
param_list <- lapply(params_list, function(sublist) {
  values <- unlist(sublist)
  numeric_values <- as.numeric(values)
  names(numeric_values) <- names(values)
  return(numeric_values)
})
MOI<- c(0,10^(-4),10^(-3),0.01,0.1,0.5,1,5,10)
ant <- c(0,0.25,0.5,1,1.5,2,3,4,5)
bu <- 10^7
the_rest <- 0
x <- length(MOI) *length(ant)
initial_dosage <- expand.grid(P=MOI*10^7,A=ant)
initial_densities <- c(Bu=bu,Bi=numeric(n),Br_A=the_rest,Br_i=numeric(n),Br_P=the_rest,Br_AP=the_rest,L=the_rest)
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
  times <- seq(0, 72, by = 0.1) # time step di circa 12 minuti
  
  out_temp <- deSolve::ode(y = initial_values, times = times,
                           func = combination, p = params)
  out <- cbind(out_temp, 
               as.data.frame(t(params)), 
               as.data.frame(t(initial_values)))
  
  return(out) 
}

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

# stop parallel
stopCluster(cl)
end_exec = Sys.time()
runtime = end_exec - start_exec
print(runtime)
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
# same 4 post processing functions
# function for post processing Br_P
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
# post processing function for Br_A
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
# apply the post processing to the unprocessed dataframe
start_exec_2 = Sys.time()
final_results_processed <- final_results %>%
  group_by(mu_max_rA,mu_max_rP,mu_max_rAP,MIC_rP,qA,qP, A.1, P.1) %>%
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
    }
    
    # Add the CFU column based on the given condition
    df <- df %>%
      mutate(CFU = case_when(
        Bu + Bi_tot + Br_A + Br_i_tot+Br_P+Br_AP <= 1 ~ 0,  # If the sum of Bu, Bi_tot, Br is less than or equal to 10^-2, set CFU to 0
        TRUE ~ Bu + Bi_tot + Br_A +Br_i_tot+Br_P+Br_AP       # Otherwise, calculate CFU as the sum
      ))
   
    
    return(df)
  }) %>%
  ungroup()

end_exec_2 = Sys.time()
run_t_2 <- end_exec_2 - start_exec_2
print(run_t_2)

# plotting function
plot_CFU <- function(df,mic_rp,ma,mp,map,qa,qp,title){
  pl<-df %>%
    filter(MIC_rP==mic_rp)%>%
    filter(mu_max_rA==ma)%>%
    filter(mu_max_rP==mp)%>%
    filter(mu_max_rAP==map)%>%
    filter(qA==qa)%>%
    filter(qP==qp) %>%
    ggplot(aes(x = time, y = CFU)) +
    geom_hline(yintercept = 1000, linetype = "dashed", color = "black",size=2)+
    geom_hline(yintercept = 1, linetype = "dashed", color = "green",size=2)+
    facet_grid(A.1~P.1,scales="fixed")+
    geom_line(linewidth=2)+
    theme_bw(base_size=60)+
    scale_y_log10()+
    scale_x_continuous(breaks=c(0,24,48,72))+
    labs(x="time (h)",y="CFU (CFU/g)")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(size = 40),  # Increase x-axis label size
          axis.title.y = element_text(size = 40))
  combined_plot <- pl
  jpeg(filename = paste0(output_dir,title), 
       width = 2000, height = 1500, res = 50)
  print(combined_plot)
  dev.off()
  return(pl)
}

plot_Res <- function(df,mic_rp,ma,mp,map,qa,qp,title){
  pl<-df %>%
    filter(MIC_rP==mic_rp)%>%
    filter(mu_max_rA==ma)%>%
    filter(mu_max_rP==mp)%>%
    filter(mu_max_rAP==map)%>%
    filter(qA==qa)%>%
    filter(qP==qp) %>%
    pivot_longer(cols = c( Br_A,Br_P,Br_AP), 
                 names_to = "Variable", 
                 values_to = "Value") %>%  # Reshape the data to long format
    ggplot(aes(x = time, y = Value, color = Variable)) +
    geom_hline(yintercept = 10^-2, linetype = "dashed", color = "green", size = 2) +
    facet_grid(A.1 ~ P.1, scales = "fixed") +  # Use facet grid if you want
    geom_line(linewidth = 2) +
    theme_bw(base_size = 60) +
    scale_y_log10() +
    scale_x_continuous(breaks = c(0, 24, 48, 72)) +
    labs(x = "Time (h)", y = "Value") +  # Update label for y-axis
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(size = 40),  # Increase x-axis label size
          axis.title.y = element_text(size = 40)) # Increase y-axis label size
  combined_plot <- pl
  jpeg(filename = paste0(output_dir,title), 
       width = 2000, height = 1500, res = 50)
  grid.draw(combined_plot)
  dev.off()
  return(pl)
}
create_df <- function(df,mic_rp,ma,mp,map,qa,qp){
  d<- df %>%
    filter(MIC_rP==mic_rp)%>%
    filter(mu_max_rA==ma)%>%
    filter(mu_max_rP==mp)%>%
    filter(mu_max_rAP==map)%>%
    filter(qA==qa)%>%
    filter(qP==qp) %>%
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
  d$BacteriaEradication <- ifelse(d$CFU_at_72==0|d$min_CFU==0,1,0)
  return(d)}

forced_df <- function(df){
  forced_df<- df
  forced_df$CFU_at_72 <- ifelse(forced_df$min_CFU==0,0,forced_df$CFU_at_72)
  return(forced_df)
}

TO <- function(df,title){
  heatmap_mild <- df %>%
    dplyr::select(P.1, A.1, BacteriaEradication) %>%
    pivot_wider(names_from = A.1, values_from = BacteriaEradication, values_fill = list(BacteriaEradication = 0))  # fill with 0 if 
  df_long_mild <- heatmap_mild %>%
    pivot_longer(cols = -P.1, names_to = "A.1", values_to = "BacteriaEradication")
  # Create the heatmap plot
  heatmap_plot_mild <- ggplot(df_long_mild, aes(x = A.1, y = factor(P.1), fill = factor(BacteriaEradication))) +
    geom_tile(color = "black", size = 0.2) +
    scale_fill_manual(values = c("pink", "darkblue")) +  # Adjust colors (white for 0, darkblue for 1)
    theme_minimal() +
    theme_minimal() +
    labs(x = "Antibiotic Dose (mg)", y = "Phage Dose (PFU/g)", fill = "Bacteria Eradication") +
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

  
