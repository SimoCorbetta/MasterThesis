#Sequential dosing 
# in this setting, we will administer one treatment after the other, varying the time delay 
# and order of administration of treatment. Model structure and parameter values are exactly the same
# the only thing that changes is the implementation of the model because rxode2 needs the explicit equation for each single compartment.
# In this setting, it would require more work to change the number of transit compartments
library(rxode2)
library(tidyverse)
library(grid)
library(doParallel)
library(tidyr)
rm(list=ls())

output_dir <- "/home/corbettas/CTSequentialDrugAdm_DoubleRes_CE_def"
if (!file.exists(output_dir)) {dir.create (output_dir)}
# sequential treatment administration
# absence of CS (MIC of phage resistant strain = MIC of WT strain)
mod_noCS <- function(burst=150) {
  # Initial conditions and parameters
  ini({
    kmax <- 3 # kmax
    MIC <- 1 # MIC of WT strain
    MIC_r <- 10 # MIC of antibiotic resistant strain
    H <- 2 # Hill coefficient
    y <- 0.03 # antibiotic decay rate
    beta <- 10^(-7.40) # linear phage adsorption rate
    b <- burst # burst size
    mu_max_u  = 0.7 # max growth rate of WT strain
    mu_max_rA = 0.7 # max growth rate of antibiotic resistant strain
    mu_max_rP = 0.7 # max growth rate of phage resistant strain
    mu_max_rAP = 0.7 # max growth rate of double resistant strain
    qA = 10^-9 # acquisition rate of antibiotic resistance
    qP = 10^-7 # acquisition rate of phage resistance
    B_max   = 10^10
    Dec_P = 0.07
    P50 = 10^7.5
    k = 10 / 0.4
  })
  
  model({
    # Variables for the differential equations
    Bi_tot <- Bi_1 + Bi_2 + Bi_3 + Bi_4 + Bi_5 + Bi_6 + Bi_7 + Bi_8 + Bi_9 + Bi_10
    Br_i_tot <- Bri_1 + Bri_2 + Bri_3 + Bri_4 + Bri_5 + Bri_6 + Bri_7 + Bri_8 + Bri_9 + Bri_10
    
    P_eff <- beta / (1 + P / P50) # phage predation (Peff)
    B_growth_u <- mu_max_u * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of uninfected WT cells
    B_growth_rA <- mu_max_rA * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of antibiotic resistant cells
    B_growth_rP <- mu_max_rP * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of phage resistant cells
    B_growth_rAP <- mu_max_rAP * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of double resistant
    Keff_rP <- kmax * ((A / MIC)^H) / (((A / MIC)^H) + ((kmax - mu_max_rP) / mu_max_rP)) # antibiotic killing of the phage resistant strain
    Keff_u <- kmax * ((A / MIC)^H) / (((A / MIC)^H) + ((kmax - mu_max_u) / mu_max_u)) # antibiotic killing of WT cells
    Keff_rA <- kmax * ((A / MIC_r)^H) / (((A / MIC_r)^H) + ((kmax - mu_max_rA) / mu_max_rA)) # antibiotic killing of the antibiotic resistant strain
    Keff_rAP <- kmax * ((A / MIC_r)^H) / (((A / MIC_r)^H) + ((kmax - mu_max_rAP) / mu_max_rAP)) # antibiotic killing of the double resistant strain
    
    # Differential equations
    d/dt(Bu) <- B_growth_u * Bu - Keff_u * Bu - P * P_eff * Bu - Bu * qA - Bu * qP
    
    # Explicit equations for infected bacteria (Bi_1 to Bi_10)
    d/dt(Bi_1) <- P_eff * Bu * P - k * Bi_1 - Keff_u * Bi_1
    d/dt(Bi_2) <- k * Bi_1 - k * Bi_2 - Keff_u * Bi_2
    d/dt(Bi_3) <- k * Bi_2 - k * Bi_3 - Keff_u * Bi_3
    d/dt(Bi_4) <- k * Bi_3 - k * Bi_4 - Keff_u * Bi_4
    d/dt(Bi_5) <- k * Bi_4 - k * Bi_5 - Keff_u * Bi_5
    d/dt(Bi_6) <- k * Bi_5 - k * Bi_6 - Keff_u * Bi_6
    d/dt(Bi_7) <- k * Bi_6 - k * Bi_7 - Keff_u * Bi_7
    d/dt(Bi_8) <- k * Bi_7 - k * Bi_8 - Keff_u * Bi_8
    d/dt(Bi_9) <- k * Bi_8 - k * Bi_9 - Keff_u * Bi_9
    d/dt(Bi_10) <- k * Bi_9 - k * Bi_10 - Keff_u * Bi_10
    
    d/dt(Br_A) <- B_growth_rA * Br_A + Bu * qA - P_eff * P * Br_A - Br_A * qP - Br_A * Keff_rA
    
    # Explicit equations for infected antibiotic-resistant bacteria (Bri_1 to Bri_10)
    d/dt(Bri_1) <- P_eff * Br_A * P - k * Bri_1 - Keff_rA * Bri_1
    d/dt(Bri_2) <- k * Bri_1 - k * Bri_2 - Keff_rA * Bri_2
    d/dt(Bri_3) <- k * Bri_2 - k * Bri_3 - Keff_rA * Bri_3
    d/dt(Bri_4) <- k * Bri_3 - k * Bri_4 - Keff_rA * Bri_4
    d/dt(Bri_5) <- k * Bri_4 - k * Bri_5 - Keff_rA * Bri_5
    d/dt(Bri_6) <- k * Bri_5 - k * Bri_6 - Keff_rA * Bri_6
    d/dt(Bri_7) <- k * Bri_6 - k * Bri_7 - Keff_rA * Bri_7
    d/dt(Bri_8) <- k * Bri_7 - k * Bri_8 - Keff_rA * Bri_8
    d/dt(Bri_9) <- k * Bri_8 - k * Bri_9 - Keff_rA * Bri_9
    d/dt(Bri_10) <- k * Bri_9 - k * Bri_10 - Keff_rA * Bri_10
    
    d/dt(Br_P) <- B_growth_rP * Br_P + Bu * qP - Keff_rP * Br_P - Br_P * qA
    d/dt(Br_AP) <- B_growth_rAP * Br_AP + qA * Br_P + qP * Br_A - Keff_rAP * Br_AP
    d/dt(P) <- b * k * (Bi_10 + Bri_10) - P_eff * (Bu + Br_A + Bi_tot + Br_i_tot) * P - Dec_P * P
    d/dt(A) <- -y * A
    
    # Initial conditions
    Bu(0) <- 10^7
    Bi_1(0) <- 0
    Bi_2(0) <- 0
    Bi_3(0) <- 0
    Bi_4(0) <- 0
    Bi_5(0) <- 0
    Bi_6(0) <- 0
    Bi_7(0) <- 0
    Bi_8(0) <- 0
    Bi_9(0) <- 0
    Bi_10(0) <- 0
    Br_A(0) <- 0
    Bri_1(0) <- 0
    Bri_2(0) <- 0
    Bri_3(0) <- 0
    Bri_4(0) <- 0
    Bri_5(0) <- 0
    Bri_6(0) <- 0
    Bri_7(0) <- 0
    Bri_8(0) <- 0
    Bri_9(0) <- 0
    Bri_10(0) <- 0
    Br_P(0) <- 0
    Br_AP(0) <- 0
    P(0) <- 0
    A(0) <- 0
  })
}

###  weak CS (MIC of phage resistant strain = 0.5 * MIC WT)
mod_weakCS <- function() {
  # Initial conditions and parameters
  ini({
    kmax <- 3 # kmax
    MIC <- 1 # MIC of WT strain
    MIC_rP <- 0.5
    MIC_r <- 10 # MIC of antibiotic resistant strain
    H <- 2 # Hill coefficient
    y <- 0.03 # antibiotic decay rate
    beta <- 10^(-7.40) # linear phage adsorption rate
    b <- 150 # burst size
    mu_max_u  = 0.7 # max growth rate of WT strain
    mu_max_rA = 0.7 # max growth rate of antibiotic resistant strain
    mu_max_rP = 0.7 # max growth rate of phage resistant strain
    mu_max_rAP = 0.7 # max growth rate of double resistant strain
    qA = 10^-9 # acquisition rate of antibiotic resistance
    qP = 10^-7 # acquisition rate of phage resistance
    B_max   = 10^10
    Dec_P = 0.07
    P50 = 10^7.5
    k = 10 / 0.4
  })
  
  model({
    # Variables for the differential equations
    Bi_tot <- Bi_1 + Bi_2 + Bi_3 + Bi_4 + Bi_5 + Bi_6 + Bi_7 + Bi_8 + Bi_9 + Bi_10
    Br_i_tot <- Bri_1 + Bri_2 + Bri_3 + Bri_4 + Bri_5 + Bri_6 + Bri_7 + Bri_8 + Bri_9 + Bri_10
    
    P_eff <- beta / (1 + P / P50) # phage predation (Peff)
    B_growth_u <- mu_max_u * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of uninfected WT cells
    B_growth_rA <- mu_max_rA * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of antibiotic resistant cells
    B_growth_rP <- mu_max_rP * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of phage resistant cells
    B_growth_rAP <- mu_max_rAP * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of double resistant
    Keff_rP <- kmax * ((A / MIC_rP)^H) / (((A / MIC_rP)^H) + ((kmax - mu_max_rP) / mu_max_rP)) # antibiotic killing of the phage resistant strain
    Keff_u <- kmax * ((A / MIC)^H) / (((A / MIC)^H) + ((kmax - mu_max_u) / mu_max_u)) # antibiotic killing of WT cells
    Keff_rA <- kmax * ((A / MIC_r)^H) / (((A / MIC_r)^H) + ((kmax - mu_max_rA) / mu_max_rA)) # antibiotic killing of the antibiotic resistant strain
    Keff_rAP <- kmax * ((A / MIC_r)^H) / (((A / MIC_r)^H) + ((kmax - mu_max_rAP) / mu_max_rAP)) # antibiotic killing of the double resistant strain
    
    # Differential equations
    d/dt(Bu) <- B_growth_u * Bu - Keff_u * Bu - P * P_eff * Bu - Bu * qA - Bu * qP
    
    # Explicit equations for infected bacteria (Bi_1 to Bi_10)
    d/dt(Bi_1) <- P_eff * Bu * P - k * Bi_1 - Keff_u * Bi_1
    d/dt(Bi_2) <- k * Bi_1 - k * Bi_2 - Keff_u * Bi_2
    d/dt(Bi_3) <- k * Bi_2 - k * Bi_3 - Keff_u * Bi_3
    d/dt(Bi_4) <- k * Bi_3 - k * Bi_4 - Keff_u * Bi_4
    d/dt(Bi_5) <- k * Bi_4 - k * Bi_5 - Keff_u * Bi_5
    d/dt(Bi_6) <- k * Bi_5 - k * Bi_6 - Keff_u * Bi_6
    d/dt(Bi_7) <- k * Bi_6 - k * Bi_7 - Keff_u * Bi_7
    d/dt(Bi_8) <- k * Bi_7 - k * Bi_8 - Keff_u * Bi_8
    d/dt(Bi_9) <- k * Bi_8 - k * Bi_9 - Keff_u * Bi_9
    d/dt(Bi_10) <- k * Bi_9 - k * Bi_10 - Keff_u * Bi_10
    
    d/dt(Br_A) <- B_growth_rA * Br_A + Bu * qA - P_eff * P * Br_A - Br_A * qP - Br_A * Keff_rA
    
    # Explicit equations for infected antibiotic-resistant bacteria (Bri_1 to Bri_10)
    d/dt(Bri_1) <- P_eff * Br_A * P - k * Bri_1 - Keff_rA * Bri_1
    d/dt(Bri_2) <- k * Bri_1 - k * Bri_2 - Keff_rA * Bri_2
    d/dt(Bri_3) <- k * Bri_2 - k * Bri_3 - Keff_rA * Bri_3
    d/dt(Bri_4) <- k * Bri_3 - k * Bri_4 - Keff_rA * Bri_4
    d/dt(Bri_5) <- k * Bri_4 - k * Bri_5 - Keff_rA * Bri_5
    d/dt(Bri_6) <- k * Bri_5 - k * Bri_6 - Keff_rA * Bri_6
    d/dt(Bri_7) <- k * Bri_6 - k * Bri_7 - Keff_rA * Bri_7
    d/dt(Bri_8) <- k * Bri_7 - k * Bri_8 - Keff_rA * Bri_8
    d/dt(Bri_9) <- k * Bri_8 - k * Bri_9 - Keff_rA * Bri_9
    d/dt(Bri_10) <- k * Bri_9 - k * Bri_10 - Keff_rA * Bri_10
    
    d/dt(Br_P) <- B_growth_rP * Br_P + Bu * qP - Keff_rP * Br_P - Br_P * qA
    d/dt(Br_AP) <- B_growth_rAP * Br_AP + qA * Br_P + qP * Br_A - Keff_rAP * Br_AP
    d/dt(P) <- b * k * (Bi_10 + Bri_10) - P_eff * (Bu + Br_A + Bi_tot + Br_i_tot) * P - Dec_P * P
    d/dt(A) <- -y * A
    
    # Initial conditions
    Bu(0) <- 10^7
    Bi_1(0) <- 0
    Bi_2(0) <- 0
    Bi_3(0) <- 0
    Bi_4(0) <- 0
    Bi_5(0) <- 0
    Bi_6(0) <- 0
    Bi_7(0) <- 0
    Bi_8(0) <- 0
    Bi_9(0) <- 0
    Bi_10(0) <- 0
    Br_A(0) <- 0
    Bri_1(0) <- 0
    Bri_2(0) <- 0
    Bri_3(0) <- 0
    Bri_4(0) <- 0
    Bri_5(0) <- 0
    Bri_6(0) <- 0
    Bri_7(0) <- 0
    Bri_8(0) <- 0
    Bri_9(0) <- 0
    Bri_10(0) <- 0
    Br_P(0) <- 0
    Br_AP(0) <- 0
    P(0) <- 0
    A(0) <- 0
  })
}
# strong CS (MIC of phage resistant strain = 0.1 * MIC WT)
mod_strongCS <- function() {
  # Initial conditions and parameters
  ini({
    kmax <- 3 # kmax
    MIC <- 1 # MIC of WT strain
    MIC_rP <- 0.1
    MIC_r <- 10 # MIC of antibiotic resistant strain
    H <- 2 # Hill coefficient
    y <- 0.03 # antibiotic decay rate
    beta <- 10^(-7.40) # linear phage adsorption rate
    b <- 150 # burst size
    mu_max_u  = 0.7 # max growth rate of WT strain
    mu_max_rA = 0.7 # max growth rate of antibiotic resistant strain
    mu_max_rP = 0.7 # max growth rate of phage resistant strain
    mu_max_rAP = 0.7 # max growth rate of double resistant strain
    qA = 10^-9 # acquisition rate of antibiotic resistance
    qP = 10^-7 # acquisition rate of phage resistance
    B_max   = 10^10
    Dec_P = 0.07
    P50 = 10^7.5
    k = 10 / 0.4
  })
  
  model({
    # Variables for the differential equations
    Bi_tot <- Bi_1 + Bi_2 + Bi_3 + Bi_4 + Bi_5 + Bi_6 + Bi_7 + Bi_8 + Bi_9 + Bi_10
    Br_i_tot <- Bri_1 + Bri_2 + Bri_3 + Bri_4 + Bri_5 + Bri_6 + Bri_7 + Bri_8 + Bri_9 + Bri_10
    
    P_eff <- beta / (1 + P / P50) # phage predation (Peff)
    B_growth_u <- mu_max_u * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of uninfected WT cells
    B_growth_rA <- mu_max_rA * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of antibiotic resistant cells
    B_growth_rP <- mu_max_rP * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of phage resistant cells
    B_growth_rAP <- mu_max_rAP * (1 - (Bu + Bi_tot + Br_A + Br_i_tot + Br_P + Br_AP) / B_max) # growth of double resistant
    Keff_rP <- kmax * ((A / MIC_rP)^H) / (((A / MIC_rP)^H) + ((kmax - mu_max_rP) / mu_max_rP)) # antibiotic killing of the phage resistant strain
    Keff_u <- kmax * ((A / MIC)^H) / (((A / MIC)^H) + ((kmax - mu_max_u) / mu_max_u)) # antibiotic killing of WT cells
    Keff_rA <- kmax * ((A / MIC_r)^H) / (((A / MIC_r)^H) + ((kmax - mu_max_rA) / mu_max_rA)) # antibiotic killing of the antibiotic resistant strain
    Keff_rAP <- kmax * ((A / MIC_r)^H) / (((A / MIC_r)^H) + ((kmax - mu_max_rAP) / mu_max_rAP)) # antibiotic killing of the double resistant strain
    
    # Differential equations
    d/dt(Bu) <- B_growth_u * Bu - Keff_u * Bu - P * P_eff * Bu - Bu * qA - Bu * qP
    
    # Explicit equations for infected bacteria (Bi_1 to Bi_10)
    d/dt(Bi_1) <- P_eff * Bu * P - k * Bi_1 - Keff_u * Bi_1
    d/dt(Bi_2) <- k * Bi_1 - k * Bi_2 - Keff_u * Bi_2
    d/dt(Bi_3) <- k * Bi_2 - k * Bi_3 - Keff_u * Bi_3
    d/dt(Bi_4) <- k * Bi_3 - k * Bi_4 - Keff_u * Bi_4
    d/dt(Bi_5) <- k * Bi_4 - k * Bi_5 - Keff_u * Bi_5
    d/dt(Bi_6) <- k * Bi_5 - k * Bi_6 - Keff_u * Bi_6
    d/dt(Bi_7) <- k * Bi_6 - k * Bi_7 - Keff_u * Bi_7
    d/dt(Bi_8) <- k * Bi_7 - k * Bi_8 - Keff_u * Bi_8
    d/dt(Bi_9) <- k * Bi_8 - k * Bi_9 - Keff_u * Bi_9
    d/dt(Bi_10) <- k * Bi_9 - k * Bi_10 - Keff_u * Bi_10
    
    d/dt(Br_A) <- B_growth_rA * Br_A + Bu * qA - P_eff * P * Br_A - Br_A * qP - Br_A * Keff_rA
    
    # Explicit equations for infected antibiotic-resistant bacteria (Bri_1 to Bri_10)
    d/dt(Bri_1) <- P_eff * Br_A * P - k * Bri_1 - Keff_rA * Bri_1
    d/dt(Bri_2) <- k * Bri_1 - k * Bri_2 - Keff_rA * Bri_2
    d/dt(Bri_3) <- k * Bri_2 - k * Bri_3 - Keff_rA * Bri_3
    d/dt(Bri_4) <- k * Bri_3 - k * Bri_4 - Keff_rA * Bri_4
    d/dt(Bri_5) <- k * Bri_4 - k * Bri_5 - Keff_rA * Bri_5
    d/dt(Bri_6) <- k * Bri_5 - k * Bri_6 - Keff_rA * Bri_6
    d/dt(Bri_7) <- k * Bri_6 - k * Bri_7 - Keff_rA * Bri_7
    d/dt(Bri_8) <- k * Bri_7 - k * Bri_8 - Keff_rA * Bri_8
    d/dt(Bri_9) <- k * Bri_8 - k * Bri_9 - Keff_rA * Bri_9
    d/dt(Bri_10) <- k * Bri_9 - k * Bri_10 - Keff_rA * Bri_10
    
    d/dt(Br_P) <- B_growth_rP * Br_P + Bu * qP - Keff_rP * Br_P - Br_P * qA
    d/dt(Br_AP) <- B_growth_rAP * Br_AP + qA * Br_P + qP * Br_A - Keff_rAP * Br_AP
    d/dt(P) <- b * k * (Bi_10 + Bri_10) - P_eff * (Bu + Br_A + Bi_tot + Br_i_tot) * P - Dec_P * P
    d/dt(A) <- -y * A
    
    # Initial conditions
    Bu(0) <- 10^7
    Bi_1(0) <- 0
    Bi_2(0) <- 0
    Bi_3(0) <- 0
    Bi_4(0) <- 0
    Bi_5(0) <- 0
    Bi_6(0) <- 0
    Bi_7(0) <- 0
    Bi_8(0) <- 0
    Bi_9(0) <- 0
    Bi_10(0) <- 0
    Br_A(0) <- 0
    Bri_1(0) <- 0
    Bri_2(0) <- 0
    Bri_3(0) <- 0
    Bri_4(0) <- 0
    Bri_5(0) <- 0
    Bri_6(0) <- 0
    Bri_7(0) <- 0
    Bri_8(0) <- 0
    Bri_9(0) <- 0
    Bri_10(0) <- 0
    Br_P(0) <- 0
    Br_AP(0) <- 0
    P(0) <- 0
    A(0) <- 0
  })
}
# same 4 post processing functions used in case of simultaneous administration
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

# function to solve the equations for a given timetable of treatment
# we get as output a dataframe where we have the time profile of the different species present
# we store also the outcome of treatment, minimum density reached, density of the resistant strain
# at the end of simulation..
dosing_function <- function(model,dosing_scheme,ADose,PDose){
  
  df_output_pp <- df_output
  df_output_pp <- post_process_Br_P(df_output_pp,30,1)
  df_output_pp <- post_process_Br_A(df_output_pp,30,1)
  df_output_pp <- post_process_Br_i_tot(df_output_pp,30,1)
  df_output_pp<- post_process_Br_AP(df_output_pp,30,1)
  df_output_pp  <- df_output_pp  %>%
    mutate(Bu = case_when(
      Bu < 1 ~ 0,  # If the sum of Bu, Bi_tot, Br is less than df_output <- model %>% rxSolve(dosing_scheme,method="lsoda")or equal to 1, set CFU to 0
      TRUE ~ Bu ))
  for (i in 1:nrow(df_output_pp)) {
    # Check if Bu, Br_A, and Br_P are 0, and Br_AP is below the threshold
    if (df_output_pp$Bu[i] == 0 && df_output_pp$Br_A[i] == 0 && df_output_pp$Br_P[i] == 0 && df_output_pp$Br_AP[i] < 1) {
      # Set Br_AP to 0 from this row till the end of the dataframe
      df_output_pp$Br_AP[i:nrow(df_output_pp)] <- 0
      break }
  }
  
  df_output_pp  <- df_output_pp  %>%
    mutate(CFU = case_when(
      Bu + Bi_tot+Br_A+Br_P+Br_AP +Br_i_tot< 1 ~ 0,  # If the sum of Bu, Bi_tot, Br is less than or equal to 10^-2, set CFU to 0
      TRUE ~ Bu + Bi_tot+Br_A+Br_P+Br_AP+ Br_i_tot ))
  df_output_pp$A.1 <- ADose
  df_output_pp$P.1 <- PDose
  df_output_pp$min_CFU <- min(df_output_pp$CFU)
  df_output_pp$CFU_at_24 <- df_output_pp$CFU[df_output_pp$time == 24][1]
  df_output_pp$CFU_at_72 <- df_output_pp$CFU[df_output_pp$time == 72][1]
  df_output_pp$Br_A_at_72 <- df_output_pp$Br_A[df_output_pp$time==72][1]
  df_output_pp$Br_P_at_72 <- df_output_pp$Br_P[df_output_pp$time==72][1]
  df_output_pp$Br_AP_at_72 <- df_output_pp$Br_AP[df_output_pp$time==72][1]
  df_output_pp$maxP <- max(df_output_pp$P)
  df_output_pp$time_to_min_CFU = df_output_pp$time[which.min(df_output_pp$CFU)]
  df_output_pp$BacteriaEradication <- ifelse(df_output_pp$CFU_at_72==0|df_output_pp$min_CFU==0,1,0)
  
  return(df_output_pp)
}
# plotting function for CFU and for all the species in the simulation
plot_CFU <- function(df,title){
  pl<-df %>% 
    ggplot(aes(x=time,y=CFU))+
    geom_line(linewidth=1)+
    #theme_bw(base_size=60)+
    geom_hline(yintercept = 1, linetype = "dashed", color = "green",size=2)+
    geom_hline(yintercept = 10^7, linetype = "dashed", color = "blue",size=2)+
    scale_y_log10(limits=c(10^-3,10^11)) +
    scale_x_continuous(breaks = c(0, 24, 48, 72)) +
    labs(x="time (h)",y="CFU")+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 25),  # Increase x-axis tick size
      axis.text.y = element_text(size = 25),  # Increase y-axis tick size
      axis.title.x = element_text(size = 40),  # Increase x-axis label size
      axis.title.y = element_text(size = 40)   # Increase y-axis label size
    )
  combined_plot <- pl
  jpeg(filename = paste0(output_dir,title), 
       width = 2000, height = 1500, res = 50)
  grid.draw(combined_plot)
  dev.off()
  return(pl)
}

plot_Spes <-function(df,title){ 
  pl<-df %>% pivot_longer(cols = c( Bu,Br_A,Br_P,Br_AP,P), 
                          names_to = "Variable", 
                          values_to = "Value") %>%  # Reshape the data to long format
    ggplot(aes(x = time, y = Value, color = Variable)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "green", size = 2) +
    geom_line(linewidth = 2) +
    theme_bw(base_size = 60) +
    scale_y_log10(limits=c(10^-3,10^11)) +
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
# we will have 6 different time delay between treatments:
# 1) 0H (simultaneous)  2) 0.2H (shorter than latency period)  3) 1H
# 4) 2H     5) 4H   6) 6H
# in one case we start with the phages, and the other we start with the antibiotic
# we choose 10 different phage antibiotic combinations:(antibiotic dose expressed as fraction of the MIC)
#1) A.1 = 1.5 & MOI = 10^-2   2)  A.1 = 1.5 & MOI = 10   3) A.1 = 1 & MOI = 10^-2
#4) A.1 = 1 & MOI = 10^-2  5) A.1 = 0.5 & MOI = 10^-2   6)  A.1 = 0.5 & MOI = 10
# these doses are the lowest effective doses in absence of CS,weak CS and strong CS in case of simultaneous administration
# 7) A.1 = 2 & MOI = 10^-4   8) A.1 = 3 & MOI = 10^-3
# 9) A.1 = 3 & MOI = 10^-3  8) A.1 = 4 & MOI = 10^-4
# these doses instead result in failure to eradicate bacteria when given in simultaneous

# we define now the treatment schemes here
t_list <- list()

ev_P_0_A_0_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 0
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
x_mod <- mod_noCS(100)
p <- x_mod %>% rxSolve(ev_P_0_A_0_P.1_5_A.1_1, method = "lsoda")
ev_P_0_A_0_P.1_5_A.1_1$timelag <- 0
ev_P_0_A_0_P.1_5_A.1_1$P.1 <- 10^5
ev_P_0_A_0_P.1_5_A.1_1$A.1 <- 1
ev_P_0_A_0_P.1_5_A.1_1$FirstTreat <- "Phages"
t_list[[length(t_list)+1]] <- ev_P_0_A_0_P.1_5_A.1_1

ev_P_0_A_halftau_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 0.2
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_5_A.1_1$timelag <- 0.2
ev_P_0_A_halftau_P.1_5_A.1_1$P.1 <- 10^5
ev_P_0_A_halftau_P.1_5_A.1_1$A.1 <- 1
ev_P_0_A_halftau_P.1_5_A.1_1$FirstTreat <- "Phages"
t_list[[length(t_list)+1]] <- ev_P_0_A_halftau_P.1_5_A.1_1

ev_P_0_A_1H_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 1
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_5_A.1_1$timelag <- 1
ev_P_0_A_1H_P.1_5_A.1_1$P.1 <- 10^5
ev_P_0_A_1H_P.1_5_A.1_1$A.1 <- 1
ev_P_0_A_1H_P.1_5_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_5_A.1_1

ev_P_0_A_2H_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 2
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_5_A.1_1$timelag <- 2
ev_P_0_A_2H_P.1_5_A.1_1$P.1 <- 10^5
ev_P_0_A_2H_P.1_5_A.1_1$A.1 <- 1
ev_P_0_A_2H_P.1_5_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <-ev_P_0_A_2H_P.1_5_A.1_1 

ev_P_0_A_4H_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 4
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_5_A.1_1$timelag <- 4
ev_P_0_A_4H_P.1_5_A.1_1$P.1 <- 10^5
ev_P_0_A_4H_P.1_5_A.1_1$A.1 <- 1
ev_P_0_A_4H_P.1_5_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <-ev_P_0_A_4H_P.1_5_A.1_1 

ev_P_0_A_6H_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 6
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_5_A.1_1$timelag <- 6
ev_P_0_A_6H_P.1_5_A.1_1$P.1 <- 10^5
ev_P_0_A_6H_P.1_5_A.1_1$A.1 <- 1
ev_P_0_A_6H_P.1_5_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <-ev_P_0_A_6H_P.1_5_A.1_1 

ev_P_0_A_0_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 0
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_0_P.1_8_A.1_1$timelag <- 0
ev_P_0_A_0_P.1_8_A.1_1$P.1 <- 10^8
ev_P_0_A_0_P.1_8_A.1_1$A.1 <- 1
ev_P_0_A_0_P.1_8_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_0_P.1_8_A.1_1

ev_P_0_A_halftau_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 0.2
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_8_A.1_1$timelag <- 0.2
ev_P_0_A_halftau_P.1_8_A.1_1$P.1 <- 10^8
ev_P_0_A_halftau_P.1_8_A.1_1$A.1 <- 1
ev_P_0_A_halftau_P.1_8_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_halftau_P.1_8_A.1_1

ev_P_0_A_1H_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 1
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_8_A.1_1$timelag <- 1
ev_P_0_A_1H_P.1_8_A.1_1$P.1 <- 10^8
ev_P_0_A_1H_P.1_8_A.1_1$A.1 <- 1
ev_P_0_A_1H_P.1_8_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_8_A.1_1

ev_P_0_A_2H_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 2
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_8_A.1_1$timelag <- 2
ev_P_0_A_2H_P.1_8_A.1_1$P.1 <- 10^8
ev_P_0_A_2H_P.1_8_A.1_1$A.1 <- 1
ev_P_0_A_2H_P.1_8_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_2H_P.1_8_A.1_1

ev_P_0_A_4H_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 4
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_8_A.1_1$timelag <- 4
ev_P_0_A_4H_P.1_8_A.1_1$P.1 <- 10^8
ev_P_0_A_4H_P.1_8_A.1_1$A.1 <- 1
ev_P_0_A_4H_P.1_8_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_4H_P.1_8_A.1_1

ev_P_0_A_6H_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 6
  et(amt = 1, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_8_A.1_1$timelag <- 6
ev_P_0_A_6H_P.1_8_A.1_1$P.1 <- 10^8
ev_P_0_A_6H_P.1_8_A.1_1$A.1 <- 1
ev_P_0_A_6H_P.1_8_A.1_1$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_6H_P.1_8_A.1_1

ev_P_0_A_0_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_0_P.1_5_A.1_0.5$timelag <- 0
ev_P_0_A_0_P.1_5_A.1_0.5$P.1 <- 10^5
ev_P_0_A_0_P.1_5_A.1_0.5$A.1 <- 0.5
ev_P_0_A_0_P.1_5_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_0_P.1_5_A.1_0.5

ev_P_0_A_halftau_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0.2
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_5_A.1_0.5$timelag <- 0.2
ev_P_0_A_halftau_P.1_5_A.1_0.5$P.1 <- 10^5
ev_P_0_A_halftau_P.1_5_A.1_0.5$A.1 <- 0.5
ev_P_0_A_halftau_P.1_5_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_halftau_P.1_5_A.1_0.5

ev_P_0_A_1H_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 1
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_5_A.1_0.5$timelag <- 1
ev_P_0_A_1H_P.1_5_A.1_0.5$P.1 <- 10^5
ev_P_0_A_1H_P.1_5_A.1_0.5$A.1 <- 0.5
ev_P_0_A_1H_P.1_5_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_5_A.1_0.5

ev_P_0_A_2H_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 2
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_5_A.1_0.5$timelag <- 2
ev_P_0_A_2H_P.1_5_A.1_0.5$P.1 <- 10^5
ev_P_0_A_2H_P.1_5_A.1_0.5$A.1 <- 0.5
ev_P_0_A_2H_P.1_5_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_2H_P.1_5_A.1_0.5

ev_P_0_A_4H_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 4
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_5_A.1_0.5$timelag <- 4
ev_P_0_A_4H_P.1_5_A.1_0.5$P.1 <- 10^5
ev_P_0_A_4H_P.1_5_A.1_0.5$A.1 <- 0.5
ev_P_0_A_4H_P.1_5_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_4H_P.1_5_A.1_0.5

ev_P_0_A_6H_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 6
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_5_A.1_0.5$timelag <- 6
ev_P_0_A_6H_P.1_5_A.1_0.5$P.1 <- 10^5
ev_P_0_A_6H_P.1_5_A.1_0.5$A.1 <- 0.5
ev_P_0_A_6H_P.1_5_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_6H_P.1_5_A.1_0.5

ev_P_0_A_0_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_0_P.1_8_A.1_0.5$timelag <- 0
ev_P_0_A_0_P.1_8_A.1_0.5$P.1 <- 10^8
ev_P_0_A_0_P.1_8_A.1_0.5$A.1 <- 0.5
ev_P_0_A_0_P.1_8_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_0_P.1_8_A.1_0.5

ev_P_0_A_halftau_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0.2
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_8_A.1_0.5$timelag <- 0.2
ev_P_0_A_halftau_P.1_8_A.1_0.5$P.1 <- 10^8
ev_P_0_A_halftau_P.1_8_A.1_0.5$A.1 <- 0.5
ev_P_0_A_halftau_P.1_8_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_halftau_P.1_8_A.1_0.5

ev_P_0_A_1H_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 1
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_8_A.1_0.5$timelag <- 1
ev_P_0_A_1H_P.1_8_A.1_0.5$P.1 <- 10^8
ev_P_0_A_1H_P.1_8_A.1_0.5$A.1 <- 0.5
ev_P_0_A_1H_P.1_8_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_8_A.1_0.5

ev_P_0_A_2H_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 2
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_8_A.1_0.5$timelag <- 2
ev_P_0_A_2H_P.1_8_A.1_0.5$P.1 <- 10^8
ev_P_0_A_2H_P.1_8_A.1_0.5$A.1 <- 0.5
ev_P_0_A_2H_P.1_8_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_2H_P.1_8_A.1_0.5

ev_P_0_A_4H_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 4
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_8_A.1_0.5$timelag <- 4
ev_P_0_A_4H_P.1_8_A.1_0.5$P.1 <- 10^8
ev_P_0_A_4H_P.1_8_A.1_0.5$A.1 <- 0.5
ev_P_0_A_4H_P.1_8_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_4H_P.1_8_A.1_0.5

ev_P_0_A_6H_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 6
  et(amt = 0.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_8_A.1_0.5$timelag <- 6
ev_P_0_A_6H_P.1_8_A.1_0.5$P.1 <- 10^8
ev_P_0_A_6H_P.1_8_A.1_0.5$A.1 <- 0.5
ev_P_0_A_6H_P.1_8_A.1_0.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_6H_P.1_8_A.1_0.5

ev_P_0_A_0_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 0
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_0_P.1_5_A.1_1.5$timelag <- 0
ev_P_0_A_0_P.1_5_A.1_1.5$P.1 <- 10^5
ev_P_0_A_0_P.1_5_A.1_1.5$A.1 <- 1.5
ev_P_0_A_0_P.1_5_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_0_P.1_5_A.1_1.5

ev_P_0_A_halftau_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 0.2
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_5_A.1_1.5$timelag <- 0.2
ev_P_0_A_halftau_P.1_5_A.1_1.5$P.1 <- 10^5
ev_P_0_A_halftau_P.1_5_A.1_1.5$A.1 <- 1.5
ev_P_0_A_halftau_P.1_5_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_halftau_P.1_5_A.1_1.5

ev_P_0_A_1H_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 1
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_5_A.1_1.5$timelag <- 1
ev_P_0_A_1H_P.1_5_A.1_1.5$P.1 <- 10^5
ev_P_0_A_1H_P.1_5_A.1_1.5$A.1 <- 1.5
ev_P_0_A_1H_P.1_5_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_5_A.1_1.5

ev_P_0_A_2H_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 2
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_5_A.1_1.5$timelag <- 2
ev_P_0_A_2H_P.1_5_A.1_1.5$P.1 <- 10^5
ev_P_0_A_2H_P.1_5_A.1_1.5$A.1 <- 1.5
ev_P_0_A_2H_P.1_5_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_2H_P.1_5_A.1_1.5

ev_P_0_A_4H_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 4
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_5_A.1_1.5$timelag <- 4
ev_P_0_A_4H_P.1_5_A.1_1.5$P.1 <- 10^5
ev_P_0_A_4H_P.1_5_A.1_1.5$A.1 <- 1.5
ev_P_0_A_4H_P.1_5_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_4H_P.1_5_A.1_1.5

ev_P_0_A_6H_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 6
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_5_A.1_1.5$timelag <- 6
ev_P_0_A_6H_P.1_5_A.1_1.5$P.1 <- 10^5
ev_P_0_A_6H_P.1_5_A.1_1.5$A.1 <- 1.5
ev_P_0_A_6H_P.1_5_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_6H_P.1_5_A.1_1.5

ev_P_0_A_0_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 0
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_0_P.1_8_A.1_1.5$timelag <- 0
ev_P_0_A_0_P.1_8_A.1_1.5$P.1 <- 10^8
ev_P_0_A_0_P.1_8_A.1_1.5$A.1 <- 1.5
ev_P_0_A_0_P.1_8_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_0_P.1_8_A.1_1.5

ev_P_0_A_halftau_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 0.2
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_8_A.1_1.5$timelag <- 0.2
ev_P_0_A_halftau_P.1_8_A.1_1.5$P.1 <- 10^8
ev_P_0_A_halftau_P.1_8_A.1_1.5$A.1 <- 1.5
ev_P_0_A_halftau_P.1_8_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_halftau_P.1_8_A.1_1.5

ev_P_0_A_1H_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 1
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_8_A.1_1.5$timelag <- 1
ev_P_0_A_1H_P.1_8_A.1_1.5$P.1 <- 10^8
ev_P_0_A_1H_P.1_8_A.1_1.5$A.1 <- 1.5
ev_P_0_A_1H_P.1_8_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_8_A.1_1.5

ev_P_0_A_2H_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 2
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_8_A.1_1.5$timelag <- 2
ev_P_0_A_2H_P.1_8_A.1_1.5$P.1 <- 10^8
ev_P_0_A_2H_P.1_8_A.1_1.5$A.1 <- 1.5
ev_P_0_A_2H_P.1_8_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_2H_P.1_8_A.1_1.5

ev_P_0_A_4H_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 4
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_8_A.1_1.5$timelag <- 4
ev_P_0_A_4H_P.1_8_A.1_1.5$P.1 <- 10^8
ev_P_0_A_4H_P.1_8_A.1_1.5$A.1 <- 1.5
ev_P_0_A_4H_P.1_8_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_4H_P.1_8_A.1_1.5

ev_P_0_A_6H_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 6
  et(amt = 1.5, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_8_A.1_1.5$timelag <- 6
ev_P_0_A_6H_P.1_8_A.1_1.5$P.1 <- 10^8
ev_P_0_A_6H_P.1_8_A.1_1.5$A.1 <- 1.5
ev_P_0_A_6H_P.1_8_A.1_1.5$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_6H_P.1_8_A.1_1.5
##### now we invert the order of treatment
ev_A_0_P_halftau_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 0.2
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_5_A.1_1$timelag <- 0.2
ev_A_0_P_halftau_P.1_5_A.1_1$P.1 <- 10^5
ev_A_0_P_halftau_P.1_5_A.1_1$A.1 <- 1
ev_A_0_P_halftau_P.1_5_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_5_A.1_1

ev_A_0_P_1H_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 1
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_5_A.1_1$timelag <- 1
ev_A_0_P_1H_P.1_5_A.1_1$P.1 <- 10^5
ev_A_0_P_1H_P.1_5_A.1_1$A.1 <- 1
ev_A_0_P_1H_P.1_5_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_5_A.1_1

ev_A_0_P_2H_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 2
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_5_A.1_1$timelag <- 2
ev_A_0_P_2H_P.1_5_A.1_1$P.1 <- 10^5
ev_A_0_P_2H_P.1_5_A.1_1$A.1 <- 1
ev_A_0_P_2H_P.1_5_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_5_A.1_1

ev_A_0_P_4H_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 4
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_5_A.1_1$timelag <- 4
ev_A_0_P_4H_P.1_5_A.1_1$P.1 <- 10^5
ev_A_0_P_4H_P.1_5_A.1_1$A.1 <- 1
ev_A_0_P_4H_P.1_5_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_5_A.1_1

ev_A_0_P_6H_P.1_5_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 6
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_5_A.1_1$timelag <- 6
ev_A_0_P_6H_P.1_5_A.1_1$P.1 <- 10^5
ev_A_0_P_6H_P.1_5_A.1_1$A.1 <- 1
ev_A_0_P_6H_P.1_5_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_5_A.1_1

ev_A_0_P_halftau_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 0.2
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_8_A.1_1$timelag <- 0.2
ev_A_0_P_halftau_P.1_8_A.1_1$P.1 <- 10^8
ev_A_0_P_halftau_P.1_8_A.1_1$A.1 <- 1
ev_A_0_P_halftau_P.1_8_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_8_A.1_1

ev_A_0_P_1H_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 1
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_8_A.1_1$timelag <- 1
ev_A_0_P_1H_P.1_8_A.1_1$P.1 <- 10^8
ev_A_0_P_1H_P.1_8_A.1_1$A.1 <- 1
ev_A_0_P_1H_P.1_8_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_8_A.1_1

ev_A_0_P_2H_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 2
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_8_A.1_1$timelag <- 2
ev_A_0_P_2H_P.1_8_A.1_1$P.1 <- 10^8
ev_A_0_P_2H_P.1_8_A.1_1$A.1 <- 1
ev_A_0_P_2H_P.1_8_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_8_A.1_1

ev_A_0_P_4H_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 4
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time =4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_8_A.1_1$timelag <- 4
ev_A_0_P_4H_P.1_8_A.1_1$P.1 <- 10^8
ev_A_0_P_4H_P.1_8_A.1_1$A.1 <- 1
ev_A_0_P_4H_P.1_8_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_8_A.1_1

ev_A_0_P_6H_P.1_8_A.1_1<- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 6
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_8_A.1_1$timelag <- 6
ev_A_0_P_6H_P.1_8_A.1_1$P.1 <- 10^8
ev_A_0_P_6H_P.1_8_A.1_1$A.1 <- 1
ev_A_0_P_6H_P.1_8_A.1_1$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_8_A.1_1

ev_A_0_P_halftau_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 0.2
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_5_A.1_0.5$timelag <- 0.2
ev_A_0_P_halftau_P.1_5_A.1_0.5$P.1 <- 10^5
ev_A_0_P_halftau_P.1_5_A.1_0.5$A.1 <- 0.5
ev_A_0_P_halftau_P.1_5_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_5_A.1_0.5

ev_A_0_P_1H_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 1
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_5_A.1_0.5$timelag <- 1
ev_A_0_P_1H_P.1_5_A.1_0.5$P.1 <- 10^5
ev_A_0_P_1H_P.1_5_A.1_0.5$A.1 <- 0.5
ev_A_0_P_1H_P.1_5_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_5_A.1_0.5

ev_A_0_P_2H_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 2
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_5_A.1_0.5$timelag <- 2
ev_A_0_P_2H_P.1_5_A.1_0.5$P.1 <- 10^5
ev_A_0_P_2H_P.1_5_A.1_0.5$A.1 <- 0.5
ev_A_0_P_2H_P.1_5_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_5_A.1_0.5

ev_A_0_P_4H_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 4
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_5_A.1_0.5$timelag <- 4
ev_A_0_P_4H_P.1_5_A.1_0.5$P.1 <- 10^5
ev_A_0_P_4H_P.1_5_A.1_0.5$A.1 <- 0.5
ev_A_0_P_4H_P.1_5_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_5_A.1_0.5

ev_A_0_P_6H_P.1_5_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 6
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_5_A.1_0.5$timelag <- 6
ev_A_0_P_6H_P.1_5_A.1_0.5$P.1 <- 10^5
ev_A_0_P_6H_P.1_5_A.1_0.5$A.1 <- 0.5
ev_A_0_P_6H_P.1_5_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_5_A.1_0.5

ev_A_0_P_halftau_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 0.2
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_8_A.1_0.5$timelag <- 0.2
ev_A_0_P_halftau_P.1_8_A.1_0.5$P.1 <- 10^8
ev_A_0_P_halftau_P.1_8_A.1_0.5$A.1 <- 0.5
ev_A_0_P_halftau_P.1_8_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_8_A.1_0.5

ev_A_0_P_1H_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 1
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_8_A.1_0.5$timelag <- 1
ev_A_0_P_1H_P.1_8_A.1_0.5$P.1 <- 10^8
ev_A_0_P_1H_P.1_8_A.1_0.5$A.1 <- 0.5
ev_A_0_P_1H_P.1_8_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_8_A.1_0.5

ev_A_0_P_2H_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 2
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_8_A.1_0.5$timelag <- 2
ev_A_0_P_2H_P.1_8_A.1_0.5$P.1 <- 10^8
ev_A_0_P_2H_P.1_8_A.1_0.5$A.1 <- 0.5
ev_A_0_P_2H_P.1_8_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_8_A.1_0.5

ev_A_0_P_4H_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 4
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_8_A.1_0.5$timelag <- 4
ev_A_0_P_4H_P.1_8_A.1_0.5$P.1 <- 10^8
ev_A_0_P_4H_P.1_8_A.1_0.5$A.1 <- 0.5
ev_A_0_P_4H_P.1_8_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_8_A.1_0.5

ev_A_0_P_6H_P.1_8_A.1_0.5<- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 6
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_8_A.1_0.5$timelag <- 6
ev_A_0_P_6H_P.1_8_A.1_0.5$P.1 <- 10^8
ev_A_0_P_6H_P.1_8_A.1_0.5$A.1 <- 0.5
ev_A_0_P_6H_P.1_8_A.1_0.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_8_A.1_0.5

ev_A_0_P_halftau_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 0.2
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_5_A.1_1.5$timelag <- 0.2
ev_A_0_P_halftau_P.1_5_A.1_1.5$P.1 <- 10^5
ev_A_0_P_halftau_P.1_5_A.1_1.5$A.1 <- 1.5
ev_A_0_P_halftau_P.1_5_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_5_A.1_1.5

ev_A_0_P_1H_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 1
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_5_A.1_1.5$timelag <- 1
ev_A_0_P_1H_P.1_5_A.1_1.5$P.1 <- 10^5
ev_A_0_P_1H_P.1_5_A.1_1.5$A.1 <- 1.5
ev_A_0_P_1H_P.1_5_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_5_A.1_1.5

ev_A_0_P_2H_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 2
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_5_A.1_1.5$timelag <- 2
ev_A_0_P_2H_P.1_5_A.1_1.5$P.1 <- 10^5
ev_A_0_P_2H_P.1_5_A.1_1.5$A.1 <- 1.5
ev_A_0_P_2H_P.1_5_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_5_A.1_1.5

ev_A_0_P_4H_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 4
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_5_A.1_1.5$timelag <- 4
ev_A_0_P_4H_P.1_5_A.1_1.5$P.1 <- 10^5
ev_A_0_P_4H_P.1_5_A.1_1.5$A.1 <- 1.5
ev_A_0_P_4H_P.1_5_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_5_A.1_1.5

ev_A_0_P_6H_P.1_5_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^5 of P at time 6
  et(amt = 10^5, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_5_A.1_1.5$timelag <- 6
ev_A_0_P_6H_P.1_5_A.1_1.5$P.1 <- 10^5
ev_A_0_P_6H_P.1_5_A.1_1.5$A.1 <- 1.5
ev_A_0_P_6H_P.1_5_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_5_A.1_1.5

ev_A_0_P_halftau_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 0.2
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_8_A.1_1.5$timelag <- 0.2
ev_A_0_P_halftau_P.1_8_A.1_1.5$P.1 <- 10^8
ev_A_0_P_halftau_P.1_8_A.1_1.5$A.1 <- 1.5
ev_A_0_P_halftau_P.1_8_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_8_A.1_1.5

ev_A_0_P_1H_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 1
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_8_A.1_1.5$timelag <- 1
ev_A_0_P_1H_P.1_8_A.1_1.5$P.1 <- 10^8
ev_A_0_P_1H_P.1_8_A.1_1.5$A.1 <- 1.5
ev_A_0_P_1H_P.1_8_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_8_A.1_1.5

ev_A_0_P_2H_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 2
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_8_A.1_1.5$timelag <- 2
ev_A_0_P_2H_P.1_8_A.1_1.5$P.1 <- 10^8
ev_A_0_P_2H_P.1_8_A.1_1.5$A.1 <- 1.5
ev_A_0_P_2H_P.1_8_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_8_A.1_1.5

ev_A_0_P_4H_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 4
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_8_A.1_1.5$timelag <- 4
ev_A_0_P_4H_P.1_8_A.1_1.5$P.1 <- 10^8
ev_A_0_P_4H_P.1_8_A.1_1.5$A.1 <- 1.5
ev_A_0_P_4H_P.1_8_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_8_A.1_1.5

ev_A_0_P_6H_P.1_8_A.1_1.5<- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^8 of P at time 6
  et(amt = 10^8, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_8_A.1_1.5$timelag <- 6
ev_A_0_P_6H_P.1_8_A.1_1.5$P.1 <- 10^8
ev_A_0_P_6H_P.1_8_A.1_1.5$A.1 <- 1.5
ev_A_0_P_6H_P.1_8_A.1_1.5$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_8_A.1_1.5

ev_P_0_A_halftau_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 0.2
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_3_A.1_3$timelag <- 0.2
ev_P_0_A_halftau_P.1_3_A.1_3$P.1 <- 10^3
ev_P_0_A_halftau_P.1_3_A.1_3$A.1 <- 3
ev_P_0_A_halftau_P.1_3_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_halftau_P.1_3_A.1_3

ev_P_0_A_1H_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 1
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_3_A.1_3$timelag <- 1
ev_P_0_A_1H_P.1_3_A.1_3$P.1 <- 10^3
ev_P_0_A_1H_P.1_3_A.1_3$A.1 <- 3
ev_P_0_A_1H_P.1_3_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_3_A.1_3

ev_P_0_A_2H_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 2
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_3_A.1_3$timelag <- 2
ev_P_0_A_2H_P.1_3_A.1_3$P.1 <- 10^3
ev_P_0_A_2H_P.1_3_A.1_3$A.1 <- 3
ev_P_0_A_2H_P.1_3_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_2H_P.1_3_A.1_3

ev_P_0_A_4H_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 4
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_3_A.1_3$timelag <- 4
ev_P_0_A_4H_P.1_3_A.1_3$P.1 <- 10^3
ev_P_0_A_4H_P.1_3_A.1_3$A.1 <- 3
ev_P_0_A_4H_P.1_3_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_4H_P.1_3_A.1_3

ev_P_0_A_6H_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 6
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_3_A.1_3$timelag <- 6
ev_P_0_A_6H_P.1_3_A.1_3$P.1 <- 10^3
ev_P_0_A_6H_P.1_3_A.1_3$A.1 <- 3
ev_P_0_A_6H_P.1_3_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_6H_P.1_3_A.1_3
# A.1 = 3 & P.1 = 10^4
ev_P_0_A_halftau_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 10^4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 0.2
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_4_A.1_3$timelag <- 0.2
ev_P_0_A_halftau_P.1_4_A.1_3$P.1 <- 10^4
ev_P_0_A_halftau_P.1_4_A.1_3$A.1 <- 3
ev_P_0_A_halftau_P.1_4_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_halftau_P.1_4_A.1_3

ev_P_0_A_1H_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 10^4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 1
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_4_A.1_3$timelag <- 1
ev_P_0_A_1H_P.1_4_A.1_3$P.1 <- 10^4
ev_P_0_A_1H_P.1_4_A.1_3$A.1 <- 3
ev_P_0_A_1H_P.1_4_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_4_A.1_3

ev_P_0_A_2H_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 10^4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 2
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_4_A.1_3$timelag <- 2
ev_P_0_A_2H_P.1_4_A.1_3$P.1 <- 10^4
ev_P_0_A_2H_P.1_4_A.1_3$A.1 <- 3
ev_P_0_A_2H_P.1_4_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_2H_P.1_4_A.1_3

ev_P_0_A_4H_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 10^4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 4
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_4_A.1_3$timelag <- 4
ev_P_0_A_4H_P.1_4_A.1_3$P.1 <- 10^4
ev_P_0_A_4H_P.1_4_A.1_3$A.1 <- 3
ev_P_0_A_4H_P.1_4_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_4H_P.1_4_A.1_3

ev_P_0_A_6H_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 10^4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 6
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_4_A.1_3$timelag <- 6
ev_P_0_A_6H_P.1_4_A.1_3$P.1 <- 10^4
ev_P_0_A_6H_P.1_4_A.1_3$A.1 <- 3
ev_P_0_A_6H_P.1_4_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_6H_P.1_4_A.1_3
# first A, A = 3 P.1 = 10^3
ev_A_0_P_halftau_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 0.2
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_3_A.1_3$timelag <- 0.2
ev_A_0_P_halftau_P.1_3_A.1_3$P.1 <- 10^3
ev_A_0_P_halftau_P.1_3_A.1_3$A.1 <- 3
ev_A_0_P_halftau_P.1_3_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_3_A.1_3

ev_A_0_P_1H_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 1
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_3_A.1_3$timelag <- 1
ev_A_0_P_1H_P.1_3_A.1_3$P.1 <- 10^3
ev_A_0_P_1H_P.1_3_A.1_3$A.1 <- 3
ev_A_0_P_1H_P.1_3_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_3_A.1_3

ev_A_0_P_2H_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 2
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_3_A.1_3$timelag <- 2
ev_A_0_P_2H_P.1_3_A.1_3$P.1 <- 10^3
ev_A_0_P_2H_P.1_3_A.1_3$A.1 <- 3
ev_A_0_P_2H_P.1_3_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_3_A.1_3

ev_A_0_P_4H_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 4
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_3_A.1_3$timelag <- 4
ev_A_0_P_4H_P.1_3_A.1_3$P.1 <- 10^3
ev_A_0_P_4H_P.1_3_A.1_3$A.1 <- 3
ev_A_0_P_4H_P.1_3_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_3_A.1_3

ev_A_0_P_6H_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 6
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_3_A.1_3$timelag <- 6
ev_A_0_P_6H_P.1_3_A.1_3$P.1 <- 10^3
ev_A_0_P_6H_P.1_3_A.1_3$A.1 <- 3
ev_A_0_P_6H_P.1_3_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_3_A.1_3
# first A, A = 3 P.1 = 10^4
ev_A_0_P_halftau_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^4 of P at time 0.2
  et(amt = 10^4, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_4_A.1_3$timelag <- 0.2
ev_A_0_P_halftau_P.1_4_A.1_3$P.1 <- 10^4
ev_A_0_P_halftau_P.1_4_A.1_3$A.1 <- 3
ev_A_0_P_halftau_P.1_4_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_4_A.1_3

ev_A_0_P_1H_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^4 of P at time 1
  et(amt = 10^4, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_4_A.1_3$timelag <- 1
ev_A_0_P_1H_P.1_4_A.1_3$P.1 <- 10^4
ev_A_0_P_1H_P.1_4_A.1_3$A.1 <- 3
ev_A_0_P_1H_P.1_4_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_4_A.1_3

ev_A_0_P_2H_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^4 of P at time 2
  et(amt = 10^4, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_4_A.1_3$timelag <- 2
ev_A_0_P_2H_P.1_4_A.1_3$P.1 <- 10^4
ev_A_0_P_2H_P.1_4_A.1_3$A.1 <- 3
ev_A_0_P_2H_P.1_4_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_4_A.1_3

ev_A_0_P_4H_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^4 of P at time 4
  et(amt = 10^4, addl = 0, ii = 24, cmt = "P", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_4_A.1_3$timelag <- 4
ev_A_0_P_4H_P.1_4_A.1_3$P.1 <- 10^4
ev_A_0_P_4H_P.1_4_A.1_3$A.1 <- 3
ev_A_0_P_4H_P.1_4_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_4_A.1_3

ev_A_0_P_6H_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^4 of P at time 6
  et(amt = 10^4, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_4_A.1_3$timelag <- 6
ev_A_0_P_6H_P.1_4_A.1_3$P.1 <- 10^4
ev_A_0_P_6H_P.1_4_A.1_3$A.1 <- 3
ev_A_0_P_6H_P.1_4_A.1_3$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_4_A.1_3
# A.1 = 2, P.1 = 10^3
ev_P_0_A_halftau_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 2 of A at time 0.2
  et(amt = 2, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_3_A.1_2$timelag <- 0.2
ev_P_0_A_halftau_P.1_3_A.1_2$P.1 <- 10^3
ev_P_0_A_halftau_P.1_3_A.1_2$A.1 <- 2
ev_P_0_A_halftau_P.1_3_A.1_2$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_halftau_P.1_3_A.1_2

ev_P_0_A_1H_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 2 of A at time 1
  et(amt = 2, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_3_A.1_2$timelag <- 1
ev_P_0_A_1H_P.1_3_A.1_2$P.1 <- 10^3
ev_P_0_A_1H_P.1_3_A.1_2$A.1 <- 2
ev_P_0_A_1H_P.1_3_A.1_2$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_3_A.1_2

ev_P_0_A_2H_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 2 of A at time 2
  et(amt = 2, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_3_A.1_2$timelag <- 2
ev_P_0_A_2H_P.1_3_A.1_2$P.1 <- 10^3
ev_P_0_A_2H_P.1_3_A.1_2$A.1 <- 2
ev_P_0_A_2H_P.1_3_A.1_2$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_2H_P.1_3_A.1_2

ev_P_0_A_4H_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 2 of A at time 4
  et(amt = 2, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_3_A.1_2$timelag <- 4
ev_P_0_A_4H_P.1_3_A.1_2$P.1 <- 10^3
ev_P_0_A_4H_P.1_3_A.1_2$A.1 <- 2
ev_P_0_A_4H_P.1_3_A.1_2$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_4H_P.1_3_A.1_2

ev_P_0_A_6H_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 2 of A at time 6
  et(amt = 2, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_3_A.1_2$timelag <- 6
ev_P_0_A_6H_P.1_3_A.1_2$P.1 <- 10^3
ev_P_0_A_6H_P.1_3_A.1_2$A.1 <- 2
ev_P_0_A_6H_P.1_3_A.1_2$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_6H_P.1_3_A.1_2
# A.1 = 4, P.1 = 10^3
ev_P_0_A_halftau_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 4 of A at time 0.2
  et(amt = 4, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_halftau_P.1_3_A.1_4$timelag <- 0.2
ev_P_0_A_halftau_P.1_3_A.1_4$P.1 <- 10^3
ev_P_0_A_halftau_P.1_3_A.1_4$A.1 <- 4
ev_P_0_A_halftau_P.1_3_A.1_4$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_halftau_P.1_3_A.1_4

ev_P_0_A_1H_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 4 of A at time 1
  et(amt = 4, addl = 0, ii = 24, cmt = "A", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_1H_P.1_3_A.1_4$timelag <- 1
ev_P_0_A_1H_P.1_3_A.1_4$P.1 <- 10^3
ev_P_0_A_1H_P.1_3_A.1_4$A.1 <- 4
ev_P_0_A_1H_P.1_3_A.1_4$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_1H_P.1_3_A.1_4

ev_P_0_A_2H_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 4 of A at time 2
  et(amt = 4, addl = 0, ii = 24, cmt = "A", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_2H_P.1_3_A.1_4$timelag <- 2
ev_P_0_A_2H_P.1_3_A.1_4$P.1 <- 10^3
ev_P_0_A_2H_P.1_3_A.1_4$A.1 <- 4
ev_P_0_A_2H_P.1_3_A.1_4$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_2H_P.1_3_A.1_4

ev_P_0_A_4H_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 4 of A at time 4
  et(amt = 4, addl = 0, ii = 24, cmt = "A", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_4H_P.1_3_A.1_4$timelag <- 4
ev_P_0_A_4H_P.1_3_A.1_4$P.1 <- 10^3
ev_P_0_A_4H_P.1_3_A.1_4$A.1 <- 4
ev_P_0_A_4H_P.1_3_A.1_4$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_4H_P.1_3_A.1_4

ev_P_0_A_6H_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 4 of A at time 6
  et(amt = 4, addl = 0, ii = 24, cmt = "A", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_0_A_6H_P.1_3_A.1_4$timelag <- 6
ev_P_0_A_6H_P.1_3_A.1_4$P.1 <- 10^3
ev_P_0_A_6H_P.1_3_A.1_4$A.1 <- 4
ev_P_0_A_6H_P.1_3_A.1_4$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_P_0_A_6H_P.1_3_A.1_4
# A.1 = 2 P.1 = 10^3, first antibiotic
ev_A_0_P_halftau_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 2,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 0.2
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_3_A.1_2$timelag <- 0.2
ev_A_0_P_halftau_P.1_3_A.1_2$P.1 <- 10^3
ev_A_0_P_halftau_P.1_3_A.1_2$A.1 <- 2
ev_A_0_P_halftau_P.1_3_A.1_2$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_3_A.1_2

ev_A_0_P_1H_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 2,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 1
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_3_A.1_2$timelag <- 1
ev_A_0_P_1H_P.1_3_A.1_2$P.1 <- 10^3
ev_A_0_P_1H_P.1_3_A.1_2$A.1 <- 2
ev_A_0_P_1H_P.1_3_A.1_2$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_3_A.1_2

ev_A_0_P_2H_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 2,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 2
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_3_A.1_2$timelag <- 2
ev_A_0_P_2H_P.1_3_A.1_2$P.1 <- 10^3
ev_A_0_P_2H_P.1_3_A.1_2$A.1 <- 2
ev_A_0_P_2H_P.1_3_A.1_2$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_3_A.1_2

ev_A_0_P_4H_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 2,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 4
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_3_A.1_2$timelag <- 4
ev_A_0_P_4H_P.1_3_A.1_2$P.1 <- 10^3
ev_A_0_P_4H_P.1_3_A.1_2$A.1 <- 2
ev_A_0_P_4H_P.1_3_A.1_2$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_3_A.1_2

ev_A_0_P_6H_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 2,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 6
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_3_A.1_2$timelag <- 6
ev_A_0_P_6H_P.1_3_A.1_2$P.1 <- 10^3
ev_A_0_P_6H_P.1_3_A.1_2$A.1 <- 2
ev_A_0_P_6H_P.1_3_A.1_2$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_3_A.1_2
# A.1 = 4, P.1 = 10^3, first A
ev_A_0_P_halftau_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 0.2
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 0.2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_halftau_P.1_3_A.1_4$timelag <- 0.2
ev_A_0_P_halftau_P.1_3_A.1_4$P.1 <- 10^3
ev_A_0_P_halftau_P.1_3_A.1_4$A.1 <- 4
ev_A_0_P_halftau_P.1_3_A.1_4$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_halftau_P.1_3_A.1_4

ev_A_0_P_1H_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 1
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 1) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_1H_P.1_3_A.1_4$timelag <- 1
ev_A_0_P_1H_P.1_3_A.1_4$P.1 <- 10^3
ev_A_0_P_1H_P.1_3_A.1_4$A.1 <- 4
ev_A_0_P_1H_P.1_3_A.1_4$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_1H_P.1_3_A.1_4

ev_A_0_P_2H_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 2
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 2) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_2H_P.1_3_A.1_4$timelag <- 2
ev_A_0_P_2H_P.1_3_A.1_4$P.1 <- 10^3
ev_A_0_P_2H_P.1_3_A.1_4$A.1 <- 4
ev_A_0_P_2H_P.1_3_A.1_4$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_2H_P.1_3_A.1_4

ev_A_0_P_4H_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 4
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 4) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_4H_P.1_3_A.1_4$timelag <- 4
ev_A_0_P_4H_P.1_3_A.1_4$P.1 <- 10^3
ev_A_0_P_4H_P.1_3_A.1_4$A.1 <- 4
ev_A_0_P_4H_P.1_3_A.1_4$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_4H_P.1_3_A.1_4

ev_A_0_P_6H_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 10^3 of P at time 6
  et(amt = 10^3, addl = 0, ii = 24, cmt = "P", evid = 1, time = 6) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_A_0_P_6H_P.1_3_A.1_4$timelag <- 6
ev_A_0_P_6H_P.1_3_A.1_4$P.1 <- 10^3
ev_A_0_P_6H_P.1_3_A.1_4$A.1 <- 4
ev_A_0_P_6H_P.1_3_A.1_4$FirstTreat <- "Antibiotic"
t_list[[lengtht_list+1]] <- ev_A_0_P_6H_P.1_3_A.1_4



ev_sim_P.1_3_A.1_3<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 0
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_sim_P.1_3_A.1_3$timelag <- 0
ev_sim_P.1_3_A.1_3$P.1 <- 10^3
ev_sim_P.1_3_A.1_3$A.1 <- 3
ev_sim_P.1_3_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_sim_P.1_3_A.1_3

ev_sim_P.1_4_A.1_3<- et( # First dose for P at time 0
  amt  = 10^4,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 3 of A at time 0
  et(amt = 3, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_sim_P.1_4_A.1_3$timelag <- 0
ev_sim_P.1_4_A.1_3$P.1 <- 10^4
ev_sim_P.1_4_A.1_3$A.1 <- 3
ev_sim_P.1_4_A.1_3$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_sim_P.1_4_A.1_3

ev_sim_P.1_3_A.1_2<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 2 of A at time 0
  et(amt = 2, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_sim_P.1_3_A.1_2$timelag <- 0
ev_sim_P.1_3_A.1_2$P.1 <- 10^3
ev_sim_P.1_3_A.1_2$A.1 <- 2
ev_sim_P.1_3_A.1_2$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_sim_P.1_3_A.1_2

ev_sim_P.1_3_A.1_4<- et( # First dose for P at time 0
  amt  = 10^3,   # mg
  addl = 0,      # No additional doses
  ii   = 24,     # Not necessary in this case since there's only one dose
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 4 of A at time 0
  et(amt = 4, addl = 0, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_sim_P.1_3_A.1_4$timelag <- 0
ev_sim_P.1_3_A.1_4$P.1 <- 10^3
ev_sim_P.1_3_A.1_4$A.1 <- 4
ev_sim_P.1_3_A.1_4$FirstTreat <- "Phages"
t_list[[lengtht_list+1]] <- ev_sim_P.1_3_A.1_4
## for each treatment combination i call the dosing function on it, plot the cfu and all bacterial species
# for each unique combination of antibiotic and phage dose
# to account for CS, we just need to call the dosing function with the weak CS or strong CS model
d_list <- rep(0,14)
for (t in t_list){
  A.1<- t$A.1[1]
  P.1<- t$P.1[1]
  if (P.1 == 10^3) {
    P.1_s <- "10^3"}
  else if (P.1==10^4) {
    P.1_s <- "10^4"}
  else if (P.1==10^5){
    P.1_s <- "10^5"
  } else {
    P.1_s <- "10^8"
  }
  time_l <- t$timelag[1]
  #print(time_l)
  first_t <- t$FirstTreat[1]
  if (first_t == "Phages") {
    simu_df <- dosing_function(mod_noCS,t,A.1,P.1)
    summary_df <- simu_df[1,40:50]
    summary_df$timelag <- time_l
    time_l <- paste0(as.character(time_l),"H")
    summary_df$First_Treat <-first_t
    summary_df$CE <- "NO CS"
    d_list <- rbind(d_list,summary_df)
    CFU_0 <- paste("/FirstP_CFU-time_A_",as.character(A.1))
    CFU_1 <- paste("_P_",P.1_s)
    CFU_2 <- paste("_","NO_CS_")
    CFU_3 <- paste(time_l,".jpg")
    cfu_title <- paste0(CFU_0,CFU_1,CFU_2,CFU_3)
    all_0 <- paste("/FirstP_AllSpes-time_A_",as.character(A.1))
    all_1 <- paste("_P_",P.1_s)
    all_2 <- paste("_","NO_CS_")
    all_3 <- paste(time_l,".jpg")
    all_title <- paste0(all_0,all_1,all_2,all_3)
    pl <- plot_CFU(simu_df,cfu_title)
    pl <- plot_Spes(simu_df,all_title)
  }
  else {
    simu_df <- dosing_function(mod_noCS,t,A.1,P.1)
    summary_df <- simu_df[1,40:50]
    summary_df$timelag <- time_l
    time_l <- paste0(as.character(time_l),"H")
    summary_df$First_Treat <-first_t
    summary_df$CE <- "NO CS"
    d_list <- rbind(d_list,summary_df)
    CFU_0 <- paste("/FirstAnt_CFU-time_A_",as.character(A.1))
    CFU_1 <- paste("_P_",P.1_s)
    CFU_2 <- paste("_","NO_CS_")
    CFU_3 <- paste(time_l,".jpg")
    cfu_title <- paste0(CFU_0,CFU_1,CFU_2,CFU_3)
    all_0 <- paste("/FirstAnt_AllSpes-time_A_",as.character(A.1))
    all_1 <- paste("_P_",P.1_s)
    all_2 <- paste("_","NO_CS_")
    all_3 <- paste(time_l,".jpg")
    all_title <- paste0(all_0,all_1,all_2,all_3)
    pl <- plot_CFU(simu_df,cfu_title)
    pl <- plot_Spes(simu_df,all_title)}}
# we store the results in a dataframe
d_list_NOCS <- d_list[d_list$CE=="NO CS",]
d_list_NOCS$P.1 <- case_when(
  d_list_NOCS$P.1 == 10^3 ~ "10^3",
  d_list_NOCS$P.1 == 10^4 ~ "10^4",
  d_list_NOCS$P.1 == 10^5 ~ "10^5",
  d_list_NOCS$P.1 == 10^8 ~ "10^8",
  TRUE ~ as.character(d_list_NOCS$P.1)
)
d_list_NOCS$Dose_Comb <- paste(d_list_NOCS$A.1,d_list_NOCS$P.1,sep="-")
d_list_NOCS_firstP <-d_list_NOCS[d_list_NOCS$First_Treat=="Phages",] 
d_list_NOCS_firstA <-d_list_NOCS[d_list_NOCS$First_Treat=="Antibiotic",]
# i have an heatmap of bacteria eradication where on the y axis i have the phage-antibiotic dose combination
# and on the x axis i have the time delay between treatments
# we have two separate cases, one when we start with phages and the other when we start with the antibiotic
heatmap_NOCS_treatmentRec <- d_list_NOCS_firstP %>%
  dplyr::select(Dose_Comb,timelag,BacteriaEradication) %>%
  pivot_wider(names_from = timelag, values_from = BacteriaEradication, values_fill = list(BacteriaEradication = 0))  # fill with 0 if 
df_long_mild <- heatmap_NOCS_treatmentRec %>%
  pivot_longer(cols = -Dose_Comb, names_to = "timelag", values_to = "BacteriaEradication")
# Create the heatmap plot
heatmap_plot_mild <- ggplot(df_long_mild, aes(x = timelag, y = factor(Dose_Comb), fill = factor(BacteriaEradication))) +
  geom_tile(color = "black", size = 0.2) +
  scale_fill_manual(values = c("0" = "pink", "1" = "darkblue")) +   # Adjust colors (white for 0, darkblue for 1)
  theme_minimal() +
  theme_minimal() +
  labs(x = "Time delay between treatments (H)", y = "Antibiotic-phage dose", fill = "Bacteria eradication") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 20),  # Increase size of x-axis labels
    axis.title.x = element_text(size = 20),  # Increase size of x-axis title
    axis.text.y = element_text(size = 14),  # Reduce size of y-axis labels (adjust to your preference)
    axis.title.y = element_text(size = 24),  # Increase size of y-axis title (adjust to your preference)
    panel.grid.major = element_line(color = "black", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "gray", size = 0.5)  # Minor grid lines
  )
jpeg(filename = paste0(output_dir,"/TO_noCS_delayTreat_firstP_allDoses.jpg"), 
     width = 1500, height = 1300, res = 150)  # Increased resolution for better quality
print(heatmap_plot_mild)
dev.off()

heatmap_NOCS_treatmentRec <- d_list_NOCS_firstA %>%
  dplyr::select(Dose_Comb,timelag,BacteriaEradication) %>%
  pivot_wider(names_from = timelag, values_from = BacteriaEradication, values_fill = list(BacteriaEradication = 0))  # fill with 0 if 
df_long_mild <- heatmap_NOCS_treatmentRec %>%
  pivot_longer(cols = -Dose_Comb, names_to = "timelag", values_to = "BacteriaEradication")
# Create the heatmap plot
heatmap_plot_mild <- ggplot(df_long_mild, aes(x = timelag, y = factor(Dose_Comb), fill = factor(BacteriaEradication))) +
  geom_tile(color = "black", size = 0.2) +
  scale_fill_manual(values = c("0" = "pink", "1" = "darkblue")) +   # Adjust colors (white for 0, darkblue for 1)
  theme_minimal() +
  theme_minimal() +
  labs(x = "Time delay between treatments (H)", y = "Antibiotic-phage dose", fill = "Bacteria eradication") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 20),  # Increase size of x-axis labels
    axis.title.x = element_text(size = 20),  # Increase size of x-axis title
    axis.text.y = element_text(size = 14),  # Reduce size of y-axis labels (adjust to your preference)
    axis.title.y = element_text(size = 24),  # Increase size of y-axis title (adjust to your preference)
    panel.grid.major = element_line(color = "black", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "gray", size = 0.5)  # Minor grid lines
  )
jpeg(filename = paste0(output_dir,"/TO_noCS_delayTreat_firstA_allDoses.jpg"), 
     width = 1500, height = 1300, res = 150)  # Increased resolution for better quality
print(heatmap_plot_mild)
dev.off()
