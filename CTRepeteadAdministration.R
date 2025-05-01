#The model structure is the same as the sequential dosing, we changed the treatment scheme
# in this case, we give our treatment for 3 days. We will still evaluate sequential vs simultaneous,
# for sequential delay we are gonna have a time delay of 2H, that has been proven in the previous section
# to not cause failure of treatment at low antibiotic dose and to restore bacteria eradication and high antibiotic dose
# Simultaneous administration:
#1) QID   2) BID    3) TID    4) infusion
# Sequential administration
#1) Phages at 0 & antibiotic at 2H (repeat for 2 more days) 
#2) Phages at 0,12 & antibiotic at 2,14H (repeat for 2 more days)
#3) Phages at 0,8,16 & antibiotic at 2,10,18H (repeat for 2 more days)
#4) Antibiotic at 0 & phages at 2H (repeat for 2 more days) 
#5) Antibiotic at 0,12 & phages at 2,14H (repeat for 2 more days)
#6) Antibiotic at 0,8,16 & phages at 2,10,18H (repeat for 2 more days)

# 6 phage-antibiotic dose:
# 1) A.1 = 1.5 & MOI = 10^-2  2) A.1 = 1.5 & MOI = 10
# 3) A.1 = 1 & MOI = 10^-2  4) A.1 = 1 & MOI = 10
# 5) A.1 = 0.5 & MOI = 10^-2  6) A.1 = 0.5 & MOI = 10
library(rxode2)
library(tidyverse)
library(grid)
library(doParallel)
library(tidyr)
rm(list=ls())
output_dir <- "/home/corbettas/IntershipScriptsCTRepeatedDrugAdm_DoubleRes_CE"
if (!file.exists(output_dir)) {dir.create (output_dir)}
# no CS
mod_noCS <- function() {
  # Initial conditions and parameters
  ini({
    kmax <- 3 # kmax
    MIC <- 1 # MIC of WT strain
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

###  weak CS
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
# strong CS
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
# same 4 post processing functions as before
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
# same dosing function as the one for sequential dosing
dosing_function <- function(model,dosing_scheme){
  df_output <- model %>% rxSolve(dosing_scheme,method="lsoda")
  df_output_pp <- df_output
  df_output_pp <- post_process_Br_P(df_output_pp,30,1)
  df_output_pp <- post_process_Br_A(df_output_pp,30,1)
  df_output_pp <- post_process_Br_i_tot(df_output_pp,30,1)
  df_output_pp<- post_process_Br_AP(df_output_pp,30,1)
  df_output_pp  <- df_output_pp  %>%
    mutate(Bu = case_when(
      Bu < 1 ~ 0,  # If the sum of Bu, Bi_tot, Br is less than or equal to 1, set CFU to 0
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
# same plotting functions
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
# we define here all the treatment schemes
t_list <- list()
ev_QD_simultaneous_A_1_P_8 <- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 0, 24, 48
  et(amt = 1, addl = 2, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_QD_simultaneous_A_1_P_8$TypeOfTreat <- "Simultaneous"
ev_QD_simultaneous_A_1_P_8$FirstTreat <- NA
ev_QD_simultaneous_A_1_P_8$P_adm_freq <- "daily"
ev_QD_simultaneous_A_1_P_8$A_adm_freq <- "daily"
ev_QD_simultaneous_A_1_P_8$timelagBetweenTreat <- 0
ev_QD_simultaneous_A_1_P_8$SinglePhageAdm<- 10^8
ev_QD_simultaneous_A_1_P_8$SingleAntAdm <- 1
ev_QD_simultaneous_A_1_P_8$Tot_A<- 3
ev_QD_simultaneous_A_1_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_QD_simultaneous_A_1_P_8

ev_BID_simultaneous_A_1_P_8 <- et( # First dose for P at time 0
  amt  = 5*10^7,   # mg
  addl = 5,      # 5 additional doses (for 12h,24h,36,48,60)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0,12, 24,36, 48,60
  et(amt = 0.5, addl = 5, ii = 12, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_BID_simultaneous_A_1_P_8$TypeOfTreat <- "Simultaneous"
ev_BID_simultaneous_A_1_P_8$FirstTreat <- NA
ev_BID_simultaneous_A_1_P_8$P_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_1_P_8$A_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_1_P_8$timelagBetweenTreat <- 0
ev_BID_simultaneous_A_1_P_8$SinglePhageAdm<- 5*10^7
ev_BID_simultaneous_A_1_P_8$SingleAntAdm <- 0.5
ev_BID_simultaneous_A_1_P_8$Tot_A<- 3
ev_BID_simultaneous_A_1_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_BID_simultaneous_A_1_P_8

ev_TID_simultaneous_A_1_P_8 <- et( # First dose for P at time 0
  amt  = 33333333,   # mg
  addl = 8,      # 8 additional doses (for 8h,16h,24,32,40,48,56,64)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.33 of A at time 0,8,16,24,32,40,48,56,64
  et(amt = 0.333, addl = 8, ii = 8, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_TID_simultaneous_A_1_P_8$TypeOfTreat <- "Simultaneous"
ev_TID_simultaneous_A_1_P_8$FirstTreat <- NA
ev_TID_simultaneous_A_1_P_8$P_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_1_P_8$A_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_1_P_8$timelagBetweenTreat <- 0
ev_TID_simultaneous_A_1_P_8$SinglePhageAdm<- 33333333
ev_TID_simultaneous_A_1_P_8$SingleAntAdm <- 0.333
ev_TID_simultaneous_A_1_P_8$Tot_A<- 3
ev_TID_simultaneous_A_1_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_TID_simultaneous_A_1_P_8

ev_infusion_A_1_P_8 <- et( # First dose for P at time 0
  amt  = 3*10^8,   # mg
  dur=72,
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 3, dur=72,cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_infusion_A_1_P_8$TypeOfTreat <- "Simultaneous"
ev_infusion_A_1_P_8$FirstTreat <- NA
ev_infusion_A_1_P_8$P_adm_freq <- "Infusion"
ev_infusion_A_1_P_8$A_adm_freq <- "Infusion"
ev_infusion_A_1_P_8$timelagBetweenTreat <- 0
ev_infusion_A_1_P_8$SinglePhageAdm<- NA
ev_infusion_A_1_P_8$SingleAntAdm <- NA
ev_infusion_A_1_P_8$Tot_A<- 3
ev_infusion_A_1_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_infusion_A_1_P_8

ev_QD_simultaneous_A_1_P_5 <- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 0, 24, 48
  et(amt = 1, addl = 2, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_QD_simultaneous_A_1_P_5$TypeOfTreat <- "Simultaneous"
ev_QD_simultaneous_A_1_P_5$FirstTreat <- NA
ev_QD_simultaneous_A_1_P_5$P_adm_freq <- "daily"
ev_QD_simultaneous_A_1_P_5$A_adm_freq <- "daily"
ev_QD_simultaneous_A_1_P_5$timelagBetweenTreat <- 0
ev_QD_simultaneous_A_1_P_5$SinglePhageAdm<- 10^5
ev_QD_simultaneous_A_1_P_5$SingleAntAdm <- 1
ev_QD_simultaneous_A_1_P_5$Tot_A<- 3
ev_QD_simultaneous_A_1_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_QD_simultaneous_A_1_P_5

ev_BID_simultaneous_A_1_P_5 <- et( # First dose for P at time 0
  amt  = 5*10^4,   # mg
  addl = 5,      # 5 additional doses (for 12h,24h,36,48,60)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0,12, 24,36, 48,60
  et(amt = 0.5, addl = 5, ii = 12, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_BID_simultaneous_A_1_P_5$TypeOfTreat <- "Simultaneous"
ev_BID_simultaneous_A_1_P_5$FirstTreat <- NA
ev_BID_simultaneous_A_1_P_5$P_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_1_P_5$A_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_1_P_5$timelagBetweenTreat <- 0
ev_BID_simultaneous_A_1_P_5$SinglePhageAdm<- 5*10^4
ev_BID_simultaneous_A_1_P_5$SingleAntAdm <- 0.5
ev_BID_simultaneous_A_1_P_5$Tot_A<- 3
ev_BID_simultaneous_A_1_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_BID_simultaneous_A_1_P_5

ev_TID_simultaneous_A_1_P_5 <- et( # First dose for P at time 0
  amt  = 33333,   # mg
  addl = 8,      # 8 additional doses (for 8h,16h,24,32,40,48,56,64)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.33 of A at time 0,8,16,24,32,40,48,56,64
  et(amt = 0.333, addl = 8, ii = 8, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_TID_simultaneous_A_1_P_5$TypeOfTreat <- "Simultaneous"
ev_TID_simultaneous_A_1_P_5$FirstTreat <- NA
ev_TID_simultaneous_A_1_P_5$P_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_1_P_5$A_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_1_P_5$timelagBetweenTreat <- 0
ev_TID_simultaneous_A_1_P_5$SinglePhageAdm<- 33333
ev_TID_simultaneous_A_1_P_5$SingleAntAdm <- 0.333
ev_TID_simultaneous_A_1_P_5$Tot_A<- 3
ev_TID_simultaneous_A_1_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_TID_simultaneous_A_1_P_5

ev_infusion_A_1_P_5 <- et( # First dose for P at time 0
  amt  = 3*10^5,   # mg
  dur=72,
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 3, dur=72,cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_infusion_A_1_P_5$TypeOfTreat <- "Simultaneous"
ev_infusion_A_1_P_5$FirstTreat <- NA
ev_infusion_A_1_P_5$P_adm_freq <- "Infusion"
ev_infusion_A_1_P_5$A_adm_freq <- "Infusion"
ev_infusion_A_1_P_5$timelagBetweenTreat <- 0
ev_infusion_A_1_P_5$SinglePhageAdm<- NA
ev_infusion_A_1_P_5$SingleAntAdm <- NA
ev_infusion_A_1_P_5$Tot_A<- 3
ev_infusion_A_1_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_infusion_A_1_P_5

ev_QD_simultaneous_A_0.5_P_8 <- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1 of A at time 0, 24, 48
  et(amt = 0.5, addl = 2, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_QD_simultaneous_A_0.5_P_8$TypeOfTreat <- "Simultaneous"
ev_QD_simultaneous_A_0.5_P_8$FirstTreat <- NA
ev_QD_simultaneous_A_0.5_P_8$P_adm_freq <- "daily"
ev_QD_simultaneous_A_0.5_P_8$A_adm_freq <- "daily"
ev_QD_simultaneous_A_0.5_P_8$timelagBetweenTreat <- 0
ev_QD_simultaneous_A_0.5_P_8$SinglePhageAdm<- 10^8
ev_QD_simultaneous_A_0.5_P_8$SingleAntAdm <- 0.5
ev_QD_simultaneous_A_0.5_P_8$Tot_A<- 1.5
ev_QD_simultaneous_A_0.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_QD_simultaneous_A_0.5_P_8

ev_BID_simultaneous_A_0.5_P_8 <- et( # First dose for P at time 0
  amt  = 5*10^7,   # mg
  addl = 5,      # 5 additional doses (for 12h,24h,36,48,60)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0,12, 24,36, 48,60
  et(amt = 0.25, addl = 5, ii = 12, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_BID_simultaneous_A_0.5_P_8$TypeOfTreat <- "Simultaneous"
ev_BID_simultaneous_A_0.5_P_8$FirstTreat <- NA
ev_BID_simultaneous_A_0.5_P_8$P_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_0.5_P_8$A_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_0.5_P_8$timelagBetweenTreat <- 0
ev_BID_simultaneous_A_0.5_P_8$SinglePhageAdm<- 5*10^7
ev_BID_simultaneous_A_0.5_P_8$SingleAntAdm <- 0.25
ev_BID_simultaneous_A_0.5_P_8$Tot_A<- 1.5
ev_BID_simultaneous_A_0.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_BID_simultaneous_A_0.5_P_8

ev_TID_simultaneous_A_0.5_P_8 <- et( # First dose for P at time 0
  amt  = 33333333,   # mg
  addl = 8,      # 8 additional doses (for 8h,16h,24,32,40,48,56,64)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.17 of A at time 0,8,16,24,32,40,48,56,64
  et(amt = 0.17, addl = 8, ii = 8, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_TID_simultaneous_A_0.5_P_8$TypeOfTreat <- "Simultaneous"
ev_TID_simultaneous_A_0.5_P_8$FirstTreat <- NA
ev_TID_simultaneous_A_0.5_P_8$P_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_0.5_P_8$A_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_0.5_P_8$timelagBetweenTreat <- 0
ev_TID_simultaneous_A_0.5_P_8$SinglePhageAdm<- 33333333
ev_TID_simultaneous_A_0.5_P_8$SingleAntAdm <- 0.17
ev_TID_simultaneous_A_0.5_P_8$Tot_A<- 1.5
ev_TID_simultaneous_A_0.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_TID_simultaneous_A_0.5_P_8

ev_infusion_A_0.5_P_8 <- et( # First dose for P at time 0
  amt  = 3*10^8,   # mg
  dur=72,
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 1.5, dur=72,cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_infusion_A_0.5_P_8$TypeOfTreat <- "Simultaneous"
ev_infusion_A_0.5_P_8$FirstTreat <- NA
ev_infusion_A_0.5_P_8$P_adm_freq <- "Infusion"
ev_infusion_A_0.5_P_8$A_adm_freq <- "Infusion"
ev_infusion_A_0.5_P_8$timelagBetweenTreat <- 0
ev_infusion_A_0.5_P_8$SinglePhageAdm<- NA
ev_infusion_A_0.5_P_8$SingleAntAdm <- NA
ev_infusion_A_0.5_P_8$Tot_A<- 1.5
ev_infusion_A_0.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_infusion_A_0.5_P_8

ev_QD_simultaneous_A_0.5_P_5 <- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0, 24, 48
  et(amt = 0.5, addl = 2, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_QD_simultaneous_A_0.5_P_5$TypeOfTreat <- "Simultaneous"
ev_QD_simultaneous_A_0.5_P_5$FirstTreat <- NA
ev_QD_simultaneous_A_0.5_P_5$P_adm_freq <- "daily"
ev_QD_simultaneous_A_0.5_P_5$A_adm_freq <- "daily"
ev_QD_simultaneous_A_0.5_P_5$timelagBetweenTreat <- 0
ev_QD_simultaneous_A_0.5_P_5$SinglePhageAdm<- 10^5
ev_QD_simultaneous_A_0.5_P_5$SingleAntAdm <- 0.5
ev_QD_simultaneous_A_0.5_P_5$Tot_A<- 1.5
ev_QD_simultaneous_A_0.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_QD_simultaneous_A_0.5_P_5

ev_BID_simultaneous_A_0.5_P_5 <- et( # First dose for P at time 0
  amt  = 5*10^4,   # mg
  addl = 5,      # 5 additional doses (for 12h,24h,36,48,60)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.25 of A at time 0,12, 24,36, 48,60
  et(amt = 0.25, addl = 5, ii = 12, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_BID_simultaneous_A_0.5_P_5$TypeOfTreat <- "Simultaneous"
ev_BID_simultaneous_A_0.5_P_5$FirstTreat <- NA
ev_BID_simultaneous_A_0.5_P_5$P_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_0.5_P_5$A_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_0.5_P_5$timelagBetweenTreat <- 0
ev_BID_simultaneous_A_0.5_P_5$SinglePhageAdm<- 5*10^4
ev_BID_simultaneous_A_0.5_P_5$SingleAntAdm <- 0.25
ev_BID_simultaneous_A_0.5_P_5$Tot_A<- 1.5
ev_BID_simultaneous_A_0.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_BID_simultaneous_A_0.5_P_5

ev_TID_simultaneous_A_0.5_P_5 <- et( # First dose for P at time 0
  amt  = 33333,   # mg
  addl = 8,      # 8 additional doses (for 8h,16h,24,32,40,48,56,64)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.17 of A at time 0,8,16,24,32,40,48,56,64
  et(amt = 0.17, addl = 8, ii = 8, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_TID_simultaneous_A_0.5_P_5$TypeOfTreat <- "Simultaneous"
ev_TID_simultaneous_A_0.5_P_5$FirstTreat <- NA
ev_TID_simultaneous_A_0.5_P_5$P_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_0.5_P_5$A_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_0.5_P_5$timelagBetweenTreat <- 0
ev_TID_simultaneous_A_0.5_P_5$SinglePhageAdm<- 33333
ev_TID_simultaneous_A_0.5_P_5$SingleAntAdm <- 0.17
ev_TID_simultaneous_A_0.5_P_5$Tot_A<- 1.5
ev_TID_simultaneous_A_0.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_TID_simultaneous_A_0.5_P_5


ev_infusion_A_0.5_P_5 <- et( # First dose for P at time 0
  amt  = 3*10^5,   # mg
  dur=72,
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 1.5, dur=72,cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_infusion_A_0.5_P_5$TypeOfTreat <- "Simultaneous"
ev_infusion_A_0.5_P_5$FirstTreat <- NA
ev_infusion_A_0.5_P_5$P_adm_freq <- "Infusion"
ev_infusion_A_0.5_P_5$A_adm_freq <- "Infusion"
ev_infusion_A_0.5_P_5$timelagBetweenTreat <- 0
ev_infusion_A_0.5_P_5$SinglePhageAdm<- NA
ev_infusion_A_0.5_P_5$SingleAntAdm <- NA
ev_infusion_A_0.5_P_5$Tot_A<- 1.5
ev_infusion_A_0.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_infusion_A_0.5_P_5

ev_QD_simultaneous_A_1.5_P_8 <- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 0, 24, 48
  et(amt = 1.5, addl = 2, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_QD_simultaneous_A_1.5_P_8$TypeOfTreat <- "Simultaneous"
ev_QD_simultaneous_A_1.5_P_8$FirstTreat <- NA
ev_QD_simultaneous_A_1.5_P_8$P_adm_freq <- "daily"
ev_QD_simultaneous_A_1.5_P_8$A_adm_freq <- "daily"
ev_QD_simultaneous_A_1.5_P_8$timelagBetweenTreat <- 0
ev_QD_simultaneous_A_1.5_P_8$SinglePhageAdm<- 10^8
ev_QD_simultaneous_A_1.5_P_8$SingleAntAdm <- 1.5
ev_QD_simultaneous_A_1.5_P_8$Tot_A<- 4.5
ev_QD_simultaneous_A_1.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_QD_simultaneous_A_1.5_P_8

ev_BID_simultaneous_A_1.5_P_8 <- et( # First dose for P at time 0
  amt  = 5*10^7,   # mg
  addl = 5,      # 5 additional doses (for 12h,24h,36,48,60)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.75 of A at time 0,12, 24,36, 48,60
  et(amt = 0.75, addl = 5, ii = 12, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_BID_simultaneous_A_1.5_P_8$TypeOfTreat <- "Simultaneous"
ev_BID_simultaneous_A_1.5_P_8$FirstTreat <- NA
ev_BID_simultaneous_A_1.5_P_8$P_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_1.5_P_8$A_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_1.5_P_8$timelagBetweenTreat <- 0
ev_BID_simultaneous_A_1.5_P_8$SinglePhageAdm<- 5*10^7
ev_BID_simultaneous_A_1.5_P_8$SingleAntAdm <- 0.75
ev_BID_simultaneous_A_1.5_P_8$Tot_A<- 4.5
ev_BID_simultaneous_A_1.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_BID_simultaneous_A_1.5_P_8

ev_TID_simultaneous_A_1.5_P_8 <- et( # First dose for P at time 0
  amt  = 33333333,   # mg
  addl = 8,      # 8 additional doses (for 8h,16h,24,32,40,48,56,64)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0,8,16,24,32,40,48,56,64
  et(amt = 0.5, addl = 8, ii = 8, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_TID_simultaneous_A_1.5_P_8$TypeOfTreat <- "Simultaneous"
ev_TID_simultaneous_A_1.5_P_8$FirstTreat <- NA
ev_TID_simultaneous_A_1.5_P_8$P_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_1.5_P_8$A_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_1.5_P_8$timelagBetweenTreat <- 0
ev_TID_simultaneous_A_1.5_P_8$SinglePhageAdm<- 33333333
ev_TID_simultaneous_A_1.5_P_8$SingleAntAdm <- 0.5
ev_TID_simultaneous_A_1.5_P_8$Tot_A<- 4.5
ev_TID_simultaneous_A_1.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_TID_simultaneous_A_1.5_P_8

ev_infusion_A_1.5_P_8 <- et( # First dose for P at time 0
  amt  = 3*10^8,   # mg
  dur=72,
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 4.5, dur=72,cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_infusion_A_1.5_P_8$TypeOfTreat <- "Simultaneous"
ev_infusion_A_1.5_P_8$FirstTreat <- NA
ev_infusion_A_1.5_P_8$P_adm_freq <- "Infusion"
ev_infusion_A_1.5_P_8$A_adm_freq <- "Infusion"
ev_infusion_A_1.5_P_8$timelagBetweenTreat <- 0
ev_infusion_A_1.5_P_8$SinglePhageAdm<- NA
ev_infusion_A_1.5_P_8$SingleAntAdm <- NA
ev_infusion_A_1.5_P_8$Tot_A<- 4.5
ev_infusion_A_1.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_infusion_A_1.5_P_8

ev_QD_simultaneous_A_1.5_P_5 <- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 1.5 of A at time 0, 24, 48
  et(amt = 1.5, addl = 2, ii = 24, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_QD_simultaneous_A_1.5_P_5$TypeOfTreat <- "Simultaneous"
ev_QD_simultaneous_A_1.5_P_5$FirstTreat <- NA
ev_QD_simultaneous_A_1.5_P_5$P_adm_freq <- "daily"
ev_QD_simultaneous_A_1.5_P_5$A_adm_freq <- "daily"
ev_QD_simultaneous_A_1.5_P_5$timelagBetweenTreat <- 0
ev_QD_simultaneous_A_1.5_P_5$SinglePhageAdm<- 10^5
ev_QD_simultaneous_A_1.5_P_5$SingleAntAdm <- 1.5
ev_QD_simultaneous_A_1.5_P_5$Tot_A<- 4.5
ev_QD_simultaneous_A_1.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_QD_simultaneous_A_1.5_P_5

ev_BID_simultaneous_A_1.5_P_5 <- et( # First dose for P at time 0
  amt  = 5*10^4,   # mg
  addl = 5,      # 5 additional doses (for 12h,24h,36,48,60)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.75 of A at time 0,12, 24,36, 48,60
  et(amt = 0.75, addl = 5, ii = 12, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_BID_simultaneous_A_1.5_P_5$TypeOfTreat <- "Simultaneous"
ev_BID_simultaneous_A_1.5_P_5$FirstTreat <- NA
ev_BID_simultaneous_A_1.5_P_5$P_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_1.5_P_5$A_adm_freq <- "twice_a_day"
ev_BID_simultaneous_A_1.5_P_5$timelagBetweenTreat <- 0
ev_BID_simultaneous_A_1.5_P_5$SinglePhageAdm<- 5*10^4
ev_BID_simultaneous_A_1.5_P_5$SingleAntAdm <- 0.75
ev_BID_simultaneous_A_1.5_P_5$Tot_A<- 4.5
ev_BID_simultaneous_A_1.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_BID_simultaneous_A_1.5_P_5

ev_TID_simultaneous_A_1.5_P_5 <- et( # First dose for P at time 0
  amt  = 33333,   # mg
  addl = 8,      # 8 additional doses (for 8h,16h,24,32,40,48,56,64)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  # Give 0.5 of A at time 0,8,16,24,32,40,48,56,64
  et(amt = 0.5, addl = 8, ii = 8, cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_TID_simultaneous_A_1.5_P_5$TypeOfTreat <- "Simultaneous"
ev_TID_simultaneous_A_1.5_P_5$FirstTreat <- NA
ev_TID_simultaneous_A_1.5_P_5$P_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_1.5_P_5$A_adm_freq <- "3_a_day"
ev_TID_simultaneous_A_1.5_P_5$timelagBetweenTreat <- 0
ev_TID_simultaneous_A_1.5_P_5$SinglePhageAdm<- 33333
ev_TID_simultaneous_A_1.5_P_5$SingleAntAdm <- 0.5
ev_TID_simultaneous_A_1.5_P_5$Tot_A<- 4.5
ev_TID_simultaneous_A_1.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_TID_simultaneous_A_1.5_P_5

ev_infusion_A_1.5_P_5 <- et( # First dose for P at time 0
  amt  = 3*10^5,   # mg
  dur=72,
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 4.5, dur=72,cmt = "A", evid = 1, time = 0) %>%
  # Define sampling every 6 minutes for both A and P
  et(seq(0, 72, by = 0.1), evid = 0)
ev_infusion_A_1.5_P_5$TypeOfTreat <- "Simultaneous"
ev_infusion_A_1.5_P_5$FirstTreat <- NA
ev_infusion_A_1.5_P_5$P_adm_freq <- "Infusion"
ev_infusion_A_1.5_P_5$A_adm_freq <- "Infusion"
ev_infusion_A_1.5_P_5$timelagBetweenTreat <- 0
ev_infusion_A_1.5_P_5$SinglePhageAdm<- NA
ev_infusion_A_1.5_P_5$SingleAntAdm <- NA
ev_infusion_A_1.5_P_5$Tot_A<- 4.5
ev_infusion_A_1.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_infusion_A_1.5_P_5
### alternated
# P at 0,24,48H A at 2,26,50
# TotA=3 & totP = 3*10^8
ev_P_once_A_once_2h_A_1_P_8 <- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 1, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 1, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 1, evid = 1, time = 50, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_1_P_8$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_1_P_8$FirstTreat <- "Phages"
ev_P_once_A_once_2h_A_1_P_8$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1_P_8$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1_P_8$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_1_P_8$SinglePhageAdm<- 10^8
ev_P_once_A_once_2h_A_1_P_8$SingleAntAdm <- 1
ev_P_once_A_once_2h_A_1_P_8$Tot_A<- 3
ev_P_once_A_once_2h_A_1_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_1_P_8

#P at 0,12,24,36,48,60 A at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_1_P_8 <- et( # First dose for P at time 0
  amt  = 5*10^7,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.5, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 14, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 38, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 62, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_1_P_8 $TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_1_P_8 $FirstTreat <- "Phages"
ev_P_twice_A_twice_2h_A_1_P_8 $P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1_P_8 $A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1_P_8$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_1_P_8$SinglePhageAdm<- 5*10^7
ev_P_twice_A_twice_2h_A_1_P_8$SingleAntAdm <- 0.5
ev_P_twice_A_twice_2h_A_1_P_8$Tot_A <- 3
ev_P_twice_A_twice_2h_A_1_P_8$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_1_P_8 

#P at 0,8,16,24,32,40,48,56,64 A at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_1_P_8 <- et( # First dose for P at time 0
  amt  = 33333333,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.333, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 10, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 18, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 34, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 42, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 58, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 66, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_1_P_8 $TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_1_P_8 $FirstTreat <- "Phages"
ev_P_thrice_A_thrice_2h_A_1_P_8 $P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1_P_8 $A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1_P_8$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_1_P_8$SinglePhageAdm<- 33333333
ev_P_thrice_A_thrice_2h_A_1_P_8$SingleAntAdm <- 0.333
ev_P_thrice_A_thrice_2h_A_1_P_8$Tot_A <- 3
ev_P_thrice_A_thrice_2h_A_1_P_8$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_1_P_8

# totA= 3 & totP = 3*10^5
ev_P_once_A_once_2h_A_1_P_5 <- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 1, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 1, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 1, evid = 1, time = 50, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_1_P_5$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_1_P_5$FirstTreat <- "Phages"
ev_P_once_A_once_2h_A_1_P_5$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1_P_5$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1_P_5$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_1_P_5$SinglePhageAdm<- 10^5
ev_P_once_A_once_2h_A_1_P_5$SingleAntAdm <- 1
ev_P_once_A_once_2h_A_1_P_5$Tot_A<- 3
ev_P_once_A_once_2h_A_1_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_1_P_5

#P at 0,12,24,36,48,60 A at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_1_P_5 <- et( # First dose for P at time 0
  amt  = 5*10^4,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.5, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 14, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 38, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 62, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_1_P_5 $TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_1_P_5 $FirstTreat <- "Phages"
ev_P_twice_A_twice_2h_A_1_P_5 $P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1_P_5 $A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1_P_5$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_1_P_5$SinglePhageAdm<- 5*10^4
ev_P_twice_A_twice_2h_A_1_P_5$SingleAntAdm <- 0.5
ev_P_twice_A_twice_2h_A_1_P_5$Tot_A <- 3
ev_P_twice_A_twice_2h_A_1_P_5$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_1_P_5 

#P at 0,8,16,24,32,40,48,56,64 A at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_1_P_5 <- et( # First dose for P at time 0
  amt  = 33333,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.333, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 10, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 18, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 34, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 42, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 58, cmt = "A") %>%
  et(amt = 0.333, evid = 1, time = 66, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_1_P_5 $TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_1_P_5 $FirstTreat <- "Phages"
ev_P_thrice_A_thrice_2h_A_1_P_5 $P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1_P_5 $A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1_P_5$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_1_P_5$SinglePhageAdm<- 33333
ev_P_thrice_A_thrice_2h_A_1_P_5$SingleAntAdm <- 0.333
ev_P_thrice_A_thrice_2h_A_1_P_5$Tot_A <- 3
ev_P_thrice_A_thrice_2h_A_1_P_5$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_1_P_5

#totA=1.5 & totP = 3*10^8
ev_P_once_A_once_2h_A_0.5_P_8 <- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.5, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 50, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_0.5_P_8$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_0.5_P_8$FirstTreat <- "Phages"
ev_P_once_A_once_2h_A_0.5_P_8$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_0.5_P_8$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_0.5_P_8$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_0.5_P_8$SinglePhageAdm<- 10^8
ev_P_once_A_once_2h_A_0.5_P_8$SingleAntAdm <- 0.5
ev_P_once_A_once_2h_A_0.5_P_8$Tot_A<- 1.5
ev_P_once_A_once_2h_A_0.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_0.5_P_8

#P at 0,12,24,36,48,60 A at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_0.5_P_8 <- et( # First dose for P at time 0
  amt  = 5*10^7,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.25, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 14, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 38, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 62, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_0.5_P_8 $TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_0.5_P_8 $FirstTreat <- "Phages"
ev_P_twice_A_twice_2h_A_0.5_P_8 $P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_0.5_P_8 $A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_0.5_P_8$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_0.5_P_8$SinglePhageAdm<- 5*10^7
ev_P_twice_A_twice_2h_A_0.5_P_8$SingleAntAdm <- 0.25
ev_P_twice_A_twice_2h_A_0.5_P_8$Tot_A <- 1.5
ev_P_twice_A_twice_2h_A_0.5_P_8$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_0.5_P_8 

#P at 0,8,16,24,32,40,48,56,64 A at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_0.5_P_8 <- et( # First dose for P at time 0
  amt  = 33333333,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.17, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 10, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 18, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 34, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 42, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 58, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 66, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_0.5_P_8 $TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_0.5_P_8 $FirstTreat <- "Phages"
ev_P_thrice_A_thrice_2h_A_0.5_P_8 $P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_0.5_P_8 $A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_0.5_P_8$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_0.5_P_8$SinglePhageAdm<- 33333333
ev_P_thrice_A_thrice_2h_A_0.5_P_8$SingleAntAdm <- 0.17
ev_P_thrice_A_thrice_2h_A_0.5_P_8$Tot_A <- 1.5
ev_P_thrice_A_thrice_2h_A_0.5_P_8$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_0.5_P_8

# totA= 1.5 & P.1 = 10^5
ev_P_once_A_once_2h_A_0.5_P_5 <- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.5, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 50, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_0.5_P_5$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_0.5_P_5$FirstTreat <- "Phages"
ev_P_once_A_once_2h_A_0.5_P_5$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_0.5_P_5$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_0.5_P_5$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_0.5_P_5$SinglePhageAdm<- 10^5
ev_P_once_A_once_2h_A_0.5_P_5$SingleAntAdm <- 0.5
ev_P_once_A_once_2h_A_0.5_P_5$Tot_A<- 1.5
ev_P_once_A_once_2h_A_0.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_0.5_P_5

#P at 0,12,24,36,48,60 A at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_0.5_P_5 <- et( # First dose for P at time 0
  amt  = 5*10^4,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.25, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 14, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 38, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.25, evid = 1, time = 62, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_0.5_P_5 $TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_0.5_P_5 $FirstTreat <- "Phages"
ev_P_twice_A_twice_2h_A_0.5_P_5 $P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_0.5_P_5 $A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_0.5_P_5$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_0.5_P_5$SinglePhageAdm<- 5*10^4
ev_P_twice_A_twice_2h_A_0.5_P_5$SingleAntAdm <- 0.25
ev_P_twice_A_twice_2h_A_0.5_P_5$Tot_A <- 1.5
ev_P_twice_A_twice_2h_A_0.5_P_5$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_0.5_P_5 

#P at 0,8,16,24,32,40,48,56,64 A at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_0.5_P_5 <- et( # First dose for P at time 0
  amt  = 33333,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.17, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 10, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 18, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 34, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 42, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 58, cmt = "A") %>%
  et(amt = 0.17, evid = 1, time = 66, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_0.5_P_5 $TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_0.5_P_5 $FirstTreat <- "Phages"
ev_P_thrice_A_thrice_2h_A_0.5_P_5 $P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_0.5_P_5 $A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_0.5_P_5$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_0.5_P_5$SinglePhageAdm<- 33333
ev_P_thrice_A_thrice_2h_A_0.5_P_5$SingleAntAdm <- 0.17
ev_P_thrice_A_thrice_2h_A_0.5_P_5$Tot_A <- 1.5
ev_P_thrice_A_thrice_2h_A_0.5_P_5$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_0.5_P_5

# totA=4.5 & totP=3*10^8
ev_P_once_A_once_2h_A_1.5_P_8 <- et( # First dose for P at time 0
  amt  = 10^8,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 1.5, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 1.5, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 1.5, evid = 1, time = 50, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_1.5_P_8$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_1.5_P_8$FirstTreat <- "Phages"
ev_P_once_A_once_2h_A_1.5_P_8$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1.5_P_8$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1.5_P_8$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_1.5_P_8$SinglePhageAdm<- 10^8
ev_P_once_A_once_2h_A_1.5_P_8$SingleAntAdm <- 1.5
ev_P_once_A_once_2h_A_1.5_P_8$Tot_A<- 4.5
ev_P_once_A_once_2h_A_1.5_P_8$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_1.5_P_8

#P at 0,12,24,36,48,60 A at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_1.5_P_8 <- et( # First dose for P at time 0
  amt  = 5*10^7,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.75, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 14, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 38, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 62, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_1.5_P_8 $TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_1.5_P_8 $FirstTreat <- "Phages"
ev_P_twice_A_twice_2h_A_1.5_P_8 $P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1.5_P_8 $A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1.5_P_8$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_1.5_P_8$SinglePhageAdm<- 5*10^7
ev_P_twice_A_twice_2h_A_1.5_P_8$SingleAntAdm <- 0.75
ev_P_twice_A_twice_2h_A_1.5_P_8$Tot_A <- 4.5
ev_P_twice_A_twice_2h_A_1.5_P_8$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_1.5_P_8 

#P at 0,8,16,24,32,40,48,56,64 A at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_1.5_P_8 <- et( # First dose for P at time 0
  amt  = 33333333,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.5, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 10, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 18, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 34, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 42, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 58, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 66, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_1.5_P_8 $TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_1.5_P_8 $FirstTreat <- "Phages"
ev_P_thrice_A_thrice_2h_A_1.5_P_8 $P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1.5_P_8 $A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1.5_P_8$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_1.5_P_8$SinglePhageAdm<- 33333333
ev_P_thrice_A_thrice_2h_A_1.5_P_8$SingleAntAdm <- 0.5
ev_P_thrice_A_thrice_2h_A_1.5_P_8$Tot_A <- 4.5
ev_P_thrice_A_thrice_2h_A_1.5_P_8$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_1.5_P_8

# totA=4.5 & totP = 3*10^5
ev_P_once_A_once_2h_A_1.5_P_5 <- et( # First dose for P at time 0
  amt  = 10^5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 1.5, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 1.5, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 1.5, evid = 1, time = 50, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_1.5_P_5$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_1.5_P_5$FirstTreat <- "Phages"
ev_P_once_A_once_2h_A_1.5_P_5$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1.5_P_5$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1.5_P_5$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_1.5_P_5$SinglePhageAdm<- 10^5
ev_P_once_A_once_2h_A_1.5_P_5$SingleAntAdm <- 1.5
ev_P_once_A_once_2h_A_1.5_P_5$Tot_A<- 4.5
ev_P_once_A_once_2h_A_1.5_P_5$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_1.5_P_5

#P at 0,12,24,36,48,60 A at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_1.5_P_5 <- et( # First dose for P at time 0
  amt  = 5*10^4,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.75, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 14, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 38, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.75, evid = 1, time = 62, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_1.5_P_5 $TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_1.5_P_5 $FirstTreat <- "Phages"
ev_P_twice_A_twice_2h_A_1.5_P_5 $P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1.5_P_5 $A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1.5_P_5$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_1.5_P_5$SinglePhageAdm<- 5*10^4
ev_P_twice_A_twice_2h_A_1.5_P_5$SingleAntAdm <- 0.75
ev_P_twice_A_twice_2h_A_1.5_P_5$Tot_A <- 4.5
ev_P_twice_A_twice_2h_A_1.5_P_5$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_1.5_P_5 

#P at 0,8,16,24,32,40,48,56,64 A at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_1.5_P_5 <- et( # First dose for P at time 0
  amt  = 33333,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "P",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 0.5, evid = 1, time = 2, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 10, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 18, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 26, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 34, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 42, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 50, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 58, cmt = "A") %>%
  et(amt = 0.5, evid = 1, time = 66, cmt = "A") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_1.5_P_5 $TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_1.5_P_5 $FirstTreat <- "Phages"
ev_P_thrice_A_thrice_2h_A_1.5_P_5 $P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1.5_P_5 $A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1.5_P_5$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_1.5_P_5$SinglePhageAdm<- 33333
ev_P_thrice_A_thrice_2h_A_1.5_P_5$SingleAntAdm <- 0.5
ev_P_thrice_A_thrice_2h_A_1.5_P_5$Tot_A <- 4.5
ev_P_thrice_A_thrice_2h_A_1.5_P_5$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_1.5_P_5

# totA= 3, totP = 3*10^8, first A
ev_P_once_A_once_2h_A_1_P_8_firstA <- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 10^8, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 10^8, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 10^8, evid = 1, time = 50, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_1_P_8_firstA$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_1_P_8_firstA$FirstTreat <- "Antibiotic"
ev_P_once_A_once_2h_A_1_P_8_firstA$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1_P_8_firstA$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1_P_8_firstA$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_1_P_8_firstA$SinglePhageAdm<- 10^8
ev_P_once_A_once_2h_A_1_P_8_firstA$SingleAntAdm <- 1
ev_P_once_A_once_2h_A_1_P_8_firstA$Tot_A<- 3
ev_P_once_A_once_2h_A_1_P_8_firstA$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_1_P_8_firstA

#A at 0,12,24,36,48,60 P at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_1_P_8_firstA <- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 5*10^7, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 14, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 38, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 62, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_1_P_8_firstA$TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_1_P_8_firstA$FirstTreat <- "Antibiotic"
ev_P_twice_A_twice_2h_A_1_P_8_firstA$P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1_P_8_firstA$A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1_P_8_firstA$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_1_P_8_firstA$SinglePhageAdm<- 5*10^7
ev_P_twice_A_twice_2h_A_1_P_8_firstA$SingleAntAdm <- 0.5
ev_P_twice_A_twice_2h_A_1_P_8_firstA$Tot_A <- 3
ev_P_twice_A_twice_2h_A_1_P_8_firstA$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_1_P_8_firstA 

#A at 0,8,16,24,32,40,48,56,64 P at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA <- et( # First dose for P at time 0
  amt  = 0.333,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 33333333, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 10, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 18, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 34, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 42, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 58, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 66, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA$TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA$FirstTreat <- "Antibiotic"
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA$P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA$A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA$SinglePhageAdm<- 33333333
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA$SingleAntAdm <- 0.333
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA$Tot_A <- 3
ev_P_thrice_A_thrice_2h_A_1_P_8_firstA$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_1_P_8_firstA

# A at 0,24,48H P at 2,26,50
# TotA=3 & totP = 3*10^5
ev_P_once_A_once_2h_A_1_P_5_firstA <- et( # First dose for P at time 0
  amt  = 1,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 10^5, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 10^5, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 10^5, evid = 1, time = 50, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_1_P_5_firstA$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_1_P_5_firstA$FirstTreat <- "Antibiotic"
ev_P_once_A_once_2h_A_1_P_5_firstA$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1_P_5_firstA$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1_P_5_firstA$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_1_P_5_firstA$SinglePhageAdm<- 10^5
ev_P_once_A_once_2h_A_1_P_5_firstA$SingleAntAdm <- 1
ev_P_once_A_once_2h_A_1_P_5_firstA$Tot_A<- 3
ev_P_once_A_once_2h_A_1_P_5_firstA$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_1_P_5_firstA

#A at 0,12,24,36,48,60 P at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_1_P_5_firstA <- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 5*10^4, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 14, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 38, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 62, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_1_P_5_firstA$TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_1_P_5_firstA$FirstTreat <- "Antibiotic"
ev_P_twice_A_twice_2h_A_1_P_5_firstA$P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1_P_5_firstA$A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1_P_5_firstA$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_1_P_5_firstA$SinglePhageAdm<- 5*10^4
ev_P_twice_A_twice_2h_A_1_P_5_firstA$SingleAntAdm <- 0.5
ev_P_twice_A_twice_2h_A_1_P_5_firstA$Tot_A <- 3
ev_P_twice_A_twice_2h_A_1_P_5_firstA$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_1_P_5_firstA 

#A at 0,8,16,24,32,40,48,56,64 P at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA <- et( # First dose for P at time 0
  amt  = 0.333,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 33333, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 10, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 18, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 34, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 42, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 58, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 66, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA$TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA$FirstTreat <- "Antibiotic"
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA$P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA$A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA$SinglePhageAdm<- 33333
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA$SingleAntAdm <- 0.333
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA$Tot_A <- 3
ev_P_thrice_A_thrice_2h_A_1_P_5_firstA$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_1_P_5_firstA


# totA=1.5 & totP = 3*10^8
ev_P_once_A_once_2h_A_0.5_P_8_firstA <- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 10^8, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 10^8, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 10^8, evid = 1, time = 50, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_0.5_P_8_firstA$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_0.5_P_8_firstA$FirstTreat <- "Antibiotic"
ev_P_once_A_once_2h_A_0.5_P_8_firstA$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_0.5_P_8_firstA$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_0.5_P_8_firstA$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_0.5_P_8_firstA$SinglePhageAdm<- 10^8
ev_P_once_A_once_2h_A_0.5_P_8_firstA$SingleAntAdm <- 0.5
ev_P_once_A_once_2h_A_0.5_P_8_firstA$Tot_A<- 1.5
ev_P_once_A_once_2h_A_0.5_P_8_firstA$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_0.5_P_8_firstA
#A at 0,12,24,36,48,60 P at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA <- et( # First dose for P at time 0
  amt  = 0.25,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 5*10^7, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 14, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 38, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 62, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA$TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA$FirstTreat <- "Antibiotic"
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA$P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA$A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA$SinglePhageAdm<- 5*10^7
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA$SingleAntAdm <- 0.25
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA$Tot_A <- 1.5
ev_P_twice_A_twice_2h_A_0.5_P_8_firstA$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_0.5_P_8_firstA 

#A at 0,8,16,24,32,40,48,56,64 P at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA <- et( # First dose for P at time 0
  amt  = 0.17,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 33333333, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 10, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 18, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 34, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 42, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 58, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 66, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA$TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA$FirstTreat <- "Antibiotic"
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA$P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA$A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA$SinglePhageAdm<- 33333333
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA$SingleAntAdm <- 0.17
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA$Tot_A <- 1.5
ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_0.5_P_8_firstA

# totA= 1.5,totP = 3*10^5, first A
ev_P_once_A_once_2h_A_0.5_P_5_firstA <- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 10^5, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 10^5, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 10^5, evid = 1, time = 50, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_0.5_P_5_firstA$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_0.5_P_5_firstA$FirstTreat <- "Antibiotic"
ev_P_once_A_once_2h_A_0.5_P_5_firstA$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_0.5_P_5_firstA$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_0.5_P_5_firstA$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_0.5_P_5_firstA$SinglePhageAdm<- 10^5
ev_P_once_A_once_2h_A_0.5_P_5_firstA$SingleAntAdm <- 0.5
ev_P_once_A_once_2h_A_0.5_P_5_firstA$Tot_A<- 1.5
ev_P_once_A_once_2h_A_0.5_P_5_firstA$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_0.5_P_5_firstA

#A at 0,12,24,36,48,60 P at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA <- et( # First dose for P at time 0
  amt  = 0.25,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 5*10^4, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 14, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 38, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 62, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA$TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA$FirstTreat <- "Antibiotic"
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA$P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA$A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA$SinglePhageAdm<- 5*10^4
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA$SingleAntAdm <- 0.25
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA$Tot_A <- 1.5
ev_P_twice_A_twice_2h_A_0.5_P_5_firstA$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_0.5_P_5_firstA 

#A at 0,8,16,24,32,40,48,56,64 P at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA <- et( # First dose for P at time 0
  amt  = 0.17,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 33333, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 10, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 18, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 34, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 42, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 58, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 66, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA$TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA$FirstTreat <- "Antibiotic"
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA$P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA$A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA$SinglePhageAdm<- 33333
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA$SingleAntAdm <- 0.17
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA$Tot_A <- 1.5
ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_0.5_P_5_firstA

# totA = 4.5, totP = 3*10^8, first antibiotic
ev_P_once_A_once_2h_A_1.5_P_8_firstA <- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 10^8, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 10^8, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 10^8, evid = 1, time = 50, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_1.5_P_8_firstA$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_1.5_P_8_firstA$FirstTreat <- "Antibiotic"
ev_P_once_A_once_2h_A_1.5_P_8_firstA$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1.5_P_8_firstA$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1.5_P_8_firstA$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_1.5_P_8_firstA$SinglePhageAdm<- 10^8
ev_P_once_A_once_2h_A_1.5_P_8_firstA$SingleAntAdm <- 1.5
ev_P_once_A_once_2h_A_1.5_P_8_firstA$Tot_A<- 4.5
ev_P_once_A_once_2h_A_1.5_P_8_firstA$Tot_P<- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_1.5_P_8_firstA

#A at 0,12,24,36,48,60 P at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA <- et( # First dose for P at time 0
  amt  = 0.75,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 5*10^7, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 14, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 38, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 5*10^7, evid = 1, time = 62, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA$TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA$FirstTreat <- "Antibiotic"
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA$P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA$A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA$SinglePhageAdm<- 5*10^7
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA$SingleAntAdm <- 0.75
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA$Tot_A <- 4.5
ev_P_twice_A_twice_2h_A_1.5_P_8_firstA$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_1.5_P_8_firstA 

#A at 0,8,16,24,32,40,48,56,64 P at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA <- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 33333333, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 10, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 18, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 34, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 42, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 58, cmt = "P") %>%
  et(amt = 33333333, evid = 1, time = 66, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA$TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA$FirstTreat <- "Antibiotic"
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA$P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA$A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA$SinglePhageAdm<- 33333333
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA$SingleAntAdm <- 0.5
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA$Tot_A <- 4.5
ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA$Tot_P <- 3*10^8
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_1.5_P_8_firstA

# totA=4.5 & totP = 3*10^5
ev_P_once_A_once_2h_A_1.5_P_5_firstA <- et( # First dose for P at time 0
  amt  = 1.5,   # mg
  addl = 2,      # Two additional doses (for 24h and 48h)
  ii   = 24,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 10^5, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 10^5, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 10^5, evid = 1, time = 50, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_once_A_once_2h_A_1.5_P_5_firstA$TypeOfTreat <- "Alternated"
ev_P_once_A_once_2h_A_1.5_P_5_firstA$FirstTreat <- "Antibiotic"
ev_P_once_A_once_2h_A_1.5_P_5_firstA$P_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1.5_P_5_firstA$A_adm_freq <- "daily"
ev_P_once_A_once_2h_A_1.5_P_5_firstA$timelagBetweenTreat <- NA
ev_P_once_A_once_2h_A_1.5_P_5_firstA$SinglePhageAdm<- 10^5
ev_P_once_A_once_2h_A_1.5_P_5_firstA$SingleAntAdm <- 1.5
ev_P_once_A_once_2h_A_1.5_P_5_firstA$Tot_A<- 4.5
ev_P_once_A_once_2h_A_1.5_P_5_firstA$Tot_P<- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_once_A_once_2h_A_1.5_P_5_firstA
#A at 0,12,24,36,48,60 P at 2,14,26,38,50,62
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA <- et( # First dose for P at time 0
  amt  = 0.75,   # mg
  addl = 5,      # Two additional doses (for 24h and 48h)
  ii   = 12,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 5*10^4, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 14, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 38, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 5*10^4, evid = 1, time = 62, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA$TypeOfTreat <- "Alternated"
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA$FirstTreat <- "Antibiotic"
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA$P_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA$A_adm_freq <- "2_a_day"
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA$timelagBetweenTreat <- NA
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA$SinglePhageAdm<- 5*10^4
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA$SingleAntAdm <- 0.75
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA$Tot_A <- 4.5
ev_P_twice_A_twice_2h_A_1.5_P_5_firstA$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_twice_A_twice_2h_A_1.5_P_5_firstA 

#A at 0,8,16,24,32,40,48,56,64 P at 2,10,18,26,34,42,50,58,66
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA <- et( # First dose for P at time 0
  amt  = 0.5,   # mg
  addl = 8,      # Two additional doses (for 24h and 48h)
  ii   = 8,     # Dosing interval (24 hours between doses)
  cmt  = "A",     # Dosing compartment (P)
  evid = 1        # Event id (evid) 1 means a dosing record
) %>%
  et(amt = 33333, evid = 1, time = 2, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 10, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 18, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 26, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 34, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 42, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 50, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 58, cmt = "P") %>%
  et(amt = 33333, evid = 1, time = 66, cmt = "P") %>%
  et(seq(0, 72, by = 0.1), evid = 0)
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA$TypeOfTreat <- "Alternated"
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA$FirstTreat <- "Antibiotic"
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA$P_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA$A_adm_freq <- "3_a_day"
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA$timelagBetweenTreat <- NA
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA$SinglePhageAdm<- 33333
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA$SingleAntAdm <- 0.5
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA$Tot_A <- 4.5
ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA$Tot_P <- 3*10^5
t_list[[lengtht_list+1]] <- ev_P_thrice_A_thrice_2h_A_1.5_P_5_firstA
# here i call the dosing function for each treatment schemes
# no CS
# i just need to call the dosing function on the weak CS and strong CS function
d_list <- rep(0,16)
for (t in t_list){
  type_of_treat = t$TypeOfTreat[1]
  if (type_of_treat=="Simultaneous"){
    freq = t$P_adm_freq[1]
    totAnt = t$Tot_A[1]
    totPhages = t$Tot_P[1]
    if (totPhages == 3*10^5){
      totPhages_string = "3e05"
    }
    else {
      totPhages_string = "3e08"
    }
    simu_df <- dosing_function(mod_noCS,t)
    #print(colnames(simu_df))
    summary_df <- simu_df[1,40:48]
    summary_df$TypeofTreat <- type_of_treat
    summary_df$FirstTreat <- NA
    summary_df$P_adm_freq <- freq
    summary_df$A_adm_freq <- freq
    summary_df$Tot_A <- totAnt
    summary_df$Tot_P <- totPhages
    summary_df$timelag <- NA
    summary_df$CE <- "NO CS"
    d_list <- rbind(d_list,summary_df)
    CFU_0 <- paste("/CFU-time_",type_of_treat)
    CFU_1 <- paste("_",freq)
    CFU_2 <- paste("_",as.character(totAnt),"_",totPhages_string,"_NO_CS.jpg")
    cfu_title <- paste0(CFU_0,CFU_1,CFU_2)
    
    all_0 <- paste("/AllSpes-time_",type_of_treat)
    all_1 <- paste("_",freq)
    all_2 <- paste("_",as.character(totAnt),"_",totPhages_string,"_NO_CS.jpg")
    all_title <- paste0(all_0,all_1,all_2)
    pl <- plot_CFU(simu_df,cfu_title)
    pl <- plot_Spes(simu_df,all_title)
  }
  else {
    first_t = t$FirstTreat[1]
    if (first_t=="Phages"){
      totAnt = t$Tot_A[1]
      totPhages = t$Tot_P[1]
      if (totPhages == 3*10^5){
        totPhages_string = "3e05"
      }
      else {
        totPhages_string = "3e08"
      }
      simu_df <- dosing_function(mod_noCS,t)
      p_adm_freq <- t$P_adm_freq[1]
      a_adm_freq <- t$A_adm_freq[1]
      time_lag <- t$timelagBetweenTreat[1]
      summary_df <- simu_df[1,40:48]
      summary_df$TypeofTreat <- type_of_treat
      summary_df$FirstTreat <- "Phages"
      summary_df$Tot_A <- totAnt
      summary_df$Tot_P <- totPhages
      summary_df$P_adm_freq <- p_adm_freq
      summary_df$A_adm_freq <- a_adm_freq
      summary_df$timelag <- time_lag
      summary_df$CE <- "NO CS"
      d_list <- rbind(d_list,summary_df)
      CFU_0 <- paste("/FirstPhages_CFU-time_","P_",p_adm_freq)
      CFU_1 <- paste("_A_",a_adm_freq)
      CFU_2 <- paste("_",as.character(totAnt),"_",totPhages_string,as.character(time_lag),"_NO_CS.jpg")
      cfu_title <- paste0(CFU_0,CFU_1,CFU_2)
      all_0 <- paste("/FirstPhages_AllSpes-time_","P_",p_adm_freq)
      all_1 <- paste("_A_",a_adm_freq)
      all_2 <- paste("_",as.character(totAnt),"_",totPhages_string,as.character(time_lag),"_NO_CS.jpg")
      all_title <- paste0(all_0,all_1,all_2)
      pl <- plot_CFU(simu_df,cfu_title)
      pl <- plot_Spes(simu_df,all_title)
    }
    else {
      totAnt = t$Tot_A[1]
      totPhages = t$Tot_P[1]
      if (totPhages == 3*10^5){
        totPhages_string = "3e05"
      }
      else {
        totPhages_string = "3e08"
      }
      simu_df <- dosing_function(mod_noCS,t)
      p_adm_freq <- t$P_adm_freq[1]
      a_adm_freq <- t$A_adm_freq[1]
      time_lag <- t$timelagBetweenTreat[1]
      summary_df <- simu_df[1,40:48]
      summary_df$TypeofTreat <- type_of_treat
      summary_df$FirstTreat <- "Antibiotic"
      summary_df$Tot_A <- totAnt
      summary_df$Tot_P <- totPhages
      summary_df$P_adm_freq <- p_adm_freq
      summary_df$A_adm_freq <- a_adm_freq
      summary_df$timelag <- time_lag
      summary_df$CE <- "NO CS"
      d_list <- rbind(d_list,summary_df)
      CFU_0 <- paste("/FirstAntibiotic_CFU-time_","P_",p_adm_freq)
      CFU_1 <- paste("_A_",a_adm_freq)
      CFU_2 <- paste("_",as.character(totAnt),"_",totPhages_string,as.character(time_lag),"_NO_CS.jpg")
      cfu_title <- paste0(CFU_0,CFU_1,CFU_2)
      all_0 <- paste("/FirstAntibiotic_AllSpes-time_","P_",p_adm_freq)
      all_1 <- paste("_A_",a_adm_freq)
      all_2 <- paste("_",as.character(totAnt),"_",totPhages_string,as.character(time_lag),"_NO_CS.jpg")
      all_title <- paste0(all_0,all_1,all_2)
      pl <- plot_CFU(simu_df,cfu_title)
      pl <- plot_Spes(simu_df,all_title)}
  }
  
}
# store the treatment success in a dataframe
d_list_noCS_TS <- d_list[d_list$BacteriaEradication==1,]