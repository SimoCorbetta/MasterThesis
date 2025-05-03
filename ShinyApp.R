# shiny app for sequential dosing, 24 H time length of the simulation and we 
# always start with phages
library(rxode2)
library(tidyverse)
library(grid)
library(tidyr)
library(shiny)
library(dplyr)
library(ggplot2)
# transit rate is equal to the number of compartments (10) divided by the latency
# Define the model function (simplified version here for the app)
# only parameter you cannot change from the linebar are the starting conditions(how many bacteria)
# 
# change model parameters here
model_function <- function(lin=10^(-7.40),burst=150,transit=25,ant_decay=0.03,phage_decay=0.07,
                           mic=1,mic_r=10,mic_rp=1,kmax=3,hill=2,bmax=10^10,p50=10^7.5,
                           qa=10^-9,qp=10^-7,g_wt=0.7,g_rp=0.7,g_ra=0.7,g_dr=0.7) {
  # Initial conditions and parameters
 
   ini({
    kmax <- kmax # kmax
    MIC <- mic # MIC of WT strain
    MIC_r <- mic_r # MIC of antibiotic resistant strain
    MIC_rP <- mic_rp # MIC of phage resistant strain
     H <- hill # Hill coefficient
    y <- ant_decay # antibiotic decay rate
    beta <- lin # linear phage adsorption rate
    b <- burst # burst size
    mu_max_u = g_wt # max growth rate of WT strain
    mu_max_rA = g_rp # max growth rate of antibiotic resistant strain
    mu_max_rP = g_ra # max growth rate of phage resistant strain
    mu_max_rAP = g_dr # max growth rate of double resistant strain
    qA = qa # acquisition rate of antibiotic resistance
    qP = qp # acquisition rate of phage resistance
    B_max = bmax
    Dec_P = phage_decay
    P50 = p50
    k = transit
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
    # change eventual initial conditions here
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
    Br_A(0) <- 0 # here
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
    Br_P(0) <- 0 # here
    Br_AP(0) <- 0 # here
    P(0) <- 0
    A(0) <- 0
  })
}
# we always start with the phages with this formulation and the time length is of 72H 
dose_adm <- function(phage_dose,antibiotic_dose,time_delay){
  dose_table<- et(
    amt = phage_dose,   # Phage dose input
    addl = 0,          # No additional doses
    ii = 24,    # Time delay input
    cmt = "P",         # Dosing compartment (P)
    evid = 1           # Event id (evid) 1 means a dosing record
  ) %>%
    et(amt = antibiotic_dose, addl = 0, ii = 24, cmt = "A", evid = 1, time = time_delay) %>%
    # to change time scale, change here 72 with the time length you prefer
    et(seq(0, 24, by = 0.1), evid = 0)
  return(dose_table)
} 



# Updated server function
# i let the user define from the control bar phage dose, antibiotic dose and time delay between phages and antibiotic
server <- function(input, output) {
  
  observeEvent(input$solve_button, {
    # Get user inputs from the UI (e.g., phage dose, antibiotic dose, time delay, lin, burst, transit, start_bact)
    phage_dose <- as.numeric(input$phage_dose)
    antibiotic_dose <- as.numeric(input$antibiotic_dose)
    time_delay <- as.numeric(input$time_delay)
    lin <- as.numeric(input$lin)    # Linear adsorption rate (lin)
    burst <- as.numeric(input$burst)  # Burst size (burst)
    transit <- as.numeric(input$transit)  # Transit rate (transit)
    ant_decay <- as.numeric(input$ant_decay) # antibiotic decay rate
    phage_decay <- as.numeric(input$phage_decay) # phage decay rate
    MIC <- as.numeric(input$mic) # MIC of the WT strain
    MIC_r <- as.numeric(input$mic_r) # MIC of the antibiotic resistant strain
    MIC_rP <- as.numeric(input$mic_rp) # MIC of the phage resistant strain
    kmax <- as.numeric(input$kmax) # kmax of the antibiotic
    hill <- as.numeric(input$hill) # hill coeffcient
    bmax <- as.numeric(input$bmax) # carrying capacity of the population
    p50 <- as.numeric(input$p50) # P50
    qa <- as.numeric(input$qa) # Rate of de-novo antibiotic resistance
    qp <- as.numeric(input$qp) # Rate of de-novo phage resistance
    g_wt <- as.numeric(input$g_wt) # growth rate of WT strain
    g_ra <- as.numeric(input$g_ra) # growth rate of antibiotic resistant strain
    g_rp <- as.numeric(input$g_rp) # growth rate of phage resistant strain
    g_dr <- as.numeric(input$g_dr) # growth rate of double resistant strain
    # Create the dosing event
    dosing_table <- dose_adm(phage_dose, antibiotic_dose, time_delay)
    
    # Solve the model with the dosing table, passing the user-defined parameters
    df_output <- model_function(lin = lin, burst = burst, transit = transit,ant_decay=ant_decay,phage_decay=phage_decay,
                                mic=MIC,mic_r=MIC_r,mic_rp=MIC_rP,kmax=kmax,hill=hill,
                                bmax=bmax,p50=p50,qa=qa,qp=qp,
                                g_wt=g_wt,g_ra=g_ra,g_rp=g_rp,g_dr=g_dr) %>% rxSolve(dosing_table, method = "lsoda")
    
    # Process and compute CFU (assuming you want the sum of all bacterial populations)
    df_output <- df_output %>%
      mutate(CFU = case_when(
        Bu + Bi_tot + Br_A + Br_P + Br_AP + Br_i_tot < 1 ~ 0,
        TRUE ~ Bu + Bi_tot + Br_A + Br_P + Br_AP + Br_i_tot
      ))
    
    # Create plot
    pl <- df_output %>% 
      ggplot(aes(x = time, y = CFU)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "green", size = 2) +
      geom_hline(yintercept = 10^7, linetype = "dashed", color = "blue", size = 2) +
      geom_hline(yintercept = 10^4, linetype = "dashed", color = "red", size = 2) +
      scale_y_log10(limits = c(10^-3, 10^11)) +
      scale_x_continuous(breaks = c(0,6,12,18, 24)) +
      labs(x = "time (h)", y = "CFU") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 35)
      )
    
    # Render plot in the UI
    output$cfu_plot <- renderPlot({
      pl
    })
  })
}


# i return the CFU plots vs time

ui <- fluidPage(
  titlePanel("Model Simulation"),
  
  sidebarLayout(
    sidebarPanel(
      # Input for Phage Dose
      numericInput("phage_dose", "Phage Dose", value = 10^3, min = 0, step = 1),
      
      # Input for Antibiotic Dose
      numericInput("antibiotic_dose", "Antibiotic Dose", value = 1, min = 0, step = 0.01),
      
      # Input for Time Delay (in hours)
      numericInput("time_delay", "Time Delay for Antibiotic (hours)", value = 1, min = 0, step = 0.01),
      
      # Input for Linear Adsorption Rate (lin)
      numericInput("lin", "Linear Adsorption Rate", value = 10^(-7.40), min = 0, step = 1e-7),
      
      # Input for Burst Size (burst)
      numericInput("burst", "Burst Size", value = 150, min = 1, step = 1),
      
      # Input for Transit Rate (transit)
      numericInput("transit", "Transit Rate", value = 25, min = 1, step = 1),
      # Input for Antibiotic decay rate (ant_decay)
      numericInput("ant_decay", "Antibiotic decay rate ", value = 0.03, min = 0.01, step = 0.01),
      # Input for Phage decay rate (phage_decay)
      numericInput("phage_decay", "Phage decay rate ", value = 0.07, min = 0.01, step = 0.01),
      # Input for MIC WT 
      numericInput("mic", "MIC of WT strain", value = 1, min = 0.01, step = 0.01),
      # Input for MIC antibiotic resistant 
      numericInput("mic_r", "MIC of antibiotic resistant strain", value = 10, min = 0.01, step = 0.01),
      # Input for MIC phage resistant 
      numericInput("mic_rp", "MIC of Phage resistant strain", value = 1, min = 0.01, step = 0.01),
      # Input for kmax 
      numericInput("kmax", "maximum killing effect constant", value = 3, min = 1, step = 0.1),
      # Input for Hill coefficient 
      numericInput("hill", "Hill Coefficient", value = 2, min = 0.5, step = 0.1),
      # Input for bmax 
      numericInput("bmax", "Carrying capacity", value = 10^10, min = 10^10, step = 1000),
      # Input for p50 
      numericInput("p50", "P50", value = 10^7.5, min = 10^6, step = 1000),
      # Input for qa 
      numericInput("qa", "Mutation rate for antibiotic resistance", value = 10^-9, min = 10^-15, step = 10^-18),
      # Input for qp 
      numericInput("qp", "Mutation rate for phage resistance", value = 10^-7, min = 10^-13, step = 10^-16),
      # Input for growth rate WT strain
      numericInput("g_wt", "max growth rate of WT strain", value = 0.7, min = 0.1, step = 0.01),
      # Input for growth rate phage resistant strain
      numericInput("g_rp", "max growth rate of phage resistant strain", value = 0.7, min = 0.1, step = 0.01),
      # Input for growth rate antibiotic resistant strain
      numericInput("g_ra", "max growth rate of antibiotic resistant strain", value = 0.7, min = 0.1, step = 0.01),
      # Input for growth rate double resistant strain
      numericInput("g_dr", "max growth rate of double resistant strain", value = 0.7, min = 0.1, step = 0.01),
      # Input for growth rate double resistant strain
      numericInput("g_dr", "max growth rate of double resistant strain", value = 0.7, min = 0.1, step = 0.01),
      # Button to solve the model and plot the results
      actionButton("solve_button", "Solve Model")
    ),
    
    mainPanel(
      # Plot output for CFU vs Time graph
      plotOutput("cfu_plot")
    )
  )
)


shinyApp(ui = ui, server = server)



 




