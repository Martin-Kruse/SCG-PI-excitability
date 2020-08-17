setwd("~/Documents/modeling/SCG_excitability_PI_metabolism/final model gKM regulation")

library(deSolve)

time <- seq(from=0, to=60, by = 0.0001)

source(file = "parameters_SCG_PI_metabolism_excitability_ode_gKM.R")

state <- c(RL_M = 0.0,
           R_M = 15.87, #15.87
           RG_GDP_M = 0.024823414763006247,
           RLG_GDP_M = 0.0,
           Gbeta_M = 0.02515331956638026,
           RLGbeta_M = 0.0,
           RGbeta_M = 1.5684107797961652E-5,
           G_GDP_M = 39.95000758156281,
           GaGTP_M = 1.9066743682731433E-4,
           Ga_GDP_M = 0.023838701032892817,
           GaGTP_PLC_M = 0,
           GaGDP_PLC_M = 0,
           PLC_M = 3.12,
           PI_M = 226975, 
           boundPIP2_M = 6464, 
           PIP_M = 4540, 
           PIP2_M = 3232, 
           Tubby_C = 0.0, # uM (set to 1.0 for simulation of Tubby-domain FRET)
           Tubby_PIP2_M = 0.0,
           DAG = 13.0,
           C1_dom_C = 0.0, #uM (set to 0.5 for simulation of C1-domain FRET)
           C1_DAG = 0.0,
           IP3_C = 0.02,
           LIBRAvIII_C = 0.0, #uM (set to 6.0 for simulation of LIBRAvIII FRET)
           IP3_LIBRAvIII_C = 0.0,
           PH_YFP_C = 0.0, #(set to 1.0 for simulation of PH-domain FRET)
           PH_YFP_PIP2 = 0.0,
           IP3_PH_YFP_C = 0.0,
           P4M_C = 0.0, #uM (set to 1.0 for simulation of P4M-domain FRET)
           P4M_PIP_M = 0.0,
           Ca_Buffer_high_aff = 9.771928, #uM
           Ca_Buffer_high_aff_Ca = 0.2280719, #uM
           Ca_Buffer_low_aff_Ca = 0.04657031, #uM
           Ca_Buffer_low_aff = 19.95343, #uM
           Ca_C = 0.002334125, #uM
           Ca_Ca_store = 361.8635, #uM
           IP3R_IP3 = 0.0,
           IP3R_IRBIT = 0.01,
           noninactiveIP3R = 0.6643741,
           inactiveIP3R = 0.3356259,
           Fura2_C = 0.0, #set to 1.0 for simulation of fura2-signals)
           Ca_Fura2_C = 0.0,
           k = 0.07769375,
           h = 0.9897161,
           m = 0.04623722,
           n = 0.003244884,
           l = 0.5065645,
           c = 1.705697e-05,
           c2 = 0.000264531,
           o = 0.001875062,
           voltage = -62.61479,
           KCNQ = 1.5314375363931854,
           KCNQPIP2M = 2.4685624636068146)
           
PI_metabolism <- function(t, state, parameters){
  with(
    as.list(c(state, parameters)),{
      
      if (t > 10 && t < 30.0) {
        oxoM_Ex <- 1.0 * 10
        gKM <- (KCNQPIP2M / KCNQPIP2M_start) * 0.000224 
        gBK <- 1.0 * 0.0008 
        gcalbar2 <- (1.0 - (RL_M / R_M_start)) * 0.00006 
      } else {
        oxoM_Ex <- 0.0 * 10
        gKM <- (KCNQPIP2M / KCNQPIP2M_start) * 0.000224 
        gBK <- 1.0 * 0.0008 
        gcalbar2 <- (1.0 - (RL_M / R_M_start)) * 0.00006
      }
      
      vtrap <- function(x,y) {
        if(x/y < 1E-6) {
          result_vtrap <- y * (1-x/y/2)
        } else {
          result_vtrap <- x / (exp(x/y) - 1)
        }
      }
      
      efun <- function(z) {
        if(z < 1E-4) {
          efun <- 1 - z/2
        } else {
          efun <- z / (exp(z) - 1)
        }
      }
      
      alpha_KM <- function(voltage) exp(1E-3 * zeta_KM * (voltage - vhalf_KM) * 9.648E4 / (8.315 * (273.16 + celsius)))
      beta_KM <- function(voltage) exp(1E-3 * zeta_KM * gm_KM * (voltage - vhalf_KM) * 9.648E4 / (8.315 * (273.16 + celsius)))
      kinf <- function() 1 / (1 + alpha_KM(voltage))
      tau_KM <- function() (beta_KM(voltage) / (q10KM * a0_KM * (1 + alpha_KM(voltage)))) / 1000
      alphah <- function(voltage) q10 * 0.128*exp((17-(voltage+65))/18)
      betah <- function(voltage) q10 * 4 / (exp((40.0 - (voltage+65))/5) + 1)
      hinf <- function() alphah(voltage) / (alphah(voltage) + betah(voltage))
      tauh <- function() (1 / (alphah(voltage) + betah(voltage))) / 1000
      vtrap_alpham <- function(voltage) vtrap((13.1-(voltage+65)), 4)
      alpham <- function(voltage) q10 * 0.32 * vtrap_alpham(voltage)
      vtrap_betam <- function(voltage) vtrap(((voltage+65)-40.1), 5)
      betam <- function(voltage) q10 * 0.28 * vtrap_betam(voltage)
      minf <- function() alpham(voltage) / (alpham(voltage) + betam(voltage))
      taum <- function() (1 / (alpham(voltage) + betam(voltage))) / 1000
      alphan <- function(voltage) exp(1E-3 * zetan * (voltage - vhalfn) * 9.648E4 / (8.315 * (273.16 + celsius)))
      betan <- function(voltage) exp(1E-3 * zetan * gmn * (voltage - vhalfn) * 9.648E4 / (8.315 * (273.16 + celsius)))
      alphal <- function(voltage) exp(1E-3 * zetal * (voltage - vhalfl) * 9.648E4 / (8.315 * (273.16 + celsius)))
      betal <- function(voltage) exp(1E-3 * zetal * gml * (voltage - vhalfl) * 9.648E4 / (8.315 * (273.16 + celsius)))
      ninf <- function() 1 / (1 + alphan(voltage))
      taun <- function() (betan(voltage) / (q10 * a0n * (1+alphan(voltage)))) / 1000
      linf <- function() 1 / (1 + alphal(voltage))
      taul <- function() (betal(voltage) / (q10 * a0l * (1 + alphal(voltage)))) / 1000
      alphac <- function(voltage) 15.69 * (-1.0 * voltage + 81.5) / (exp((-1.0 * voltage + 81.5) / 10.0)-1.0)
      betac <- function(voltage) 0.29 * exp(-1.0 * voltage / 10.86)
      alphac2 <-function(voltage) 0.055 * (-27.01 - voltage) / (exp((-27.01 - voltage) / 3.8) -1)
      betac2 <- function(voltage) 0.94 * exp((-63.01 - voltage) / 17)
      tauc <- function() ((1 / (tfa * (alphac(voltage) + betac(voltage)))) / 1000) * 10
      cinf <- function() alphac(voltage) / (alphac(voltage) + betac(voltage))
      tauc2 <- function() ((1 / (tfa2 * (alphac2(voltage) + betac2(voltage)))) / 1000) * 10
      cinf2 <- function() alphac2(voltage) / (alphac2(voltage) + betac2(voltage))
      alphaBK <- function(Ca_C, voltage) ((Ca_C / 1E3) * abar) / ((Ca_C / 1E3) + (k1 * exp(-2*d1*FARADAY*voltage/R/(273.15 + celsius))))
      betaBK <- function(Ca_C, voltage) bbar / (1 + (Ca_C / 1E3)/(k2 * exp(-2*d2*FARADAY*voltage/R/(273.15 + celsius))))
      tauo <- function() (1/(alphaBK(Ca_C, voltage) + betaBK(Ca_C, voltage))) / 1000
      oinf <- function() alphaBK(Ca_C, voltage) / (alphaBK(Ca_C, voltage) + betaBK(Ca_C, voltage))
      
      IKM <- function(k, voltage) gKM * k^st * (voltage - EK)
      IKM_basal <- function(k, voltage) gKMbasal * k^st * (voltage - EK)
      IL <- function(voltage) gL * (voltage - EL)
      IBK <- function(o, voltage) gBK * o^st * (voltage - EK)
      IBK_basal <- function(o, voltage) gBKbasal * o^st * (voltage - EK)
      IKDR <- function(n, l, voltage) (n^3) * l * gKDR * (voltage - EK)
      INa <- function(m, h, voltage) (m^3) * h * gNA * (voltage - ENa)
      ICa <- function(c, Ca_C, voltage) (gcalbar * c * c * (ki / (ki + (Ca_C / 1E3)))) * (-f * (1.0 - ((Ca_C / 1E3)/cao) * exp(voltage / f)) * efun(voltage / f))
      ICa2 <- function(c2, Ca_C, voltage) (gcalbar2 * c2 * (ki / (ki + (Ca_C / 1E3)))) * (-f2 * (1.0 - ((Ca_C / 1E3)/cao) * exp(voltage / f2)) * efun(voltage / f2))
      I_Na_Ca_Ex <- function(voltage, Ca_C) -1.0 * K2f_ex * (nai^3 * cao * exp(E_1 * voltage) - nao^3 * (Ca_C / 1E3) * exp(-E_2 * voltage))
      I_ATPase <- function(Ca_C) K2f_ATPase * (f_ATPase * (Ca_C / 1E3) / (f_ATPase * (Ca_C / 1E3) + b_ATPase))
      
      I_Na_Ca_Ex_plot <- I_Na_Ca_Ex(voltage, Ca_C) * 41
      ICa_plot <- ICa(c, Ca_C, voltage) * 41
      ICa2_plot <- ICa2(c2, Ca_C, voltage) * 41
      I_ATPase_plot <- I_ATPase(Ca_C) * 41
      IKDR_plot <- IKDR(n, l, voltage) * 41
      IKM_plot <- IKM(k, voltage) * 41
      IKM_basal_plot <- IKM_basal(k, voltage) * 41
      INa_plot <- INa(m, h, voltage) * 41
      IBK_plot <- IBK(o, voltage) * 41
      IBK_basal_plot <- IBK_basal(o, voltage) * 41
      IL_plot <- IL(voltage) * 41
      
      dRL_M <- ((KfL1 * oxoM_Ex * R_M) - (KrL1 * RL_M)) - ((KfG2 * RL_M * G_GDP_M) - (KrG2 * RLG_GDP_M)) - ((KfG2 * RL_M * Gbeta_M) - (KrG2 * RLGbeta_M))
      dR_M <- - ((KfL1 * oxoM_Ex * R_M) - (KrL1 * RL_M)) - ((KfG1 * G_GDP_M * R_M) - (KrG1 * RG_GDP_M)) - ((KfG1 * Gbeta_M * R_M) - (KrG1 * RGbeta_M))
      dRG_GDP_M <- ((KfG1 * G_GDP_M * R_M) - (KrG1 * RG_GDP_M)) - ((KfL2 * oxoM_Ex * RG_GDP_M) - (KrL2 * RLG_GDP_M)) - ((Kf_NE_RG * RG_GDP_M) - (Kr_NE_RG * GaGTP_M * RGbeta_M))
      dRLG_GDP_M <- ((KfG2 * RL_M * G_GDP_M) - (KrG2 * RLG_GDP_M)) + ((KfL2 * oxoM_Ex * RG_GDP_M) - (KrL2 * RLG_GDP_M)) - (Kf_NE_RLG * RLG_GDP_M)
      dGbeta_M <- ((Kf_NE_G * G_GDP_M) - (Kr_NE_G * Gbeta_M * GaGTP_M)) - ((KfG2 * RL_M * Gbeta_M) - (KrG2 * RLGbeta_M)) - ((KfG1 * Gbeta_M * R_M) - (KrG1 * RGbeta_M)) - (Kf_reconstitution * Gbeta_M * Ga_GDP_M)
      dRLGbeta_M <- ((KfG2 * RL_M * Gbeta_M) - (KrG2 * RLGbeta_M)) + (Kf_NE_RLG * RLG_GDP_M) + ((KfL2 * RGbeta_M * oxoM_Ex) - (KrL2 * RLGbeta_M))
      dRGbeta_M <- ((KfG1 * Gbeta_M * R_M) - (KrG1 * RGbeta_M)) + ((Kf_NE_RG * RG_GDP_M) - (Kr_NE_RG * GaGTP_M * RGbeta_M)) - ((KfL2 * RGbeta_M * oxoM_Ex) - (KrL2 * RLGbeta_M))
      dG_GDP_M <- - ((KfG1 * G_GDP_M * R_M) - (KrG1 * RG_GDP_M)) - ((KfG2 * RL_M * G_GDP_M) - (KrG2 * RLG_GDP_M)) - ((Kf_NE_G * G_GDP_M) - (Kr_NE_G * Gbeta_M * GaGTP_M)) + (Kf_reconstitution * Gbeta_M * Ga_GDP_M)
      dGaGTP_M <- ((Kf_NE_RG * RG_GDP_M) - (Kr_NE_RG * GaGTP_M * RGbeta_M)) + ((Kf_NE_G * G_GDP_M) - (Kr_NE_G * Gbeta_M * GaGTP_M)) + (Kf_NE_RLG * RLG_GDP_M) - (Kf_GTPase_Ga * GaGTP_M) - ((Kf_PLC_assoc * PLC_M * GaGTP_M) - (Kr_PLC_assoc * GaGTP_PLC_M))
      dGa_GDP_M <- (Kf_GTPase_Ga * GaGTP_M) - (Kf_reconstitution * Gbeta_M * Ga_GDP_M) + ((Kf_PLC_diss * GaGDP_PLC_M) - (Kr_PLC_diss * PLC_M * Ga_GDP_M))
      dGaGTP_PLC_M <- ((Kf_PLC_assoc * PLC_M * GaGTP_M) - (Kr_PLC_assoc * GaGTP_PLC_M)) + (Kf_NE_GaP * GaGDP_PLC_M) - (Kf_GTPase_GaP * GaGTP_PLC_M)
      dGaGDP_PLC_M <- (Kf_GTPase_GaP * GaGTP_PLC_M) - (Kf_NE_GaP * GaGDP_PLC_M) - ((Kf_PLC_diss * GaGDP_PLC_M) - (Kr_PLC_diss * PLC_M * Ga_GDP_M))
      dPLC_M <- - ((Kf_PLC_assoc * PLC_M * GaGTP_M) - (Kr_PLC_assoc * GaGTP_PLC_M)) + ((Kf_PLC_diss * GaGDP_PLC_M) - (Kr_PLC_diss * PLC_M * Ga_GDP_M))
      dPI_M <- 0
      dboundPIP2_M <- (((-1.0 + foldPIP2) * speed_buffering) * PIP2_M) - (speed_buffering * boundPIP2_M)
      dPIP_M <- (k4_rest * foldPIP2 * PI_M) - (k4P * PIP_M) - (k5k_rest * foldPIP2 * PIP_M) + (k5P_rest * foldPIP2 * PIP2_M) - (PIP_M * PLC_efficiency_PIP * (PLC_basal + K_PLC * GaGTP_PLC_M)) - ((speed_P4M_PIP * P4M_C * PIP_M) - (speed_P4M_PIP * KD_P4M_PIP * P4M_PIP_M))
      dPIP2_M <- (k5k_rest * foldPIP2 * PIP_M) - (k5P_rest * foldPIP2 * PIP2_M) -dboundPIP2_M - ((speed_Tubby_PIP2 * Tubby_C * PIP2_M) - (speed_Tubby_PIP2 * KD_Tubby_PIP2 * Tubby_PIP2_M)) - (PIP2_M * (PLC_basal + K_PLC * foldPIP2 * GaGTP_PLC_M)) - ((speed_PH_PIP2 * PH_YFP_C * PIP2_M) - (KD_PH_PIP2 * speed_PH_PIP2 * PH_YFP_PIP2)) - ((speed_PIP2_KCNQ * KCNQ * PIP2_M) - (speed_PIP2_KCNQ * KD_PIP2_KCNQ * KCNQPIP2M))
      dTubby_C <- (-1.0 * (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_C * ((speed_Tubby_PIP2 * Tubby_C * PIP2_M) - (speed_Tubby_PIP2 * KD_Tubby_PIP2 * Tubby_PIP2_M))))
      dTubby_PIP2_M <- ((speed_Tubby_PIP2 * Tubby_C * PIP2_M) - (speed_Tubby_PIP2 * KD_Tubby_PIP2 * Tubby_PIP2_M))
      dDAG <- (PIP2_M * (PLC_basal + K_PLC * foldPIP2 * GaGTP_PLC_M)) + (PIP_M * PLC_efficiency_PIP * (PLC_basal + K_PLC * GaGTP_PLC_M)) - (Kf_DAGPase * DAG) - ((Kf_C1DAG * DAG * C1_dom_C) - (KD_DAG_C1 * speed_DAG_C1 * C1_DAG))
      dC1_dom_C <- (-1.0 * (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_C * ((speed_DAG_C1 * C1_dom_C * DAG) - (speed_DAG_C1 * KD_DAG_C1 * C1_DAG))))
      dC1_DAG <- ((Kf_C1DAG * DAG * C1_dom_C) - (KD_DAG_C1 * speed_DAG_C1 * C1_DAG))
      dIP3_C <- (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_C * (PIP2_M * (PLC_basal + K_PLC * foldPIP2 * GaGTP_PLC_M))) - (K_IP3ase * IP3_C) - ((speed_LIBRAvIII * IP3_C * LIBRAvIII_C) - (speed_LIBRAvIII * KD_LIBRAvIII * IP3_LIBRAvIII_C)) - ((speed_PH_IP3 * PH_YFP_C * IP3_C) - (KD_PH_IP3 * speed_PH_IP3 * IP3_PH_YFP_C))
      dLIBRAvIII_C <- -1.0 * ((speed_LIBRAvIII * IP3_C * LIBRAvIII_C) - (speed_LIBRAvIII * KD_LIBRAvIII * IP3_LIBRAvIII_C))
      dIP3_LIBRAvIII_C <- (((speed_LIBRAvIII * LIBRAvIII_C) * IP3_C) - (speed_LIBRAvIII * KD_LIBRAvIII * IP3_LIBRAvIII_C))
      dPH_YFP_C <- -1.0 * (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_C * ((speed_PH_PIP2 * PH_YFP_C * PIP2_M) - (KD_PH_PIP2 * speed_PH_PIP2 * PH_YFP_PIP2))) - ((speed_PH_IP3 * PH_YFP_C * IP3_C) - (KD_PH_IP3 * speed_PH_IP3 * IP3_PH_YFP_C))
      dPH_YFP_PIP2 <- ((speed_PH_PIP2 * PH_YFP_C * PIP2_M) - (KD_PH_PIP2 * speed_PH_PIP2 * PH_YFP_PIP2))
      dIP3_PH_YFP_C <- (((speed_PH_IP3 * PH_YFP_C) * IP3_C) - (speed_PH_IP3 * KD_PH_IP3 * IP3_PH_YFP_C))
      dP4M_C <- (-1.0 * (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_C * ((speed_P4M_PIP * P4M_C * PIP_M) - (speed_P4M_PIP * KD_P4M_PIP * P4M_PIP_M))))
      dP4M_PIP_M <- ((speed_P4M_PIP * P4M_C * PIP_M) - (speed_P4M_PIP * KD_P4M_PIP * P4M_PIP_M))
      dCa_Buffer_high_aff <- -1.0 * ((speed_Ca_Buffer_high_aff * Ca_Buffer_high_aff * Ca_C) - (speed_Ca_Buffer_high_aff * KD_Ca_Buffer_high_aff * Ca_Buffer_high_aff_Ca))
      dCa_Buffer_high_aff_Ca <- ((speed_Ca_Buffer_high_aff * Ca_Buffer_high_aff * Ca_C) - (speed_Ca_Buffer_high_aff * KD_Ca_Buffer_high_aff * Ca_Buffer_high_aff_Ca))
      dCa_Buffer_low_aff_Ca <- ((speed_Ca_Buffer_low_aff * Ca_Buffer_low_aff * Ca_C) - (speed_Ca_Buffer_low_aff * KD_Ca_Buffer_low_aff * Ca_Buffer_low_aff_Ca))
      dCa_Buffer_low_aff <- -1.0 * ((speed_Ca_Buffer_low_aff * Ca_Buffer_low_aff * Ca_C) - (speed_Ca_Buffer_low_aff * KD_Ca_Buffer_low_aff * Ca_Buffer_low_aff_Ca))
      #dCa_C <- ((-1.0 * (KFlux_ERM_C * (-1.0 * r_KVH * vL * (1.0 - (Ca_C / Ca_Ca_store)) * r_KVH)) - (KFlux_ERM_C * ((r_KVH * vP * (Ca_C)^2) / ((kP)^2 + (Ca_C)^2))) - (KFlux_ERM_C * (-1.0 * r_KVH * Jmax2 * (IP3R_IP3 / (IP3R_IP3 + IP3R_IRBIT))^4 * (1.0 - (Ca_C / Ca_Ca_store)) * ((noninactiveIP3R / (noninactiveIP3R + inactiveIP3R)) * IP3_C * Ca_C / ((IP3_C + KD_IP3_IP3R) * (Ca_C + KD_act_Ca_IP3R)))^3)))) - ((speed_Ca_Buffer_low_aff * Ca_Buffer_low_aff * Ca_C) - (speed_Ca_Buffer_low_aff * KD_Ca_Buffer_low_aff * Ca_Buffer_low_aff_Ca)) - ((speed_Ca_Buffer_high_aff * Ca_Buffer_high_aff * Ca_C) - (speed_Ca_Buffer_high_aff * KD_Ca_Buffer_high_aff * Ca_Buffer_high_aff_Ca)) - ((speed_Ca_Fura2 * Fura2_C * Ca_C) - (speed_Ca_Fura2 * KD_Fura2 * Ca_Fura2_C))
      dCa_C <- (((((-ICa(c, Ca_C, voltage)) - I_Na_Ca_Ex(voltage, Ca_C) - ICa2(c2, Ca_C, voltage) - I_ATPase(Ca_C)) * 3.141593 * diam * diam / (2 * FARADAY_2)) * 1E3) / dt) - ((speed_Ca_Buffer_low_aff * Ca_Buffer_low_aff * Ca_C) - (speed_Ca_Buffer_low_aff * KD_Ca_Buffer_low_aff * Ca_Buffer_low_aff_Ca)) - ((speed_Ca_Buffer_high_aff * Ca_Buffer_high_aff * Ca_C) - (speed_Ca_Buffer_high_aff * KD_Ca_Buffer_high_aff * Ca_Buffer_high_aff_Ca))
      dCa_Ca_store <- ((KFlux_ERM_Ca_store * (-1.0 * r_KVH * vL * (1.0 - (Ca_C / Ca_Ca_store)) * r_KVH)) + (KFlux_ERM_Ca_store * ((r_KVH * vP * (Ca_C)^2) / ((kP)^2 + (Ca_C)^2))) + (KFlux_ERM_Ca_store * (-1.0 * r_KVH * Jmax2 * (IP3R_IP3 / (IP3R_IP3 + IP3R_IRBIT))^4 * (1.0 - (Ca_C / Ca_Ca_store)) * ((noninactiveIP3R / (noninactiveIP3R + inactiveIP3R)) * IP3_C * Ca_C / ((IP3_C + KD_IP3_IP3R) * (Ca_C + KD_act_Ca_IP3R)))^3)))
      dIP3R_IP3 <- ((speed_IP3_IRBIT * KD_IP3_IP3R_IRBIT * IP3_C * IP3R_IRBIT) - (speed_IP3_IRBIT * KD_IP3R_IRBIT * IP3R_IP3))
      dIP3R_IRBIT <- -1.0 * ((speed_IP3_IRBIT * KD_IP3_IP3R_IRBIT * IP3_C * IP3R_IRBIT) - (speed_IP3_IRBIT * KD_IP3R_IRBIT * IP3R_IP3))
      dnoninactiveIP3R <- -1.0 * (-1.0 * (KD_inh_Ca_IP3R - (Ca_C + KD_inh_Ca_IP3R) * (noninactiveIP3R / (noninactiveIP3R + inactiveIP3R))) * Kon_Ca_IP3R)
      dinactiveIP3R <- (-1.0 * (KD_inh_Ca_IP3R - (Ca_C + KD_inh_Ca_IP3R) * (noninactiveIP3R / (noninactiveIP3R + inactiveIP3R))) * Kon_Ca_IP3R)
      dFura2_C <- -1.0 * ((speed_Ca_Fura2 * Fura2_C * Ca_C) - (speed_Ca_Fura2 * KD_Fura2 * Ca_Fura2_C))
      dCa_Fura2_C <- ((speed_Ca_Fura2 * Fura2_C * Ca_C) - (speed_Ca_Fura2 * KD_Fura2 * Ca_Fura2_C))
      dk <- ((kinf() - k) / tau_KM())
      dh <- ((hinf() - h) / tauh())
      dm <- ((minf() - m) / taum())
      dn <- ((ninf() - n) / taun())
      dl <- ((linf() - l) / taul())
      dc <- ((cinf() - c) / tauc())
      dc2 <- ((cinf2() - c2) / tauc2())
      do <- ((oinf() - o) / tauo())
      dvoltage <- (-IKM(k, voltage) - IL(voltage) - IBK(o, voltage) - IKDR(n, l, voltage) - INa(m, h, voltage) - ICa(c, Ca_C, voltage) - ICa(c2, Ca_C, voltage) - I_ATPase(Ca_C) - I_Na_Ca_Ex(voltage, Ca_C) - IKM_basal(k, voltage) - IBK_basal(o, voltage))/ C
      dKCNQ <- -1.0 * ((speed_PIP2_KCNQ * KCNQ * PIP2_M) - (speed_PIP2_KCNQ * KD_PIP2_KCNQ * KCNQPIP2M))
      dKCNQPIP2M <- ((speed_PIP2_KCNQ * KCNQ * PIP2_M) - (speed_PIP2_KCNQ * KD_PIP2_KCNQ * KCNQPIP2M))
      return(list(c(dRL_M, dR_M, dRG_GDP_M, dRLG_GDP_M, dGbeta_M, dRLGbeta_M, dRGbeta_M, dG_GDP_M, dGaGTP_M, dGa_GDP_M, dGaGTP_PLC_M, dGaGDP_PLC_M, dPLC_M, dPI_M, dboundPIP2_M, dPIP_M, dPIP2_M, dTubby_C, dTubby_PIP2_M, dDAG, dC1_dom_C, dC1_DAG, dIP3_C, dLIBRAvIII_C, dIP3_LIBRAvIII_C, dPH_YFP_C, dPH_YFP_PIP2, dIP3_PH_YFP_C, dP4M_C, dP4M_PIP_M, dCa_Buffer_high_aff, dCa_Buffer_high_aff_Ca, dCa_Buffer_low_aff_Ca, dCa_Buffer_low_aff, dCa_C, dCa_Ca_store, dIP3R_IP3, dIP3R_IRBIT, dnoninactiveIP3R, dinactiveIP3R, dFura2_C, dCa_Fura2_C, dk, dh, dm, dn, dl, dc, dc2, do, dvoltage, dKCNQ, dKCNQPIP2M), I_Na_Ca_Ex_plot = I_Na_Ca_Ex_plot, ICa_plot = ICa_plot, ICa2_plot = ICa2_plot, I_ATPase_plot = I_ATPase_plot, IKDR_plot = IKDR_plot, IKM_plot = IKM_plot, INa_plot = INa_plot, IBK_plot = IBK_plot, IL_plot = IL_plot, IKM_basal_plot = IKM_basal_plot, IBK_basal_plot = IBK_basal_plot, gKM = gKM, gcalbar2 = gcalbar2, gBK = gBK))
    }
  )
}

out <- ode(y = state, times = time, func = PI_metabolism, parms = parameters)

result <- as.data.frame(out)
