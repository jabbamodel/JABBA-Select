#-------------------------------------------------------------------
# Stochastic Age-Structured Simulator (runs from asmsim2jabba_data.R)  
# Based on the R package CCSRA (Thorson & Cope 2015)
# Modified by Henning Winker
# henning.winker@gmail.com
#-------------------------------------------------------------------


OM_jabba_Fn <-function (F_method, Nyears, AgeMax, SigmaR, M, F1, W_a, S_a, 
                        Mat_a, h, SB0, Frate, Fequil, SigmaF, Fdynamics = "Endogenous", 
                        Fmax = NA, Ncomp_per_year, SurveyCV, Recruitment_dynamics = "Stationary", 
                        Regime_multiplier = NA, dS_y = 0, q = q) 
{
  Cw_t = SB_t = F_t = Bexploit_t = rep(NA, Nyears)
  Zn_at = Dn_at = Cn_at = N_at = matrix(NA, nrow = AgeMax + 
                                          1, ncol = Nyears)
  if (Fdynamics == "Ramp") 
    Framp_t = c(rampup = seq(F1, Fmax, length = floor(Nyears/2)), 
                peak = rep(Fmax, floor((Nyears - floor(Nyears/2))/2)), 
                managed = rep(Fmax/3, Nyears - floor(Nyears/2) - 
                                floor((Nyears - floor(Nyears/2))/2)))
  if (Recruitment_dynamics == "Stationary") {
    Regime_multiplier = rep(0, Nyears + AgeMax)
  }
  if (Recruitment_dynamics == "Regime") {
    if (is.na(Regime_multiplier) || length(Regime_multiplier) != 
        (Nyears + AgeMax)) 
      stop("Please specify `RegimeRmultiplier` for regime shift dynamics")
  }
  RecDev = rnorm(Nyears + AgeMax, mean = 0, sd = SigmaR) + 
    Regime_multiplier
  RecMult = exp(RecDev - SigmaR^2/2)
  if (Fdynamics == "Endogenous") 
    F_t[1] = F1
  if (Fdynamics == "Ramp") 
    F_t[1] = Framp_t[1]
  N_at[, 1] = R0 * exp(-M * 0:AgeMax) * RecMult[(AgeMax + 1):1]
  SB_t[1] = sum(N_at[, 1] * W_a * Mat_a)
  Bexploit_t[1] = sum(N_at[, 1] * W_a * S_a[,1])
  if (F_method == -1 | F_method == 1) {
    Zn_at[, 1] = N_at[, 1] * (1 - exp(-M - F_t[1] * S_a[,1]))
    Dn_at[, 1] = Zn_at[, 1] * (M)/(M + F_t[1] * S_a[,1])
    Cn_at[, 1] = Zn_at[, 1] * (F_t[1] * S_a[,1])/(M + F_t[1] * 
                                                    S_a[,1])
  }
  if (F_method == -2 | F_method == 2) {
    Dn_at[, 1] = N_at[, 1] * (1 - exp(-M/2))
    Cn_at[, 1] = N_at[, 1] * exp(-M/2) * F_t[1] * S_a[,1]
    Dn_at[, 1] = Dn_at[, 1] + N_at[, 1] * exp(-M/2) * (1 - 
                                                         F_t[1] * S_a[,1]) * (1 - exp(-M/2))
    Zn_at[, 1] = Dn_at[, 1] + Cn_at[, 1]
  }
  for (YearI in 2:Nyears) {
    #><> Add selectivity change index
    
    if (Fdynamics == "Endogenous") {
      F_t[YearI] = F_t[YearI - 1] * (SB_t[YearI - 1]/(Fequil * 
                                                        SB0))^Frate * exp(rnorm(1, mean = -SigmaF^2/2, 
                                                                                sd = SigmaF))
      if (F_t[YearI] > 0.95 & (F_method == -2 | F_method == 
                               2)) 
        F_t[YearI] = 0.95
    }
    if (Fdynamics == "Ramp") {
      F_t[YearI] = Framp_t[YearI]
    }
    N_at[-1, YearI] = N_at[-(AgeMax + 1), YearI - 1] - Zn_at[-(AgeMax + 
                                                                 1), YearI - 1]
    SB_t[YearI] = sum((N_at[, YearI] * W_a * Mat_a)[-1])
    N_at[1, YearI] = 4 * h * R0 * SB_t[YearI]/(SB0 * (1 - 
                                                        h) + SB_t[YearI] * (5 * h - 1)) * RecMult[AgeMax + 
                                                                                                    YearI]
    Bexploit_t[YearI] = sum(N_at[, YearI] * W_a * S_a[,YearI])
    if (F_method == -1 | F_method == 1) {
      Zn_at[, YearI] = N_at[, YearI] * (1 - exp(-M - F_t[YearI] * 
                                                  S_a[,YearI]))
      Dn_at[, YearI] = Zn_at[, YearI] * (M)/(M + F_t[YearI] * 
                                               S_a[,YearI])
      Cn_at[, YearI] = Zn_at[, YearI] * (F_t[YearI] * S_a[,YearI])/(M + 
                                                                      F_t[YearI] * S_a[,YearI])
    }
    if (F_method == -2 | F_method == 2) {
      Dn_at[, YearI] = N_at[, YearI] * (1 - exp(-M/2))
      Cn_at[, YearI] = N_at[, YearI] * exp(-M/2) * F_t[YearI] * 
        S_a[,YearI]
      Dn_at[, YearI] = Dn_at[, YearI] + N_at[, YearI] * 
        exp(-M/2) * (1 - F_t[YearI] * S_a[,YearI]) * (1 - exp(-M/2))
      Zn_at[, YearI] = Dn_at[, YearI] + Cn_at[, YearI]
    }
  }
  Cw_t = (W_a %*% Cn_at)[1, ]
  AgeComp_at = array(0, dim = dim(N_at))
  for (YearI in 1:Nyears) {
    AgeComp_at[, YearI] = rmultinom(n = 1, size = Ncomp_per_year, 
                                    prob = Cn_at[, YearI])[, 1]
  }
  # add catchability q
  Index_t = cbind(Bexploit_t* q * exp(rnorm(Nyears, mean = 0, 
                                            sd = SurveyCV)), SurveyCV*0.5)
  DataList = list(Cw_t = Cw_t, SB_t = SB_t, F_t = F_t, N_at = N_at, 
                  AgeComp_at = AgeComp_at, Index_t = Index_t, RecDev = RecDev, 
                  RecMult = RecMult, Bexploit_t = Bexploit_t, survey_q = q, 
                  Cn_at = Cn_at, Regime_multiplier = Regime_multiplier,S_a=S_a)
  return(DataList)
}



aspm_Input_Fn <- function (Version , Method, K_prior, M_prior, h_prior, D_prior, 
                           SigmaR_prior, AgeComp_at, Index_t, Cw_t, W_a, Mat_a, RecDev_biasadj, rec_method = "dev",R0,q) 
{
  Nyears = ncol(AgeComp_at)
  MaxAge = nrow(AgeComp_at) - 1
  h_alpha = ((h_prior[1] - 0.2)/0.8) * (((h_prior[1] - 0.2)/0.8) * 
                                          (1 - ((h_prior[1] - 0.2)/0.8))/h_prior[2]^2 - 1)
  h_beta = (1 - ((h_prior[1] - 0.2)/0.8)) * (((h_prior[1] - 
                                                 0.2)/0.8) * (1 - ((h_prior[1] - 0.2)/0.8))/h_prior[2]^2 - 
                                               1)
  #ln_R0_prior = c(10, 30, R0, 999, 999, 1)
  ln_R0_prior = c(1, 3, R0, 999, 999, 1)
  
  F_t_prior = c(0, 3, 0.1, 999, 0.1, 1, Nyears)
  h_prior = c(0.2, 1, ifelse(Method == "CC", 0.9999, h_prior[1]), 
              h_alpha, h_beta, 1)
  K_prior = c(0, 1, K_prior[1], K_prior[1], K_prior[2], 4)
  
  M_prior = c(0, 1, M_prior[1], M_prior[1], M_prior[2], 4)
  D_prior = c(D_prior[1], D_prior[2], ifelse(Method == "SRA", 
                                             1, 0))
  S_a = S_a
  SigmaR_prior = c(0, 1, 0.6, SigmaR_prior[1], SigmaR_prior[2], 
                   -1)
  RecDev_prior = c(-3, 3, 0, NA, 999, 5)
  if (rec_method == "dev") 
    RecDev_prior[4] = 1
  if (rec_method == "ln_R0_ratio") 
    RecDev_prior[4] = 2
  CatchCV = 0.01
  Data = list(Nyears = Nyears, AgeMax = AgeMax, F_method = F_method, 
              CatchCV = CatchCV, ln_R0_prior = ln_R0_prior,K_prior = K_prior, M_prior = M_prior, 
              h_prior = h_prior, F_t_prior = F_t_prior, D_prior = D_prior, 
              SigmaR_prior = SigmaR_prior,RecDev_prior = RecDev_prior, RecDev_biasadj = RecDev_biasadj, 
              Cw_t = Cw_t, W_a = W_a, Mat_a = Mat_a, AgeComp_at = AgeComp_at, 
              Index_t = Index_t,S_a=S_a)
  #if (Version %in% c("CCSRA_v7")) 
  # add extra value to S50 and slope
  Parameters = list(ln_R0 = ln_R0_prior[3], ln_M = log(M_prior[3]), 
                    input_h = qlogis((h_prior[3] - 0.2)/0.8), ln_SigmaR = log(SigmaR_prior[4]), 
                    Survey_par = c(log(q), log(0.0001)), ln_F_t_input = log(rep(0.1, 
                                                                                Nyears)), Rec_par = rep(0, AgeMax + Nyears))
  Map = list()
  if (Method == "CC") {
    Map[["ln_F_t_input"]] = factor(rep(1, length(Parameters$ln_F_t_input)))
    Map[["ln_R0"]] = factor(NA)
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["input_h"]] = factor(NA)
    Map[["Survey_par"]] = factor(rep(NA, 2))
  }
  if (Method == "SRA") {
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["Survey_par"]] = factor(rep(NA, 2))
  }
  if (Method == "CCSRA") {
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["Survey_par"]] = factor(rep(NA, 2))
  }
  if (Method == "AS") {
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["Survey_par"]] = factor(c(1, NA))
  }
  if (Method == "ASSP") {
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["Survey_par"]] = factor(c(1, NA))
  }
  Random = NULL
  if (Method == "SRA") {
    Random = c(Random, "ln_F_t_input")
  }
  if (Method == "CCSRA") {
    Random = c(Random, "ln_F_t_input")
  }
  if (Method == "AS") {
    Random = c(Random, "ln_F_t_input")
  }
  if (Method == "ASSP") {
    Random = c(Random, "ln_F_t_input")
  }
  if ("RecDev_hat" %in% names(Parameters)) {
    Random = union(Random, "RecDev_hat")
  }
  if ("Rec_par" %in% names(Parameters)) {
    Random = union(Random, "Rec_par")
  }
  InputList = list(Data = Data, Parameters = Parameters, Random = Random, 
                   Map = Map)
  return(InputList)
}


Calc_estimates <-function (Obj, SD = NULL) 
{
  if (length(Obj$env$random) > 0) {
    par = Obj$env$last.par[-Obj$env$random]
  }else {
    par = Obj$env$last.par
  }
  Data = Obj$env$data
  Report = Obj$report()
  ParHat = Obj$env$parList(par)
  TB_t = as.vector(Data$W_a %*% Report$N_at)
  Yield_Fn = function(Fmean, Return_type = "Yield") {
    Data_new = Data
    Data_new[["Nyears"]] = 1000
    Data_new[["RecDev_biasadj"]] = rep(0, Data_new[["Nyears"]] + 
                                         Data_new[["AgeMax"]])
    Data_new[["Cw_t"]] = rep(1, Data_new[["Nyears"]])
    Data_new[["AgeComp_at"]] = matrix(0, nrow = Data_new[["AgeMax"]] + 
                                        1, ncol = Data_new[["Nyears"]])
    Data_new[["Index_t"]] = cbind(NA, rep(1, Data_new[["Nyears"]]))
    # ADD selectivity matrix here
    Data_new[["S_a"]] =  matrix(rep(Report$S_a[,ncol(Report$S_a)], Data_new[["Nyears"]]),nrow = Data_new[["AgeMax"]] +1, ncol = Data_new[["Nyears"]])
    
    
    ParHat_new = ParHat
    ParHat_new[["ln_F_t_input"]] = rep(log(Fmean + 1e-10), 
                                       Data_new[["Nyears"]])
    
    
    if ("RecDev_hat" %in% names(ParHat)) 
      ParHat_new[["RecDev_hat"]] = rep(0, Data_new[["Nyears"]] + 
                                         Data_new[["AgeMax"]])
    if ("Rec_par" %in% names(ParHat)) 
      ParHat_new[["Rec_par"]] = rep(0, Data_new[["Nyears"]] + 
                                      Data_new[["AgeMax"]])
    Obj_new = MakeADFun(data = Data_new, parameters = ParHat_new, 
                        inner.control = list(maxit = 1000))
    Report_new = Obj_new$report()
    if (Return_type == "Yield") 
      Return = rev(Report_new$Cw_t_hat)[1]
    if (Return_type == "Report") 
      Return = Report_new
    return(Return)
  }
  Fmsy = optimize(f = Yield_Fn, interval = c(0, 3), maximum = TRUE)$maximum
  Report_msy = Yield_Fn(Fmean = Fmsy, Return_type = "Report")
  Report_0 = Yield_Fn(Fmean = 0, Return_type = "Report")
  TBmsy = rev(as.vector(Data$W_a %*% Report_msy$N_at))[1]
  SBmsy = rev(Report_msy$SB_t)[1]
  MSY = rev(Report_msy$Cw_t_hat)[1]
  TB0 = rev(as.vector(Data$W_a %*% Report_0$N_at))[1]
  SB0 = rev(Report_0$SB_t)[1]
  Return = list(Fmsy = Fmsy, SB0 = SB0, TB0 = TB0, TB_t = TB_t, 
                D_t = Report$D_t, SB_t = Report$SB_t, F_t =  Report$F_t, MSY = MSY, TBmsy = TBmsy, 
                SBmsy = SBmsy)
  if (!is.null(SD)) {
    Summ = summary(SD)
    if ("Est. (bias.correct)" %in% colnames(Summ)) {
      Summ[, "Estimate"] = ifelse(is.na(Summ[, "Est. (bias.correct)"]), 
                                  Summ[, "Estimate"], Summ[, "Est. (bias.correct)"])
    }
    Return$ln_D_t = Summ[which(rownames(Summ) == "ln_D_t"), 
                         ]
    Return$ln_SB_t = Summ[which(rownames(Summ) == "ln_SB_t"), 
                          ]
    Return$D_t = Summ[which(rownames(Summ) == "D_t"), ]
    Return$SB_t = Summ[which(rownames(Summ) == "SB_t"), ]
    Return$Rec_t = Summ[which(rownames(Summ) == "Rec_t"),]
    Return$F_t = Summ[which(rownames(Summ) == "F_t"),]
  }
  return(Return)
}


