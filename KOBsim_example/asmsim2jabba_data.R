##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><> ><><><
## asmsim2jabba_data.R
## Stock Assessment Simulation Input File generator for JABBA-SELECT 
## written by Henning Winker; henning.winker@gmail.com
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

# Delete all objects
rm(list=ls())
gc()
# Install required packages if missing
list.of.packages <- c("gplots", "coda","rjags","R2jags","fitdistrplus","reshape","mvtnorm","scales","devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# Load Packages
library(gplots);library(coda);library(rjags);library(R2jags);library("fitdistrplus");library(reshape);
library(mvtnorm);library(scales);library(devtools)

if("ThorsonUtilities" %in% installed.packages()[,"Package"]==FALSE)install_github("James-Thorson/utilities")
library(ThorsonUtilities)

#---------------------------------------------------------------------
# Set Working directory file where to store the results
File = "C:/Work/Research/GitHub/JABBA-SELECT/KOBsim_example"
# Set Assessment
assessment = "KOBsim"
# Create directory
assessment.dir = paste0(File,"/",assessment) 
dir.create(assessment.dir,showWarnings = FALSE)

# Get adjusted function based CCSRA (Thorson and Cope 2015)
source(paste0(File,"/OM_ccsra2jabba_Fn.R"))

# Set up stock parameters    
AgeMax = 20
Nyears = 40
Ncomp_per_year = 1e2 # not used
SurveyCV = 0.2 # CV on CPUE
Nrep = 1

Recruitment_dynamics = c("Stationary", "Regime")[1]
RegimeMult = c( -0.33, 0.33 )  # Additive effect on rec-devs for regimes 1 and 2
SimSeed = ceiling(runif(1,min=0,max=1e6))
if(Recruitment_dynamics=="Stationary") RegimeMult[] = 0

#========================
# Biological parameters
#========================
SigmaR = 0.6 # Recruitment variation
# VBGF
K = 0.115
Linf = 1372
t0 = -0.815
# Weight-Length
W_alpha = 0.000006
W_beta = 3.07
# Age-at-Maturity
Amat = 3
# Natural Mortality
M = 0.18
# Steepness
h = 0.8
R0 = 1.5 # Unfished mean recruitment 

# Selectivity parameters 
SL50 = S50 = c(300,500) #><> two value option for change point in selectivity
Sslope = 1/(SL50 * 0.05)
dS_y = 26 # year of selectivity change (input control)
q = 0.05 # ><> catchability
# Effort/Fishing dynamics  parameters
F1 = 0.01
Fequil = 0.17 #changed
Frate = 0.14
SigmaF = 0.15

# Derive Popdyns
L_a = Linf *(1- exp(-K*(0:AgeMax-t0)))
W_a = W_alpha * L_a^W_beta      # In grams
S_a = matrix(NA, nrow = AgeMax +1, ncol = Nyears)
  for(YearI in 1:Nyears){
    i = ifelse(YearI<dS_y,1,2)
    S_a[,YearI] = 1 / (1 + exp( -Sslope[i] * (L_a - SL50[i]) ))
  }
  Mat_a = ifelse(0:AgeMax>Amat,1,0)
  h = 0.8
  SB0 = sum( R0 * exp(-M * 0:AgeMax) * W_a * Mat_a )
  SBPR0 = SB0 / R0
  
  # Set seed
  set.seed( SimSeed + runif(123))

  # Decide on regimes (No Regime shift)
  Regime_multiplier = RegimeMult[sample(1:length(RegimeMult),replace=FALSE)]
  Regime_multiplier = c( rep(Regime_multiplier[1], each=floor(AgeMax+Nyears/2)), rep(Regime_multiplier[2], each=ceiling(Nyears/2)) )
    
    condition=1
    while(condition<2){
    DataList_Original = OM_jabba_Fn( F_method=1, Nyears=Nyears, AgeMax=AgeMax, SigmaR=SigmaR, M=M, F1=F1, W_a=W_a, S_a=S_a, Mat_a=Mat_a, h=h, SB0=SB0, Frate=Frate, Fequil=Fequil, SigmaF=SigmaF, Ncomp_per_year=Ncomp_per_year, SurveyCV=SurveyCV, Recruitment_dynamics=Recruitment_dynamics, Regime_multiplier=Regime_multiplier,dS_y = dS_y, q=q)
    if(min(DataList_Original$SB_t)>1000) condition = condition+1 # add condition
    }
    DataList_Original$SB0 = SB0
    save( DataList_Original, file=paste0(File,"/DataList_Original.RData") )
  

    # Plot time series
    ThorsonUtilities::save_fig( file=paste0(File,"/Dynamics"), width=4, height=4 )
    par( yaxs="i", xaxs="i", mgp=c(2,0.5,0), tck=-0.02, mar=c(4,2,1,1)  )
    Mat_tz = cbind(DataList_Original[['SB_t']]/SB0,DataList_Original[['F_t']],DataList_Original[['Cw_t']]/max(DataList_Original[['Cw_t']]),DataList_Original[['N_at']][1,]/R0,exp(DataList_Original[["Regime_multiplier"]][-c(1:AgeMax)]))
    matplot( Mat_tz, type="l", col=c("black","red","blue","darkgreen","lightgreen"), lty="solid", xlab="Year", ylab="", ylim=c(0,max(c(1.2,Mat_tz))) )
    legend("top", bty="n", fill=c("black","red","blue","darkgreen","lightgreen"), legend=c("Depletion","F_t","Relative catch","R_t","Regime_t"), ncol=2 )
    dev.off()

  
    
    DataList = DataList_Original
    
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
    # Create JABBA-SELECT input files
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
    
    
    
    # For GitHub Example add one more index with higher CV and different q
    DataList$Index2_t = cbind(DataList$Bexploit_t* 0.02 * exp(rnorm(Nyears, mean = 0, 
                                                        sd = 0.35)), 0.2)
    # Create Catch, CPUE and SE input files
    cpue = se = matrix(0,Nyears,4)
    catch = matrix(0,Nyears,2)
    for(y in 1:Nyears){
    if(y<dS_y){
      cpue[y,] = c(DataList$Index_t[y,1],NA,DataList$Index2_t[y,1],NA)   
      se[y,] = c(DataList$Index_t[y,2],NA,DataList$Index2_t[y,2],NA)
      catch[y,] = c(DataList$Cw_t[y],0)
    } else {
      cpue[y,] = c(NA,DataList$Index_t[y,1],NA,DataList$Index2_t[y,1])   
      se[y,] = c(NA,DataList$Index_t[y,2],NA,DataList$Index2_t[y,2])
      catch[y,] = c(0,DataList$Cw_t[y])
    }}
    
    cpue = data.frame(Year=1:Nyears,cpue)
    se = data.frame(Year=1:Nyears,se)
    catch = data.frame(Year=1:Nyears,catch)
    colnames(cpue) <- c("Year",paste0("A1.S",c("1","2")),paste0("A2.S",c("1","2")))
    colnames(se) <- names(cpue) 
    colnames(catch) <- c("Year",paste0("S",1:2))
    names(catch)
    names(cpue)
    
    # Create Select Table
    select = data.frame(Fleet.id=1:6,Fleet=c(names(catch[,-1]),names(cpue[,-1])),Selectivity=rep(1:2,3),CPUE.units=1,Catch=c(rep(TRUE,2),rep(FALSE,4)),CPUE=c(rep(FALSE,2),rep(TRUE,4)),q=c(0,0,1,1,2,2))
    # Create Selex (Selectivity) Table
    
    
    # JABBA-Selex function
    jabba.selex <- function(pars,dat){
      L = dat[,1]
      sel = dat[,2]/max(dat[,2])
      SL50	= pars[1]
      SL95	=  pars[2]
      SL.desc =  pars[3] 	
      CV.desc	=  pars[4]
      min.desc = 	 pars[5]
      psel_a = 1/(1+exp(-log(19)*(L-SL50)/(SL95-SL50)))
      psel_b = dnorm(L,SL.desc,CV.desc*SL.desc)/max(dnorm(L,SL.desc,CV.desc*SL.desc))
      psel_c = 1+(min.desc-1)*(psel_b-1)/-1
      psel = ifelse(L<SL.desc,psel_a,psel_c)
      resids = sel-psel
      return(list(ll=sum(resids^2),results=data.frame(L=L,obs=sel,fit=psel,logis=psel_a,halfnorm=psel_b,height=psel_c)))
    }
    # Likelihood
    jsel.ll = function(pars,dat){
      jabba.selex(pars,dat)$ll
    } 
    
    sel.pos = c(1,Nyears) # early and late selectivity
    Seldat = data.frame(L_a,S1=S_a[,1],S2=S_a[,Nyears])
    
    selex = data.frame(Parameter=c("SL50","SL95","SL.desc","CV.desc","min.desc"), S1=NA,S2=NA)
    
    par(mfrow=c(1,1))
    for(i in 1:2){
    peak = which((Seldat[,i+1])>0.95*max(Seldat[,i+1]))
    # get inits of jabba.selex parameters
    pars = c(SL50	= Seldat[peak[1],1]*0.9,SL95=Seldat[peak[1],1]*0.99,SL.desc=Seldat[max(peak),1],CV.desc=0.2,min.desc=0.001)
    # Minimize
    jsel.est = optim(pars, fn = jsel.ll,method="L-BFGS-B",lower=10^-3,upper=max(Seldat$L_a*1.4), dat = Seldat[,c(1,i+1)], hessian = TRUE)
    jsel.out = jabba.selex(jsel.est$par,dat=data.frame(L=seq(100,max(Seldat[,1]),1),pr=0.5))$results 
    if(i==1)plot(Seldat[,1] ,Seldat[,2],pch=1,ylab="Selectivity",xlab="Length",type="n")
    points(Seldat[,1] ,Seldat[,1+i],pch=16, col=2+i)
    lines(jsel.out$L,jsel.out$fit,lwd=2,col=2+i)
    selex[1:2,i+1] = jsel.est$par[1:2]  # only needs SL50 and SL95 for logistic 
    }
    
    # Create CSV files
    write.csv(catch,paste0(assessment.dir,"/catch",assessment,".csv"),row.names = FALSE)
    write.csv(cpue,paste0(assessment.dir,"/cpue",assessment,".csv"),row.names = FALSE)
    write.csv(se,paste0(assessment.dir,"/se",assessment,".csv"),row.names = FALSE)
    write.csv(selex,paste0(assessment.dir,"/selex",assessment,".csv"),row.names = FALSE)
    write.csv(select,paste0(assessment.dir,"/select",assessment,".csv"),row.names = FALSE)
    
    