##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
## Prime File for SWOss3 JABBA-SELECT example
## written by Henning Winker
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Delete all objects
rm(list=ls())
gc()
# Install required packages if missing
list.of.packages <- c("gplots", "coda","rjags","R2jags","fitdistrplus","reshape","mvtnorm","scales","reshape2","r4ss")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load Packages
library(gplots);library(coda);library(rjags);library(R2jags);library("fitdistrplus");library(reshape);
library(mvtnorm);library(scales);library(reshape2);library(r4ss)


#---------------------------------------------------------------------
# Set Working directory file where to store the results
File = "C:/Work/Research/GitHub/JABBA-SELECT"
# Set working directory for JABBA R source code
JABBA.file = "C:/Work/Research/GitHub/JABBA-SELECT"
# Set Assessment
assessment = "SWOss3"
# Version
version = "v1.1"

#-------------------------------------------
# Load ss2jabba rdata for JABBA-Select
#-------------------------------------------

load(paste0(File,"/",assessment,"/ss4js_",assessment,".rdata"),verbose=T)

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Graphic, Output, Saving (.RData) settings 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
KOBE.plot = TRUE # Produces JABBA Kobe plot 
KOBE.type = c("ICCAT","IOTC")[2] # ICCAT uses 3 colors; IOTC 4 (incl. orange) 
Biplot= TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
SP.plot = c("standard","phase")[2] # Produces standard or 'Kobe phase' SP plot  
save.trajectories =TRUE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
harvest.label = c("Hmsy","Fmsy")[2] # choose label preference H/Hmsy versus Fmsy
CPUE.plot= TRUE # Runs state-tool to produce "alligned" multi-CPUE plot  
meanCPUE = FALSE # Uses averaged CPUE from state-space tool instead of individual indices  
Projection = TRUE # Use Projections: requires to define TACs vectors 
save.projections = TRUE # saves projection posteriors as .RData object 
catch.metric = "(t)" # Define catch input metric e.g. (tons) "000 t" etc 
Reproduce.seed = FALSE # If FALSE a random seed assigned to each run, if TRUE set.seed(123)
runASEM=TRUE 
jabba2FRL = TRUE
# Save entire posterior as .RData object
save.all = FALSE # (if TRUE, a very large R object of entire posterior is saved)  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


# Choose Sceanrio name for creating a seperate folder
# Scenario 1: SELECT option TRUE
# Scenario 2: Pella-Tomlinson - ignoring selectivity

Scenarios = c("SWOselect","SWOpella")
s=1

for(s in 1:length(Scenarios)){
Scenario = Scenarios[s] 
if(s==2) SELECT = FALSE  

 #--------------------------------------------------
 # Read csv files
 #--------------------------------------------------
 # Load assessment data
 catch = read.csv(paste0(File,"/",assessment,"/catch",assessment,".csv"))
 cpue = read.csv(paste0(File,"/",assessment,"/cpue",assessment,".csv"))#

 # Use SEs from csv file for abudance indices (TRUE/FALSE)
 SE.I = FALSE
 
 if(SE.I ==TRUE){
  se = read.csv(paste0(File,"/",assessment,"/se",assessment,".csv"))
 }
 
 
 # Read select csv
 select = read.csv(paste0(File,"/",assessment,"/select",assessment,".csv"))
 
 # Read selex (selectivity) csv
 selex = read.csv(paste0(File,"/",assessment,"/selex",assessment,".csv"))
 
 names(cpue)
 names(catch)

  #---------------------------------------
  # option to exclude CPUE time series 
  #---------------------------------------
  
  
  #------------------------------------------------------
  # Option use mean CPUE from state-space cpue averaging
  #-----------------------------------------------------
  # Produce CPUE plot average plot
  CPUE.plot = TRUE
  
  #------------------------------------------------------
  # mean and CV and sd for unfished spawning biomass SB0
  #------------------------------------------------------
  mu.SB0 = 200000; CV.SB0 = 2; sd.SB0=sqrt(log(CV.SB0^2+1)) 
  
  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  # Set the initial depletion prior SB1/SB0 
  # To be converted into a lognormal prior (with upper bound at 1.1)
  # psi.prior = "lnorm"
  # or to be converted into a Beta prior
  # psi.prior = "beta"
  
  psi.prior= "lnorm"
  # specify as mean and CV 
  mu.psi = 1
  CV.psi = 0.05
  
  P_bound = c(0.03,1.2)
  #--------------------------------------------------------------
  # Determine estimation for catchability q 
  #--------------------------------------------------------------
  # Assign q to abundance
  sets.q = 1:(ncol(cpue)-1)  # here 1: South Early+Recent, 2: South-East Early+Recent  
  
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Observation Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  
  #To Estimate additional observation variance set sigma.est = TRUE
  sigma.est = c(TRUE,FALSE)[2]
  
  # Assign common variance estimates to abundance indices
  # Here it assumed that same flagged fleets have the same additional observation variance
  sets.var = c(1,2,2,3,3,3,4,5,6)
  
  # As option for data-weighing
  # minimum fixed observation error for each variance set (optional choose 1 value for both)
  fixed.obsE = c(0.3) # Important if SE.I is not availble
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Process Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  #Estimate set sigma.proc == True
  sigma.proc = TRUE
  proc.dev.all= 1979 # start process error in 1979 as in ss3 assessment
  #------------------------------------------
  if(sigma.proc == TRUE){
    proc.type = c("igamma","lnorm")[1] # choose 1: inverse-gamma or 2: lognormal
    
    if(proc.type=="lnorm"){
    pr.proc = c(log(0.08),0.2) # Option for lognormal process error
    }
    if(proc.type=="igamma"){
      pr.proc = c(0.001,0.001) # Option for inverse-gamma prior
      pr.proc = c(0.001,0.001)
      gamma.check = 1/rgamma(1000,pr.proc[1],pr.proc[2]) # Process error check
    # check mean process error + CV
    mu.proc = sqrt(mean(gamma.check)); CV.proc = sd(sqrt(gamma.check))/mean(sqrt(gamma.check))
    # check CV
    round(c(mu.proc,CV.proc),3)
    quantile(sqrt(gamma.check),c(0.1,0.9))
    }  
    }else{
    sigma.proc = 0.1 #IF Fixed (sigma.est = FALSE): typicallly 0.05-0.15 (see Ono et al. 2012)
  
    }
  
  #--------------------------------------------
  
  
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Prior specification for Model 5: JABBA-SELECT
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # The following section is ignored if: Model < 5
    
    # Option to set SBmsy reference point (e.g. SBmsy/SB0 = 0.4)
    # If SBmsy_SB0 = NULL the reference SBmsy_SB0 that produces MSY will be used 
    SBmsy_SB0 = NULL   # Standard Reference setting for SA Linefish assessment reference points
    
  #---------------------------------------------------------------
  # STOCK PARAMETERS for prior generation Hmsy as a fuction of r
  #---------------------------------------------------------------
    
    minage <- 0  																						
    maxage <- ss4js$stock.pars$Amax
    PlusGroup = TRUE
    # Number of sexes order Female and Males (this model is sex-structured)
    nsexes = 2
    # VBGF parameters
    Linf <- c(ss4js$stock.pars$vbgf_F[1],ss4js$stock.pars$vbgf_M[1])
    kappa <- c(ss4js$stock.pars$vbgf_F[2],ss4js$stock.pars$vbgf_M[2])
    t0 <- c(ss4js$stock.pars$vbgf_F[3],ss4js$stock.pars$vbgf_M[3])
    
    #Length-weight
    aW <- exp(c(ss4js$stock.pars$LW_F[1],ss4js$stock.pars$LW_F[1])) 																						
    bW <- c(ss4js$stock.pars$LW_F[2],ss4js$stock.pars$LW_F[2])
    
    # Maturity 
    maturity = c(ss4js$stock.pars$mat,0)  #c(Lm50,Lm95,L=0/age=1) # two values are taken as length-based: Lm50 and Lm95 
    # Reduce maturity to approximate mean length-at-first capture
    if(s==3) maturity = c(100,ss4js$stock.pars$mat[2],0)  #c(Lm50,Lm95,L=0/age=1) # two values are taken as length-based: Lm50 and Lm95 
    
    # Natural mortality estimate (affects Hmsy)
    M = 0.2
    CV.M = 0.2 # ss3  M "fixed"
    
    # steepness B&H SSR (effects Hmsy and determines SBmsy/SB0)
    h = 0.8
    CV.h = 0.1
    
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Optional: Do TAC Projections
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  Projection = TRUE # Switch on by Projection = TRUE 
  
  # Check final year catch
  apply(catch[,-1],1,sum,na.rm=TRUE)

  # Set range for alternative TAC projections
  TACs = seq(9000,12500,500) #example
  
  # Set year of first TAC implementation
  imp.yr = 2018
  # Set number of projections years
  pyrs = 10
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Execute model and produce output
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # MCMC settings
  ni <- 25000 # Number of iterations
  nt <- 2 # Steps saved
  nb <- 3000 # Burn-in
  nc <- 1 # number of chains
  nsaved = (ni-nb)/nt*nc
  cat(paste0("\n","- Run Model","\n"))
  source(paste0(JABBA.file,"/JABBA_SELECT",version,".r"))
}

#-----------------------------------------------------
# Create FRL type object to compare JABBA-Select model
#-----------------------------------------------------

Bs = Fs = BBmsys = FFmsys = BKs = SPs =Bi =NULL
mods = list(timeseries = NULL,refpts=NULL,pfunc=NULL,curves=NULL,dgs=NULL)
for(s in 1:length(Scenarios)){
  Scenario = Scenarios[s]
  Mod.names = c("JS")
  output.dir = paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Output")
  load(paste0(output.dir,"/jabbaout_",assessment,"_",Scenario,".RData"),verbose=T)
  load(paste0(output.dir,"/jb.",Scenario,".RData"),verbose=T)
  
  mods$timeseries = rbind(mods$timeseries,jb$timeseries) 
  mods$refpts = rbind(mods$refpts,jb$refpts) 
  mods$pfunc = rbind(mods$pfunc,jb$pfunc) 
  mods$curves = rbind(mods$curves,jb$curves) 
  jb$dgs$xval = ifelse(jb$dgs$tails<jb$dgs$year,TRUE,FALSE)
  mods$dgs = rbind(mods$dgs,jb$dgs) 
  Bs = cbind(Bs,ifelse(Mod.names==rep("JS",nrow(jabba_out$Stock.trj)),jabba_out$Stock.trj$SBt.Median,jabba_out$Stock.trj$Bt.Median)) 
  
  Fs = cbind(Fs,jabba_out$Stock.trj$Ft.Median) 
  BBmsys = cbind(BBmsys,ifelse(Mod.names==rep("JS",nrow(jabba_out$Stock.trj)),jabba_out$Stock.trj$SBt_SBmsy.Median,jabba_out$Stock.trj$Bt_Bmsy.Median)) 
  FFmsys = cbind(FFmsys,jabba_out$Stock.trj$Ft_Fmsy.Median) 
  BKs = cbind(BKs,ifelse(Mod.names==rep("JS",nrow(jabba_out$Stock.trj)),jabba_out$Stock.trj$SBt_SB0.Median,jabba_out$Stock.trj$Bt_K.Median)) 
  SPs =cbind(SPs,jabba_out$SPphase$SP)
  Bi = cbind(Bi,jabba_out$SPphase$SB_i)
}
save(mods,file=paste0(File,"/",assessment,"/BaseDiags.rdata"))

Labels = Scenarios 
# 3x2 
Par = list(mfrow=c(3,2),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/Compare_",assessment,".png"), width = 6, height = 8, 
    res = 200, units = "in")
par(Par)
#><> plot B
years = order(unique(mods$timeseries$year))
y = Bs
plot(years,years,type="n",ylim=c(0,max(y)),ylab="Biomass (t)",xlab="Year")
for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=i,lwd=2,lty=1)
}
legend("bottomleft",paste(Labels),col=1:length(Scenarios),bty="n",cex=0.75,pt.cex=0.7,lwd=c(2,rep(2,length(Labels))))

y = Fs
plot(years,years,type="n",ylim=c(0,max(y)),ylab="Fishing Mortality F",xlab="Year")
for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=i,lwd=2,lty=1)
}

y = BKs
plot(years,years,type="n",ylim=c(0,1),ylab="B/K",xlab="Year")
for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=i,lwd=2,lty=1)
}
abline(h=0.4,lty=2)

# Plot SP
plot(years,years,type="n",ylim=c(0,max(SPs)),xlim=c(0,max(Bi)),ylab="Surplus Production (t)",xlab="Biomass (t)")
for(i in 1:length(Scenarios)){
  lines(Bi[,i],SPs[,i],col=i,lwd=1,lty=5)
  points(mean(Bi[SPs[,i]==max(SPs[,i]),i]),max(SPs[,i]),col=i,pch=16,cex=1.2)
}

# Plot BtoBmsy
y = BBmsys
plot(years,years,type="n",ylim=c(0,max(y)),ylab=expression(paste(B/B[MSY])),xlab="Year")
for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=i,lwd=2,lty=1)
}
abline(h=1,lty=2)

#><> Plot FtoFmsy
y = FFmsys
plot(years,years,type="n",ylim=c(0,max(y)),ylab=expression(paste(F/F[MSY])),xlab="Year")
for(i in 1:length(Scenarios)){
  lines(years,y[,i],col=i,lwd=2,lty=1)
}
abline(h=1,lty=2)
dev.off()


