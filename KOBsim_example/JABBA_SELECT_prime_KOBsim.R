##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
## Prime File for KOBsim JABBA-SELECT example
## written by Henning Winker
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

# Delete all objects
rm(list=ls())
gc()
# Install required packages if missing
list.of.packages <- c("gplots", "coda","rjags","R2jags","fitdistrplus","reshape","mvtnorm","scales")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load Packages
library(gplots);library(coda);library(rjags);library(R2jags);library("fitdistrplus");library(reshape);
library(mvtnorm);library(scales)

#---------------------------------------------------------------------
# Set Working directory file where to store the results
#---------------------------------------------------------------------
File = "C:/Work/Research/GitHub/JABBA-SELECT/KOBsim_example"
# Set working directory for JABBA-Select R source code
JABBA.file = "C:/Work/Research/GitHub/JABBA-SELECT"
# JABBA-Select Version
version = "v1.1"
# Set Assessment file: assement folder within File that includes .csv input files
assessment = "KOBSim"


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Graphic, Output, Saving (.RData) settings 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
SELECT = TRUE
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

Scenarios = c("SELECT","Pella")
s=1

for(s in 1:length(Scenarios)){
  Scenario = Scenarios[s] 
  if(s==2) SELECT = FALSE # Reduce JABBA-Select to a Bayesian Pella-Tomlinson
  
  
  #--------------------------------------------------
  # Read csv files
  #--------------------------------------------------
  # Load assessment data
  catch = read.csv(paste0(File,"/",assessment,"/catch",assessment,".csv"))
  cpue = read.csv(paste0(File,"/",assessment,"/cpue",assessment,".csv"))#
  
  # Use SEs from csv file for abudance indices (TRUE/FALSE)
  SE.I = TRUE
  if(SE.I ==TRUE){
    se = read.csv(paste0(File,"/",assessment,"/se",assessment,".csv"))
  }
  names(cpue)
  names(catch)
 
  # Read select csv
  select = read.csv(paste0(File,"/",assessment,"/select",assessment,".csv"))
  
  # NOTE: 
  ## only unique SL50 values (no replicates) 
  ## Selectivity SL50 must be sufficiently different (+-5%) between "fleets" to seperate r 
  
  
  
  # Read selex (selectivity) csv
  selex = read.csv(paste0(File,"/",assessment,"/selex",assessment,".csv"))
  

  #---------------------------------------
  # option to exclude CPUE time series 
  #---------------------------------------
  # Not used here
  
  #------------------------------------------------------
  # Option use mean CPUE from state-space cpue averaging
  #-----------------------------------------------------
  
  CPUE.plot = TRUE
    
  #------------------------------------------------
  # mean and CV and sd for unfished biomass K (SB0)
  #------------------------------------------------
  mu.SB0 = 30000; CV.SB0 = 2; sd.SB0=sqrt(log(CV.SB0^2+1)) 
  SB0.pr = c(mu.SB0,CV.SB0) # Use CV as input
  
  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  # Set the initial depletion prior SB/SB0 
  # To be converted into a lognormal prior (with upper bound at 1.1)
    
    psi.prior= c("lognormal","beta")[1]
    # specify as mean and CV 
    mu.psi = 1 
    CV.psi = 0.3 
    
 
    #--------------------------------------------
    
    
    # Soft Penalty to control lower and upper SB_t/SB0
    P_bound = c(0.0001,1.5) # default if not specified
    

    
    #---------------------------------------------------------------
    # STOCK PARAMETERS for prior generation Hmsy as a function of r
    #---------------------------------------------------------------
    # Age
    minage <- 0  																						
    maxage <- 20
    plusgroup = c(FALSE,TRUE)[1] 
    #Growth paramters (von Bertalnffy)
    Linf <- 1372
    kappa <- 0.115
    t0 <-  -0.815
    
    #Length-weight
    aW <- 0.000006 																						
    bW <- 3.07
    
    # Maturity 
    maturity = c(2.5,2.6,1) # knife edge age-3 # a single value is taken as age (knife-edge) 
    #c(Lm50,Lm95,0) # two values are taken as length-based: Lm50 and Lm95 
    
    # Natural mortality estimate (affects Hmsy)
    M = 0.18
    CV.M = 0.25 
    
    # steepness B&H SSR (effects Hmsy and determines SBmsy/SB0)
    h = 0.8
    CV.h = 0.15
    
   
    #--------------------------------------------------------------
    # Determine estimation for catchability q and observation error 
    #--------------------------------------------------------------
    
    # Assign q to CPUE from select file
    sets.q = select$q[select$CPUE]
    if(SELECT==F) sets.q = 1:length(select$q[select$CPUE])
    
    # Assign additional observation variances to indices
    sets.var = select$q[select$CPUE] # assigns indices with same obervation variance 
    
    #To Estimate additional observation variance set sigma.est = TRUE
    sigma.est = c(TRUE,FALSE)[1]
    # As option for data-weighing
    # minimum additional observation error for each variance set (optional choose 1 value for both)
    fixed.obsE = 0.01
    
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
    # Process Error
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
    
    #------------------------------------------
    #Estimate set sigma.proc == True
    sigma.proc = TRUE
    proc.dev.all=1 # start process error in year 1
    #------------------------------------------
    if(sigma.proc == TRUE){
      proc.type = c("igamma","lnorm")[1] # choose 1: inverse-gamma or 2: lognormal
      
      if(proc.type=="lnorm"){
        pr.proc = c(log(0.08),0.2) # Option for lognormal process error
      }
      if(proc.type=="igamma"){
        #pr.proc = c(0.001,0.001) # Option for inverse-gamma prior
        #pr.proc = c(10,0.1)
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
   
    
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
    # Option to set Bmsy reference point (e.g. SBmsy/SB0 = 0.4)
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
    
    
    # If SBmsy_SB0 = NULL the model implicite estimate (SBmsy_SB0 where MSY) will be used 
    SBmsy_SB0 = 0.4  # Standard Reference setting for SA Linefish assessment reference points
   
    
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
    # Optional: Do TAC Projections
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
    Projection = TRUE # Switch on by Projection = TRUE 
    
    # Check final year catch
    curC =  apply(catch[,-1],1,sum,na.rm=TRUE)[nrow(catch)]
    
    # Set range for alternative TAC projections
    TACs = ceiling(seq(0.5,1.2,0.1)*curC) #example
    
    # Set year of first TAC implementation
    imp.yr = 41
    
    # Set number of projections years
    pyrs = 10
    
    
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
    # Execute model and produce output
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
    
    # MCMC settings
    ni <- 14000 # Number of iterations
    nt <- 2 # Steps saved
    nb <- 4000 # Burn-in
    nc <- 3 # number of chains
    nsaved = (ni-nb)/nt*nc
    cat(paste0("\n","- Run Model","\n"))
    
    dir.create(paste0(File,"/",assessment),showWarnings = F)
    
    # RUN JABBA-SELECT
    source(paste0(JABBA.file,"/JABBA_SELECT",version,".r"))
    

#-------------------------------------------       
#    ADDED CODE TO 
#    Compare with True Simulated Scenarios
#------------------------------------------    
if("ThorsonUtilities" %in% installed.packages()[,"Package"]==FALSE)install_github("James-Thorson/utilities")
library(ThorsonUtilities)
load(paste0(File,"/DataList_Original.rdata"),verbose=T)

ThorsonUtilities::save_fig( file=paste0(File,"/",Scenarios[s]), width=4, height=6 )
par( mfrow=c(2,1), yaxs="i", xaxs="i", mgp=c(2,0.5,0), tck=-0.02, mar=c(2,4,1,0)  )
# Depletion
mu.y = apply(posteriors$P,2,quantile,c(0.025,0.5,0.975))
ylim = c(0, max(mu.y))
cord.x <- c(years,rev(years))
cord.y <- c(mu.y[1,],rev(mu.y[3,]))
plot(years, mu.y[2,], type="n", col=c(1), lty="solid",lwd=1, xlab="Year", ylab="SB_t/SB0", ylim=ylim) 
polygon(cord.x,cord.y,col='grey',border=0,lty=2,grey(0.5,0.2))
lines(DataList_Original$SB_t/DataList_Original$SB0,col=2) 
lines(mu.y[2,],col=1) 
legend("topright", bty="n", fill=c("black","red"), legend=c("Estimated","True") )
# SB_t
mu.y = apply(posteriors$SB,2,quantile,c(0.025,0.5,0.975))
ylim = c(0, max(mu.y))
cord.x <- c(years,rev(years))
cord.y <- c(mu.y[1,],rev(mu.y[3,]))
plot(years, mu.y[2,], type="n", col=c(1), lty="solid",lwd=1, xlab="Year", ylab="SB_t", ylim=ylim) 
polygon(cord.x,cord.y,col='grey',border=0,lty=2,grey(0.5,0.2))
lines(DataList_Original$SB_t,col=2) 
lines(mu.y[2,],col=1) 
legend("topright", bty="n", fill=c("black","red"), legend=c("Estimated","True") )
dev.off()
}
    
    
    
