##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><
## JABBA-SELECT: JABBA extension to integrate Life History and Selectivity Distortion
## #Executeable to JABBA-SELECT 
## written by Henning Winker
## henning.winker@gmail.com
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
Mod.names = "JS" 
cat(paste0("\n","- Run Model ",Mod.names,"\n"))
#setwd(paste(File))
dir.create(paste0(File,"/",assessment,"/",Scenario,"_",Mod.names),showWarnings = F)
dir.create(paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Input"),showWarnings = F)
input.dir = paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Input")


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Define objects to make sure they exist
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
if(exists("proc.type")==FALSE) proc.type = "igamma"  # Produces JABBA Kobe plot 
if(exists("pr.proc")==FALSE) igamma = c(4,0.01)  # Produces JABBA Kobe plot 
if(exists("Plim")==FALSE) Plim = 0  # Produces JABBA Kobe plot 
if(exists("KOBE.plot")==FALSE) KOBE.plot = TRUE # Produces JABBA Kobe plot 
if(exists("KOBE.type")==FALSE) KOBE.type = c("ICCAT","IOTC")[2] # ICCAT uses 3 colors; IOTC 4 (incl. orange) 
if(exists("SP.plot")==FALSE) SP.plot = c("standard","phase")[2] # Produces standard or 'Kobe phase' SP plot  
if(exists("Biplot")==FALSE) Biplot= TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
if(exists("save.trajectories")==FALSE) save.trajectories =FALSE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
if(exists("catch.metric")==FALSE) catch.metric = "(t)" # Runs state-tool to produce "alligned" multi-CPUE plot  
if(exists("harvest.label")==FALSE) harvest.label  = c("Hmsy","Fmsy")[1] # choose label preference H/Hmsy versus Fmsy
if(exists("CPUE.plot")==FALSE) CPUE.plot= TRUE # Runs state-tool to produce "alligned" multi-CPUE plot  
if(exists("meanCPUE")==FALSE) meanCPUE = TRUE # Uses averaged CPUE from state-space tool instead of individual indices  
if(exists("Projection")==FALSE) Projection = FALSE # Use Projections: requires to define TACs vectors 
if(exists("TACint")==FALSE) TACint = mean(apply(catch[(nrow(catch)-2):nrow(catch),-1],1,sum,na.rm=TRUE)) # use mean catch from last years
if(exists("save.projections")==FALSE) save.projections = FALSE# saves projection posteriors as .RData object 
if(exists("imp.yr")==FALSE) imp.yr = max(catch[,1])+1
if(exists("Reproduce.seed")==FALSE) Reproduce.seed = FALSE # If FALSE a random seed assigned to each run (default)
if(exists("P_bound")==FALSE) P_bound = c(0.02,1.3)  # Soft penalty bounds for P 
if(exists("q_bounds")==FALSE) q_bounds= c(10^-30,1000) # Defines lower and upper bounds for q 
if(exists("sigmaobs_bound")==FALSE) sigmaobs_bound = 1 # Adds an upper bound to the observation variance  
if(exists("sigmaproc_bound")==FALSE) sigmaproc_bound = 0.2 # Adds an upper bound to the process variance  
if(exists("SB0_bounds")==FALSE) SB0_bounds= c(0.01,10^10) # Defines lower and upper bounds for q 
if(exists("runASEM")==FALSE) runASEM = TRUE # Run Monte-Carlo ASEM
if(exists("PlusGroup")==FALSE) PlusGroup = FALSE # if TRUE add PlusGroup to ASEM
if(exists("init.values")==FALSE) init.values = FALSE # Option to fix init values of q's and SB0
if(exists("nsexes")==FALSE) nsexes = 1
if(exists("SELECT")==FALSE) SELECT = TRUE
if(exists("save.all")==FALSE) save.all = FALSE #  
if(exists("SBmsy_SB0")==FALSE) SBmsy_SB0 = NULL
if(exists("proc.dev.all")==FALSE) proc.dev.all=TRUE
refB = ifelse(is.null(SBmsy_SB0)==FALSE,paste(SBmsy_SB0*100),"MSY") 
if(exists("p")==FALSE) p=0
if(exists("GFUN")==FALSE) GFUN = 1 # Plan to implement alternative growth function
if(exists("jabba2FRL")==FALSE) jabba2FRL = FALSE
if(sigma.proc==FALSE){
pr.proc=c(2,2)  
}

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


#-------------------------
# Prepare input data
#-------------------------
indices = names(cpue)[2:ncol(cpue)]
n.indices = max(length(indices),1)
catches = names(catch)[2:ncol(catch)]
n.catches = length(catches)

years=catch[,1]
styr = min(years)
endyr = max(years)
n.years = length(years)
styr.cpue = min(cpue[,1])
styr.I = styr.cpue-styr+1 


# Convert input data to matrices
conv.cpue = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
CPUE=matrix(conv.cpue,nrow=n.years,ncol=n.indices)

if(SE.I==FALSE){
  se = cpue  
  conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
  se2 = matrix(ifelse(fixed.obsE>0,fixed.obsE^2,10^-10),n.years,n.indices)#/2
} else{
  conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(se[,-1])))
  #conv.se = sqrt(conv.se^2+fixed.obsE^2) 
  se2 = matrix(ifelse(is.na(conv.se),0.3^2,conv.se)^2,n.years,n.indices)+fixed.obsE^2#/2
}

conv.catch = 0.00001+as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.catches),styr.I-1,n.catches),as.matrix(catch[,-1])))
Catch=matrix(conv.catch,nrow=n.years,ncol=n.catches)
Catch[is.na(Catch)] = 0.00001 # Replace any NA by zero 

# Total Catch
TC = apply(Catch,1,sum)

# hindcast option
if(exists("tails")==FALSE) tails = max(years)
if(tails<max(years)){
cpue.raw = read.csv(paste0(File,"/",assessment,"/cpue",assessment,".csv"))  
cpue.raw = cpue.raw[which(cpue.raw$Yr%in%years),which(colnames(cpue.raw)%in%colnames(cpue))]   
} else {
  cpue.raw = cpue
}

#---------------------
# Index color palette
#---------------------
jabba.colors = as.character(c("#e6194b", "#0082c8","#3cb44b",
                                 "#f58231", "#911eb4",
                                  "#46f0f0", "#f032e6", "#d2f53c",
                                  "#fabebe", "#008080","#e6beff", "#aa6e28","#ffe119",rainbow(12)[seq(1,12,3)],rainbow(12)[seq(2,12,3)],rainbow(12)[seq(3,12,3)]))



#####################################################################################

# Plot Catch

cat(paste0("\n","- Plot Catch in Input subfolder","\n"))


Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0), tck = -0.02,cex=0.8)
png(file = paste0(input.dir,"/Catches_",assessment,".png"), width = 7, height = 5, 
    res = 200, units = "in")
par(Par)
catch[is.na(catch)]=0
xpol = c(years,rev(years))
plot(years,TC,ylim=c(0,max(TC)),ylab=paste0("Catch ",catch.metric),xlab="Year",type="n")
Cacc = rep(0,length(catch[,1]))
for(i in 2:ncol(catch)){
polygon(xpol,c(Cacc,rev(Cacc+catch[,i])),col=jabba.colors[i])
Cacc=Cacc+catch[,i]
lines(catch[,1],Cacc,lty=(1),lwd=1,col=1)
}
legend("topleft",paste(names(catch)[2:ncol(catch)]),pt.cex=1.5,pch=22,pt.bg=jabba.colors[-1],bty="n",cex=1)
dev.off()

#--------------------------
# Capture Stock Parameters
#--------------------------
stockpars = list(minage=minage,maxage=maxage,nsexes=nsexes,PlusGroup=PlusGroup,
                 Growth= paste0("VBGF (Linf, k, t0)"), vb.female= c(Linf[1],kappa[1],t0[1]),vb.male= c(Linf[1],kappa[1],t0[1]),
                 LW = paste0("Length-Weight relationship (a,b)"),
                 lw.female = as.numeric(c(aW[1],bW[1])),
                 lw.male = as.numeric(c(aW[2],bW[2])),
                 Maturation= paste0("Maturity@", ifelse(maturity[3]==1,"Age","Length"),"(50%, 95%)"),
                 maturity = maturity[1:2],
                 NatMortality = "Natural Mortality (M): Mean, CV",
                 M = c(M,CV.M),
                 SSR = "B&H steepness (h): mean, CV",
                 h = c(h,CV.h)
                 )


#------------------------------------------------
# Selectivity will determine changes in r (Fmsy)
#------------------------------------------------
# Selectivity SL50 must be sufficiently different (+-5%) to seperate r 
# only unique SL50 values (no replicates) 

# Selectivity SL50 must be sufficiently different (+-5%) between "fleets" to seperate r 
SL50 <- as.numeric(selex[1,-1]) 
SL95 <- as.numeric(selex[2,-1])  # If unknown set to 0.05*SL50 ~ knife-edge

# Define point where descening limb starts (set Linf for logistic)
SL.desc <-as.numeric(selex[3,-1])  # mean of half-normal 
# Define rate of decreasing selectivity
CV.desc <- as.numeric(selex[4,-1]) # CV of half-normal 
# Define minimum descending limp between 0 and 1
min.desc =as.numeric(selex[5,-1]) 

# number of different Hmsy (r) priors 
nSel = length(SL50)

# Translate NAs 
SL.desc[is.na(SL.desc)] =rep(Linf,nSel)
CV.desc[is.na(CV.desc)] =rep(0.1,nSel)
min.desc[is.na(min.desc)] =rep(0.001,nSel)

#---------------------------------------------------------------
# Compile fleet/index settings from selex and select input files
#--------------------------------------------------------------- 

# Assign Selectivity to abundance indices
sets.I = select$Selectivity[select$CPUE] 

# Assign Selectivity to catch series
sets.C = select$Selectivity[select$Catch] # here 1: South, 2: South-East, 3: Trawl

# Define if index is in numbers: 0 or biomass: 1 
I.unit =  select$CPUE.units[select$CPUE]

# Reduce JABBA-Select to Pella-Tomlison Model
if(SELECT==FALSE){nSel = 1; nsexes=1; sets.C = rep(1,length(sets.C));sets.I = rep(1,length(sets.I));I.unit = rep(1,length(I.unit))}

#--------------------
# Set seed
#--------------------
if(Reproduce.seed==FALSE){
  get_seed = ceiling(runif(1,min=0,max=1e6)) } else {get_seed = 123}
set.seed(get_seed)  

#---------------------------------------------------------------------------
# CPUE run State-Space model for averaging CPUE
#---------------------------------------------------------------------------
if(CPUE.plot==TRUE){ 
  cat(paste0("\n","><> Run State-Space CPUE averaging tool","\n"))
  #find first time-series with first CPUE
  q1.y = c(1:n.years)[is.na(apply(CPUE,1,mean,na.rm=TRUE))==FALSE][1] #first year with CPUE
  q1.I = which.max(CPUE[q1.y,])
  
  qs = c(q1.I,c(1:(ncol(cpue)-1))[-q1.I])
  
  
  sink("cpueAVG.jags")
  cat("
      model {
      
      # Prior specifications  
      eps <- 0.0000000000001 # small constant    
      
      iq[1] ~ dgamma(1000,1000)
      q[1] <-  pow(iq[1],-1)
      logq[1] <- log(1)
      for(i in 2:nI){
      iq[i] ~ dgamma(0.001,0.001)
      q[i] <- pow(iq[i],-1)
      logq[i] <-  log(q[i])
      }
      
      
      ")
  
  if(sigma.proc==TRUE){
    cat("
        # Process variance
        isigma2 <- isigma2.est 
        sigma2 <- pow(isigma2,-1)
        sigma <- sqrt(sigma2)
        fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
        ",append=TRUE)  
  }else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
             sigma2 <- pow(isigma2,-1)
             sigma <- sqrt(sigma2)
             
             ",append=TRUE)}
  
  if(sigma.est==TRUE){
    cat("
        # Obsevation variance
        # Observation error
        itau2~ dgamma(0.001,0.001)
        tau2 <- 1/itau2
        
        
        for(i in 1:nI)
        {
        for(t in 1:N)
        {
        var.obs[t,i] <- SE2[t,i]+tau2
        ivar.obs[t,i] <- 1/var.obs[t,i]
        # note total observation error (TOE)     
        TOE[t,i] <- sqrt(var.obs[t,i])
        
        }}
        ",append=TRUE)  
  }else{ cat(" 
      # Obsevation variance
             # Observation error
             itau2~ dgamma(2,2)
             tau2 <- 1/itau2
             
             
             for(i in 1:nI)
             {
             for(t in 1:N)
             {
             var.obs[t,i] <- SE2[t,i] # drop tau2
             fake.tau[t,i] <- tau2
             
             ivar.obs[t,i] <- 1/var.obs[t,i]
             # note total observation error (TOE)     
             TOE[t,i] <- sqrt(var.obs[t,i])
             
             }}
             
             ",append=TRUE)}
  
  # Run rest of code  
  cat("  
      # Process variance prior
      isigma2.est ~ dgamma(0.001,0.001)
      
      
      # Priors and constraints
      logY.est[1] ~ dnorm(logY1, 1)       # Prior for initial population size
      
      mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
      
      # Likelihood
      # State process
      for (t in 1:(N-1)){
      r[t] ~ dnorm(mean.r, isigma2)
      logY.est[t+1] <- logY.est[t] + r[t] }
      
      # Observation process
      for (t in 1:N) {
      for(i in 1:nI){
      y[t,i] ~ dnorm(logY.est[t]+logq[i], ivar.obs[t,i])
      }}
      
      # Population sizes on real scale
      for (t in 1:N) {
      Y.est[t] <- exp(logY.est[t])
      }
      
  } 
      ",fill = TRUE)
  sink()
  
  
  
  q.init = 1
  mCPUE = as.matrix(CPUE[q1.y:n.years,qs])
  mSE2 = as.matrix(se2[q1.y:n.years,qs])
  if(n.indices>1) for(i in 2:n.indices){q.init[i] = mean(mCPUE[,i],na.rm=TRUE)/mean(mCPUE[,1],na.rm=TRUE)}
  # Bundle data
  jags.data <- list(y = log(mCPUE),SE2=mSE2, logY1 = log(mCPUE[1,1]), N = length(q1.y:n.years),nI=n.indices,sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc))
  
  # Initial values
  inits <- function(){list(isigma2.est=runif(1,20,100), itau2=runif(1,80,200), mean.r = rnorm(1),iq = 1/q.init)}
  
  # Parameters monitored
  parameters <- c("mean.r", "sigma","r", "Y.est","q")
  
  
  # Call JAGS from R (BRT 3 min)
  mod.cpue <- jags(jags.data, inits, parameters, "cpueAVG.jags", n.chains = nc, n.thin = max(nt,2), n.iter = max(ni/5,10000), n.burnin = nb/10)
  
  
  cat(paste0("\n","><> Plot State-Space CPUE fits  in Input subfolder <><","\n"))
  # get individual trends
  fitted <- lower <- upper <- NULL
  cpue.yrs = years[q1.y:n.years]
  
  for (t in 1:nrow(mCPUE)){
    fitted[t] <- median(mod.cpue$BUGSoutput$sims.list$Y.est[,t])
    lower[t] <- quantile(mod.cpue$BUGSoutput$sims.list$Y.est[,t], 0.025)
    upper[t] <- quantile(mod.cpue$BUGSoutput$sims.list$Y.est[,t], 0.975)}
  
  
  q.adj = apply(mod.cpue$BUGSoutput$sims.list$q,2,median)
  
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(input.dir,"/CPUE_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  u.ylim = NULL
  for(i in 1:n.indices){ u.ylim = c(u.ylim,exp(log(mCPUE[,i]/q.adj[i])+1.96*sqrt(mSE2[,i])))}  
  ylim = c(0,max(u.ylim,na.rm=TRUE))
  plot(0, 0, ylim = ylim, xlim = range(cpue.yrs), ylab = "Expected CPUE", xlab = "Year", col = "black", type = "n")
  legend("topright",paste(indices),lwd=2,col=(jabba.colors)[1:n.indices],bty="n")
  polygon(x = c(cpue.yrs,rev(cpue.yrs)), y = c(lower,rev(upper)), col = "gray", border = "gray90")
  
  for(i in 1:n.indices)
  {
    shift = runif(1,-0.1,0.1)
    cols=jabba.colors[qs[i]]
    plotCI(cpue.yrs+shift,mCPUE[,i]/q.adj[i],ui=exp(log(mCPUE[,i]/q.adj[i])+1.96*sqrt(mSE2[,i])),li=exp(log(mCPUE[,i]/q.adj[i])-1.96*sqrt(mSE2[,i])),add=TRUE,col= cols,pt.bg = cols,pch=21,gap=0)
    lines(cpue.yrs+shift,mCPUE[,i]/q.adj[i], col = cols,lwd=2)
    points(cpue.yrs+shift,mCPUE[,i]/q.adj[i], bg = cols,pch=21)
  }
  lines(cpue.yrs,fitted,lwd=2)
  
  dev.off()
  
  logSE = apply(log(mod.cpue$BUGSoutput$sims.list$Y.est),2,sd)
  
  
  if(nrow(mCPUE)<n.years) {
    fitted = c(rep(NA,q1.y-1),fitted)
    logSE = c(rep(0.2,q1.y-1),logSE)
  }    
  avgCPUE = data.frame(Year=years,CPUE= fitted,logSE=logSE)
  
  write.csv(avgCPUE,paste0(input.dir,"/avgCPUE_",assessment,"_",Scenario,".csv"))
  
  if(meanCPUE==TRUE){
    cat(paste0("\n","><> Use average CPUE as input for JABBA <><","\n"))
    
    CPUE = as.matrix(avgCPUE[,2]) 
    cpue.check = cpue[,-1]
    cpue.check[is.na(cpue[,-1])]=0
    CPUE[,1] = ifelse(apply(cpue.check,1,sum)==0,rep(NA,length(CPUE[,1])),CPUE[,1])
    se2 =  as.matrix(avgCPUE[,3]^2)     
    n.indices=1
    indices = "All"
    sets.q =1
    sets.var =1
  }
  
  }


#--------------------------------------------------------------------------
# END of CPUE State-Space tool
#--------------------------------------------------------------------------



#################################################################
# FUNCTIONS
#################################################################


#--------------------------------------------------
# Function to get beta prior parameters
#--------------------------------------------------
get_beta <- function(mu,CV,Min=0,Prior="x"){
  a = seq(0.0001,5000,0.001)
  b= (a-mu*a)/mu
  s2 = a*b/((a+b)^2*(a+b+1))
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter 
  b = (a-mu*a)/mu
  x = seq(Min,1,0.001)  
  pdf = dbeta(x,a,b)  
  plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n",ylim=c(0,max(pdf*1.2)))
  polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  return(c(a,b))
}

get_gamma <- function(mu,CV,Prior="x"){
  a = seq(0.0001,5000,0.0001)
  b = a/mu
  s2 = (a/b^2)
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter b
  b = a/mu
  x = sort(rgamma(1000,a,b))  
  pdf = dgamma(x,a,b)  
  plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
  polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  return(c(a,b))
}


# Bias corrected lognormal
plot_lnorm <- function(mu,CV,Prior="x"){
  sdev= sqrt(log(CV^2+1))
  rand.pr = rlnorm(1000,log(mu)-0.5*sdev^2,sdev)
  x = seq(min(rand.pr),quantile(rand.pr,0.995),max(rand.pr/500))  
  pdf = dlnorm(x,log(mu)-0.5*sdev^2,sdev)  
  plot(x,pdf,type="l",xlim=range(x),xlab=paste(Prior),ylab="",yaxt="n")
  polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  return(c(exp(log(mu)-0.5*sdev^2),sdev))
}



#------------------------------------
# Function kobeJabba for FLR
#------------------------------------
kobeJabba<-function(x,minyear=1){
  
  out=cbind(melt(x[,,2]),c(x[,,3]))
  names(out)=c("iter","year","stock","harvest")
  out$year=out$year+minyear-1
  out}

#-------------------------------------------------
# Function kobeJabbaProj for projections with FLR
#-------------------------------------------------
kobeJabbaProj<-function(x,minyear=1,tac=NULL){
  
  out=cbind(melt(x[,,,2]),c(x[,,,3]))
  names(out)=c("iter","year","tac","stock","harvest")
  out$year=out$year+minyear-1
  
  out}

#---------------------------------------------------------------------------
# Equilibrium Age-structured Function:
# Computes Spawner Biomass (SB), Yield and Exploitable Biomass (EB) 
# with Beverton and Holt spawner-recruitment relationship (S-R) 
#---------------------------------------------------------------------------

# Input paramters
# Growth: Linf, K, t0 (Female,Male)
# Length-weight: aW (Female), bW (Male)
# Maturity: Lm50, dM (logistic)
# Longivity: tmax
# Harvest rate: H (Catch/Exploitable Biomass)
# Selectivity: SL50, dS (logistic)
# steepness: z
# Natural Mortality: M
# Harvest rate: H (Catch/Exploitable Biomass)


get_ASEM <- function(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars,h,M,F_i,R0,PlusGroup=FALSE,nsexes=1,SELECT=TRUE){
  age = minage:maxage
  nages = length(age)
  ntilda_f = mat.or.vec(length(F_i),nages)  
  ntilda_m = mat.or.vec(length(F_i),nages)  
  
  
  #Length-at-age
  L_f = Linf[1]*(1-exp(-kappa[1]*(age-t0[1])))
  L_m = Linf[nsexes]*(1-exp(-kappa[nsexes]*(age-t0[nsexes])))
  
  if(t0[1]>=0) L_f[1] = 0.1    
  if(t0[nsexes]>=0) L_f[2] = 0.1    
    
  # Weight-at-age
  W_f = aW[1]*L_f^bW[1] 
  W_m = aW[nsexes]*L_f^bW[nsexes] 
  
  
  # Maturity-at-age
  if(maturity[3]==1){
    mat = 1/(1+exp(-log(19)*(age-maturity[1])/(maturity[2]-maturity[1])))
    } else {
    mat = 1/(1+exp(-log(19)*(L_f-maturity[1])/(maturity[2]-maturity[1])))
    }
  
  #Selectivity-at-age (logistic/Domeshaped)
  if(SELECT==TRUE){
  sel_a_f = 1/(1+exp(-log(19)*(L_f-selpars[1])/(selpars[2]-selpars[1])))
  sel_b_f = dnorm(L_f,selpars[3],selpars[3]*selpars[4])/max(dnorm(L_f,selpars[3],selpars[3]*selpars[4]))
  sel_c_f = 1+(selpars[5]-1)*(sel_b_f-1)/(-1)
  sel_f = as.numeric(ifelse(rep(selpars[3],length(L_f))>=L_f,sel_a_f,sel_c_f))
  sel_a_m = 1/(1+exp(-log(19)*(L_m-selpars[1])/(selpars[2]-selpars[1])))
  sel_b_m = dnorm(L_m,selpars[3],selpars[3]*selpars[4])/max(dnorm(L_m,selpars[3],selpars[3]*selpars[4]))
  sel_c_m = 1+(selpars[5]-1)*(sel_b_m-1)/(-1)
  sel_m = as.numeric(ifelse(rep(selpars[3],length(L_m))>=L_m,sel_a_m,sel_c_m))
  } else {
    sel_f = sel_m = mat 
  }
  
  # Length-based for plotting  
  Li = seq(0,max(L_m,L_f),1)
  LW=cbind(aW[1]*Li^bW[1],aW[nsexes]*Li^bW[nsexes])
  
  if(maturity[3]==1){
    aLi = t0[1]-1/kappa[1]*log(1-Li/Linf[1])
    Lmat = 1/(1+exp(-log(19)*(aLi-maturity[1])/(maturity[2]-maturity[1])))
    } else {
    Lmat = 1/(1+exp(-log(19)*(Li-maturity[1])/(maturity[2]-maturity[1])))
    }
  if(SELECT==TRUE){
  Lsel_a = 1/(1+exp(-log(19)*(Li-selpars[1])/(selpars[2]-selpars[1])))
  Lsel_b = dnorm(Li,selpars[3],selpars[3]*selpars[4])/max(dnorm(Li,selpars[3],selpars[3]*selpars[4]))
  Lsel_c = 1+(selpars[5]-1)*(Lsel_b-1)/(-1)
  Lsel = as.numeric(ifelse(rep(selpars[3],length(Li))>=Li,Lsel_a,Lsel_c))
  } else {
  Lsel = Lmat  
  }
  Ldat = data.frame(Li,LW,Lmat,Lsel) 
  Adat = data.frame(age,L_f,W_f,mat,L_m,W_m,sel_m)
  
  # compute unfished Spawning biomass per recruit (SBR0)
  n0_f = 0.5*exp(-(M*(0:(nages-1))))
  if(PlusGroup==TRUE){
    n0_f[nages] <- n0_f[nages]/(1-exp(-M))
  }  
  # Only female spawning if 2 sexes
  SBR0 = sum(n0_f*W_f*mat)*(2/nsexes)
  
  # compute Spawning Biomass per rectuit as a function of F (SBR)
  
  for (t in 1:nages)
  {
    
    if(t==1) ntilda_f[,t] <- 0.5
    if(t>1) ntilda_f[,t] <- ntilda_f[,t-1]*exp(-(M+sel_f[t-1]*F_i)) 
    if(PlusGroup==FALSE){
      if(t==nages) ntilda_f[,t] <- ntilda_f[,t]  
    } else {
      if(t==nages) ntilda_f[,t] <- ntilda_f[,t]/(1-exp(-(M+sel_f[t]*F_i))) 
    }
  }
  
  for (t in 1:nages)
  {
    
    if(t==1) ntilda_m[,t] <- 0.5
    if(t>1) ntilda_m[,t] <- ntilda_m[,t-1]*exp(-(M+sel_m[t-1]*F_i)) 
    if(PlusGroup==FALSE){
      if(t==nages) ntilda_m[,t] <- ntilda_m[,t]  
    } else {
      if(t==nages) ntilda_m[,t] <- ntilda_m[,t]/(1-exp(-(M+sel_m[t]*F_i))) 
    }
  }
  
  
  ntilda_f[ntilda_f < 0] = 0 
  ntilda_m[ntilda_m < 0] = 0 
  
  
  #get spawner biomass per recruit
  SBR= apply(ntilda_f %*% diag(W_f) %*% diag(mat),1,sum)*(2/nsexes)
  #get Yield per recruit
  
  # get Z matrix              
  Z_f = M+t(sel_f%*%t(F_i))
  Z_m = M+t(sel_m%*%t(F_i))
  
  # Yield-per-recruit matrix
  YPR = apply(ntilda_f %*% diag(W_f) %*% diag(sel_f)*F_i/Z_f*(1-exp(-Z_f)),1,sum)
  YPR = YPR +  apply(ntilda_m %*% diag(W_m) %*% diag(sel_m)*F_i/Z_m*(1-exp(-Z_m)),1,sum)
  
  #get Exploitable Biomass per recruit
  EBR = apply(ntilda_f %*% diag(W_f) %*% diag(sel_f),1,sum)
  EBR = EBR + apply(ntilda_m %*% diag(W_m) %*% diag(sel_m),1,sum)
  #get Exploitable Numbers per recruit
  ENR = apply(ntilda_f %*% diag(sel_f),1,sum)
  ENR = ENR + apply(ntilda_m %*% diag(sel_m),1,sum)
  
  
  # equilibrium recruitment
  RF = R0*(4*h*SBR-(1-h)*SBR0)/(SBR*(5*h-1))
  recruits = ifelse(RF<0,0,RF)
  
  # get spawner biomass depletion
  SBtoSB0 = SBR*recruits/SBR0
  
  
  # return Spawning Biomass as fraction of unfished levels
  return(list(SBtoSB0=SBR*recruits/SBR0,Yield=YPR*recruits, EB=EBR*recruits,SB=SBR*recruits,EN=ENR*recruits,Ldat=Ldat,Adat=Adat))
  
}

#-----------------------------------------------------------  
# Section for JABBA-SELECT 
#-----------------------------------------------------------
# Plot stock functions
if(SELECT==TRUE){selpars = rbind(SL50,SL95,SL.desc,CV.desc,min.desc)}else{
  selpars = rbind(Linf/3,Linf/2,SL.desc,CV.desc,min.desc)} # Unused place holder if SELECT == FALSE

stock =  get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,1],h=h,M=M,0,R0=1,PlusGroup,nsexes,SELECT)
Par = list(mfrow=c(2,2),mai=c(0.45,0.45,0,.15),omi = c(0.1,0.1,0.1,0) + 0.1,mgp=c(2,1,0), tck = -0.02,cex=0.7)
png(file = paste0(input.dir,"/StockFunctions_",assessment,".png"), width = 7, height = 6, 
    res = 200, units = "in")
cols=jabba.colors
par(Par)
plot(stock$Adat[,1],stock$Adat[,2],type="l",lwd=2,col=1,ylim=c(0,max(Linf)),xlab="Age",ylab="Length")
if(nsexes==2){lines(stock$Adat[,1],stock$Adat[,2],col=cols[1],lwd=2)
  lines(stock$Adat[,1],stock$Adat[,5],col=cols[2],lwd=2)}
plot(stock$Ldat[,1],stock$Ldat[,2],type="l",lwd=2,col=1,ylim=c(0,max(stock$Ldat[,2])),xlab="Length",ylab="Weigth")
if(nsexes==2) {lines(stock$Ldat[,1],stock$Ldat[,2],col=cols[1],lwd=2)
  lines(stock$Ldat[,1],stock$Ldat[,3],col=cols[2],lwd=2)}

plot(stock$Adat[,1],stock$Adat[,3],type="l",lwd=2,col=1,ylim=c(0,max(stock$Ldat[,2])),xlab="Length",ylab="Weigth")
if(nsexes==2){ lines(stock$Adat[,1],stock$Adat[,3],col=cols[1],lwd=2)
               lines(stock$Adat[,1],stock$Adat[,6],col=cols[2],lwd=2)}

plot(stock$Ldat[,1],stock$Ldat[,4],type="l",lwd=2,col=cols[1],ylim=c(0,max(1)),xlab="Length",ylab="Proportion")
for(j in 1:nSel){
stock =  get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,j],h=h,M=M,0,R0=1,PlusGroup,nsexes,SELECT)
lines(stock$Ldat[,1],stock$Ldat[,5],col=cols[j+1],lwd=2)
#points(stock$Adat[,2],stock$Adat[,4],pch=15,cex=0.8,col=cols[j+4])
}
lines(stock$Ldat[,1],stock$Ldat[,4],col=cols[1],lwd=2)
legend(ifelse(min(selpars[1,])>max(Linf)-max(selpars[1,]),"topleft","right"),c("Mature",paste0("Sel",1:nSel)),col=cols[c(1,2:(1+nSel))],bty="n",lwd=2,cex=0.9)
dev.off()

#-----------------------------------------------------------
# Generate Priors 
#-----------------------------------------------------------


Par = list(mfrow=c(2,2),mai=c(0.6,0.1,0,.1),omi = c(0.1,0.1,0.1,0) + 0.1,mgp=c(2.2,1,0), tck = -0.02,cex=0.75)
png(file = paste0(input.dir,"/Prior_hM_",assessment,".png"), width = 7.2, height = 7., 
    res = 200, units = "in")
par(Par)
SB0.pr = plot_lnorm(mu=mu.SB0,CV=CV.SB0,Prior=paste0("Prior B(",years[1],")/B0"))


if(psi.prior=="beta"){
  psi.pr = get_beta(mu=mu.psi,CV=CV.psi,Min=0,Prior=paste0("Prior B(",years[1],")/B0"))} else {
  psi.pr = plot_lnorm(mu=mu.psi,CV=CV.psi,Prior=paste0("Prior B(",years[1],")/B0"))  
  
}

h.pr = get_beta(mu=h,CV=CV.h,Min=0.2,Prior="Steepness h")
M.pr = get_gamma(M,CV.M,Prior="M")
mtext(paste("Density"), side=2, outer=T, at=0.5,line=0.4,cex=0.8)
dev.off()

# number of simulations for Monte-Carlo 
nsim =1000

# Admit uncertainty for M
sim.M = rgamma(nsim,M.pr[1],M.pr[2])
# Random variates of z
sim.h = rbeta(nsim, h.pr[1], h.pr[2])
sim.h=ifelse(sim.h<0.21,0.21,sim.h) # truncate at 0.2
if(runASEM==TRUE){

# matrices for storage
EBmsy.sim = matrix(rep(0,nsim*nSel),nsim,nSel)
Fmsy = matrix(rep(0,nsim*nSel),nsim,nSel)
Hmsy.sim = matrix(rep(0,nsim*nSel),nsim,nSel)
MSY.sim = matrix(rep(0,nsim*nSel),nsim,nSel)
SBmsy.sim = matrix(rep(0,nsim*nSel),nsim,nSel)
shape.sim = matrix(rep(0,nsim*nSel),nsim,nSel)
SPmax = matrix(rep(0,nsim*nSel),nsim,nSel)
cat(paste0("\n","- Generate informative priors for Hmsy and m","\n"))
ci = rep(c(1,rep(0,(round(nsim/50,1)-1))),nsim)

cat(paste0("\n","- Running Monte-Carlo with ",nsim," runs\n"))

for(i in 1:nsim)
{
  
  
  if(i==1) cat(paste("\n","|"))
  if(ci[i]==1) cat("*")
  if(i==nsim) cat(paste("|","\n"))

  for(j in 1:nSel)
  {
    
    
    F_i = (exp(seq(log(M*0.1),log(M*5),0.01)))
    # Find Hmsy
    depletion =  get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,j],h=sim.h[i],M=sim.M[i],F_i,R0=1,PlusGroup,nsexes,SELECT)
    getYield = depletion$Yield
    SPmax[i,j] = ifelse(depletion$Yield[length(depletion$Yield)]<max(depletion$Yield),1,0)
    
    # Hmsy at maximum yield
    Fmsy[i,j] = ifelse(SPmax[i,j]==0,1,F_i[getYield==max(getYield)])   
    
    pop_Fmsy = get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,j],h=sim.h[i],M=sim.M[i],F_i=Fmsy[i,j],R0=1,PlusGroup,nsexes,SELECT)
    EBmsy.sim[i,j]  = pop_Fmsy$EB
    MSY.sim[i,j] =   pop_Fmsy$Yield
    SBmsy.sim[i,j] =   pop_Fmsy$SB
    Hmsy.sim[i,j] = MSY.sim[i,j]/SBmsy.sim[i,j] 
    SBmsySB0 = ifelse(SPmax[i,j]==0,0.37,pop_Fmsy$SBtoSB0)
    
    rshape = 0
    
    # Find shape for  SBmsytoK 
    rshape = seq(0.1,5,0.001)
    
    check.shape =((rshape)^(-1/(rshape-1))-SBmsySB0)^2
      
    shape.sim[i,j] = rshape[check.shape==min(check.shape)]
    
   
    
  } # end of 2st loop
} # end of sim

}

# Determine non-linear relationship between SB and EB as a function of P = SB/K


# Determine correlation prior for Hmsy and m
Hmsy_m_sim =matrix(rep(0,nsim*nSel),nsim,2)
Hmsy_m_sim[,1] = (MSY.sim[,1]/apply(SBmsy.sim,1,mean))
Hmsy_m_sim[,2] = as.numeric(apply(shape.sim,1,mean))
# Take viable pair that attained max SP
if(length(which(SPmax[,1]==0))>0){
Hmsy_m_sim[which(SPmax[,1]==0),] = exp(rmvnorm(length(which(SPmax[,1]==0)),log(apply(Hmsy_m_sim[which(SPmax[,1]==1),],2,median)),cov(log(Hmsy_m_sim[which(SPmax[,1]==1),]))))
}

# if more than ne fishery (selectivity) express Hmsy as ratio
if(nSel>1){
ratio.Hmsy = matrix(rep(0,nsim*(nSel-1)),nsim,nSel-1)
dHmsy.pr=NULL
for(j in 2:nSel) {
ratio.Hmsy [,j-1] = (MSY.sim[,j]/MSY.sim[,1])   
dHmsy.pr = cbind(dHmsy.pr,as.numeric(fitdist(ratio.Hmsy[,j-1], "gamma",method="mme")$estimate))# conditioned on 1 period for linefish selectivity
}}

# get JAGS dmnrom prior
mu_prod_prior = apply(log(Hmsy_m_sim),2,mean) 
cov_prod_prior = cov(log(Hmsy_m_sim)) 
ln_Hmsy_m_devs = rmvnorm(nsim*10 ,mean = mu_prod_prior,sigma = cov_prod_prior)
Hmsy_m_devs = exp(ln_Hmsy_m_devs)




Par = list(mfrow=c(round((nSel+1)/2+.33,0)+1,2),mai=c(0.6,0.3,0,.15),omi = c(0.2,0.2,0.2,0) + 0.1,mgp=c(2,1,0), tck = -0.02,cex=1)
png(file = paste0(input.dir,"/PriorsInput",assessment,".png"), width = 8, height = 2.5*round((nSel+1)/2+2.5,0), 
    res = 200, units = "in")
par(Par)
M.pr = get_gamma(M,CV.M,Prior="M'")
legend("topright",c("gamma"),col=c("grey",rgb(0,0,1,0.5)),pch=15,bty="n",cex=1,pt.cex=1.5)

h.pr = get_beta(mu=h,CV=CV.h,Min=0.2,Prior="h'")
legend("topright",c("beta"),col=c("grey",rgb(0,0,1,0.5)),pch=15,bty="n",cex=1,pt.cex=1.5)

d=stats::density(Hmsy_m_devs[,1],adjust=2)
plot(d,type="l",xlim=range(d$x),xlab=bquote(H[MSY*","*S1]),yaxt="n",ylim=c(0,max(d$y)*1.25),main="")
polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$y))),col="grey")
d=stats::density(Hmsy_m_sim[,1],adjust=2)
polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$y))),col=rgb(0,0,1,0.3))
legend("topright",c("Monte-Carlo","MVN"),col=c("grey",rgb(0,0,1,0.5)),pch=15,bty="n",cex=1,pt.cex=1.5)
d=stats::density(Hmsy_m_devs[,2],adjust=2)
plot(d,type="l",xlim=range(d$x),xlab="shape m",yaxt="n",main="")
polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$y))),col="grey")    
d=stats::density(Hmsy_m_sim[,2],adjust=2)
polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$y))),col=rgb(0,0,1,0.3))
legend("topright",c("Monte-Carlo","MVN"),col=c("grey",rgb(0,0,1,0.5)),pch=15,bty="n",cex=1,pt.cex=1.5)
mtext(paste("Density"), side=2, outer=T, at=0.5,line=0.5,cex=1.1)

if(nSel>1){
# Express next priors as ratio to Hmsy1
  for(j in 2:nSel){ 
  d=stats::density(rgamma(5000,dHmsy.pr[1,j-1],dHmsy.pr[2,j-1]),adjust=2)
  plot(d,type="l",xlim=range(c(0.9,1.2,d$x)),xlab=bquote(H[MSY*",S"*.(j)] ~ "/" ~H[MSY*","*S1]),yaxt="n",main="")
  polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$y))),col="grey")    
  d=stats::density(ratio.Hmsy[,j-1],adjust=2)
  polygon(c(d$x,rev(d$x)),c(d$y,rep(0,length(d$y))),col=rgb(0,0,1,0.3))
  abline(v=1,lty=2)
  legend("topright",c("Monte-Carlo","gamma"),col=c("grey",rgb(0,0,1,0.5)),pch=15,bty="n",cex=1,pt.cex=1.5)
  }  
  }
  dev.off()  
    
   

# Hmsy_vs_m
# Prior correlation
Par = list(mfrow=c(1,2),mai=c(0.45,0.45,0,.15),omi = c(0.1,0.1,0.1,0) + 0.1,mgp=c(2.2,1,0), tck = -0.02,cex=0.7)
png(file = paste0(input.dir,"/Cor_m_Hmsy_",assessment,".png"), width = 7, height = 3, 
    res = 200, units = "in")
par(Par)
plot(log(Hmsy_m_devs[,1]),log(Hmsy_m_devs[,2]),type="n",ylab="log_m'",xlab=bquote(log_H[MSY*","*S1]),pch=19,col=grey(0.5,0.4))
points(log(Hmsy_m_devs[,1]),log(Hmsy_m_devs[,2]),col=rgb(0,0,1,0.2),pch=19)
points(log(Hmsy_m_sim[,1]),log(Hmsy_m_sim[,2]),pch=21,bg=grey(0.6,1),cex=1.2)
legend("topright",c("Monte-Carlo","MVN approximation"),col=c(1,4),pch=19,cex=1.1,bty="n")

plot(Hmsy_m_sim[,1],Hmsy_m_sim[,2],ylab="m",xlab=bquote(H[MSY*","*S1]),pch=19,col=rgb(0.,0.,1,0.5),type="n")
points(Hmsy_m_devs[,1],Hmsy_m_devs[,2],col=rgb(0,0,1,0.2),pch=19)
points(Hmsy_m_sim[,1],Hmsy_m_sim[,2],pch=21,bg=grey(0.6,1),cex=1.2)
legend("topright",c("Monte-Carlo","MVN approximation"),col=c(1,4),pch=19,cex=1.1,bty="n")
dev.off()

Par = list(mfrow=c(1,1),mai=c(0.55,0.3,0,.15),omi = c(0.0,0.2,0.2,0) + 0.1,mgp=c(2.5,1,0), tck = -0.02,cex=0.8)
png(file = paste0(input.dir,"/EBtoSB",assessment,".png"), width = 5.5, height = 4, 
    res = 200, units = "in")
par(Par)

ylim = c(0,5)
xlim=c(0,1)
cols = jabba.colors# rainbow(nSel+1)
# Find sets.EB
nI = length(sets.I)
sets.EB = selEB = unitEB = qEB= check.q =NULL 
k=l=1
for(j in 1:nI){
sets.EB[j] = which((paste0(sets.I,I.unit))[j]==unique(paste0(sets.I,I.unit)) )
if(length(which(sets.EB%in%sets.EB[j]))==1){
selEB[k] = sets.I[j] 
unitEB[k] = I.unit[j]
k = k+1
}
}

nEB = length(unique(sets.EB))
coefT = matrix(0,5,nEB)
ymax=ymin=NULL
for(j in 1:nEB)
{
  depletion = get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,selEB[j]],h=h,M=M,F_i,R0=1,PlusGroup,nsexes,SELECT)
  y = (depletion$EB/(depletion$SB+0.000001))
  y = y[depletion$SBtoSB0>0.1]
  ymax = max(y,ymax)
  ymin = min(y[y>0.01],ymin)
  } 

plot(1,1,type="n",xlab=expression(SB/SB[0]),ylab="",ylim=c(min(ymin,0.7),max(ymax,1.3)),xlim=c(0,0.9))
EBSB = NULL
for(j in 1:nEB){
  depletion = get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,selEB[j]],h=h,M=M,F_i,R0=1,PlusGroup,nsexes,SELECT)
  y = (depletion$EB/depletion$SB)
  if(unitEB[j]!=1){ y = depletion$EN/depletion$SB*depletion$EB[1]/depletion$EN[1]
  #y = y*depletion$EB[1]/depletion$EN[1]
  }
  x = depletion$SBtoSB0
  
  dat = data.frame(x,y)
  dat=subset(dat,x>0.1 & x<0.9)  
  EBSB[j] = mean(dat$y) 
  
  #dat$y=dat$y/max(dat$y)
  P1 = min(dat$x)
  P2 = max(dat$x)
  
  
  if(coef(lm(y~x,data=dat))[1] < 0.97 | coef(lm(y~x,data=dat))[1] > 1.03 & abs(coef(lm(y~x,data=dat))[2]) >0.001){ 
  fit = nls((y)~v1+(v2-v1)*((1-exp(-v3*((x)-P1)))/(1-exp(-v3*(P2-P1)))),data=dat,start=list(v1=dat[nrow(dat),2],v2=dat[1,2],v3=1))
  coefT[,j] = c(as.numeric(coef(fit)),P1,P2)
  } else {
  fit =  lm(y~x,data=dat)  
  coefT[,j] = c(0.999999,0.99999,0.0000001,0.1,0.9)
  }
  if(unitEB[j]==1){
  points(dat$x,dat$y,cex=0.5)
  lines(dat$x,(predict(fit)),col=cols[j+1],lwd=2)} else {
    points(dat$x,dat$y,cex=0.5)
    lines(dat$x,(predict(fit)),col=cols[j+1],lwd=2)  
    }
  
  abline(h=1,lty=2,col=cols[1])
}
legend("topright",c("Maturity",paste0("Sel",selEB[1:nEB],".",ifelse(unitEB[1:nEB]==rep(1,nEB),"Bio","N"))),lwd=c(1,rep(2,100)),lty=c(2,rep(1,100)),col=cols[1:(j+1)],bty="n")
mtext(paste("EB/SB"), side=2, outer=T, at=0.6,line=0.6,cex=0.8)
#  mtext(paste("B/K"), side=1, outer=T, at=0.5,line=1,cex=0.9)

dev.off()

# Plot MSY
Par = list(mfrow=c(1,1),mai=c(0.6,0.3,0,.15),omi = c(0.1,0.2,0.2,0) + 0.1,mgp=c(2.5,1,0), tck = -0.02,cex=0.8)
png(file = paste0(input.dir,"/Prodution",assessment,".png"), width = 6, height = 5, 
    res = 200, units = "in")
par(Par)
colsp = cols[-1] # jabba.colors[c(4,5,6,7,2,11,12,3,8)]  

j = 1

ymax = NULL
  for(j in 1:nSel){
  Y.it = get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,j],h=h,M=M,F_i,R0=1,PlusGroup,nsexes,SELECT)$Yield
  ymax = c(ymax,max(Y.it))  
  }  
  ymax = max(ymax)  
  plot(Y.it,Y.it,type="n",ylab="Yield",xlab="SB/SB0",lwd=2,yaxt="n",xlim=c(0,1),ylim=c(0,ymax))
  SBR0 = get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,j],h=h,M=M,0,R0=1,PlusGroup,nsexes,SELECT)$SB
  P = seq(0.0001,1,0.001) 
  B = P*SBR0
  
  for(j in 1:nSel){
  SB_K = get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,j],h=h,M=M,F_i,R0=1,PlusGroup,nsexes,SELECT)$SBtoSB0
  Y.it = get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,j],h=h,M=M,F_i,R0=1,PlusGroup,nsexes,SELECT)$Yield
  dens = stats::density(Y.it)
  lines(SB_K,Y.it,col=colsp[j],lwd=2)
  Fmsy = F_i[Y.it==max(Y.it)]   
  pop_SP = get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,j],h=h,M=M,F_i=Fmsy,R0=1,PlusGroup,nsexes,SELECT)
  MSY = pop_SP$Yield
  SBmsy = pop_SP$SB
  Hmsy = MSY/SBmsy 
  SBmsySB0 = SBmsy/get_ASEM(Linf,kappa,t0,aW,bW,minage,maxage,maturity,selpars[,j],h=h,M=M,F_i=0,R0=1,PlusGroup,nsexes,SELECT)$SB
  rshape = 0
  # Find shape for  SBmsytoK 
  rshape = seq(0.1,5,0.001)
  check.shape =((rshape)^(-1/(rshape-1))-SBmsySB0)^2
  shape = rshape[check.shape==min(check.shape)]
  SP = Hmsy/(1-1/shape)*B*(1-(B/SBR0)^(shape-1))
  lines(P,SP,lty=3,col=colsp[j],lwd=2)
  #PMSP = P[SP==max(SP)]
  PMSY = SB_K[Y.it==max(Y.it)]
  lines(rep(PMSY,2),c(0,max(SP)),lty=1,col=colsp[j])
  }
  mtext(paste("Yield"), side=2, outer=T, at=0.6,line=0.5,cex=0.9)
  legend("topright",c("ASPM","SPM",paste("Sel",1:nSel)),col=c(1,1,colsp),lwd=2,lty=c(1,3,rep(1,nSel)),bty="n")  
  dev.off()

 
  
  # PRIORS
  SB0.prior = c(round(SB0.pr[1],0),sqrt(log(CV.SB0^2+1)))
  Psi.prior = c(mu.psi,CV.psi)  
  Hmsy_m.prior = NULL
  for(i in 1:2) Hmsy_m.prior = cbind(Hmsy_m.prior ,round(c(mean(Hmsy_m_devs[,i]),sd(Hmsy_m_devs[,i])/mean(Hmsy_m_devs[,i])),3)) 
  M.prior = c(M,CV.M)
  h.prior = c(h, CV.h)
  Priors =rbind(SB0.prior,Psi.prior,t(Hmsy_m.prior),M.prior,h.prior)
  row.names(Priors) = c("SB0","Psi","Hmsy1","shape m","M","h")
  colnames(Priors) = c("Mean","CV")                          
  write.csv(Priors,paste0(input.dir,"/Priors",assessment,".S",s,".csv"))
 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
#---------------------------------------------------------------
# Setup JABBA-SELECT
#---------------------------------------------------------------
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
# remove scientific numbers
#options(scipen=999)
##########################################################################


# starting values
nq = length(unique(sets.q))
nvar = length(unique(sets.var))
nI = ncol(CPUE)


#----------------------------------------------------------
# Setup TAC projection
#---------------------------------------------------------
if(Projection==TRUE){
  nTAC = length(TACs)
  TAC = mat.or.vec(pyrs,nTAC)
  yr.last = max(years) # assessment year  
  
  for(i in 1:nTAC){
    TAC[,i] = c(rep(TACint,imp.yr-yr.last-1),rep(TACs[i],pyrs-(imp.yr-yr.last-1)))  
  }
  
}else{
  nTAC = 1  
  TAC = TC[n.years] #  
  pyrs = 1
}


#stI = c(1:n.years)[is.na(apply(CPUE,1,mean,na.rm=TRUE))==FALSE][1] #first year with CPUE
stI = ifelse(proc.dev.all==TRUE,1, which(proc.dev.all==years)) #first year with CPUE

inits <- function(){list(SB0= rlnorm(1,log(mu.SB0),0.3),q = (aggregate(EBSB[sets.EB]~sets.q,FUN=mean)[,2]*runif(nq,min(CPUE,na.rm=T)/max(apply(Catch,1,sum),na.rm=T),mean(CPUE,na.rm=T)/max(apply(Catch,1,sum)))),  proc.est=ifelse(proc.type=="igamma",runif(1,30,100),1/runif(1,20,100)), itau2=runif(nvar,80,200), psi=mu.psi,log.prod.pars=apply(rmvnorm(50 ,mean = mu_prod_prior,sigma = cov_prod_prior),2,mean))}

if(init.values==TRUE){
  inits <- function(){list(SB0= init.SB0,q = init.q, proc.est=ifelse(proc.type=="igamma",runif(1,30,100),1/runif(1,20,100)), itau2=runif(nvar,80,200), psi=mu.psi,log.prod.pars=apply(rmvnorm(500 ,mean = mu_prod_prior,sigma = cov_prod_prior),2,mean))}
}

if(sigma.proc!=TRUE) pr.proc = c(log(0.07),0.1)

surplus.dat = list(N=n.years,C=Catch, TC = TC,I=CPUE,SE2=se2,mu.prod.prior = apply(log(Hmsy_m_sim),2,mean),cov.prod.prior = cov(log(Hmsy_m_sim)) ,psi.pr=psi.pr,SB0.pr = SB0.pr,coefT=coefT,nSel=nSel,nq=nq,nC=length(sets.C),nI = n.indices,nvar=nvar,
                   sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc),sets.var=sets.var, sets.q=sets.q,sets.EB=sets.EB,sets.C=sets.C,
                   SBmsytoSB0 = ifelse(is.null(SBmsy_SB0)==TRUE,0,SBmsy_SB0),nTAC=nTAC,pyrs=pyrs,TAC=TAC,pr.proc = pr.proc,stI=stI,P_bound=P_bound,pen.bk = rep(0,n.years),q_bounds=q_bounds,sigmaobs_bound=sigmaobs_bound,sigmaproc_bound=sigmaproc_bound,SB0_bounds=SB0_bounds)


# add additional Hmsy ratio priors if nSel > 1 
if(nSel>1) surplus.dat$dHmsy.pr =as.matrix(dHmsy.pr) 


JABBA_SELECT = "JABBA_SELECT.jags"

# PARAMETERS TO MONITOR
params <- c("SB0","r", "q", "psi","sigma2", "tau2","m","Hmsy","SBmsy", "MSY","Href","SBref","Yref", "BtoBmsy","HtoHmsy","CPUE","Ihat","Proc.Dev","P","SB","Hmsy.y","prP","prBtoBmsy","prHtoHmsy","TOE","H","EB")

cat(paste0("\n","><> RUN ",Mod.names," model for ",assessment," ",Scenario," in JAGS <><","\n","\n"))




#--------------------------
# Capture Settings
#--------------------------

Settings = surplus.dat
Settings$Stock =  stockpars 
Settings$Model.type = Mod.names
Settings$proc.dev.all = proc.dev.all
Settings$Bmsy_target = ifelse(is.null(SBmsy_SB0),"Bmsy",SBmsy_SB0)
Settings$Do.Projection = Projection
Settings$TAC.implementation = imp.yr
Settings$catch.metric = catch.metric
Settings$harvest.label = harvest.label
Settings$Run.CPUE.avg.tool = CPUE.plot  
Settings$Use.avg.CPUE = meanCPUE
Settings$Specify.init.values = init.values
Settings$save.trajectories = save.trajectories
Settings$save.large.posterior.object = save.all 
Settings$Seed = get_seed
Settings$MCMC.ni = ni
Settings$MCMC.saved.steps = nt
Settings$MCMC.burnin = nb
Settings$MCMC.Chains = nc
Settings$MCMC.ni = ni
Settings$MCMC.saved = nsaved 
capture.output( Settings, file=paste0(input.dir,"/Settings.txt"))
#-------------------------------------------------------------------




# For parallel processing
# JAGS MODEL
sink("JABBA_SELECT.jags")
cat("
    
    model {
    
    # Prior specifications  

      
    eps <- 0.000001    
    
    #Catchability
    for(i in 1:nq){   
    #iq[i] ~ dgamma(0.001,0.001)
    #q[i] <- pow(iq[i],-1)
    q[i] ~ dunif(q_bounds[1],q_bounds[2])
    
    }  
    
    
    # Biomass depletion at the start
    ")

    if(psi.prior =="beta"){
    cat("
      psi ~ dbeta(psi.pr[1],psi.pr[2])
      ",append=TRUE)
    } else {
    cat("
      psi ~ dlnorm(log(psi.pr[1]),pow(psi.pr[2],-2)) #I(0.1,1.1)    
      ",append=TRUE)   # bias correct or not?
      }

    if(sigma.proc==TRUE){
    if(proc.type=="igamma"){
     cat("
      proc.est ~ dgamma(pr.proc[1],pr.proc[2])
      isigma2 <- proc.est 
      sigma2 <- pow(isigma2,-1)
      sigma <- sqrt(sigma2)
      fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
         ",append=TRUE)}else{
      cat("
      proc.est ~ dlnorm(pr.proc[1],pr.proc[2])
      sigma <- proc.est 
      sigma2 <- pow(sigma,2)
      isigma2 <- pow(sigma,-2)
      fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
          
      ",append=TRUE)}}else{ 
      cat(" 
      proc.est ~ dlnorm(pr.proc[1],pr.proc[2])
      isigma2 <- pow(sigma.fixed+eps,-2) 
      sigma2 <- pow(isigma2,-1)
      sigma <- sqrt(sigma2)
           
      ",append=TRUE)}

if(sigma.est==TRUE){
  cat("
      # Obsevation variance
      for(i in 1:nvar)
      {
      # Observation error
      itau2[i]~ dgamma(0.001,0.001)
      tau2[i] <- 1/itau2[i]
      }
      
      for(i in 1:nI)
      {
      for(t in 1:N)
      {
      var.obs[t,i] <- SE2[t,i]+tau2[sets.var[i]] 
      ivar.obs[t,i] <- 1/var.obs[t,i]
      # note total observation error (TOE)     
      TOE[t,i] <- sqrt(var.obs[t,i]) # Total observation variance
      
      }}
      ",append=TRUE)  
}else{ cat(" 
      # Obsevation variance
           for(i in 1:nvar)
           {
           # Observation error
           itau2[i]~ dgamma(4,0.01)
           tau2[i] <- 1/itau2[i]
           }
           
           for(i in 1:nI)
           {
           for(t in 1:N)
           {
           var.obs[t,i] <- SE2[t,i] # drop tau2
           fake.tau[t,i] <- tau2[sets.var[i]]
           
           ivar.obs[t,i] <- 1/var.obs[t,i]
           # note total observation error (TOE)     
           TOE[t,i] <- sqrt(var.obs[t,i])
           
           }}
           
           ",append=TRUE)}

    # Run rest of code  
      cat("  

    # Carrying Capacity SB0
    SB0 ~ dlnorm(log(SB0.pr[1]),pow(SB0.pr[2], -2))
    
    
    # get correlated production priors Hmsy1 and m
    log.prod.pars[1:2] ~ dmnorm(mu.prod.prior[1:2],inverse(cov.prod.prior[1:2,1:2]))
    
    
    Hmsy[1] <- exp(log.prod.pars[1])
    # informative priors for m
    m <- exp(log.prod.pars[2])
   ",append=TRUE)
    
    if(nSel>1){
    cat("
          
      for(i in 1:(nSel-1)){
        dHmsy[i] ~ dgamma(dHmsy.pr[1,i],dHmsy.pr[2,i])
        Hmsy[i+1] <- Hmsy[1]*dHmsy[i] 
         }
        
          ",append=TRUE)}
      
      # Run rest of code  
      cat("  
          
    # where r =Hmsy*(m-1)/(1-1/m), and m is the shape
    for(i in 1:nSel){
    r[i] <-  Hmsy[i]*(m-1)/(1-1/m)
    }
    
    #  Catch multipliers for Hmsy
    for(t in 1:N){
    for(i in 1:nC){
    rf[t,i] <- C[t,i]/TC[t]*Hmsy[sets.C[i]] 
    }}
    
    for(t in 1:N){
    Hmsy.y[t] <- sum(rf[t,])
    }
    
    #Process equation in SB
    Pmean[1] <- log(psi)
    #iPV[1] <- ifelse(1<(stI),10000,isigma2) # inverse process variance (was 10000)
    #P[1] ~ dlnorm(Pmean[1],iPV[1]) # set to small noise instead of isigma2
    iPV[1] <- ifelse(1<(stI),10000,isigma2) # inverse process variance
    lnbias[1] <- ifelse(1<(stI),0,-0.5*sigma2)
    
    P[1] ~ dlnorm(Pmean[1]+lnbias[1],iPV[1]) # Added bias correction
    

    penB[1]  <- ifelse(P[1]<(P_bound[1]),log(SB0*P[1])-log(SB0*P_bound[1]),ifelse(P[1]>P_bound[2],log(SB0*P[1])-log(SB0*(P_bound[2])),0)) # penalty if Pmean is outside viable biomass
    

    # Process equation
    for (t in 2:N) 
    {
    Pmean[t] <- log(max(P[t-1] +  Hmsy.y/(1-1/m)*P[t-1]*(1-pow(P[t-1],m-1)) - TC[t-1]/SB0,0.01))
    #iPV[t] <- ifelse(t<(stI),10000,isigma2) # inverse process variance
    #P[t] ~ dlnorm(Pmean[t],iPV[t])
    iPV[t] <- ifelse(t<(stI),10000,isigma2) # inverse process variance
    lnbias[t] <- ifelse(t<(stI),0,-0.5*sigma2)
    P[t] ~ dlnorm(Pmean[t]+lnbias[t],iPV[t]) # Added bias correction
    penB[t]  <- ifelse(P[t]<P_bound[1],log(SB0*P[t])-log(SB0*P_bound[1]),ifelse(P[t]>P_bound[2],log(SB0*P[t])-log(SB0*P_bound[2]),0)) # penalty if Pmean is outside viable biomass
    }
    
    # Process error deviation 
    for(t in 1:N){
    Proc.Dev[t] <- log(P[t]*SB0)-log(exp(Pmean[t])*SB0)} 
    
    # Enforce penaly
    for(t in 1:N){
    pen.bk[t] ~ dnorm(penB[t],1000) # enforce penalty with CV = 0.1
    }
    
  

    for (t in 1:N) 
    { 
    SB[t] <- SB0*P[t]    
    H[t] <- TC[t]/SB[t] 
    }
    
    # Observation equation in related to EB
    
    for(i in 1:nI)
    {
    for (t in 1:N) 
    { 
    EBtoSBmean[t,i] <- log(max(coefT[1,sets.EB[i]]+(coefT[2,sets.EB[i]]-coefT[1,sets.EB[i]])*((1-exp(-coefT[3,sets.EB[i]]*(P[t]-coefT[4,sets.EB[i]])))/(1-exp(-coefT[3,sets.EB[i]]*(coefT[5,sets.EB[i]]-coefT[4,sets.EB[i]])))),0.01)) 
    EBtoSB[t,i] <- exp(EBtoSBmean[t,i])
    Imean[t,i] <- log(q[sets.q[i]]*P[t]*SB0*EBtoSB[t,i]);
    I[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]));
    EB[t,i] <- P[t]*SB0*EBtoSB[t,i]   # Exploitable biomass
    CPUE[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]))#q[sets.q[i]]*P[t]*SB0*EBtoSB[t,i]
    Ihat[t,i]  <- exp(Imean[t,i])
    }}
    
    
    #Management quantaties
    SBmsy_SB0 <- (m)^(-1/(m-1))
    SBref_SB0 <- ifelse(SBmsytoSB0>0,SBmsytoSB0,(m)^(-1/(m-1)))
    
    SBmsy <- SBmsy_SB0*SB0
    for(i in 1:nSel){
    MSY[i] <- SBmsy*Hmsy[i]; # Fix
    }
    
    SBref <- SBref_SB0*SB0
    
          


    for(i in 1:nSel){
    Yref[i] <- SB0*(Hmsy[i]/(1-1/m)*SBref_SB0*(1-SBref_SB0^(m-1)))
    Href[i] <- Yref[i]/SBref
    }
     
    
    #  Catch multipliers for Href.y
    for(t in 1:N){
    for(i in 1:nC){
    rf.ref[t,i] <- C[t,i]/TC[t]*Href[sets.C[i]] 
    }}
          
    for(t in 1:N){
    Href.y[t] <- sum(rf.ref[t,])
    }    
      
    
    for (t in 1:N){
    # use x y to put them towards the end of the alphabetically sorted  mcmc object
    BtoBmsy[t] <- SB[t]/SBref
    HtoHmsy[t] <- H[t]/(Href.y[t]) 
    }

    
    # Enforce soft penalty on K if < SB0_bounds >  
    SB0.pen ~ dnorm(penSB0,1000) # enforce penalty 
          penSB0  <- ifelse(SB0<(SB0_bounds[1]),log(SB0)-log(SB0_bounds[1]),ifelse(SB0>SB0_bounds[2],log(SB0)-log(SB0_bounds[2]),0)) # penalty if Pmean is outside viable biomass
          
          
    # Enforce soft penalty on process deviance if sigma.proc > 0.2 
    proc.pen ~ dnorm(penProc,1000) # enforce penalty 
    penProc  <- ifelse(sigma>sigmaproc_bound,log(sigma)-log(sigmaproc_bound),0) 
          
          
    # Enforce soft penalty on observation error if sigma.obs > sigma_bound 
    for(i in 1:nvar){
    obs.pen[i] ~ dnorm(penObs[i],1000) # enforce penalty 
    penObs[i]  <- ifelse(pow(tau2[i],0.5)>sigmaobs_bound,log(pow(tau2[i],0.5))-log(sigmaobs_bound),0) 
    }
          
    ", append=TRUE)
      
      # PROJECTION
      if(Projection==TRUE){
        cat("
            for(i in 1:nTAC){
            # Project first year into the future
            prPmean[1,i] <- log(max(P[N] +  Hmsy.y[N]/(1-1/m)*P[N]*(1-pow(P[N],m-1)) - TC[N]/SB0,0.005))
            prP[1,i] ~ dlnorm(prPmean[1,i],isigma2) 
            # Project all following years
            for(t in 2:pyrs){
            prPmean[t,i] <- log(max(prP[t-1,i] +  Hmsy.y[N]/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1)) - TAC[t-1,i]/SB0,0.001))
            # process error (as monte-carlo simular)
            prP[t,i] ~ dlnorm(prPmean[t,i],isigma2)}
            for(t in 1:pyrs){
            prB[t,i] <- prP[t,i]*SB0
            prH[t,i] <- TAC[t,i]/prB[t,i]
            prHtoHmsy[t,i] <- prH[t,i]/Hmsy.y[N]
            prBtoBmsy[t,i] <- prB[t,i]/SBmsy
            }}  
            ",append=TRUE)} else {
              cat("
                  #Prevent error for unused input    
                  fakeTAC <-  TAC
                  fakepyrs <- pyrs 
                  fakenTAC <- nTAC    
                  prHtoHmsy <- 1
                  prP <- 1 
                  prBtoBmsy <- 1    
                  ", append=TRUE)}  
      
      cat("
          
      } # END OF MODEL
          ",append=TRUE,fill = TRUE)
      sink()
      
      
      
ptm <- proc.time()
      
mod <- jags(surplus.dat, inits,params,paste(JABBA_SELECT), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)  # adapt is burn-in

proc.time() - ptm
save.time = proc.time() - ptm

cat(paste0("\n",paste0("><> Scenario ",Scenario," completed in ",as.integer(save.time[3]/60)," min and ",round((save.time[3]/60-as.integer(save.time[3]/60))*100)," sec <><","\n")))

cat(paste0("\n","><> Produce results output for ",assessment," ",Scenario," <><","\n"))

# if run with library(rjags)
posteriors = mod$BUGSoutput$sims.list


#-----------------------------------------------------------
# <><<><<><<><<><<><<>< Outputs ><>><>><>><>><>><>><>><>><>
#-----------------------------------------------------------
output.dir = paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Output")
dir.create(output.dir, showWarnings = FALSE)

# run some mcmc convergence tests
par.dat= data.frame(posteriors[params[c(1,8,3:7)]])
geweke = geweke.diag(data.frame(par.dat))
pvalues <- 2*pnorm(-abs(geweke$z))
pvalues

heidle = heidel.diag(data.frame(par.dat))
heidle
# postrior means + 95% BCIs
#Model  parameter
apply(par.dat,2,quantile,c(0.025,0.5,0.975))

man.dat = data.frame(posteriors[params[8:13]])
#Management quantaties
apply(man.dat,2,quantile,c(0.025,0.5,0.975))

# Depletion
Depletion = posteriors$P[,c(1,n.years)]
colnames(Depletion) = c(paste0("P",years[1]),paste0("P",years[n.years]))

H_Hmsy.cur = posteriors$HtoHmsy[,c(n.years)]
B_Bmsy.cur = posteriors$BtoBmsy[,c(n.years)]


man.dat = data.frame(man.dat,Depletion,B_Bmsy.cur,H_Hmsy.cur)

results = round(t(cbind(apply(par.dat,2,quantile,c(0.025,0.5,0.975)))),3)

results = data.frame(Median = results[,2],LCI=results[,1],UCI=results[,3],Geweke.p=round(pvalues,3),Heidelberger.p=round(heidle[,3],3))

ref.points = round(t(cbind(apply(man.dat,2,quantile,c(0.025,0.5,0.975)))),3)

ref.points = data.frame(Median = ref.points[,2],LCI=ref.points[,1],UCI=ref.points[,3])

# get number of parameters
npar = length(par.dat)

#-------------------------------------------------------------------------
# Save parameters, results table and current status posterior in csv files
#-------------------------------------------------------------------------

# Safe posteriors
if(save.all==TRUE) save(posteriors,file=paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_posteriors.Rdata"))


# Save model estimates and convergence p-values
write.csv(data.frame(results),paste0(output.dir,"/Estimates_",assessment,"_",Scenario,".csv"))

# Make standard results table with parameter estimates and reference points
Table = rbind(data.frame(results)[c("SB0","m"),1:3],data.frame(ref.points))  
write.csv(Table,paste0(output.dir,"/Results_",assessment,"_",Scenario,".csv"))

#Save posterior of recent assessment year (KOBE posterior)
write.csv(data.frame(BtoBmsy=B_Bmsy.cur,FtoFmsy=H_Hmsy.cur),paste0(output.dir,"/Status_posterior",assessment,".csv"))  

#-----------------------------------------------
# Stock trajectories
#-----------------------------------------------

Bt = posteriors$SB
Ht = posteriors$H
Bt_Bmsy = posteriors$BtoBmsy
Ht_Hmsy = posteriors$HtoHmsy
Bt_K = posteriors$P
Stock_trj = cbind(t(apply(Bt,2,quantile,c(0.5,0.025,0.975))),
                  t(apply(Ht,2,quantile,c(0.5,0.025,0.975))),
                  t(apply(Bt_Bmsy,2,quantile,c(0.5,0.025,0.975))),
                  t(apply(Ht_Hmsy,2,quantile,c(0.5,0.025,0.975))),t(apply(Bt_K,2,quantile,c(0.5,0.025,0.975))))


colnames(Stock_trj) = paste0(rep(c("SBt",ifelse(harvest.label=="Hmsy","Ht","Ft"),"SBt_SBmsy",ifelse(harvest.label=="Hmsy","Ht_Hmsy","Ft_Fmsy"),"SBt_SB0"),each=3),rep(c(".Median",".LCI95%",".UCI95%"),4))
rownames(Stock_trj) = years


# Save results
write.csv(Stock_trj,paste0(output.dir,"/Stock_trj.csv"))



if(save.trajectories==TRUE){
  cat(paste0("\n","><> Saving Posteriors of FRP trajectories <><","\n"))
  
  # FRP trajectories
  trajectories = array(NA,c(nsaved,n.years,3))
  trajectories[,,1] = posteriors$P 
  trajectories[,,2] = posteriors$BtoBmsy 
  trajectories[,,3] = posteriors$HtoHmsy
  
  kb=kobeJabba(trajectories,years[1])
  save(kb,file=paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_trajectories.Rdata"))
  
}


#------------------
# Goodness-of-Fit
#------------------
DIC =round(mod$BUGSoutput$DIC,1)

  # get residuals
  Resids = NULL
  for(i in 1:n.indices){
    Resids =rbind(Resids,log(CPUE[,i])-log(apply(posteriors$CPUE[,,i],2,quantile,c(0.5))))   
  }
  
  # Standardized Residuals
  StResid = NULL
  for(i in 1:n.indices){
    StResid =rbind(StResid,log(CPUE[,i]/apply(posteriors$CPUE[,,i],2,quantile,c(0.5)))/
                     apply(posteriors$TOE[,,i],2,quantile,c(0.5))+0.5*apply(posteriors$TOE[,,i],2,quantile,c(0.5)))        
  }
  
  Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
  DF = Nobs-npar
  RMSE = round(100*sqrt(sum(Resids^2,na.rm =TRUE)/DF),1)
  SDNR = round(sqrt(sum(StResid^2,na.rm =TRUE)/(Nobs-1)),2)
  Crit.value = (qchisq(.95, df=(Nobs-1))/(Nobs-1))^0.5
  # Produce statistice describing the Goodness of the Fit

GOF = data.frame(Statistic = c("N","p","DF","SDNR","RMSE","DIC"),Value = c(Nobs,npar,DF,SDNR,RMSE,DIC))
write.csv(GOF,paste0(output.dir,"/GoodnessFit_",assessment,"_",Scenario,".csv"))

# Save Obs,Fit,Residuals
jabba.res = NULL
for(i in 1:n.indices){
  
  Yr = years
  Yr = min(Yr):max(Yr)
  yr = Yr-min(years)+1
  
  exp.i = apply(posteriors$CPUE[,is.na(cpue.raw[,i+1])==F,i],2,quantile,c(0.5))
  obs.i = cpue.raw[is.na(cpue.raw[,i+1])==F,i+1]
  sigma.obs.i = (apply(posteriors$TOE[,is.na(cpue.raw[,i+1])==F,i],2,quantile,c(0.5)))
  
  yr.i = Yr[is.na(cpue.raw[,i+1])==F]
  jabba.res = rbind(jabba.res,data.frame(scenario=Scenario,name=names(cpue)[i+1],year=yr.i,obs=obs.i,hat=exp.i,sigma.obs=sigma.obs.i,residual=log(obs.i)-log(exp.i),tails=tails))  
}

#----------------
# Total Landings
#----------------
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

png(file = paste0(output.dir,"/Landings_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

cord.x <- c(years,rev(years))
y<-rep(0,length(years))
plot(years,(TC),type="l",ylim=c(0,max(TC)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ",catch.metric),main="")
polygon(cord.x,c(TC,rev(y)),col="gray",border=1,lty=1)
dev.off()

#------------------------------
# Plot Posteriors
#------------------------------


sel.par = c(1,8,7,4,3,6,5)

out=data.frame(posteriors[params[sel.par]])
if(nSel>1) out=out[,-c(3:(3+nSel-2))]

node_id = names(out)


#Posteriors
Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Posteriors_",assessment,"_",Scenario,".png"),width  = 8, height = 2.5*round(length(node_id)/3,0), 
    res = 200, units = "in")
par(Par)

#par(mfrow=c(4,2),oma=c(0,1,1,0), mar=c(4,4,1,1))

for(i in 1:length(node_id)){
  
  post.par = as.numeric(unlist(out[paste(node_id[i])]))
  if(unlist(strsplit(node_id[i],"[.]"))[1]=="tau2"){
    post.par=sqrt(post.par)  
    node_id[i] =ifelse(length(unique(sets.var))==1,paste0("sigma.obs"),paste0("sigma.obs",unlist(strsplit(node_id[i],"[.]"))[2]))
  }
  
  
  if(i==1){
    
    rpr = rlnorm(10000,log(SB0.pr[1]),SB0.pr[2]) 
    pdf = stats::density(post.par,adjust=2)  
    prior = dlnorm(sort(rpr),log(SB0.pr[1]),SB0.pr[2])   
    plot(pdf,type="l",ylim=range(prior,pdf$y),xlim=range(c(pdf$x,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=expression(SB[0]),ylab="",xaxs="i",yaxs="i",main="")
    
    polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
    polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
    legend('right',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
    PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
    PPVM = round(mean(post.par)/mean(rpr),3)
    legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")  
    
  }  
  
  if(i<4 & i>1){
    
    prior = stats::density(Hmsy_m_devs[,i-1],adjust=2)
    pdf = stats::density(post.par,adjust=2)  
    
    plot(pdf,type="l",ylim=c(min(prior$y,pdf$y),max(prior$y,pdf$y)), xlim=range(c(pdf$x,quantile(prior$x,c(0.0001,0.96)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
    
    polygon(c(prior$x,rev(prior$x)),c(prior$y,rep(0,length(prior$y))),col=gray(0.4,1))
    polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
    PPVR = round((sd(post.par)/mean(post.par))^2/(sd(Hmsy_m_devs[,i-1])/mean(Hmsy_m_devs[,i-1]))^2,3)  
    PPVM = round(mean(post.par)/mean(Hmsy_m_devs[,i-1]),3)
    legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
  }        
  
  if(i==4){
    if(psi.prior=="beta"){
      #parm = fitdist(post.par[post.par<1 & post.par>0.01], "beta")$estimate
      rpr = rbeta(10000,(psi.pr[1]),psi.pr[2]) 
      pdf = stats::density(post.par)  
      prior = stats::dbeta(sort(rpr),psi.pr[1],psi.pr[2])   
    } else {
      rpr = rlnorm(10000,log(Psi.prior[1]),Psi.prior[2]) 
      pdf = stats::density(post.par)  
      prior = dlnorm(sort(rpr),log(Psi.prior[1])-0.5*Psi.prior[2]^2,Psi.prior[2])
    }
    plot(pdf,type="l",ylim=range(quantile(c(prior,pdf$y,c(0,0.95)))),xlim=range(c(0.5,post.par,pdf$x,rpr)),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
    polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
    polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
    PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
    PPVM = round(mean(post.par)/mean(rpr),3)
    legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")  
  }        
  
  if(i>4){
    if(i<length(node_id)){
      pdf = stats::density((post.par),adjust=2)  
      plot(pdf,type="l",xlim=range(0,(post.par)),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
    }
    
    if(i==length(node_id)){
      if(sigma.proc!=TRUE & i==length(node_id)) {
        plot(1,1,type="n",xlim=range(0,0.15^2),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")  
        abline(v=sigma.proc^2,lwd=2)
      } else {
        if(proc.type=="igamma"){
          if(i==length(node_id)& pr.proc[1]>0.9){
            rpr = sqrt(1/rgamma(10000,pr.proc[1],pr.proc[2]))
            prior = stats::density(rpr,adjust=2)
            pdf = stats::density(sqrt(post.par),adjust=2)  
            plot(pdf,type="l",ylim=c(0,max(prior$y,pdf$y)),xlim=range(0,sqrt(post.par)),yaxt="n",xlab="sigma.proc",ylab="",xaxs="i",yaxs="i",main="")
            polygon(c(prior$x,rev(prior$x)),c(prior$y,rep(0,length(prior$y))),col=gray(0.4,1))
            polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
            PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
            PPVM = round(mean(sqrt(post.par))/mean(sqrt(rpr)),3)
            legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")  
            
          }else{
            pdf = stats::density(sqrt(post.par),adjust=2)  
            plot(pdf,type="l",ylim=c(0,max(pdf$y)),xlim=range(0,sqrt(post.par)),yaxt="n",xlab="sigma.proc",ylab="",xaxs="i",yaxs="i",main="")
            polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
            
            
          }} 
        if(proc.type=="lnorm"){ 
          rpr = rlnorm(10000,pr.proc[1],pr.proc[2])
          pdf = stats::density(sqrt(post.par),adjust=2)  
          prior = stats::density(rpr,adjust=2)
          plot(pdf,type="l",ylim=c(0,max(prior$y,pdf$y)),xlim=range(0,sqrt(post.par)),yaxt="n",xlab="sigma.proc",ylab="",xaxs="i",yaxs="i",main="")
          polygon(c(prior$x,rev(prior$x)),c(prior$y,rep(0,length(prior$y))),col=gray(0.4,1))
          polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
          
          PPVR = round((sd(log(post.par)))^2/(sd(log(rpr)))^2,3)  
          PPVM = round(mean((post.par))/mean((rpr)),3)
          legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")  
          
          
        }
      }
    }}
  
  
}          

mtext(paste("Density"), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
dev.off() 

out=data.frame(posteriors[params[sel.par]])
if(nSel>1) out=out[,-c(3:(3+nSel-2))]
node_id = names(out)

#-----------------------------
# MCMC chains of posteriors
#-----------------------------

Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/MCMC_",assessment,"_",Scenario,".png"), width = 8, height = 2.5*round(length(node_id)/3,0), 
    res = 200, units = "in")
par(Par)
for(i in 1:length(node_id)){
  
  cols=rainbow(10)
  post.par = as.numeric(unlist(out[paste(node_id[i])]))
  plot(out[,i],xlab=paste(node_id[i]),ylab="",type="l",col=cols[7])
  lines(rep(mean(out[,i]),length(out[,i])),col=2,lwd=2)   
}
dev.off()

#-----------------------------------------------------------
# <><<><<><<>< Produce JABBA Model Diagnostics ><>><>><>><>
#-----------------------------------------------------------
cat(paste0("\n","><> Producing JABBA Model Fit Diagnostics <><","\n"))


# extract predicted CPUE + CIs

N = n.years
series = 1:n.indices

check.yrs = apply(CPUE,1,sum,na.rm=TRUE)
cpue.yrs = years[check.yrs>0]

#CPUE FITS
Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Fits_",assessment,"_",Scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
    res = 200, units = "in")
par(Par)


for(i in 1:n.indices){
  
  # set observed vs predicted CPUE
  #par(mfrow=c(1,1))
  Yr = years
  Yr = min(Yr):max(Yr)
  yr = Yr-min(years)+1
  
  fit = apply(posteriors$CPUE[,,i],2,quantile,c(0.025,0.5,0.975))
  fit.hat = apply(posteriors$Ihat[,,i],2,quantile,c(0.025,0.5,0.975))
  mufit = mean(fit[2,])
  fit = fit/mufit
  fit.hat = fit.hat/mufit
  
  cpue.i = CPUE[is.na(CPUE[,i])==F,i]
  yr.i = Yr[is.na(CPUE[,i])==F]
  se.i = sqrt(se2[is.na(CPUE[,i])==F,(i)])
  
  ylim = c(min(fit*0.9,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit*1.05,exp(log(cpue.i)+1.96*se.i)/mufit))
  
  cord.x <- c(Yr,rev(Yr))
  cord.y <- c(fit[1,yr],rev(fit[3,yr]))
  cord.yhat <- c(fit.hat[1,yr],rev(fit.hat[3,yr]))
  # Plot Observed vs predicted CPUE
  plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(years),type='n',xaxt="n",yaxt="n")
  axis(1,labels=TRUE,cex=0.8)
  axis(2,labels=TRUE,cex=0.8)
  polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
  polygon(cord.x,cord.yhat,col=grey(0.3,0.5),border=grey(0.3,0.5),lty=2)
  
  lines(Yr,fit[2,yr],lwd=2,col=1)
  if(SE.I ==TRUE | max(se2)>0.01){ plotCI(yr.i,cpue.i/mufit,ui=exp(log(cpue.i)+1.96*se.i)/mufit,li=exp(log(cpue.i)-1.96*se.i)/mufit,add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
    points(yr.i,cpue.i/mufit,pch=21,xaxt="n",yaxt="n",bg="white")}
  
  legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
}

mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
dev.off()

#log CPUE FITS
Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/logFits_",assessment,"_",Scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
    res = 200, units = "in")
par(Par)


for(i in 1:n.indices){
  
  # set observed vs predicted CPUE
  #par(mfrow=c(1,1))
  Yr = years
  Yr = min(Yr):max(Yr)
  yr = Yr-min(years)+1
  
  fit = apply(posteriors$CPUE[,,i],2,quantile,c(0.025,0.5,0.975))
  mufit = mean(fit[2,])
  fit = fit/mufit
  cpue.i = CPUE[is.na(CPUE[,i])==F,i]
  yr.i = Yr[is.na(CPUE[,i])==F]
  se.i = sqrt(se2[is.na(CPUE[,i])==F,(i)])
  
  ylim = log(c(min(fit[,yr[is.na(CPUE[,i])==F]]*0.8,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit[,yr[is.na(CPUE[,i])==F]]*1.3,exp(log(cpue.i)+1.96*se.i)/mufit)))
  
  cord.x <- c(Yr,rev(Yr))
  cord.y <- log(c(fit[1,yr],rev(fit[3,yr])))
  
  # Plot Observed vs predicted CPUE
  plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
  axis(1,labels=TRUE,cex=0.8)
  axis(2,labels=TRUE,cex=0.8)
  #polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
  
  
  lines(Yr,log(fit[2,yr]),lwd=2,col=4)
  if(SE.I ==TRUE | max(se2)>0.01){ plotCI(yr.i,log(cpue.i/mufit),ui=log(exp(log(cpue.i)+1.96*se.i)/mufit),li=log(exp(log(cpue.i)-1.96*se.i)/mufit),add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
    points(yr.i,log(cpue.i/mufit),pch=21,xaxt="n",yaxt="n",bg="white")}
  legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
}

mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
mtext(paste("Log Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
dev.off()



# JABBA-residual plot
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Residuals_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

plot(Yr,Yr,type = "n",ylim=ifelse(rep(max(Resids,na.rm = T),2)>0.9,range(1.2*Resids,na.rm = T),range(c(-1.3,1.2))),xlim=range(cpue.yrs),ylab="log residuals",xlab="Year")
boxplot(Resids,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
abline(h=0,lty=2)
positions=runif(n.indices,-0.2,0.2)

for(i in 1:n.indices){
  for(t in 1:n.years){
    lines(rep((Yr+positions[i])[t],2),c(0,Resids[i,t]),col=jabba.colors[i])}
  points(Yr+positions[i],Resids[i,],col=1,pch=21,bg=jabba.colors[i])}
mean.res = apply(Resids,2,mean,na.rm =TRUE)
smooth.res = predict(loess(mean.res~Yr),data.frame(Yr=cpue.yrs))
lines(cpue.yrs,smooth.res,lwd=2)
DIC =round(mod$BUGSoutput$DIC,1)
# get degree of freedom
Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
DF = Nobs-npar

RMSE = round(100*sqrt(sum(Resids^2,na.rm =TRUE)/DF),1)

legend('topright',c(paste0("RMSE = ",RMSE,"%")),bty="n")
legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,pt.cex=1.1,cex=0.75,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba.colors[series],1),lwd=c(rep(-1,n.indices),2))

dev.off()

#Save Residuals 
Res.CPUE = data.frame(Resids)
row.names(Res.CPUE) = indices   
colnames(Res.CPUE) = paste(Yr)
write.csv(Res.CPUE,paste0(output.dir,"/ResCPUE_",assessment,"_",Scenario,".csv"))

#---------------------------------------
# Stadardized Residuals
#--------------------------------------
Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/StandardizedResids_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)


plot(Yr,Yr,type = "n",ylim=c(min(-1,-1.2*max(abs(StResid),na.rm = T)),max(1,1.2*max(abs(StResid),na.rm = T))),xlim=range(cpue.yrs),ylab="Standardized residuals",xlab="Year")
boxplot(StResid,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
abline(h=0,lty=2)
positions=runif(n.indices,-0.2,0.2)

for(i in 1:n.indices){
  for(t in 1:n.years){
    lines(rep((Yr+positions[i])[t],2),c(0,StResid[i,t]),col=jabba.colors[i])}
  points(Yr+positions[i],StResid[i,],col=1,pch=21,bg=jabba.colors[i])}
mean.res = apply(StResid,2,mean,na.rm =TRUE)
smooth.res = predict(loess(mean.res~Yr),data.frame(Yr=cpue.yrs))
lines(cpue.yrs,smooth.res,lwd=2)
DIC =round(mod$BUGSoutput$DIC,1)
SDNR = round(sqrt(sum(StResid^2,na.rm =TRUE)/(Nobs-1)),2)
Crit.value = (qchisq(.95, df=(Nobs-1))/(Nobs-1))^0.5
legend('topright',c(paste0("SDNR = ",SDNR,"(",round(Crit.value,2),")")),bty="n")
legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,cex=0.75,pt.cex=1.1,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba.colors[series],1),lwd=c(rep(-1,n.indices),2))


dev.off()

#Save standardized Residuals 
StRes.CPUE = data.frame(StResid)
row.names(Res.CPUE) = indices   
colnames(Res.CPUE) = paste(Yr)
write.csv(Res.CPUE,paste0(output.dir,"/StResCPUE_",assessment,"_",Scenario,".csv"))


#------------------------------
# Plot process error deviation
#------------------------------

proc.dev = apply(posteriors$Proc.Dev,2,quantile,c(0.025,0.5,0.975))

Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)

png(file = paste0(output.dir,"/ProcDev_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

ylim = c(min(-0.22,proc.dev),max(0.22,proc.dev))#range(proc.dev)*1.1
cord.x <- c(years,rev(years))
cord.y <- c(proc.dev[1,],rev(proc.dev[3,]))
# Process Error
plot(years,proc.dev[2,],ylab="Process Error Deviates",xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,proc.dev[2,],lwd=2)
lines(years,rep(0,length(years)),lty=5)


dev.off()

procE.dev = data.frame(Scenario,Yr=years,mu=proc.dev[2,],lci=proc.dev[1,],uci=proc.dev[3,])

#-----------------------------------------------------------
# <><<><<><<><<>< JABBA Management Plots ><>><>><>><>><>><>
#-----------------------------------------------------------
cat(paste0("\n","><> Producing Management Plots <><","\n"))


#-----------------------
# Plot Biomass B_t
#-----------------------

B_t = posteriors$SB
mu.B = apply(B_t,2,quantile,c(0.025,0.5,0.975))

Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Biomass_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
    res = 200, units = "in")
par(Par)

ylim = c(0, max(mu.B ))
cord.x <- c(years,rev(years))
cord.y <- c(mu.B [1,],rev(mu.B [3,]))

# B_t
plot(years,mu.B[2,],ylab=paste0("Spawning Biomass ",catch.metric),xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,mu.B[2,],lwd=2,col=1)
lines(years,rep(mean(posteriors$SBmsy),length(years)),lty=5)
text((max(years)-min(years))/30+years[1],mean(posteriors$SBmsy)*1.11,expression(paste(B[MSY])))
dev.off()



#-----------------------
# Plot B/Bmsy and H/Hmsy
#-----------------------

HtoHmsy = posteriors$HtoHmsy
BtoBmsy = posteriors$BtoBmsy

mu.f = apply(HtoHmsy,2,quantile,c(0.025,0.5,0.975))
mu.b = apply(BtoBmsy,2,quantile,c(0.025,0.5,0.975))

f = HtoHmsy[,N]
b = BtoBmsy[,N]

Par = list(mfrow=c(1,2),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/TrendMSY_",assessment,"_",Scenario,".png"), width = 7, height = 3, 
    res = 200, units = "in")
par(Par)

ylim = c(0, max(mu.f))
cord.x <- c(years,rev(years))
cord.y <- c(mu.f[1,],rev(mu.f[3,]))

# H/Hmsy
plot(years,mu.f[2,],ylab=bquote(H/H[.(refB)]),xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,mu.f[2,],lwd=2,col=1)
lines(years,rep(1,length(years)),lty=5)

ylim = c(0, max(mu.b,1.1))

cord.x <- c(years,rev(years))
cord.y <- c(mu.b[1,],rev(mu.b[3,]))
plot(years,mu.b[2,],ylab=bquote(SB/SB[.(refB)]),xlab="Year",ylim=ylim,type="n")
polygon(cord.x,cord.y,col='grey',border=0,lty=2)
lines(years,mu.b[2,],lwd=2,col=1)
lines(years,rep(1,length(years)),lty=5)
dev.off()

  #-----------------------------------------
  # Produce JABBA SP-phase plot
  #-----------------------------------------
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/SPphase_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  m = median(posteriors$m)
  Bit = seq(1,median(posteriors$SB0),median(posteriors$SB0)/500)
  Cmsy = Bit*median(posteriors$Hmsy)
  B = apply(posteriors$SB,2,mean)
  Hmsy.sp = median(posteriors$Hmsy.y[,n.years]) 
  SB0.sp = median(posteriors$SB0)
  
  SP = Hmsy.sp/(1-1/m)*Bit*(1-(Bit/SB0.sp)^(m-1))
  
  #SP = median(posteriors$r)/(m-1)*Bit*(1-(Bit/median(posteriors$SB0))^(m-1))  
  Bmsy.sp = median(posteriors$SBmsy)
  MSY.sp = quantile(posteriors$SBmsy*posteriors$Hmsy.y[,n.years],c(0.025,0.5,0.975))
  green.x = c(max(Bit,B),max(Bit,B),Bmsy.sp,Bmsy.sp,max(Bit))
  green.y = c(Bmsy.sp,0,0,max(SP),max(Cmsy))
  
  
  red.x = c(0,0,Bmsy.sp,Bmsy.sp,0)
  red.y = c(SB0.sp,0,max(SP),SB0.sp,SB0.sp)
  plot(Bit,SP,type = "n",ylim=c(0,max(c(max(TC,na.rm=T)*1.05,max(MSY.sp*1.1)))),xlim=c(0,max(Bit,B)),ylab="Surplus Production (t)",xlab="Spawning Biomass (t)",xaxs="i",yaxs="i")
  rect(0,0,SB0.sp*1.1,SB0.sp*1.1,col="green",border=0)
  rect(0,0,SB0.sp,SB0.sp,col="yellow",border=0)
  if(KOBE.type!="ICCAT") rect(0,max(SP),SB0.sp,SB0.sp,col="orange",border=0)
  polygon(green.x,green.y,border = 0,col="green")
  polygon(red.x,red.y,border = 0,col="red")
  
  ry.sp = Bit[Bit<=Bmsy.sp]
  for(i in 1:length(ry.sp)){
    
    lines(rep(Bit[i],2),c(Cmsy[i],SP[i]),col=ifelse(i %% 2== 0,"yellow","red"),lty=3)  
    #i = i+1
  }
  
  gy.sp = Bit[Bit>Bmsy.sp]
  for(i in (length(ry.sp)+1):length(Bit)){
    
    #lines(rep(Bit[i],2),c(max(SP),Cmsy[i]),col=ifelse(i %% 2== 0,ifelse(KOBE.type=="ICCAT","yellow","orange"),"green"),lty=3)  
    #i = i+1
  }
  
  
  polygon(c(-10000,10^7,10^7,-10000),c(rep(MSY.sp[1],2),rep(MSY.sp[3],2)),border = FALSE,col=rgb(0,0,1,0.4))
  lines(Bit,SP,col=4,lwd=2)
  lines(B,apply(Catch,1,sum),lty=1,lwd=1)
  points(B,apply(Catch,1,sum),cex=0.8,pch=16)
  #lines(Bit,Cmsy,col=1,lwd=1,lty=2)
  N=n.years
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(B[sel.yr],apply(Catch,1,sum)[sel.yr],col= 1,pch=c(22,21,24),bg="white",cex=1.7)
  abline(h=max(SP),col=4,lty=5)
  sel.years =years[sel.yr]
  lines(rep(median(posteriors$SBmsy),2),c(-1000,max(SP)),lty=2,col=4)
  
  legend('topright', 
         c(expression(SB[MSY]),"MSY","SP","Catch",paste(sel.years)), 
         lty=c(2,5,1,1,1,1,1),pch=c(-1,-1,-1,16,22,21,24),pt.bg=c(0,0,0,0,rep("white",3)), 
         col=c(4,4,4,rep(1,4)),lwd=c(1,1,2,1,1,1),cex=0.8,pt.cex=c(-1,-1,-1,0.5,rep(1.3,3)),bty="n")
  
  dev.off()

  # save results
  SPphase = data.frame(Scenario,SB_i=round(Bit,1),SP=round(SP,1),Hmsy=round(Hmsy.sp,3),r=round(Hmsy.sp*(m-1)/(1-1/m),3),m=round(m,3),MSY=round(as.numeric(MSY.sp[2]),1),SB0=round(SB0.sp,1),Cmsy=round(Cmsy,1))
  




######################################
# Prepare Kobe
######################################
# extract vectors BtoBmsy and FtoFmsy

if(KOBE.plot==TRUE){
  # prepare biplot
  mu.y = apply(HtoHmsy,2,quantile,c(0.5))
  mu.x = apply(BtoBmsy,2,quantile,c(0.5))
  y = HtoHmsy[,N]
  x = BtoBmsy[,N]
  
  
  f<-HtoHmsy[,N]
  b<-BtoBmsy[,N]
  # fit kernel function
  # fit kernel function
  kernelF <- ci2d(b,f,nbins=151,factor=1.5,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,xlab= ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])))
  
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Kobe_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  #Create plot
  plot(1000,1000,type="b", xlim=c(0,max(1/SBmsySB0,mu.b[2,]) +0.05), ylim=c(0,max(apply(HtoHmsy,2,quantile,c(0.5)),quantile(f,0.85),2.)),lty=3,ylab=bquote(H/H[.(refB)]),xlab=bquote(SB/SB[.(refB)]),xaxs="i",yaxs="i")
  c1 <- c(-1,100)
  c2 <- c(1,1)
  
  # extract interval information from ci2d object
  # and fill areas using the polygon function
  zb2 = c(0,1)
  zf2  = c(1,100)
  zb1 = c(1,100)
  zf1  = c(0,1)
  polygon(c(zb1,rev(zb1)),c(0,0,1,1),col="green",border=0)
  polygon(c(zb2,rev(zb2)),c(0,0,1,1),col="yellow",border=0)
  polygon(c(1,100,100,1),c(1,1,100,100),col=ifelse(KOBE.type=="ICCAT","yellow","orange"),border=0)
  polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)
  
  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
  points(mu.b[2,],mu.f[2,],pch=16,cex=1)
  
  
  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(mu.b[2,],mu.f[2,], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.b[2,sel.yr],mu.f[2,sel.yr],col=
           1,pch=c(22,21,24),bg="white",cex=1.9)
  
  # Get Propability
  Pr.green = sum(ifelse(b>1 & f<1,1,0))/length(b)*100
  Pr.red = sum(ifelse(b<1 & f>1,1,0))/length(b)*100
  
  if(KOBE.type=="ICCAT"){               
    Pr.yellow = (sum(ifelse(b<1 & f<1,1,0))+sum(ifelse(b>1 & f>1,1,0)))/length(b)*100} else {
      Pr.yellow = sum(ifelse(b<1 & f<1,1,0))/length(b)*100
      Pr.orange = sum(ifelse(b>1 & f>1,1,0))/length(b)*100
    }
  
  
  sel.years = c(years[sel.yr])
  ## Add legend
  if(KOBE.type=="ICCAT"){
    legend('topright', 
           c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.green),1),"%")), 
           lty=c(1,1,1,rep(-1,7)),pch=c(22,21,24,rep(22,7)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","green"), 
           col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.1,3)),bty="n")
  }else{
    legend('topright', 
           c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.orange,Pr.green),1),"%")), 
           lty=c(1,1,1,rep(-1,8)),pch=c(22,21,24,rep(22,8)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","orange","green"), 
           col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.2,4)),bty="n")  
    
  }
  dev.off()
}



if(Biplot==TRUE){
  
  #---------------------------------------------------------
  # Produce 'post-modern' biplot (see Quinn and Collie 2005)
  #---------------------------------------------------------
  
  # read ftarget,bthreshold
  ftarget<-0.8
  bthreshold<-0.2
  
  # fit kernel function
  kernelF <- ci2d(f,b,nbins=151,factor=2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,ylab= ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])))
  
  
  Par = list(mfrow=c(1,1),mai=c(0.2,0.15,0,.15),omi = c(0.3,0.25,0.2,0) + 0.1, mgp =c(3,1,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Biplot_",assessment,"_",Scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  #Create plot
  plot(1000,1000,type="b", ylim=c(0,1/SBmsySB0+0.05), xlim=c(0,max(apply(HtoHmsy,2,quantile,c(0.5)),quantile(f,0.85),2.)),lty=3,xaxs="i",yaxs="i")
  
  # and fill areas using the polygon function
  fint = seq(0.001,100,0.01)
  #Zone X
  xb=bthreshold+(1.0-bthreshold)/ftarget*fint
  xf =  ifelse(xb>1,0.8,fint)
  polygon(c(0,0,xf),c(max(xb),bthreshold,xb),col="green")
  zb = bthreshold+(1.0-bthreshold)*fint
  zf  = ifelse(zb>1,1,fint) 
  polygon(c(zf,rep(max(fint),2),rep(0,2)),c(zb,max(zb),0,0,bthreshold),col="red")
  
  polygon(c(xf,rev(zf)),c(xb,rev(zb)),col="yellow")
  
  c1 <- c(-1,100)
  c2 <- c(1,1)
  
  # extract interval information from ci2d object
  # and fill areas using the polygon function
  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
  points(mu.f[2,],mu.b[2,],pch=16,cex=1)
  
  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(mu.f[2,],mu.b[2,], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.f[2,sel.yr],mu.b[2,sel.yr],col=
           1,pch=c(22,21,24),bg="white",cex=1.9)
  
  
  sel.years = years[sel.yr]
  ## Add legend
  legend('topright', 
         c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I."), 
         lty=c(1,1,1,-1,-1,-1),pch=c(22,21,24,22,22,22),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4"), 
         col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,4),1.7,1.7,1.7),bty="n")
  
  
  
  Zone  = NULL
  Status = NULL
  X  = 0.15
  Y = 0
  Z = -0.15
  
  for(i  in 1:length(f))
  {
    if(b[i]>1.0){
      if(f[i]<ftarget){
        Zone[i]<-X
      } else if (f[i]>1.0){
        Zone[i]<-Z
      } else {
        Zone[i]<-Y
      }
    } else {
      if(b[i]>bthreshold+(1.0-bthreshold)/ftarget*f[i]){
        Zone[i]<-X
      } else if(b[i]<bthreshold+(1.0-bthreshold)*f[i]){
        Zone[i]<-Z
      } else {
        Zone[i]<-Y
      }
      
      
    }}
  
  perGreen = round(length(Zone[Zone==0.15])/length(Zone)*100,1) 
  perYellow = round(length(Zone[Zone==0])/length(Zone)*100,1) 
  perRed = round(length(Zone[Zone==-0.15])/length(Zone)*100,1)
  
  mtext(bquote(SB/SB[.(refB)]), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
  mtext(bquote(H/H[.(refB)]), side=1, outer=TRUE, at=0.5,line=1,cex=0.9)
  
  text(0.65,1/SBmsySB0,paste0(perGreen,"%"))
  text(0.9,1/SBmsySB0,paste0(perYellow,"%"))
  text(1.2,1/SBmsySB0,paste0(perRed,"%"))
  
  dev.off()
  
}



#--------------------------------------
# Plot projections and safe posteriors
#--------------------------------------
if(Projection ==TRUE){
  cat(paste0("\n","><> Producing Future TAC Projections <><","\n"))
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Projections_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  
  proj.yrs =  years[n.years]:(years[n.years]+pyrs)
  # Dims 1: saved MCMC,2: Years, 3:alternatic TACs, 4: P, H/Hmsy, B/Bmsy
  projections = array(NA,c(nsaved,length(proj.yrs),nTAC,3))
  for(i in 1:nTAC){
    projections[,,i,1] = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
  }
  
  for(i in 1:nTAC){
    projections[,,i,2] = cbind(posteriors$BtoBmsy[,(n.years):n.years],posteriors$prBtoBmsy[,,i])
  }
  for(i in 1:nTAC){
    projections[,,i,3] = cbind(posteriors$HtoHmsy[,(n.years):n.years],posteriors$prHtoHmsy[,,i])
  }
  
  kjp = kobeJabbaProj(projections,proj.yrs[1])
  
  save(kjp,file=paste0(output.dir,"/",assessment,"_",Mod.names,"_",Scenario,"_projections.Rdata"))
  
  # Change here for ICCAT Bmsy plot
  Traj = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,nTAC])
  
  #plot(proj.yrs,apply(Traj,2,mean),ylim=c(0,1),xlim=c(min(proj.yrs),max(proj.yrs)+length(proj.yrs)*0.2),type="n",ylab="Biomass depletion (B/K)",xlab="Projection Years")
  plot(proj.yrs,apply(Traj,2,mean),ylim=c(0,1),xlim=c(min(proj.yrs),max(proj.yrs)),type="n",ylab="Biomass depletion (B/K)",xlab="Projection Years")
  
  cols = rev(seq(0.4,0.9,0.5/nTAC))
  plot.order = (1:nTAC)
  for(j in 1:(nTAC)){
    i =  plot.order[j]
    Traj = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
    polygon(c(proj.yrs,rev(proj.yrs)),c(apply(Traj,2,quantile,0.05),rev(apply(Traj,2,quantile,0.95))),col=grey(cols[i],1),border=NA)
  }
  for(j in 1:(nTAC)){
    i =  plot.order[j]
    Traj = cbind(posteriors$P[,(n.years):n.years],posteriors$prP[,,i])
    lines(proj.yrs,apply(Traj,2,median),col=rev(jabba.colors[j]),lwd=2)
  }
  lines(proj.yrs[1:(imp.yr-yr.last+1)],apply(Traj,2,median)[1:(imp.yr-yr.last+1)],col=1,lwd=2)
  BmsyK=(m)^(-1/(m-1))
  abline(h=BmsyK,lty=2,lwd=2)
  legend("topleft",paste(TACs,"(t)"),col=(jabba.colors[1:nTAC]),lwd=2,cex=0.8)   
  
  dev.off()    
  
}
#---------------------------------------------------------------------------------
# Save core results
#---------------------------------------------------------------------------------
rownames(Stock_trj) = 1:nrow(Stock_trj)
Stock.trj=data.frame(Scenario,Yr=years,Total.Catch=TC,data.frame(Stock_trj))
jabba_out = list(Data=surplus.dat,Pars=results,Estimates=Table,Stock.trj=Stock.trj,Fits = jabba.res,ProcErrDev=procE.dev, SPphase=SPphase)
save(jabba_out, file=paste0(output.dir,"/jabbaout_",assessment,"_",Scenario,".RData"))  
cat(paste0("\n","Scenario ",Mod.names,"_",Scenario," - DONE!","\n"))
#---------------------------------------------------------------------------------
#compile jabba2FRL
if(jabba2FRL==TRUE){
#---------------------------------------------------------------------------------
# Define elements
timeseries = data.frame(factor=assessment,level=Scenario,year=years,area=1,season=1,biomass=NA,ssb=Stock.trj$SBt.Median,rec=NA,catch=TC,proc.dev=proc.dev[2,])
refpts = data.frame(factor=assessment,level=Scenario,quant = c("hat","var"), k=c(median(posteriors$SB0),var(posteriors$SB0)),bmsy=c(median(posteriors$SBmsy),var(posteriors$SBmsy)),
                    fmsy=c(median(posteriors$Hmsy.y[,length(years)]),var(posteriors$Hmsy.y[,length(years)])),msy=c(median(posteriors$SBmsy*posteriors$Hmsy.y[,length(years)]),var(posteriors$SBmsy*posteriors$Hmsy.y[,length(years)]))) 
                                                                             
pfunc = data.frame(factor=assessment,level=Scenario, k=median(posteriors$SB0),r=median(posteriors$Hmsy.y[,length(years)])*(m-1)/(1-1/m),p = m-1,shape=median(posteriors$SBmsy)/median(posteriors$SB0),m)          
curves=data.frame(factor=assessment,level=Scenario,ssb=SPphase$SB_i,yield=SPphase$SP)                   
dgs = data.frame(factor=assessment,level=Scenario,name=jabba.res$name,year=jabba.res$year,season=1,obs=jabba.res$obs,hat=jabba.res$hat,residual=jabba.res$residual,tails=jabba.res$tails)
jb = list(timeseries=timeseries,refpts=refpts,pfunc=pfunc,ts2=NULL,curves=curves,dgs=dgs)
save(jb,file=paste0(output.dir,"/jb.",Scenario,".Rdata"))
}
#---------------------------------------------------------------------------------

# END OF CODE