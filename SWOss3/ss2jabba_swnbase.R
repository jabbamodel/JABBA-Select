#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# ss2jabba extracts all required JABBA-Select input data from a ss3 repile
# and create the ss4js rdata object 
# Written by Henning Winker, 2019
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

#---------------------------------------------------------------------
# Set Working directory file where to store the results
File = "C:/Work/Research/GitHub/JABBA-SELECT"
# Set Assessment
assessment = "SWOss3"

# Install required packages if missing
list.of.packages <- c("reshape2","r4ss")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load Packages
library(reshape2);library(r4ss)



#----------------------------------------
# Extract basic data from swnbase
#----------------------------------------

load(paste0(File,"/",assessment,"/swnbase.Rdata"),verbose=TRUE)
# Rename
ss3rep = swnbase

# Extract/Plot catch matrix
Par = list(mfrow=c(1,1),mai=c(0.5,0.5,0.1,.1),omi = c(0.1,0.1,0.1,0.1) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/",assessment,"_catch.png"), width = 6, height =4, 
    res = 200, units = "in")
par(Par)
catches = SSplotCatch(ss3rep,subplots = 2)$totcatchmat
catches = catches[,c(ncol(catches),1:(ncol(catches)-1))]
dev.off()

# Stock parameters by sex
#Females
biodat = ss3rep$endgrowth[swnbase$endgrowth$Sex==1,]
vbgf_F = coef(nls(Len_Beg~Linf*(1-exp(-K*(Age-t0))),start=list(Linf=max(biodat$Len_Beg),K=0.2,t0=-0.5),data=biodat))
LW_F = lm(log(Wt_Beg)~log(Len_Beg),biodat)$coef
Lm50 = biodat$Len_Beg[biodat$Age_Mat==0.5]   
Lm95 = Lm50*1.05          
# Males =
biodat = swnbase$endgrowth[swnbase$endgrowth$Sex==2,]
vbgf_M = coef(nls(Len_Beg~Linf*(1-exp(-K*(Age-t0))),start=list(Linf=max(biodat$Len_Beg),K=0.2,t0=-0.5),data=biodat))
LW_M = lm(log(Wt_Beg)~log(Len_Beg),biodat)$coef

# summarized stock parameters
stock.pars = list(vbgf_F=vbgf_F,vbgf_M=vbgf_M,mat=c(Lm50,Lm95),LW_F=LW_F,LW_M=LW_M,Amax=max(biodat$Age))


# Get selectivties 
sels = ss3rep$sizeselex
sfl = unique(sels$Fleet)
Sel = as.numeric(names(sels[,6:ncol(sels)]))
for(i in 1:length(sfl)){
  Sel = cbind(Sel,as.numeric(t(sels[sels$Fleet==sfl[i],6:ncol(sels)][1,])))
}
datS = data.frame(Sel)  
colnames(datS)=  c("L",paste(unlist(ss3rep$definitions[2,-1])))
names(datS)
# remove indices Age-1 to Age-4 (not suitable due to lags)
datS = datS[,-c(13:16)]
# replace selectivity for Age-5 assuming sel = maturity
datS[,13] = 1/(1+exp(-log(19)*(datS$L-Lm50)/(Lm95-Lm50)))

# Find unique selectivities
fleets=names(datS[,-1])
flpos = which(paste(unlist(ss3rep$definitions[2,-1]))%in%fleets)
CPUE.units = ss3rep$survey_units[flpos] # numbers vs biomass 
nL = nrow(datS)
selcheck =  as.numeric(datS[ceiling(nL/3),-1])
selunique = unique(selcheck) 
nSel = length(selunique)
# Create select table
select = NULL 
for(i in 1:nSel){
  col.id = which(selunique[i]==selcheck)
  select = rbind(select,data.frame(Fleet.id=col.id,Fleet=paste(fleets[col.id]),Selectivity=i,CPUE.units=CPUE.units[col.id]))  
}
select= select[order(select$Fleet.id),]
select$Catch = select$Fleet%in%names(catches[,-1])

# Extract unique selectivies
get.sel = NULL
for(i in 1:nSel){
  get.sel =c(get.sel,which(selunique[i]==selcheck)[1])
}  
ss3sel = datS[,c(1,(get.sel+1))]


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

selex.pars = NULL
# convert ss2jabba.selex
for(i in 1:nSel){
  
  # check maxima
  peak = which(max(ss3sel[,i+1])==ss3sel[,i+1])
  # get inits of jabba.selex parameters
  pars = c(SL50	= ss3sel[peak[1],1]*0.8,SL95=ss3sel[peak[1],1]*0.95,SL.desc=ss3sel[max(peak),1],CV.desc=0.2,min.desc=0.001)
  # Minimize
  jsel.est = optim(pars, fn = jsel.ll,method="L-BFGS-B",lower=10^-3,upper=max(ss3sel$L*1.1), dat = ss3sel[,c(1,i+1)], hessian = TRUE)
  #Results
  jsel.out = jabba.selex(jsel.est$par,dat=ss3sel[,c(1,i+1)])$results 
  selex.pars = cbind(selex.pars,round(jsel.est$par,4))
  # Plot
  if(i==1){
    Par = list(mfrow=c(round(nSel/2+0.01,0),ifelse(nSel==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(File,"/",assessment,"/SELEX_",assessment,".png"), width = 7, height = ifelse(nSel==1,5,ifelse(nSel==2,3.,2.5))*round(nSel/2+0.01,0), 
        res = 200, units = "in")
    par(Par)}
  # Plot
  plot(jsel.out$L ,jsel.out$obs,type="p",pch=1,ylab="Selectivity",xlab="Length")
  lines(jsel.out$L,jsel.out$logis,lwd=2,col=4)
  if(jsel.est$par[3]<0.98*max(jsel.out$L))lines(jsel.out$L,jsel.out$halfnorm,lwd=2,col=3)
  if(jsel.est$par[3]<0.98*max(jsel.out$L)) lines(jsel.out$L,jsel.out$height,lwd=2,col=7)
  lines(jsel.out$L,jsel.out$fit,lwd=2,col=2)
  if(i==1)legend("right",c("ss3","Fit","Logistic","Half-Normal","Height"),pch=c(1,rep(-1,4)),lwd=c(-1,rep(2,4)),col=c(1,2,4,3,7),cex=0.8,bty="n")
  if(i == nSel){
    mtext(paste("Length"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
    mtext(paste("Selectivity"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
    dev.off()
  }
} # End of selectivity loop
colnames(selex.pars) = paste0("S",1:nSel)




# Get CPUE
CPUE.agg = aggregate(Obs~Yr+Fleet_name+Fleet, ss3rep$cpue,mean)
fl = unique(CPUE.agg$Fleet)
flname = unique(CPUE.agg$Fleet_name)
nI = length(fl)
# Add early catch years
CPUE.agg = rbind(CPUE.agg,data.frame(Yr=catches[,1],Fleet_name="Z",Fleet=1000,Obs=NA)) 
cpues = dcast(CPUE.agg,Yr~Fleet_name,value.var = "Obs")
# Reorganize
cpues=cpues[,c("Yr",flname)]
cpues = cpues[,which(names(cpues)%in%c("Yr",fleets))]
# Remove SPN (use ESP Age5+)
cpues = cpues[,-2]
Inames = names(cpues)
nI = length(Inames)
Par = list(mfrow=c(round(nI/2+0.01,0),ifelse(nI==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(File,"/",assessment,"/Indices_",assessment,".png"), width = 7, height = ifelse(nI==1,5,ifelse(nI==2,3.,2.5))*round(nI/2+0.01,0), 
    res = 200, units = "in")
par(Par)
for(i in 1:nI){
  plot(Obs~Yr,data=CPUE.agg[CPUE.agg$Fleet==fl[i],],type="l",ylim=c(0,max(CPUE.agg[CPUE.agg$Fleet==fl[i],4])*1.1),col=4,lwd=2)  
  points(Obs~Yr,data=CPUE.agg[CPUE.agg$Fleet==fl[i],],cex=0.8)  
  legend("top",paste0(flname[i]),bty="n",cex=0.9,y.intersp = -0.2)
}
mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
mtext(paste("CPUE"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
dev.off()
# add cpue.se
se.agg = aggregate(SE~Yr+Fleet_name+Fleet, swnbase$cpue,mean)
se.agg = rbind(se.agg,data.frame(Yr=cpues[,1],Fleet_name="Z",Fleet=1000,SE=NA)) 
cpue.se = dcast(se.agg,Yr~Fleet_name,value.var = "SE")
cpue.se = cpue.se[,names(cpues)]

# Add to seltab if Fleet produced CPUE
select$CPUE = select$Fleet%in%Inames[-1]
select$q = 0
select$q[select$CPUE] = 1:length(select$q[select$CPUE]) 
selex = data.frame(Parameter=c("SL50","SL95","SL.desc","CV.desc","min.desc"),selex.pars)
# write catch,cpue,se, select and selex .csv file
write.csv(catches,paste0(File,"/",assessment,"/catch",assessment,".csv"),row.names = F)
write.csv(cpues,paste0(File,"/",assessment,"/cpue",assessment,".csv"),row.names = F)
write.csv(cpue.se,paste0(File,"/",assessment,"/se",assessment,".csv"),row.names = F)
write.csv(select,paste0(File,"/",assessment,"/select",assessment,".csv"),row.names = F)
write.csv(selex,paste0(File,"/",assessment,"/selex",assessment,".csv"),row.names = F)

# Create rdata object with all information for JABBA-Selecy
ss4js <- list(catch=catches,cpue=cpues,cpue.se=cpue.se,select=select,selex.pars=selex.pars,
                   stock.pars = stock.pars)

save(ss4js,file=paste0(File,"/",assessment,"/ss4js_",assessment,".rdata"))

