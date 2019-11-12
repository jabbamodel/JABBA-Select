JABBA-Select Model Execution Vignette
================
Written by Henning Winker <br>
henning.winker@gmail.com

### Getting started

A JABBA-Select model is desinged to be set up in a R 'Prime' file, from which the JABBA-Select source code [`JABBA_SELECTv1.1.R`](https://github.com/jabbamodel/JABBA-Select/blob/master/JABBA_SELECTv1.1.R) via the `source()` function. 

This tutorial explains the main segments of the Prime file setup using the simulation example for South African silver kob [`JABBA_SELECT_prime_KOBsim.R`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/JABBA_SELECT_prime_KOBsim.R) presented in [Winker et al. (2020)](https://www.sciencedirect.com/science/article/pii/S0165783619302103). 

JABBA-Select requires the installation of [R](https://cran.r-project.org/) and [JAGS](https://sourceforge.net/projects/mcmc-jags/) and the following R packages that are installed within R if not available already. 



``` r
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
## Prime File for KOBsim JABBA-SELECT example
## written by Henning Winker - Cape Town
## henning.winker@gmail.com
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

```

### Simulation example of Silver Kob 

The simulation example of Silver Kob is used to illustrate the basic steps of fitting a JABBA-Select model using the Prime file [`JABBA_SELECT_prime_KOBsim.R`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/JABBA_SELECT_prime_KOBsim.R). 

The population dynamics of silver kob are those assumed in Winker et al. (2019). The availbe data are catch and catch-per-unit effort (CPUE) time series over a period of 40 years. A simple logistic selectivity function is assumed, which is described by the lengths where 50% and 95% of fish are retained as catch. Length-at-50%-selectivity increased after year 25 to simulate an increase in miminum size limit regulations.  Here, two CPUE indices are simulated, which are associated with different standard errors catchability q (different scale) to represent, for example, CPUE indices from two areas that are sampled with varying precision (Winker et al. 2019).


### Input files

JABBA-Select requires four to five input comma-separated value files (.csv). 

#####<i>1.  Catch </i>
The [`catch`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) input file contains the time series of year and catch by weight, seperated fisheries that operate with selectivity patterns. Missing catch years or catch values are not allowed. 

##### <i>2. CPUE </i>
JABBA-SELECT is formulated to accommodate abundance indices from multiple sources (i.e., fleets) in a single [`cpue`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) file, which contains all considered abundance indices. The first column of the [`cpue`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) input is year, which must match the range of years provided in the [`catch`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) file. If a CPUE index is affected by a change selectivity, it has be seperately captured in two columns. In contrast to the [`catch`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) input, missing abundance index values are allowed, such that different abundance indices may correspond to smaller portions of the catch time series. Optionally, an additional [`se`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/seKOBsim.csv) input can be passed onto JABBA, containing standard error estimates associated with the abundance indices on a log scale. The [`se`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/seKOBsim.csv) input file is structurally identical to the [`cpue`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) input. Alternatively, this feature can be used to apply different weighting to individual abundance indices by assigning varying coefficients of variation (CV) to each time series. If such weighting is implemented, it is advised that the CV chosen for each indexed year approximates the observed standard error on the log scale, such that the data weights are congruent with expectations as to how well the model should fit these data [(Winker et al. 2018)](https://www.sciencedirect.com/science/article/pii/S0165783618300845).

##### <i> 3. Selectivity functions </i>
JABBA-Select aims to account for the effects of different selectivity patterns on the stock's surplus production and arising distortions between spawing biomass (SB) and exploitable biomass (EB) when fitting the CPUE data. The selectivity functions under consideration are summarized in the [`selex`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selexKOBsim.csv) .csv file. The JABBA selex function provides the option to specify a 2-parameter logistic as well as a 5-parameter piece-wise dome-shaped selectivity curve, with a logistic function for the ascending limb and the descending limb described by the mean and CV of a half-normal distribution (Huynh et al., 2018). 
<br><br>
The KOBSim example only focuses on the simple 2-parameter logistic function, which requires the lengths where 50% and 95% of fish are retained as catch. The early and recent selectivity functions are referenced a S1 and S2, respectively. The first column in the [`selex`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selexKOBsim.csv) file lists the 5 parameters SL50, SL95 (anscending logistic), SL.desc (mean of double-normal), CV.desc (rate of descent) and min.desc. The next columns provide the parameter values. If the latter three parameters (SL.desc, CV.desc, min.desc) are left blank or provided as `NA`, the selex function automically reduces to a logistic.

#####<i>4. Assigning Selectivity to Catch and CPUE </i>

JABBA-Select requires to assign a selection function to each Catch and CPUE time series provided in the form of columns in the [`catch`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) and the
[`cpue`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) input files, respectively. The [`select`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selectKOBsim.csv) input file is designed to facilitate this. <br><br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/Figures/select.png" width="700">

Column `Fleet.ID` specifies an index of unique numbers for each Catch and CPUE time series, where Catch series must be listed first, followed by CPUE indices. Column `Fleet` specifies the name labels for of each Catch and CPUE series corresponding column names chosen by the user in the [`catch`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) and the
[`cpue`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) input files.         


<br><br><br>
The four (or five) csv files are named such that the type classifier `catch`, `cpue` and `se`, `selex` and `select` is combined with the `assessment` name. In this example, we define `assessment = "KOBSim"`, so that (i) [`catchKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv), (ii) [`cpueKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv), (iii) [`seKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/seKOBsim.csv), (iv) [`selexKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selexKOBsim.csv) and (v) [`selectKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selectKOBsim.csv).


All input files have to be saved in a folder that is named after the `assessment`, here `/KOBSim`.

In the Prime file:
`File =` requires the path where the assessment folder is located
`JABBA =` requires the path where the JABBA model `JABBAv1.1.R` is located
`version =` determines the JABBA model version
`assessment =` assignes the assessment folder name

``` r
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

```

### Basic settings

JABBA provides various Graphic, Output, Saving options that can be specified in the prime file. If not specified, JABBA will automatically use the default settings as specified on top of the [`JABBAv1.1.R`](https://github.com/jabbamodel/JABBA/blob/master/JABBAv1.1.R) code.

``` r
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
# Save entire posterior as .RData object
save.all = FALSE # (if TRUE, a very large R object of entire posterior is saved)  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
```

### Looping through scenarios

JABBA makes it easy to run alternative scenarios in a loop. For this purpose, a unique name has to be assgined to each scenario. For example, a generic option to do this for 10 alternative scenarios is:

`Scenarios = c(paste0("Scenario",1:10))`

but individual names may be specified as well, e.g. `Scenarios = c("Run_high_r","Run_medium_r","Run_low_r")`. JABBA automatically creates a folder for each scenario, including the `Input` and `output` subfolders. Note that the type of model `_Schaefer`, `_Fox` and `_Pella` is automatically added to the name of scenario folder.

``` r
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Optional: Note Scenarios
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# S1: Model including Brazil1 
# S2: Model excluding Brazil1
# S3: Base-case Model with time blocks on ESP and JPN 
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Specify Scenario name for output file names
Scenarios = c(paste0("Scenario",1:3)) 

# Execute multiple JABBA runs in loop 
for(s in 1:3){
  Scenario = Scenarios[s] 
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Suplus Production model specifications
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # Choose model type: 
  # 1: Schaefer
  # 2: Fox  
  # 3: Pella-Tomlinsson  
  
  Model = c(3,3,3,3,3)[s] 
  Mod.names = c("Schaefer","Fox","Pella")[Model]
  
  # Depensation opiton:
  # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25) 
  # Choose Plim = 0 to reduce to conventional Schaefer, Fox, Pella models 
  Plim = 0
  
  # Required specification for Pella-Tomlinson (Model = 3)
  BmsyK = 0.4 # Set Surplus Production curve inflection point
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>  
```

### Read input data files

Here, the user can specify if the `se` file is available by setting \`SE.I = TRUE' We advice to always check if all csv files were correctly read-in.

``` r
#--------------------------------------------------
# Read csv files
#--------------------------------------------------
  
  # Use SEs from csv file for abudance indices (TRUE/FALSE)
  SE.I = TRUE
  
  # Load assessment data
  catch = read.csv(paste0(File,"/",assessment,"/catch",assessment,".csv"))
  cpue = read.csv(paste0(File,"/",assessment,"/cpue",assessment,".csv"))#
  
  if(SE.I ==TRUE){
    se =  read.csv(paste0(File,"/",assessment,"/se",assessment,".csv"))
  }
  
  head(cpue)
  head(catch)
```

### Input data manipulation

This section is optional. For the SWO\_SA example, we use this section to manipulate the input `cpue` and `se` data for scenarios 1-3. For sceanarios 2-3, we remove the BRA\_LL1 series. For the base-case scenario 3, we also use seperate time series (time-blocks) for Japan (JP\_LL1 & JP\_LL2) and Spain (SP\_LL1 & SP\_LL2), which we merge back into one index for JP\_1 and on for SP\_LL1 for scenarios 1-2. If the `se` input is provided, the same munipulations must be applied as for `cpue`.Again, always check if the manilupations were correctly applied.

JABBA provides provides the option to use a single averaged CPUE index instead of the individual abundance indices (see *2.5.1. State-Space model for averaging of abundance indices* in Winker et al, 2018). This feature can be activated by setting `meanCPUE = TRUE`.

``` r
  #--------------------------------------------------
  # option to exclude CPUE time series or catch year
  #--------------------------------------------------
  
  if(s<3){ # Combine SPAIN and JAPAN
    cpue[,4] = apply(cpue[,4:5],1,mean,na.rm=TRUE)
    cpue[,6] = apply(cpue[,6:7],1,mean,na.rm=TRUE)
    
    cpue = cpue[,-c(5,7)] 
    se[,4] = apply(se[,4:5],1,mean,na.rm=TRUE)
    se[,6] = apply(se[,6:7],1,mean,na.rm=TRUE)
    
    se = se[,-c(5,7)]
    cpue[!is.finite(cpue[,4]),4]=NA
    cpue[!is.finite(cpue[,5]),5]=NA
    se[!is.finite(se[,4]),4]=NA
    se[!is.finite(se[,5]),5]=NA
    
  }
  
  # Remove BrazilI
  if(s>1){
    cpue = cpue[,-c(2)]
    se = se[,-c(2)]
  }

  
  names(cpue)
  ncol(catch)
  ncol(cpue)

  #------------------------------------------------------
  # Option use mean CPUE from state-space cpue averaging
  #-----------------------------------------------------
  meanCPUE = FALSE
```

### Prior and Process variance settings

Please see section *2.3.2. Prior specification* in Winker et al. (2008) for details.

Most prior settings provide more than one option. For example, if the prior for K is meant to be specified as a lognormal prior set `K.dist = c("lnorm","range")[1]`, whereas for a range set `K.dist = c("lnorm","range")[2]`. If the prior for K is specified as lognormal, e.g. `K.prior = c(200000,1)`, it requires the untransformed mean K and the assumed CV. If the prior for K is specified as range, it requires the assumed minum and maximum values, e.g. `K.prior = c(15000,1500000)`.

The r prior provides an additional option, in that it can be specified as a generic resiliance category *Very low, Low, Medium* or *High*, such as provided by [FishBase](www.FishBase.org). This requires specifying `K.dist = c("lnorm","range")[2]` (i.e. as a range) and then setting the `K.prior` equal to one of the above reliance categories, e.g. `K.prior = "Low"`.

``` r
  #------------------------------------------------
  # Prior for unfished biomass K
  #------------------------------------------------
  # The option are: 
  # a) Specify as a lognormal prior with mean and CV 
  # b) Specify as range to be converted into lognormal prior
  
  K.dist = c("lnorm","range")[1]
  
  # if lnorm use mean and CV; if range use lower,upper bound
  K.prior = c(200000,1) 
  
  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  # Set the initial depletion prior B1/K 
  # To be converted into a lognormal prior (with upper bound at 1.1)
  
  psi.dist= c("lnorm","beta")[1]
  # specify as mean and CV 
  psi.prior = c(1,0.25) 
  
  #----------------------------------------------------
  # Determine r prior
  #----------------------------------------------------
  # The option are: 
  # a) Specifying a lognormal prior 
  # b) Specifying a resiliance category after Froese et al. (2017; CMSY)
  # Resilience: "Very low", "Low", "Medium", High" (requires r.range = TRUE)
  
  # use [1] lognormal(mean,stdev) or [2] range (min,max) or
  r.dist = c("lnorm","range")[1] 
  
  r.prior = c(0.42,0.37) 
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Process Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  #Estimate set sigma.proc == True
  sigma.proc = TRUE
  # Determines if process error deviation are estimated for all years (TRUE)  
  # or only from the point the first abundance index becomes available (FALSE)
  proc.dev.all = FALSE 
  #------------------------------------------
  if(sigma.proc == TRUE){
    igamma = c(4,0.01) #specify inv-gamma parameters
    
    # Process error check
    gamma.check = 1/rgamma(1000,igamma[1],igamma[2])
    # check mean process error + CV
    mu.proc = sqrt(mean(gamma.check)); CV.proc =    sd(sqrt(gamma.check))/mean(sqrt(gamma.check))
    
    # check CV
    round(c(mu.proc,CV.proc),3)
    quantile(sqrt(gamma.check),c(0.1,0.9))
  }else{
    sigma.proc = 0.07 #IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
  }
```

Both catchability *q* and the estimable observation variance *σ*<sub>*e**s**t*, *i*</sub><sup>2</sup> can be specified to be estimated: (1) for each CPUE index, (2) in groups or (3) as the same quantatity for all indices. For (1), simply provide a vector of unique integer in order for each index, e.g. `sets.q = 1:(ncol(cpue)-1)`. For (2), `set.q =` can be specified by grouping similar indices, e.g. `set.q = c(1,1,2,2,3)`. For (3), simply provide the indentifier 1 for all indices, e.g. `sets.q = rep(1,ncol(cpue)-1)`. The exact same principles apply for assigning *σ*<sub>*e**s**t*, *i*</sub><sup>2</sup> to individual indices *i*, i.e. `sets.var = 1:(ncol(cpue)-1)` for case (1).

``` r
  #--------------------------------------------------------------
  # Determine estimation for catchability q and observation error 
  #--------------------------------------------------------------
  # Assign q to CPUE
  sets.q = 1:(ncol(cpue)-1) 
  
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Observation Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  
  #To Estimate additional observation variance set sigma.add = TRUE
  sigma.est = TRUE
  
  # Series
  sets.var = 1:(ncol(cpue)-1) # estimate individual additional variace
  
  # As option for data-weighing
  # minimum fixed observation error for each variance set (optional choose 1 value for both)
  fixed.obsE = c(0.2) # Important if SE.I is not availble
  
  # Total observation error: TOE = sqrt(SE^2+sigma.est^2+fixed.obsE^2)
  #--------------------------------------------
```

### Projections under constant Total Allowable Catch (TAC)

JABBA enables projections under constant catch scenarios. JABBA automatically compares the difference between the last assessment year and the present year. The difference between these years is projected forward under the *current* catch, which could, for example, be determined based on updated catch inofrmation `TACint = 10058` or by assuming an average catch based on the three most recent assessment years, such that `TACint = mean(catch[nrow(catch)-3,2]:catch[nrow(catch),2])`. All prjected posteriors can be saved a `_projections.Rdata` object, which can be easily passed on JABBAgoesFLR.R for further processing, including the production of Kobe projection matrices.

``` r
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Optional: Do TAC Projections
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  Projection = TRUE # Switch on by Projection = TRUE 
  
  # Check final year catch 
  catch[nrow(catch),]
  
  # Set range for alternative TAC projections
  TACs = seq(10000,18000,1000) #example
  
  # Intermitted TAC to get to current year
  #TACint = mean(catch[nrow(catch)-3,2]:catch[nrow(catch),2]) # avg last 3 years
  TACint = 10058 # Catch for 2016
  # Set number of projections years
  pyrs = 10
```

### MCMC setting and Model execution

``` r
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Execute model and produce output
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  
  # MCMC settings
  ni <- 30000 # Number of iterations
  nt <- 5 # Steps saved
  nb <- 5000 # Burn-in
  nc <- 2 # number of chains
  nsaved = (ni-nb)/nt*nc # MUST be an integer
  
  # Run model (JABBA model file, must be in the same working directory)
  source(paste0(JABBA.file,"/JABBA",version,".R")) 
  
  }# THE END
```
