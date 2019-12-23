JABBA-Select Basic User Guide
==============================

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

The population dynamics of silver kob are those assumed in Winker et al. (2019). The availbe data are catch and catch-per-unit effort (CPUE) time series over a period of 40 years. A simple logistic selectivity function is assumed, which is described by the lengths where 50% and 95% of fish are retained as catch. Length-at-50%-selectivity increased after year 25 to simulate an increase in miminum size limit regulations.  Here, two CPUE indices are simulated, which are associated with different standard errors catchability q (different scale) to represent, for example, CPUE indices from two areas that are sampled with varying precision ([Winker et al. 2020)](https://www.sciencedirect.com/science/article/pii/S0165783619302103).


### Input files

JABBA-Select requires four to five input comma-separated value files (.csv). <br>
The four (or five) csv files are named such that the type classifier `catch`, `cpue` and `se`, `selex` and `select` is combined with the `assessment` name. In this example, we define `assessment = "KOBSim"`, so that (i) [`catchKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv), (ii) [`cpueKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv), (iii) [`seKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/seKOBsim.csv), (iv) [`selexKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selexKOBsim.csv) and (v) [`selectKOBSim.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selectKOBsim.csv). All file must be located in the assessment folder, i.e. `/KOBSim`.
 


##### <i>1.  Catch </i>
The [`catch.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) input file contains the time series of year and catch by weight, seperated fisheries that operate with selectivity patterns. Missing catch years or catch values are not allowed. In this simple example of [`JABBA_SELECT_prime_KOBsim.R`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/JABBA_SELECT_prime_KOBsim.R), fishing selectivity changed for the entire fleet between year 25 and 26.
<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/SELECT_JS/Input/Catches_KOBSim.png" width="900" >
<br>

In the more complex Atlantic Swordfish example [SWOss3](https://github.com/jabbamodel/JABBA-Select/tree/master/SWOss3) the catches seperate into eleven fleets, fishing with different selectivity pattern (see below). 
<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/SWOss3/SWOselect_JS/Input/Catches_SWOss3.png" width="900" >
<br>



##### <i>2. CPUE </i>
JABBA-SELECT is formulated to accommodate abundance indices from multiple sources (i.e., fleets) in a single [`cpue.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) file, which contains all considered abundance indices. The first column of the [`cpue`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) input is year, which must match the range of years provided in the [`catch.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) file. If a CPUE index is affected by a change selectivity, it has be seperately captured in two columns. In contrast to the [`catch.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) input, missing abundance index values are allowed, such that different abundance indices may correspond to smaller portions of the catch time series. Optionally, an additional [`se.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/seKOBsim.csv) input can be passed onto JABBA, containing standard error estimates associated with the abundance indices on a log scale. The [`se.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/seKOBsim.csv) input file is structurally identical to the [`cpue.csv`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) input. Alternatively, this feature can be used to apply different weighting to individual abundance indices by assigning varying coefficients of variation (CV) to each time series. If such weighting is implemented, it is advised that the CV chosen for each indexed year approximates the observed standard error on the log scale, such that the data weights are congruent with expectations as to how well the model should fit these data [(Winker et al. 2018)](https://www.sciencedirect.com/science/article/pii/S0165783618300845).

<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/SELECT_JS/Input/CPUE_KOBSim_SELECT.png" width="900" >
<br>


##### <i> 3. Selectivity functions </i>
JABBA-Select aims to account for the effects of different selectivity patterns on the stock's surplus production and arising distortions between spawing biomass (SB) and exploitable biomass (EB) when fitting the CPUE data. The selectivity functions under consideration are summarized in the [`selex`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selexKOBsim.csv) .csv file. The JABBA selex function provides the option to specify a 2-parameter logistic as well as a 5-parameter piece-wise dome-shaped selectivity curve, with a logistic function for the ascending limb and the descending limb described by the mean and CV of a half-normal distribution (Huynh et al., 2018). 
<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/SWOss3/SELEX_SWOss3.png">
<br>
<br>
The KOBSim example only considers on the simple 2-parameter logistic function, which requires the lengths where 50% and 95% of fish are retained as catch. The early and recent selectivity functions are referenced a S1 and S2, respectively. The first column in the [`selex`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selexKOBsim.csv) file lists the 5 parameters SL50, SL95 (anscending logistic), SL.desc (mean of double-normal), CV.desc (rate of descent) and min.desc. The next columns provide the parameter values. If the latter three parameters (SL.desc, CV.desc, min.desc) are left blank or provided as `NA`, the selex function automically reduces to a logistic. 
<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/Figures/selex.png" width="900" height="200">


##### <i> 4. Assigning Selectivity to Catch and CPUE </i>

JABBA-Select requires to assign a selection function to each Catch and CPUE time series provided in the form of columns in the [`catch`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) and the
[`cpue`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) input files, respectively. The [`select`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selectKOBsim.csv) input file is designed to facilitate this. <br><br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/Figures/select.png" width="900" height="200">

Column `Fleet.ID` specifies an index of unique numbers for each Catch and CPUE time series, where Catch series must be listed first, followed by CPUE indices. Column `Fleet` specifies the name labels for of each Catch and CPUE series corresponding column names chosen by the user in the [`catch`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/catchKOBsim.csv) and the
[`cpue`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/cpueKOBsim.csv) input files. The column `Selectivity` assigns the respective `Fleet` to a selectivity function provided in the [`selex`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selexKOBsim.csv) file. In this example, the catches of fleet S1 and S2 are assigned to the first and second listed selectivity functions (S1 & S2) in in the [`selex`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/selexKOBsim.csv) by inputting 1 and 2, respectively. The CPUE indices (A1.S1, A1.S2, A2.S1, A2.S2) are then assigned in the same way. Column `CPUE.units` specify if CPUE is in weight = 1 or number = 0. Note that catch must be in weight. Column `Catch` specifies whether the time series represents catch, if TRUE, or cpue, if FALSE. Finally, column `q` assigns a estimable catchability coefficients each cpue index or 0 to catch. 


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

All input files have to be saved in a folder that is named after the `assessment`, here `/KOBSim`.

In the Prime file:
`File =` requires the path where the assessment folder is located
`JABBA =` requires the path where the JABBA model `JABBA_SELECTv1.1.R` is located
`version =` determines the `JABBA_SELECT` version
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

JABBA-Select provides various Graphic, Output, Saving options that can be specified in the prime file. If not specified, JABBA will automatically use the default settings as specified on top of the [`JABBA-SELECTv1.1.R`](https://github.com/jabbamodel/JABBA-Select/blob/master/JABBA-SELECTv1.1.R) code.

``` r
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
```

### Looping through scenarios

JABBA-Select makes it easy to run alternative scenarios in a loop. For this purpose, a unique name has to be assgined to each scenario. For example, a generic option to do this for 10 alternative scenarios is:

`Scenarios = c(paste0("Scenario",1:10))`

but individual names may be specified as well, e.g. `Scenarios = c("Run_high_h","Run_medium_h","Run_low_h")`. In this example the Prime file is set up to run the full SELECT model, followed by a reduced Pella-Tomlinson model, which does not account for selectivity by assuming that exploitable biomass equals spawning biomass (SB = EB). 
``` r
Scenarios = c("SELECT","Pella")
s=1

for(s in 1:length(Scenarios)){
  Scenario = Scenarios[s] 
  if(s==2) SELECT = FALSE # Reduce JABBA-Select to a Bayesian Pella-Tomlinson
```
JABBA-Select automatically creates a folder for each scenario, including the [`Input`](https://github.com/jabbamodel/JABBA-Select/tree/master/KOBsim_example/KOBsim/SELECT_JS/Input) and [`Output`](https://github.com/jabbamodel/JABBA-Select/tree/master/KOBsim_example/KOBsim/SELECT_JS/Output) subfolders, where all default plots, .csv and .rdata files are saved automatically depending on the choices of Basic Settings (see above). 
<br>

### Read input data files

Here, the user can specify if the `se` file is available by setting \`SE.I = TRUE' We advice to always check if all csv files were correctly read-in.

``` r

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
  
  # Read selex (selectivity) csv
  selex = read.csv(paste0(File,"/",assessment,"/selex",assessment,".csv"))
  
```

### Biomass priors 

JABBA-Select requires two priors on Spawning Biomass (*SB*). The first is a lognormal for the unfished *SB*<sub>0</sub>, which should be typically speficified as vaguely as possible around some plauseble "guess". To this we suggest a CV of 100-200%. To minimize the influence of the choice and improve modelperformance, it adviced that the Posterior median should be always located to right of the highest density peak of the *SB*<sub>0</sub> prior distribution. A very small Prior-Posterior-Variance-Ratio (PPVR) is indicative for the estimate being minimally informaned by the choice of prior.

<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/SELECT_JS/Output/Posteriors_KOBSim_SELECT.png" width="900">
<br>

The second is a prior to inform the the model about initial spawning biomass level at the start of the catch time series `psi` = *SB*<sub>y=1</sub>/*SB*<sub>0</sub>, which can be implemented as either lognormal `psi.prior= c("lognormal","beta")[1]` or beta `psi.prior= c("lognormal","beta")[2]` distribution by specifying the mean and CV of the respective distributions. Here we assumed that the inital *SB* at the start of the catch time series was varying around the unfished *SB*<sub>0</sub> by assuming a lognormal with a mean *SB*<sub>y=1</sub>/*SB*<sub>0</sub> = 1 with a CV = 30% by setting `mu.psi = 1; CV.psi = 0.3`.  

JABBA-Select provides imposes “soft” boundary penalty on *SB* if P<sub>y</sub> = *SB*<sub>y</sub>/*SB*<sub>0</sub>, decreases below or is either equal to or greater than user-defined minimum and maximum values, respectively (default: min=0.02, max=1.3). The idea is that the likelihood is increasingly penalized the further that P<sub>y</sub> diverges from the soft boundaries, thereby improving mixing behaviour of MCMC chains. These constraints can be manually relaxed by setting e.g. `P_bound = c(0.0001,1.5)` if the model runs sufficiently stable.



``` r

 #------------------------------------------------
  # mean and CV and sd for unfished biomass K (SB0)
  #------------------------------------------------
  mu.SB0 = 30000; CV.SB0 = 2; sd.SB0=sqrt(log(CV.SB0^2+1)) 
  SB0.pr = c(mu.SB0,sd.SB0)
    
  
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
    P_bound = c(0.0001,1.5) 

```
### Life History Parameters 

JABBA-Select requires basic life history parameters as input the [Age-Structured Equilibbrium Model (ASEM) function](link here), describing growth, weight-lenght, maturation, longivity, natural mortality and the spawning-recruitment relationship recruitment 

``` r
    #---------------------------------------------------------------
    # STOCK PARAMETERS for prior generation Hmsy as a fuction of r
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
    CV.h = 0.1
    
```
Growth is decribed by the von Bertalnffy Growth Function (VBGF). The maximum age can be treated a plus group `plusgroup = c(FALSE,TRUE)[2]` or as *de facto* maximum age  `plusgroup = c(FALSE,TRUE)[1]`. Maturity as logistic a logistic function that can be either specified  by lengths where 50% (Lm50) and 95% (Lm95%) maturity is attained `maturity = c(Lm50,Lm95,0)` or the ages where 50% (Am50) and 95% (Am95%) maturity is attained `maturity = c(Am50,Am95,1)`, which is selected by choosing 0 for length, and 1 for age as the third value when specifying `maturity`. Uncertainty can be addmitted by specifying CVs for two key paramters natural mortality *M* (assumed to be age-independent) and the steepness parameter *h* of the Beverton and Holt Spawner-Recruitment Relationship.

<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/SELECT_JS/Input/Prior_hM_KOBSim.png" width="900">
<br>

The resulting relationships of length-at-age, weight-at-length, weight-at-age and maturity-at-length relative selectivity-at-length are shown in the `StockFunctions_` plot that is saved in the [`Input`](https://github.com/jabbamodel/JABBA-Select/tree/master/KOBsim_example/KOBsim/SELECT_JS/Input) folder.

<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/SELECT_JS/Input/StockFunctions_KOBSim.png" width="900">
<br>

In addition, JABBA-Select permits to specify growth and length-weight functions to be sex specific by adding `nsexes = 2` (`nsexes = 1` is default). This then requires spefifyinh two values for each paramter: the first for females and the second for males (see our more complext worked example for North Atlantic swordfish [SWOss3](https://github.com/jabbamodel/JABBA-Select/tree/master/SWOss3)). 

``` r
#---------------------------------------------------------------
# STOCK PARAMETERS for prior generation Hmsy as a fuction of r
#---------------------------------------------------------------
minage <- 0  																						
maxage <- 25
PlusGroup = TRUE
# Number of sexes order Female and Males (this model is sex-structured)
nsexes = 2
# VBGF parameters
Linf <- c(290.1,214.2)
kappa <- c(0.147,0.266)
t0 <- c(-1.163, -0.6)
```

If `nsexes = 2`, the maturity function and spawning biomass (SB) is speficific to females.

<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/SWOss3/SWOselect_JS/Input/StockFunctions_SWOss3.png" width="900">
<br>

Based on the Life History parameters and given uncertainty about *M* and *h*, JABBA-SELECT generates a Multivariate Normal (MVN) Prior for *H*<sub>MSY</sub> and the shape *m* of the surplus production curve for selectivity function *s*. 


<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/SELECT_JS/Input/Cor_m_Hmsy_KOBSim.png" width="900">

<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/SELECT_JS/Input/ProdutionKOBSim.png" width="900">
<br>

### Catchability, Observation Error and Process Error settings

Like JABBA, JABBA-SELECT allows the separation of the observation variance into three components: (1) the squared externally estimable observation error (`SE`), that is read-in via [`se`](https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/seKOBsim.csv) csv file, (2) a fixed additional input variance denoted in the R code as `fixed.obsE`, and (3) estimable variance, which is envoked if `sigma.est = TRUE`, where the default prior option for assumes an uninformative inverse-gamma distribution with both gamma scaling parameters set to 0.001. This variance can be estimated individually for each abudance index <i> i </i> for goups of indices or as single quantity common to all indices. All three variance components are additive and can be switched on or off in any combination to provide flexible data-weighting options. Please also see section  *2.3.2. Prior specification* in [Winker et al. (2018)](https://www.sciencedirect.com/science/article/pii/S0165783618300845) for additional details.

The estimable observation variance *σ*<sub>*e**s**t*, *i*</sub><sup>2</sup> can be specified to be estimated: (1) for each CPUE index, (2) in groups or (3) as the same quantatity for all indices. For (1), simply provide a vector of unique integer in order for each index, e.g. `sets.var = select$q[select$CPUE]`. For (2), `set.var =` can be specified by grouping similar indices, e.g. `sets.var = c(1,1,2,2,3)`. For (3), simply provide the indentifier 1 for all indices, e.g. `sets.var = rep(1,ncol(cpue)-1)`. The exact same principles apply for changing the assigned *q*<sub> *i*</sub> for index *i*. For example for option (1), one can simply specify `sets.var = 1:length(select$q[select$CPUE])`. <br>

In this example index A1.S1 and A1.S2 are for example assumed to have a common q = 1, but different selectivity function due to the change in size limits, whereas A1.S1 and A2.S1, which were indipendtly generated from e.g. two regions, are assigned common selectivity functions but different *q*'s. 	 

``` r

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

```
The process variation on log(SB) can be chosen to be estimable `sigma.proc = TRUE` or fixed `sigma.proc = FALSE`. If `sigma.proc = TRUE`, the user has the option to specify a prior for the estimable process variance as inverse-gamma distribution `proc.type = c("igamma","lnorm")[1]` (see [Millar & Meyer 2000](https://www.jstor.org/stable/2680768?seq=1)) or a process error prior (not variance) as lognormal distribution `proc.type = c("igamma","lnorm")[1]`, which can be for example derived of simulation studies ([Winker 2018](http://webcms.uct.ac.za/sites/default/files/image_tool/images/302/pub/2018/IWS2018/Line_Fish/MARAM_IWS2018_Linefish_P3%20-%20SPM_ProcessErrors.docx)). The defaul setting for inverse-gamma is an uninformative approximation of the Jeffrey's prior with `pr.proc = c(0.001,0.001)`. It is probably more intuitive forumalte informative lognormal priors for the process error `pr.proc = c(log(mu),CV)`, where range between 0.05 (long generation times) and 0.15 (short generation times) are thought covering to capture the plausibe ranges (see e.g. [Winker 2018](http://webcms.uct.ac.za/sites/default/files/image_tool/images/302/pub/2018/IWS2018/Line_Fish/MARAM_IWS2018_Linefish_P3%20-%20SPM_ProcessErrors.docx)). A concservative starting value for an informative inverse-gamma prior would `pr.proc = c(4,0.01)` (c.f. [Millar & Meyer 2000](https://www.jstor.org/stable/2680768?seq=1); [(Winker et al. 2018)](https://www.sciencedirect.com/science/article/pii/S0165783618300845)). The provided R code in this example Prime file includes some simple checks for approximating the mean and CV given the choice of scale and shape parameters `pr.proc = c(shape,scale)` for the inverse-gamma distribution. 

If the choice is to fix the process error by setting `sigma.proc = FALSE`, the expected mean specified by e.g. `sigma.proc = 0.1`. Note that short and/or noisy CPUE time series may not permit to estimate the process error reliably due to lack of contrast and information in the data. In such cases, it is recommended to initially develop the JABBA-Select with fixed process error and possibly explore option to use prior once the model runs stable.

``` r
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

```
### Biological Reference Points

By default the stock the stock status reference points are the harvest rate that is required to attain MSY, *H*<sub>MSY</sub> and the corresponding spawning biomass at MSY, *SB*<sub>MSY</sub> by using default `SBmsy_SB0 = NULL`.  


In addition JABBA-Select provides the user option to specify a target *SB*/*SB*<sub>0</sub>. In this example, we adopted *SB*<sub>40</sub>40 = 0.4 × *SB*<sub>0</sub> and the corresponding *H*<sub>40</sub> as the target reference for *SB*, by setting `SBmsy_SB0 = 0.4`. 

<br>
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/SELECT_JS/Output/TrendMSY_KOBSim_SELECT.png" width="900">
<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/KOBsim_example/KOBsim/SELECT_JS/Output/Kobe_KOBSim_SELECT.png" width="900">
<br>

### Projections under constant Total Allowable Catch (TAC)

JABBA-Select enables projections under constant catch scenarios, where weighted average of the annual sustainable
harvest rate H<sub>MSY,y</sub> (see Eq. 11 in [Winker et al. (2020)](https://www.sciencedirect.com/science/article/pii/S0165783619302103)) is determined by relative catch by fleet in the terminal year to account for selectivity induced changes of the stock’s surplus production during the projection phase. Note that projecting relative catch by fleet can be relatively easily customized by adding an additional year with values of to be projected catches by fleet to the `catch` input csv and another year with otherwise empty cells to the  `cpue ` and `se` input files.   

JABBA-SELECT automatically compares the difference between the last assessment year and the quota implemtation year `imp.yr`. The difference between these years is projected forward under the *current* catch. Alternative fixed catches can specied as `TACs`. In this, example, we project over a sequence from 50% to 120% of the current catch `curC` in steps of 10% increaments. The `TACs`are implemented from from year 41 and the projection horizon is set to 10 years as `pyrs = 10`. The projected posteriors can be saved as `_projections.Rdata` object by setting `save.projections = TRUE` under the *Basic settings* section.

``` r
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
    
```

### MCMC setting and Model execution

``` r
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
 

```
