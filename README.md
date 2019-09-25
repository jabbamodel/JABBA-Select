---
title: "JABBA Select - Overview"
author: "Henning Winker"
date: "Cape Town, 2019"
output: html_document
---

<br />

## JABBA-SELECT
### Incorporating life history and fisheries’ selectivity into surplus production models
The materials in this repository present the JABBA-Select stock assessment model (Winker et al. 2019). JABBA-Select textends the Bayesian state-space surplus production model JABBA [(Winker et al. 2018)](https://www.sciencedirect.com/science/article/pii/S0165783618300845) to account for selectivity-induced distortion of abundance indices and impacts on stock productivity. Like JABBA, JABBA-Select is implemented in JAGS, called from the statistical programming environment R. JABBA-Select retains the core features of a basic JABBA modelling framework (Winker et al., 2018), including its modular coding structure, a suite of options to fix or estimate process and observation variance components and inbuilt graphics to illustrate model fit diagnostics and stock status results. 

Central to our approach is the integration of prior information from spawning biomass- and yield-per-recruit models with integrated Beverton-Holt spawner recruitment relationship (BH-SRR) into JABBA-Select, which we subsequently refer to as age-structured equilibrium model (ASEM).To directly link the generalized three parameter SPM by Pella and Tomlinson (1969) to the ASEM, we assume that surplus production is a function of spawning biomass and then express surplus production as a function of our formulation of <i>H<sub>MSY</sub></i> instead of the intrinsic rate of population increase, so that:    


<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/Figures/SPM_Eq.1.png" height="70">     (Equation 1)

where <i>SB<sub>0</sub></i> is the unfished biomass and m is a shape parameter that determines at which <i>SB/SB<sub>0</sub></i> ratio maximum surplus production is attained. The functional links between the ASEM and Pella-Tomlinson SPM are illustrated in Fig. 1, which provides a means to translate typical input parameters of age-structured models into the key SPM parameters <i>r</i> and <i>m</i>. 

<br />
<br />

<img src="https://github.com/jabbamodel/JABBA-Select/blob/master/Figures/Fig2_Schematic.png" width="600">

<i> Fig. 1. Schematic of functional relationships between the productivity parameter r and the shape parameter of the surplus production function and the Age-Structured Equilibrium Model (ASEM; i.e. yield- and spawning biomass-per-recruit models with integrated spawner recruitment relationship). Numbers in boxes denote the sequence of deriving deviates of r and m from life history and selectivity parameter inputs into the ASEM. </i>

<br />

JABBA-Select has four novel elements compared to conventional Surplus Production Models:

+ The model uses the expression of harvest rate at MSY (<i>H<sub>MSY</sub></i>), which we define here as <i>H<sub>MSY</sub> = MSY /SB<sub>MSY</sub> </i>, as a surrogate for the intrinsic rate of population increase <i>r</i>, and derives the shape parameter <i>m</i> of the surplus production curve as a function of SBMSY/SB0. This provides a means to generate prior distributions of likely values of <i>H<sub>MSY</sub></i> and <i>m</i> from the ASEM using life history parameters and fishery-selectivity inputs (Fig. 2a)
+ The parameter H<sub>MSY<sub>s</sub></sub> is specific to fishing operations that fish with selectivity <i>s</i> and can be adjusted to account for selectivity-induced changes in the overall year-specific stock productivity <i>H<sub>MSY<sub>y</sub></sub></i> as well as on the abundance indices (Fig. 3b).
+ The model separates between exploitable biomass EBs and spawning biomass SB; the former is used to fit indices given selectivity s, and the latter to predict surplus production. The parameters used to describe the ratio of <i>EB<sub>s,y</sub></i> and <i>SB<sub>y</sub></i>, as a function of spawning biomass depletion relative to average unfished levels are inferred from the ASEM (Fig. 1c)
+ The model accounts for the underlying correlation structure between generated values <i>H<sub>MSY</sub></i> and <i>m</i> through the formulation of a multivariate normal (MVN) prior, which allows for estimating both parameters jointly within the model (Fig. 1d).   

![Figure 2](https://github.com/jabbamodel/JABBA-Select/blob/master/Figures/Fig1_4elements.png)
<i> Fig. 2.  Illustration of the four novel elements of JABBA-Select based on the stock parameters for silver kob: (a) Comparison of the functional forms of the yield curves produced from the Age-Structured Equilibrium Model (ASEM) with the approximation by the JABBA-Select surplus production function (Eq. 1) as function spawning biomass depletion SB / SB<sub>0</sub>, using the life history parameter input values and a range of length-at-50%-selectivity values; (b) JABBA-Select model estimates  of time-varying productivity parameters of H<sub>MSY<sub>y</sub></sub>, (c) ASEM-derived selectivity-dependent distortion in the exploitable biomass (EB) relative to the spawning biomass (SB) over a wide a range of SB<sub>0</sub> iterations, with the dashed line denoting the increase in minimum size limit for line-caught silver kob and the remainder of variations attributed to variations in the relative catch contribution the of inshore trawl; and (d) Multivariate normal (MVN) approximation of log⁡(H<sub>MSY<sub>f,s,k</sub></sub> and log(m<sub>f,s,k</sub>) random deviates generated from the ASEM via Monte-Carlo simulations </i>.

<br />

JABBA-Select is able to accommodate multiple catch time series, changes in selectivity within each fishery (e.g. due to gear
regulations), and can be simultaneously fitted to multiple abundance indices with varying selectivity. JABBA-Select reduces to a conventional Pella-Tomlinson model (JABBA-PTM) if SB = EB and then estimates a single HMSY independent of selectivity. This JABBA-PTM can be evoked via an implemented “user option” that sets all selectivity functions associated with the catch time series and abundance indices equal to the asymptotic maturity curve parameterization.



<b> This Repository includes:

+ The JABBA-Select source code [`JABBA_SELECTv1.1.R`](https://github.com/JABBAmodel/JABBA-SELECT/blob/master/JABBA_SELECTv1.2beta.R), which is called from an assessment specific "Prime" file.

+ A basic [simulation example] (https://github.com/Henning-Winker/JABBA-SELECTbeta/blob/master/KOBsim_example) [`JABBA_SELECT_prime_KOBsim.R`](https://github.com/Henning-Winker/JABBA-SELECTbeta/blob/master/KOBsim_example/JABBA_SELECT_prime_KOBsim.R). New data can be generated from a stochastic age-structured simulation model from  [`asmsim2jabba_data.R`](https://github.com/Henning-Winker/JABBA-SELECTbeta/blob/master/KOBsim_example/asmsim2jabba_data.R), which calls the function [`OM_ccsra2jabba_Fn.R`](https://github.com/Henning-Winker/JABBA-SELECTbeta/blob/master/KOBsim_example/OM_ccsra2jabba_Fn.R) that has been adopted from the R package [CCSRA](https://github.com/James-Thorson/CCSRA).

+ A complex worked example that replicates the Stock Synthesis (ss3) base-case scenario of the 2017 [ICCAT assessment of North Atlantic swordfish](https://www.iccat.int/Documents/Meetings/Docs/2017_ATL_SWO_ASS_REP_ENG.pdf) with JABBA-Select.   

</b>


**Reference**

[Winker, H., Carvalho, F., Kapur, M. (2018) <U>JABBA: Just Another Bayesian Biomass Assessment.</U> *Fisheries Research* **204**: 275-288.](https://www.sciencedirect.com/science/article/pii/S0165783618300845)   

