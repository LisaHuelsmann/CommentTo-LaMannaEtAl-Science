# Data simulation and CNDD analyses
Lisa Hülsmann & Florian Hartig  
February 15 2018  



## Goal of the analysis

The goal of the simulation null model was to be able to create synthetic data with and without CNDD, varying also a range of other processes and community characteristics that could affect the analysis (see list below), to evaluate the behaviour of the Ricker and the offset-power model under different conditions. 

CNDD and the correlation of CNDD with abundance were assessed for

*  synthetic data created from a simulation model with the settings explained before, and
*  data from the Barro Colorado Island (BCI, 7th survey) plot ([*1*](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/tree/master/code#1-r-condit-s-lao-r-p%C3%A9rez-s-b-dolins-r-b-foster-s-p-hubbell-barro-colorado-forest-census-plot-data-2012-version-center-for-tropical-forest-science-databases-2012)).

*Note that we cannot provide the original dataset in this repository due to the data sharing agreement of the BCI data.* 

It must be requested from the data owners via http://ctfs.si.edu/webatlas/datasets/bci/


## Simulation settings
We used randomly simulated data with varying process and community characteristics. Counts of recruits and adults were generated for 10x10m quadrats assuming different CNDD, dispersal, adult survival, habitat specificity, species abundance and adult proportion. Counts of adults and recruits in each quadrat were randomly drawn from a Poisson distribution with mean values defined by species abundance, the proportion of adults and the dynamic processes as follows.

#### CNDD
Conspecific negative density dependence was implemented by reducing local recruits as a function of local adults following an exponential survival function (see [code](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/blob/81fcea17e852a8792407676473f762d679235ac5/code/functions_data_simulation.R#L65-L67)). To maintain the recruit-to-adult-ratio, we randomly placed the removed recruits to all quadrats (see [code](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/blob/81fcea17e852a8792407676473f762d679235ac5/code/functions_data_simulation.R#L70-L71)). We vary CNDD between 0 (i.e. no CNDD) and -2 (very strong CNDD) using the same value for all species.

#### Dispersal
Dispersal was defined as the percentage of local recruits that originate from conspecific adults in the meta-community; the remaining recruits originate from conspecific adults in the same quadrat (see [code](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/blob/1a0e902c91dd0cc21f1fcbf144f8d6b686842cc0/code/functions_data_simulation.R#L24-L47)). Dispersal was varied between 0 (i.e. local dispersal) and 1 (i.e. global dispersal). LaManna *et al*. assumed that only 10% of recruits derive from outside the quadrats which seems to produce unrealistically low dispersal considering empirical values for primary dispersal ([*2*](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/tree/master/code#2-h-c-muller-landau-s-j-wright-o-calder-r-condit-s-p-hubbell-interspecific-variation-in-primary-seed-dispersal-in-a-tropical-forest-j-ecol-96-653-667-2008)). A uniform dispersal kernel with only 10m radius would predict approximately 52% of all seeds to fall outside a 10x10m quadrat.

#### Adult survival
To account for the time that a seed needs to develop a DBH above the threshold of 1cm, we randomly removed 1-survival of the adults after generating the recruits (see [code](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/blob/81fcea17e852a8792407676473f762d679235ac5/code/functions_data_simulation.R#L58-L61)). To obtain the same density of adults after mortality, we increased the initial abundance of adults correspondingly (see [code](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/blob/81fcea17e852a8792407676473f762d679235ac5/code/functions_data_simulation.R#L18)). We vary adult survival between 0.1 (i.e. observed recruits and adults are strongly decoupled, which corresponds to a species with very fast dynamics) and 1 (i.e. no temporal decoupling, which corresponds to the instantaneous dispersal of recruits).

#### Habitat specificity
Niche effects were considered by allowing a species to occupy only a proportion of the quadrats (see [code](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/blob/81fcea17e852a8792407676473f762d679235ac5/code/functions_data_simulation.R#L10)). With increasing specificity, the local density of adults and recruits at suitable quadrats increases, although species abundance at the plot scale remains equal. We varied habitat specificity between 0 (i.e. no niche effects) and 0.9 (i.e. strong niche effects).

#### Species richness and abundances
We generated species abundances for plots with increasing species richness using an exponential decay function (see [code](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/blob/81fcea17e852a8792407676473f762d679235ac5/code/functions_data_simulation.R#L112-L116)). We varied species richness between 13 and 447 species to reproduce a latitudinal gradient between temperate forests and the tropics.

#### Adult proportion
The proportion of individuals that are considered as adults was attempted at 0.2 by LaManna *et al*., but realized proportions varied due to a stepwise reduction of the DBH limit for recruit identification (i.e. 10, 5 or 2cm). On average, 0.4 of the individuals were adults in the BCI data. We varied adult proportion between 0.05 and 0.75 (see [code](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/blob/81fcea17e852a8792407676473f762d679235ac5/code/functions_data_simulation.R#L20-L21)).

## CNDD Analyses
Parameter estimation followed the methodology by LaManna *et al*., i.e. using the Ricker and the offset-power models, including weighted regressions for the CNDD-abundance correlation, and models accounting for heterospecific adult and recruit density. 

Note that we show abundance as the number of trees per hectare (N), since the CNDD bias is mathematically linked to N, but N and basal area were strongly correlated in LaManna *et al*..

## References

###### 1: R. Condit, S. Lao, R. Pérez, S. B. Dolins, R. B. Foster, S. P. Hubbell, Barro Colorado Forest Census Plot Data, 2012 Version. Center for Tropical Forest Science Databases (2012).
###### 2: H. C. Muller-Landau, S. J. Wright, O. Calder, R. Condit, S. P. Hubbell, Interspecific Variation in Primary Seed Dispersal in a Tropical Forest. *J Ecol* **96**, 653-667 (2008).

