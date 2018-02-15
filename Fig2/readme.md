# Fig. 2
Lisa Hülsmann & Florian Hartig  
February 15 2018  


In Fig. 2, we show how CNDD estimated with the Ricker and the offset-power model varies with abundance for real data from the tropical BCI plot and for simulated data with varying spatial association between adults and recruits. In the second row of the plot, we depict how bias arises in the Ricker and the offset-power model. Details on the simulations can be found [here](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/tree/master/code). 

*Note that we cannot provide the original dataset in this repository due to the data sharing agreement of the BCI data.* 

It must be requested from the data owners via http://ctfs.si.edu/webatlas/datasets/bci/



## Load packages, functions and BCI data


```r
library(sads)
library(gnm)
library(MASS)
library(lme4)
library(mgcv)

source("../code/functions_data_simulation.R")
source("../code/functions_analyses.R")

try(load(file = "../data/bciCounts.Rdata"))

set.seed(123)
```



## Run simulations and analyses

Run analyses with simulated data using four different dispersal settings

```r
fig2 = runAnalyses(abundanceDist = bci*12, numQuadrats= 5000, dispersal = c(0, 0.4, 0.8, Inf))
```

```
## [1] "running 0"
## [1] "running 0.4"
## [1] "running 0.8"
## [1] "running Inf"
```

Run analyses with real counts from the BCI data

```r
if(exists("bciCounts")){
  fig2[[5]] = runAnalysis(externalData = bciCounts)
  names(fig2)[5] = "BCI"
}
```

Run analyses with simulated data using the specifications by LaManna *et al*.

```r
fig2[[6]] = runAnalysis(abundanceDist=bci*12, numQuadrats= 5000, dispersal = 0.1, theta = 1, adultSurvival = "LaManna")
names(fig2)[6] = "LaManna"
```


## Fig. 2


```r
gr <- seq(0.8, 0, length.out = 4)
col = c(rgb(gr, gr, gr, 0.25),  rgb(0, 160, 176, alpha = 0.25*255, max = 255), rgb(235, 104, 65, alpha = 0.25*255, max = 255))
col_full = c(rgb(gr, gr, gr, 0.8),  rgb(0, 160, 176, max = 255), rgb(235, 104, 65, max = 255))
pch = c(rep(16,4), 17,15)


oldpar <- par(mfrow = c(2,2), oma = c(3,3,1,1), mar = c(4,4,2,2), las=1, mgp=c(2.5, 0.5, 0), lwd=75/75, cex.axis=3/3, cex.lab=1, font.lab=2, tcl=(-0.2))

for(j in 1:6) {
  if(!is.null(fig2[[j]])){
    x <- log10(fig2[[j]]$abundance)[fig2[[j]]$RickerLaManna$CNDD > -10]
    y <- fig2[[j]]$RickerLaManna$CNDD[fig2[[j]]$RickerLaManna$CNDD > -10]
    if (j == 1) {
      plot(x, y, main = "Ricker model", xlim = c(-0.4, 3.4),  ylim = c(-8,2), col = col[1], pch=pch[1], xlab = "Species abundance log10(N/ha)", ylab = "CNDD")
    } else {
      points(x, y, col = col[j], cex=1, pch=pch[j])    
    }
    smo <- gam(y ~ s(x, k = 3))
    smo_fit <- predict(smo, newdata = data.frame(x=sort(unique(x))), type = "response")
    lines(sort(unique(x)), smo_fit, col=col_full[j], lwd=2)
  }
}


for(j in 1:6) {
   if(!is.null(fig2[[j]])){
    x <- log10(fig2[[j]]$abundance)
    y <- fig2[[j]]$PowerLaManna$CNDD
    if (j == 1) {
      plot(x, y, main = "Offset-power model", xlim = c(-0.4, 3.4),  ylim = c(-0.5,2.5), col = col[1], pch=pch[1], cex=1, xlab = "Species abundance log10(N/ha)", ylab = "CNDD")
      meanA = seq(min(fig2[[j]]$abundanceA), max(fig2[[j]]$abundanceA), length.out = 20)/10000
      meanB = seq(min(fig2[[j]]$abundanceB), max(fig2[[j]]$abundanceB), length.out = 20)/10000
      estCNDD = log10(meanB + 1)/log10(meanA + 1)
      lines(log10(meanA*50000), estCNDD)
    } else {
      points(x, y, col = col[j], cex=1, pch=pch[j])    
    }
    smo <- gam(y ~ s(x, k=3))
    smo_fit <- predict(smo, newdata = data.frame(x=sort(unique(x))), type = "response")
    lines(sort(unique(x)), smo_fit, col=col_full[j], lwd=2)
   }
}


for(j in 1:6) {
  if(!is.null(fig2[[j]])){
    x <- fig2[[j]]$zerosOnes[fig2[[j]]$RickerLaManna$CNDD > -10]
    y <- fig2[[j]]$RickerLaManna$CNDD[fig2[[j]]$RickerLaManna$CNDD > -10]
    if (j == 1) {
      plot(x, y, main = "",  xlim = c(0,1), ylim = c(-8,2), col = col[1], pch=pch[1], cex=1, xlab = "Proportion recruits without adults", ylab = "CNDD")
    } else {
      points(x, y, col = col[j], cex=1, pch=pch[j])    
    }
    if (j != 1) {
      smo <- gam(y ~ s(x, k=3))
      smo_fit <- predict(smo, newdata = data.frame(x=sort(unique(x))), type = "response")
      lines(sort(unique(x)), smo_fit, col=col_full[j], lwd=2)
  }
  }
}


for(j in 1:6) {
  if(!is.null(fig2[[j]])){
    x <- log10(fig2[[j]]$abundance)
    y <- fig2[[j]]$PowerRandomR$CNDD
    if (j == 1) {
      plot(x, y, main = "", xlim = c(-0.4, 3.4),  ylim = c(-0.5,2.5), col = col[1], pch=pch[1], cex=1, xlab = "Species abundance log10(N/ha)", ylab = "CNDD")
    } else {
      points(x, y, col = col[j], cex=1, pch=pch[j])    
    }
    smo <- gam(y ~ s(x, k=3))
    smo_fit <- predict(smo, newdata = data.frame(x=sort(unique(x))), type = "response")
    lines(sort(unique(x)), smo_fit, col=col_full[j], lwd=2)
  }
}
 
mtext(c("A", "B", "C", "D"), side = rep(3,4), c(-2,-2,-25, -25), at = c(0,0.5, 0, 0.5), outer = T, cex = 1.2, font=2)
```

![](readme_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


**Fig. 2** Estimated CNDD vs. abundance (Log10 N/ha) per species from the Ricker model (**A**) and the offset-power model (**B**) at the 10x10m scale. Gray circles result from data simulated randomly without CNDD and changing spatial association between adults and recruits, ranging from perfect spatial coupling (lightest gray) to no spatial association (black). Blue triangles depict CNDD estimates for the tropical Barro Colorado Island (BCI) forest plot; orange squares the simulation model without CNDD used in the appendix of LaManna *et al*.. For the Ricker model, CNDD bias is highly correlated with the proportion of corrected adult counts (**C**). Fitting a species-specific recruit-to-adult ratio in the offset-power model removes the CNDD-abundance correlation (**D**). For details on the simulation settings see [here](https://github.com/LisaHuelsmann/CommentTo-LaMannaEtAl-Science/tree/master/code). 






