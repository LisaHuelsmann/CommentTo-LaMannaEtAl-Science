
# Code by Lisa Hülsmann and Florian Hartig
# for a technical comment to LaManna et al 2017 Science




createSpeciesData <- function(numQuadrats, abundance, dispersal, theta = NULL, adultSurvival = 1, CNDD = 0, gaps = 1, Aratio = 0.2, ABsd = 0, suitability = 1) {
  
  numQuadratsTrue = ifelse(suitability==1, numQuadrats, round(suitability * numQuadrats))

  if (adultSurvival == "LaManna"){
    adultSurvival = 1
    LaManna = T
  } else LaManna = F
  
  # local adult abundance (considering that some will die later on)
  abundance = abundance / numQuadratsTrue * Aratio * (1/adultSurvival) # does not consider ABsd...
  
  ABprod = (1-Aratio) / Aratio
  ABratio = rnorm(1, ABprod * adultSurvival * 1/ gaps, sd = ABsd)  
  
  if(is.null(theta)){
    if(dispersal == Inf | dispersal == 1){
      A = rpois(numQuadratsTrue, abundance)
      B = rpois(numQuadratsTrue, ABratio * abundance)    
    }else if(dispersal == 0){
      A = rpois(numQuadratsTrue, abundance)
      B = rpois(numQuadratsTrue, ABratio * A)
    }else if(dispersal > 0){
      A = rpois(numQuadratsTrue, abundance)
      temp = ABratio * A
      temp2 = ABratio * A * (1-dispersal)
      B = rpois(numQuadratsTrue, temp2 + sum(temp)*dispersal/numQuadratsTrue)
    }else stop("wrong argument")    
  }else{
    if(dispersal == Inf | dispersal == 1){
      A = rnbinom(numQuadratsTrue, size = theta, mu=abundance)
      B = rnbinom(numQuadratsTrue, size = theta, mu=ABratio * abundance)    
    }else if(dispersal == 0){
      A = rnbinom(numQuadratsTrue, size = theta, mu=abundance)
      B = rnbinom(numQuadratsTrue, size = theta, mu=ABratio * A)
    }else if(dispersal > 0){
      A = rnbinom(numQuadratsTrue, size = theta, mu=abundance)
      temp = ABratio* A
      temp2 = ABratio* A * (1-dispersal)
      B = rnbinom(numQuadratsTrue, size = theta, mu=temp2 + sum(temp)*dispersal/numQuadratsTrue)
    }else stop("wrong argument")  
    
  }
  
  # this will kill all recruits that are not in a gap
  if(gaps < 1){
    sel = sample.int(numQuadratsTrue, numQuadratsTrue * (1-gaps))
    B[sel] = 0
  }
  
  if(LaManna == T) A = rpois(length(A), A)
  else if(adultSurvival < 1){
    A = rbinom(length(A), size = A, prob = adultSurvival)
  } 

  
  if (CNDD < 0) {
    survival = exp(CNDD * A)
    B0 = B
    B = rbinom(length(B), size = B, prob = survival)  
    
    # to maintain abundance, we add randomly placed individuals 
    dens = sum(B0-B) / numQuadratsTrue
    B = B + rpois(length(B), dens)
  }
  
  dat = data.frame(A=rep(0, numQuadrats), B=rep(0, numQuadrats))
  dat[sample.int(n = numQuadrats, size = numQuadratsTrue), ] <- data.frame(A, B)
  return(dat)
}




#' @param abundanceDist named vector with abundance per species on the entire plot
createPlotData <- function(abundanceDist , numQuadrats , dispersal, theta = NULL, adultSurvival = 1, CNDD = 0, gaps = 1, Aratio = 0.2, ABsd = 0, suitability = 1){
  
  if(is.null(names(abundanceDist))) names(abundanceDist) = 1:length(abundanceDist)
  
  nObs = length(abundanceDist) * numQuadrats

  temp = data.frame(species = rep(names(abundanceDist), each = numQuadrats), 
                    quadrat = rep(1:numQuadrats, times = length(abundanceDist)), 
                    A = rep(NA, nObs), B = rep(NA, nObs))
  for(i in 1:length(abundanceDist)){
    counts = abundanceDist[i]
    temp[((i-1) * numQuadrats +1) :(i*numQuadrats), c("A", "B")] = createSpeciesData(abundance = counts, numQuadrats = numQuadrats, dispersal = dispersal, theta = theta, adultSurvival = adultSurvival, CNDD = CNDD, gaps = gaps, Aratio = Aratio, ABsd = ABsd, suitability = suitability)
  }    
  
  adult_all = rep(by(temp$A, temp$quadrat, sum), length(abundanceDist))
  recruit_all = rep(by(temp$B, temp$quadrat, sum), length(abundanceDist))
  
  # calculate heterospecific density
  temp$b <- recruit_all - temp$B
  temp$a <- adult_all - temp$A
  
  return(temp)
}




createAbundance <- function(species, total = 130000){
  
  abundance = rep(1, species)
  
  decay = 20 /(species + 100)
  prob = dexp(1:species, rate = decay)
  
  abundance = abundance + as.vector(rmultinom(1, total-species, prob = prob))
  
  names(abundance) = 1:species
  return(abundance)
}








