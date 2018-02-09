
# Code by Lisa Hülsmann and Florian Hartig
# for a technical comment to LaManna et al 2017 Science





Ricker <- function(A, a, b) {
  list(predictors = list(r = 1, CNDD = 1, HNDDa = 1, HNDDr = 1),
       variables = list(substitute(A), substitute(a), substitute(b)),
       term = function(predictors, variables) {
         pred <- paste("(", variables[1],")*exp(", predictors[1],
                       ")*exp(", predictors[2], "*", variables[1],
                       ")*exp(", predictors[3], "*", variables[2],
                       ")*exp(", predictors[4], "*", variables[3],")+0.0001",
                       sep = "")
       }
  )
}
class(Ricker) <- "nonlin"



getRicker <- function(data, correctZero = F){
  if(correctZero ==T) data$A[data$B > 0 & data$A == 0] = 0.1
  fit = try(gnm(B ~ -1 + Ricker(A, a, b) , family = poisson(link = "identity"), data = data, verbose = F), silent = T)
  out = tryCatch(c(coef(fit), seCNDD = se(fit, checkEstimability = F)[2, 2]), error=function(e) rep(NA, 5))
  if(is.null(out)) out = rep(NA, 5)
  return(out)
}


# extract random effects
se.ranef=function (object) {
  se.bygroup <- ranef(object, condVar = TRUE)
  n.groupings <- length(se.bygroup)
  for (m in 1:n.groupings) {
    vars.m <- attr(se.bygroup[[m]], "postVar")
    K <- dim(vars.m)[1]
    J <- dim(vars.m)[3]
    se.bygroup[[m]] <- array(NA, c(J, K))
    for (j in 1:J) {
      se.bygroup[[m]][j, ] <- sqrt(diag(as.matrix(vars.m[, , j])))
    }
    names.full <- dimnames(se.bygroup$species)
    dimnames(se.bygroup[[m]]) <- list(names.full[[1]], names.full[[2]])
  }
  return(se.bygroup)
}



getPowerOffsetMixed <- function(data, randomR = F){
  if(randomR == F){
    
    fit = try(lmer(log(B+1) ~ log(A+1) + (log(A+1)-1|species) + scale(a) + scale(b), data = data), silent = T)
    ran = tryCatch(ranef(fit), error=function(e) data.frame(species = rep(NA, length(unique(data$species)))))
    fixed = tryCatch(fixef(fit), error=function(e) rep(NA, 4))
    seCNDD = se.ranef(fit)$species[, ncol(ran$species)]
    
    out = data.frame(r = fixed[1], CNDD = fixed[2] + ran$species, HNDDa = fixed[3], HNDDr = fixed[4], seCNDD = seCNDD,
                     row.names = rownames(ran$species))
    colnames(out)[2] = "CNDD"
  }
  if(randomR == T){
    fit = try(lmer(log(B+1) ~ log(A+1) + (log(A+1)|species) + scale(a) + scale(b), data = data, control = lmerControl(optimizer = "bobyqa")), silent = T)
    ran = tryCatch(ranef(fit), error=function(e) list(species = list('log(A + 1)'=rep(NA, length(unique(data$species))),
                                                                     '(Intercept)'=rep(NA, length(unique(data$species))))))
    fixed = tryCatch(fixef(fit), error=function(e) rep(NA, 4))
    seCNDD = se.ranef(fit)$species[, ncol(ran$species)]
    
    out = data.frame(r = fixed[1] + ran$species$`(Intercept)`, CNDD = fixed[2] + ran$species$`log(A + 1)`, HNDDa = fixed[3], HNDDr = fixed[4], seCNDD = seCNDD,
                     row.names = rownames(ran$species))
  }
  return(out)
}


convertMat <- function(x){
  nam = names(x)
  x = matrix(unlist(x), ncol = 5, byrow = T) 
  rownames(x) = nam
  colnames(x) = c("r", "CNDD", "HNDDa", "HNDDr", "seCNDD")
  return(as.data.frame(x))
}



#' @which which model to be used
runAnalysis <- function(which="all", abundanceDist = NULL, externalData = NULL, numQuadrats = NULL, dispersal = Inf, remove = T, theta = NULL, adultSurvival = 1, CNDD = 0, gaps = 1, Aratio = 0.2, ABsd = 0, suitability = 1){
  
  if(! is.null(externalData) & ! is.null(abundanceDist)) stop()
  if(! is.null(externalData))  {
    data = externalData
    abundanceDist = as.vector(by(data, data$species, function(x) sum(x$A + x$B)))
    names(abundanceDist) = levels(data$species)
    numQuadrats = nrow(data) / length(unique(data$species))
  } else  {
    data = createPlotData(abundanceDist, numQuadrats = numQuadrats, dispersal = dispersal, theta = theta, adultSurvival = adultSurvival, CNDD = CNDD, gaps = gaps, Aratio = Aratio, ABsd = ABsd, suitability = suitability)
  }
  
  out = list()
  
  out$data = data
  
  out$abundance = as.vector(by(data, data$species, function(x) sum(x$A) + sum(x$B)))/numQuadrats*100
  out$abundanceA = as.vector(by(data, data$species, function(x) sum(x$A)))/numQuadrats*100
  out$abundanceB = as.vector(by(data, data$species, function(x) sum(x$B)))/numQuadrats*100
  
  if(remove == T){
    incl = levels(data$species)[!(out$abundanceA*numQuadrats/100 < 10 | out$abundanceB*numQuadrats/100 < 10)]
    data = data[data$species %in% incl, ]
    data$species = droplevels(data$species)
    out$dataReduced = data
    
    out$abundance = as.vector(by(data, data$species, function(x) sum(x$A) + sum(x$B)))/numQuadrats*100
    out$abundanceA = as.vector(by(data, data$species, function(x) sum(x$A)))/numQuadrats*100
    out$abundanceB = as.vector(by(data, data$species, function(x) sum(x$B)))/numQuadrats*100
  }  
  
  out$zerosOnes= as.vector(by(data, data$species, function(x) sum(x$B > 0 & x$A == 0 ) / sum(x$B > 0)))
  out$zeros = as.vector(by(data, data$species, function(x) sum(x$B == 0 | x$A == 0) / length(x$B)))
  out$zerosA = as.vector(by(data, data$species, function(x) sum(x$A == 0 ) / length(x$B)))
  out$zerosB = as.vector(by(data, data$species, function(x) sum(x$B == 0 ) / length(x$B)))
  
  if (which == "all") {
    out$RickerLaManna = convertMat(by(data, data$species, getRicker, correctZero =T))
    out$RickerRemoved = convertMat(by(data, data$species, getRicker))
    out$PowerLaManna = getPowerOffsetMixed(data, randomR = F)
    out$PowerRandomR = getPowerOffsetMixed(data, randomR = T)
  } else {
    out$RickerLaManna = convertMat(by(data, data$species, getRicker, correctZero =T))
  }
  
  out$meanCNDD = weighted.mean(x = out$RickerLaManna$CNDD[out$RickerLaManna$CNDD > -15], 
                               w = 1/out$RickerLaManna$seCNDD[out$RickerLaManna$CNDD > -15], na.rm=T)
  out$lmfit = lm(out$RickerLaManna$CNDD[out$RickerLaManna$CNDD > -15] ~ log10(out$abundance[out$RickerLaManna$CNDD > -15]), 
                 weights = 1/out$RickerLaManna$seCNDD[out$RickerLaManna$CNDD > -15])
  
  
  return(out)
}



runAnalyses <- function(which = "all", abundanceDist = NULL, numQuadrats = NULL, externalData = NULL, dispersal = c(0, 0.2, Inf), theta = NULL, adultSurvival = 1, CNDD = 0, gaps = 1, Aratio = 0.2, ABsd = 0, suitability = 1){
  res = list()
  for(i in 1:length(dispersal)){
    print(paste("running" , dispersal[i]))
    res[[i]] = runAnalysis(which=which, abundanceDist = abundanceDist, externalData = externalData, numQuadrats = numQuadrats, dispersal = dispersal[i], theta = theta, adultSurvival = adultSurvival, CNDD = CNDD, gaps = gaps, Aratio = Aratio, ABsd = ABsd, suitability = suitability)
  } 
  
  names(res) = dispersal
  return(res)
}


runAnalysesAbundance <- function(which = "all", abundanceDist = NULL, species = NULL, numQuadrats = 2500, externalData = NULL, dispersal , theta = NULL, adultSurvival, CNDD, gaps = 1, Aratio = 0.2, ABsd = 0, suitability = 1){
  
  
  if(is.null(abundanceDist)) abundanceDist = sapply(species, createAbundance, simplify = F)
  
  res = list()
  pars = matrix(ncol = 3, nrow = length(abundanceDist))
  
  for(i in 1:length(abundanceDist)){
    res[[i]] = runAnalysis(which = which, abundanceDist = abundanceDist[[i]], externalData = externalData, numQuadrats = numQuadrats, dispersal = dispersal[i], theta = theta, adultSurvival = adultSurvival[i], CNDD = CNDD[i], gaps = gaps[i], Aratio = Aratio[i], ABsd = ABsd [i], suitability = suitability[i])
    pars[i,1] = res[[i]]$meanCNDD
    pars[i,2:3] = coef(res[[i]]$lmfit)
  } 
  colnames(pars) = c("mean", "intercept", "slope")
  return(list(summary = pars, full = res))
}






