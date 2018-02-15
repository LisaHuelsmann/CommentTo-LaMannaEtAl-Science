# Preparation of BCI data
Lisa Hülsmann & Florian Hartig  
February 15 2018  


This code was used to generate counts of conspecific and heterospecific adults and recruits from BCI data.

The data were requested at http://ctfs.si.edu/webatlas/datasets/bci/. In accordance with LaManna *et al*., we used the 7th survey at BCI from http://ctfs.si.edu/ctfsrep. We believe that this is the same version of the BCI data that was used in LaManna *et al*..

*Note that we cannot provide the original dataset in this repository due to the data sharing agreement of the BCI data.* 

It must be requested from the data owners via http://ctfs.si.edu/webatlas/datasets/bci/



## Load packages



```r
library(spatstat)
library(dplyr)
```


## Read BCI data



```r
dat <- read.csv("bci7.txt", sep="\t", header=T, stringsAsFactors = F)
```


## Modify structure of the dataset

Change column names and coding (for applicability with code by LaManna *et al*.)

```r
dat$gx <- dat$PX
dat$gy <- dat$PY
dat$dbh <- dat$DBH
dat$species <- toupper(dat$Mnemonic)
```

Remove multiple stems (following LaManna *et al*.)

```r
dat <- dat[order(dat$Quadrat, dat$Tag, -dat$dbh), ]
dat <- dat[!duplicated(dat$Tag), ]
```

Remove Anaxagorea panamensis, a species which reproduces clonally (following LaManna *et al*.)

```r
dat <- dat[dat$species != "ANAXPA", ]
```

Select living trees

```r
dat$status <- ifelse(dat$Status == "alive", "A", "U")
dat_liv <- dat[dat$status=="A" & !is.na(dat$dbh), ]
```


## Functions to generate counts

The counts are named as follows:

* Conspecific adults *A* and recruits *B*
* Heterospecific adults *a* and recruits *b*


Generate tree status 

```r
assign_recruit <- function(dat_liv) {
  dat_liv$adultcutoff <- 100
  
  # recruitment based on rules by LaManna et al.
  dat_liv$type1 <- ifelse(dat_liv$dbh < 100, "recruit", "adult")
  dat_check <- dat_liv %>%
    group_by(species) %>% 
    summarise(perc_adult = sum(type1=="adult")/length(type1)) %>% 
    filter(perc_adult < 0.2)
  dat_liv$type1[dat_liv$species %in% dat_check$species] <- ifelse(dat_liv$dbh[dat_liv$species %in% dat_check$species] <= 50, "recruit", "adult")
  dat_liv$adultcutoff[dat_liv$species %in% dat_check$species] <- 50
  
  dat_check <- dat_liv %>%
    group_by(species) %>% 
    summarise(perc_adult = sum(type1=="adult")/length(type1)) %>% 
    filter(perc_adult < 0.2)
  dat_liv$type1[dat_liv$species %in% dat_check$species] <- ifelse(dat_liv$dbh[dat_liv$species %in% dat_check$species] <= 20, "recruit", "adult")
  dat_liv$adultcutoff[dat_liv$species %in% dat_check$species] <- 20
  
  return(dat_liv)
}
```

Counts of all adults and recruits per quadrat to get heterospecific as the difference between adult_all - A. Recruits, respectively.

```r
count_quadrat_all <- function(dat_liv, res) {
  dat.ppp <- ppp(dat_liv$gx, dat_liv$gy, window = owin(c(0, 1000), c(0, 500)), marks = dat_liv$type1)
  # plot(density(dat.ppp), cex=0.1)
  split_recruit <- split.ppp(dat.ppp, as.factor(dat.ppp$marks))
  adult <- split_recruit$adult
  recruit <- split_recruit$recruit
  
  count_adult <- as.numeric(quadratcount(adult, xbreaks = seq(0, 1000, res), ybreaks = seq(0, 500, res)))
  count_recruit <- as.numeric(quadratcount(recruit, xbreaks = seq(0, 1000, res), ybreaks = seq(0, 500, res)))  
  
  count_all <- data.frame(recruit_all = count_recruit, adult_all = count_adult)
  return(count_all)
}
```

Count adults and recruits per species

```r
count_quadrat_sp <- function(act.sp, dat_liv, count_all, res) {
  dat.sp <- dat_liv[dat_liv$species == act.sp, ]
  dat.ppp <- ppp(dat.sp$gx, dat.sp$gy, window = owin(c(0, 1000), c(0, 500)), marks = dat.sp$type1)
  # plot(density(dat.ppp), cex=0.1)
  split_recruit <- split.ppp(dat.ppp, as.factor(dat.ppp$marks))
  adult <- split_recruit$adult
  recruit <- split_recruit$recruit
  
  if (is.null(recruit)) count_recruit <- rep(0, 500*1000/(res*res)) else {
    count_recruit <- as.numeric(quadratcount(recruit, xbreaks = seq(0, 1000, res), ybreaks = seq(0, 500, res)))  
  }
  if (is.null(adult)) count_adult <- rep(0, 500*1000/(res*res)) else {
    count_adult <- as.numeric(quadratcount(adult, xbreaks = seq(0, 1000, res), ybreaks = seq(0, 500, res)))
  }
  count_sp <- data.frame(B = count_recruit, A = count_adult)
  count_sp$b <- count_all$recruit_all - count_sp$B
  count_sp$a <- count_all$adult_all - count_sp$A
  count_sp$species <- act.sp
  
  return(count_sp)
}
```


## Generate count data for BCI data


Assign recruit and adult

```r
dat_liv <- assign_recruit(dat_liv)
```

Count trees in quadrats

```r
count_all <- count_quadrat_all(dat_liv, res=10)  # all trees per quadrat
```

```
## Warning: data contain duplicated points
```

```r
species <- sort(unique(dat_liv$species))
count_sp <- lapply(species, count_quadrat_sp, dat_liv = dat_liv, count_all = count_all, res = 10)  # per species while a = adult_all - A
count_sp <- do.call(rbind.data.frame, count_sp)
```

Clean count (following LaManna *et al*.)

```r
quadrat_count <- count_sp %>%
  group_by(species) %>%
  summarise(Nquadrat_B = sum(B > 0),
            Nquadrat_A = sum(A > 0),
            Nquadrat_A_and_B = sum(B > 0 & A > 0)) %>%
  filter(Nquadrat_B < 10 | Nquadrat_A < 10 | Nquadrat_A_and_B == 0)

bciCounts <- count_sp[!count_sp$species %in% quadrat_count$species, ]
bciCounts <- bciCounts[, c("species", "A", "B", "a", "b")]
```


## Save BCI counts of conspecific and heterospecific adults and recruits


```r
save(bciCounts, file = "bciCounts.Rdata")
```



