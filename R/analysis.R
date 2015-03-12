
devtools::load_all("C:/work/repos/Faskally/rtemp")

setwd("C:/work/repos/papers/rtemp_paper")


# get basis transformation

write.table(subset(stream, site == 10), file = "Burn10.csv", sep = ",", row.names = FALSE)
write.table(subset(stream, site == 11), file = "Burn11.csv", sep = ",", row.names = FALSE)

# set up a fourier basis of dimension 48 and order 6 
hourbasis <- create.fourier.basis(c(0, 24), 47, 24)
# and a harmonic accelerator penalty... penalises deviations from sinusiodal curves
penalty <- vec2Lfd(c(0, (2*pi/24)^2, 0), c(0, 24)) 

# would be good to use the same bases for all experiments...

for (choice in 1:5) {

cat("choice:", choice, "\n"); flush.console()

if (choice == 1) {
  years <- 2003:2009
  site2 <- 2
  site10 <- 10
} else if (choice == 2) {
  years <- 2003:2009
  site2 <- 2
  site10 <- 11
} else if (choice == 3) {
  years <- 1988:1994
  site2 <- 11
  site10 <- 10  
} else if (choice == 4) {
  years <- 1996:2002
  site2 <- 10
  site10 <- 11  
} else if (choice == 5) {
  years <- 1993:1995
  site2 <- 11
  site10 <- 7  
}


# select a few days data and add in decimal hour
dat <- with(subset(stream, site == site2 & year %in% years), 
            data.frame(temp = Original.Value, 
                       dhour = hour + min/60, 
                       month = factor(month.abb[mon + 1], levels = month.abb), 
                       yday = yday, 
                       year = year,  
                       site = site)) 
dat <- subset(dat, yday < 365) # remove last day of leap year

getcoef <- function(lambda, .yday, .year = 2009, .site = site2) {
  sdat <- subset(dat, yday == .yday & year == .year & site == .site)
  #if (nrow(sdat) < 23) return (NULL)
  if (sum(unique(round(sdat $ dhour)) %in% 0:24) < 24) return(NULL)  
  c(coef(smooth.basis(sdat $ dhour, sdat $ temp, fdPar(hourbasis, penalty, lambda))))
}

yday <- rep(0:365, length(years))
year <- rep(years, each = 366)

ck <- mapply(getcoef, .yday = yday, .year = year, MoreArgs = list(lambda = exp(5)))
whichNULL <- sapply(ck, is.null)
ck <- simplify2array(ck[!whichNULL])

info2 <-data.frame(yday = yday[!whichNULL], year = year[!whichNULL])

tempfd2 <- fd(ck, hourbasis)
pc2 <- pca.fd(tempfd2, nharm = 4) # is this sensitive to missing values?...

par(mfrow = c(2,2))
plot(pc2)

info2[pc2 $ harmonics $ fdnames[[2]]] <- as.data.frame(pc2 $ scores)



# select a few days data and add in decimal hour
dat <- with(subset(stream, site == site10 & year %in% years), 
            data.frame(temp = Original.Value, 
                       dhour = hour + min/60, 
                       month = factor(month.abb[mon + 1], levels = month.abb), 
                       yday = yday, 
                       year = year, 
                       site = site))  
dat <- subset(dat, yday < 365) # remove last day of leap year

getcoef <- function(lambda, .yday, .year = 2009, .site = site10) {
  sdat <- subset(dat, yday == .yday & year == .year & site == .site)
  #if (nrow(sdat) < 23) return (NULL)
  if (sum(unique(round(sdat $ dhour)) %in% 0:24) < 24) return(NULL)  
  c(coef(smooth.basis(sdat $ dhour, sdat $ temp, fdPar(hourbasis, penalty, lambda))))
  #c(coef(smooth.basis(sdat $ dhour, sdat $ temp, pc2 $ harmonics)))
}

yday <- rep(0:365, length(years))
year <- rep(years, each = 366)

ck <- mapply(getcoef, .yday = yday, .year = year, MoreArgs = list(lambda = exp(5)))
whichNULL <- sapply(ck, is.null)
ck <- simplify2array(ck[!whichNULL])

info10 <-data.frame(yday = yday[!whichNULL], year = year[!whichNULL])

tempfd10 <- fd(ck, hourbasis)
pc10 <- pca.fd(tempfd10, nharm = 4) # is this sensitive to missing values?...

par(mfrow = c(2,2))
plot(pc10)

info10[pc10 $ harmonics $ fdnames[[2]]] <- as.data.frame(pc10 $ scores)

info <- merge(info2, info10, by = c("yday", "year"), suffixes = c(".2", ".10"))
info $ day <- with(info, yday + (year - 2003) * 365)
info $ month <- month.abb[strptime(paste(info $ yday + 1, info $ year), format = "%j %Y") $ mon + 1]
info $ season <- with(info, ifelse(month %in% c("Nov", "Dec", "Jan", "Feb"), "Winter", 
                                   ifelse(month %in% c("Mar", "Apr"), "Spring",
                                          ifelse(month %in% c("May", "Jun", "Jul", "Aug"), "Summer", "Autumn"))))
info $ season <- factor(info $ season, levels = c("Winter", "Spring", "Summer", "Autumn"))
info <- info[order(info $ day),]  

names(info) <- gsub("[.]2", ".x", names(info))
names(info) <- gsub("[.]10", ".y", names(info))


}







# load required libraries
library(lattice)
library(fda)
#library(reshape)
library(mgcv)
library(MASS)
library(nlme)


# get the data for each site:


# get full data set

choices <- data.frame(burn = c(10, 11, 7, 10, 11), year = c(1989, 1997, 1994, 2004, 2004))


for (ichoices in 1:5) {
  
  fname <- paste0("Z:/Loch_Ard_Felling/BHS/data/B", choices $ burn[ichoices], "_", choices $ year[ichoices], ".rda")
  load(fname)
  
  #load("Z:/Loch_Ard_Felling/BHS/data/B10_2004.rda")
  #load("Z:/Loch_Ard_Felling/BHS/data/B11_2004.rda")
  #load("Z:/Loch_Ard_Felling/BHS/data/B7_1994.rda")
  #load("Z:/Loch_Ard_Felling/BHS/data/B10_1989.rda")
  #load("Z:/Loch_Ard_Felling/BHS/data/B11_1997.rda")
  
  # sort out day
  info $ day <- info $ day - min(info $ day) + info $ yday[1]
  
  
  #plot.pcs(info, main = "All data") 
  
  #xyplot(PC1.y ~ yday | factor(year), data = info, as.table = TRUE)
  
  decays <- seq(0, 1, length = 11)
  gcv1 <- sapply(decays, fitDecay, gcv = TRUE, data = within(info, {y = PC1.y}))
  gcv2 <- sapply(decays, fitDecay, gcv = TRUE, data = within(info, {y = PC2.y}))
  gcv3 <- sapply(decays, fitDecay, gcv = TRUE, data = within(info, {y = PC3.y}))
  gcv4 <- sapply(decays, fitDecay, gcv = TRUE, data = within(info, {y = PC4.y}))
  gcvs <- gcv1 + gcv2 + gcv3 + gcv4
  
  decays2 <- seq(decays[max(1, which.min(gcvs) - 1)], decays[max(1, which.min(gcvs) + 1)], length = 11)
  gcv1 <- sapply(decays2, fitDecay, gcv = TRUE, data = within(info, {y = PC1.y}))
  gcv2 <- sapply(decays2, fitDecay, gcv = TRUE, data = within(info, {y = PC2.y}))
  gcv3 <- sapply(decays2, fitDecay, gcv = TRUE, data = within(info, {y = PC3.y}))
  gcv4 <- sapply(decays2, fitDecay, gcv = TRUE, data = within(info, {y = PC4.y}))
  gcvs2 <- gcv1 + gcv2 + gcv3 + gcv4
  
  gcvs <- c(gcvs, gcvs2)
  decays <- c(decays, decays2)
  
  plot(decays, gcvs)
  
  decay <- decays[which.min(gcvs)]
  full1 <- fitDecay(decay, data = within(info, {y = PC1.y}))
  full2 <- fitDecay(decay, data = within(info, {y = PC2.y}))
  full3 <- fitDecay(decay, data = within(info, {y = PC3.y}))
  full4 <- fitDecay(decay, data = within(info, {y = PC4.y}))
  
  info $ fell <- feffect(info $ day, decay = decay, type = 1)
  xyplot(fell ~ day, data = info)
  
  X <- eval.fd(seq(0, 24, length = 96), pc10 $ harmonics)
  Beta <- rbind(0, getFelledPC(full1), getFelledPC(full2), getFelledPC(full3), getFelledPC(full4))
  
  Beta.all <- rbind(1, getFelledPC(full1, which = NULL), getFelledPC(full2, which = NULL), getFelledPC(full3, which = NULL), getFelledPC(full4, which = NULL))
  
  Xmu <- eval.fd(seq(0, 24, length = 96), pc10 $ mean)
  X <- cbind(Xmu, X)
  
  
  fits <- X %*% Beta
  fits.all <- X %*% Beta.all
  
  # we can compute daily summaries of felling effects
  info $ dailyMax <- apply(fits, 2, max)
  info $ dailyMax.all <- apply(fits.all, 2, max)
  
  
  # what about confidence intervals...
  n <- 1000
  sim1 <- getFelledPC(full1, sim = TRUE, n = n)
  sim2 <- getFelledPC(full2, sim = TRUE, n = n)
  sim3 <- getFelledPC(full3, sim = TRUE, n = n)
  sim4 <- getFelledPC(full4, sim = TRUE, n = n)
  
  simfits <- sapply(1:n, function(i) X %*% rbind(0, sim1[,i], sim2[,i], sim3[,i], sim4[,i]))
  dim(simfits) <- c(nrow(X), ncol(Beta), n)
  
  # we can compute daily summaries of felling effects
  simdailyMax <- apply(simfits, 2:3, max)
  info $ dailyMax.cil <- apply(simdailyMax, 1, quantile, .025)
  info $ dailyMax.ciu <- apply(simdailyMax, 1, quantile, .975)
  
  
  p <- xyplot(dailyMax + dailyMax.cil  + dailyMax.ciu ~ yday | factor(year), data = info, 
         pch = 19, cex = .7, as.table = TRUE, grid = TRUE, type = c("p"),
         main = paste0(choices $ burn[ichoices], "_", choices $ year[ichoices]))
  
  plot(p)

  info $ event <- do.call(paste, c(choices[ichoices,], list(sep = "_"))) 
  
  
  vname <- do.call(paste, c(choices[ichoices,], list(sep = "_"))) 
  assign(paste0("info", vname), info, envir = .GlobalEnv)
  #assign(paste0("simfits", vname), simfits, envir = .GlobalEnv)
  
  fname <- paste0("Z:/Loch_Ard_Felling/BHS/data/B", choices $ burn[ichoices], "_", choices $ year[ichoices], "fits.rda")
  save(list = paste0(c("info"), vname),  file = fname)
  
}


info.all <- rbind(info10_1989, info10_2004, info11_1997, info11_2004, info7_1994)


xyplot(dailyMax ~ yday | year * event, data = info.all)
fname <- paste0("Z:/Loch_Ard_Felling/BHS/data/B", choices $ burn[ichoices], "_", choices $ year[ichoices], "fits.rda")
save(info.all,  file = "Z:/Loch_Ard_Felling/BHS/data/infoall.rda")


































