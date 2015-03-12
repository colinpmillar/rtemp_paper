
devtools::load_all("C:/work/repos/Faskally/rtemp")

setwd("C:/work/repos/papers/rtemp_paper")



# Set up a problem
years <- 2003:2009
control_site <- 2
felled_site <- 10
felled_year <- 2004

# select a few days data and add in decimal hour
dat <- with(subset(stream, site %in% c(control_site, felled_site) & year %in% years), 
            data.frame(temp = Original.Value, 
                       dhour = hour + min/60, 
                       month = factor(month.abb[mon + 1], levels = month.abb), 
                       yday = yday, 
                       year = year,  
                       site = site)) 
# 
dat <- cleanData(dat)

# set up a new set of bases
reddat <- reduceData(dat[1:200,])


# a squiz at the PCs
par(mfrow = c(2,2))
plot(attr(reddat, "pc"))

    



  # sort out day... not sure why...
  reddat $ day <- reddat $ day - min(reddat $ day) + reddat $ yday[1]
  
  #plot_pcs(reddat, main = "All data") 
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


































