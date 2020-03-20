error
rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/analyze_burgmer/")

#TODO: think about adding multiple species?
#TODO: think about leave-one-out cross-validation at the rep level

## set-up
# load packages
require(rEDM); require(BayesianTools)
require(viridis)
source("../pts_r_package/pttstability/R/bayesfun.R")
source("../pts_r_package/pttstability/R/fake_data.R")
source("../pts_r_package/pttstability/R/logit_funs.R")
source("../pts_r_package/pttstability/R/particlefilter.R")
source("../pts_r_package/pttstability/R/pttstability_man.R")

if(FALSE) {
  # load data
  datfull<-read.table("data/burgmer phyto biovol.csv", header=TRUE, sep=";", stringsAsFactors = FALSE)
  datfull<-datfull[rowSums(!is.na(datfull[,-c(1:3)]),na.rm=T)!=0,]
  
  # get dates
  datedat<-data.frame(date=datfull$date, year=NA, month=NA, day=NA)
  datedat$year<-as.numeric(substr(datedat$date, 7, 10))
  datedat$month<-as.numeric(substr(datedat$date, 4, 5))
  datedat$day<-as.numeric(substr(datedat$date, 1, 2))
  datedat$time<-suppressWarnings(as.POSIXlt(datedat$date, format = "%d.%m.%Y")$yday)+(as.numeric(datedat$year==2006)*365)
  datedat$time<-datedat$time-min(datedat$time)+1
  
  # get treatments
  trtmat<-data.frame(trt=as.character(sort(unique(datfull$treatment))),
                     temperature=NA, variability=NA, predator=NA)
  trtmat$temperature<-substr(trtmat$trt, 1, 1)
  trtmat$variability<-substr(trtmat$trt, 2, 2)
  trtmat$predator<-substr(trtmat$trt, 3, 3)
  spps<-c(4:22)
  collst<-viridis(length(spps))
  
  # make time series for each treatmetn and species
  dat<-data.frame(treatment=datfull$treatment, number=datfull$number, time=datedat$time, datfull[,-c(1:3)])
  dat<-dat[order(dat$treatment, dat$number, dat$time),]
  
  # make plots for each
  par(mfrow=c(4,2), mar=c(4,4,2,2))
  libmat<-NULL
  datnum<-1:nrow(dat)
  doplot<-TRUE
  
  for(i in 1:nrow(trtmat)) {
    ps1<-which(dat$treatment==trtmat$trt[i])
    replst<-sort(unique(dat$number[ps1]))
    
    if(doplot) {
      plot(range(dat$time), range(dat[,spps])+1, xlab="days", ylab="abundance",
           type="n", log="y", axes=F, main=trtmat$trt[i])
      axis(1)
      axis(2, at=(10^(0:8)), labels = c(0, (10^(1:8))), las=2)
      box(); abline(h=1, lty=3)
    }
      
    for(j in 1:length(replst)) {
      ps2<-which(dat$number[ps1]==replst[j])
      if(doplot) {
        matlines(dat[ps1,][ps2,]$time, dat[ps1,][ps2,][,spps]+1, type="b", lty=1, col=collst, pch=1, cex=0.5)
      }
      
      libmat<-rbind(libmat, data.frame(trt=trtmat$trt[i], rep=replst[j], start=min(datnum[ps1][ps2]), end=max(datnum[ps1][ps2])))
    }
  }
  
  ## run particle filter
  trtuse<-"LSA"
  libuse<-as.matrix(libmat[libmat$trt==trtuse,c("start", "end")])
  yps<-which(dat$treatment==trtuse)
  y<-dat[yps,c("Chlamydomonas.terricola", "Ankistrodesmus.sp.","Coccoid.Cyanobacterium")]
  libuse_y<-libuse-min(libuse)+1
  y<-apply(y, 2, function(x) x/sd(x,na.rm=T))
  
  par(mfrow=c(1,1))
  matplot(dat$time[yps], y, xlab="time", ylab="abundance", type="p", col=c(1,2,4), pch=1)
  
  # get EDM parameters
  tuse<-numeric(ncol(y))
  euse<-numeric(ncol(y))
  for(i in 1:ncol(y)) {
    sout<-s_map(y[,i], E=2:4, silent = TRUE, lib = libuse_y)
    tuse[i]<-sout$theta[which.max(sout$rho)]
    euse[i]<-sout$E[which.max(sout$rho)]
  }
  spred<-NULL
  smap_coefs<-NULL
  for(i in 1:ncol(y)) {
    spred[[i]]<-s_map(y[,i], E=euse[i], theta=tuse[i], silent = TRUE, lib = libuse_y, stats_only = FALSE, save_smap_coefficients = TRUE)
    plot(spred[[i]]$model_output[[1]]$obs, pmax(0, spred[[i]]$model_output[[1]]$pred), xlab="obs", ylab="pred"); abline(a=0, b=1, lty=2); abline(h=0, v=0, lty=3)
    print(cor(spred[[i]]$model_output[[1]]$obs, pmax(0, spred[[i]]$model_output[[1]]$pred), use="complete")^2)
    
    smap_coefs[[i]]<-spred[[i]]$smap_coefficients[[1]]
  }
  
  # set priors
  minvUSE_edm<-c(log(0.01), rep(c(log(0.01), log(0.01)), ncol(y)))
  maxvUSE_edm<-c(log(1), rep(c(log(1), log(3)), ncol(y)))
  
  #density, sampler, and prior functions for EDM function
  density_fun_USE_edm<-function(param) density_fun0(param = param, minv = minvUSE_edm, maxv=maxvUSE_edm)
  sampler_fun_USE_edm<-function(x) sampler_fun0(n = 1, minv = minvUSE_edm, maxv=maxvUSE_edm)
  prior_edm <- createPrior(density = density_fun_USE_edm, sampler = sampler_fun_USE_edm,
                           lower = minvUSE_edm, upper = maxvUSE_edm)
  ## Run filter
  niter<-1e4 #number of steps for the MCMC sampler
  N<-2e3 #number of particles
  Euse<-euse #number of embedding dimensions
  
  #likelihood and bayesian set-ups for EDM functions
  likelihood_EDM_piecewise<-function(x) {
    LLtot<-0
    
    for(ii in 1:ncol(y)) {
      xuse<-x[c(1, 2:3+2*(ii-1))]
      tuse_edm<-tuse[ii]
      
      for(i in 1:nrow(libuse_y)) {
        ysegment<-y[libuse_y[i,1]:libuse_y[i,2],ii]
        smap_coefs_segment<-smap_coefs[[ii]][libuse_y[i,1]:libuse_y[i,2],]
        
        LLtot<-LLtot+likelihood0(param = xuse, y=ysegment, parseparam = function(x) parseparam0(x, colparam=c(logit(1e-6), log(0.1))),
                  detfun = EDMfun0, edmdat = list(E=Euse[ii], theta=tuse_edm, smp_cf=smap_coefs_segment), N = N)
      }
    }
    
    return(LLtot)
  }
  
  
  particleFilterLL_piecewise<-function(paramtot, N) {
    pfout<-NULL
    for(ii in 1:ncol(y)) {
      param<-paramtot[c(1, 2:3+2*(ii-1))]
      
      pars<-parseparam0(param, colparam=c(logit(1e-6), log(0.1)))
      tuse_edm<-tuse[ii]
      
      pfout[[ii]]<-data.frame(rep=NA, Nest=rep(NA, length(y[,ii])), Nsd=NA)
      for(i in 1:nrow(libuse_y)) {
        ysegment<-y[libuse_y[i,1]:libuse_y[i,2],ii]
        smap_coefs_segment<-smap_coefs[[ii]][libuse_y[i,1]:libuse_y[i,2],]
        
        tmp<-particleFilterLL(ysegment, pars, N=N, detfun = EDMfun0,
                              edmdat = list(E=Euse[ii], theta=tuse_edm, smp_cf=smap_coefs_segment),
                              dotraceback = TRUE)
        pfout[[ii]]$Nest[libuse_y[i,1]:libuse_y[i,2]]<-tmp$Nest
        pfout[[ii]]$Nsd[libuse_y[i,1]:libuse_y[i,2]]<-tmp$Nest
        pfout[[ii]]$rep[libuse_y[i,1]:libuse_y[i,2]]<-i
      }
    }
    
    return(pfout)
  }
  
  bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM_piecewise, prior = prior_edm)
  
  #run MCMC optimization
  out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=20))
  plot(out_EDM, start = floor(niter/5))
  correlationPlot(out_EDM, start = floor(niter/5))
  
  #get trajectories
  smp_EDM<-getSample(out_EDM, start=floor(niter/5))
  simout<-array(dim=c(nrow(y), nrow(smp_EDM), ncol(y)))
  for(i in 1:nrow(smp_EDM)) {
    tmp<-particleFilterLL_piecewise(smp_EDM[i,], N=N)
    simout[,i,]<-sapply(tmp, function(x) x$Nest)
    
    if(i/10 == floor(i/10)) {
      print(round(i/nrow(smp_EDM),3))
    }
  }
  save.image("tmp_n.rda", version = 2)
} else {
  load("tmp_n.rda")
}

#plot outputs
for(ii in 1:ncol(y)) {
  par(mar=c(4,4,2,1))
  m<-rbind(c(1,1,2,5),
           c(1,1,2,5),
           c(1,1,3,5),
           c(1,1,3,6),
           c(1,1,4,6),
           c(1,1,4,6))
  layout(m)
  plot(range(dat$time[yps]), c(0, 5), xlab="time", ylab="abundance", type="n", xaxs="i", main="Corrected")
  for(i in 1:nrow(libuse_y)) {
    tmp<-t(apply(simout[libuse_y[i,1]:libuse_y[i,2],,ii], 1, function(x) quantile(x, probs = pnorm(c(-2,2)))))
    tm<-dat$time[yps][libuse_y[i,1]:libuse_y[i,2]]
    polygon(c(tm, rev(tm)), c(tmp[,1], rev(tmp[,2])), col=adjustcolor(collst[i], alpha.f = 0.8))
  }
  abline(h=0, lty=3)
  abline(v=(dat$time[yps][libuse_y[1,1]:libuse_y[1,2]])[Euse[ii]], lty=3)
  
  hist(exp(smp_EDM[,1]), main="Obs. Error_B0", xlab="obs", breaks = 20, probability = TRUE)
  abline(v=mean(exp(smp_EDM[,1])), lwd=2, col="blue")
  hist(exp(smp_EDM[,2+2*(ii-1)]), main="Proc. Noise_B0", xlab="proc", breaks = 20, probability = TRUE)
  abline(v=mean(exp(smp_EDM[,2+2*(ii-1)])), lwd=2, col="red")
  hist(exp(smp_EDM[,3+2*(ii-1)]), main="Proc. Noise_B1", xlab="proc", breaks = 20, probability = TRUE)
  abline(v=mean(exp(smp_EDM[,3+2*(ii-1)])), lwd=2, col="red")
  
  plot(range(dat$time[yps]), c(0, 5), xlab="time", ylab="abundance", type="n", xaxs="i", main="Raw")
  for(i in 1:nrow(libuse_y)) {
    lines(dat$time[yps][libuse_y[i,1]:libuse_y[i,2]], y[libuse_y[i,1]:libuse_y[i,2],ii], col=adjustcolor(collst[i], alpha.f = 0.8))
  }
  tmp<-t(apply(simout[,,ii], 1, function(x) quantile(x, probs = pnorm(c(-2,0,2)))))
  abline(h=0, lty=3)
  
  
  plot(tmp[,2], y, xlab="Corrected", ylab="Raw", col=adjustcolor(1, alpha.f = 0.2))
  segments(tmp[,1], y, tmp[,3], y[,ii], lend=2, col=adjustcolor(1, alpha.f = 0.2))
  obsest<-pmax(0.01, y[,ii]*mean(exp(smp_EDM[,1])))
  segments(tmp[,2], y[,ii]-obsest, tmp[,2], y[,ii]+obsest, lend=2, col=adjustcolor(1, alpha.f = 0.2))
  abline(a=0, b=1, lty=2, lwd=2)
  abline(h=0, v=0, lty=3)
}
