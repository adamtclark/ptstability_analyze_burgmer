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
  #y<-dat$Chlamydomonas.terricola[yps]
  y<-dat$Fragilaria.capucina[yps]
  libuse_y<-libuse-min(libuse)+1
  y<-y/sd(y)
  
  par(mfrow=c(1,1))
  plot(dat$time[yps], y, xlab="time", ylab="abundance", type="p")
  
  # get EDM parameters
  sout<-s_map(y, E=2:4, silent = TRUE, lib = libuse_y)
  tuse<-sout$theta[which.max(sout$rho)]
  euse<-sout$E[which.max(sout$rho)]
  
  spred<-s_map(y, E=euse, theta=tuse, silent = TRUE, lib = libuse_y, stats_only = FALSE, save_smap_coefficients = TRUE)
  plot(spred$model_output[[1]]$obs, pmax(0, spred$model_output[[1]]$pred), xlab="obs", ylab="pred"); abline(a=0, b=1, lty=2); abline(h=0, v=0, lty=3)
  cor(spred$model_output[[1]]$obs, pmax(0, spred$model_output[[1]]$pred), use="complete")^2
  
  # set priors
  #minvUSE_edm<-c(-5, -5, -5, -5, -5) #minimum interval for obs, proc, and theta
  #maxvUSE_edm<-c(0, 0, 0, 1, 2) #maximum interval for obs, proc, and theta
  minvUSE_edm<-c(log(0.01), log(0.01), log(0.01))
  maxvUSE_edm<-c(log(1), log(1), log(3))
  
  #density, sampler, and prior functions for EDM function
  density_fun_USE_edm<-function(param) density_fun0(param = param, minv = minvUSE_edm, maxv=maxvUSE_edm)
  sampler_fun_USE_edm<-function(x) sampler_fun0(n = 1, minv = minvUSE_edm, maxv=maxvUSE_edm)
  prior_edm <- createPrior(density = density_fun_USE_edm, sampler = sampler_fun_USE_edm,
                           lower = minvUSE_edm, upper = maxvUSE_edm)
  ## Run filter
  niter<-1e4 #number of steps for the MCMC sampler
  N<-2e3 #number of particles
  Euse<-euse #number of embedding dimensions
  
  smap_coefs<-spred$smap_coefficients[[1]]
  
  #likelihood and bayesian set-ups for EDM functions
  likelihood_EDM_piecewise<-function(x) {
    xuse<-x#[1:(length(x)-1)]
    tuse_edm<-tuse#exp(x[length(x)])
    
    #smap_out<-s_map(y, E=Euse, theta=tuse_edm, silent = TRUE, lib = libuse_y, save_smap_coefficients = TRUE)
    #smap_coefs<-smap_out$smap_coefficients[[1]]
    
    LLtot<-0
    
    for(i in 1:nrow(libuse_y)) {
      ysegment<-y[libuse_y[i,1]:libuse_y[i,2]]
      smap_coefs_segment<-smap_coefs[libuse_y[i,1]:libuse_y[i,2],]
      
      LLtot<-LLtot+likelihood0(param = xuse, y=ysegment, parseparam = function(x) parseparam0(x, colparam=c(logit(1e-6), log(0.1))),
                detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse_edm, smp_cf=smap_coefs_segment, ytot=y), N = N)
    }
    
    return(LLtot)
  }
  
  
  particleFilterLL_piecewise<-function(param, N) {
    pars<-parseparam0(param, colparam=c(logit(1e-6), log(0.1)))
      #parseparam0(param[1:(length(param)-1)], colparam=c(logit(1e-6), log(0.1)))
    tuse_edm<-tuse#exp(param[length(param)])
    
    #smap_out<-s_map(y, E=Euse, theta=tuse_edm, silent = TRUE, lib = libuse_y, save_smap_coefficients = TRUE)
    #smap_coefs<-smap_out$smap_coefficients[[1]]
    
    pfout<-data.frame(rep=NA, Nest=rep(NA, length(y)), Nsd=NA, Nsmp=NA, Nsmp_noproc=NA)
    for(i in 1:nrow(libuse_y)) {
      ysegment<-y[libuse_y[i,1]:libuse_y[i,2]]
      smap_coefs_segment<-smap_coefs[libuse_y[i,1]:libuse_y[i,2],]
      
      tmp<-particleFilterLL(ysegment, pars, N=N, detfun = EDMfun0,
                            edmdat = list(E=Euse, theta=tuse_edm, smp_cf=smap_coefs_segment),
                            dotraceback = TRUE, fulltraceback = TRUE)
      pfout$Nest[libuse_y[i,1]:libuse_y[i,2]]<-tmp$Nest
      pfout$Nsd[libuse_y[i,1]:libuse_y[i,2]]<-tmp$Nsd
      pfout$rep[libuse_y[i,1]:libuse_y[i,2]]<-i
      pfout$Nsmp[libuse_y[i,1]:libuse_y[i,2]]<-c(indexsort(tmp$fulltracemat, tmp$fulltraceindex, nsmp=1))
      pfout$Nsmp_noproc[libuse_y[i,1]:libuse_y[i,2]]<-c(indexsort(tmp$fulltracemat_noproc, tmp$fulltraceindex, nsmp=1))
    }
    
    return(pfout)
  }
  
  
  bayesianSetup_EDM <- createBayesianSetup(likelihood = likelihood_EDM_piecewise, prior = prior_edm)
  
  #run MCMC optimization
  out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=20))
  plot(out_EDM, start = floor(niter/5))
  correlationPlot(out_EDM, start = floor(niter/5))
  
  #get trajectories
  smp_EDM<-(getSample(out_EDM, start=floor(niter/5)))
  simout<-matrix(nrow=length(y), ncol=nrow(smp_EDM))
  sdout<-matrix(nrow=length(y), ncol=nrow(smp_EDM))
  simout_smp<-matrix(nrow=length(y), ncol=nrow(smp_EDM))
  simout_smp_noproc<-matrix(nrow=length(y), ncol=nrow(smp_EDM))
  for(i in 1:nrow(smp_EDM)) {
    tmp<-particleFilterLL_piecewise(smp_EDM[i,], N=N)
    simout[,i]<-tmp$Nest
    sdout[,i]<-tmp$Nsd
    simout_smp[,i]<-tmp$Nsmp
    simout_smp_noproc[,i]<-tmp$Nsmp_noproc
    
    if(i/10 == floor(i/10)) {
      print(round(i/nrow(smp_EDM),3))
    }
  }
  #save.image("tmp3.rda", version = 2)
  #save.image("tmp2.rda", version = 2)
  #save.image("tmp.rda", version = 2)
} else {
  load("tmp3.rda")
  #load("tmp2.rda")
  #load("tmp.rda")
}

#plot outputs
par(mar=c(4,4,2,1))
m<-rbind(c(1,1,2,5),
         c(1,1,2,5),
         c(1,1,3,5),
         c(1,1,3,6),
         c(1,1,4,6),
         c(1,1,4,6))
layout(m)
ymax<-6
fullerror<-FALSE
tracedparticles<-TRUE

plot(range(dat$time[yps]), c(0, ymax), xlab="time", ylab="abundance", type="n", xaxs="i", main="Corrected")

for(i in 1:nrow(libuse_y)) {
  if(fullerror) {
    tmpsmp[]<-pmax(0, rnorm(n = length(simout[libuse_y[i,1]:libuse_y[i,2],]),
                  mean = simout[libuse_y[i,1]:libuse_y[i,2],],
                  sd = sdout[libuse_y[i,1]:libuse_y[i,2],]))
  } else {
    if(tracedparticles) {
      tmpsmp<-simout_smp[libuse_y[i,1]:libuse_y[i,2],]
      psps<-which(is.na(tmpsmp[]))
      tmpsmp[psps]<-simout[libuse_y[i,1]:libuse_y[i,2],][psps]
    } else {
      tmpsmp<-simout[libuse_y[i,1]:libuse_y[i,2],]
    }
  }
  tmp<-t(apply(tmpsmp, 1, function(x) quantile(x, probs = pnorm(c(-2:2)))))
  tm<-dat$time[yps][libuse_y[i,1]:libuse_y[i,2]]
  polygon(c(tm, rev(tm)), c(tmp[,2], rev(tmp[,4])), col=adjustcolor(collst[i], alpha.f = 0.8))
  i<-i+1
}
abline(h=0, lty=3)
abline(v=(dat$time[yps][libuse_y[1,1]:libuse_y[1,2]])[Euse], lty=3)

hist(exp(smp_EDM[,1]), main="Obs. Error_B0", xlab="obs", breaks = 20, probability = TRUE)
abline(v=mean(exp(smp_EDM[,1])), lwd=2, col="blue")
hist(exp(smp_EDM[,2]), main="Proc. Noise_B0", xlab="proc", breaks = 20, probability = TRUE)
abline(v=mean(exp(smp_EDM[,2])), lwd=2, col="blue")
hist(exp(smp_EDM[,3]), main="Proc. Noise_B1", xlab="proc", breaks = 20, probability = TRUE)
abline(v=mean(exp(smp_EDM[,3])), lwd=2, col="blue")

plot(range(dat$time[yps]), c(0, ymax), xlab="time", ylab="abundance", type="n", xaxs="i", main="Raw")
for(i in 1:nrow(libuse_y)) {
  lines(dat$time[yps][libuse_y[i,1]:libuse_y[i,2]], y[libuse_y[i,1]:libuse_y[i,2]], col=adjustcolor(collst[i], alpha.f = 0.8))
}
abline(h=0, lty=3)

if(fullerror) {
  tmp<-matrix(nrow=nrow(simout), ncol=ncol(simout))
  tmp[]<-rnorm(length(simout), simout, sdout)
  tmp<-t(apply(tmp, 1, function(x) quantile(x, probs = pnorm(c(-1,0,1)))))
} else {
  if(tracedparticles) {
    tmp<-t(apply(simout_smp, 1, function(x) quantile(x, probs = pnorm(c(-1,0,1)), na.rm=T)))
  } else {
    tmp<-t(apply(simout, 1, function(x) quantile(x, probs = pnorm(c(-1,0,1)))))
  }
}

plot(tmp[,2], y, xlab="Corrected", ylab="Raw", col=adjustcolor(1, alpha.f = 0.2), xlim=c(0, ymax), ylim=c(0, ymax))
segments(tmp[,1], y, tmp[,3], y, lend=2, col=adjustcolor(1, alpha.f = 0.2))
obsest<-pmax(0.01, y*mean(exp(smp_EDM[,1])))
segments(tmp[,2], y-obsest, tmp[,2], y+obsest, lend=2, col=adjustcolor(1, alpha.f = 0.2))
abline(a=0, b=1, lty=2, lwd=2)
abline(h=0, v=0, lty=3)






#Use inverse process function
#Make sure to use full set of outcomes to calculate p-values. 

#Think about last time point...

plot(c(0,1), c(0, 1), ylab="pmor[obs]", xlab="pmor[sim]", type="n")
for(i in 1:nrow(libuse_y)) {
  tmpsmp<-simout_smp_noproc[libuse_y[i,1]:libuse_y[i,2],]
  tmpsmp_proc<-simout_smp[libuse_y[i,1]:libuse_y[i,2],]
  
  plot(rowMeans(tmpsmp_proc))
  plot(rowMeans(simout[libuse_y[i,1]:libuse_y[i,2],]))
  
  tmp<-apply(tmpsmp, 1, function(x) mean(procfun0(c(mean(smp_EDM[,2]), mean(smp_EDM[,3])),x,inverse = TRUE)))
  prod(1-tmp[!is.na(tmp)]) #pr of extinction
  mean(tmpsmp_proc[17,]==0) #actual extinction
  #need to think about effect of minsd?
  #maybe to think about adding 2 obs params again?... Or reducing?
  tmp<-(exp(mean(smp_EDM[,1]))*tmpsmp_proc[17,])
  hist(tmp[tmp<0.1], breaks = 20)
  hist(pmax(tmp[tmp<0.1], 0.01), breaks = 20)
  
  x<-0.1999774
  pnorm(x, sqrt(exp(mean(smp_EDM[,2]))*x^exp(mean(smp_EDM[,3]))))
  
  std_tmp<-exp(mean(smp_EDM[,1]))*x
  dnorm((0-x)/std_tmp,log=TRUE)-log(std_tmp)
  pnorm(-0/(0.01))
  
  
  
  apply(tmpsmp_proc, 1, function(x) length(unique(x)))
  
  std_tmp<-sqrt(exp(rep(smp_EDM[,2], each=nrow(tmpsmp)))*tmpsmp^exp(
    rep(smp_EDM[,3], each=nrow(tmpsmp))
  ))
  
  tm<-dat$time[yps][libuse_y[i,1]:libuse_y[i,2]]
  
  if(FALSE) {
    #full sets of particles
    pr_died_cv<-pnorm(0, tmpsmp, std_tmp)
    pr_died_cv[tmpsmp==0 | std_tmp==0]<-NA
    
    pr_died<-tmpsmp_proc[-1,]==0 & tmpsmp_proc[-nrow(tmpsmp_proc),]>0 & is.finite(tmpsmp_proc[-1,]) & is.finite(tmpsmp_proc[-nrow(tmpsmp_proc),])
    pr_died[!(tmpsmp_proc[-nrow(tmpsmp_proc),]>0 & is.finite(tmpsmp_proc[-1,]) & is.finite(tmpsmp_proc[-nrow(tmpsmp_proc),]))]<-NA
    
    
    cum_pr_died_cv<-1-apply(pr_died_cv, 2, function(x) {x[is.na(x)]<-0; cumprod(1-x)})
    cum_pr_died<-1-apply(pr_died, 2, function(x) {x[is.na(x)]<-0; cumprod(1-x)})
    
    #BINOMIAL MODEL!!
    plot(cum_pr_died_cv[-nrow(cum_pr_died_cv),], as.numeric(pr_died))
    
    plot(rowMeans(cum_pr_died_cv[-nrow(cum_pr_died_cv),]), rowMeans(cum_pr_died))
  }
  
  if(FALSE) {
    #plot states
    tmp<-t(apply(tmpsmp, 1, function(x) quantile(x, probs = pnorm(c(-2:2)), na.rm=T)))
    plot(tm, tmp[,3])
    
    tmp_proc<-t(apply(tmpsmp_proc, 1, function(x) quantile(x, probs = pnorm(c(-2:2)), na.rm=T)))
    plot(tm, tmp_proc[,3])
  }
  if(FALSE) {
    #plot CV
    plot(tmpsmp, std_tmp)
    abline(a=0, b=1, lty=2, col=2, lwd=2)
    
    plot(tmpsmp, std_tmp/tmpsmp)
  }
  
  #probabilities
  pr_died<-rowSums(tmpsmp_proc[-1,]==0 & tmpsmp_proc[-nrow(tmpsmp_proc),]>0 & is.finite(tmpsmp_proc[-1,]) & is.finite(tmpsmp_proc[-nrow(tmpsmp_proc),]))/
    rowSums(tmpsmp_proc[-nrow(tmpsmp_proc),]>0 & is.finite(tmpsmp_proc[-1,]) & is.finite(tmpsmp_proc[-nrow(tmpsmp_proc),]))
  
  pr_died_cv<-pnorm(0, tmpsmp, std_tmp)
  pr_died_cv[tmpsmp==0 | std_tmp==0]<-NA
  pr_died_cv_mean<-apply(pr_died_cv, 1, function(x) mean(x[is.finite(x)]))
  if(FALSE) {
    plot(pr_died_cv_mean[-length(pr_died_cv_mean)], pr_died)
    
    plot(cumprod((1-pr_died_cv_mean[-length(pr_died_cv_mean)])[-c(1:3)]),
         cumprod((1-pr_died)[-c(1:3)]))
  }
  points(1-cumprod((1-pr_died_cv_mean[-length(pr_died_cv_mean)])[-c(1:3)]),
         1-cumprod((1-pr_died)[-c(1:3)]),
         col=adjustcolor(collst[i], alpha.f = 0.8))
  #cv_int<-t(apply(std_tmp/tmpsmp, 1, function(x) quantile(x, probs = pnorm(c(-2:2)), na.rm=T)))
  #polygon(c(tm, rev(tm)), c(cv_int[,2], rev(cv_int[,4])), col=adjustcolor(collst[i], alpha.f = 0.8))
  #lines(c(tm), c(cv_int[,3]), col=adjustcolor(collst[i], alpha.f = 0.8), lwd=1.5)
  #i<-i+1
}
abline(h=c(0, 1), lty=3)
