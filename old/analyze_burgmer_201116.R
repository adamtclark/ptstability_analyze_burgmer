#!/usr/bin/env Rscript

if(FALSE) {
  error
  rm(list=ls())
  setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/analyze_burgmer/")
}
  
## set-up
nclus <- 32
nrep <- 5

# load packages
if(length(dir("/cl_tmp/clarka/Rpkg/"))>0) {
  .libPaths("/cl_tmp/clarka/Rpkg/")
}
require(rEDM)
require(BayesianTools)
require(viridis)

source("../pts_r_package/pttstability/R/bayesfun.R")
source("../pts_r_package/pttstability/R/fake_data.R")
source("../pts_r_package/pttstability/R/logit_funs.R")
source("../pts_r_package/pttstability/R/particlefilter.R")

# load treatments
trtsplst<-read.csv("data/trtmat.csv", stringsAsFactors = FALSE)

commArgin<-commandArgs(1)
if(length(commArgin)==0) {
  commArgin<-round(runif(1)*1e6)
  commArg_ps<-69
} else {
  commArg_ps<-as.numeric(commArgin)
}
print(commArg_ps)

commArgin_kp = commArg_ps

for(iclu in 1:nrep) {
  commArgin = commArgin_kp + (iclu - 1)*nclus
  
  try({
  if(commArgin <= 152) {
    commArg_ps = commArgin
    
    simname<-paste(c(trtsplst[commArg_ps,]), collapse = "_")
    
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
    libmat<-NULL
    datnum<-1:nrow(dat)
    doplot<-FALSE
    if(doplot) {
      par(mfrow=c(4,2), mar=c(4,4,2,2))
    }
    
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
    trtuse<-as.character(trtsplst[commArg_ps,2])
    libuse<-as.matrix(libmat[libmat$trt==trtuse,c("start", "end")])
    yps<-which(dat$treatment==trtuse)
    y<-dat[,as.character(trtsplst[commArg_ps,1])][yps]
    libuse_y<-libuse-min(libuse)+1
    y<-y/sd(y)
    
    #par(mfrow=c(1,1))
    #plot(dat$time[yps], y, xlab="time", ylab="abundance", type="p")
    
    # get EDM parameters
    sout<-NULL
    for(E in 2:4) {
      sout<-rbind(sout, s_map(y, E=E, silent = TRUE, lib = libuse_y))
    }
    tuse<-sout$theta[which.max(sout$rho)]
    euse<-sout$E[which.max(sout$rho)]
    
    spred<-s_map(y, E=euse, theta=tuse, silent = TRUE, lib = libuse_y, stats_only = FALSE, save_smap_coefficients = TRUE)
    #plot(spred$model_output[[1]]$Obs, pmax(0, spred$model_output[[1]]$Predictions), xlab="Obs", ylab="Pred"); abline(a=0, b=1, lty=2); abline(h=0, v=0, lty=3)
    edm_r2<-cor(spred$model_output[[1]]$Obs, pmax(0, spred$model_output[[1]]$Predictions), use="complete")^2
    edm_mae<-mean(abs(spred$model_output[[1]]$Obs-pmax(0, spred$model_output[[1]]$Predictions)),na.rm=T)
    edm_mte<-mean(abs(spred$model_output[[1]]$Obs-mean(spred$model_output[[1]]$Obs, na.rm=T)), na.rm=T)
    
    
    # set priors
    minvUSE_edm<-c(log(0.01), log(0.01), log(0.01))
    maxvUSE_edm<-c(log(2), log(2), log(3))
    
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
    smap_coefs <- process_scof(smap_coefs)
    
    #likelihood and bayesian set-ups for EDM functions
    likelihood_EDM_piecewise<-function(x) {
      xuse<-x
      tuse_edm<-tuse
      
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
      tuse_edm<-tuse
      
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
    out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=200))
    #plot(out_EDM, start = floor(niter/5))
    #correlationPlot(out_EDM, start = floor(niter/5))
    
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
    
    # calculate pmor per timestep
    xtedm<-simout_smp_noproc
    stdedm<-sqrt(rep(exp(smp_EDM[,2]), each=nrow(xtedm))*
                   xtedm^rep(exp(smp_EDM[,3]), each=nrow(xtedm)))
    
    pmedm_analy<-matrix(nrow=ncol(xtedm), ncol=nrow(libuse_y)+1)
    
    for(i in 1:ncol(xtedm)) {
      for(j in 1:nrow(libuse_y)) {
        ps<-which(!is.na(xtedm[,i]) & xtedm[,i]>0)
        ps<-ps[ps%in%c(libuse_y[j,1]:libuse_y[j,2])]
        if(length(ps)>0) {
          pmedm_analy[i,j]<-sum(pnorm(0, xtedm[ps,i], stdedm[ps,i]))/length(ps)
        } else {
          pmedm_analy[i,j]<-0
        }
      }
      
      ps<-(!is.na(xtedm[,i]) & xtedm[,i]>0)
      pmedm_analy[i,nrow(libuse_y)+1]<-sum(pnorm(0, xtedm[ps,i], stdedm[ps,i]))/sum(ps)
    }
    
    
    save(list = c("out_EDM", "smp_EDM", "simout", "sdout", "simout_smp", "simout_smp_noproc", "pmedm_analy"),
         file = paste("datout/", simname, ".rda", sep=""), version = 2)
  }
  })
}
