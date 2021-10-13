#!/usr/bin/env Rscript

if(FALSE) {
  error
  rm(list=ls())
  setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/analyze_burgmer/")
}
  
# load packages
require(rEDM)
require(BayesianTools)
#require(pttstability)

if(length(dir("/cl_tmp/clarka/Rpkg/"))>0) {
  .libPaths("/cl_tmp/clarka/Rpkg/")
}

source("../pts_r_package/pttstability/R/bayesfun.R")
source("../pts_r_package/pttstability/R/fake_data.R")
source("../pts_r_package/pttstability/R/logit_funs.R")
source("../pts_r_package/pttstability/R/particlefilter.R")

# load treatments
trtsplst<-read.csv("data/trtmat.csv", stringsAsFactors = FALSE)
zero_cutoff = 0

commArgin<-commandArgs(1)
if(length(commArgin)==0) {
  commArgin<-sample(nrow(trtsplst),1)
  commArg_ps<-commArgin
} else {
  commArg_ps<-as.numeric(commArgin)
}

#for(iclu in c(9:16)) {#c(13,15)) {
#  commArg_ps = iclu
  
  simname<-paste(c(trtsplst[commArg_ps,]), collapse = "_")
  
  # load data
  datfull<-read.csv2("data/burgmer phyto biovol_casestudy.csv")
  datfull<-datfull[!is.na(datfull$Chlamydomonas.terricola),]
  
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
  
  # make time series for each treatmetn and species
  dat<-data.frame(treatment=datfull$treatment, number=datfull$number, time=datedat$time, Chlamydomonas.terricola = datfull[,-c(1:3)])
  dat<-dat[order(dat$treatment, dat$number, dat$time),]
  
  # make plots for each
  libmat<-NULL
  datnum<-1:nrow(dat)
  
  for(i in 1:nrow(trtmat)) {
    ps1<-which(dat$treatment==trtmat$trt[i])
    replst<-sort(unique(dat$number[ps1]))
    
   
    for(j in 1:length(replst)) {
      ps2<-which(dat$number[ps1]==replst[j])
      
      libmat<-rbind(libmat, data.frame(trt=trtmat$trt[i], rep=replst[j], start=min(datnum[ps1][ps2]), end=max(datnum[ps1][ps2])))
    }
  }
  
  ## run particle filter
  trtuse<-as.character(trtsplst[commArg_ps,2])
  libuse<-as.matrix(libmat[libmat$trt==trtuse,c("start", "end")])
  yps<-which(dat$treatment==trtuse)
  y<-dat[,as.character(trtsplst[commArg_ps,1])][yps]
  libuse_y<-libuse-min(libuse)+1
  y<-y/mean(y)
  
  
  # get EDM parameters
  sout<-NULL
  for(E in 2:4) {
    sout<-rbind(sout, s_map(y, E=E, silent = TRUE, lib = libuse_y))
  }
  tuse<-sout$theta[which.max(sout$rho)]
  euse<-sout$E[which.max(sout$rho)]
  
  spred<-s_map(y, E=euse, theta=tuse, silent = TRUE, lib = libuse_y, stats_only = FALSE, save_smap_coefficients = TRUE)
  edm_r2<-cor(spred$model_output[[1]]$Obs, pmax(0, spred$model_output[[1]]$Predictions), use="complete")^2
  edm_mae<-mean(abs(spred$model_output[[1]]$Obs-pmax(0, spred$model_output[[1]]$Predictions)),na.rm=T)
  edm_mte<-mean(abs(spred$model_output[[1]]$Obs-mean(spred$model_output[[1]]$Obs, na.rm=T)), na.rm=T)
  
  
  # set priors
  minvUSE_edm<-c(log(0.001), log(0.001))
  maxvUSE_edm<-c(log(2), log(2))
  
  #density, sampler, and prior functions for EDM function
  density_fun_USE_edm<-function(param) density_fun0(param = param, minv = minvUSE_edm, maxv=maxvUSE_edm)
  sampler_fun_USE_edm<-function(x) sampler_fun0(n = 1, minv = minvUSE_edm, maxv=maxvUSE_edm)
  prior_edm <- createPrior(density = density_fun_USE_edm, sampler = sampler_fun_USE_edm,
                           lower = minvUSE_edm, upper = maxvUSE_edm)
  ## Run filter
  niter<-6e3 #number of steps for the MCMC sampler
  N<-1e3 #number of particles
  Euse<-euse #number of embedding dimensions
  
  smap_coefs<-spred$smap_coefficients[[1]]
  smap_coefs <- process_scof(smap_coefs)
  
  #likelihood and bayesian set-ups for EDM functions
  likelihood_EDM_piecewise<-function(x, lowerbound = -999, maxNuse = 512000) {
    # piecewise likelihood function with automatic N selection for each time series chunk
    xuse<-x
    tuse_edm<-tuse
    
    LLtot<-0
    
    for(i in 1:nrow(libuse_y)) {
      if(!is.na(LLtot)) {
        ysegment<-y[libuse_y[i,1]:libuse_y[i,2]]
        smap_coefs_segment<-smap_coefs[libuse_y[i,1]:libuse_y[i,2],]
        
        Nuse = N
        LLtmp = -Inf
        while(LLtmp <=lowerbound & Nuse <= maxNuse) {
          LLtmp = likelihood0(param = xuse, y=ysegment, parseparam = function(x) parseparam0(x, colparam=c(logit(1e-6), log(0.1))),
                              detfun = EDMfun0, edmdat = list(E=Euse, theta=tuse_edm, smp_cf=smap_coefs_segment, ytot=y), N = Nuse, lowerbound = lowerbound)
          Nuse = 2*Nuse
          #print(paste(i, Nuse))
        }
        if(Nuse <= maxNuse) {
          LLtot<-LLtot+LLtmp
        } else {
          LLtot<-NA
        }
      }
    }
    if(is.na(LLtot)) {
      LLtot = lowerbound
    }
    
    return(sum(LLtot))
  }
  
  particleFilterLL_piecewise<-function(param, N, lowerbound = -999, maxNuse = 512000) {
    # piecewise particle filter function with automatic N selection for each time series chunk
    pars<-parseparam0(param, colparam=c(logit(1e-6), log(0.1)))
    tuse_edm<-tuse
    
    pfout<-data.frame(rep=NA, Nest=rep(NA, length(y)), Nsd=NA, Nsmp=NA, Nsmp_noproc=NA)
    for(i in 1:nrow(libuse_y)) {
      ysegment<-y[libuse_y[i,1]:libuse_y[i,2]]
      smap_coefs_segment<-smap_coefs[libuse_y[i,1]:libuse_y[i,2],]
      
      Nuse = N
      LLtmp = -Inf
      while(LLtmp <=lowerbound & Nuse <= maxNuse) {
        tmp<-particleFilterLL(ysegment, pars, N=N, detfun = EDMfun0,
                                                 edmdat = list(E=Euse, theta=tuse_edm, smp_cf=smap_coefs_segment),
                                                 dotraceback = TRUE, fulltraceback = TRUE)
        LLtmp = tmp$LL
        Nuse = 2*Nuse
      }
      
      
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
  out_EDM <- runMCMC(bayesianSetup = bayesianSetup_EDM, settings = list(iterations=niter, consoleUpdates=1))
  #summary(out_EDM)
  #gelmanDiagnostics(out_EDM)
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
  stdedm<-exp(smp_EDM[,2])
  
  pmedm_analy<-matrix(nrow=ncol(xtedm), ncol=nrow(libuse_y)+1)
  
  for(i in 1:ncol(xtedm)) {
    for(j in 1:nrow(libuse_y)) {
      ps<-which(!is.na(xtedm[,i]) & xtedm[,i]>zero_cutoff)
      ps<-ps[ps%in%c(libuse_y[j,1]:libuse_y[j,2])]
      if(length(ps)>0) {
        pmedm_analy[i,j]<-sum(pnorm(0, xtedm[ps,i], stdedm[i]))/length(ps)
      } else {
        pmedm_analy[i,j]<-0
      }
    }
    
    ps<-(!is.na(xtedm[,i]) & xtedm[,i]>zero_cutoff)
    pmedm_analy[i,nrow(libuse_y)+1]<-sum(pnorm(0, xtedm[ps,i], stdedm[i]))/sum(ps)
  }
  
  
  save(list = c("out_EDM", "smp_EDM", "simout", "sdout", "simout_smp", "simout_smp_noproc", "pmedm_analy"),
       file = paste("datout/", simname, ".rda", sep=""), version = 2)
#}
