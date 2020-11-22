error
rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/analyze_burgmer/")

## set-up
# load packages
require(rEDM)
require(BayesianTools)
require(viridis)

source("../pts_r_package/pttstability/R/bayesfun.R")
source("../pts_r_package/pttstability/R/fake_data.R")
source("../pts_r_package/pttstability/R/logit_funs.R")
source("../pts_r_package/pttstability/R/particlefilter.R")

# load treatments
trtsplst<-read.csv("data/trtmat.csv", stringsAsFactors = FALSE)

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

for(i in 1:nrow(trtmat)) {
  ps1<-which(dat$treatment==trtmat$trt[i])
  replst<-sort(unique(dat$number[ps1]))
  
  for(j in 1:length(replst)) {
    ps2<-which(dat$number[ps1]==replst[j])
    libmat<-rbind(libmat, data.frame(trt=trtmat$trt[i], rep=replst[j], start=min(datnum[ps1][ps2]), end=max(datnum[ps1][ps2])))
  }
}

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

dlst<-c("Chlamydomonas.terricola_LSA.rda", "Chlamydomonas.terricola_LSP.rda")

#set up plotting
ymax<-5.4
axcx<-1.6; axln<-(-1.2); axadj<-0.02
cbxfun<-function(dw1=0.05, dw2=0.17) {
  xp<-par("usr")[1:2]
  yp<-par("usr")[3:4]
  dxs<-(xp[2]-xp[1])*dw1
  dyx<-(yp[2]-yp[1])*dw2
  
  polygon(c(xp[1]+dxs, xp[1]+dxs, xp[1]-dxs, xp[1]-dxs),
          c(yp[2]-dyx, yp[2]+dyx, yp[2]+dyx, yp[2]-dyx), col="white", border=NA)
}

#plot outputs
fullerror<-FALSE
tracedparticles<-TRUE
dev.off()

pdf("plotout/burgmer_examples_201116.pdf", width=8, height=6, colormodel = "cmyk", useDingbats = FALSE)
par(mar=c(3.5,3.5,2,1), oma=c(1,3,0,0), mfrow=c(3,3))

for(ii in 1:length(dlst)) {
  if(ii==1) {
    dcolps<-6
  } else {
    dcolps<-3
  }
  
  load(paste("datout/", dlst[ii], sep=""))
  simname<-gsub(".rda", "", dlst[ii])
  tmpsn<-strsplit(simname, "_", fixed=T)[[1]]
  
  commArg_ps<-which((trtsplst[,1]==tmpsn[1]) & (trtsplst[,2]==tmpsn[2]))
  
  #set-up dat
  trtuse<-as.character(trtsplst[commArg_ps,2])
  libuse<-as.matrix(libmat[libmat$trt==trtuse,c("start", "end")])
  yps<-which(dat$treatment==trtuse)
  y<-dat[,as.character(trtsplst[commArg_ps,1])][yps]
  libuse_y<-libuse-min(libuse)+1
  y<-y/sd(y)
  
  if(ii==1) {
    collst<-c(viridis(nrow(libuse_y)), "black")
  } else {
    collst<-c(magma(nrow(libuse_y)), "black")
  }
  
  
  # get EDM parameters
  sout = NULL
  for(E in 2:4) {
    sout<-rbind(sout, s_map(y, E=E, silent = TRUE, lib = libuse_y))
  }
  tuse<-sout$theta[which.max(sout$rho)]
  euse<-sout$E[which.max(sout$rho)]
  Euse<-euse
  
  spred<-s_map(y, E=euse, theta=tuse, silent = TRUE, lib = libuse_y, stats_only = FALSE, save_smap_coefficients = TRUE)
  edm_r2<-cor(spred$model_output[[1]]$obs, pmax(0, spred$model_output[[1]]$pred), use="complete")^2
  edm_mae<-mean(abs(spred$model_output[[1]]$obs-pmax(0, spred$model_output[[1]]$pred)),na.rm=T)
  edm_mte<-mean(abs(spred$model_output[[1]]$obs-mean(spred$model_output[[1]]$obs, na.rm=T)), na.rm=T)
  
  if(ii==1) {
    par(mfg=c(1,1,3,3))
    luse<-"a."
  } else {
    par(mfg=c(1,2,3,3))
    luse<-"d."
  }
  par(mar=c(3.5,1.2,2,1.2))
  plot(c((dat$time[yps][libuse_y[1,1]:libuse_y[1,2]])[Euse], max(dat$time[yps])), c(0, ymax), xlab="", ylab="", type="n", xaxs="i")
  
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
    #lines(tm, y[libuse_y[i,1]:libuse_y[i,2]], col=collst[i], lwd=1.5)
    polygon(c(tm, rev(tm)), c(tmp[,2], rev(tmp[,4])), col=adjustcolor(collst[i], alpha.f = 0.5))
  }
  abline(h=0, lty=3)
  mtext("Time (days)", side = 1, line = 2.6)
  if(ii==1) {
    mtext("Estimated Abundance", side = 2, line = 2.6)
  }
  cbxfun()
  box()
  title(luse, line = axln, cex.main=axcx, adj=axadj)
  
  #extinction
  tmdat_ci<-array(dim=c(median(libuse_y[,2]-libuse_y[,1]+1),
                        nrow(libuse_y)+1, 3))
  
  for(i in 1:nrow(libuse_y)) {
    pstmp<-libuse_y[i,1]:libuse_y[i,2]
    tmdat_ci[(nrow(tmdat_ci)-length(pstmp)+1):nrow(tmdat_ci),i,2]<-rowMeans(simout_smp[pstmp,]==0, na.rm=T)
  }
  tmdat_ci[,dim(tmdat_ci)[2],2]<-rowMeans(tmdat_ci[,1:(dim(tmdat_ci)[2]-1),2], na.rm=T)
  
  nrep<-matrix(nrow=dim(tmdat_ci)[1], ncol=dim(tmdat_ci)[2],
               data=c(rep(1, dim(tmdat_ci)[1]*(dim(tmdat_ci)[2]-1)),
                      rep((dim(tmdat_ci)[2]-1), dim(tmdat_ci)[1])))
  tmdat_ci[,,1]<-pmax(0, tmdat_ci[,,2]-sqrt(tmdat_ci[,,2]*(1-tmdat_ci[,,2])/nrep))
  tmdat_ci[,,3]<-pmin(1, tmdat_ci[,,2]+sqrt(tmdat_ci[,,2]*(1-tmdat_ci[,,2])/nrep))
  
  if(ii==1) {
    par(mfg=c(2,1,3,3))
    luse<-"b."
  } else {
    par(mfg=c(2,2,3,3))
    luse<-"e."
  }
  
  plot(c((dat$time[yps][libuse_y[1,1]:libuse_y[1,2]])[Euse], max(dat$time[yps])), c(0, 1), xlab="", ylab="", type="n", xaxs="i")
  for(i in 1:(nrow(libuse_y)+1)) {
    if(i<=nrow(libuse_y)) {
      tm<-dat$time[yps][libuse_y[i,1]:libuse_y[i,2]]
    } else {
      tm<-dat$time[yps][libuse_y[i-1,1]:libuse_y[i-1,2]]
    }
    if(length(tm)<nrow(tmdat_ci)) {
      tm<-c(rep(0, nrow(tmdat_ci)-length(tm)), tm)
    }
    
    if(i<=nrow(libuse_y)) {
      lines(tm, tmdat_ci[,i,2], col=collst[i], lwd=1.5)
    } else {
      polygon(c(tm, rev(tm)), c(tmdat_ci[,i,1], rev(tmdat_ci[,i,3])), col=adjustcolor(collst[i], alpha.f = 0.5))
    }
  }
  abline(h=c(0,1), lty=3)
  mtext("Time (days)", side = 1, line = 2.6)
  if(ii==1) {
    mtext(expression(paste("Estimated Mortality, ", Pr[mor])), side = 2, line = 2.6)
  }
  cbxfun()
  box()
  title(luse, line = axln, cex.main=axcx, adj=axadj)
  
  
  #parameters
  par(mar=c(3.5,3.5,2,1.2))
  lfuse<-exp
  par(mfg=c(1,3,3,3))
  plot(density(lfuse(smp_EDM[,1]), bw = diff(range(lfuse(smp_EDM[,1])))/20,
               from=lfuse(minvUSE_edm[1]), to=lfuse(maxvUSE_edm[1])),
       main="", xlab="", ylab="", lwd=1.5, axes=F, lty=1, col=collst[dcolps])
  abline(v=mean(lfuse(smp_EDM[,1])), lwd=2, lty=2, col=collst[dcolps])
  if(ii==1) {
    axis(2)
    axis(1)
    box()
    abline(h=0, lty=3)
    abline(v=c(lfuse(minvUSE_edm[1]), lfuse(maxvUSE_edm[1])), lty=3)
    mtext(expression(paste("Obs. Error Slope, ", beta[obs])), side = 1, line = 2.6)
    mtext("Density", side = 2, line = 2.6)
    cbxfun()
    box()
    title("g.", line = axln, cex.main=axcx, adj=axadj)
  }
  
  par(mfg=c(2,3,3,3))
  plot(density(lfuse(smp_EDM[,2]), bw = diff(range(lfuse(smp_EDM[,2])))/20,
               from=lfuse(minvUSE_edm[2]), to=lfuse(maxvUSE_edm[2])),
       main="", xlab="", ylab="", lwd=1.5, axes=F, lty=1, col=collst[dcolps])
  abline(v=mean(lfuse(smp_EDM[,2])), lwd=2, lty=2, col=collst[dcolps])
  if(ii==1) {
    axis(2)
    axis(1)
    box()
    abline(h=0, lty=3)
    abline(v=c(lfuse(minvUSE_edm[2]), lfuse(maxvUSE_edm[2])), lty=3)
    mtext(expression(paste("Proc. Noise Intercept, ", beta[proc[0]])), side = 1, line = 2.6)
    mtext("Density", side = 2, line = 2.6)
    cbxfun()
    box()
    title("h.", line = axln, cex.main=axcx, adj=axadj)
  }
  
  par(mfg=c(3,3,3,3))
  plot(density(lfuse(smp_EDM[,3]), bw = diff(range(lfuse(smp_EDM[,3])))/20,
               from=lfuse(minvUSE_edm[3]), to=lfuse(maxvUSE_edm[3])),
       main="", xlab="", ylab="", lwd=1.5, axes=F, lty=1, col=collst[dcolps])
  abline(v=mean(lfuse(smp_EDM[,3])), lwd=2, lty=2, col=collst[dcolps])
  if(ii==1) {
    axis(2)
    axis(1)
    box()
    abline(h=0, lty=3)
    abline(v=c(lfuse(minvUSE_edm[3]), lfuse(maxvUSE_edm[3])), lty=3)
    mtext(expression(paste("Proc. Noise Slope, ", beta[proc[1]])), side = 1, line = 2.6)
    mtext("Density", side = 2, line = 2.6)
    cbxfun()
    box()
    title("i.", line = axln, cex.main=axcx, adj=axadj)
  }
  
  
  #time to extinction
  par(mar=c(3.5,1.2,2,1.2))
  if(ii==1) {
    par(mfg=c(3,1,3,3))
    luse<-"c."
  } else {
    par(mfg=c(3,2,3,3))
    luse<-"f."
  }
  
  stepsize<-mean(diff(dat$time[yps][libuse_y[1,1]:libuse_y[1,2]]))
  
  text_extrap<-(1/pmedm_analy)*stepsize
  
  #plot outputs
  dnsrng<-c(0, max(text_extrap))
  dns<-NULL
  
  for(i in 1:ncol(text_extrap)) {
    dns[[i]]<-density(text_extrap[,i], from = dnsrng[1])
  }
  
  yrng<-c(0, max(sapply(dns, function(x) max(x$y))))
  
  xmx<-2000#pmin(365*20, quantile(text_extrap, 0.99))
  plot(c(dnsrng[1], xmx), yrng, type="n", xlab="", ylab="")
  if(ii==1) {
    collst2<-c(adjustcolor(viridis(ncol(text_extrap)-1), alpha.f = 0.6), "black")
  } else {
    collst2<-c(adjustcolor(magma(ncol(text_extrap)-1), alpha.f = 0.6), "black")
  }
  lwdlst<-c(rep(1, ncol(text_extrap)-1), 2)
  ltylst<-c(rep(1, ncol(text_extrap)-1), 1)
  for(i in 1:ncol(text_extrap)) {
    lines(dns[[i]]$x, dns[[i]]$y, col=collst2[i], lwd=lwdlst[i], lty=ltylst[i])
  }
  abline(h=0, v=c(0,1), lty=3)
  mormu<-mean(text_extrap[,ncol(text_extrap)],na.rm=T)
  abline(v=mormu, lwd=2, lty=2)
  
  mtext("Time to Extinction (days)", side = 1, line = 2.6)
  if(ii==1) {
    mtext("Density", side = 2, line = 2.6)
  }
  
  pmedm_tot<-text_extrap
  
  cbxfun()
  box()
  title(luse, line = axln, cex.main=axcx, adj=axadj)
}  

dev.off()
