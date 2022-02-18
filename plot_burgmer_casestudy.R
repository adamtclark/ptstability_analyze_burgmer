error
rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/analyze_burgmer/")

## set-up
# load packages
require(rEDM)
require(BayesianTools)
require(viridis)
require(pttstability)

# load treatments
trtsplst<-read.csv("data/trtmat.csv", stringsAsFactors = FALSE)

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
spps<-c(4:22)
collst<-viridis(length(spps))

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

#likelihood and bayesian set-ups for EDM functions
likelihood_EDM_piecewise<-function(x, lowerbound = -999, maxNuse = 512000) {
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

dlst<-c("Chlamydomonas.terricola_HSP_long.rda", "Chlamydomonas.terricola_LSP_long.rda")

#set up plotting
ymax<-5
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

pdf("plotout/burgmer_examples_220218.pdf", width=6, height=8, colormodel = "cmyk", useDingbats = FALSE)
nm= 2
psm = c(rep(1,nm), 2, rep(3:4, each=nm))
m = cbind(c(1,7,2,4)[psm], c(3,7,6,5)[psm])
layout(m)
par(oma=c(1,1.5,0,0))

for(ii in 1:length(dlst)) {
  par(mar=c(3.5,3.5,2,1))
  if(ii==1) {
    dcolps<-6
  } else {
    dcolps<-3
  }
  
  load(paste("datout/", dlst[ii], sep=""))
  #summary(out_EDM)
  #gelmanDiagnostics(out_EDM, start = floor(niter/5))
  #plot(out_EDM, start = floor(niter/5))
  #correlationPlot(out_EDM, start = floor(niter/5))
  
  simname<-gsub(".rda", "", dlst[ii])
  tmpsn<-strsplit(simname, "_", fixed=T)[[1]]
  
  commArg_ps<-which((trtsplst[,1]==tmpsn[1]) & (trtsplst[,2]==tmpsn[2]))
  
  #set-up dat
  trtuse<-as.character(trtsplst[commArg_ps,2])
  libuse<-as.matrix(libmat[libmat$trt==trtuse,c("start", "end")])
  yps<-which(dat$treatment==trtuse)
  y<-dat[,as.character(trtsplst[commArg_ps,1])][yps]
  libuse_y<-libuse-min(libuse)+1
  meany = mean(y)
  y<-y/meany
  # convert from um^3/ml to mm^3/dL
  if(ii == 1) {
    conversion_factor = meany*(((1e-6)^3)/((1e-3)^3))*(1e3/1)*(1/10)
  } else {
    conversion_factor = c(conversion_factor, meany*(((1e-6)^3)/((1e-3)^3))*(1e3/1)*(1/10))
  }
  
  if(ii==1) {
    collst<-c(viridis(nrow(libuse_y)), "black")
  } else {
    collst<-c(magma(nrow(libuse_y)), "black")
    mcol2 = c(viridis(nrow(libuse_y)), "black")[6]
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
  edm_r2<-cor(spred$model_output[[1]]$Obs, pmax(0, spred$model_output[[1]]$Predictions), use="complete")^2
  edm_mae<-mean(abs(spred$model_output[[1]]$Obs-pmax(0, spred$model_output[[1]]$Predictions)),na.rm=T)
  edm_mte<-mean(abs(spred$model_output[[1]]$Obs-mean(spred$model_output[[1]]$Obs, na.rm=T)), na.rm=T)
  
  if(ii==1) {
    luse<-"a."
  } else {
    luse<-"b."
    axadj<-0.06
  }
  #par(mar=c(3.5,1.2,2,1.2))
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

    polygon(c(tm, rev(tm)), c(tmp[,2], rev(tmp[,4]))*conversion_factor[ii], col=adjustcolor(collst[i], alpha.f = 0.5))
  }
  abline(h=0, lty=3)
  mtext("Time (days)", side = 1, line = 2.6)
  if(ii==1) {
    mtext(expression(paste("Est. Abund., mm"^3, "dL"^-1)), side = 2, line = 2.6)
  }
  if(ii==1) {
    cbxfun()
  }
  box()
  title(luse, line = axln, cex.main=axcx, adj=axadj)
  axadj<-0.02
  
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
  
  #parameters
  if(ii==1) {
    tmp1 = smp_EDM[,1]
    tmp2 = minvUSE_edm[1]
    tmp3 = maxvUSE_edm[1]
  }
  
  if(ii==2) {
    #par(mar=c(3.5,3.5,2,1.2))
    lfuse<-exp
    plot(density(lfuse(smp_EDM[,1]), bw = diff(range(lfuse(smp_EDM[,1])))/20,
                 from=lfuse(minvUSE_edm[1]), to=lfuse(maxvUSE_edm[1])),
         main="", xlab="", ylab="", lwd=1.5, axes=F, lty=1, col=collst[dcolps],
         ylim=c(0,5.6))
    abline(v=mean(lfuse(smp_EDM[,1])), lwd=2, lty=2, col=collst[dcolps])
    
    tmpd = density(lfuse(tmp1), bw = diff(range(lfuse(tmp1)))/20,
                   from=lfuse(tmp2), to=lfuse(tmp3))
    lines(tmpd$x, tmpd$y,
         lwd=1.5, lty=1, col=mcol2)
    abline(v=mean(lfuse(tmp1)), lwd=2, lty=2, col=mcol2)
    
    axis(2)
    axis(1)
    box()
    abline(h=0, lty=3)
    abline(v=c(lfuse(minvUSE_edm[1]), lfuse(maxvUSE_edm[1])), lty=3)
    mtext(expression(paste("Obs. Error, ", sigma[italic(O)])), side = 1, line = 2.6)
    mtext("Density", side = 2, line = 2.6)
    cbxfun()
    box()
    title("f.", line = axln, cex.main=axcx, adj=axadj)
  }
  
  
  if(ii==1) {
    tmp11 = smp_EDM[,2]
    tmp22 = minvUSE_edm[2]
    tmp33 = maxvUSE_edm[2]
  }
  
  if(ii==2) {
    plot(density(lfuse(smp_EDM[,2])*conversion_factor[ii], bw = diff(range(lfuse(smp_EDM[,2])*conversion_factor[ii]))/20,
                 from=lfuse(minvUSE_edm[2])*conversion_factor[ii], to=lfuse(maxvUSE_edm[2])*conversion_factor[ii]),
         main="", xlab="", ylab="", lwd=1.5, axes=F, lty=1, col=collst[dcolps],
         ylim=c(0,10))
    abline(v=mean(lfuse(smp_EDM[,2])*conversion_factor[ii]), lwd=2, lty=2, col=collst[dcolps])
    
    
    tmpd2 = density(lfuse(tmp11)*conversion_factor[2], bw = diff(range(lfuse(tmp11)*conversion_factor[2]))/20,
                   from=lfuse(tmp22)*conversion_factor[2], to=lfuse(tmp33)*conversion_factor[2])
    lines(tmpd2$x, tmpd2$y,
        lwd=1.5, lty=1, col=mcol2)
    abline(v=mean(lfuse(tmp11)*conversion_factor[2]), lwd=2, lty=2, col=mcol2)
    
    axis(2)
    axis(1)
    box()
    abline(h=0, lty=3)
    abline(v=c(lfuse(minvUSE_edm[2]), lfuse(maxvUSE_edm[2]))*conversion_factor[1], lty=3)
    mtext(expression(paste("Proc. Noise, ", sigma[italic(P)])), side = 1, line = 2.6)
    cbxfun()
    box()
    title("g.", line = axln, cex.main=axcx, adj=axadj)
  }
  
  #time to extinction
  #par(mar=c(3.5,1.2,2,1.2))
  if(ii==1) {
    luse<-"d."
  } else {
    luse<-"e."
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
  
  xmx<-1000#pmin(365*20, quantile(text_extrap, 0.99))
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
  
  tmpsim = rowMeans(simout)
  tmpsim = tmpsim[tmpsim>1e-4]
  bardattmp = cbind(c(mean(mean(exp(smp_EDM[,1]))*tmpsim)^2,
                   (mean(exp(smp_EDM[,2])))^2,
                   NA))
  bardattmp[3] = var(y,na.rm=TRUE)-sum(bardattmp[1:2])
  
  if(ii == 1) {
    bardat = NULL
  }
  bardat = cbind(bardat, bardattmp)
  
  if(ii == 2) {
    par(mar=c(3,5.5,1,13.5))
    acol = function(col, a=c(0.9, 0.4, 0)) {
      c(adjustcolor(col, a[1]),
        adjustcolor(col, a[2]),
        adjustcolor(col, a[3]))
    }
    
    tmp = barplot(bardat[,2:1], axes=F, names.arg = rep("", ncol(bardat)),
                  col = cbind(acol("black")), horiz = TRUE)
    axis(1)
    axis(2, at = tmp, labels = c("Low Temp.", "High Temp."), cex.axis=1.2, las=2)
    title("c.", line = axln+1.3, cex.main=axcx, adj=axadj)
    mtext("Standardised Temporal Variance", 1, line = 2.4)
    #box()
    legend(2.1, 2.5, legend = c("Observation Error", "Process Noise", "Deterministic Variation"),
           fill = acol("black"),
           border = 1,
           bty="n", xpd = NA, cex=1.2)
  }
}  

dev.off()
