error
rm(list=ls())
setwd("~/Dropbox/Projects/041_Powerscaling_stability/src/analyze_burgmer/")

#identify successful runs
trtsplst<-read.csv("data/trtmat.csv", stringsAsFactors = FALSE)
hpcres<-dir("datout/")
hpcres<-hpcres[grep(".csv", hpcres)]
hpcres<-unique(gsub("_summary", "", gsub(".csv", "", gsub(".rda", "", hpcres))))
hpcres<-t(sapply(strsplit(hpcres, "_", fixed=TRUE), cbind))
hpcres<-data.frame(hpcres)
colnames(hpcres)<-c("species", "trt")

trtsplst$success<-paste(trtsplst$species, trtsplst$trt)%in%c(paste(hpcres$species, hpcres$trt))

#note cases of zero abundance
trtsplst$fzero<-FALSE
trtsplst$fzero_E4<-FALSE
datfull<-read.table("data/burgmer phyto biovol.csv", header=TRUE, sep=";", stringsAsFactors = FALSE)
datfull<-datfull[rowSums(!is.na(datfull[,-c(1:3)]),na.rm=T)!=0,]

datedat<-data.frame(date=datfull$date, year=NA, month=NA, day=NA)
datedat$year<-as.numeric(substr(datedat$date, 7, 10))
datedat$month<-as.numeric(substr(datedat$date, 4, 5))
datedat$day<-as.numeric(substr(datedat$date, 1, 2))
datedat$time<-suppressWarnings(as.POSIXlt(datedat$date, format = "%d.%m.%Y")$yday)+(as.numeric(datedat$year==2006)*365)
datedat$time<-datedat$time-min(datedat$time)+1

datfull$time<-datedat$time[match(datfull$date, datedat$date)]
tmlst<-sort(unique(datfull$time))

for(i in 1:nrow(trtsplst)) {
  ps<-which(datfull$treatment==trtsplst$trt[i])
  trtsplst$fzero[i]<-mean(datfull[ps,trtsplst$species[i]]<=0,na.rm=T)
  
  ps2<-which(datfull$treatment==trtsplst$trt[i] & (datfull$time>tmlst[4]))
  trtsplst$fzero_E4[i]<-mean(datfull[ps2,trtsplst$species[i]]<=0,na.rm=T)
}

#note - only cases without non-zero elements beyond E4 work
hist(trtsplst$fzero_E4[!trtsplst$success])
hist(trtsplst$fzero_E4[trtsplst$success])

success_table<-table(trtsplst$species[trtsplst$success],trtsplst$trt[trtsplst$success])

#populate simulation summaries
trtsplst$Euse<-NA
trtsplst$tuse<-NA
trtsplst$edm_r2<-NA
trtsplst$edm_mae<-NA
trtsplst$edm_mte<-NA
trtsplst$mormu<-NA
trtsplst$morsd<-NA

for(i in 1:nrow(trtsplst)) {
  if(trtsplst$success[i]) {
    tmp<-read.csv(paste("datout/", trtsplst$species[i], "_", trtsplst$trt[i], "_summary.csv", sep=""), stringsAsFactors = FALSE)
    trtsplst[i,6:12]<-tmp[3:9]
  }
}


#keep only fully successful species
trtlst_use<-trtsplst[trtsplst$species%in%c(row.names(success_table[rowSums(success_table)==8,])),]
trtlst_use$E1<-(1-trtlst_use$edm_mae/trtlst_use$edm_mte)

tmp<-trtlst_use[,c("species", "trt", "E1", "mormu", "morsd")]
tmp[,3:5]<-round(tmp[,3:5],3)


tmp


#In any case, focus on these examples:
#Coccoid.Cyanobacterium LOW as example?
#Chlamydomonas.terricola HI
#Chlamydomonas.terricola probably best (HSP vs. HVP; HSP vs. LSP; LSA vs. LSP)


