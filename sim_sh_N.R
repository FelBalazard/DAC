#!/usr/bin/Rscript
load("env_power.RData")
args<-commandArgs(TRUE)
library(AUC)
K=2/1000
eff=3
shift=print(as.numeric(args[1]))
N=as.numeric(args[2])
enveff<-c(rep(eff,floor(1000/(1+eff))),rep(1,1000-floor(1000/(1+eff))))
nsim=100000
psim<-numeric(nsim)
evi<-function(x){
  return(log(x/(1-x)))
}
genmod<-ISISWT
genmod[,2]<-evi(ISISWT[,2])
genmod[,2]<-genmod[,2]+shift*genmod[,1]
calibration.model<-glm(genmod[,1]~genmod[,2],family="binomial")
genmod[,2]<-evi(calibration.model$fitted.values)
aucmod=auc(roc(genmod[,2],factor(genmod[,1])))
for(sim in 1:nsim){
  set.seed(sim)
  riskenv<-sample(log(enveff),floor(N/K),replace=T)
  riskgen<-c(sample(genmod[genmod[,1]==0,2],floor(N/K)-N,replace=T),sample(genmod[genmod[,1]==1,2],N,replace=T))
  risk=riskenv+riskgen
  prob<-evi(runif(floor(N/K)))
  T1D=which(risk-prob>=sort(risk-prob,decreasing = T)[N])
  if (lm(riskgen[T1D]~riskenv[T1D])$coefficients[2]<0){
    psim[sim]<-summary(lm(riskgen[T1D]~(riskenv[T1D])))$coefficients[2,4]/2
  } else {
    psim[sim]<-1-summary(lm(riskgen[T1D]~(riskenv[T1D])))$coefficients[2,4]/2
  }
}
result<-paste(K,N,aucmod,sum(psim<0.05)/nsim)
system(paste("echo ",result," >>result_auc"))
