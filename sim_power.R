#!/usr/bin/Rscript
load("env_power.RData")
args<-commandArgs(TRUE)
name=args[1]
K=2/1000
nsim=100000
N=sum(!is.na(ISISenvgen[,name]))
evi<-function(x){
  return(log(x/(1-x)))
}

eff=min(c(result_cc[name,"size_m"],result_cc[name,"size_pr"]))
enveff<-c()
env<-c()
lev<-length(table(ISISenvgen[,name]))
for(k in 1:lev){
  alpha=eff^((k-1)/(lev-1))
  enveff<-c(enveff, rep(alpha,floor((1+K*(alpha-1))/alpha*table(ISISenvgen[,name])[[k]])))
  env<-c(env,rep((k-1)/(lev-1),floor((1+K*(alpha-1))/alpha*table(ISISenvgen[,name])[[k]])))
}

psim<-numeric(nsim)
for(sim in 1:nsim){
  set.seed(sim)
  envindex<-sample.int(length(env),floor(N/K),replace=T)
  riskenv<-log(enveff[envindex])
  riskgen<-evi(c(sample(ISISWT[ISISWT[,1]==0,2],floor(N/K)-N,replace=T),sample(ISISWT[ISISWT[,1]==1,2],N,replace=T)))
  risk=riskenv+riskgen
  prob<-evi(runif(floor(N/K)))
  T1D=which(risk-prob>=sort(risk-prob,decreasing = T)[N])
  psim[sim]<-summary(lm(riskgen[T1D]~(env[envindex[T1D]])))$coefficients[2,4]
  if (summary(lm(riskgen[T1D]~(env[envindex[T1D]])))$coefficients[2,1]>0){
    psim[sim]<-psim[sim]/2
  } else {
    psim[sim]<-1-psim[sim]/2
  }
}

result<-paste(N,name,sum(psim<0.05)/nsim)
system(paste("echo ",result," >>result_power"))
