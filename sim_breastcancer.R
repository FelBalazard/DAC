evi<-function(x){
  return(log(x/(1-x)))
}
print(12)
#alcohol brca1
a=chisq.test(as.table(matrix(c(74,173,190,656),2,2)))
K=0.12
RR=1.32
eff=eff=RR*K/(1-RR*K)/K*(1-K)
BRCAeff=0.6/(1-0.6)/K*(1-K)
N=53000
N1=1093
enveff<-c(rep(eff,floor(846/RR)),rep(1,1093-846))
geneff<-c(rep(BRCAeff,1),rep(1,999))
nsim=100000
CPR<-numeric(nsim)
Nbrca<-numeric(nsim)
Nalc<-numeric(nsim)

for(sim in 1:nsim){
  set.seed(sim)
  riskenv<-sample(log(enveff),floor(N/K),replace=T)
  riskgen<-sample(log(geneff),floor(N/K),replace=T)
  risk=riskenv+riskgen
  prob<-evi(runif(floor(N/K)))
  BreastCancer=(risk-prob>=sort(risk-prob,decreasing = T)[N])
  StudySelect=union(which(riskgen>0 & BreastCancer),sample(which(riskgen==0 & BreastCancer),N1-sum(riskgen>0 & BreastCancer)))
  t=table(riskenv[StudySelect],riskgen[StudySelect])
  CPR[sim]=t[1,1]*t[2,2]/t[1,2]/t[2,1]
  #if(sim %%50==0){print(table(riskenv[StudySelect],riskgen[StudySelect]))}
  Nbrca[sim]=sum(riskgen[StudySelect]>0)
  Nalc[sim]=sum(riskenv[StudySelect]>0)
}
print(mean(Nalc))
print(190+656)
print(mean(Nbrca))
print(74+190)
alc1=numeric(3)
alc1[1]=median(CPR)
alc1[2]=quantile(CPR,0.025)
alc1[3]=quantile(CPR,0.975)
print(alc1)

#age at menarche and brca2
b=chisq.test(as.table(matrix(c(35,182,66,205),2,2)))
RR=1.05^-4
K=0.12
eff=RR*K/(1-RR*K)/K*(1-K)
BRCAeff=0.6/(1-0.6)/K*(1-K)
N=10000
N1=488
enveff<-c(rep(eff,floor(271/RR)),rep(1,N1-271))
geneff<-c(rep(BRCAeff,2),rep(1,998))
nsim=100000
CPR<-numeric(nsim)
Nbrca<-numeric(nsim)
Nearly<-numeric(nsim)

for(sim in 1:nsim){
  set.seed(sim)
  riskenv<-sample(log(enveff),floor(N/K),replace=T)
  riskgen<-sample(log(geneff),floor(N/K),replace=T)
  risk=riskenv+riskgen
  prob<-evi(runif(floor(N/K)))
  BreastCancer=(risk-prob>=sort(risk-prob,decreasing = T)[N])
  StudySelect=union(which(riskgen>0 & BreastCancer),sample(which(riskgen==0 & BreastCancer),N1-sum(riskgen>0 & BreastCancer)))
  t=table(riskenv[StudySelect],riskgen[StudySelect])
  #if(sim %%50==0){print(table(riskenv[StudySelect],riskgen[StudySelect]))}
  CPR[sim]=1/(t[1,1]*t[2,2]/t[1,2]/t[2,1])
  Nbrca[sim]=sum(riskgen[StudySelect]>0)
  Nearly[sim]=sum(riskenv[StudySelect]==0)
}
print(mean(Nbrca))
print(35+66)
print(mean(Nearly))
print(35+182)
menarc2=numeric(3)
menarc2[1]=median(CPR)
menarc2[2]=quantile(CPR,0.025)
menarc2[3]=quantile(CPR,0.975)
print(menarc2)

#parity brca1
c=chisq.test(as.table(matrix(c(81,281,63,175),2,2)))
K=0.12
R1=K*0.84
R0=K*1.29
eff=R1*(1-R0)/R0/(1-R1)
BRCAeff=0.6/(1-0.6)/K*(1-K)
N=29000
N1=600
enveff<-c(rep(eff,floor(238*R0/R1)),rep(1,362))
geneff<-c(rep(BRCAeff,1),rep(1,999))
nsim=100000
CPR<-numeric(nsim)
Nbrca<-numeric(nsim)
Nulliparous<-numeric(nsim)

for(sim in 1:nsim){
  set.seed(sim)
  riskenv<-sample(log(enveff),floor(N/K),replace=T)
  riskgen<-sample(log(geneff),floor(N/K),replace=T)
  risk=riskenv+riskgen
  prob<-evi(runif(floor(N/K)))
  BreastCancer=(risk-prob>=sort(risk-prob,decreasing = T)[N])
  StudySelect=union(which(riskgen>0 & BreastCancer),sample(which(riskgen==0 & BreastCancer),N1-sum(riskgen>0 & BreastCancer)))
  t=table(riskenv[StudySelect],riskgen[StudySelect])
  CPR[sim]=1/(t[1,1]*t[2,2]/t[1,2]/t[2,1])
  #if(sim %%50==0){print(table(riskenv[StudySelect],riskgen[StudySelect]))}
  Nbrca[sim]=sum(riskgen[StudySelect]>0)
  Nulliparous[sim]=sum(riskenv[StudySelect]==0)
}
print(mean(Nulliparous))
print(81+281)
print(mean(Nbrca))
print(81+63)
parity1=numeric(3)
parity1[1]=median(CPR)
parity1[2]=quantile(CPR,0.025)
parity1[3]=quantile(CPR,0.975)
print(parity1)

#parity brca2
c=chisq.test(as.table(matrix(c(48,281,78,299),2,2)))
K=0.12
R1=K
R0=K*1.29
eff=R1*(1-R0)/R0/(1-R1)
BRCAeff=0.6/(1-0.6)/K*(1-K)
N=13000
N1=706
enveff<-c(rep(eff,floor(377*R0/R1)),rep(1,700-377))
geneff<-c(rep(BRCAeff,2),rep(1,998))
nsim=100000
CPR<-numeric(nsim)
Nbrca<-numeric(nsim)
Nulliparous<-numeric(nsim)

for(sim in 1:nsim){
  set.seed(sim)
  riskenv<-sample(log(enveff),floor(N/K),replace=T)
  riskgen<-sample(log(geneff),floor(N/K),replace=T)
  risk=riskenv+riskgen
  prob<-evi(runif(floor(N/K)))
  BreastCancer=(risk-prob>=sort(risk-prob,decreasing = T)[N])
  StudySelect=union(which(riskgen>0 & BreastCancer),sample(which(riskgen==0 & BreastCancer),N1-sum(riskgen>0 & BreastCancer)))
  t=table(riskenv[StudySelect],riskgen[StudySelect])
  CPR[sim]=1/(t[1,1]*t[2,2]/t[1,2]/t[2,1])
  #if(sim %%50==0){print(table(riskenv[StudySelect],riskgen[StudySelect]))}
  Nbrca[sim]=sum(riskgen[StudySelect]>0)
  Nulliparous[sim]=sum(riskenv[StudySelect]==0)
}
print(mean(Nulliparous))
print(48+281)
print(mean(Nbrca))
print(48+78)
parity2=numeric(3)
parity2[1]=median(CPR)
parity2[2]=quantile(CPR,0.025)
parity2[3]=quantile(CPR,0.975)
print(parity2)
