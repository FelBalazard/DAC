#Script used to plot figure 3.
result_auc <- read.table("~/Documents/publications/case-only gene environment/result/result_auc", quote="\"", comment.char="")
result_K <- read.table("~/Documents/publications/case-only gene environment/result/result", quote="\"", comment.char="")
par(mfrow=c(1,2),mar=c(5,4.5,1,.5))
result_K_m=matrix(rep(0,36),12,3)
i=1
for(K in c(0.002,0.006,0.01)){
  for(N in c(500,1500,3000,5000)){
    result_K_m[i,1]=K
    result_K_m[i,2]=N
    result_K_m[i,3]=mean(result_K[result_K[,1]==K & result_K[,2]==N,5])
    i=i+1}
}

plot(result_K_m[,2],result_K_m[,3],ylim=c(0,1),pch=3,xlab="Sample size (case-only)",ylab="Power")
lines(result_K_m[1:4,2],result_K_m[1:4,3],lty=as.numeric(factor(result_K_m[1:4,1])))
lines(result_K_m[5:8,2],result_K_m[5:8,3],lty=as.numeric(factor(result_K_m[,1])[5:8]))
lines(result_K_m[9:12,2],result_K_m[9:12,3],lty=as.numeric(factor(result_K_m[,1])[9:12]))
legend(500,1,bty='n',lty=c(3,2,1),c("K=1%","K=0.6%","K=0.2%"))

result_auc_m=matrix(rep(0,48),16,3)
i=1
auc=0.86
for(N in c(500,1500,3000,5000)){
  result_auc_m[i,1]=auc
  result_auc_m[i,2]=N
  result_auc_m[i,3]=result_K_m[result_K_m[,1]==0.002 & result_K_m[,2]==N,3]
  i=i+1
}
for(auc in as.numeric(levels(factor(result_auc[,3])))){
  for(N in c(500,1500,3000,5000)){
    result_auc_m[i,1]=auc
    result_auc_m[i,2]=N
    result_auc_m[i,3]=mean(result_auc[result_auc[,3]==auc & result_auc[,2]==N,5])
    i=i+1}
}
plot(result_auc_m[,2],result_auc_m[,3],ylim=c(0,1),pch=3,xlab="Sample size (case-only)",ylab=" ")
lines(result_auc_m[1:4,2],result_auc_m[1:4,3],lty=as.numeric(factor(result_auc_m[1:4,3])))
lines(result_auc_m[5:8,2],result_auc_m[5:8,3],lty=as.numeric(factor(result_auc_m[,3])[5:8]))
lines(result_auc_m[9:12,2],result_auc_m[9:12,3],lty=as.numeric(factor(result_auc_m[,3])[9:12]))
lines(result_auc_m[13:16,2],result_auc_m[13:16,3],lty=as.numeric(factor(result_auc_m[,3])[13:16]))
legend(500,1,bty='n',lty=c(4,3,2,1),c("AUC=0.92","AUC=0.90","AUC=0.88","AUC=0.86"))

