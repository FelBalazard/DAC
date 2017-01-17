#Script used to plot figure 3.
result_auc <- read.table("~/Documents/publications/valid env markers gen risk/results/result_auc", quote="\"", comment.char="")
result_K <- read.table("~/Documents/publications/valid env markers gen risk/results/result", quote="\"", comment.char="")
par(mfrow=c(1,2),mar=c(5,4.5,1,.5))
plot(result_K[,2],result_K[,4],ylim=c(0,1),pch=3,xlab="Sample size (case-only)",ylab="Power")
lines(result_K[1:4,2],result_K[1:4,4],lty=as.numeric(factor(result_K[1:4,1])))
lines(result_K[5:8,2],result_K[5:8,4],lty=as.numeric(factor(result_K[,1])[5:8]))
lines(result_K[9:12,2],result_K[9:12,4],lty=as.numeric(factor(result_K[,1])[9:12]))
legend(500,1,bty='n',lty=c(3,2,1),c("K=1%","K=0.6%","K=0.2%"))

plot(result_auc[,2],result_auc[,4],ylim=c(0,1),pch=3,xlab="Sample size (case-only)",ylab=" ")
lines(result_auc[1:4,2],result_auc[1:4,4],lty=as.numeric(factor(result_auc[1:4,3])))
lines(result_auc[5:8,2],result_auc[5:8,4],lty=as.numeric(factor(result_auc[,3])[5:8]))
lines(result_auc[9:12,2],result_auc[9:12,4],lty=as.numeric(factor(result_auc[,3])[9:12]))
lines(result_auc[13:16,2],result_auc[13:16,4],lty=as.numeric(factor(result_auc[,3])[13:16]))
legend(500,1,bty='n',lty=c(4,3,2,1),c("AUC=0.92","AUC=0.90","AUC=0.88","AUC=0.86"))

