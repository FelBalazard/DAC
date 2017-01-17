# geneticrisk.R
#This script is the main analysis script of the paper. It starts by training a genetic risk estimation
#on the Wellcome Trust T1D data. This estimator is then used to predict risk of T1D on Cases of Non Autoimmune diseases
#(CNAD) of the Wellcome trust and on the patients of the Isis-Diab study. At this point, SNPs have 
#already been QCed for the training set and selected to be associated with T1D (p<10e-5)

#The following datasets are in Plink format. WT controls refer to CNAD.
ISIS <- read.table("~/Documents/G-E/ISIS.raw", quote="\"", comment.char="",header=T)
training <- read.table("~/Documents/G-E/training.raw", quote="\"", comment.char="",header=T)
WTcontrols<-read.table("~/Documents/G-E/WTcontrols.raw", quote="\"", comment.char="",header=T)

#The following are used to check that no difference in the quality of the data remains because of
#the different chips used. Results not shown.
mid_id<-read.table("~/Documents/G-E/DIAB_PB_610-1.fam", quote="\"", comment.char="",header=F)
OMNI5_id<-read.table("~/Documents/G-E/DIAB_PB_OMNI5-1.fam", quote="\"", comment.char="",header=F)

#The following is to have access to additional information, in particular, ethnic information.
identite <- read.csv("~/Documents/G-E/identite", header=T, sep=";")
identite<-identite[-1,]
idIID <-as.character(read.csv("~/Documents/G-E/idIID",header=F)[,])
idgen<-identite[na.omit(match(ISIS$IID,idIID)),]

#The frq datasets are used to check which allele is coded 1 in the previous .raw datasets. 
#Everything is flipped to correspond to the training dataset.
ISIS_frq <- read.table("~/Documents/G-E/ISIS.frq", quote="\"", comment.char="",header=T)
training_frq<-read.table("~/Documents/G-E/training.frq", quote="\"", comment.char="",header=T)
WTcontrols_frq<-read.table("~/Documents/G-E/WTcontrols.frq", quote="\"", comment.char="",header=T)
ATCG=c("A","T","C","G")
names(ATCG)<-c("T","A","G","C")
same<-(training_frq[,3]==as.character(ISIS_frq[,3]))*(training_frq[,4]==as.character(ISIS_frq[,4]))
complementary<-(training_frq[,3]==as.character(ATCG[as.character(ISIS_frq[,3])]))*(training_frq[,4]==as.character(ATCG[as.character(ISIS_frq[,4])]))
flipISIS=2*(same+complementary)-1
same<-(training_frq[,3]==as.character(WTcontrols_frq[,3]))*(training_frq[,4]==as.character(WTcontrols_frq[,4]))
complementary<-(training_frq[,3]==as.character(ATCG[as.character(WTcontrols_frq[,3])]))*(training_frq[,4]==as.character(ATCG[as.character(WTcontrols_frq[,4])]))
flipWTcontrols= 2*(same+complementary)-1

snps<-seq(7,ncol(ISIS))
for (i in which(flipISIS==-1)){
  ISIS[snps[i]]=2-ISIS[snps[i]]
}
names(ISIS)=names(training)

for (i in which(flipWTcontrols==-1)){
  WTcontrols[snps[i]]=2-WTcontrols[snps[i]]
}
names(WTcontrols)=names(training)

#This is the last piece of QC: missingness rate<5% in the ISIS data.
snp_miss_I<-apply(is.na(ISIS),2,sum)
snp_miss_tr<-apply(is.na(training),2,sum)
badsnp<-union(which(snp_miss_I>nrow(ISIS)/20),which(snp_miss_tr>nrow(training)/20))
ISIS<-ISIS[,-badsnp]
training<-training[,-badsnp]
WTcontrols<-WTcontrols[,-badsnp]

snps<-seq(7,ncol(ISIS))

#The following completes the remaining missing data by sampling from the training set.
set.seed(12)#12 is picked as seed as it is an excellent number.
for (i in snps){
  training[is.na(training[,i]),i]=sample(training[!is.na(training[,i]),i],sum(is.na(training[,i])),replace = FALSE)
  ISIS[is.na(ISIS[,i]),i]=sample(training[,i],sum(is.na(ISIS[,i])),replace = FALSE)
  WTcontrols[is.na(WTcontrols[,i]),i]=sample(training[,i],sum(is.na(WTcontrols[,i])),replace = FALSE)
}
#phenotype 0= control 1=case
WTcontrols$PHENOTYPE=0
training$PHENOTYPE=(training$PHENOTYPE+9)/10
ISIS$PHENOTYPE=1

#The genetic risk estimation is trained here and then used to predict on the Isis Diab and CNAD
library(e1071)
Pheno= as.factor(training$PHENOTYPE)
svm.model <- svm(Pheno ~ ., data = training[,snps],probability=TRUE)
ISIS.pred  <- predict(svm.model,ISIS[,snps],probability=TRUE)
WTcontrols.pred<-predict(svm.model,WTcontrols[,snps],probability=T)
#Non-Europeans are excluded from the ISIS data at this point.
ISISWT<-rbind(cbind(ISIS[idgen[,38]=="Europe","PHENOTYPE"],attributes(ISIS.pred)$probabilities[idgen[,38]=="Europe",2]),cbind(WTcontrols$PHENOTYPE,attributes(WTcontrols.pred)$probabilities[,2]))

#This is the logit function but I prefer to call it the evidence function following Jaynes's book
#Probability theory: the logic of science. 
evi<-function(x){
  return(log(x/(1-x)))
}

#Calibration of the genetic risk estimation.
calibration.model<-glm(ISISWT[,1]~evi(ISISWT[,2]),family="binomial")
ISISWT[,2]<-calibration.model$fitted.values

calibration.plot <- function(obs, pred, bins=10) {
  #  Plots a calibration plot of a set of predictions from a classifier
  #
  # Args:
  #   obs: Vector of true labels. Should be binary (0 or 1)
  #   pred: Vector of predictions of each observation from the classifier. Should be real
  #       number
  #   bins: The number of bins to use in the calibration plot
  require(plyr)
  library(Hmisc)
  
  bin.pred <- cut(pred, bins)
  
  k <- ldply(levels(bin.pred), function(x) {
    idx <- x == bin.pred
    c(sum(obs[idx]) / length(obs[idx]), mean(pred[idx]))
  })
  
  is.nan.idx <- !is.nan(k$V2)
  k <- k[is.nan.idx,]  
  plot(k$V2, k$V1, xlim=c(0,1), ylim=c(0,1),pch=3,cex.lab=1.3, xlab="Mean Prediction", ylab="Observed Fraction", type="o")
  lines(c(0,1),c(0,1), col="grey")
}

#Calibration plot of figure 2
calibration.plot(ISISWT[,1],ISISWT[,2])
#Density plot of figure 2
densityplot(~evi(ISISWT[,2]),groups=factor(ISISWT[,1]),xlim=c(-8,4),plot.points=F,
            ylab=list(label="Density",cex=1.4),xlab=list(label="Genetic risk (logit)",cex=1.4),scales = list(cex=1),lty=c(1,2),col=c(1,2),
            trellis.par.set(superpose.line = list(
              col = rep(
                c(1,2)
                ,14),
              lwd = rep( 2, 28),
              lty = rep( c(1,2), 14) )
            ),
            key = list(corner=c(0,1), 
                       lines=list( 
                         col = trellis.par.get()$superpose.line$col[1:2], 
                         lwd = trellis.par.get()$superpose.line$lwd[1:2], 
                         lty = trellis.par.get()$superpose.line$lty[1:2]
                       ),
                       text=list(c("CNAD","Isis-Diab T1D"))))

library(AUC)
auc(roc(ISISWT[,2],factor(ISISWT[,1])))
#ROCurve of figure 2
plot(roc(ISISWT[,2],factor(ISISWT[,1])),cex.lab=1.4)
text(0.4,0.7,"AUC=0.86",cex=1.4)

#ISISgen contains two identifiers for patients and calibrated genetic risk estimate.
ISISgen<-cbind(ISIS[idgen[,38]=="Europe",1:2],ISISWT[1:sum(idgen[,38]=="Europe"),2])
#AGE####
attach(idgen[idgen[,38]=="Europe",])
age=1/365.25*as.numeric(as.Date(Date.du.diagnostic,format="%Y-%m-%d")-as.Date(Date.de.naissance,format="%Y-%m-%d"))
detach(idgen[idgen[,38]=="Europe",])
summary(lm(age~evi(ISISgen[,3]) ))#First analysis announced.
confint.lm(lm(age~evi(ISISgen[,3])))#Confidence interval for estimate of effect size of association
ISISgen<-cbind(ISISgen,age)


#ENV####
#The following is to homogenize and merge the large questionnaire and the short questionnaire.
envmarkers<-read.csv("~/Documents/G-E/envmarkers",header=F)#The list of markers deemed significant in case-control study.
ISISenv<-read.csv("~/Documents/G-E/ISISall",header=T)
IIDenv<-read.csv("~/Documents/G-E/IIDall",header=F)
petitenv=read.csv("~/Documents/G-E/petitquest.csv")
petitenv=petitenv[petitenv$Type.Questionnaire=="PATIENT",]
for(name in names(petitenv)[6:ncol(petitenv)]){
  if(length(levels(factor(petitenv[,name])))!=2){
    petitenv[,name]<-(max(petitenv[,name],na.rm=T)-petitenv[,name])/(max(petitenv[,name],na.rm=T)-min(petitenv[,name],na.rm=T))
  }
}
#For the few variables ordered in a reasonable manner, we reverse the ordering back.
for(name in c("X434","X435")){
  petitenv[,name]<-1-petitenv[,name]
}
petitenv$X456.1[petitenv$X456==2]=2
petitenv$X456.2[petitenv$X456==2]=2
petitenvmarkers=as.character(envmarkers[as.character(envmarkers[,1])%in%names(petitenv),1])
names(ISISenv)[1:2]=c("Date.de.naissance","Date.de.diagnostic")

ISISenv$IID=IIDenv[,1]
Merged.env=rbind(petitenv,ISISenv[,which(names(ISISenv)%in%names(petitenv))])
ISISgenenv<-ISISgen[na.omit(match(Merged.env$IID,ISISgen$IID)),]
ISISenvgen<-Merged.env[Merged.env$IID%in%ISISgen$IID,]
for(name in names(ISISenvgen)[6:ncol(ISISenvgen)]){
  if(length(levels(factor(ISISenvgen[,name])))==2){
    ISISenvgen[,name]=as.numeric(ISISenvgen[,name]==1)
  }
}

#The following is the main analysis of the environment against genetic risk. The default is used 
#for the gam function of the mgcv package.
library(mgcv)
p<-numeric(length(names(ISISenvgen))-5)
names(p)=names(ISISenvgen)[6:54]
for(name in names(ISISenvgen)[6:54]){
  if (var(ISISenvgen[,name],na.rm=T)<0.001){
    p[name]=1
  }else{
    p[name]<-summary(gam(ISISenvgen[,name]~evi(ISISgenenv[,3])+s(ISISgenenv$age)))$p.pv[2]
  }
}

#One-sided test
for( name in petitenvmarkers){
  if (summary(gam(evi(ISISgenenv[,3])~ISISenvgen[,name]+s(ISISgenenv$age)))$p.coeff[2]>0){
    p[name]<-p[name]/2
  } else {
    p[name]<-1-p[name]/2
  }
}

#Used to have the effect size estimated in the case-control study.
result_cc <- read.csv("~/Documents/publications/isis_diab/resultat.csv")
row.names(result_cc)<-result_cc[,1]

#Environment used for simulations. Simulations are parallelized.
rm(list=as.character(ls()[!ls()%in%c("ISISenvgen","ISISWT","result_cc")]))
save.image("~/Documents/G-E/env_power.RData")