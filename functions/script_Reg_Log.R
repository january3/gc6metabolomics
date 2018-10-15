rm(list = ls())
# library(knitr)
# setwd("D:/Projekt_GC6/manuscript_21_naturemed/includes/")
# purl("calculations.rmd")
setwd("D:/Projekt_GC6/Data_log_reg/")
load("Reg_Log_Data.RData")
# File  Reg_Log_Data.RData contain variable cutoff=5 and rgR

####### Reg.Log to all at Proximate state
#data prepare
data<-rgR$targets[rgR$targets$TIME_TO_TB < cutoff,c("group","SITE","TIME_TO_TB","SEX","SAMPLE_TYPE","CASE_GROUP")]
data[,(dim(data)[2]+1):(dim(data)[2]+dim(rgR$E)[1])]<-as.data.frame(t(rgR$E[,rgR$targets$TIME_TO_TB < cutoff]))
#regression with output calc
res<-matrix(nrow=dim(rgR$genes)[1],ncol=8)
for (i in 1:dim(rgR$genes)[1]){
tmp<-data[, c(1,i+6)] 
regress<-glm(data$group ~. ,data = tmp, family = binomial(logit))
res[i,c(1,2,3)]<-summary(regress)$coefficients[2,c(1,2,4)]# coef for fea
res[i,4]<-1-pchisq(regress$null.deviance-regress$deviance, regress$df.null-regress$df.residual) #LRT test
res[i,6]<-exp(res[i,1]) #OR
res[i,7]<-exp(res[i,1]-(qnorm(0.975)*res[i,2])) #OR CI down
res[i,8]<-exp(res[i,1]+(qnorm(0.975)*res[i,2])) #OR CI up
}
res[,5]<-p.adjust(res[,4],method = "BH")

#######Reg.Log to all at ALL state
#data prepare
data<-rgR$targets[,c("group","SITE","TIME_TO_TB","SEX","SAMPLE_TYPE","CASE_GROUP")]
data[,(dim(data)[2]+1):(dim(data)[2]+dim(rgR$E)[1])]<-as.data.frame(t(rgR$E))
#regression with output calc
res.2<-matrix(nrow=dim(rgR$genes)[1],ncol=8)
for (i in 1:dim(rgR$genes)[1]){
  tmp<-data[, c(1,i+6)] 
  regress<-glm(data$group ~. ,data = tmp, family = binomial(logit))
  res.2[i,c(1,2,3)]<-summary(regress)$coefficients[2,c(1,2,4)]# coef for fea
  res.2[i,4]<-1-pchisq(regress$null.deviance-regress$deviance, regress$df.null-regress$df.residual) #LRT test
  res.2[i,6]<-exp(res.2[i,1]) #OR
  res.2[i,7]<-exp(res.2[i,1]-(qnorm(0.975)*res[i,2])) #OR CI down
  res.2[i,8]<-exp(res.2[i,1]+(qnorm(0.975)*res[i,2])) #OR CI up
}
res.2[,5]<-p.adjust(res.2[,4],method = "BH")


# Results bind and selection
colnames(res)<-colnames(res.2)<-c("Beta","Beta.SE","Beta.pval","LRT.pval","LRT.qval","OR","OR.CI.low","OR.CI.up")
res.all<-cbind(res,res.2)
rownames(res.all)<-rgR$genes$BIOCHEMICAL
sig<-res.all[which(res[,5]<0.01 | res.2[,5]<0.01),] # extract of relevant fea



