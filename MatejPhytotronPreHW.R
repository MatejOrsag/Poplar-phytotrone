#clear memory
rm(list=ls())

#Load Libraries
library(lme4)
library(lmerTest)
library(MuMIn)
library(car)
library(MASS)

#Read data
gdat<-read.delim('Inputs/gasometry_M.txt',header=T,sep=';')
bdat<-read.delim('Inputs/biometric.csv',header=T,sep=',')

#Generate Day of Experiment for gas exchange
DOX<-rep(NA,dim(gdat)[1])
DOX[which(gdat$Date=='29.7.2014')]<-7
DOX[which(gdat$Date=='1.8.2014')]<-10
DOX[which(gdat$Date=='5.8.2014')]<-14
DOX[which(gdat$Date=='12.8.2014')]<-21
DOX[which(gdat$Date=='19.8.2014')]<-28
DOX[which(gdat$Date=='26.8.2014')]<-35
DOX[which(gdat$Date=='2.9.2014')]<-42
#plot(DOX,type='l')#visual check

#get preHW data
preHW<-gdat[which(DOX<15),]

#assemble data for photosynthesis, Conductance, WUE
Tree<-factor(preHW$Tree)
Ttr<-factor(paste(preHW$Temp,preHW$VPD,sep=''))
CO2<-factor(preHW$CO2)
CO2[which(CO2=='C')]<-'AC'
CO2<-factor(CO2)
Wtr<-factor(preHW$WD)
Asat<-preHW$Photo
Cond<-preHW$Cond
WUE<-preHW$WUE
Day<-DOX[DOX<15]
preHWdat<-data.frame(cbind(Tree,Ttr,Wtr,CO2,Day,Asat,Cond,WUE))
for(i in 1:4) preHWdat[,i]<-factor(preHWdat[,i]) #make sure factors are assigned as such
#Check Normality
#qqp(preHWdat$Asat,'norm')

#fit LMM
fitasat<-lmer(Asat~(Ttr)+(Wtr)+(CO2)+Wtr:CO2+Wtr:CO2+(1|Ttr/Wtr/Tree),data=preHWdat,REML=F,na.action="na.fail")
print(dredge(fitasat,rank='AICc'))


print(dredge(fitasat,rank='BIC'))

finalasat<-lmer(Asat~(Wtr)+(CO2)+(1|Ttr/Wtr/Tree),data=preHWdat,REML=F,na.action="na.fail")
print(summary(finalasat))
print(coef(finalasat))
print(r.squaredGLMM(finalasat))

tiff(filename="Outputs/PreHWdataAsat.tif", width=1240, height=560, units='px', type='cairo')
plot(Ttr:CO2:Wtr,Asat,xaxt='n')
points(Ttr:CO2:Wtr,fitted(finalasat),col=2,cex=2)
par(cex.axis=.6)
axis(1,cex=.7,labels=unique(Ttr:CO2:Wtr),at=c(1:2,5:12))
dev.off()

#fit LMM for Conductance
fitasat<-lmer(Cond~(Ttr)+(Wtr)+(CO2)+Wtr:CO2+Wtr:CO2+(1|Ttr/Wtr/Tree),data=preHWdat,REML=F,na.action="na.fail")
print(dredge(fitasat,rank='AICc'))


print(dredge(fitasat,rank='BIC'))

finalasat<-lmer(Cond~(Wtr)+(CO2)+(1|Ttr/Wtr/Tree),data=preHWdat,REML=F,na.action="na.fail")
print(summary(finalasat))
print(coef(finalasat))
print(r.squaredGLMM(finalasat))

tiff(filename="Outputs/PreHWdataCond.tif", width=1240, height=560, units='px', type='cairo')
plot(Ttr:CO2:Wtr,Cond,xaxt='n')
points(Ttr:CO2:Wtr,fitted(finalasat),col=2,cex=2)
par(cex.axis=.6)
axis(1,cex=.7,labels=unique(Ttr:CO2:Wtr),at=c(1:2,5:12))
dev.off()

#fit LMM for WUE
preHWdat2<-preHWdat
preHWdat2<-preHWdat2[Cond>0.001,]

fitasat<-lmer(WUE~(Ttr)+(Wtr)+(CO2)+Wtr:CO2+Wtr:CO2+(1|Ttr/Wtr/Tree),data=preHWdat2,REML=F,na.action="na.fail")
print(dredge(fitasat,rank='AICc'))


print(dredge(fitasat,rank='BIC'))

finalasat<-lmer(WUE~(Wtr)+(CO2)+(Wtr):(CO2)+(1|Ttr/Wtr/Tree),data=preHWdat2,REML=F,na.action="na.fail")
print(summary(finalasat))
print(coef(finalasat))
print(r.squaredGLMM(finalasat))

tiff(filename="Outputs/PreHWdataWUE.tif", width=1240, height=560, units='px', type='cairo')
plot(preHWdat2$Ttr:preHWdat2$CO2:preHWdat2$Wtr,preHWdat2$WUE,xaxt='n')
points(preHWdat2$Ttr:preHWdat2$CO2:preHWdat2$Wtr,fitted(finalasat),col=2,cex=2)
par(cex.axis=.6)
axis(1,cex=.7,labels=unique(Ttr:CO2:Wtr),at=c(1:2,5:12))
dev.off()
