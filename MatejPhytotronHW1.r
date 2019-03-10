#clear memory
rm(list=ls())

#Load Libraries
library(lme4)
library(lmerTest)
library(MuMIn)
library(car)
library(MASS)

#Read data
gdat<-read.delim('Inputs/gasometry_M.csv',header=T,sep=';')
bdat<-read.delim('Inputs/biometric.csv',header=T,sep=',')
postHW1Asat<-read.csv('Outputs/preHWpredAsat.csv')
postHW1Cond<-read.csv('Outputs/preHWpredCond.csv')

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

#get post HW1 data
postHW1<-gdat[which(DOX>15&DOX<34),]


#assemble data for photosynthesis, Conductance, WUE
Tree<-factor(postHW1$Tree)
Ttr<-factor(paste(postHW1$Temp,postHW1$VPD,sep=''))
CO2<-factor(postHW1$CO2)
CO2[which(CO2=='C')]<-'AC'
CO2<-factor(CO2)
Wtr<-factor(postHW1$WD)
Asat<-postHW1$Photo
Cond<-postHW1$Cond
WUE<-postHW1$WUE
Day<-DOX[DOX>15&DOX<34]
postHW1dat<-data.frame(Tree,Ttr,Wtr,CO2,Day,Asat,Cond,WUE)
for(i in 1:4) postHW1dat[,i]<-factor(postHW1dat[,i]) #make sure factors are assigned as such

#get predicted Asat and Cond
TreeID<-Tree:(Wtr:Ttr)
PredCondTree<-numeric()
PredAsatTree<-numeric()
for(i in 1:dim(postHW1dat)[1]){
  PredAsatTree<-c(PredAsatTree,postHW1Asat[which(postHW1Asat[,1]==TreeID[i]),2])
  PredCondTree<-c(PredCondTree,postHW1Cond[which(postHW1Cond[,1]==TreeID[i]),2])
}

postHW1dat<-data.frame(postHW1dat,PredAsatTree,PredCondTree)


#fit LMM for Asat
fitasat<-lmer(Asat~(Ttr)+(Wtr)+(CO2)+Wtr:CO2+Wtr:CO2+PredAsatTree+(1|Ttr/Wtr/Tree),postHW1dat,REML=F,na.action="na.fail")
print(dredge(fitasat,rank='AICc'))


print(dredge(fitasat,rank='BIC'))

finalasat<-lmer(Asat~(Wtr)+(1|Ttr/Wtr/Tree),postHW1dat,REML=F,na.action="na.fail")
print(summary(finalasat))
print(coef(finalasat))
print(r.squaredGLMM(finalasat))

tiff(filename="Outputs/PostHW1Asat.tif", width=1240, height=560, units='px', type='cairo')
plot(Ttr:CO2:Wtr,Asat,xaxt='n')
predAsat<-fitted(finalasat)
points(Ttr:CO2:Wtr,predAsat,col=2,cex=2)
#write.csv(predAsat,'Outputs/PreHWpredAsat.csv',row.names=F)
par(cex.axis=.6)
axis(1,cex=.7,labels=unique(Ttr:CO2:Wtr),at=c(1:2,5:12))
dev.off()

#fit LMM for Conductance
fitasat<-lmer(Cond~(Ttr)+(Wtr)+(CO2)+Wtr:CO2+Wtr:CO2+PredCondTree+(1|Ttr/Wtr/Tree),postHW1dat,REML=F,na.action="na.fail")
print(dredge(fitasat,rank='AICc'))


print(dredge(fitasat,rank='BIC'))

finalasat<-lmer(Cond~PredCondTree+(1|Ttr/Wtr/Tree),postHW1dat,REML=F,na.action="na.fail")
print(summary(finalasat))
print(coef(finalasat))
print(r.squaredGLMM(finalasat))

tiff(filename="Outputs/PostHW1Cond.tif", width=1240, height=560, units='px', type='cairo')
plot(Ttr:CO2:Wtr,Cond,xaxt='n')
predCond<-fitted(finalasat)
points(Ttr:CO2:Wtr,predCond,col=2,cex=2)
par(cex.axis=.6)
axis(1,cex=.7,labels=unique(Ttr:CO2:Wtr),at=c(1:2,5:12))
dev.off()


tiff(filename="Outputs/PostHW1CondbyTree.tif", width=560, height=560, units='px', type='cairo')
xx<-which(Wtr=='WET')
plot(postHW1dat$PredCondTree[xx],postHW1dat$Cond[xx],type='p',col=4,xlim=c(0,0.6),ylim=c(0,0.6),xlab='PreHW1 Cond',ylab='PostHW1 Cond')
points(postHW1dat$PredCondTree[-xx],postHW1dat$Cond[-xx],type='p',col=2)
legend('topright',c('Wet','Dry'),pch=1,col=c(4,2))
dev.off()
