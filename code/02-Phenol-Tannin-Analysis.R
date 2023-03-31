rm(list =ls())

library(ggplot2)
library (ordinal)
library(glmm)
library(lme4)
library(lmerTest)
library(optimx)
library(plyr)
require(mdthemes)
require(readr)
require(gridExtra)


##Phenol and tannin analysis

PhenolsTannins <- read_csv("PhenolsTannins.csv")
PhenolsTannins<-PhenolsTannins[PhenolsTannins$Temp!='Varying',]

PhenolsTannins$CT_prop_dm<-PhenolsTannins$CT_percent_dm/100
PhenolsTannins$Phenols_prop_dm<-PhenolsTannins$Phenol_percent_dm/100

PhenolsTanninsD2<-PhenolsTannins[PhenolsTannins$SampleDay<3,]
PhenolsTanninsD5<-PhenolsTannins[PhenolsTannins$SampleDay>3,]

#phenol results
#36 hours
m1<-lmer(Phenol_percent_dm~Laby*Temp+(1|ID),data=PhenolsTanninsD2, na.action = na.omit)
summary(m1)


#5 days
m2<-lmer(Phenol_percent_dm~Laby*Temp+(1|ID),data=PhenolsTanninsD5, na.action = na.omit)
summary(m2)

#Condensed Tannin results
#36 hours
m1<-lmer(CT_percent_dm~Laby*Temp+(1|ID),data=PhenolsTanninsD2, na.action = na.omit)
summary(m1)


#5 days
m2<-lmer(CT_percent_dm~Laby*Temp+(1|ID),data=PhenolsTanninsD5, na.action = na.omit)
summary(m2)

#graph results
PhenolsTanninsSummary<-ddply(PhenolsTannins, .(Temp, Laby, SampleDay, Treatment), summarize, meanCT=mean(CT_percent_dm, na.rm=TRUE), sdCT=sd(CT_percent_dm, na.rm=TRUE), meanPhenols=mean(Phenol_percent_dm, na.rm=TRUE), sdPhenols=sd(Phenol_percent_dm, na.rm=TRUE))
PhenolsTanninsSummary$SampleDay<-ifelse(PhenolsTanninsSummary$SampleDay=='5', '5 days post-exposure', '2 days post-exposure') 


p2<-ggplot(PhenolsTanninsSummary, aes(x=Laby, y=meanCT, fill=Temp), color="black")+geom_col(position=position_dodge(), color="black")+
  geom_errorbar(aes(ymin=meanCT-sdCT, ymax=meanCT+sdCT), width=.2,
                position=position_dodge(.9), color="black") +
  scale_fill_manual(values=c('blue','red'), name="Temperature", labels=c("11 °C", "18 °C"))
p2<-p2+facet_grid(.~SampleDay)+theme_bw()+ylab("Condensed Tannins \n(% Dry Mass)")+labs(x=expression(paste(italic("Labyrinthula"), " treatment")))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position = c(0.9, 0.85), legend.title = element_blank())+scale_y_continuous(expand = c(0, 0), limits=c(NA, 0.6))


p1<-ggplot(PhenolsTanninsSummary, aes(x=Laby, y=meanPhenols, fill=Temp), color="black")+geom_col(position=position_dodge(), color="black")+
  geom_errorbar(aes(ymin=meanPhenols-sdPhenols, ymax=meanPhenols+sdPhenols), width=.2,
                position=position_dodge(.9), color="black") +
  scale_fill_manual(values=c('blue','red'), name="Temperature", labels=c("11 °C", "18 °C"))+ylab("Phenols\n(% Dry Mass)")
p1<-p1+facet_grid(.~SampleDay)+theme_bw()+theme_bw() +
  labs(x=expression(paste(italic("Labyrinthula"), " treatment")))+ylab("Total Phenolic Compounds \n(% Dry Mass)")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(legend.position="none")+scale_y_continuous(expand = c(0, 0), limits=c(NA, 1.9))

grid.arrange(p2, p1, ncol=1)


