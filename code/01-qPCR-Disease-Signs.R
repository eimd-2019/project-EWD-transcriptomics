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


#upload data
#setwd()
Adult36 <- read.delim("36hourAdult.txt")
Adult36$Temp<-ifelse(Adult36$Temp=='Cold', '11 °C', '18 °C')
Adult5<-read.delim("5dayAdult.txt")
Adult5$Temp<-ifelse(Adult5$Temp=='Cold', '11 °C', '18 °C')



#####After 36 hours, impact of treatments on disease presence and pathogen load
m3a<-glmer(diseased~Laby*Temp+(1|ID),data=Adult36, na.action = na.omit, family=binomial, nAGQ=0, glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
summary(m3a)

m3b<-lmer(log10cells~Laby*Temp+(1|ID),na.action = na.omit, data=Adult36)
summary(m3b)



###Graph after 36hours
#Disease signs
Adult36$predicteddisease <- predict(m3a, re.form=NA, type="response")

Adult36Sum<-ddply(Adult36, .(Temp, Laby), summarize, meanDisease=mean(predicteddisease), seDisease=sd(predicteddisease)/sqrt(length(predicteddisease)))


p<-ggplot(Adult36Sum, aes(x=Laby, y=meanDisease, fill=Temp))+geom_col(position = "dodge")+
  geom_errorbar(aes(ymin=meanDisease-seDisease, ymax=meanDisease+seDisease), width=.2,
                position=position_dodge(.9)) +scale_fill_manual(values=c("blue", "red"))+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p+labs(x=expression(paste(italic("Labyrinthula"), "treatment")), y=expression(paste("Predicted prevalence of disease signs")))

#################################
##After 5 days, impact on disease presence 
Adult5<-Adult5[!is.na(Adult5$Diseased),]
m3c<-glmer(Diseased~Laby*Temp+(1|ID),data=Adult5, nAGQ=0, na.action = na.omit,family=binomial)
summary(m3c)


Adult5$predicteddisease <- predict(m3c, re.form=NA, type="response")

Adult5Sum<-ddply(Adult5, .(Temp, Laby), summarize, meanDisease=mean(predicteddisease), seDisease=sd(predicteddisease)/sqrt(length(predicteddisease)))

#graph
p<-ggplot(Adult5Sum, aes(x=Laby, y=meanDisease, fill=Temp))+geom_col(position = "dodge")+
  geom_errorbar(aes(ymin=meanDisease-seDisease, ymax=meanDisease+seDisease), width=.2,
                position=position_dodge(.9)) +scale_fill_manual(values=c("blue", "red"))+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p+labs(x=expression(paste(italic("Labyrinthula"), "treatment")), y=expression(paste("Predicted prevalence of disease signs")))



#Graph disease sign data for days 2 and 5 together

Adult36Sum<-ddply(Adult36, .(Temp, Laby, ID), summarize, meanDisease=sum(diseased, na.rm=TRUE)/length(diseased))
Adult36Sum$Day<-'Day 2'

Adult5Sum<-ddply(Adult5, .(Temp, Laby, ID), summarize, meanDisease=sum(Diseased, na.rm=TRUE)/length(Diseased))
Adult5Sum$Day<-'Day 5'

AdultSum<-rbind(Adult36Sum, Adult5Sum)

AdultSum2<-ddply(AdultSum, .(Temp, Laby, Day), summarize, meanDisease2=mean(meanDisease, na.rm=TRUE), seDisease=sd(meanDisease)/2)
AdultSum2$Day<-ifelse(AdultSum2$Day=='Day 2', '2 days post-exposure', '5 days post-exposure')

p<-ggplot(AdultSum2, aes(x=Laby, y=meanDisease2, fill=Temp))+geom_col(position = "dodge")+
  geom_errorbar(aes(ymin=meanDisease2-seDisease, ymax=meanDisease2+seDisease), width=.2, position=position_dodge(.9))+
  scale_fill_manual(values=c("blue", "red"))+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p+facet_grid(.~Day)+labs(x=expression(paste(italic("Labyrinthula"), " treatment")), y=expression(paste("Prevalence of disease signs")))+scale_y_continuous(limits = c(0,1.1), expand = c(0, 0))+ theme(legend.title = element_blank())+ theme(legend.position = c(0.1, 0.8))



#qPCR results


#qPCR results graph
Adult36$predictedlog10load<- predict(m3b, re.form=NA, type="response")
Adult36Sum<-ddply(Adult36, .(Temp, Laby), summarize, meanpredictedlog10=mean(predictedlog10load), meanlog10=log10(mean(cells)), selog10cells=(sd(log10(cells)))/sqrt(length(cells)))
Adult36Sum$Temp<-ifelse(Adult36Sum$Temp=='11 °C', '11 °C', '18 °C')

p<-ggplot(Adult36Sum, aes(x=Laby, y=meanlog10, fill=Temp))+geom_col(position = "dodge")+
  geom_errorbar(aes(ymin=meanlog10-selog10cells, ymax=meanlog10+selog10cells), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c("blue", "red"))+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p+labs(x=expression(paste(italic("Labyrinthula"), " treatment")), y=expression(paste("log"[10]* " (", italic('Labyrinthula '), 'cells/ mg eelgrass tissue)')))+scale_y_continuous(limits = c(0,5), expand = c(0, 0))+ theme(legend.title = element_blank())+ theme(legend.position = c(0.2, 0.8))