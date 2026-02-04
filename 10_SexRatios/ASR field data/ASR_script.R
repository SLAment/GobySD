#tHIS CODE IS FOR THE GOBY ADULT SEX RATIO MANUSCRIPT
#Analysing Adult sex ratio data from the field

#setwd("C://Users//ivainm//Working_Document//DYNAMAR//popGen//data_analysis//sexID")
setwd("C://Users//ivmar9032//Documents//DYNAMAR local//Manuscripts//sex_det//data_analysis")


# Package import----

library(effects)
library(car)
library(ggplot2)
library(dplyr) 
library(tidyverse)
library(reshape2) 
library(lubridate)
library(lme4)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#____________________----
# 2022 ASR Field data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##Import OSR data ---- 
#with field observations of population census and additional info on the transects
OSRdata<-read.csv("./sexID/OSRdata.csv")
OSRdata$area<-as.factor(OSRdata$area)
OSRdata$area<- factor(OSRdata$area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
OSRdata$timepoint2<-as.factor(OSRdata$timepoint2)
levels(OSRdata$timepoint2)<-c("Early season (1)" , "Late season (3)", "Mid-season (2)")
OSRdata$timepoint2<- factor(OSRdata$timepoint2, levels=c("Early season (1)" ,  "Mid-season (2)", "Late season (3)"))
OSRdata$date<-as.Date(OSRdata$date)
OSRdata$date<-OSRdata$date %m-% years(2)
#Extract the count of males and females and format in long (optional)
#longdensexf<-OSRdata[,c("area", "transectlength","Nf")]
#longdensexm<-OSRdata[,c("area", "transectlength","Nm")]
#colnames(longdensexf)<-c("area","transectlength","N")
#colnames(longdensexm)<-c("area","transectlength","N")
#longdensex<-rbind(cbind(longdensexf,sex="f"),cbind(longdensexm,sex="m"))
#longdensex$sex<-as.factor(longdensex$sex)
#Add latitude info
#longdensex<-merge(longdensex,Coord[,1:2], by.x="area", by.y="Area")
#head(longdensex)

str(OSRdata)
sum(OSRdata$Census)

## Import catch data----
#Data from catch of fish and phenotyping and land, to be compared to the census data

catch<-read.csv("./sexID/sexrpheno.csv")
catch$area<-as.factor(catch$area)
catch$area<- factor(catch$area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
catch$timepoint<-as.factor(catch$timepoint)
levels(catch$timepoint)<-c("Early season (1)" , "Mid-season (2)","Late season (3)")
str(catch)

catch$N<-catch$f1+catch$f2+catch$f3+catch$m+catch$u
sum(catch$N)
## STATISTICS---------
#Here we analyse the 2022 data
#The idea is that we have sites/sublocations with temporal trajectories.
#We use the census data, Adult Sex Ratio ASR is the response variable
#date is a continuous predictor. We have to go quadratic for the temporal trend

###Prepare data----
str(OSRdata)
OSRdata$location<-as.factor(OSRdata$location)
OSRdata$yearday<-yday(OSRdata$date)
OSRdata$yearday<-scale(OSRdata$yearday)
summary(OSRdata$yearday)
#Only keep locations that have 3 time points
table(OSRdata$location,OSRdata$timepoint)
OSRdata_cut<-OSRdata[OSRdata$location!="Litj-sorroy" &
                       OSRdata$location!="Saltkjaerholmane" ,]

### Model choice----
#Area as a fixed effect so we can see what happens by pop, location as random effect if needed
#Try different random effect structures and compare
modelASR1<-glm(cbind(Nm,Nf)~(yearday+I(yearday^2))*(area), data=OSRdata, family="binomial")
modelASR2<-glmer(cbind(Nm,Nf)~(yearday+I(yearday^2))*(area)+(yearday|location), data=OSRdata, family="binomial")
modelASR3<-glmer(cbind(Nm,Nf)~(yearday+I(yearday^2))*(area)+(1|location), data=OSRdata, family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(modelASR1)
summary(modelASR2)
summary(modelASR3)
#The model with random slopes has better support AIC/BIC
Anova(modelASR2, type=3)
plot(effect("yearday:area",x.var="yearday",modelASR2), multiline=T)

### model fit for plot ----
fittime<-as.data.frame(effect("yearday:area",modelASR2,x.var="yearday", xlevels=20,se=list(compute=T,level=0.95)))
fittime$trueday<-as.Date(mean(yday(OSRdata$date))+fittime$yearday*sd(scale(OSRdata$date,scale=F)), origin = "2022-01-01")
#I would like to avoid fitting a line where there is no data
#get first sampling date for each loc
startdates<- OSRdata %>%
  group_by(area) %>%
  summarise(start=min(date))
fittimecut<-fittime[fittime$area=="Kristineberg" & fittime$trueday> startdates[startdates$area=="Kristineberg",]$start-5,]
for (i in 2:6) {fittimecut<-rbind(fittimecut, fittime[fittime$area==levels(startdates$area)[i] & fittime$trueday> startdates[startdates$area==levels(fittime$area)[i],]$start-5,])}



## PLOTS ----
Areacolors<-c("red3",  "brown1",  "darkorchid4",  "mediumpurple1",  "blue3",  "steelblue1")

### FIGURE 7 ---- 
#(without face wrap, with facet wrap for supplementary )
#Final plot with data and model fit
ggplot(fittimecut,aes(color=fct_rev(area),x=trueday,y=fit))+
  geom_point(data=OSRdata,aes(x=date,y=ASR,fill=fct_rev(area), size=Census),shape=21,alpha=0.8)+
  geom_line(size=2, alpha=0.5)+
  geom_hline(yintercept=0.5,alpha=0.6)+
  geom_ribbon(data=fittimecut, aes( ymin=lower, ymax=upper, fill=area), color=NA,alpha=0.1)+
  scale_fill_manual("Population",values=rev( Areacolors))+
  scale_color_manual("Population",values=rev(Areacolors))+
  scale_size_continuous("Sample size")+
  #facet_wrap(~fct_rev(area))+
  labs(x="Date, year 2022", y="Adult Sex Ratio", title="")+
  #theme(text=element_text(size=15))+
  theme_bw(base_size=14)
ggsave("..//figures//ASR_Fig6.pdf", width = 8, height = 5)

#suppl figure for 2022 ASR
ggplot(fittimecut,aes(color=fct_rev(area),x=trueday,y=fit))+
  geom_point(data=OSRdata,aes(x=date,y=ASR,fill=fct_rev(area), size=Census),shape=21,alpha=0.8)+
  geom_line(size=2, alpha=0.5)+
  geom_hline(yintercept=0.5,alpha=0.6)+
  geom_ribbon(data=fittimecut, aes( ymin=lower, ymax=upper, fill=area), color=NA,alpha=0.1)+
scale_fill_manual("Population",values=rev( Areacolors))+
scale_color_manual("Population",values=rev(Areacolors))+
  scale_size_continuous("Sample size")+
 facet_wrap(~fct_rev(area))+
  labs(x="Date, year 2022", y="Adult Sex Ratio", title="")+
  #theme(text=element_text(size=15))+
  theme_bw(base_size=14)
ggsave("..//figures//ASR_Suppl.pdf", width = 8, height = 5)


#Catch data, maybe enough to show that they correlate with census data?
str(OSRdata)
str(catch)
catch$location<-as.factor(catch$location)
OSRdata$timepoint<-as.factor(OSRdata$timepoint)
catch$fem<-catch$f1+catch$f2+catch$f3
levels(catch$timepoint)<-c("1","2","3")
censcatch<-merge(catch[,c("fem","m","location","timepoint","area")],
                 OSRdata[,c("Nm","Nf","location","timepoint","area","date")],
                 by=c("location","timepoint","area"))

## Census vs Catch Figure----
#plot showing correlation of Observation and catch
ggplot(censcatch,aes(x=Nm/(Nm+Nf),y=m/(m+fem),fill=fct_rev(area)))+
  geom_point(shape=21,aes(size=m+fem))+
  scale_fill_manual("Population",values=rev( Areacolors))+
  scale_size_continuous("Catch sample size")+
  geom_hline(yintercept=0.5)+geom_vline(xintercept=0.5)+
  theme(text=element_text(size=15))+
  ylim(0,0.8)+xlim(0,0.8)+
  theme_bw(base_size=14)+
  labs(x="ASR - Observation transect", y="ASR - Catch data")
ggsave("..//figures//SUppl_Obs_vs_Catch.pdf", width =7, height = 5)

cor.test(censcatch$Nm/(censcatch$Nm+censcatch$Nf),censcatch$m/(censcatch$m+censcatch$fem))

#plot(effect("length:timepoint:sex",modcondsize,x.var="timepoint" ), multiline=T)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#____________________----
# "Historical" data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Import and merge various datasets ----

# Bergen data from AC and Karen 2019 (beach seine)
Bergen2019<-read.delim("./sexID/Bergen_2019.txt", header=T)
Bergen2019$date<-paste(as.character(Bergen2019$Date),as.character(Bergen2019$Year),sep=".")
Bergen2019$date<-as.Date(Bergen2019$date,"%d.%m.%Y")
str(Bergen2019)
Bergen2019

#Kberg and Bergen data from Elisabet 2007
KB2007<-read.delim("./sexID/K_B_2007.txt", header=T)
KB2007$date<-paste(KB2007$day,KB2007$month,KB2007$year, sep="/")
KB2007$date<-as.Date(KB2007$date,"%d/%B/%Y")
#Collapse the repeated observations into sums?
SumF<-aggregate(KB2007$Tot.F, by=list(site=KB2007$site, date=KB2007$date), FUN=sum)
SumM<-aggregate(KB2007$Tot.M, by=list(site=KB2007$site, date=KB2007$date), FUN=sum)
colnames(SumF)[3]<-"Females"
Sums<-cbind(SumF,Males=SumM[,3])
Sums
KB7_s<-merge(x=Sums,y=KB2007[,c(1,2,8)],by=c("site","date"),all.x=F)
KB7_s<-unique(KB7_s)
str(KB7_s)
KB7_s
KB2007

#Kberg data from 2010
K2010<-read.delim("./sexID/K_2010.txt", header=T)
K2010$date<-as.Date(K2010$Date,"%d/%m/%Y")
#Collapse the repeated observations into sums?
SumF<-aggregate(K2010$Sum.F, by=list(site=K2010$Locality, date=K2010$date), FUN=sum)
SumM<-aggregate(K2010$Sum.M, by=list(site=K2010$Locality, date=K2010$date), FUN=sum)
colnames(SumF)[3]<-"Females"
Sums<-cbind(SumF,Males=SumM[,3])
K2010_s<-cbind(Sums,Location="Kristineberg")
str(K2010_s)

#Kberg data from 2008
K2008<-read.delim("./sexID/K_2008.txt", header=T)
K2008$date<-as.Date(K2008$date,"%d/%m/%Y")
#Collapse the repeated observations into sums?
SumF<-aggregate(K2008$Sum.F, by=list(site=K2008$Locality, date=K2008$date), FUN=sum)
SumM<-aggregate(K2008$Sum.M, by=list(site=K2008$Locality, date=K2008$date), FUN=sum)
colnames(SumF)[3]<-"Females"
Sums<-cbind(SumF,Males=SumM[,3])
K2008_s<-cbind(Sums,Location="Kristineberg")
K2008_s

#Tvarminne 2008 (Finland)
F08<-read.delim("./sexID/Tvarminne_2008.txt", header=T)
F08$date<-as.Date(F08$Date,"%d/%m/%Y")
#Collapse the repeated observations into sums?
SumF<-aggregate(F08$F.total, by=list(site=F08$Locality, date=F08$date), FUN=sum)
SumM<-aggregate(F08$M.total, by=list(site=F08$Locality, date=F08$date), FUN=sum)
SumF
colnames(SumF)[3]<-"Females"
Sums<-cbind(SumF,Males=SumM[,3])
F08_s<-cbind(Sums,Location="Tvarminne")
str(F08_s)
F08_s

#Beach Seine data
BS<-read.delim("./sexID/beach_seines.txt", header=T)
BS$date<-paste(BS$Date,BS$Year, sep="-")
BS$date<-as.Date(BS$date,"%d-%b-%Y")
BS<-na.omit(BS)
BS<-BS[BS$Method=="Beach_seine",]
#Collapse the repeated observations into sums?
SumF<-aggregate(BS$Females, by=list(Site=BS$Site, date=BS$date), FUN=sum)
SumM<-aggregate(BS$Males, by=list(Site=BS$Site, date=BS$date), FUN=sum)
SumF
colnames(SumF)[3]<-"Females"
Sums<-cbind(SumF,Males=SumM[,3])
Sums
BS_s<-merge(Sums,BS[,c(1,2,4,9)], by=c("Site","date"))
BS_s<-unique(BS_s)
BS_s[BS_s$Site=="Kvalen_N" | BS_s$Site=="Kvalen_S" | BS_s$Site=="Kvalen_S_and_N",]$Site<-"Kvalen"
BS_s[BS_s$Site=="Langskar-Granbusken",]$Site<-"Granbusken"
BS_s

## Format for plot ----
#Setup some code to extract sex ratio data from several historical datasets and format it
Histo1<-data.frame(Males=Bergen2019$Males,Females=Bergen2019$Femals,Area="Bergen",Year="2019",Date=Bergen2019$date, Method="Beach Seine", Location=Bergen2019$Location) #Bergen 2019 data
Histo2<-data.frame(Males=KB7_s$Males, Females=KB7_s$Females,Area=KB7_s$location,Year="2007", Date=KB7_s$date, Method="Observation", Location=KB7_s$site) #Kberg and Bergen 2007 data
Histo3<-data.frame(Males=K2010_s$Males,Females=K2010_s$Females,Area=K2010_s$Location, Year="2010",Date=K2010_s$date, Method="Observation",Location=K2010_s$site) #Kberg 2010 data
Histo4<-data.frame(Males=F08_s$Males, Females=F08_s$Females, Area="Tvarminne",Year="2008", Date=F08_s$date, Method="Observation",Location=F08_s$site) #Findland 2008 data
Histo5<-data.frame(Males=K2008_s$Males,Females=K2008_s$Females,Area="Kristineberg", Year="2008", Date=K2008_s$date, Method="Observation",Location=K2008_s$site) #Kberg 2008
Histo6<-data.frame(Males=BS_s$Males,Females=BS_s$Females, Area=BS_s$Location, Year=BS_s$Year, Date=BS_s$date, Method="Beach Seine",Location=BS_s$Site)
Histo<-rbind(Histo1,Histo2, Histo3, Histo4, Histo5, Histo6)
str(Histo)

Histo$N<-Histo$Females+Histo$Males
Histo$Area_Year<-paste(Histo$Area,Histo$Year,sep = "-")
Histo$pm<-Histo$Males/Histo$N
Histo$date_short<-as.Date(format(Histo$Date,"%m/%d"),"%m/%d")
Histo$Area<-as.factor(Histo$Area)
Histo$Location<-as.factor(Histo$Location)
Histo$Area<- factor(Histo$Area, levels=c('Kristineberg', 'Bergen','Tvarminne'))
#Filter by sample size? or just don't calculate error bars if the normal approximation does not apply
Histo$check<-Histo$pm*Histo$N #check if Normal approximation is going to be ok. All values should be >5
Histo$CI<-NA
Histo[Histo$check>5,]$CI<-1.96*sqrt(Histo[Histo$check>5,]$pm*(1-Histo[Histo$check>5,]$pm)/Histo[Histo$check>5,]$N) #confidence interval using the normal approximation, justified by large sample size and high p
Histo$Year<-as.integer(Histo$Year)

str(Histo)
summary(Histo$N)

## Plot historical data ----
### FIGURE 8 ----
ggplot(Histo, aes(y=pm, x=date_short, fill=(Year), group=N))+
  geom_point(aes( size=N, shape=Method,fill=(Year), ),alpha=0.8,  position = position_jitterdodge(jitter.width=0.3, dodge.width=0.1))+
  geom_hline(yintercept=0.5, alpha=0.5)+
  geom_line(aes(group=factor(Year):Location),linewidth=1.5, alpha=0.2, color="orange")+
  geom_errorbar(width=0.4,alpha=0.3,aes(ymin=pm-CI,ymax=pm+CI),position = position_dodge(width=0.3))+
  facet_wrap(~Area)+
  scale_shape_manual("Sampling method",values=c(22,24))+
  scale_fill_gradient("Sampling year",low="darkblue", high="skyblue")+
  scale_size_continuous("Sample size", range = c(0.5,5), breaks=c(10,100,200,1000))+
  ylim(0,0.8)+
  theme_bw()+
  labs(y="Proportion of males", x="Time of the year")
ggsave("..//figures//Historical_ASR.pdf", width =9, height = 4)


### more exhaustive figure but a bit too large ----
ggplot(Histo, aes(y=pm, x=date_short, fill=Method, group=N))+
  geom_point(aes( size=N, shape=Method ),alpha=0.5,  position = position_jitterdodge(jitter.width=0.3, dodge.width=0.1))+
  geom_line(aes(group=Location),size=1, alpha=0.3)+
  geom_hline(yintercept=0.5, alpha=0.5)+
  geom_errorbar(width=0.4,alpha=0.3,aes(ymin=pm-CI,ymax=pm+CI),position = position_dodge(width=0.3))+
  facet_grid(rows=vars(Area),cols=vars(Year))+
  scale_shape_manual("Sampling method",values=c(22,24))+
  scale_fill_manual("Sampling method",values=c("orange","navyblue"))+
  scale_size_continuous("Sample size", range = c(0.5,5), breaks=c(10,100,200,1000))+
  ylim(0,1)+
  labs(y="Proportion of males", x="Location - year")


  