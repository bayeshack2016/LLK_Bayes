setwd("~/Desktop/Bayes")

######aggregate satellite (got from Alex, Star, Jeff)
sat<-read.csv("bayes_viirs_centiles.csv",header=T)
sat<-sat[sat$year=="2014",] ##only select year 2014
sat$id<-paste(sat$GEOID,sat$year,sat$month,sep='-') #create unique identifier

######read in mental disorder data
setwd("~/Desktop/Bayes/HHS")
mental<-read.delim("mental.csv",header=T,sep="\t")
mental$month<-as.integer(substr(mental$Month.Code,6,7))
mental$id<-paste(mental$County.Code,mental$Year,mental$month,sep='-') #create unique identifier

######merge two data by id
mergedata<-merge(mental,sat,by="id",all=F)


######further merge with SES data
edu<-read.csv("education.csv",header=T)
pov<-read.csv("PovertyEstimates.csv",header = T)
ump<-read.csv("Unemployment.csv",header = T)
pop<-read.csv("PopulationEstimates.csv",header=T)

edu<-edu[edu$FIPS.Code!=0,]
pov<-pov[pov$FIPStxt!=0,]

names(edu)[names(edu)=="FIPS.Code"] <- "FIPS"
names(ump)[names(ump)=="FIPS_Code"] <- "FIPS"
names(pov)[names(pov)=="FIPStxt"] <- "FIPS"

ses<-Reduce(function(x, y) merge(x, y, by="FIPS",all=F), list(edu, pop, ump,pov))
names(mergedata)[names(mergedata)=="GEOID"] <- "FIPS"
finaldata<-merge(mergedata,ses,by="FIPS",all.x=F)

######After preprocessing, we store cleaned data in finaldata.csv
data<-read.csv("finaldata.csv",header=T)


######run superlearner
llibrary(SuperLearner)

set.seed(1127)

index<-sample(seq_len(nrow(new)),size=floor(0.75*nrow(new)))

train<-new[index,] ##75% training

test<-new[-index,] ##25% test

SL.library<-c('SL.ridge','SL.rpartPrune','SL.polymars','SL.mean',
              'SL.randomForest','SL.gam','SL.glm','SL.glmnet',
              'SL.svm','SL.gbm','SL.nnet')

X<-train[c("month","pop","per_povall","median_income","ump_rate","R_DOMESTIC_MIG_2014",
           "R_INTERNATIONAL_MIG_2014","R_NATURAL_INC_2014","R_death_2014","per_bechelor",
           "per_college","per_highsch","per_lesshighsch","sum","lat","lng")]   

SL.out<-SuperLearner(Y=train$Deaths,X=X,SL.library = SL.library,family = 'gaussian',cvControl = list(V=5))

######model performance on test set
testX<-test[c("month","pop","per_povall","median_income","ump_rate","R_DOMESTIC_MIG_2014",
              "R_INTERNATIONAL_MIG_2014","R_NATURAL_INC_2014","R_death_2014","per_bechelor",
              "per_college","per_highsch","per_lesshighsch","sum","lat","lng")]   

pred<-predict.SuperLearner(SL.out,newdata=testX,onlySL = TRUE)

