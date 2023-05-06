# LMC experiment 2023

#Lysiphlebus wing measurements

#Intraclass correlation coefficients for all 3 measures (wing length, A-B and A-C)
#testing repeatability

#Based on Koo & Li 2016
#ICC < 0.5 = poor
#ICC 0.5-0.75 = moderate
#ICC 0.75-0.9 = good
#ICC > 0.9 = excellent

library(irr)
library(tidyr)

#read in data
wingL <- read.delim("IOR-ICC.txt")

LengthICC <- cbind(wingL$Wing.length, wingL$Wing.length.1)
LengthICC <- as.data.frame(LengthICC)
LengthICC <- LengthICC %>% drop_na()
LengthICC$V1 <- as.numeric(LengthICC$V1)
LengthICC$V2 <- as.numeric(LengthICC$V2)
View(LengthICC)
icc(LengthICC, model="twoway", type="agreement")

#Single Score Intraclass Correlation

#Subjects = 67 
#Raters = 2 
#ICC(A,1) = 0.915 (= excellent based on Koo & Li 2016)

#F-Test, H0: r0 = 0 ; H1: r0 > 0 
#F(66,54.1) = 24 , p = 1.55e-24 

#95%-Confidence Interval for ICC Population Values:
#  0.861 < ICC < 0.948
ABICC <- cbind(wingL$A.B.length, wingL$A.B.length.1)
ABICC <- as.data.frame(ABICC)
ABICC <- ABICC %>% drop_na()
ABICC$V1 <- as.numeric(ABICC$V1)
ABICC$V2 <- as.numeric(ABICC$V2)
View(ABICC)
icc(ABICC, model="twoway", type="agreement")

#Single Score Intraclass Correlation

#Subjects = 66 
#Raters = 2 
#ICC(A,1) = 0.902 = excellent

#F-Test, H0: r0 = 0 ; H1: r0 > 0 
#F(65,65.5) = 19.7 , p = 4.72e-26 

#95%-Confidence Interval for ICC Population Values:
#  0.845 < ICC < 0.939

BCICC <- cbind(wingL$B.C.length, wingL$B.C.length.1)
BCICC <- as.data.frame(BCICC)
BCICC <- BCICC %>% drop_na()
BCICC$V1 <- as.numeric(BCICC$V1)
BCICC$V2 <- as.numeric(BCICC$V2)
View(BCICC)
icc(BCICC, model="twoway", type="agreement")

#Single Score Intraclass Correlation

#Subjects = 66 
#Raters = 2 
#ICC(A,1) = 0.836 = good

#F-Test, H0: r0 = 0 ; H1: r0 > 0 
#F(65,65) = 11 , p = 1.01e-18 

#95%-Confidence Interval for ICC Population Values:
#  0.745 < ICC < 0.896

#check - any difference between left and right wing 


#------------------------------------------------------------------------------

#data wrangling - need to combine sex ratio data and wing length and amalgamate any duplicates

Raw.SR <- read.delim("Raw.SR.txt")

Raw.SR$Treat <- as.factor(Raw.SR$Treatment....females.)
Raw.SR$Dayspp <- as.factor(Raw.SR$Days.pp.1)
Raw.SR$Tot <- Raw.SR$Males + Raw.SR$Females

F.WL <- read.delim("F.WL.txt")
M.WL <- read.delim("M.WL.txt")

detach("package:plyr", unload = TRUE) #need to do this or the next bit doesn't work
library(dplyr)

SR.dat <- Raw.SR  %>%
  group_by(Tube.tub) %>%
  summarise(m = sum(Males), f = sum(Females), tot = sum(Tot))

Treat.dat <- Raw.SR %>%
  group_by(Tube.tub) %>%
  summarise(Treatment = names(which.max(table(Treat))))

Time.dat <- Raw.SR %>%
  group_by(Tube.tub) %>%
  summarise(Time = names(which.max(table(Days.pp.1))))

Total.tubes <- Raw.SR %>%
  group_by(Tube.tub) %>%
  summarise(Tube_ID = names(which.max(table(Tube_number))))

Total.wasps <- Raw.SR %>%
  group_by(Tube_number) %>%
  summarise(Total = sum(Tot))

colnames(Total.wasps)<- c('Tube_ID', 'Total')

All.SR <- SR.dat %>%
  full_join(Treat.dat, by = "Tube.tub") %>%
  full_join(Time.dat, by = "Tube.tub") %>%
  full_join(Total.tubes, by = "Tube.tub") %>%
  mutate_if(is.numeric, function(x) replace_na(x, -1))

#need to summarise totals and overall sex ratios for each replicate 

Fem.WL <- F.WL %>%
  group_by(Tube.tub) %>%
  summarise(Winglength = mean(Wing.length))

Male.WL <- M.WL %>%
  group_by(Tube.tub) %>%
  summarise(Winglength = mean(Wing.length))

combine <- Male.WL %>%
  full_join(Fem.WL, by = "Tube.tub") %>%
  mutate_if(is.numeric, function(x) replace_na(x, -1))

SR.WL <- All.SR %>%
  full_join(combine, by = "Tube.tub") %>%
  mutate_if(is.numeric, function(x) replace_na(x, -1))

SR.WL <- SR.WL[complete.cases(SR.WL), ]

SR.WL$M.Winglength <- na_if(SR.WL$Winglength.x, -1)
SR.WL$F.Winglength <- na_if(SR.WL$Winglength.y, -1)

SR.WL <- subset(SR.WL, m>=0 & f>=0)

Males <- SR.WL[, c('Tube.tub', 'm', 'f', 'Treatment', 'Time', 'M.Winglength')]
colnames(Males)<- c('Tube.tub', 'm', 'f', 'Treatment', 'Time', 'Winglength')

Females <- SR.WL[, c('Tube.tub', 'm', 'f', 'Treatment', 'Time', 'F.Winglength')]
colnames(Females)<- c('Tube.tub', 'm', 'f', 'Treatment', 'Time', 'Winglength')


Males$sex <- 1
Females$sex <- 0

library(data.table)

M.F.WL <- bind_rows(Males,Females)

#Questions
library(lme4)
library(glmmTMB)
library(DHARMa)
library(car)
library(binom)
library(ggplot2)

###### Sex ratio

# (1) Is sex ratio consistent across treatment and time?

#remove tube 13 (only 1 male wasp emerged - mix up with 18 so can't use)

SR.WL <- SR.WL[-c(445), ]

SR.WL$Timecont <- as.numeric(SR.WL$Time)
SexRatioMod <- glmer(cbind(m,f) ~ Treatment*Timecont + (1|Tube_ID), family =  binomial(), data=SR.WL)
SimulationOutputsr <- simulateResiduals(SexRatioMod, plot = T)

summary(SexRatioMod)
car::Anova(SexRatioMod)

#Only effect is time - sex ratio becomes female biased over time
males <- SR.WL %>%
  group_by(Treatment) %>%
  summarise( 
    sum=sum(m))
total<- SR.WL %>%
  group_by(Treatment) %>%
  summarise( 
    sum=sum(tot))

library(binom)
total$binomial_CI <- binom.confint(males$sum, total$sum, conf.level = 0.95, methods = "logit")

mean(total$binomial_CI$mean)

mean(total$binomial_CI$mean - total$binomial_CI$lower)


#plot means
SR.WL$SR <- SR.WL$m/SR.WL$tot
SR.WL$Timecont <- as.numeric(SR.WL$Time)

level_order <- c('1', '5', '10') 
cols <- c("1" = "red", "5" = "yellow", "10" = "blue")
Time.Treat.SR <- ggplot(SR.WL, aes(Timecont, SR, colour=factor(Treatment, level = level_order))) +
  geom_point() + geom_smooth(method = "lm") +
  xlab("Emergence time (days") + ylab("Emergence sex ratio (proportion males)") + labs(colour = "Treatment") + scale_colour_manual(values = cols)
Time.Treat.SR <- Time.Treat.SR + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                               size="12", color="black"), legend.position="top")

# repeat using cumulative sums of males and females

SR.WL$csm <- ave(SR.WL$m, SR.WL$Tube_ID, FUN=cumsum)
SR.WL$csf <- ave(SR.WL$f, SR.WL$Tube_ID, FUN=cumsum)
SR.WL$cstot <- ave(SR.WL$tot, SR.WL$Tube_ID, FUN=cumsum)


cummales <- SR.WL %>%
  group_by(Treatment) %>%
  summarise( 
    sum=sum(csm))
cumtotal<- SR.WL %>%
  group_by(Treatment) %>%
  summarise( 
    sum=sum(cstot))

library(binom)
cumtotal$binomial_CI <- binom.confint(cummales$sum, cumtotal$sum, conf.level = 0.95, methods = "logit")

mean(cumtotal$binomial_CI$mean)
mean(cumtotal$binomial_CI$mean - cumtotal$binomial_CI$lower)


CumSexRatioMod <- glmer(cbind(csm,csf) ~ Treatment + Timecont + (1|Tube_ID), family =  binomial(), data=SR.WL)
SimulationOutputsr <- simulateResiduals(CumSexRatioMod, plot = T)

summary(CumSexRatioMod)
car::Anova(CumSexRatioMod)

#Only effect still time

#plot means
SR.WL$cumSR <- SR.WL$csm/SR.WL$cstot
SR.WL$Timecont <- as.numeric(SR.WL$Time)


level_order <- c('1', '5', '10') 
cols <- c("1" = "red", "5" = "yellow", "10" = "blue")
Time.Treat.cumSR <- ggplot(SR.WL, aes(Timecont, cumSR, colour=factor(Treatment, level = level_order))) +
  geom_point() + geom_smooth(method = "lm") +
  xlab("Emergence time (days") + ylab("Cumulative sex ratio (proportion males)") + labs(colour = "Treatment") + scale_colour_manual(values = cols)
Time.Treat.cumSR <- Time.Treat.cumSR + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                                 size="12", color="black"), legend.position="top")
# same - just more female biased over time and no differences between groups

# (2) Is the overall sex ratio related to total number of wasps that emerged (in each rep)?

Sum.SR.dat <- SR.WL  %>%
  group_by(Tube_ID) %>%
  summarise(m = sum(m), f = sum(f), tot = sum(tot))

Sum.Treat.dat <- SR.WL %>%
  group_by(Tube_ID) %>%
  summarise(Treatment = names(which.max(table(Treatment))))

Sum.SR <- Sum.SR.dat %>%
  full_join(Sum.Treat.dat, by = "Tube_ID") 

SumSexRatioMod <- glm(cbind(m,f) ~ Treatment*tot, family =  binomial(), data=Sum.SR)
SimulationOutputsr <- simulateResiduals(SumSexRatioMod, plot = T)

summary(SumSexRatioMod)
car::Anova(SumSexRatioMod)

#no effect of offspring number or treatment (or interaction) on overall sex ratio

males<- Sum.SR %>%
  group_by(Treatment) %>%
  summarise( 
    sum=sum(m))
tot.wasps <- Sum.SR %>%
  group_by(Treatment) %>%
  summarise( 
    sum=sum(tot))

level_order <- c('1', '5', '10') 
cols <- c("1" = "red", "5" = "yellow", "10" = "blue")
Tot.Treat.SR <- ggplot(SR.WL, aes(tot, SR, colour=factor(Treatment, level = level_order))) +
  geom_point() + geom_smooth(method = "lm") +
  xlab("Total wasps emerging") + ylab("Sex ratio") + labs(colour = "Treatment") + scale_colour_manual(values = cols)
Tot.Treat.SR <- Tot.Treat.SR + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                             size="12", color="black"), legend.position="top")

# (3) is there an effect of treatment on the number of wasps that emerge?

TotMod <- glm(tot ~ Treatment, data=Sum.SR)
SimulationOutputsr <- simulateResiduals(TotMod, plot = T)
summary(TotMod)
car::Anova(TotMod)

#sig diff between 1 and 5 and 1 and 10
#test whether 5 and 10 differ

library(plyr)
Sum.SR$Treat.recode <- revalue(Sum.SR$Treatment, c("1" = "10", "10" = "1"))
Sum.SR$Treat.recode <- as.character(Sum.SR$Treat.recode)
ReTotMod <- glm(tot ~ Treat.recode, data=Sum.SR)
SimulationOutputsr <- simulateResiduals(ReTotMod, plot = T)
summary(ReTotMod)
car::Anova(ReTotMod)
detach("package:plyr", unload = TRUE)

#no difference between 10 and 5 (p = 0.7)


Totalemerge <- Sum.SR %>%
  group_by(Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(tot),
    sd=sd(tot),
    max=max(tot),
    min=min(tot),
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

Totalemerge$mean[[1]]
Totalemerge$ic[[1]]

Tot.Treat.fig <-  ggplot(Totalemerge, aes(x=factor(Treatment, level = level_order), y=mean, fill = Treatment)) + 
  geom_bar(position=position_dodge(), stat="identity", fill=cols, colour="black") +
  geom_errorbar( aes(ymin=mean-ic, ymax=mean+ic), width=0.2) +
  xlab("Treatment") + ylab("total wasps emerged") 
Tot.Treat.fig <- Tot.Treat.fig + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                                   size="12", color="black"), legend.position="top")

#no difference between 5 and 10 in terms of mean number of wasps - 
# didn't count number of nymphs, possibly not enough? Or something else - saturation point?
# Max for 1 female is 94, on average 58; max for 10 females is 205, avg is 98
# max for 5 females is 179, avg is 193 


# (4) is the sequence of emergence the same across treatments?

#GAMM for tot - do males and females show different emergence patterns?

Females <- SR.WL[, c('Tube.tub', 'f', 'Tube_ID', 'Treatment', 'Timecont')]
colnames(Females)<- c('Tube.tub', 'n', 'Tube_ID', 'Treatment', 'Time')
Females$Sex <- "F"


Males <- SR.WL[, c('Tube.tub', 'm', 'Tube_ID', 'Treatment', 'Timecont')]
colnames(Males)<- c('Tube.tub', 'n', 'Tube_ID', 'Treatment', 'Time')
Males$Sex <- "M"

Long.SR.WL <- rbind(Females, Males)

library(mgcv)

LMC.emergence <- gamm(n ~ s(Time, fx = FALSE, k=-1, bs = "ds"), random = list(Tube_ID=~1), 
              , family = poisson(link = "log"), data = Long.SR.WL)
Long.SR.WL$Sex <- as.factor(Long.SR.WL$Sex)

LMC.sex.emergence <- gamm(n ~ s(Time, fx = FALSE, k=-1, bs = "ds", by = Sex), random = list(Tube_ID=~1), 
                          , family = poisson(link = "log"), method = "REML", data = Long.SR.WL)

Long.SR.WL$Treatment <- as.factor(Long.SR.WL$Treatment)

LMC.Treat.emergence <- gamm(n ~ s(Time, fx = FALSE, k=-1, bs = "ds", by = Treatment), random = list(Tube_ID=~1), 
                          , family = poisson(link = "log"), method = "REML", data = Long.SR.WL)


gam.check(LMC.emergence$gam)
gam.check(LMC.sex.emergence$gam)
gam.check(LMC.Treat.emergence$gam)

anova(LMC.emergence$lme, LMC.sex.emergence$lme, LMC.Treat.emergence$lme)
#improve when you put sex in and when you put treatment in

summary(LMC.emergence$gam)
anova(LMC.emergence$gam)
plot(LMC.emergence$gam,pages=1)

summary(LMC.sex.emergence$gam)
anova(LMC.sex.emergence$gam)
plot(LMC.sex.emergence$gam,pages=1)

summary(LMC.Treat.emergence$gam)
anova(LMC.Treat.emergence$gam)
plot(LMC.Treat.emergence$gam,pages=1)

#Treatment by sex

#Females

LMC.Females.emergence <- gamm(n ~ s(Time, fx = FALSE, k=-1, bs = "ds"), random = list(Tube_ID=~1), 
                            family = poisson(link = "log"), method = "REML", data = Females)

Females$Treatment <- as.factor(Females$Treatment)

LMC.Females.Treat.emergence <- gamm(n ~ s(Time, fx = FALSE, k=-1, bs = "ds", by = Treatment), random = list(Tube_ID=~1), 
                              family = poisson(link = "log"), method = "REML", data = Females)

summary(LMC.Females.emergence$gam)
anova(LMC.Females.Treat.emergence$gam)
plot(LMC.Females.emergence$gam,pages=1)
plot(LMC.Females.Treat.emergence$gam,pages=1)


#adding in Treatment reduces AIC 

#Males

Males$Treatment <- as.factor(Males$Treatment)

LMC.Males.emergence <- gamm(n ~ s(Time, fx = FALSE, k=-1, bs = "ds"), random = list(Tube_ID=~1), 
                              family = poisson(link = "log"), method = "REML", data = Males)

LMC.Males.Treat.emergence <- gamm(n ~ s(Time, fx = FALSE, k=-1, bs = "ds", by = Treatment), random = list(Tube_ID=~1), 
                                   family = poisson(link = "log"), method = "REML", data = Males)

summary(LMC.Males.emergence$gam)
anova(LMC.Males.Treat.emergence$gam)
plot(LMC.Males.emergence$gam,pages=1)
plot(LMC.Males.Treat.emergence$gam,pages=1)

# what does this actually tell me? Emergence patterns are different but how quantitatively?
#SSM not necessary as we have all the data (i.e. no hidden structure that needs estimating)

# do males emerge earlier than females?
  #does this differ over treatments

# recode data

data <- Females
stretch_females <- 
  data.frame(data[rep(seq_len(dim(data)[1]), data$n), c(1,3:6), drop = FALSE], row.names=NULL) 
stretch_females$n <- 1

data <- Males
stretch_males <- 
  data.frame(data[rep(seq_len(dim(data)[1]), data$n), c(1,3:6), drop = FALSE], row.names=NULL) 
stretch_males$n <- 1                            

stretch_data <- rbind(stretch_males, stretch_females)

#use cox model with mixed effects

library(coxme)
Time.mod.Treat.Sex.int<- coxme(Surv(Time) ~ Treatment*Sex + (1|Tube_ID), stretch_data)
Time.mod.Treat.Sex<- coxme(Surv(Time) ~ Treatment +Sex + (1|Tube_ID), stretch_data)
Time.mod.Treat <- coxme(Surv(Time) ~ Treatment + (1|Tube_ID), stretch_data)
Time.mod.Sex <- coxme(Surv(Time) ~ Sex + (1|Tube_ID), stretch_data)
Time.mod.Null<- coxme(Surv(Time) ~ 1+ (1|Tube_ID), stretch_data)


anova(Time.mod.Null, Time.mod.Sex, Time.mod.Treat, Time.mod.Treat.Sex, Time.mod.Treat.Sex.int)
#lowest LL includes treat and sex (interaction doesn't make a difference)

summary(Time.mod.Treat.Sex)
car::Anova(Time.mod.Treat.Sex)

#No interaction effect but both treatment and sex are related to emergence time
#5 and 10 sig diff than 1, recode - are they different from each other?

library(plyr)
stretch_data$Treatment.recode <- revalue(stretch_data$Treatment, c("1" = "10", "10" = "1"))
stretch_data$Treatment.recode <- as.character(stretch_data$Treatment.recode)
detach("package:plyr", unload = TRUE)

Time.mod.Treat.Sex.recode <- coxme(Surv(Time) ~ Treatment.recode +
         Sex + (1|Tube_ID), stretch_data)
                                    
summary(Time.mod.Treat.Sex.recode)
car::Anova(Time.mod.Treat.Sex.recode)

#no sig diff in emergence time of both sexes between 5 and 10

#check with female and male data only

Time.mod.Treat.fem <- coxme(Surv(Time) ~ Treatment + (1|Tube_ID), stretch_females)
summary(Time.mod.Treat.fem)
car::Anova(Time.mod.Treat.fem)

library(plyr)
stretch_females$Treatment.recode <- revalue(stretch_females$Treatment, c( "1" = "10", "10" = "1"))
stretch_females$Treatment.recode <- as.character(stretch_females$Treatment.recode)

Time.mod.Treat.Fem.recode <- coxme(Surv(Time) ~ Treatment.recode +
                                     (1|Tube_ID), stretch_females)
summary(Time.mod.Treat.Fem.recode)
car::Anova(Time.mod.Treat.Fem.recode)

#for females - emergence date sig for 1 foundress, 5 and 10 the same

Time.mod.Treat.male<- coxme(Surv(Time) ~ Treatment + (1|Tube_ID), stretch_males)
summary(Time.mod.Treat.male)
car::Anova(Time.mod.Treat.male)

stretch_males$Treatment.recode <- revalue(stretch_males$Treatment, c("1" = "10", "10" = "1"))
stretch_males$Treatment.recode <- as.character(stretch_males$Treatment.recode)
detach("package:plyr", unload = TRUE)

Time.mod.Treat.male.recode <- coxme(Surv(Time) ~ Treatment.recode +
                                     (1|Tube_ID), stretch_males)
summary(Time.mod.Treat.male.recode)
car::Anova(Time.mod.Treat.male.recode)

#for males only 1 and 10 sig differ (1 and 5 and 10 and 5 the same) 
#average earlier emergence when only 1 foundress

#plot emergence of males and females over time (all treatments)

Sex.emergence.fig <- ggplot(data = Long.SR.WL, aes(x = Time, y = n, colour=Sex)) +
  geom_smooth(method = "loess", alpha=0.2)
Sex.emergence.fig <- Sex.emergence.fig + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                                         size="12", color="black"), legend.position="top")
Sex.emergence.fig + labs(colour = "Treatment") + xlab("Emergence time (days)") + ylab("# emerging wasps")


#plot emergence of females in different treatments

level_order <- c('1', '5', '10') 
Fem.emergence.fig <- ggplot(data = Females, aes(x = Time, y = n, colour=factor(Treatment, level = level_order))) +
  geom_smooth(method = "loess", alpha=0.2) + scale_y_continuous(limits = c(0, NA))
Fem.emergence.fig <- Fem.emergence.fig + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                                   size="12", color="black"), legend.position="top")
Fem.emergence.fig + labs(colour = "Treatment") + xlab("Emergence time (days)") + ylab("# emerging female wasps")

Male.emergence.fig <- ggplot(data = Males, aes(x = Time, y = n, colour=factor(Treatment, level = level_order))) +
  geom_smooth(method = "loess", alpha=0.2) + scale_y_continuous(limits = c(0, NA))
Male.emergence.fig <- Male.emergence.fig + theme_classic() + theme(axis.title=element_text(face="bold.italic",
size="12", color="black"), legend.position="top")
Male.emergence.fig + labs(colour = "Treatment") + xlab("Emergence time (days)") + ylab("# emerging male wasps")


#bar plot of means for all treatments and by sex

Long.emerge.sum <-  stretch_data %>%
  group_by(Treatment, Sex) %>%
summarise(
  n=n(),
  mean=mean(Time),
  sd=sd(Time),
  max=max(Time),
  min=min(Time),
) %>%
  mutate( se=sd/sqrt(n)) %>%
  mutate(  ic=se * qt((1-0.05)/2 + .5, n-1))

  
level_order <- c('1', '5', '10') 
Full.emergence.fig <- ggplot(data = Long.emerge.sum,
       aes(x = Sex, y = mean, fill = factor(Treatment, level = level_order))) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-ic, ymax=mean+ic), width=.2,
                position=position_dodge(.9)) 
Full.emergence.fig <- Full.emergence.fig + coord_cartesian(ylim=c(12,13.5))
Full.emergence.fig <- Full.emergence.fig + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                                       size="12", color="black"), legend.position="top")
Full.emergence.fig  + ylab("Emergence time (days)") + xlab("Sex") + labs(fill= "Treatment") + scale_colour_manual(values = cols)

female.emerge.mean <- mean(Long.emerge.sum$mean[[1]],Long.emerge.sum$mean[[3]],Long.emerge.sum$mean[[5]])
female.emerge.ic <- mean(Long.emerge.sum$ic[[1]],Long.emerge.sum$ic[[3]],Long.emerge.sum$ic[[5]])

male.emerge.mean <- mean(Long.emerge.sum$mean[[2]],Long.emerge.sum$mean[[4]],Long.emerge.sum$mean[[6]])
male.emerge.ic <- mean(Long.emerge.sum$ic[[2]],Long.emerge.sum$ic[[4]],Long.emerge.sum$ic[[6]])

oneemerge.mean <- mean(Long.emerge.sum$mean[[1]],Long.emerge.sum$mean[[2]])
oneemerge.ic <- mean(Long.emerge.sum$ic[[1]],Long.emerge.sum$ic[[2]])

fiveemerge.mean <- mean(Long.emerge.sum$mean[[5]],Long.emerge.sum$mean[[6]])
fiveemerge.ic <- mean(Long.emerge.sum$ic[[5]],Long.emerge.sum$ic[[6]])

tenemerge.mean <- mean(Long.emerge.sum$mean[[3]],Long.emerge.sum$mean[[4]])
tenemerge.ic <- mean(Long.emerge.sum$ic[[3]],Long.emerge.sum$ic[[4]])

#Earlier average emergence for females and males when only 1 foundress

#test using glm

# WING LENGTH

Males <- SR.WL[, c('Tube.tub', 'Tube_ID','Treatment', 'Time', 'M.Winglength')]
colnames(Males)<- c('Tube.tub', 'Tube_ID', 'Treatment', 'Time', 'M.Winglength')

Females <- SR.WL[, c('Tube.tub', 'Tube_ID', 'Treatment', 'Time', 'F.Winglength')]
colnames(Females)<- c('Tube.tub', 'Tube_ID', 'Treatment', 'Time', 'F.Winglength')


Males$sex <- 1
Females$sex <- 0

library(data.table)

M.F.WL <- bind_rows(Males,Females)
M.F.WL$Timecont<- as.numeric(M.F.WL$Time)

M.F.WL.clean <-
M.F.WL %>% drop_na(Winglength)


Wing.length.mod <- lmer(Winglength ~ sex*Treatment + Treatment*Timecont
                        + Timecont*sex + Treatment*Timecont*sex +
                          (1|Tube_ID), data = M.F.WL.clean)

SimulationOutputsr <- simulateResiduals(Wing.length.mod, plot = T)

summary(Wing.length.mod)
car::Anova(Wing.length.mod)

library(plyr)
M.F.WL.clean$sex <- as.character(M.F.WL.clean$sex)
M.F.WL.clean$Sex <- revalue(M.F.WL.clean$sex, c( 
                           "0" = "F", "1" = "M"))
detach("package:plyr", unload = TRUE)


#(1) Is there a size difference between males and females?

#yes - males bigger


# (2) Is there an overall difference in size of wasps that emerge from high, medium 
#or low LMC treatments?

#no 


# (3) Is there an interaction between treatment and wasp sex?

#no

#plotting treatment and sex

#

level_order <- c('1', '5', '10') 
Full.WL.fig <- ggplot(data = WL.sum,
                             aes(x = Sex, y = mean, fill = factor(Treatment, level = level_order))) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-ic, ymax=mean+ic), width=.2,
                position=position_dodge(.9)) 
Full.WL.fig <- Full.WL.fig + coord_cartesian(ylim=c(1500,1700))
Full.WL.fig <- Full.WL.fig + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                                           size="12", color="black"), legend.position="top")
Full.WL.fig  + ylab("Winglength (µM)") + xlab("Sex") + labs(fill= "Treatment") + scale_colour_manual(values = cols)



# (4) Is wasp size correlated with emergence time?

#yes - smaller over time - males and females


Time.WL.fig <- ggplot(data = M.F.WL.clean, aes(x = Timecont,  y = Winglength, colour=Sex)) +
  geom_smooth(method = "loess", alpha=0.2)
Time.WL.fig <- Time.WL.fig + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                                        size="12", color="black"), legend.position="top")
Time.WL.fig + labs(colour = "Sex") + xlab("Emergence time (days)") + ylab("Winglength (µM)")


#what about sexual dimorphism?
#male size-female size

Males <- SR.WL[, c('Treatment', 'Time', 'M.Winglength')]
colnames(Males)<- c('Treatment', 'Time', 'M.Winglength')

Females <- SR.WL[, c('Treatment', 'Time', 'F.Winglength')]
colnames(Males)<- c('Treatment', 'Time', 'F.Winglength')


#sexual dimorphism in wing length doesn't change over time or with treatment

SD.WL <- SR.WL %>%
  group_by(Time, Treatment) %>%
  summarise(SD = M.Winglength-F.Winglength)

SD.WL <- SD.WL[complete.cases(SD.WL), ]

Mean.SD.WL <- SD.WL %>%
  group_by(Time, Treatment) %>% 
  summarise(meanSD=mean(SD))

Mean.SD.WL$Timecont <- as.numeric(Mean.SD.WL$Time)

SD.WL.mod <- glm(meanSD ~ Timecont*Treatment, data = Mean.SD.WL)
hist(SD.WL.mod$residuals) #looks ok

summary(SD.WL.mod)
car::Anova(SD.WL.mod)

#no change in SD over time or with treatment

level_order <- c('1', '5', '10') 
Time.WL.SD.fig <- ggplot(data = Mean.SD.WL, aes(x = Timecont,  y = meanSD, colour=factor(Treatment, level = level_order))) +
  geom_smooth(method = "loess", alpha=0.2)
Time.WL.SD.fig <- Time.WL.SD.fig + theme_classic() + theme(axis.title=element_text(face="bold.italic", 
                                                                             size="12", color="black"), legend.position="top")
Time.WL.SD.fig + labs(colour = "Treatment") + xlab("Emergence time (days)") + ylab("Sexual dimorphism in Winglength (µM)")



#does male wing size correlate with sex ratio in the replicate?
#if being big is an advantage expect bigger males in more male biased conditions?
#could also be the opposite

SR.WL.comp <- SR.WL[complete.cases(SR.WL), ]

SR.WL.by.tube_ID <- SR.WL.comp  %>%
  group_by(Tube_ID, Treatment) %>%
  summarise(m = sum(m), f = sum(f), tot = sum(tot), M.WL = mean(M.Winglength), F.WL = mean(F.Winglength))

SR.WL.by.tube_ID$SR <- SR.WL.by.tube_ID$m/SR.WL.by.tube_ID$tot

M.SR.WL.null <- glm(M.WL ~ 1, family = gaussian(link = "identity"),
                   data = SR.WL.by.tube_ID)
M.SR.WL.SR <- glm(M.WL ~ SR, family = gaussian(link = "identity"),
                    data = SR.WL.by.tube_ID)
M.SR.WL.tot <- glm(M.WL ~ tot, family = gaussian(link = "identity"),
                  data = SR.WL.by.tube_ID)
M.SR.WL.m<- glm(M.WL ~ m, family = gaussian(link = "identity"),
                  data = SR.WL.by.tube_ID)
M.SR.WL.f <- glm(M.WL ~ f, family = gaussian(link = "identity"),
                  data = SR.WL.by.tube_ID)
#not enough df to estimate effects of male number, female number, total or sex ratio
# pick best model

AIC(M.SR.WL.null, M.SR.WL.SR, M.SR.WL.tot, M.SR.WL.m, M.SR.WL.f)
#lowest AIC is null - no effect

#female wing length

F.SR.WL.null <- glm(F.WL ~ 1, family = gaussian(link = "identity"),
                    data = SR.WL.by.tube_ID)
F.SR.WL.SR <- glm(F.WL ~ SR, family = gaussian(link = "identity"),
                  data = SR.WL.by.tube_ID)
F.SR.WL.tot <- glm(F.WL ~ tot, family = gaussian(link = "identity"),
                   data = SR.WL.by.tube_ID)
F.SR.WL.m<- glm(F.WL ~ m, family = gaussian(link = "identity"),
                data = SR.WL.by.tube_ID)
F.SR.WL.f <- glm(F.WL ~ f, family = gaussian(link = "identity"),
                 data = SR.WL.by.tube_ID)

AIC(F.SR.WL.null, F.SR.WL.SR, F.SR.WL.tot, F.SR.WL.m, F.SR.WL.f)

#null is best model

SR.WL.by.tube_ID$SD <- SR.WL.by.tube_ID$M.WL - SR.WL.by.tube_ID$F.WL

SD.SR.WL.null <- glm(SD ~ 1, family = gaussian(link = "identity"),
                    data = SR.WL.by.tube_ID)
SD.SR.WL.SR <- glm(SD ~ SR, family = gaussian(link = "identity"),
                  data = SR.WL.by.tube_ID)
SD.SR.WL.tot <- glm(SD ~ tot, family = gaussian(link = "identity"),
                   data = SR.WL.by.tube_ID)
SD.SR.WL.m<- glm(SD ~ m, family = gaussian(link = "identity"),
                data = SR.WL.by.tube_ID)
SD.SR.WL.f <- glm(SD ~ f, family = gaussian(link = "identity"),
                 data = SR.WL.by.tube_ID)

AIC(SD.SR.WL.null, SD.SR.WL.SR, SD.SR.WL.tot, SD.SR.WL.m, SD.SR.WL.f)

### Summary of results

#Emergence sex ratio
 # Becomes more female biased over time. No effect of number of foundresses.

#Cumulative sex ratio 
   # same pattern as emergence sex ratio. 

#No effect of the total number of wasps that emerged or the number of foundresses 
#on the sex ratio

#more wasps in total emerge when 5 or 10 foundresses compared to 1, but no diff between 5 and 10
  #suggests saturation point reached

# For emergence time - GAMMs suggest adding sex and treatment to emergence time data both improve 
   #model fit 
#can't get sex and treatment in the model - but separate univariate models show differences 
#emergence patterns for males and females - females increase in numbers more rapidly
#at the start compared to males but no obvious qualitative patterns 
  #males also no obvious patterns

#also used a Cox model with random effects (coxme) to test emergence sequence - 
#best model included treatment and sex b

#plotted figure of average females across all groups - males emerged eariler, 
#individuals in treatment 1 emerged earlier than 5 and 10, no diff between 5 and 10
#males emerged earlier than females

#winglength - individuals that emerged earlier are larger than later emergers. 
#Males are larger than females but no effect of treatment on wing length
#sexual dimorphism in winglength is the same over time and across treatments
#winglength not related to the sex ratio or number offspring - not a sex allocation
#strategy to optimise fitness by producing larger males when more competition


#useful tidyverse tutorial
#https://www.youtube.com/watch?v=sVISY_27znA&list=PLLxj8fULvXwGOf8uHlL4Tr62oXSB5k_in&index=7



