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

Raw.SR <- read.delim("~Raw.SR.txt")

Raw.SR$Treat <- as.factor(Raw.SR$Treatment....females.)
Raw.SR$Dayspp <- as.factor(Raw.SR$Days.pp.1)

F.WL <- read.delim("~F.WL.txt")
M.WL <- read.delim("~M.WL.txt")


SR.dat <- Raw.SR  %>%
  group_by(Tube.tub) %>%
  summarise(m = sum(Males), f = sum(Females))

Treat.dat <- Raw.SR %>%
  group_by(Tube.tub) %>%
  summarise(Treatment = names(which.max(table(Treat))))

Time.dat <- Raw.SR %>%
  group_by(Tube.tub) %>%
  summarise(Time = names(which.max(table(Days.pp.1))))

All.SR <- SR.dat %>%
  full_join(Treat.dat, by = "Tube.tub") %>%
  full_join(Time.dat, by = "Tube.tub") %>%
  mutate_if(is.numeric, function(x) replace_na(x, -1))

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

  

#Questions
library(lme4)

# (1) Is there a size difference between males and females?

# (2) Is there an overall difference in size of wasps that emerge from high, medium 
#or low LMC treatments?

# (3) Is there an interaction between treatment and wasp sex?

# (4) Is wasp size correlated with emergence time?

# include random effect of tube ID 

