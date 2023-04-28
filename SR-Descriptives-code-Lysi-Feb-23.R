#This bit of code calculates cumulative sums of males and females over time
#across 3 treatments in the LMC experiment. 
#No stats, just numbers so far

library(dplyr)
library(ggplot2)
library(Rmisc)
library(lme4)
library(Hmisc)


LMC.Lysi <- read.delim("~LMC-clean-15.2.23.txt")
LMC.Lysi$Treatment <- as.character(LMC.Lysi$Treatment....females.)
LMC.Lysi <- LMC.Lysi[!is.na(LMC.Lysi$Treatment),]

#get cumulative sums of males and females per tube over time

cumfem <- aggregate(Females ~ Tube_number, LMC.Lysi, cumsum)
                                   
LMC_cumul_sums <- 
  LMC.Lysi %>% 
  arrange(as.numeric(Hours.post.parasitism)) %>%
  group_by(Tube_number) %>%
  mutate(cumul_males = cumsum(Males),
         cumul_females = cumsum(Females)) 


#get cumulative means of males and females per treatment over time

males <- group.CI(cumul_males ~ Treatment*Hours.post.parasitism, data = LMC_cumul_sums, ci = 0.95)
females <- group.CI(cumul_females ~ Treatment*Hours.post.parasitism, data = LMC_cumul_sums, ci = 0.95)

males_females_list <- list(males, females)
males_females <- Reduce(function(x,y) merge(x, y, all=TRUE), males_females_list)

write.table(males_females, "~/Documents/Work in progress/Lysiphlebus /LMC/Data and code/males_females.txt")

#I had to modify this in excel to add males and females to the same column, probably a way to do it in R though
cumul_means <- read.delim("~/Documents/Work in progress/Lysiphlebus /LMC/Data and code/cumul_means.txt")
cumul_means <- cumul_means[!is.na(cumul_means$Treatment),]

#use this to make some nice graphs (subset by treatment or it gets a bit messy)
T1.means <- subset(cumul_means, Treatment == 1)
T1.males.females.plot <- ggplot(T1.means, aes(x = Hours.post.parasitism, y = cumul_males.mean, color = Sex)) +
  geom_point() + geom_smooth() + aes(ymin = 0, ymax = 70) +  xlab("Time (days post parasitism)") + ylab("Number of individuals emerged")
T1.males.females.plot

T5.means <- subset(cumul_means, Treatment == 5)
T5.males.females.plot <- ggplot(T5.means, aes(x = Hours.post.parasitism, y = cumul_males.mean, color = Sex)) +
  geom_point() + geom_smooth() + aes(ymin = 0, ymax = 70)  +  xlab("Time (days post parasitism)") + ylab("Number of individuals emerged")
T5.males.females.plot

T10.means <- subset(cumul_means, Treatment == 10)
T10.males.females.plot <- ggplot(T10.means, aes(x = Hours.post.parasitism, y = cumul_males.mean, color = Sex)) +
  geom_point() + geom_smooth() + aes(ymin = 0, ymax = 70)  +  xlab("Time (days post parasitism)") + ylab("Number of individuals emerged")
T10.males.females.plot

#Sex ratio over time per treatment
Cumul_dat <- as.data.frame(LMC_cumul_sums)
Cumul_dat$cumul_Tot <- Cumul_dat$cumul_males + Cumul_dat$cumul_females
Cumul_dat$SR <- Cumul_dat$cumul_males/(Cumul_dat$cumul_Tot)
Cumul_dat$Tot.Fdisp <- Cumul_dat$cumul_males + Cumul_dat$Females
Cumul_dat$SR.Fdisp <- Cumul_dat$cumul_males/(Cumul_dat$Tot.Fdisp)
Cumul_dat$Tot <- Cumul_dat$Males + Cumul_dat$Females
Cumul_dat$SR.disp <- Cumul_dat$Males/(Cumul_dat$Tot)


#this plot is the cumulative of males and females 
No.disp.SR.plot <-  ggplot(Cumul_dat, aes(x = Hours.post.parasitism, y = SR, color = Treatment)) +
  geom_point() + geom_smooth() + aes(ymin = 0, ymax = 1)  +  xlab("Time (days post parasitism)") + ylab("Sex ratio (proportion males)")
No.disp.SR.plot <- SR.plot + theme_classic()
No.disp.SR.plot + theme(axis.title=element_text(face="bold.italic",
                                      size="12", color="black"), legend.position="top") + scale_color_discrete(name = "# foundresses", limits = c("1", "5", "10"),
                                                                                                               labels= c("1", "5", "10"))
#this plot assumes females disperse and males stay (i.e. cumulative males, not cumulative females)

F.disp.SR.plot <-  ggplot(Cumul_dat, aes(x = Hours.post.parasitism, y = SR.Fdisp, color = Treatment)) +
  geom_point() + geom_smooth() + aes(ymin = 0, ymax = 1)  +  xlab("Time (days post parasitism)") + ylab("Sex ratio (proportion males)")
F.disp.SR.plot <- F.disp.SR.plot + theme_classic()
F.disp.SR.plot + theme(axis.title=element_text(face="bold.italic",
                                        size="12", color="black"), legend.position="top") + scale_color_discrete(name = "# foundresses", limits = c("1", "5", "10"),
                                                                                                               labels= c("1", "5", "10"))

#this plot assumes males and females disperse (not cumulative sums for either sex)

All.disp.SR.plot <-  ggplot(Cumul_dat, aes(x = Hours.post.parasitism, y = SR.disp, color = Treatment)) +
  geom_point() + geom_smooth() + aes(ymin = 0, ymax = 1)  +  xlab("Time (days post parasitism)") + ylab("Sex ratio (proportion males)")
All.disp.SR.plot <- All.disp.SR.plot + theme_classic()
All.disp.SR.plot + theme(axis.title=element_text(face="bold.italic",
                                               size="12", color="black"), legend.position="top") + scale_color_discrete(name = "# foundresses", limits = c("1", "5", "10"),
                                                                                                                        labels= c("1", "5", "10"))

#quick and dirty LMC test - does the cumulative sex ratio change with foundress number?

SR.cum.mod <- glmer(cbind(cumul_males, cumul_females) ~ Treatment*Hours.post.parasitism + (1|Tube_number), family = binomial(link = "logit"), data = Cumul_dat)
simulationOutput <- simulateResiduals(fittedModel = SR.cum.mod, plot = T, use.u = T)
#bit dodgy - probs ok?

summary(SR.cum.mod)
car::Anova(SR.cum.mod)

# what about non-cumulative?
SR.disp.mod <- glmer(cbind(Males, Females) ~ Treatment*Hours.post.parasitism + (1|Tube_number), family = binomial(link = "logit"), data = Cumul_dat)
simulationOutput <- simulateResiduals(fittedModel = SR.disp.mod, plot = T, use.u = T)

summary(SR.disp.mod)
car::Anova(SR.disp.mod)

# sex ratio becomes more female biased over time (cumulative and dispersal)

# if we assume females disperse:
SR.F.disp.mod <- glmer(cbind(cumul_males, Females) ~ Treatment*Hours.post.parasitism + (1|Tube_number), family = binomial(link = "logit"), data = Cumul_dat)
simulationOutput <- simulateResiduals(fittedModel = SR.F.disp.mod, plot = T, use.u = T)

summary(SR.F.disp.mod)
car::Anova(SR.F.disp.mod)

#no effect of treatment on cumulative SR, disp SR or F disp SR
#only sig effect is time

#what about total SR

LMC_Tube_sums <- 
  LMC.Lysi %>%
    group_by(Treatment, Tube_number) %>%
  summarise(
    count_males = sum(Males),
    count_females =sum(Females),
  )

SR.Tot.mod <- glm(cbind(count_males, count_females) ~ Treatment,  family = binomial(link = "logit"), data = LMC_Tube_sums)
summary(SR.Tot.mod)
car::Anova(SR.Tot.mod)


#No effect of treatment on overall sex ratio - Lysiphlebus aren't allocating sex according to LMC
#Sex ratio varies the most based on timing

#get binomial confidence intervals and means

LMC_Treat_sums <- 
  LMC.Lysi %>%
  group_by(Treatment) %>%
  summarise(
    count_males = sum(Males),
    count_females =sum(Females),
  )


LMC_Treat_sums$Total <- LMC_Treat_sums$count_males + LMC_Treat_sums$count_females
Biconfint_overall_SR <- binconf(x=LMC_Treat_sums$count_males, n=LMC_Treat_sums$Total, alpha=.05)
Bi.means <- round(Biconfint_overall_SR, digits = 3)
Bi.means <- as.data.frame(Bi.means)
Treatment <- c(1,10,5)

Bi.means$Treatment <- as.character(Treatment)

Overall.SR.plot <-  ggplot(Bi.means, aes(x = Treatment, y = PointEst)) + 
  geom_bar(position=position_dodge(), stat="identity")  + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) + 
                xlab("# foundresses") + ylab("Overall sex ratio (proportion males)")

positions <- c("1", "5", "10")
Overall.SR.plot <- Overall.SR.plot + scale_x_discrete(limits = positions) + theme_classic()
Overall.SR.plot + theme(axis.title=element_text(face="bold.italic",
                                                 size="12", color="black"), legend.position="top") 

#sort out Y tick marks                                                                                                                          
