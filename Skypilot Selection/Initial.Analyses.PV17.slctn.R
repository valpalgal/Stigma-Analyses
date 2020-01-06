setwd("~/Desktop/Skypilot 2017")
library(ggplot2)
library(plyr)
library(nlme)
library(car)
library(emmeans)
library(multcompView)

fruit <- read.csv("Skypilot 2017 (MASTER) Tab3 08-12-19.csv", header=T) # all data included (final)
fruit$OtherSeedcount[is.na(fruit$OtherSeedcount)] <- 0
fruit$Total.SeedsPerPlant <- fruit$No.TrmtSeeds+fruit$OtherSeedcount

cols <- c("Site", "Habitat", "Treatment", "Plant.Number", "Morph", "No.Pollinated.Total", "Total.RepUnits", "Prct.Pollinated",
          "Avg.SepL",	"Avg.TubeL",	"Avg.LobeL",	"Avg.LobW",	"Avg.CorollaFlare", "No.SeedsPerFruit", "Fruit.Maturity.Rank", 
          "Total.SeedsPerFruit", "Total.SeedsPerPlant")
fruit.sm <- subset(fruit, select=cols) # remove raw data columns

fruit.sm <- subset(fruit.sm, Avg.CorollaFlare>0) # remove those without flower measurements
fruit.sm$Prct.Pollinated <- as.numeric(as.character(fruit.sm$Prct.Pollinated)) # change Prct.Pollinated to numeric (it was a factor)
fruit.sm <- subset(fruit.sm, !(Treatment=="Supplementation" & Prct.Pollinated < 0.6666)) # remove those supplement flowers with <66.66% fruits pollinated
fruit.sm <- subset(fruit.sm, Fruit.Maturity.Rank < 3) # remove fruit maturity 3 (includes all exlusion as well) 
fruit.sm$CorollaL <- fruit.sm$Avg.TubeL + fruit.sm$Avg.LobeL
max(fruit.sm$Total.SeedsPerPlant)
fruit.sm$RelativeFitness <- fruit.sm$Total.SeedsPerPlant/138

# further subsets created for the analyses below
fruit.sw <- subset(fruit.sm, Morph=="Sweet") # subst to sweet morph only
fruit.sw.np <- subset(fruit.sw, Site %in% c("Niwot Ridge", "Penn Mountain"))
fruit.sw.nat <- subset(fruit.sw, Treatment=="Natural")
fruit.sw.sup <- subset(fruit.sw, Treatment=="Supplementation")
fruit.nat <- subset(fruit.sm, Treatment=="Natural") # subset to natural treatment only
fruit.sup <- subset(fruit.sm, Treatment=="Supplementation") # subset to supplement treatment only
fruit.np <- subset(fruit.sm, Site %in% c("Niwot Ridge", "Penn Mountain"))

flwr <- subset(fruit, Avg.CorollaFlare>0) # remove those without flower measurements
flwr.nat <- subset(flwr, Treatment=="Natural") # subset to natural treatment only

## exploration of replication
ddply(fruit.sm, .(Site, Habitat, Treatment, Morph), summarise, N=length(Total.SeedsPerPlant))
ddply(fruit.sm, .(Site, Habitat, Treatment), summarise, N=length(Total.SeedsPerPlant)) # 
ddply(fruit.nat, .(Site, Habitat, Treatment, Morph), summarise, N=length(Total.SeedsPerPlant)) # skunky not found @ NR & no NAT skunky @ CP...Need to check notes, but Jake has them

ddply(flwr, .(Site, Habitat, Treatment, Morph), summarise, N=length(Avg.CorollaFlare))
ddply(flwr, .(Site, Habitat, Morph), summarise, N=length(Avg.CorollaFlare))

## Galen flower size data (code not correct yet) ####
## we have flower size and fruit counts from Candi (1985?), but we haven't used them yet. Could be interesting to explore...
## not on individual plant basis, so would have to look at the effects at the site level...
galen.flwr <- read.csv("Galen data Tab2 3-22-18.csv", header=T)
galen.flwr$Avg.CorollaFlare <- colMeans(galen.flwr[17:19,], na.rm=T) # not sure why this isn't working; they are all numeric but getting error that "'x' must be numeric"
ddply(galen.flwr, .(Site, Habitat, Morph), summarise, N=length(Total.SeedsPerPlant))


# SUPPLEMENTATION EXPERIMENT ####
## can't use cumberland pass (no supplmement seeds were successfully collected...triple checking that)
## marginal interaction effect: total seeds per plant increased with supplementation of skunky, but not sweet, flowers
  ### Sweet has this weird relationship where sup seedset goes up with flower size. 
  ### Jake swears size differences shouldn't have influenced the effectiveness of the treatment

## sweet morph only
### no effect of habitat => remove
### no treatment effect, but Candi saw this in some years - might it relate to weaking selection? 
qqnorm(fruit.sw.np$Total.SeedsPerPlant)
qqline(fruit.sw.np$Total.SeedsPerPlant)
hist(fruit.sw.np$Total.SeedsPerPlant)
qqnorm(sqrt(fruit.sw.np$Total.SeedsPerPlant))
qqline(sqrt(fruit.sw.np$Total.SeedsPerPlant))
hist(sqrt(fruit.sw.np$Total.SeedsPerPlant))

fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment, random=~1|Site, data=fruit.sw.np, method="ML", na.action=na.omit)  
anova.lme(fit)
# extract least squares mean which accounts for random effect
marginal = lsmeans(fit, ~ Treatment)
lsm = multcomp::cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
# back transform from sqrt
lsm$mean <- lsm$lsmean^2
lsm$SE.low <- (lsm$lsmean-lsm$SE)^2
lsm$SE.high <- (lsm$lsmean+lsm$SE)^2
lsm
# Create a point graph with the least squares mean (NOT BACK TRANSFORMED) and SE bars
qplot(x = Treatment, y = mean, data = lsm) +
  geom_errorbar(aes(ymin  = SE.low, ymax  = SE.high, width = 0.15)) +
  ylab("Total seeds per plant") +
  xlab("Treatment") +
  ylim(20,60) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18)) 
# plot (not controlling for random effects)
ggplot(fruit.sw.np, aes(Treatment,Total.SeedsPerPlant)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  facet_wrap(~Morph) +
  scale_y_continuous(trans="sqrt") +
  ylab("Seeds per plant") +
  xlab("Morph") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(fruit.sw.np, aes(Treatment,Total.SeedsPerPlant)) +
  geom_boxplot() + 
  scale_y_continuous(trans="sqrt") +
  ylab("Seeds per plant") +
  xlab("Morph") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

# with habitat effects
fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment, random=~1|Site/Habitat, data=fruit.sw.np, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lm(sqrt(Total.SeedsPerPlant)~Habitat, data=fruit.sw.np)
anova.lme(fit)
fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment*Habitat, random=~1|Site, data=fruit.sw.np, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lm(sqrt(Total.SeedsPerPlant)~Site, data=fruit.sw.np)
anova.lme(fit)


# both morphs
# marginal interaction effect: total seeds per plant increased with supplementation of skunky, but not sweet, flowers
## Sweet has this weird relationship where sup seedset goes up with flower size. 
## Jake swears size differences shouldn't have influenced the effectiveness of the treatment
qqnorm(fruit.np$Total.SeedsPerPlant)
qqline(fruit.np$Total.SeedsPerPlant)
hist(fruit.np$Total.SeedsPerPlant)
qqnorm(sqrt(fruit.np$Total.SeedsPerPlant))
qqline(sqrt(fruit.np$Total.SeedsPerPlant))
hist(sqrt(fruit.np$Total.SeedsPerPlant))

fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment*Morph, random=~1|Site, data=fruit.np, method="ML", na.action=na.omit)  
anova.lme(fit)
# extract least squares mean which accounts for random effect
marginal = lsmeans(fit, ~ Treatment*Morph)
lsm = cld(marginal,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")
#back transform from sqrt
lsm$mean <- lsm$lsmean^2
lsm$SE.low <- (lsm$lsmean-lsm$SE)^2
lsm$SE.high <- (lsm$lsmean+lsm$SE)^2
lsm
# Create a point graph with the least squares mean (NOT BACK TRANSFORMED) and SE bars
qplot(x = Treatment, y = mean, data = lsm) +
  geom_errorbar(aes(ymin  = SE.low, ymax  = SE.high, width = 0.15)) +
  ylab("Total seeds per plant") +
  xlab("Treatment") +
  ylim(20,80) +
  facet_wrap(~Morph) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18)) 
# plot (not controlling for random effects)
ggplot(fruit.sm, aes(Treatment,Total.SeedsPerPlant)) +
  geom_boxplot() + 
  scale_y_continuous(trans="sqrt") +
  ylab("Seeds per plant") +
  facet_wrap(~Morph) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

# with habitat effect
fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment*Habitat, random=~1|Site, data=fruit.np, method="ML", na.action=na.omit)  
anova.lme(fit)



  

#### SELECTION sweet flowers only ####
# would experience the strongest selection and sk flowers biased by habitat (almost exclusively in the krummholz)
# including Cumberland pass - replication is low but not biased

#### Seeds per plant
# without supplementation treatment (treatment not evenly administered and missing CP - more important to include multiple sites)
qqnorm(fruit.sw.nat$Total.SeedsPerPlant)
qqline(fruit.sw.nat$Total.SeedsPerPlant)
hist(fruit.sw.nat$Total.SeedsPerPlant)
qqnorm(sqrt(fruit.sw.nat$Total.SeedsPerPlant))
qqline(sqrt(fruit.sw.nat$Total.SeedsPerPlant))
hist(sqrt(fruit.sw.nat$Total.SeedsPerPlant))

# Site
anova(lm(sqrt(Total.SeedsPerPlant)~Site, data=fruit.sw.nat, na.action=na.omit))
TukeyHSD(aov(sqrt(Total.SeedsPerPlant)~Site, data=fruit.sw.nat, na.action=na.omit))
ggplot(fruit.sw.nat, aes(x=Site,y=Total.SeedsPerPlant)) +
  geom_boxplot() + 
  scale_y_sqrt()+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

# corolla flare 
fit<-lme(sqrt(Total.SeedsPerPlant)~Avg.CorollaFlare, random=~1|Site, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
ggplot(fruit.sw.nat, aes(x=Avg.CorollaFlare,y=Total.SeedsPerPlant, color=Site)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  ylab("Seeds per plant") +
  xlab("Corolla flare") +  
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20), legend.position = "top")
fit<-lm(sqrt(Total.SeedsPerPlant)~Avg.CorollaFlare*Site, data=fruit.sw.nat, na.action=na.omit)  
Anova(fit, type="III")
# tube length
fit<-lme(sqrt(Total.SeedsPerPlant)~Avg.TubeL, random=~1|Site, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
ggplot(fruit.sw.nat, aes(x=Avg.TubeL,y=Total.SeedsPerPlant, color=Site)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  ylab("Seeds per plant") +
  xlab("Tube length") +  
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20), legend.position = "top")
fit<-lm(sqrt(Total.SeedsPerPlant)~Avg.TubeL*Site, data=fruit.sw.nat, na.action=na.omit)  
Anova(fit, type="III")
# flower measurements combined via pca ####
## flower measurements combined via pca (tube length)
cor(fruit.sw.nat$Avg.CorollaFlare,fruit.sw.nat$Avg.TubeL)
cor.test(fruit.sw.nat$Avg.CorollaFlare,fruit.sw.nat$Avg.TubeL)
plot(fruit.sw.nat$Avg.CorollaFlare,fruit.sw.nat$Avg.TubeL)
measures <- c("Avg.CorollaFlare", "Avg.TubeL")
flwr.size<-as.matrix(subset(fruit.sw.nat, select=measures))
pc<-princomp(flwr.size)
pc<-prcomp(flwr.size,center=TRUE)
print(pc)
plot(pc, type="l")
summary(pc)
pc.extracted <- as.data.frame(predict(pc, newdata=fruit.sw.nat))
fruit.sw.nat$PC1 <- pc.extracted$PC1

## flower measurment combined via pca (corolla length)
fruit.sw.nat$CorollaL <- fruit.sw.nat$Avg.TubeL + fruit.sw.nat$Avg.LobeL
cor(fruit.sw.nat$Avg.CorollaFlare,fruit.sw.nat$CorollaL)
cor.test(fruit.sw.nat$Avg.CorollaFlare,fruit.sw.nat$CorollaL)
plot(fruit.sw.nat$Avg.CorollaFlare,fruit.sw.nat$CorollaL)

measures <- c("CorollaL", "Avg.CorollaFlare")
flwr.size<-as.matrix(subset(fruit.sw.nat, select=measures))
pc<-prcomp(flwr.size,center=TRUE)
print(pc)
plot(pc, type="l")
summary(pc)
pc.extracted <- as.data.frame(predict(pc, newdata=fruit.sw.nat))
fruit.sw.nat$PC1.cl <- pc.extracted$PC1

## flower measurements combined via pca (all 4 traits - same as pressed) ####
measures <- c("Avg.CorollaFlare", "Avg.TubeL", "Avg.LobeL", "Avg.LobW")
flwr.size<-as.matrix(subset(fruit.sw.nat, select=measures))
pc<-prcomp(flwr.size,center=TRUE)
print(pc)
plot(pc, type="l")
summary(pc)
pc.extracted <- as.data.frame(predict(pc, newdata=fruit.sw.nat))
#INVERTED PC1 when added to the df, because all loadings are negative. This way, larger PC denotes a larger flower. 
fruit.sw.nat$PC1.all <- -pc.extracted$PC1 

## total seedset per plant (annual fecundity)
fit<-lme(sqrt(Total.SeedsPerPlant)~PC1.all, random=~1|Site, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
ggplot(fruit.sw.nat, aes(x=PC1,y=Total.SeedsPerPlant, color=Site)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  ylab("Seeds per plant") +
  xlab("PC1") +  
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20), legend.position = "top")

fit<-lme(sqrt(Total.SeedsPerPlant)~PC1.cl, random=~1|Site, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)

## relative fecundity (standardized by max annual fecundity) ####
max(fruit.sw.nat$Total.SeedsPerPlant)
fruit.sw.nat$RelativeFitness <- fruit.sw.nat$Total.SeedsPerPlant/138

qqnorm(fruit.sw.nat$RelativeFitness)
qqline(fruit.sw.nat$RelativeFitness)
hist(fruit.sw.nat$RelativeFitness)
qqnorm(sqrt(fruit.sw.nat$RelativeFitness))
qqline(sqrt(fruit.sw.nat$RelativeFitness))
hist(sqrt(fruit.sw.nat$RelativeFitness))

fit<-lme(sqrt(RelativeFitness)~PC1, random=~1|Site, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
ggplot(fruit.sw.nat, aes(x=PC1,y=RelativeFitness)) +
  geom_point(size=2) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("black")) +
  ylab("Relative fitness") +
  xlab("PC1") +  
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit<-lme(sqrt(RelativeFitness)~PC1.cl, random=~1|Site, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
ggplot(fruit.sw.nat, aes(x=PC1.cl,y=RelativeFitness)) +
  geom_point(size=2) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("black")) +
  xlim(-1,1) +
  ylab("Relative fitness") +
  xlab("PC1") +  
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit<-lme(sqrt(RelativeFitness)~Avg.CorollaFlare, random=~1|Site, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
ggplot(fruit.sw.nat, aes(x=CorollaL,y=RelativeFitness)) +
  geom_point(size=2) + 
  scale_y_sqrt()+
  ylab("Relative fitness") +
  xlab("Corolla flare (cm)") +  
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit<-lme(sqrt(RelativeFitness)~CorollaL, random=~1|Site, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
ggplot(fruit.sw.nat, aes(x=CorollaL,y=RelativeFitness)) +
  geom_point(size=2) + 
  scale_y_sqrt()+
  ylab("Relative fitness") +
  xlab("Corolla length (cm)") +  
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit<-lme(sqrt(RelativeFitness)~Avg.TubeL, random=~1|Site, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
ggplot(fruit.sw.nat, aes(x=Avg.TubeL,y=RelativeFitness)) +
  geom_point(size=2) + 
  scale_y_sqrt()+
  ylab("Relative fitness") +
  xlab("Flower tube length (cm)") +  
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))



# with supplementation treatment (treatment was more effective for larger flowers...Jake swears flwr size didn't affect the trmt)####
## total seedset per plant
qqnorm(fruit.sw$Total.SeedsPerPlant)
qqline(fruit.sw$Total.SeedsPerPlant)
hist(fruit.sw$Total.SeedsPerPlant)
qqnorm(sqrt(fruit.sw$Total.SeedsPerPlant))
qqline(sqrt(fruit.sw$Total.SeedsPerPlant))
hist(sqrt(fruit.sw$Total.SeedsPerPlant))

fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment*Avg.CorollaFlare, random=~1|Site/Habitat, data=fruit.sw, method="ML", na.action=na.omit)  
anova.lme(fit)
summary(fit)
ggplot(fruit.sw, aes(x=Avg.CorollaFlare,y=Total.SeedsPerPlant, color=Treatment, shape=Treatment)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
fit<-lme(sqrt(Total.SeedsPerPlant)~Avg.CorollaFlare, random=~1|Site/Habitat, data=fruit.sw.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lme(sqrt(Total.SeedsPerPlant)~Avg.CorollaFlare, random=~1|Site/Habitat, data=fruit.sw.sup, method="ML", na.action=na.omit)  
anova.lme(fit)

fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment*Avg.TubeL, random=~1|Site/Habitat, data=fruit.sw, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment*CorollaL, random=~1|Site/Habitat, data=fruit.sw, method="ML", na.action=na.omit)  
anova.lme(fit)

## relative fitness ####
fit<-lme(sqrt(RelativeFitness)~Treatment*Avg.CorollaFlare, random=~1|Site, data=fruit.sw, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lme(sqrt(RelativeFitness)~Treatment*Avg.CorollaFlare, random=~1|Site, data=fruit.sw, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lme(sqrt(RelativeFitness)~Treatment*Avg.CorollaFlare, random=~1|Site, data=fruit.sw, method="ML", na.action=na.omit)  
anova.lme(fit)

## seeds per fruit ####
qqnorm(fruit.sw$No.SeedsPerFruit)
qqline(fruit.sw$No.SeedsPerFruit)
hist(fruit.sw$No.SeedsPerFruit)
qqnorm(sqrt(fruit.sw$No.SeedsPerFruit))
qqline(sqrt(fruit.sw$No.SeedsPerFruit))
hist(sqrt(fruit.sw$No.SeedsPerFruit))

fit<-lme(sqrt(No.SeedsPerFruit)~Treatment*Avg.CorollaFlare, random=~1|Site/Habitat, data=fruit.sw, method="ML", na.action=na.omit)  
anova.lme(fit)
summary(fit)

ggplot(fruit.sw, aes(x=Avg.CorollaFlare,y=No.SeedsPerFruit, color=Treatment, shape=Treatment)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))


fit<-lme(sqrt(No.SeedsPerFruit)~Treatment*Avg.TubeL, random=~1|Site/Habitat, data=fruit.sw, method="ML", na.action=na.omit)  
anova.lme(fit)
summary(fit)

fit<-lme(sqrt(No.SeedsPerFruit)~Treatment*Avg.LobW, random=~1|Site/Habitat, data=fruit.sw, method="ML")  
anova.lme(fit)
summary(fit)

fit<-lme(sqrt(No.SeedsPerFruit)~Treatment*Avg.LobeL, random=~1|Site/Habitat, data=fruit.sw, method="ML")  
anova.lme(fit)
summary(fit)




#### both morphs: natural treatment only ####
#### Seeds per plant
# without supplementation treatment
qqnorm(fruit.nat$Total.SeedsPerPlant)
qqline(fruit.nat$Total.SeedsPerPlant)
hist(fruit.nat$Total.SeedsPerPlant)
qqnorm(sqrt(fruit.nat$Total.SeedsPerPlant))
qqline(sqrt(fruit.nat$Total.SeedsPerPlant))
hist(sqrt(fruit.nat$Total.SeedsPerPlant))
#site
anova(lm(sqrt(Total.SeedsPerPlant)~Morph, data=fruit.nat, na.action=na.omit))
anova(lm(sqrt(Total.SeedsPerPlant)~Site, data=fruit.nat, na.action=na.omit))
TukeyHSD(aov(sqrt(Total.SeedsPerPlant)~Site, data=fruit.nat, na.action=na.omit))
ggplot(fruit.nat, aes(x=Site,y=No.SeedsPerFruit)) +
  geom_boxplot() + 
  scale_y_sqrt()+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
#corolla flare
fit<-lme(sqrt(Total.SeedsPerPlant)~Avg.CorollaFlare*Morph, random=~1|Site, data=fruit.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
#tube length
fit<-lme(sqrt(Total.SeedsPerPlant)~Avg.TubeL*Morph, random=~1|Site, data=fruit.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
# flower measurements combined via pca (all 4 traits - same as pressed)
measures <- c("Avg.CorollaFlare", "Avg.TubeL", "Avg.LobeL", "Avg.LobW")
flwr.size<-as.matrix(subset(fruit.nat, select=measures))
pc<-prcomp(flwr.size,center=TRUE)
print(pc)
plot(pc, type="l")
summary(pc)
pc.extracted <- as.data.frame(predict(pc, newdata=fruit.nat))
fruit.nat$PC1.all <- pc.extracted$PC1 
# annual fecundity
fit<-lme(sqrt(Total.SeedsPerPlant)~PC1.all*Morph, random=~1|Site, data=fruit.nat, method="ML", na.action=na.omit)  
anova.lme(fit)

ggplot(fruit.nat, aes(x=Avg.CorollaFlare,y=No.SeedsPerFruit, color=Morph, shape=Treatment)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

#### all morphs: site, habitat and morph are random effects...####
#### seeds per fruit ####
qqnorm(fruit.sm$No.SeedsPerFruit)
qqline(fruit.sm$No.SeedsPerFruit)
hist((fruit.sm$No.SeedsPerFruit))
hist(sqrt(fruit.sm$No.SeedsPerFruit))
qqnorm(sqrt(fruit.sm$No.SeedsPerFruit))
qqline(sqrt(fruit.sm$No.SeedsPerFruit))

fit<-lme(sqrt(No.SeedsPerFruit)~Treatment*Avg.CorollaFlare, random=~1|Site/Habitat/Morph, data=fruit.sm, method="ML", na.action=na.omit)  
anova.lme(fit)
summary(fit) # sign. dif by treatment but difficult to interpret figures, if transformed looks like both go up with corolla flare and only intercept is different...
fruit.nat <- subset(fruit.sm, Treatment=="Natural")
fruit.sup <- subset(fruit.sm, Treatment=="Supplementation")

fit<-lme(sqrt(No.SeedsPerFruit)~Avg.CorollaFlare, random=~1|Site/Habitat, data=fruit.nat, method="ML")
anova.lme(fit)
fit<-lme(sqrt(No.SeedsPerFruit)~Avg.CorollaFlare, random=~1|Site/Habitat, data=fruit.sup, method="ML")
anova.lme(fit)

fit<-lme(sqrt(No.SeedsPerFruit)~Treatment*Avg.TubeL, random=~1|Site/Habitat/Morph, data=fruit.sm, method="ML")  
anova.lme(fit)
summary(fit)

fit<-lme(sqrt(No.SeedsPerFruit)~Avg.LobeL, random=~1|Site/Habitat/Morph, data=fruit.nat, method="ML")  
anova.lme(fit)
summary(fit)


ggplot(fruit.sm, aes(x=Avg.CorollaFlare,y=No.SeedsPerFruit, color=Treatment, shape=Treatment)) +
  geom_point(size=4) + 
  scale_colour_manual(values=c("blue", "red")) +
  xlim(0.8,2.4) +
  ylim(0,15) +
  geom_smooth(method='lm', aes(linetype=Treatment), fill=NA) +
  ylab("Number of seeds per fruit") +
  xlab("Corolla Flare (cm)") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
ggplot(fruit.sm, aes(x=Avg.CorollaFlare,y=No.SeedsPerFruit, color=Treatment, shape=Treatment)) +
  geom_point(size=4) + 
  scale_y_continuous(limits=c(0,15), trans="sqrt") +
  xlim(0.8,2.4) +
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', aes(linetype=Treatment), fill=NA) +
  ylab("Number of seeds per fruit") +
  xlab("Corolla Flare (cm)") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

ggplot(fruit.sm, aes(x=Avg.TubeL,y=No.SeedsPerFruit, color=Treatment, shape=Treatment)) +
  geom_point(size=4) + 
  scale_y_continuous(limits=c(0,15), trans="sqrt") +
  scale_colour_manual(values=c("blue", "red")) +
  ylab("Number of seeds per fruit") +
  xlab("Tube Length (cm)") +  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

ggplot(fruit.sm, aes(x=Avg.LobeL,y=No.SeedsPerFruit, color=Treatment, shape=Treatment)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

#fruit.sm$fit <- predict(fit)
#ggplot(fruit.sm,aes(Avg.CorollaFlare, No.SeedsPerFruit, group=Treatment, col=Treatment, shape=Treatment )) + 
#  scale_y_sqrt()+
#  geom_smooth(method='lm', aes(y=fit, lty=Treatment), size=0.8) +
#  geom_point() + 
#  geom_hline(yintercept=0, linetype="dashed") +
#  theme_bw()
# fit

### add tube length and lobe length to get corolla length
fruit.sm$CorollaL <- fruit.sm$Avg.TubeL + fruit.sm$Avg.TubeL
fit<-lme(sqrt(No.SeedsPerFruit)~Treatment*CorollaL, random=~1|Site/Habitat/Morph, data=fruit.sm, method="ML")  
anova.lme(fit)
summary(fit)

#### seeds per plant
qqnorm(fruit.sm$Total.SeedsPerPlant)
qqline(fruit.sm$Total.SeedsPerPlant)
hist(fruit.sm$Total.SeedsPerPlant)
qqnorm(sqrt(fruit.sm$Total.SeedsPerPlant))
qqline(sqrt(fruit.sm$Total.SeedsPerPlant))
hist(sqrt(fruit.sm$Total.SeedsPerPlant))

fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment*Avg.CorollaFlare, random=~1|Site/Habitat/Morph, data=fruit.sm, method="ML", na.action=na.omit)  
anova.lme(fit)
summary(fit) 

fit<-lme(sqrt(Total.SeedsPerPlant)~Treatment*Avg.TubeL, random=~1|Site/Habitat/Morph, data=fruit.sm, method="ML", na.action=na.omit)  
anova.lme(fit)
summary(fit)

#### FLOWER SIZE by habitat & morph ####
# no difference in flower size by scent morph or habitat => sweet and skunky have the same size flowers
# flowers are cumberland pass are smaller than niwot ridge & penn mountain. 
##  This parallels differences in seed set, but when we test for a direct relationship via regression (see below), there is no difference
##  Perhaps it has more to do with a weather effect, since selection doesn't vary among sites??? [stilll need to verify this last statement]
qqnorm(flwr.nat$Avg.CorollaFlare)
qqline(flwr.nat$Avg.CorollaFlare)
hist(flwr.nat$Avg.CorollaFlare)

fit<-lm(Avg.CorollaFlare~Site, data=flwr.nat, na.action=na.omit)  
anova(fit)
TukeyHSD(aov(Avg.CorollaFlare~Site, data=flwr.nat, na.action=na.omit))
ggplot(flwr.nat, aes(x=Site,y=Avg.CorollaFlare)) +
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
fit<-lme(Avg.CorollaFlare~Site, random=~1|Morph, data=flwr.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lme(Avg.CorollaFlare~Habitat*Morph, random=~1|Site, data=flwr.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lm(Avg.CorollaFlare~Site*Morph, data=flwr.nat, na.action=na.omit)  
Anova(fit,type="III") 
fit<-lm(Avg.CorollaFlare~Site*Morph*Habitat, data=flwr.nat, na.action=na.omit)  
anova(fit) # unbalanced, so cannot run type 3. 

qqnorm(flwr.nat$Avg.TubeL)
qqline(flwr.nat$Avg.TubeL)
hist(flwr.nat$Avg.TubeL)
fit<-lm(Avg.TubeL~Site, data=flwr.nat, na.action=na.omit)  
anova(fit)
fit<-lme(Avg.TubeL~Habitat*Morph, random=~1|Site, data=flwr.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lm(Avg.TubeL~Site*Morph, data=flwr.nat, na.action=na.omit)  
Anova(fit,type="III") 
fit<-lm(Avg.TubeL~Site*Morph*Habitat, data=flwr.nat, na.action=na.omit)  
anova(fit)# unbalanced, so cannot run type 3. 

measures <- c("Avg.CorollaFlare", "Avg.TubeL", "Avg.LobeL", "Avg.LobW")
flwr.size<-as.matrix(subset(flwr.nat, select=measures))
pc<-prcomp(flwr.size,center=TRUE)
print(pc)
plot(pc, type="l")
summary(pc)
pc.extracted <- as.data.frame(predict(pc, newdata=flwr.nat))
flwr.nat$PC1.all <- pc.extracted$PC1 
qqnorm(flwr.nat$PC1.all)
qqline(flwr.nat$PC1.all)
fit<-lm(PC1.all~Site, data=flwr.nat, na.action=na.omit)  
anova(fit)
ggplot(flwr.nat, aes(x=Site,y=PC1.all)) +
  geom_boxplot() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
fit<-lme(PC1.all~Habitat*Morph, random=~1|Site, data=flwr.nat, method="ML", na.action=na.omit)  
anova.lme(fit)
fit<-lm(PC1.all~Site*Morph, data=flwr.nat, na.action=na.omit)  
Anova(fit,type="III") 
