setwd("~/Desktop/Skypilot 2017")
library(ggplot2)
library(plyr)
library(nlme)
library(car)
library(emmeans)
library(multcompView)
library(multcomp)

stigma <- read.csv("Skypilot 2017 (MASTER) Tab2 08-12-19__Purity added.csv")
cols <- c("Specimen.No.", "Year", "Era", "Species", "Scent.Morph", "Date", "Location", "Altitude..m.",
          "Habitat..NMS.", "Corolla.Flare_1", "Sepal.Length_1", "Tube.Length_1", "Lobe.Length_1", "Lobe.Width_1", 
          "No.PVgrains_stigma", "No.PVgrains_off", "No.OtherGrains", "No.OtherTypes")
stigma.sm <- subset(stigma, select = cols)

#-----# summarize data & add columns where appropriate#-----#
ddply(stigma.sm, .(Era,Location, Altitude..m., Scent.Morph), summarise, con=mean(No.PVgrains_stigma, na.rm=T), het=mean(No.OtherGrains, na.rm=T))

stigma.sm$No.TotalPVgrains <- stigma.sm$No.PVgrains_stigma + stigma.sm$No.PVgrains_off
stigma.sm$Purity.tot <- stigma.sm$No.TotalPVgrains/(stigma.sm$No.TotalPVgrains + stigma.sm$No.OtherGrains)
stigma.sm$Purity.tot <- as.numeric(as.character(stigma.sm$Purity.tot))
stigma.sm$Purity.on <- stigma.sm$No.PVgrains_stigma/(stigma.sm$No.PVgrains_stigma + stigma.sm$No.OtherGrains)
stigma.sm$Corolla.Length <- stigma.sm$Tube.Length_1 + stigma.sm$Lobe.Length_1

#-----# subsets for analysis #-----#
# subset to complete cases
stigma.sm <- subset(stigma.sm, stigma.sm$Corolla.Flare_1 > 0 & stigma.sm$Sepal.Length_1 > 0 & stigma.sm$Tube.Length_1 > 0 & 
                      stigma.sm$Lobe.Length_1 >0 & stigma.sm$Lobe.Width_1 >0) #eliminates all flowers without complete measurements
# subset by morph
stigma.sw <- subset(stigma.sm, Scent.Morph =="Sweet") #subset to sweet morph
stigma.sk <- subset(stigma.sm, Scent.Morph =="Skunky") #subset to skunky morph
nohab <- subset(stigma.sm, select = -c(stigma$Habitat..NMS.)) #subset to eliminate habitat
stigma.new <- subset(stigma.sm, Era == "present") #selects for new era
stigma.sw.new <- subset(stigma.sw, Era == "present") #selects for new sweet
stigma.sw$No.Totalgrains <- stigma.sw$No.TotalPVgrains+stigma.sw$No.OtherGrains
stigma.sk.new <- subset(stigma.sk, Era == "present") #selects for new skunky
stigma.sk$No.Totalgrains <- stigma.sk$No.TotalPVgrains+stigma.sk$No.OtherGrains

#-----# Exploration of Replication #-----#
ddply(stigma.sm, .(Scent.Morph, Location, Habitat..NMS., Era), summarise, N.TL=length(Tube.Length_1[!is.na(Tube.Length_1)]), 
      N.St=length(No.PVgrains_stigma[!is.na(No.PVgrains_stigma)]))
ddply(stigma.sm, .(Location, Habitat..NMS., Scent.Morph), summarise, N=length(No.TotalPVgrains))
ddply(stigma.sm, .(Location, Habitat..NMS.), summarise, N= length(No.TotalPVgrains))
ddply(stigma.sw, .(Location, Habitat..NMS.), summarise, N=length(No.TotalPVgrains))
ddply(stigma.sk, .(Location, Habitat..NMS.), summarise, N=length(No.TotalPVgrains))

#-----# flower measurment combined via pca #-----# WITHOUT FLOWER NUMBER!!!!!
#SWEET
measures <- c("Corolla.Flare_1", "Corolla.Length","Lobe.Width_1")
sw.size<-as.matrix(subset(stigma.sw.new, select=measures))
sw.size<-sw.size[complete.cases(sw.size),]
pc<-princomp(sw.size)
print(pc)
plot(pc, type="l")
summary(pc)
loadings(pc)  
pc.extracted <- as.data.frame(predict(pc, newdata=sw.size))
#Val: looking at the loadings() function, I think that all floral measurements are positively correlated with pc1 when we only look at 2017 data
stigma.sw.new$PC1 <- pc.extracted$Comp.1 # not negatively correlated with pc1, so I did not reverse it 
stigma.sw.new$PC2 <- pc.extracted$Comp.2 # only tube length is negatively correlated with pc2
plot(stigma.sw.new$PC1, stigma.sw.new$PC2, ylim=c(-1,1), xlim=c(-1.5,1.5), ylab="PC2 (18.4%)", xlab="PC1 (77.1%)")

#SKUNKY
measures <- c("Corolla.Flare_1", "Corolla.Length","Lobe.Width_1")
sk.size<-as.matrix(subset(stigma.sk.new, select=measures))
sk.size<-sk.size[complete.cases(sk.size),]
pc<-princomp(sk.size)
print(pc)
plot(pc, type="l")
summary(pc)
loadings(pc) 
pc.extracted <- as.data.frame(predict(pc, newdata=sk.size))
stigma.sk.new$PC1 <- -pc.extracted$Comp.1 # all floral measurements negatively correclated with pc1, so I reversed it to be more intuitive
stigma.sk.new$PC2 <- pc.extracted$Comp.2
plot(stigma.sk.new$PC1, stigma.sk.new$PC2, ylim=c(-1,1), xlim=c(-1.5,1.5), ylab="PC2 (24.4%)", xlab="PC1 (71.4%)")

####----BOTH MORPHS 2017----####
#Total PVgrains distribution
qqnorm(stigma.new$No.TotalPVgrains)
qqline(stigma.new$No.TotalPVgrains)
hist(stigma.new$No.TotalPVgrains)
qqnorm(sqrt(stigma.new$No.TotalPVgrains))
qqline(sqrt(stigma.new$No.TotalPVgrains))
hist(sqrt(stigma.new$No.TotalPVgrains))
qqnorm(log1p(stigma.new$No.TotalPVgrains))
qqline(log1p(stigma.new$No.TotalPVgrains))
hist(log1p(stigma.new$No.TotalPVgrains))

#Effects of corolla flare*morph on Total P.v. pollen
fit <- lme(sqrt(No.TotalPVgrains) ~ Corolla.Flare_1*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of tube length*morph on Total P.v. pollen
fit <- lme(sqrt(No.TotalPVgrains) ~ Tube.Length_1*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of Corolla length*morph on Total P.v. pollen
fit <- lme(sqrt(No.TotalPVgrains) ~ Corolla.Length*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Other Grains distribution
qqnorm(stigma.new$No.OtherGrains)
qqline(stigma.new$No.OtherGrains)
hist(stigma.new$No.OtherGrains)
qqnorm(sqrt(stigma.new$No.OtherGrains))
qqline(sqrt(stigma.new$No.OtherGrains))
hist(sqrt(stigma.new$No.OtherGrains))
qqnorm(log1p(stigma.new$No.OtherGrains))
qqline(log1p(stigma.new$No.OtherGrains))
hist(log1p(stigma.new$No.OtherGrains))


#Effects of Corolla Flare * Morph on Other Pollen
fit <- lme(log1p(No.OtherGrains) ~ Corolla.Flare_1*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of tube length*morph on Other pollen
fit <- lme(log1p(No.OtherGrains) ~ Tube.Length_1*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of Corolla length*morph on Other pollen
fit <- lme(sqrt(No.OtherGrains) ~ Corolla.Length*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Val: skunky morph seems to pick up a lot more foreign pollen
#Getting message: "'stat_bindot()' using 'bins=30'.  Pick better value with 'binwidth'" 
#what does this mean?
ggplot(stigma.new, aes(Scent.Morph, No.OtherGrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(breaks=c(0,500,1000,2000,4000), trans="log1p") +
  ylab("No. Other Grains") +
  xlab("Morph") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))


#Total Purity Distribution
#Even with sqrt or log1p, not getting distributions that approx. normal here
x = stigma.new$Purity.tot
qqnorm(x)
qqline(x)
hist(x)
qqnorm(sqrt(x))
qqline(sqrt(x))
hist(sqrt(x))
qqnorm(log1p(x))
qqline(log1p(x))
hist(log1p(x))


#Effects of Corolla Flare*morph on Total Purity
fit <- lme(log1p(Purity.tot) ~ Corolla.Flare_1*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of tube length*morph on Total Purity
fit <- lme(log1p(Purity.tot) ~ Tube.Length_1*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of Corolla length*morph on Total Purity
fit <- lme(sqrt(Purity.tot) ~ Corolla.Length*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Val: Skunky seems to have lower total pollen purity
ggplot(stigma.new, aes(Scent.Morph, Purity.tot)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(breaks=c(0,500,1000,2000,4000), trans="log1p") +
  ylab("Total Pollen Purity") +
  xlab("Morph") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

#Distribution of No. Other Pollen Types
x = stigma.new$No.OtherTypes
qqnorm(x)
qqline(x)
hist(x)
qqnorm(sqrt(x))
qqline(sqrt(x))
hist(sqrt(x))
qqnorm(log1p(x))
qqline(log1p(x))
hist(log1p(x))

#Effects of Corolla Flare * No. Other Pollen Types
fit <- lme(log1p(No.OtherTypes) ~ Corolla.Flare_1*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of tube length*morph on Other pollen
fit <- lme(log1p(No.OtherTypes) ~ Tube.Length_1*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of Corolla length*morph on No. Other Pollen Types
fit <- lme(sqrt(No.OtherTypes) ~ Corolla.Length*Scent.Morph, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)



#---# Sweet Only #---#
#Total PVgrains distribution
qqnorm(stigma.sw.new$No.TotalPVgrains)
qqline(stigma.sw.new$No.TotalPVgrains)
hist(stigma.sw.new$No.TotalPVgrains)
qqnorm(sqrt(stigma.sw.new$No.TotalPVgrains))
qqline(sqrt(stigma.sw.new$No.TotalPVgrains))
hist(sqrt(stigma.sw.new$No.TotalPVgrains))
qqnorm(log1p(stigma.sw.new$No.TotalPVgrains))
qqline(log1p(stigma.sw.new$No.TotalPVgrains))
hist(log1p(stigma.sw.new$No.TotalPVgrains))

#Effects of PC1 on No. Total P.v. Pollen

fit <- lme(sqrt(No.TotalPVgrains) ~ PC1, random = ~1|Habitat..NMS., 
           data = stigma.sw.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Other Grains distribution
qqnorm(stigma.sw.new$No.OtherGrains)
qqline(stigma.sw.new$No.OtherGrains)
hist(stigma.sw.new$No.OtherGrains)
qqnorm(sqrt(stigma.sw.new$No.OtherGrains))
qqline(sqrt(stigma.sw.new$No.OtherGrains))
hist(sqrt(stigma.sw.new$No.OtherGrains))
qqnorm(log1p(stigma.sw.new$No.OtherGrains))
qqline(log1p(stigma.sw.new$No.OtherGrains))
hist(log1p(stigma.sw.new$No.OtherGrains))

#Effects of Corolla Flare on Other Pollen
fit <- lme(log1p(No.OtherGrains) ~ Corolla.Flare_1, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of tube length on Other pollen
fit <- lme(log1p(No.OtherGrains) ~ Tube.Length_1, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)

#Effects of Corolla length on Other pollen
fit <- lme(sqrt(No.OtherGrains) ~ Corolla.Length, random = ~1|Habitat..NMS., 
           data = stigma.new, method = "ML", na.action = na.omit)
anova.lme(fit)




