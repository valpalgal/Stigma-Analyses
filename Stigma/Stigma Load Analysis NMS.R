setwd("~/Documents/,Research/1 Skypilot Flwr Evo/data")

library(ggplot2)
library(plyr)
library(nlme)

stigma <- read.csv("Skypilot 2017 (MASTER) Pressed Skypilots (past & pres) Tab2 9-18-18.csv", header=TRUE) # all data included (final)

cols <- c("Specimen.No.", "Year", "Era", "Scent.Morph", "Date", "Location", "Altitude..m.", "Habitat..NMS.",
          "Corolla.Flare_1", "Sepal.Length_1", "Tube.Length_1", "Lobe.Length_1", "Lobe.Width_1", "No.PVgrains_stigma", 
          "No.PVgrains_off", "No.OtherGrains", "No.OtherTypes", "CP.No.PVgrains_stigma", "CP.No.PVgrains_off",
          "CP.No.OtherGrains", "CP.No.OtherTypes", "CP.Purity", "No.TotalPVgrains", "CP.No.TotalPVgrains")
stigma.sm <- subset(stigma, select = cols) #sets to raw data columns, I did not save purity, because it is calculated incorrectly

#-----# summarize data & add columns where appropriate#-----#
ddply(stigma.sm, .(Era,Location, Altitude..m., Scent.Morph), summarise, con=mean(No.PVgrains_stigma, na.rm=T), het=mean(No.OtherGrains, na.rm=T))

# stigma.sm <- subset(stigma.sm, na.rm = TRUE) #Removes all NA flower measurements/missing flowers
# stigma.sm <- subset(stigma.sm, (stigma$No.PVgrains_stigma+stigma$No.PVgrains_off)>0 | 
#                    (stigma$CP.No.PVgrains_stigma+stigma$CP.No.PVgrains_off)>0 | stigma$No.OtherGrains>0 | 
#                    stigma$CP.No.OtherGrains>0) #removes any stigmas with no PV or other pollen, to find purity
stigma.sm$No.TotalPVgrains <- stigma.sm$No.PVgrains_stigma + stigma.sm$No.PVgrains_off
stigma.sm$CP.TotalPVgrains <- stigma.sm$CP.No.PVgrains_stigma + stigma.sm$CP.No.PVgrains_off
stigma.sm$Purity.tot <- stigma.sm$No.TotalPVgrains/(stigma.sm$No.TotalPVgrains + stigma.sm$No.OtherGrains)
stigma.sm$Purity.tot <- as.numeric(as.character(stigma.sm$Purity.tot))
stigma.sm$Purity.on <- stigma.sm$No.PVgrains_stigma/(stigma.sm$No.PVgrains_stigma + stigma.sm$No.OtherGrains)
# stigma.sm$CP.Purity <- (stigma.sm$CP.No.PVgrains_stigma.sm + stigma.sm$CP.No.PVgrains_off)/stigma.sm$CP.No.OtherGrains #calculating purity of ours and Carlys data
# stigma.sm$Purity <- as.numeric(as.character(stigma.sm$Purity))
stigma.sm$Corolla.Length <- stigma.sm$Tube.Length_1 + stigma.sm$Lobe.Length_1

#-----# subsets for analysis #-----#
# subset to complete cases
stigma.sm <- subset(stigma.sm, stigma.sm$Corolla.Flare_1 > 0 & stigma.sm$Sepal.Length_1 > 0 & stigma.sm$Tube.Length_1 > 0 & 
                      stigma.sm$Lobe.Length_1 >0 & stigma.sm$Lobe.Width_1 >0) #eliminates all flowers without complete measurements
# subset by morph
stigma.sw <- subset(stigma.sm, Scent.Morph =="Sweet") #subset to sweet morph
stigma.sk <- subset(stigma.sm, Scent.Morph =="Skunky") #subset to skunky morph

# nohab <- subset(stigma.sm, select = -c(stigma$Habitat..NMS.)) #subset to eliminate habitat
# flwr.sw <- subset(flwr, Scent.Morph %in% "Sweet") #subset for sweet flower measurements
# stigma.old <- subset(stigma.sm, Era == "past") #selects for old era
# stigma.new <- subset(stigma.sm, Era == "present") #selects for new era
# stigma.sw.old <- subset(stigma.sw, Era == "past") #selects for old sweet
# stigma.sw.new <- subset(stigma.sw, Era == "present") #selects for new sweet
# stigma.sw$No.Totalgrains <- stigma.sw$No.TotalPVgrains+stigma.sw$No.OtherGrains

#-----# Exploration of Replication #-----#
ddply(stigma.sm, .(Scent.Morph, Location, Habitat..NMS., Era), summarise, N.TL=length(Tube.Length_1[!is.na(Tube.Length_1)]), 
      N.St=length(No.PVgrains_stigma[!is.na(No.PVgrains_stigma)]))
ddply(stigma.sm, .(Location, Habitat..NMS., Scent.Morph), summarise, N=length(No.TotalPVgrains))
ddply(stigma.sm, .(Location, Habitat..NMS.), summarise, N= length(No.TotalPVgrains))
ddply(stigma.sw, .(Location, Habitat..NMS.), summarise, N=length(No.TotalPVgrains))

#-----# flower measurment combined via pca #-----# WITHOUT FLOWER NUMBER!!!!!
#SWEET
measures <- c("Corolla.Flare_1", "Corolla.Length","Lobe.Width_1")
sw.size<-as.matrix(subset(stigma.sw, select=measures))
sw.size<-sw.size[complete.cases(sw.size),]
pc<-princomp(sw.size)
print(pc)
plot(pc, type="l")
summary(pc)
loadings(pc)
pc.extracted <- as.data.frame(predict(pc, newdata=stigma.sw))
stigma.sw$PC1 <- -pc.extracted$Comp.1 # all floral measurements negatively correclated with pc1, so I reversed it to be more intuitive
stigma.sw$PC2 <- pc.extracted$Comp.2 # only tube length is negatively correlated with pc2
plot(stigma.sw$PC1, stigma.sw$PC2, ylim=c(-1,1), xlim=c(-1.5,1.5), ylab="PC2 (26.2%)", xlab="PC1 (68.0%)")


####----BOTH MORPHS----####
#Total PVgrains distribution
qqnorm(stigma.sm$No.TotalPVgrains)
qqline(stigma.sm$No.TotalPVgrains)
hist(stigma.sm$No.TotalPVgrains)
qqnorm(sqrt(stigma.sm$No.TotalPVgrains))
qqline(sqrt(stigma.sm$No.TotalPVgrains))
hist(sqrt(stigma.sm$No.TotalPVgrains))
qqnorm(log1p(stigma.sm$No.TotalPVgrains))
qqline(log1p(stigma.sm$No.TotalPVgrains))
hist(log1p(stigma.sm$No.TotalPVgrains))

fit <- lme(log1p(No.TotalPVgrains) ~ No.OtherGrains*Era*Scent.Morph, random = ~1|Habitat..NMS., 
             data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(Era,No.TotalPVgrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(breaks=c(0,500,1000,2000,4000), trans="log1p") +
  ylab("Total onspecific grains") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sm, aes(x=No.OtherGrains,y=No.TotalPVgrains, shape=Era, color=Era)) +
  geom_point(size=4) + 
  scale_y_continuous(breaks=c(0,500,1000,2000,4000), trans="log1p") +
  scale_colour_manual(values=c("black")) +
  ylab("PV grains per flower") +
  xlab("Other grains per flower") + 
  scale_colour_manual(values=c("blue", "red")) +
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(No.TotalPVgrains) ~ No.OtherTypes*Era, random = ~1|Habitat..NMS.,
           data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(No.OtherTypes, No.TotalPVgrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  facet_wrap(~Era) +
  scale_y_continuous(trans="log1p") +
  ylab("PV grains per plant") +
  xlab("No. Other Pollen Types") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sm, aes(x=No.OtherTypes,y=No.TotalPVgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(No.TotalPVgrains) ~ No.OtherGrains*No.OtherTypes, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(No.OtherGrains,No.TotalPVgrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  facet_wrap(~No.OtherTypes) +
  scale_y_continuous(trans="log") +
  ylab("PV grains per plant") +
  xlab("No. Other Types") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sm, aes(x=No.OtherGrains,y=No.TotalPVgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(No.TotalPVgrains) ~ Era*Tube.Length_1, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(Tube.Length_1,No.TotalPVgrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  facet_wrap(~Era) +
  scale_y_continuous(trans="log") +
  ylab("PV grains per plant") +
  xlab("Tube Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sm, aes(x=Tube.Length_1,y=No.TotalPVgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(No.TotalPVgrains) ~ Era*Corolla.Flare_1, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(x=Corolla.Flare_1,y=No.TotalPVgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
ggplot(stigma.sm, aes(x=Corolla.Flare_1,y=No.TotalPVgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  xlim(0.5,2)+
  scale_y_continuous(breaks=c(0,150,300,600,1200,2400), trans="log10") +
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  ylab("Total conspecific pollen") + 
  xlab("Corolla flare (cm)") +    theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(No.TotalPVgrains) ~ Location*Era, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(Location,No.TotalPVgrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  facet_wrap(~Era) +
  scale_y_continuous(trans="log") +
  ylab("PV grains per plant") +
  xlab("Location") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

#Total Purity distribution
qqnorm(stigma.sm$Purity)
qqline(stigma.sm$Purity)
hist(stigma.sm$Purity)
qqnorm(sqrt(stigma.sm$Purity))
qqline(sqrt(stigma.sm$Purity))
hist(sqrt(stigma.sm$Purity))
qqnorm(log(stigma.sm$Purity))
qqline(log(stigma.sm$Purity))
hist(log(stigma.sm$Purity))

fit <- lme(log(Purity) ~ Era*Tube.Length_1, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(Era,Purity)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="log") +
  ylab("Purity") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sm, aes(x=Tube.Length_1,y=Purity, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(Purity) ~ Era*Corolla.Flare_1, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(x=Corolla.Flare_1,y=Purity, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(Purity) ~ Scent.Morph*Era, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(Era,Purity)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="log") +
  facet_wrap(~Scent.Morph) +
  ylab("Purity") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

fit <- lme(log(Purity) ~ Location*Era, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sm, aes(Era,Purity)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="log") +
  facet_wrap(~Location) +
  ylab("Purity") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

####----SWEET MORPH----####
#-----# Total PV grains #-----#
qqnorm(stigma.sw$No.TotalPVgrains)
qqline(stigma.sw$No.TotalPVgrains)
hist(stigma.sw$No.TotalPVgrains)
qqnorm(sqrt(stigma.sw$No.TotalPVgrains))
qqline(sqrt(stigma.sw$No.TotalPVgrains))
hist(sqrt(stigma.sw$No.TotalPVgrains))
qqnorm(log1p(stigma.sw$No.TotalPVgrains))
qqline(log1p(stigma.sw$No.TotalPVgrains))
hist(log(stigma.sw$No.TotalPVgrains))

fit <- lme(sqrt(No.TotalPVgrains) ~ No.OtherGrains, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)

fit <- lme(sqrt(No.TotalPVgrains) ~ Era*Tube.Length_1, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era, No.TotalPVgrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(breaks=c(0,100,200,400,600,800),trans="sqrt") +
  ylab("Total conspecific pollen grains") +
  xlab("Era") +
  geom_hline(yintercept=0) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=Tube.Length_1,y=No.TotalPVgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(sqrt(No.TotalPVgrains) ~ Era*Corolla.Flare_1, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era, No.TotalPVgrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(breaks=c(0,100,200,400,600,800),trans="sqrt") +
  ylab("Total conspecific pollen grains") +
  xlab("Era") +
  geom_hline(yintercept=0) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=Corolla.Flare_1,y=No.TotalPVgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(sqrt(No.TotalPVgrains) ~ Era*PC1, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era, No.TotalPVgrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(breaks=c(0,100,200,400,600,800),trans="sqrt") +
  ylab("Total conspecific pollen grains") +
  xlab("Era") +
  geom_hline(yintercept=0) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=PC1,y=No.TotalPVgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

#purity
qqnorm(stigma.sw$Purity.tot)
qqline(stigma.sw$Purity.tot)
hist(stigma.sw$Purity.tot)

fit <- lme((Purity.tot) ~ Era*Tube.Length_1, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era, Purity.tot)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  ylab("Purity") +
  xlab("Era") +
  geom_hline(yintercept=0) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=Tube.Length_1,y=Purity.tot, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme((Purity.tot) ~ Era*Corolla.Flare_1, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era, Purity.tot)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  ylab("Purity") +
  xlab("Era") +
  geom_hline(yintercept=0) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=Corolla.Flare_1,y=Purity.tot, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme((Purity.tot) ~ Era*PC1, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era, Purity.tot)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  ylab("Purity") +
  xlab("Era") +
  geom_hline(yintercept=0) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=PC1,y=Purity.tot, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))



#-----# PV grains on stigma #-----#
qqnorm(stigma.sw$No.PVgrains_stigma)
qqline(stigma.sw$No.PVgrains_stigma)
hist(stigma.sw$No.PVgrains_stigma)
qqnorm(sqrt(stigma.sw$No.PVgrains_stigma))
qqline(sqrt(stigma.sw$No.PVgrains_stigma))
hist(sqrt(stigma.sw$No.PVgrains_stigma))
qqnorm(log1p(stigma.sw$No.PVgrains_stigma))
qqline(log1p(stigma.sw$No.PVgrains_stigma))
hist(log(stigma.sw$No.PVgrains_stigma))

fit <- lme(sqrt(No.PVgrains_stigma) ~ No.OtherGrains, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)

fit <- lme(sqrt(No.PVgrains_stigma) ~ Era*Tube.Length_1, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era, No.PVgrains_stigma)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(breaks=c(0,100,200,400,600,800),trans="sqrt") +
  ylab("Total conspecific pollen grains") +
  xlab("Era") +
  geom_hline(yintercept=0) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=Tube.Length_1,y=No.PVgrains_stigma, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(sqrt(No.PVgrains_stigma) ~ Era*Corolla.Flare_1, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era, No.PVgrains_stigma)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(breaks=c(0,100,200,400,600,800),trans="sqrt") +
  ylab("Total conspecific pollen grains") +
  xlab("Era") +
  geom_hline(yintercept=0) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=Corolla.Flare_1,y=No.PVgrains_stigma, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(sqrt(No.PVgrains_stigma) ~ Era*PC1, random = ~1|Location/Habitat..NMS.,
           data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era, No.PVgrains_stigma)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(breaks=c(0,100,200,400,600,800),trans="sqrt") +
  ylab("Total conspecific pollen grains") +
  xlab("Era") +
  geom_hline(yintercept=0) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=PC1,y=No.PVgrains_stigma, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))



#Total grains, sweet
qqnorm(stigma.sw$No.Totalgrains)
qqline(stigma.sw$No.Totalgrains)
hist(stigma.sw$No.Totalgrains)
qqnorm(sqrt(stigma.sw$No.Totalgrains))
qqline(sqrt(stigma.sw$No.Totalgrains))
hist(sqrt(stigma.sw$No.Totalgrains))
qqnorm(log(stigma.sw$No.Totalgrains))
qqline(log(stigma.sw$No.Totalgrains))
hist(log(stigma.sw$No.Totalgrains))

fit <- lme(log(No.Totalgrains) ~ Era*No.OtherTypes, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.OtherTypes,y=No.Totalgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(No.Totalgrains) ~ Era, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era,No.Totalgrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="sqrt") +
  ylab("log(Total grains per stigma)") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=No.OtherTypes,y=No.Totalgrains, shape=Era, color=Era)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("black")) +
  ylab("Total Pollen Per Stigma") +
  xlab("log(Other Pollen Types)") + 
  scale_colour_manual(values=c("blue", "red")) +
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))


##No other grains vs no. other types
qqnorm(stigma.sw$No.OtherGrains)
qqline(stigma.sw$No.OtherGrains)
hist(stigma.sw$No.OtherGrains)
qqnorm(sqrt(stigma.sw$No.OtherGrains))
qqline(sqrt(stigma.sw$No.OtherGrains))
hist(sqrt(stigma.sw$No.OtherGrains))
qqnorm(log(stigma.sw$No.OtherGrains))
qqline(log(stigma.sw$No.OtherGrains))
hist(log(stigma.sw$No.OtherGrains))

fit <- lme(sqrt(No.OtherGrains) ~ No.OtherTypes, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
fit <- lme(sqrt(No.OtherGrains) ~ Era, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era,No.OtherGrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="sqrt") +
  ylab("sqrt(Total other grains per stigma)") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=No.OtherTypes,y=No.Totalgrains, shape=Era, color=Era)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("black")) +
  ylab("Other Pollen Per Stigma") +
  xlab("Number Other Types") + 
  scale_colour_manual(values=c("blue", "red")) +
  geom_hline(yintercept=0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))



##Do Tuckey's post hoc
#Sweet Purity

qqnorm(stigma.sw$Purity)
qqline(stigma.sw$Purity)
hist(stigma.sw$Purity)
qqnorm(sqrt(stigma.sw$Purity))
qqline(sqrt(stigma.sw$Purity))
hist(sqrt(stigma.sw$Purity))
qqnorm(log(stigma.sw$Purity))
qqline(log(stigma.sw$Purity))
hist(log(stigma.sw$Purity))

fit <- lme(log(Purity) ~ Era*Tube.Length_1, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era,Purity)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="log") +
  ylab("Purity") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=Tube.Length_1,y=Purity, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(Purity) ~ Era*No.TotalPVgrains, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.TotalPVgrains,y=Purity, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(Purity) ~ Era*No.PVgrains_stigma, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.PVgrains_stigma,y=Purity, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(Purity) ~ Era*Corolla.Flare_1, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=Corolla.Flare_1,y=Purity, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(Purity) ~ Scent.Morph*Era, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era,Purity)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="log") +
  facet_wrap(~Scent.Morph) +
  ylab("Purity") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

fit <- lme(log(Purity) ~ Location*Era, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era,Purity)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="log") +
  facet_wrap(~Location) +
  ylab("Purity") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

#Flower size, sweet morph only
qqnorm(stigma.sw$Tube.Length_1)
qqline(stigma.sw$Tube.Length_1)
hist(stigma.sw$Tube.Length_1)
qqnorm(sqrt(stigma.sw$Tube.Length_1))
qqline(sqrt(stigma.sw$Tube.Length_1))
hist(sqrt(stigma.sw$Tube.Length_1))
qqnorm(log(stigma.sw$Tube.Length_1))
qqline(log(stigma.sw$Tube.Length_1))
hist(log(stigma.sw$Tube.Length_1))

fit <- lme(Tube.Length_1 ~ Era*No.TotalPVgrains, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era,Tube.Length_1)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="log") +
  ylab("Tube Length") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))
ggplot(stigma.sw, aes(x=Tube.Length_1,y=No.TotalPVgrains, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(Tube.Length_1 ~ Era*Purity, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=Tube.Length_1,y=Purity, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(Tube.Length_1 ~ Era*No.PVgrains_stigma, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=Tube.Length_1,y=No.PVgrains_stigma, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(Tube.Length_1 ~ Era*No.OtherTypes, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.OtherTypes,y=Tube.Length_1, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(Tube.Length_1 ~ No.OtherGrains*No.PVgrains_off*No.OtherTypes, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.PVgrains_off,y=Tube.Length_1)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(Tube.Length_1 ~ Purity*No.OtherTypes, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=Tube.Length_1,y=Purity)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
ggplot(stigma.sw, aes(x=No.OtherTypes,y=Tube.Length_1)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
p <- plot_ly(x = stigma.sw$No.OtherTypes, y = stigma.sw$Tube.Length_1, z = stigma.sw$Purity) add_surface() ##I think my version of R is too old to run this package...only partially installed?













#Sweet other pollen by era and tube length 
qqnorm(stigma.sw$No.OtherGrains)
qqline(stigma.sw$No.OtherGrains)
hist(stigma.sw$No.OtherGrains)
qqnorm(sqrt(stigma.sw$No.OtherGrains))
qqline(sqrt(stigma.sw$No.OtherGrains))
hist(sqrt(stigma.sw$No.OtherGrains))
qqnorm(log(stigma.sw$No.OtherGrains))
qqline(log(stigma.sw$No.OtherGrains))
hist(log(stigma.sw$No.OtherGrains))

fit <- lme(sqrt(No.OtherGrains) ~ Era, random = ~1|Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Era,No.OtherGrains)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="sqrt") +
  ylab("No. Other Grains") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

fit <- lme(sqrt(No.OtherGrains) ~ Tube.Length_1, random = ~1|Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.OtherGrains,y=Tube.Length_1)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

#PV on stigma
qqnorm(stigma.sm$No.PVgrains_stigma)
qqline(stigma.sm$No.PVgrains_stigma)
hist(stigma.sm$No.PVgrains_stigma)
qqnorm(sqrt(stigma.sm$No.PVgrains_stigma))
qqline(sqrt(stigma.sm$No.PVgrains_stigma))
hist(sqrt(stigma.sm$No.PVgrains_stigma))
qqnorm(log(stigma.sm$No.PVgrains_stigma))
qqline(log(stigma.sm$No.PVgrains_stigma))
hist(log(stigma.sm$No.PVgrains_stigma))

fit <- lme(sqrt(No.PVgrains_stigma) ~ Era, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)
fit <- lme(sqrt(No.PVgrains_stigma) ~ Tube.Length_1, random = ~1|Habitat..NMS., data = stigma.sm, method = "ML", na.action = na.omit)
anova.lme(fit)










#Sweet pv on vs tube length
ggplot(stigma.sw, aes(x=Purity,y=No.PVgrains_stigma, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

#Sweet PV grains on stigma vs tube length
ggplot(stigma.sw, aes(x=Tube.Length_1,y=No.PVgrains_stigma, color=Era, shape=Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))



#####ADDITIONAL ANALYSIS AFTER MEETNG W/ DR. S
###PV on stigma analysis

qqnorm(stigma.sw$No.PVgrains_off)
qqline(stigma.sw$No.PVgrains_off)
hist(stigma.sw$No.PVgrains_off)
qqnorm(sqrt(stigma.sw$No.PVgrains_off))
qqline(sqrt(stigma.sw$No.PVgrains_off))
hist(sqrt(stigma.sw$No.PVgrains_off))
qqnorm(log(stigma.sw$No.PVgrains_off))
qqline(log(stigma.sw$No.PVgrains_off))
hist(log(stigma.sw$No.PVgrains_off))

fit <- lme(sqrt(No.PVgrains_off) ~ Era*No.OtherGrains, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.OtherGrains,y=No.PVgrains_off, color=Era, shape = Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
ggplot(stigma.sw, aes(Era, No.PVgrains_off)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="sqrt") +
  ylab("sqrt(PV grains off stigma)") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

fit <- lme(sqrt(No.PVgrains_off) ~ Era*Tube.Length_1, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=Tube.Length_1,y=No.PVgrains_off, color=Era, shape = Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_sqrt()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(sqrt(No.PVgrains_off) ~ Era*Tube.Length_1*No.OtherGrains, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)

fit <- lme(sqrt(No.PVgrains_off) ~ Era*Location, random = ~1|Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Location, No.PVgrains_stigma)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="sqrt") +
  ylab("sqrt(PV grains on stigma)") +
  xlab("Location") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

##PV off stigma analysis
qqnorm(stigma.sw$No.PVgrains_stigma)
qqline(stigma.sw$No.PVgrains_stigma)
hist(stigma.sw$No.PVgrains_stigma)
qqnorm(sqrt(stigma.sw$No.PVgrains_stigma))
qqline(sqrt(stigma.sw$No.PVgrains_stigma))
hist(sqrt(stigma.sw$No.PVgrains_stigma))
qqnorm(log(stigma.sw$No.PVgrains_stigma))
qqline(log(stigma.sw$No.PVgrains_stigma))
hist(log(stigma.sw$No.PVgrains_stigma))

fit <- lme(sqrt(No.PVgrains_stigma) ~ Era*No.OtherGrains, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.OtherGrains,y=No.PVgrains_stigma, color=Era, shape = Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
ggplot(stigma.sw, aes(No.PVgrains_stigma, Era)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="sqrt") +
  ylab("PV grains per plant") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

fit <- lme(sqrt(No.PVgrains_stigma) ~ Era*Tube.Length_1, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=Tube.Length_1,y=No.PVgrains_stigma, color=Era, shape = Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(sqrt(No.PVgrains_stigma) ~ Era*Tube.Length_1*No.OtherGrains, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)

fit <- lme(sqrt(No.PVgrains_stigma) ~ Era*Location, random = ~1|Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(Location, No.PVgrains_off)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="sqrt") +
  ylab("sqrt(PV grains off stigma)") +
  xlab("Location") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

##Purity
qqnorm(stigma.sw$Purity)
qqline(stigma.sw$Purity)
hist(stigma.sw$Purity)
qqnorm(sqrt(stigma.sw$Purity))
qqline(sqrt(stigma.sw$Purity))
hist(sqrt(stigma.sw$Purity))
qqnorm(log(stigma.sw$Purity))
qqline(log(stigma.sw$Purity))
hist(log(stigma.sw$Purity))

fit <- lme(log(Purity) ~ Era*Tube.Length_1, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=Tube.Length_1,y=No.PVgrains_stigma, color=Era, shape = Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
ggplot(stigma.sw, aes(Era, Purity)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="log") +
  ylab("Purity") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

fit <- lme(log(Purity) ~ No.Totalgrains*Era, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.Totalgrains,y=Purity, color=Era, shape = Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))

fit <- lme(log(Purity) ~ Era*Tube.Length_1*No.PVgrains_stigma*No.OtherTypes, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.PVgrains_stigma,y=Purity, color=Era, shape = Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))


##Types of pollen

qqnorm(stigma.sw$No.OtherTypes)
qqline(stigma.sw$No.OtherTypes)
hist(stigma.sw$No.OtherTypes)
qqnorm(sqrt(stigma.sw$No.OtherTypes))
qqline(sqrt(stigma.sw$No.OtherTypes))
hist(sqrt(stigma.sw$No.OtherTypes))
qqnorm(log(stigma.sw$No.OtherTypes))
qqline(log(stigma.sw$No.OtherTypes))
hist(log(stigma.sw$No.OtherTypes))

fit <- lme(sqrt(No.OtherTypes) ~ Era*Tube.Length_1, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.OtherTypes,y=Tube.Length_1, color=Era, shape = Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
ggplot(stigma.sw, aes(Era, No.OtherTypes)) +
  geom_violin(fill='gray') +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)  +
  scale_y_continuous(trans="sqrt") +
  ylab("sqrt(No. Other Types)") +
  xlab("Era") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=18))

fit <- lme(sqrt(No.OtherTypes) ~ Era*Tube.Length_1*No.OtherGrains, random = ~1|Location/Habitat..NMS., data = stigma.sw, method = "ML", na.action = na.omit)
anova.lme(fit)
ggplot(stigma.sw, aes(x=No.OtherTypes,y=No.OtherGrains, color=Era, shape = Era)) +
  scale_shape_manual(values=c(16,17)) +
  geom_point(size=4) + 
  scale_y_log10()+
  scale_colour_manual(values=c("blue", "red")) +
  geom_smooth(method='lm', fill=NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=20))
