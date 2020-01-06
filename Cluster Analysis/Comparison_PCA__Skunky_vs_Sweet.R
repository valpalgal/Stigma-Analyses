setwd("~/Desktop/Skypilot 2017")
library(ggplot2)
library(plyr)
library(nlme)
library(car)
library(emmeans)
library(multcompView)
library(multcomp)
library(tidyverse)  
library(FactoMineR)

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

stigma.new.dmfa <- subset(stigma.sm, select = c("Scent.Morph", "Corolla.Flare_1", "Sepal.Length_1", 
                                                "Tube.Length_1", "Lobe.Length_1", "Lobe.Width_1"))


#-----# Dual Multiple Factor Analysis #-----#
# DMFA is part of the multitable PCA family - it analyzes several sets of observations measured on the 
# same set of variables, provides a set of common factor scores (compromise) and projects each of the original 
# data sets onto the compromise to analyze commonalities and discrepancies

DMFA(stigma.new.dmfa, num.fact = 1, scale.unit = TRUE, ncp = 5, 
     quanti.sup = NULL, quali.sup = NULL, graph = TRUE, axes=c(1,2))