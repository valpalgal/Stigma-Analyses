setwd("~/Desktop/Skypilot 2017")
library(ggplot2)
library(plyr)
library(nlme)
library(car)
library(emmeans)
library(multcompView)
library(multcomp)
library(tidyverse)  
library(cluster)    
library(factoextra) 

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



#-----# cluster analysis using k-clustering#-----#
#In k-clustering, we choose the number of clusters

stigma.new.flwrSizeOnly <- subset(stigma.new, select = c("Corolla.Flare_1", "Sepal.Length_1", "Tube.Length_1", "Lobe.Length_1", "Lobe.Width_1"))

str(stigma.new.flwrSizeOnly)

stigma.new.flwrSizeScaled <- scale(stigma.new.flwrSizeOnly)

distance <- get_dist(stigma.new.flwrSizeScaled)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

k2 <- kmeans(stigma.new.flwrSizeScaled, centers = 2, nstart = 25)
str(k2)
k2
fviz_cluster(k2, data = stigma.new.flwrSizeScaled)

k2 <- kmeans(stigma.new.flwrSizeScaled, centers = 4, nstart = 25)
fviz_cluster(k2, data = stigma.new.flwrSizeScaled)

wss <- function(k) {
  kmeans(stigma.new.flwrSizeScaled, k, nstart = 10 )$tot.withinss
}
k.values <- 1:15
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

viz_nbclust(stigma.new.flwrSizeScaled, kmeans, method = "wss")


#-----# cluster analysis using hierarchical clustering #-----#
#This method allows us to form a dendrogram of clusters (we don't have to choose
#the number of clusters at the outset)

#Compare different types of hierarchical clustering methods 
#to establish which has the largest cluster structure value
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(stigma.new.flwrSizeScaled, method = x)$ac
}

map_dbl(m, ac)

hc3 <- agnes(stigma.new.flwrSizeScaled, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 


# Cut tree into 2 groups
sub_grp <- cutree(hc3, k = 2)

# Number of members in each cluster
table(sub_grp)


fviz_cluster(list(data = stigma.new.flwrSizeScaled, cluster = sub_grp))




