#Biominerals Compost Pot Study
#Soil Chemistry statistics

#Jakub Papik; jkbppk@gmail.com , Jakub.Papik@vscht.cz  
# University of Chemistry and Technology, Prague 


library(here)
library(vegan)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

data <- read.table(file ="chemieTAILINGS_16S2.csv",sep =";", dec=".", header=T, row.names=1 )

dataph <- select(data, c(treatment, pH))
dataorganic <- select(data, c(treatment, Organic.carbon))
datainoorganic <- select(data, c(treatment, Inorganic.carbon))
datanitrogen <- select(data, c(treatment, Total.nitrogen))
dataBa <- select(data, c(treatment, Ba))
dataCu <- select(data, c(treatment, Cu))
dataCo <- select(data, c(treatment, Co))
dataMg <- select(data, c(treatment, Mg))
dataNi <- select(data, c(treatment, Ni))
dataPb <- select(data, c(treatment, Pb))
dataSb <- select(data, c(treatment, Sb))
dataCd <- select(data, c(treatment, Cd))
dataAs <- select(data, c(treatment, As))
dataMn <- select(data, c(treatment, Mn))
dataTl <- select(data, c(treatment, Tl))
dataZn <- select(data, c(treatment, Zn))

#statistical testing
m1 <- lm(datainoorganic$Inorganic.carbon~treatment, data=data)
shapiro.test(residuals(m1)) 
pairwise.wilcox.test(datainoorganic$Inorganic.carbon, data$treatment, p.adjust.method="fdr")

m1 <- lm(dataorganic$Organic.carbon~treatment, data=dataorganic)
shapiro.test(residuals(m1)) 
kruskal.test(dataorganic$Organic.carbon~dataorganic$treatment) # 
pairwise.wilcox.test(dataorganic$Organic.carbon, dataorganic$treatment, p.adjust.method="fdr")

m1 <- lm(datanitrogen$Total.nitrogen~treatment, data=datanitrogen)
shapiro.test(residuals(m1)) 
kruskal.test(datanitrogen$Total.nitrogen~datanitrogen$treatment) # 
pairwise.wilcox.test(datanitrogen$Total.nitrogen, datanitrogen$treatment, p.adjust.method="fdr", correct=FALSE)
warnings()

m1 <- lm(dataCu$Cu~treatment, data=dataCu)
shapiro.test(residuals(m1)) 
kruskal.test(dataCu$Cu~dataCu$treatment) # 
pairwise.wilcox.test(dataCu$Cu, dataBa$treatment, p.adjust.method="fdr")

m1 <- lm(dataCo$Co~treatment, data=dataCo)
shapiro.test(residuals(m1)) 
kruskal.test(dataCo$Co~dataCo$treatment) # 
pairwise.wilcox.test(dataCo$Co, dataBa$treatment, p.adjust.method="fdr")

m1 <- lm(dataMg$Mg~treatment, data=dataMg)
shapiro.test(residuals(m1)) 
kruskal.test(dataMg$Mg~dataMg$treatment) # 
pairwise.wilcox.test(dataMg$Mg, dataMg$treatment, p.adjust.method="fdr")

m1 <- lm(dataPb$Pb~treatment, data=dataPb)
shapiro.test(residuals(m1)) 
kruskal.test(dataPb$Pb~dataPb$treatment) # 
pairwise.wilcox.test(dataPb$Pb, dataPb$treatment, p.adjust.method="fdr")

m1 <- lm(dataSb$Sb~treatment, data=dataSb)
shapiro.test(residuals(m1)) 
kruskal.test(dataSb$Sb~dataSb$treatment) # 
pairwise.wilcox.test(dataSb$Sb, dataSb$treatment, p.adjust.method="fdr")

m1 <- lm(dataAs$As~treatment, data=dataAs)
shapiro.test(residuals(m1)) 
kruskal.test(dataAs$As~dataAs$treatment) # 
pairwise.wilcox.test(dataAs$As, dataAs$treatment, p.adjust.method="fdr")

m1 <- lm(dataCd$Cd~treatment, data=dataCd)
shapiro.test(residuals(m1)) 
kruskal.test(dataCd$Cd~dataCd$treatment) # 
pairwise.wilcox.test(dataCd$Cd, dataCd$treatment, p.adjust.method="fdr")

m1 <- lm(dataMn$Mn~treatment, data=dataMn)
shapiro.test(residuals(m1)) 
kruskal.test(dataMn$Mn~dataMn$treatment) # 
pairwise.wilcox.test(dataMn$Mn, dataMn$treatment, p.adjust.method="fdr")

m1 <- lm(data$Zn$Zn~treatment, data=dataZn)
shapiro.test(residuals(m1)) 
kruskal.test(dataZn$Zn~dataZn$treatment) # 
pairwise.wilcox.test(dataZn$Zn, dataZn$treatment, p.adjust.method="fdr")

m1 <- lm(dataph$pH~treatment, data=data)
shapiro.test(residuals(m1)) 
kruskal.test(dataph$pH~dataph$treatment)
pairwise.wilcox.test(dataph$pH, data$treatment, p.adjust.method="fdr")

