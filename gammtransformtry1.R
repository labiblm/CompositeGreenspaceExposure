

library(rgdal)
library(sf)
library(spdep)
library(spatialreg)
library(ggplot2)


#Data file
CompositeLSOA <- read_sf(dsn = ".", layer = "CompositefulldataF")

#Neighboirhoods and weights objects
compositeLSOA_nbq <- poly2nb(CompositeLSOA, queen=TRUE) #Queen’s neighborhood
summary(compositeLSOA_nbq)

compositeLSOA_nbq_w <- nb2listw(compositeLSOA_nbq) #Queen’s neighborhood wights
summary(compositeLSOA_nbq_w)

YrLLGamma <- dgamma(CompositeLSOA$Yrpotlife, shape=.83, rate= 0.0213)



ComLSOA_LAGGAMM <- lagsarlm (YrLLGamma ~ Compositem + SHDiv + Income_Sco + Education_ + Crime_Scor + no2_mean_a, 
                             data = CompositeLSOA, compositeLSOA_nbq_w)

summary(ComLSOA_LAGGAMM)

bptest.sarlm (ComLSOA_LAGGAMM)

moran.mc(ComLSOA_LAGGAMM$residuals, compositeLSOA_nbq_w, 999)