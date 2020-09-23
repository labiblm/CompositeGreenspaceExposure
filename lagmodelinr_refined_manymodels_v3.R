
#Modelling the relationship between exposure and health, inequality 
library(car)
library(MASS)
library(sf)
library(spdep)
library(spatialreg)
library(ggplot2)
library(sp)
library(RColorBrewer)
library(tmap)
library(tidyverse)
library(broom)
library(DescTools)
library(mediation)
library(mvtnorm)
library(sandwich)

library(jtools)

set.seed(123)



#Data file
CompositeLSOA <- read_sf(dsn = ".", layer = "CompositefulldataF")

plot(CompositeLSOA)

YPLL <- CompositeLSOA$Yrpotlife
CompositeLSOA$YPLL <- YPLL
summary (CompositeLSOA$Yrpotlife)

tm_shape(CompositeLSOA) + 
  tm_fill("YPLL", title = "YPLL per 1000 people", palette = "-RdBu") +
  tm_borders(alpha = 0.1) +
  tm_layout(main.title = "YPLL Distribution", main.title.size = 0.9 ,
            legend.position = c("left", "bottom"), legend.title.size = 0.8) + 
  tm_scale_bar(breaks= 10, text.size = 0.7,position=c("right", "bottom")) +
  tm_compass(type="4star", position=c("left", "top"), show.labels = 1)


#dependent names
#Avilability: AvailComN
#Accessibility: AccessS
#Visibility: GVILSOA
# Composite: Compositem

#Comilldis
#Yrpotlife



#Neighboirhoods and weights objects
compositeLSOA_nbq <- poly2nb(CompositeLSOA, queen=TRUE) #Queen’s neighborhood
summary(compositeLSOA_nbq)

compositeLSOA_nbq_w <- nb2listw(compositeLSOA_nbq, style="W") #Queen’s neighborhood wights
summary(compositeLSOA_nbq_w)


########################################################################################
#---------------------Non-spatial models for all LSOA----------------------------------#

#---------Dependent Yrpotlife------------#
#Model for green availability

modelLSOA_NDVI_yUnadj <- glm(Yrpotlife~ NDVIm, 
                             family = Gamma(link = "identity"), data= CompositeLSOA)

summary(modelLSOA_NDVI_yUnadj, Nagelkerke=T)


modelLSOA_NDVI_y <- glm(Yrpotlife~ NDVIm + Income_Sco + Crime_Scor +SHDiv, 
                        family = Gamma (link = "identity"), data= CompositeLSOA)

summary(modelLSOA_NDVI_y, Nagelkerke=T)

moran.mc(modelLSOA_NDVI_y$residuals, compositeLSOA_nbq_w, 999)



modelLSOA_GAVI_yUnadj <- glm(Yrpotlife~ AvailComN, 
                        family = Gamma(link = "identity"), data= CompositeLSOA)

summary(modelLSOA_GAVI_yUnadj, Nagelkerke=T)

PseudoR2(modelLSOA_GAVI_yUnadj, c("McFadden", "Nagel"))


modelLSOA_GAVI_y <- glm(Yrpotlife~ AvailComN + Income_Sco + Crime_Scor +SHDiv, 
                      family = Gamma (link = "identity"), data= CompositeLSOA)

summary(modelLSOA_GAVI_y)

summ(modelLSOA_GAVI_y)

PseudoR2(modelLSOA_GAVI_y, c("McFadden", "Nagel"))

moran.mc(modelLSOA_GAVI_y$residuals, compositeLSOA_nbq_w, 999)

LMtestGAVI <- lm.LMtests (modelLSOA_GAVI_y, compositeLSOA_nbq_w, test=c("LMerr", "LMlag","SARMA"))

summary(LMtestGAVI)



#Model for green accessibility
modelLSOA_GAC_y <- glm(Yrpotlife~ AccessS + Income_Sco + Crime_Scor + SHDiv, 
                     family = Gamma(link = "identity"), data= CompositeLSOA)

summary(modelLSOA_GAC_y)

PseudoR2(modelLSOA_GAC_y, c("McFadden", "Nagel"))

moran.mc(modelLSOA_GAC_y$residuals, compositeLSOA_nbq_w, 999)

LMtestGAC <- lm.LMtests (modelLSOA_GAC_y, compositeLSOA_nbq_w, test=c("LMerr", "LMlag","SARMA"))

summary(LMtestGAC)


#Model for green visibility
modelLSOA_VGVI_y <- glm(Yrpotlife~  Income_Sco + Crime_Scor + SHDiv + GVILSOA, 
                      family = Gamma(link = "identity"), data= CompositeLSOA)

summary(modelLSOA_VGVI_y)

PseudoR2(modelLSOA_VGVI_y, c("McFadden", "Nagel"))

moran.mc(modelLSOA_VGVI_y$residuals, compositeLSOA_nbq_w, 999)

LMtestVGVI <- lm.LMtests (modelLSOA_VGVI_y, compositeLSOA_nbq_w, test=c("LMerr", "LMlag","SARMA"))

summary(LMtestVGVI)


#Model for green composite greenness index
#+ pm10_mean_

modelLSOA_CEGI_y <- glm(Yrpotlife~ Income_Sco + Crime_Scor + SHDiv+ Compositem, 
                      family = Gamma(link = "identity"), data= CompositeLSOA)

summary(modelLSOA_CEGI_y)

PseudoR2(modelLSOA_CEGI_y, c("McFadden", "Nagel"))


LMtestCGEI <- lm.LMtests (modelLSOA_CEGI_y, compositeLSOA_nbq_w, test=c("LMerr", "LMlag","SARMA"))

summary(LMtestCGEI)



plot(modelLSOA_CEGI_y$residuals)

vif(modelLSOA_CEGI_y)

summ(modelLSOA_CEGI_y)

moran.mc(modelLSOA_CEGI_y$residuals, compositeLSOA_nbq_w, 999)


CompositeLSOA$resmodelLSOA_CEGI_y <- residuals(modelLSOA_CEGI_y)


tm_shape(CompositeLSOA) + 
  tm_fill("resmodelLSOA_CEGI_y", title = "GLM Residuals", palette = "-RdBu") +
  tm_borders(alpha = 0.1) +
  tm_layout(main.title = "GLM model Residuals", main.title.size = 0.9 ,
            legend.position = c("left", "bottom"), legend.title.size = 0.8)



######################################################################################
#--------------------------SPATIAL MODELS For all LSOA-------------------------------#
#INDEPENDENT
# Avilability: AvailComN
# Accessibility: AccessS
# Visibility: GVILSOA
# Composite: Compositem

# Possible Confounding adjustment
# Income Deprivation: Income_Sco
# Air pollution: pm10_mean_
# Eduction Deprivation: Education_ (highly multi-colinar with Income)
# Crime Score: Crime_Scor (Not high VIF, but has effect of model residual)
# Barriers to Housing: Barriers_H (Somewhat multicolinar)
# Diversity: SHDiv
#CompositeLSOA$Crime_Scor
CompositeLSOA$Comilldis

#---------Dependent Yrpotlife------------#


#Model for NDVI
Spatial_LAG_NDVI_yUnadj <- sacsarlm (Yrpotlife ~ NDVIm, 
                                     data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_NDVI_yUnadj, Nagelkerke=T)

moran.mc(Spatial_LAG_NDVI_yUnadj$residuals, compositeLSOA_nbq_w, 999)

#adjusted
Spatial_LAG_NDVI_y <- sacsarlm (Yrpotlife ~ NDVIm + Income_Sco + Crime_Scor + SHDiv, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_NDVI_y, Nagelkerke=T)

moran.mc(Spatial_LAG_NDVI_y$residuals, compositeLSOA_nbq_w, 999)


#Model for green availability
Spatial_LAG_GAVI_yUnadj <- sacsarlm (Yrpotlife ~ AvailComN, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_GAVI_yUnadj, Nagelkerke=T)

moran.mc(Spatial_LAG_GAVI_yUnadj$residuals, compositeLSOA_nbq_w, 999)


#adjusted
Spatial_LAG_GAVI_y <- sacsarlm (Yrpotlife ~ Income_Sco + Crime_Scor + SHDiv + AvailComN, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_GAVI_y, Nagelkerke=T)

moran.mc(Spatial_LAG_GAVI_y$residuals, compositeLSOA_nbq_w, 999)



#Model for green accessibility

Spatial_LAG_GAC_yUnadj  <- sacsarlm (Yrpotlife ~ AccessS, 
                               data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_GAC_yUnadj , Nagelkerke=T)

moran.mc(Spatial_LAG_GAC_yUnadj$residuals, compositeLSOA_nbq_w, 999)



Spatial_LAG_GAC_y <- sacsarlm (Yrpotlife ~ Income_Sco + Crime_Scor + SHDiv + AccessS, 
                               data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_GAC_y, Nagelkerke=T)

moran.mc(Spatial_LAG_GAC_y$residuals, compositeLSOA_nbq_w, 999)



#Model for green visibility

Spatial_LAG_VGVI_yUnadj <- sacsarlm (Yrpotlife ~ GVILSOA, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_VGVI_yUnadj, Nagelkerke=T)

moran.mc(Spatial_LAG_VGVI_yUnadj$residuals, compositeLSOA_nbq_w, 999)



Spatial_LAG_VGVI_y <- sacsarlm (Yrpotlife ~ Income_Sco + Crime_Scor + SHDiv + GVILSOA, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_VGVI_y, Nagelkerke=T)

moran.mc(Spatial_LAG_VGVI_y$residuals, compositeLSOA_nbq_w, 999)



#Model for green composite greenness index
#Specification and Estimation of Spatial Autoregressive Models/ Spatial autoregressive combined model  (both side take account of lag)
Spatial_LAG_CEGI_yUnadj <- sacsarlm (Yrpotlife ~ Compositem, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_CEGI_yUnadj, Nagelkerke=T)

moran.mc(Spatial_LAG_CEGI_yUnadj$residuals, compositeLSOA_nbq_w, 999)




Spatial_LAG_CEGI_y <- sacsarlm (Yrpotlife ~ Income_Sco + Crime_Scor+ SHDiv + Compositem, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_CEGI_y, Nagelkerke=T)

moran.mc(Spatial_LAG_CEGI_y$residuals, compositeLSOA_nbq_w, 999)



#errorsarlm#################################################################


Spatial_LAG_CEGI_ySE <- errorsarlm (Yrpotlife ~ Income_Sco + Crime_Scor+ SHDiv + AvailComN, 
                                data = CompositeLSOA, compositeLSOA_nbq_w)
summary(Spatial_LAG_CEGI_ySE, Nagelkerke=T)

moran.mc(Spatial_LAG_CEGI_ySE$residuals, compositeLSOA_nbq_w, 999)





#test air poll####################################################################

Spatial_LAG_CEGI_yIn <- sacsarlm (Yrpotlife ~ Income_Sco + Crime_Scor+ SHDiv + pm10_mean_ + Compositem, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_CEGI_yIn, Nagelkerke=T)

moran.mc(Spatial_LAG_CEGI_yIn$residuals, compositeLSOA_nbq_w, 999)



Spatial_LAG_GAVI_yPM <- sacsarlm (Yrpotlife ~ Income_Sco + Crime_Scor + SHDiv + pm10_mean_ + AvailComN, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_GAVI_yPM, Nagelkerke=T)

moran.mc(Spatial_LAG_GAVI_yPM$residuals, compositeLSOA_nbq_w, 999)



Spatial_LAG_GAC_yPM <- sacsarlm (Yrpotlife ~ Income_Sco + Crime_Scor + SHDiv + pm10_mean_ + AccessS, 
                                  data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_GAC_yPM, Nagelkerke=T)

moran.mc(Spatial_LAG_GAC_yPM$residuals, compositeLSOA_nbq_w, 999)



Spatial_LAG_VGVI_yPM <- sacsarlm (Yrpotlife ~ Income_Sco + Crime_Scor + SHDiv + pm10_mean_ + GVILSOA, 
                                 data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_VGVI_yPM, Nagelkerke=T)

moran.mc(Spatial_LAG_VGVI_yPM$residuals, compositeLSOA_nbq_w, 999)





W <- as(compositeLSOA_nbq_w, "CsparseMatrix")
trMC <- trW(W, type="MC")
im<-impacts(Spatial_LAG_CEGI_y, tr=trMC, R= 999)
summary(im,  zstats=T)

Spatial_LAG_CEGI_y$residuals <- residuals(Spatial_LAG_CEGI_y)

plot(Spatial_LAG_CEGI_y$residuals)

bptest.sarlm (Spatial_LAG_CEGI_y)

set.seed(123)
moran.mc(Spatial_LAG_CEGI_y$residuals, compositeLSOA_nbq_w, 999)


CompositeLSOA$Spatial_LAG_CEGI__y_cres <- residuals(Spatial_LAG_CEGI_y)

tm_shape(CompositeLSOA) + 
  tm_fill("Spatial_LAG_CEGI__y_cres", title = "Residuals", palette = "-RdBu") +
  tm_borders(alpha = 0.1) +
  tm_layout(main.title = "SAC lag model Residuals", main.title.size = 0.9,
            legend.position = c("left", "bottom"), legend.title.size = 0.8)




#######################Air Pollution Sensitivity###########################

pm10_mean_








#######################testmediation################################################


greenPM10 <- sacsarlm (pm10_mean_ ~ Compositem + Income_Sco + Crime_Scor + SHDiv, 
                       data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")

summary(greenPM10, Nagelkerke=T)


Spatial_LAG_CEGI_yGLM <- sacsarlm (YPLL ~ Compositem + Income_Sco + Crime_Scor + SHDiv + pm10_mean_, 
                                data = CompositeLSOA, compositeLSOA_nbq_w, type = "sac")
summary(Spatial_LAG_CEGI_yGLM, Nagelkerke=T)


med.output1 <- mediate (greenPM10, Spatial_LAG_CEGI_yGLM, treat = "Compositem", mediator = "pm10_mean_", robustSE = TRUE, sims = 1000)

summary(med.output1)

plot (med.output1)



greenPM102 <- glm (pm10_mean_ ~ Compositem + Income_Sco + Crime_Scor + SHDiv, 
                 family = Gamma(link = "identity"), data = CompositeLSOA)

summary(greenPM102, Nagelkerke=T)

Spatial_LAG_CEGI_y2 <- glm (Yrpotlife ~ Compositem + Income_Sco + Crime_Scor+ SHDiv + pm10_mean_, 
                            family = Gamma(link = "identity"), data = CompositeLSOA)
summary(Spatial_LAG_CEGI_y2, Nagelkerke=T)

med.output2 <- mediate(greenPM102, Spatial_LAG_CEGI_y2, treat = "Compositem", mediator = "pm10_mean_", robustSE = TRUE, sims = 1000)

summary(med.output2)

plot(med.output2)


###################SubGroups###################
# Avilability: AvailComN
# Accessibility: AccessS
# Visibility: GVILSOA
# Composite: Compositem


#---------------IMD_1_2 Analyses--------------#

#Data

LSOA_IMD_1_2 <- read_sf(dsn = ".", layer = "Quintile_group_5")

#nh
LSOA_IMD_1_2_nbq <- poly2nb(LSOA_IMD_1_2, queen=TRUE) #Queen’s neighborhood

LSOA_IMD_1_2_nbq_w <- nb2listw(LSOA_IMD_1_2_nbq, style = "W", zero.policy = TRUE) #Queen’s neighborhood wights


#non-spatial model

model_IMD1_2 <- glm(Yrpotlife~ Compositem + Income_Sco+ Crime_Scor+ SHDiv, 
              family = Gamma(link = "identity"), data= LSOA_IMD_1_2)

summary(model_IMD1_2)

moran.mc(model_IMD1_2$residuals, LSOA_IMD_1_2_nbq_w, zero.policy = TRUE, 999)



#spatial model

LSOA_IMD_1_2Sav <- sacsarlm (Yrpotlife ~AvailComN + Income_Sco+ Crime_Scor+ SHDiv, 
                           data = LSOA_IMD_1_2, LSOA_IMD_1_2_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_1_2Sav, Nagelkerke=T)


LSOA_IMD_1_2Sac <- sacsarlm (Yrpotlife ~AccessS + Income_Sco+ Crime_Scor+ SHDiv, 
                             data = LSOA_IMD_1_2, LSOA_IMD_1_2_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_1_2Sac, Nagelkerke=T)


LSOA_IMD_1_2Sv <- sacsarlm (Yrpotlife ~GVILSOA + Income_Sco+ Crime_Scor+ SHDiv, 
                             data = LSOA_IMD_1_2, LSOA_IMD_1_2_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_1_2Sv, Nagelkerke=T)


LSOA_IMD_1_2Sua <- sacsarlm (Yrpotlife ~Compositem, 
                           data = LSOA_IMD_1_2, LSOA_IMD_1_2_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_1_2Sua, Nagelkerke=T)


LSOA_IMD_1_2S <- sacsarlm (Yrpotlife ~Compositem + Income_Sco+ Crime_Scor+ SHDiv, 
                          data = LSOA_IMD_1_2, LSOA_IMD_1_2_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_1_2S, Nagelkerke=T)


W <- as(LSOA_IMD_1_2_nbq_w, "CsparseMatrix")
trMC <- trW(W, type="MC")
IMD_1_2im <-impacts(LSOA_IMD_1_2S, tr=trMC, R= 999)
summary(IMD_1_2im,  zstats=T)

LSOA_IMD_1_2S$residuals <- residuals(LSOA_IMD_1_2S)

plot(LSOA_IMD_1_2S$residuals)

bptest.sarlm (LSOA_IMD_1_2S)


moran.mc(LSOA_IMD_1_2S$residuals, LSOA_IMD_1_2_nbq_w, zero.policy = TRUE, 999)


LSOA_IMD_1_2$LSOA_IMD_1_2Sres <- residuals(LSOA_IMD_1_2S)

tm_shape(LSOA_IMD_1_2) + 
  tm_fill("LSOA_IMD_1_2Sres", title = "Residuals",  style = "fixed", breaks = spa_breaks, palette = "-RdBu") +
  tm_borders(alpha = 0.1) +
  tm_layout(main.title = "SARAR lag model Residuals", main.title.size = 0.9,
            legend.position = c("left", "bottom"), legend.title.size = 0.8)





#---------------IMD_3_4 Analyses--------------#

#Data

LSOA_IMD_3_4 <- read_sf(dsn = ".", layer = "Quintile_group_4")

plot(LSOA_IMD_3_4)

#nh
LSOA_IMD_3_4_nbq <- poly2nb(LSOA_IMD_3_4, queen=TRUE) #Queen’s neighborhood

LSOA_IMD_3_4_nbq_w <- nb2listw(LSOA_IMD_3_4_nbq, style = "W", zero.policy = TRUE) #Queen’s neighborhood wights
summary(LSOA_IMD_3_4_nbq_w)

#non-spatial model

model_IMD_3_4 <- glm(Yrpotlife~  Compositem + Income_Sco + Crime_Scor + SHDiv, 
                    family = Gamma(link = "identity"), data= LSOA_IMD_3_4)

summary(model_IMD_3_4)

plot(model_IMD_3_4)

moran.mc(model_IMD_3_4$residuals, LSOA_IMD_3_4_nbq_w, zero.policy = TRUE, 999)



#spatial model

LSOA_IMD_3_4Sav <- sacsarlm (Yrpotlife ~  AvailComN + Income_Sco+ Crime_Scor+ SHDiv, 
                           data = LSOA_IMD_3_4, LSOA_IMD_3_4_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_3_4Sav, Nagelkerke=T)


LSOA_IMD_3_4Sac <- sacsarlm (Yrpotlife ~  AccessS + Income_Sco+ Crime_Scor+ SHDiv, 
                             data = LSOA_IMD_3_4, LSOA_IMD_3_4_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_3_4Sac, Nagelkerke=T)

LSOA_IMD_3_4Sv <- sacsarlm (Yrpotlife ~  GVILSOA + Income_Sco+ Crime_Scor+ SHDiv, 
                             data = LSOA_IMD_3_4, LSOA_IMD_3_4_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_3_4Sv, Nagelkerke=T)


LSOA_IMD_3_4Sun <- sacsarlm (Yrpotlife ~  Compositem, 
                           data = LSOA_IMD_3_4, LSOA_IMD_3_4_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_3_4Sun, Nagelkerke=T)


LSOA_IMD_3_4S <- sacsarlm (Yrpotlife ~  Compositem + Income_Sco+ Crime_Scor+ SHDiv, 
                           data = LSOA_IMD_3_4, LSOA_IMD_3_4_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_3_4S, Nagelkerke=T)

W <- as(LSOA_IMD_3_4_nbq_w, "CsparseMatrix")
trMC <- trW(W, type="MC")
IMD_3_4im <-impacts(LSOA_IMD_3_4S, tr=trMC, R= 999)
summary(IMD_1_2im,  zstats=T)

LSOA_IMD_3_4S$residuals <- residuals(LSOA_IMD_3_4S)

plot(LSOA_IMD_3_4S$residuals)

bptest.sarlm (LSOA_IMD_3_4S)


moran.mc(LSOA_IMD_3_4S$residuals, LSOA_IMD_3_4_nbq_w, zero.policy = TRUE, 999)


LSOA_IMD_3_4$LSOA_IMD_3_4Sres <- residuals(LSOA_IMD_3_4S)

tm_shape(LSOA_IMD_3_4) + 
  tm_fill("LSOA_IMD_3_4Sres", title = "Residuals",  style = "fixed", breaks = spa_breaks, palette = "-RdBu") +
  tm_borders(alpha = 0.1) +
  tm_layout(main.title = "SARAR lag model Residuals", main.title.size = 0.9,
            legend.position = c("left", "bottom"), legend.title.size = 0.8)



#---------------IMD_5_6 Analyses--------------#

#Data

LSOA_IMD_5_6 <- read_sf(dsn = ".", layer = "Quintile_group_3")

plot(LSOA_IMD_5_6)

#nh
LSOA_IMD_5_6_nbq <- poly2nb(LSOA_IMD_5_6, queen=TRUE) #Queen’s neighborhood

LSOA_IMD_5_6_nbq_w <- nb2listw(LSOA_IMD_5_6_nbq, style = "W", zero.policy = TRUE) #Queen’s neighborhood wights

summary(LSOA_IMD_5_6_nbq_w)

#non-spatial model

model_IMD_5_6 <- glm(Yrpotlife~  Compositem + Income_Sco + pm10_mean_ + SHDiv, 
                     family = Gamma(link = "identity"), data= LSOA_IMD_5_6)

summary(model_IMD_5_6)

plot(model_IMD_5_6)

moran.mc(model_IMD_5_6$residuals, LSOA_IMD_5_6_nbq_w, zero.policy = TRUE, 999)



#spatial model
LSOA_IMD_5_6Sav <- sacsarlm (Yrpotlife ~ AvailComN + Income_Sco+ Crime_Scor+ SHDiv, 
                           data = LSOA_IMD_5_6, LSOA_IMD_5_6_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_5_6Sav, Nagelkerke=T)

LSOA_IMD_5_6Sac <- sacsarlm (Yrpotlife ~ AccessS + Income_Sco+ Crime_Scor+ SHDiv, 
                             data = LSOA_IMD_5_6, LSOA_IMD_5_6_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_5_6Sac, Nagelkerke=T)


LSOA_IMD_5_6Sac <- sacsarlm (Yrpotlife ~ GVILSOA + Income_Sco+ Crime_Scor+ SHDiv, 
                             data = LSOA_IMD_5_6, LSOA_IMD_5_6_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_5_6Sac, Nagelkerke=T)



LSOA_IMD_5_6Sun <- sacsarlm (Yrpotlife ~ Compositem, 
                           data = LSOA_IMD_5_6, LSOA_IMD_5_6_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_5_6Sun, Nagelkerke=T)


LSOA_IMD_5_6S <- sacsarlm (Yrpotlife ~ Compositem + Income_Sco+ Crime_Scor+ SHDiv, 
                           data = LSOA_IMD_5_6, LSOA_IMD_5_6_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_5_6S, Nagelkerke=T)

W <- as(LSOA_IMD_5_6_nbq_w, "CsparseMatrix")
trMC <- trW(W, type="MC")
IMD_5_6im <-impacts(LSOA_IMD_5_6S, tr=trMC, R= 999)
summary(IMD_5_6im,  zstats=T)

LSOA_IMD_5_6S$residuals <- residuals(LSOA_IMD_5_6S)

plot(LSOA_IMD_5_6S$residuals)

bptest.sarlm (LSOA_IMD_5_6S)


moran.mc(LSOA_IMD_5_6S$residuals, LSOA_IMD_5_6_nbq_w, zero.policy = TRUE, 999)


LSOA_IMD_5_6$LSOA_IMD_5_6Sres <- residuals(LSOA_IMD_5_6S)

tm_shape(LSOA_IMD_5_6) + 
  tm_fill("LSOA_IMD_5_6Sres", title = "Residuals",  style = "fixed", breaks = spa_breaks, palette = "-RdBu") +
  tm_borders(alpha = 0.1) +
  tm_layout(main.title = "SARAR lag model Residuals", main.title.size = 0.9,
            legend.position = c("left", "bottom"), legend.title.size = 0.8)



#---------------IMD_7_8 Analyses--------------#

#Data

LSOA_IMD_7_8 <- read_sf(dsn = ".", layer = "Quintile_group_2")

plot(LSOA_IMD_7_8)

#nh
LSOA_IMD_7_8_nbq <- poly2nb(LSOA_IMD_7_8, queen=TRUE) #Queen’s neighborhood

LSOA_IMD_7_8_nbq_w <- nb2listw(LSOA_IMD_7_8_nbq, style = "W", zero.policy = TRUE) #Queen’s neighborhood wights
summary(LSOA_IMD_7_8_nbq_w)

#non-spatial model

model_IMD_7_8 <- glm(Yrpotlife~  AccessS + Income_Sco + Crime_Scor + SHDiv, 
                     family = Gamma(link = "identity"), data= LSOA_IMD_7_8)

summary(model_IMD_7_8)

plot(model_IMD_7_8)

moran.mc(model_IMD_7_8$residuals, LSOA_IMD_7_8_nbq_w, zero.policy = TRUE, 999)



#spatial model
LSOA_IMD_7_8Sav <- sacsarlm (Yrpotlife ~  AvailComN + Income_Sco+ Crime_Scor+ SHDiv, 
                           data = LSOA_IMD_7_8, LSOA_IMD_7_8_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_7_8Sav, Nagelkerke=T)


LSOA_IMD_7_8Sac <- sacsarlm (Yrpotlife ~  AccessS + Income_Sco + Crime_Scor+ SHDiv, 
                             data = LSOA_IMD_7_8, LSOA_IMD_7_8_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_7_8Sac, Nagelkerke=T)



LSOA_IMD_7_8Sv <- sacsarlm (Yrpotlife ~  GVILSOA + Income_Sco+ Crime_Scor+ SHDiv, 
                             data = LSOA_IMD_7_8, LSOA_IMD_7_8_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_7_8Sv, Nagelkerke=T)



LSOA_IMD_7_8Sun <- sacsarlm (Yrpotlife ~  Compositem, 
                           data = LSOA_IMD_7_8, LSOA_IMD_7_8_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_7_8Sun, Nagelkerke=T)



LSOA_IMD_7_8S <- sacsarlm (Yrpotlife ~  Compositem + Income_Sco+ Crime_Scor+ SHDiv, 
                           data = LSOA_IMD_7_8, LSOA_IMD_7_8_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_7_8S, Nagelkerke=T)

W <- as(LSOA_IMD_7_8_nbq_w, "CsparseMatrix")
trMC <- trW(W, type="MC")
IMD_7_8im <-impacts(LSOA_IMD_7_8S, tr=trMC, R= 999)
summary(IMD_7_8im,  zstats=T)

LSOA_IMD_7_8S$residuals <- residuals(LSOA_IMD_7_8S)

plot(LSOA_IMD_7_8S$residuals)

bptest.sarlm (LSOA_IMD_7_8S)


moran.mc(LSOA_IMD_7_8S$residuals, LSOA_IMD_7_8_nbq_w, zero.policy = TRUE, 999)


LSOA_IMD_7_8$LSOA_IMD_7_8Sres <- residuals(LSOA_IMD_7_8S)

tm_shape(LSOA_IMD_7_8) + 
  tm_fill("LSOA_IMD_7_8Sres", title = "Residuals",  style = "fixed", breaks = spa_breaks, palette = "-RdBu") +
  tm_borders(alpha = 0.1) +
  tm_layout(main.title = "SARAR lag model Residuals", main.title.size = 0.9,
            legend.position = c("left", "bottom"), legend.title.size = 0.8)



#---------------IMD_9_10_Analyses--------------#

#Data

LSOA_IMD_9_10 <- read_sf(dsn = ".", layer = "Quintile_group_1")

plot(LSOA_IMD_9_10)

#nh
LSOA_IMD_9_10_nbq <- poly2nb(LSOA_IMD_9_10, queen=TRUE) #Queen’s neighborhood

LSOA_IMD_9_10_nbq_w <- nb2listw(LSOA_IMD_9_10_nbq, style = "W", zero.policy = TRUE) #Queen’s neighborhood wights


#non-spatial model

model_IMD_9_10 <- glm(Yrpotlife~  Compositem + Income_Sco + pm10_mean_ + SHDiv, 
                    family = Gamma(link = "identity"), data= LSOA_IMD_9_10)

summary(model_IMD_9_10)

vif(model_IMD_9_10)

moran.mc(model_IMD_9_10$residuals, LSOA_IMD_9_10_nbq_w, zero.policy = TRUE, 999)



#spatial model

LSOA_IMD_9_10Sav <- sacsarlm (Yrpotlife ~  AvailComN + Income_Sco+ Crime_Scor+ SHDiv, 
                           data = LSOA_IMD_9_10, LSOA_IMD_9_10_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_9_10Sav, Nagelkerke=T)


LSOA_IMD_9_10Sac <- sacsarlm (Yrpotlife ~  AccessS + Income_Sco+ Crime_Scor+ SHDiv, 
                            data = LSOA_IMD_9_10, LSOA_IMD_9_10_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_9_10Sac, Nagelkerke=T)


LSOA_IMD_9_10Sv <- sacsarlm (Yrpotlife ~  GVILSOA + Income_Sco+ Crime_Scor+ SHDiv, 
                              data = LSOA_IMD_9_10, LSOA_IMD_9_10_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_9_10Sv, Nagelkerke=T)



LSOA_IMD_9_10Sua <- sacsarlm (Yrpotlife ~  Compositem, 
                            data = LSOA_IMD_9_10, LSOA_IMD_9_10_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_9_10Sua, Nagelkerke=T)


LSOA_IMD_9_10S <- sacsarlm (Yrpotlife ~  Compositem + Income_Sco+ Crime_Scor+ SHDiv, 
                             data = LSOA_IMD_9_10, LSOA_IMD_9_10_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_9_10S, Nagelkerke=T)


W <- as(LSOA_IMD_9_10_nbq_w, "CsparseMatrix")
trMC <- trW(W, type="MC")
IMD_9_10im <-impacts(LSOA_IMD_9_10S, tr=trMC, R= 999)
summary(IMD_9_10im,  zstats=T)

LSOA_IMD_9_10S$residuals <- residuals(LSOA_IMD_9_10S)

plot(LSOA_IMD_9_10S$residuals)

bptest.sarlm (LSOA_IMD_9_10S)


moran.mc(LSOA_IMD_9_10S$residuals, LSOA_IMD_9_10_nbq_w, zero.policy = TRUE, 999)


LSOA_IMD_9_10$LSOA_IMD_9_10Sres <- residuals(LSOA_IMD_9_10S)

tm_shape(LSOA_IMD_9_10) + 
  tm_fill("LSOA_IMD_9_10Sres", title = "Residuals",  style = "fixed", breaks = spa_breaks, palette = "-RdBu") +
  tm_borders(alpha = 0.1) +
  tm_layout(main.title = "SARAR lag model Residuals", main.title.size = 0.9,
            legend.position = c("left", "bottom"), legend.title.size = 0.8)






#INDEPENDENT
# Avilability: AvailComN
# Accessibility: AccessS
# Visibility: GVILSOA
# Composite: Compositem








#IMD_Tercile

#MostDeprived


LSOA_IMDMostDeprived <- read_sf(dsn = ".", layer = "IMDtercile_mostDeprive")

plot(LSOA_IMDMostDeprived )

#nh
LSOA_IMDMostDeprived_nbq <- poly2nb(LSOA_IMDMostDeprived, queen=TRUE) #Queen’s neighborhood

LSOA_IMDMostDeprived_nbq_w <- nb2listw(LSOA_IMDMostDeprived_nbq, style = "W", zero.policy = TRUE) #Queen’s neighborhood wights


LSOA_IMD_IMDMostDeprivedS <- sacsarlm (Yrpotlife ~ Compositem + Income_Sco + no2_mean_a, 
                            data = LSOA_IMDMostDeprived, LSOA_IMDMostDeprived_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMD_IMDMostDeprivedS, Nagelkerke=T)




#ModDeprived


LSOA_IMDModDeprived <- read_sf(dsn = ".", layer = "IMDtercile_ModDeprive")

plot(LSOA_IMDModDeprived)

#nh
LSOA_IMDModDeprived_nbq <- poly2nb(LSOA_IMDModDeprived, queen=TRUE) #Queen’s neighborhood

LSOA_IMDModDeprived_nbq_w <- nb2listw(LSOA_IMDModDeprived_nbq, style = "W", zero.policy = TRUE) #Queen’s neighborhood wights


LSOA_IMDModDeprivedS <- sacsarlm (Yrpotlife ~ Compositem + Income_Sco + no2_mean_a, 
                                       data = LSOA_IMDModDeprived, LSOA_IMDModDeprived_nbq_w, type = "sac", zero.policy = TRUE)
summary(LSOA_IMDModDeprivedS, Nagelkerke=T)



#ANOVA for IMD
#Avilability: AvailComN
#Accessibility: AccessS
#Visibility: GVILSOA
# Composite: Compositem

cor(CompositeLSOA$Compositem, CompositeLSOA$IMD_Decile, method = c("spearman"))

cor.test(CompositeLSOA$Compositem, CompositeLSOA$IMD_Decile, 
         method = "spearman")




Comanova <- aov(Compositem ~ IMD_Decile, data = CompositeLSOA)

summary(Comanova)

TukeyHSD(aov(CompositeLSOA$Compositem ~ as.factor(CompositeLSOA$IMD_Decile)))




Availanova <- aov(AvailComN ~ IMD_Decile, data = CompositeLSOA)

summary(Availanova)

TukeyHSD(aov(CompositeLSOA$AvailComN ~ as.factor(CompositeLSOA$IMD_Decile)))


Accessanova <- aov(AccessS ~ IMD_Decile, data = CompositeLSOA)

summary(Accessanova)

TukeyHSD(aov(CompositeLSOA$AccessS ~ as.factor(CompositeLSOA$IMD_Decile)))


Visanova <- aov(GVILSOA ~ IMD_Decile, data = CompositeLSOA)

summary(Visanova)

TukeyHSD(aov(CompositeLSOA$GVILSOA ~ as.factor(CompositeLSOA$IMD_Decile)))




#exporting images

tiff("test.tiff", units="in", width=10, height=10, res=300)

ggplot(CompositeLSOA) +
  aes(x = Compositem, y = Yrpotlife) +
  xlab("CGEI") + ylab("YPLL") +
  geom_point(size = 1, colour = "#0c4c8a") +
  geom_smooth(span = 0.5) +
  theme_minimal() +
  facet_wrap(vars(IMD_Decile))
dev.off()







#tes GEE

library(geepack)
library(gee)
library(spind)

modelGee <- geeglm (Yrpotlife~ Compositem + Income_Sco + no2_mean_a, id = UNID,
                 family = gaussian("identity"), data= CompositeLSOA, corstr = "exchangeable")
summary(modelGee)

moran.mc(modelGee$residuals, compositeLSOA_nbq_w, 999)

