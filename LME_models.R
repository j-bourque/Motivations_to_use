## Motivations for cannabis and alcohol use in adolescents prone to psychosis 
## Bourque J, Livet A, Afzali MH, Castellanos-Ryan N, Conrod PJ 


## LME models for the second and third aim of the manuscript ##
## AIM 2: Exploring the relationships between year-to-year changes in discriminant motives to use and a continuous 
## measure of psychosis risk

## AIM 3: Exploring the relationships between year-to-year changes in discriminant motives to use and other mental
## health outcomes

#####################
## LOAD PACKAGES ####
#####################
library(dplyr)
library(caret)
library(foreign)
library(lme4)
library(data.table)
library(car)
library(MASS)
library(fitdistrplus)
library(lmerTest)
library(optimx)
library(numDeriv)
library(MuMIn)

#####################
## ORGANIZE DATA ####
#####################
# At this point the dataframe should still be in a wide format in order to create the between-person, concurrent
# within-person, and lagged within-person effects.

## 1. CREATE BETWEEN-PERSON VARIABLES OR MEANS ##
# Mean of all timepoints for an individual on a specific variable
mean.n   <- function(coventure_all, n) {
  means <- apply(as.matrix(coventure_all), 1, mean, na.rm = TRUE)
  nvalid <- apply(as.matrix(coventure_all), 1, function(coventure_all) sum(!is.na(coventure_all)))
  ifelse(nvalid >= n, means, NA)
}

#Creating the between-person variables of the the coping for depression and social motives subscales from the
#Modified Drinking Motives Questionnaire - Revised (DMQ). We chose a minimum of 2 timepoints with available data to
#create these between-person variables.
coventure_all$DMQ_dep_mean <- mean.n(coventure_all[c(29,34,39,44,49)], 2)
coventure_all$DMQ_soc_mean <- mean.n(coventure_all[c(26,31,36,41,46)], 2)

#Creating the between-person variables of the the coping for depression and social motives subscales from the
#Cannabis Motives Questionnaire (CMQ). We chose a minimum of 2 timepoints with available data to create these 
#between-person variables.
coventure_all$CMQ_dep_mean <- mean.n(coventure_all[c(54,59,64,69,74)], 2)
coventure_all$CMQ_soc_mean <- mean.n(coventure_all[c(51,56,61,66,71)], 2)


## 2. PERSON-MEAN CENTERING OF THE TIME-VARYING VARIABLES - CREATING CONCURRENT WITHIN-SUBJECT VARIABLES ####
coventure_all = coventure_all %>%
  mutate(DMQ_depm_Y1 = DMQ_dep_Y1-DMQ_dep_mean) %>%
  mutate(DMQ_depm_Y2 = DMQ_dep_Y2-DMQ_dep_mean) %>%
  mutate(DMQ_depm_Y3 = DMQ_dep_Y3-DMQ_dep_mean) %>%
  mutate(DMQ_depm_Y4 = DMQ_dep_Y4-DMQ_dep_mean) %>%
  mutate(DMQ_depm_Y5 = DMQ_dep_Y5-DMQ_dep_mean) %>%
  
  mutate(DMQ_socm_Y1 = DMQ_soc_Y1-DMQ_soc_mean) %>%
  mutate(DMQ_socm_Y2 = DMQ_soc_Y2-DMQ_soc_mean) %>%
  mutate(DMQ_socm_Y3 = DMQ_soc_Y3-DMQ_soc_mean) %>%
  mutate(DMQ_socm_Y4 = DMQ_soc_Y4-DMQ_soc_mean) %>%
  mutate(DMQ_socm_Y5 = DMQ_soc_Y5-DMQ_soc_mean)

coventure_all = coventure_all %>%
  mutate(CMQ_depm_Y1 = CMQ_dep_Y1-CMQ_dep_mean) %>%
  mutate(CMQ_depm_Y2 = CMQ_dep_Y2-CMQ_dep_mean) %>%
  mutate(CMQ_depm_Y3 = CMQ_dep_Y3-CMQ_dep_mean) %>%
  mutate(CMQ_depm_Y4 = CMQ_dep_Y4-CMQ_dep_mean) %>%
  mutate(CMQ_depm_Y5 = CMQ_dep_Y5-CMQ_dep_mean) %>%
  
  mutate(CMQ_socm_Y1 = CMQ_soc_Y1-CMQ_soc_mean) %>%
  mutate(CMQ_socm_Y2 = CMQ_soc_Y2-CMQ_soc_mean) %>%
  mutate(CMQ_socm_Y3 = CMQ_soc_Y3-CMQ_soc_mean) %>%
  mutate(CMQ_socm_Y4 = CMQ_soc_Y4-CMQ_soc_mean) %>%
  mutate(CMQ_socm_Y5 = CMQ_soc_Y5-CMQ_soc_mean)

## 3. RESHAPE DATA INTO LONG FORMAT ##
coventure.long <- reshape(coventure_all, varying=c("PSYCHOTIC_Scoretotal_Y1","PSYCHOTIC_Scoretotal_Y2",
                                                   "PSYCHOTIC_Scoretotal_Y3","PSYCHOTIC_Scoretotal_Y4",
                                                   "PSYCHOTIC_Scoretotal_Y5",
                                                   "SDQ_int_Y1","SDQ_int_Y2","SDQ_int_Y3","SDQ_int_Y4",
                                                   "SDQ_int_Y5",
                                                   "SDQ_ext_Y1","SDQ_ext_Y2","SDQ_ext_Y3","SDQ_ext_Y4",
                                                   "SDQ_ext_Y5",
                                                   "DEPAPO_ECHELLE_TOTAL_Y1","DEPAPO_ECHELLE_TOTAL_Y2",
                                                   "DEPAPO_ECHELLE_TOTAL_Y3","DEPAPO_ECHELLE_TOTAL_Y4",
                                                   "DEPAPO_ECHELLE_TOTAL_Y5",
                                                   "DMQ_depm_Y1","DMQ_depm_Y2","DMQ_depm_Y3","DMQ_depm_Y4",
                                                   "DMQ_depm_Y5","DMQ_socm_Y1","DMQ_socm_Y2","DMQ_socm_Y3",
                                                   "DMQ_socm_Y4","DMQ_socm_Y5",
                                                   "CMQ_depm_Y1","CMQ_depm_Y2","CMQ_depm_Y3","CMQ_depm_Y4",
                                                   "CMQ_depm_Y5","CMQ_socm_Y1","CMQ_socm_Y2","CMQ_socm_Y3",
                                                   "CMQ_socm_Y4","CMQ_socm_Y5"), idvar = "ID", direction="long", 
                          sep="_Y", timevar="time")

## 4. CREATE THE LAGGED WITHIN-SUBJECT VARIABLES ##
coventure.long_sorted <- coventure.long[order(coventure.long$ID) , ]

lg <- function(x)c(NA, x[1:(length(x)-1)])
coventure.long.ddt <- data.table(coventure.long_sorted)
coventure.long.ddt[,lag_DMQ_dep :=lg(DMQ_depm), by=c("ID")]
coventure.long.ddt[,lag_DMQ_soc :=lg(DMQ_socm), by=c("ID")]

coventure.long.ddt[,lag_CMQ_dep :=lg(CMQ_depm), by=c("ID")]
coventure.long.ddt[,lag_CMQ_soc :=lg(CMQ_socm), by=c("ID")]

coventure.long_lagged <- as.data.frame(coventure.long.ddt)

###################################################################
## BUILD LME MODELS - AIM 2: PREDICT CONTINUOUS PSYCHOSIS RISK ####
###################################################################
# A first alcohol motives model to visually inspect residual plots of the outcome variable

# Abbreviations: DEM_01_Y1 is sex (M or F), AgeAtTesting_Y1 is age at baseline, DMQ_soc_mean is the time-invariant
# between-person predictor of the social motive for alcohol use, DMQ_dep_mean is the time-invariant between-person 
# predictor of the coping for depression motive for alcohol use, DMQ_socm is the time-varying concurrent within-
# person predictor of the social motive for alcohol use, DMQ_depm is the time-varying concurrent within-person 
# predictor of the coping for depression motive for alcohol use, lag_DMQ_soc is the time-varying lagged within-
# person predictor of the social motive for alcohol use, lag_DMQ_dep is the time-varying lagged within-person 
# predictor of the coping for depression motive for alcohol use

PLE_DMQ.lme <- lmer(PSYCHOTIC_Scoretotal ~ DEM_01_Y1+AgeAtTesting_Y1+DMQ_soc_mean+DMQ_dep_mean+
                      DMQ_socm+DMQ_depm+lag_DMQ_dep+lag_DMQ_soc+
                      (1|ID), data=coventure.long_lagged, REML=T)

plot(PLE_DMQ.lme)

# A first cannabis motives model to visually inspect residual plots of the outcome variable
PLE_CMQ.lme <- lmer(PSYCHOTIC_Scoretotal ~ DEM_01_Y1+AgeAtTesting_Y1+CMQ_soc_mean+CMQ_dep_mean+
                      CMQ_socm+CMQ_depm+lag_CMQ_dep+lag_CMQ_soc+
                      (1|ID), data=coventure.long_lagged, REML=T)

plot(PLE_CMQ.lme)

# Try Box-Cox transformation of dependent variable and use lmer. For the Box-Cox transformation to work, we can't
# have 0 values on a dependent variable. Consequently, just add a constant. And remove NAs from the dependent 
# variable
coventure.long_lagged_no_NA <- coventure.long_lagged %>%
  mutate(PSYCHOTIC_Scoretotal_C = PSYCHOTIC_Scoretotal+0.1) %>%
  filter(!is.na(PSYCHOTIC_Scoretotal_C))
  
PLE_BC <- caret::BoxCoxTrans(coventure.long_lagged_no_NA$PSYCHOTIC_Scoretotal_C)
print(PLE_BC)

coventure.long_lagged_no_NA<-cbind(coventure.long_lagged_no_NA, PLEs_BC=predict(PLE_BC,coventure.long_lagged_no_NA$PSYCHOTIC_Scoretotal_C))


# A second alcohol motives model with the Box-Cox transformation as the dependent variable
PLE_DMQ_bc.lme <- lmer(PLEs_BC ~ DEM_01_Y1+AgeAtTesting_Y1+DMQ_soc_mean+DMQ_dep_mean+
                      DMQ_socm+DMQ_depm+lag_DMQ_dep+lag_DMQ_soc+
                      (1|ID), data=coventure.long_lagged_no_NA, REML=T)

plot(PLE_DMQ_bc.lme) ## the residuals do not seem to deviate from normality
qqnorm(resid(PLE_DMQ_bc.lme)) 
summary(PLE_DMQ_bc.lme) ## Results reported in the manuscript

# A second cannabis motives model with the Box-Cox transformation as the dependent variable
PLE_CMQ_bc.lme <- lmer(PLEs_BC ~ DEM_01_Y1+AgeAtTesting_Y1+CMQ_soc_mean+CMQ_dep_mean+
                         CMQ_socm+CMQ_depm+lag_CMQ_dep+lag_CMQ_soc+
                         (1|ID), data=coventure.long_lagged_no_NA, REML=T)

plot(PLE_CMQ_bc.lme) ## the residuals do not seem to deviate from normality
qqnorm(resid(PLE_CMQ_bc.lme)) 
summary(PLE_CMQ_bc.lme)  ## Results reported in the manuscript

#################################################################
## BUILD LME MODELS - AIM 3a: PREDICT INTERNALIZING SYMPTOMS ####
#################################################################
# A first alcohol motives model to visually inspect residual plots of the outcome variable

SDQ_int_DMQ.lme <- lmer(SDQ_int ~ DEM_01_Y1+AgeAtTesting_Y1+DMQ_soc_mean+DMQ_dep_mean+
                      DMQ_socm+DMQ_depm+lag_DMQ_dep+lag_DMQ_soc+
                      (1|ID), data=coventure.long_lagged, REML=T)

plot(SDQ_int_DMQ.lme)
qqnorm(resid(SDQ_int_DMQ.lme))
summary(SDQ_int_DMQ.lme)

# A first cannabis motives model to visually inspect residual plots of the outcome variables
SDQ_int_CMQ.lme <- lmer(SDQ_int ~ DEM_01_Y1+AgeAtTesting_Y1+CMQ_soc_mean+CMQ_dep_mean+
                      CMQ_socm+CMQ_depm+lag_CMQ_dep+lag_CMQ_soc+
                      (1|ID), data=coventure.long_lagged, REML=T)

plot(SDQ_int_CMQ.lme)
qqnorm(resid(SDQ_int_CMQ.lme))
summary(SDQ_int_CMQ.lme)

# Try Box-Cox transformation of dependent variable and use lmer. For the Box-Cox transformation to work, we can't
# have 0 values on a dependent variable. Consequently, just add a constant. And remove NAs from the dependent 
# variable
test_SDQ <- coventure.long_lagged %>%
  mutate(SDQ_int_C = SDQ_int+0.1) %>%
  filter(!is.na(SDQ_int_C))

SDQ_int_BC <- caret::BoxCoxTrans(test_SDQ$SDQ_int_C)
print(SDQ_int_BC)

test_SDQ<-cbind(test_SDQ, SDQ_int_bc=predict(SDQ_int_BC,test_SDQ$SDQ_int_C))

# A second alcohol motives model with the Box-Cox transformation as the dependent variable
SDQ_int_DMQ_bc.lme <- lmer(SDQ_int_bc ~ DEM_01_Y1+AgeAtTesting_Y1+DMQ_soc_mean+DMQ_dep_mean+
                          DMQ_socm+DMQ_depm+lag_DMQ_dep+lag_DMQ_soc+
                          (1|ID), data=test_SDQ, REML=T)

plot(SDQ_int_DMQ_bc.lme) ## the residuals do not seem to deviate from normality
qqnorm(resid(SDQ_int_DMQ_bc.lme)) 
summary(SDQ_int_DMQ_bc.lme) ## Results reported in the manuscript

# A second cannabis motives model with the Box-Cox transformation as the dependent variable
SDQ_int_CMQ_bc.lme <- lmer(SDQ_int_bc ~ DEM_01_Y1+AgeAtTesting_Y1+CMQ_soc_mean+CMQ_dep_mean+
                             CMQ_socm+CMQ_depm+lag_CMQ_dep+lag_CMQ_soc+
                             (1|ID), data=test_SDQ, REML=T)

plot(SDQ_int_CMQ_bc.lme) ## the residuals do not seem to deviate from normality
qqnorm(resid(SDQ_int_CMQ_bc.lme))  
summary(SDQ_int_CMQ_bc.lme) ## Results reported in the manuscript

#################################################################
## BUILD LME MODELS - AIM 3b: PREDICT EXTERNALIZING SYMPTOMS ####
#################################################################
# A first alcohol motives model to visually inspect residual plots of the outcome variable

SDQ_ext_DMQ.lme <- lmer(SDQ_ext ~ DEM_01_Y1+AgeAtTesting_Y1+DMQ_soc_mean+DMQ_dep_mean+
                          DMQ_socm+DMQ_depm+lag_DMQ_dep+lag_DMQ_soc+
                          (1|ID), data=coventure.long_lagged, REML=T)

plot(SDQ_ext_DMQ.lme) ## the residuals do not deviate from normality, no need to transform the dependent variable
qqnorm(resid(SDQ_ext_DMQ.lme)) 
summary(SDQ_ext_DMQ.lme) ## Results reported in the mansucript


# A first cannabis motives model to visually inspect residual plots of the outcome variable
SDQ_ext_CMQ.lme <- lmer(SDQ_ext ~ DEM_01_Y1+AgeAtTesting_Y1+CMQ_soc_mean+CMQ_dep_mean+
                          CMQ_socm+CMQ_depm+lag_CMQ_dep+lag_CMQ_soc+
                          (1|ID), data=coventure.long_lagged, REML=T)

plot(SDQ_ext_CMQ.lme) ## the residuals do not deviate from normality, no need to transform the dependent variable
qqnorm(resid(SDQ_ext_CMQ.lme))
summary(SDQ_ext_CMQ.lme) ## Results reported in the mansucript

#################################################################
## BUILD LME MODELS - AIM 3c: PREDICT SUBSTANCE USE PROBLEMS ####
#################################################################
# A first alcohol motives model to visually inspect residual plots of the outcome variable
SU_DMQ.lme <- lmer(SU ~ DEM_01_Y1+AgeAtTesting_Y1+DMQ_soc_mean+DMQ_dep_mean+
                        DMQ_socm+DMQ_depm+lag_DMQ_dep+lag_DMQ_soc+
                        (1|ID), data=coventure.long_lagged, REML=T)

plot(SU_DMQ.lme) 
qqnorm(resid(SU_DMQ.lme))
summary(SU_DMQ.lme)

# A first cannabis motives model to visually inspect residual plots of the outcome variable
SU_CMQ.lme <- lmer(SU ~ DEM_01_Y1+AgeAtTesting_Y1+CMQ_soc_mean+CMQ_dep_mean+
                     CMQ_socm+CMQ_depm+lag_CMQ_dep+lag_CMQ_soc+
                     (1|ID), data=coventure.long_lagged, REML=T)

plot(SU_CMQ.lme) 
qqnorm(resid(SU_CMQ.lme))
summary(SU_CMQ.lme)

# Try Box-Cox transformation of dependent variable and use lmer. For the Box-Cox transformation to work, we can't
# have 0 values on a dependent variable. Consequently, just add a constant. And remove NAs from the dependent 
# variable
test_SU <- coventure.long_lagged %>%
  mutate(SU_C=DEPAPO_ECHELLE_TOTAL+0.1) %>%
  filter(!is.na(SU_C))

SU_BC <- caret::BoxCoxTrans(test_SU$SU_C)
print(SU_BC)

test_SU<-cbind(test_SU, SU_bc=predict(SU_BC,test_SU$SU_C))

# A second alcohol motives model with the Box-Cox transformation as the dependent variable
SU_DMQ_bc.lme <- lmer(SU_bc ~ DEM_01_Y1+AgeAtTesting_Y1+DMQ_soc_mean+DMQ_dep_mean+
                             DMQ_socm+DMQ_depm+lag_DMQ_dep+lag_DMQ_soc+
                             (1|ID), data=test_SU, REML=T)

plot(SU_DMQ_bc.lme) ## the residuals do not seem to deviate from normality
qqnorm(resid(SU_DMQ_bc.lme)) 
summary(SU_DMQ_bc.lme) ## Results reported in the manuscript

# A second cannabis motives model with the Box-Cox transformation as the dependent variable
SU_CMQ_bc.lme <- lmer(SU_bc ~ DEM_01_Y1+AgeAtTesting_Y1+CMQ_soc_mean+CMQ_dep_mean+
                        CMQ_socm+CMQ_depm+lag_CMQ_dep+lag_CMQ_soc+
                        (1|ID), data=test_SU, REML=T)

plot(SU_CMQ_bc.lme) ## the residuals do not seem to deviate from normality
qqnorm(resid(SU_CMQ_bc.lme)) 
summary(SU_CMQ_bc.lme) ## Results reported in the manuscript


######################################################################
## SENSITIVITY ANALYSES : ADDING SCHOOLS AS ANOTHER RANDOM EFFECT ####
######################################################################
## EXAMPLE FOR AIM 2: PREDICT THE CONTINUOUS PSYCHOSIS RISK ##

# LME for alcohol motives
PLE_school_DMQ_bc.lme <- lmer(PLEs_BC ~ DEM_01_Y1+AgeAtTesting_Y1+DMQ_soc_mean+DMQ_dep_mean+
                         DMQ_socm+DMQ_depm+lag_DMQ_dep+lag_DMQ_soc+
                         (1|School_Y5)+(1|ID:School_Y5), data=coventure.long_lagged_no_NA, REML=T)

plot(PLE_school_DMQ_bc.lme) ## the residuals do not seem to deviate from normality
qqnorm(resid(PLE_school_DMQ_bc.lme))
summary(PLE_school_DMQ_bc.lme) # Results reported in the supplementary materials

# LME for cannabis motives
PLE_school_CMQ_bc.lme <- lmer(PLEs_BC ~ DEM_01_Y1+AgeAtTesting_Y1+CMQ_soc_mean+CMQ_dep_mean+
                         CMQ_socm+CMQ_depm+lag_CMQ_dep+lag_CMQ_soc+
                         (1|School_Y5)+(1|ID:School_Y5), data=coventure.long_lagged_no_NA, REML=T)

plot(PLE_school_CMQ_bc.lme) ## the residuals do not seem to deviate from normality
qqnorm(resid(PLE_school_CMQ_bc.lme)) 
summary(PLE_school_CMQ_bc.lme) # Results reported in the supplementary materials

