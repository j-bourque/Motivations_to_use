## Motivations for cannabis and alcohol use in adolescents prone to psychosis ##
## Bourque J, Livet A, Afzali MH, Castellanos-Ryan N, Conrod PJ ##


## GLMM model to explore the association between using to cope for depression and the continuous psychosis risk 
## while also controlling for the association between using to cope for depression and other mental health outcomes.

## To do so, the dependent variable is now frequency of using alcohol or cannabis to cope with depression at every
## assessments, and the independent variables are the continuous psychosis risk, internalizing and externalizing
## symptoms, as well as substance use problems. 

#####################
## LOAD PACKAGES ####
#####################
library(dplyr)
library(caret)
library(foreign)
library(lme4)
library(data.table)
library(nlme)
library(car)
library(MASS)
library(fitdistrplus)
library(predictmeans)
library(lmerTest)
library(optimx)
library(numDeriv)
library(MuMIn)

#####################
## ORGANIZE DATA ####
#####################
# At this point the dataframe should still be in a wide format in order to create the between-person effects. 
# For reasons of non-convergence of GLMM models when including between-person, concurrent within-person, and lagged
# within-person effects of the 4 different vulnerability to psychopathology (12 predictors + 2 covariates); and 
# considering that the most consistent effects reported in aims 2 and 3 are at the between-level, we chose to only
# inlcude between-person effects in the GLMM models.

## 1. CREATE BETWEEN-PERSON VARIABLES OR MEANS ##
#Mean of all timepoints for an individual on a specific variable
mean.n   <- function(coventure_all, n) {
  means <- apply(as.matrix(coventure_all), 1, mean, na.rm = TRUE)
  nvalid <- apply(as.matrix(coventure_all), 1, function(coventure_all) sum(!is.na(coventure_all)))
  ifelse(nvalid >= n, means, NA)
}

# Creating the between-person variables of the different vulnerability to mental health problems. We chose a 
# minimum of 2 timepoints with available data to create these between-person variables.
coventure_all$PLE_mean <- mean.n(coventure_all[c(1:5)], 2)
coventure_all$SDQ_int_mean <- mean.n(coventure_all[c(11,14,17,20,23)], 2)
coventure_all$SDQ_ext_mean <- mean.n(coventure_all[c(12,15,18,21,24)], 2)
coventure_all$SU_mean <- mean.n(coventure_all[c(77:81)], 2)


## 2. RESHAPE DATA INTO LONG FORMAT ##
coventure.long <- reshape(coventure_all, varying=c("DMQ_dep_Y1","DMQ_dep_Y2","DMQ_dep_Y3","DMQ_dep_Y4","DMQ_dep_Y5",
                                                   "CMQ_dep_Y1","CMQ_dep_Y2","CMQ_dep_Y3","CMQ_dep_Y4","CMQ_dep_Y5"), 
                          idvar = "ID", direction="long", 
                          sep="_Y", timevar="time")


#######################
## BUILD GLMM MODELS ##
#######################
# Considering the dependent variables here consisted of non-normal data (range 0 to 4 with and excess of 0), and 
# that Box-Cox transformation did not improve so much the normality of residuals, we chose to transform the dependent
# variables into a binary variable (0 = 0, 1 = score > 0)

# Predicting using alcohol to cope for depression

# Abbreviations: DMQ_dep_bin is the frequency of using alcohol to cope for depression at every assessment dichotomized
# (0 or 1), DEM_01_Y1 is sex (M or F), AgeAtTesting_Y1 is age at the baseline assessment, PLE_mean is the psychosis
# risk between-person predictor, SDQ_int_mean is the internalizing symptoms between-person predictor, SDQ_ext_mean 
# is the externalizing symptoms between-person predictor, SU_mean is the substance use problems between-person
# predictor.

DMQ_dep_bet.glme <- glmer(DMQ_dep_bin ~ DEM_01_Y1+AgeAtTesting_Y1+PLE_mean+SDQ_int_mean+SDQ_ext_mean+SU_mean+
                            (1|ID), family=binomial(link="logit"), data=coventure.long_lagged,
                          control=glmerControl(optimizer="bobyqa",
                                               optCtrl=list(maxfun=2e5)))

summary(DMQ_dep_bet.glme)
se_DMQ<-sqrt(diag(vcov(DMQ_dep_bet.glme)))
#Get CI instead of SE
(tab_DMQ <- cbind(Est = fixef(DMQ_dep_bet.glme), LL = fixef(DMQ_dep_bet.glme) - 1.96 * se_DMQ, UL = fixef(DMQ_dep_bet.glme) 
                  + 1.96 * se_DMQ))
#Get odds ratios instead of coefficients
exp(tab_DMQ) # results reported in manuscript

# Predicting using cannabis to cope for depression
CMQ_dep_bet.glme <- glmer(CMQ_dep_bin ~ DEM_01_Y1+AgeAtTesting_Y1+PLE_mean+SDQ_int_mean+SDQ_ext_mean+SU_mean+
                            (1|ID), family=binomial(link="logit"), data=coventure.long_lagged)

summary(CMQ_dep_bet.glme)
se_CMQ<-sqrt(diag(vcov(CMQ_dep_bet.glme)))
#Get CI instead of SE
(tab_CMQ <- cbind(Est = fixef(CMQ_dep_bet.glme), LL = fixef(CMQ_dep_bet.glme) - 1.96 * se_CMQ, UL = fixef(CMQ_dep_bet.glme) 
                  + 1.96 * se_CMQ))
#Get odds ratios instead of coefficients
exp(tab_CMQ) # results reported in manuscript


