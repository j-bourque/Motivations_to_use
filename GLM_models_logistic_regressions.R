## Motivations for cannabis and alcohol use in adolescents prone to psychosis 
## Bourque J, Livet A, Afzali MH, Castellanos-Ryan N, Conrod PJ 


## GLM models (logistic regressions) for the first aim of the manuscript ##
## AIM 1: Investigate whether specific motives to use cannabis and alcohol discriminate youths at risk for psychosis
## from their low risk peers. This cross-sectional analysis used data from the last assessment considering the
## highest prevalence of cannabis and alcohol users were observed at the last assessment (i.e., 16.8 years old)

#####################
## LOAD PACKAGES ####
#####################
library(nnet)
library(foreign)
library(car)
library(ICC)
library(psych)
library(miceadds)
library(sandwich)

########################
## BUILD GLM MODELS ####
########################
# Logistic regression predicting at-risk from low risk youths using all 5 subscales from the Modified Drinking
# Motives Questionnaire - Revised

# Abbreviations: PSYCHOTIC_Abn_Y5 is the grouping variable (at-risk vs. low risk) based on a previously validated 
# cutoff (please refer to manuscript), DMQ_soc_Y5 is the reported frequency of using alcohol for social motives at 
# the last assessment, DMQ_enh_Y5 is the reported frequency of using alcohol for enhancement motives at the last 
# assessment, DMQ_conf_Y5 is the reported frequency of using alcohol for conformity motives at the last assessment,
# DMQ_anx_Y5 is the reported frequency of using alcohol to cope with anxiety at the last assessment, DMQ_dep_Y5 is
# the reported frequency of using alcohol to cope with depression at the last assessment, AgeAtTesting_Y5, is age at
# last assessment, DEM_01_Y5 is sex (M or F).

motives_alc_Y5.glm = glm(PSYCHOTIC_Abn_Y5 ~ DMQ_soc_Y5 + DMQ_enh_Y5 + DMQ_conf_Y5 + DMQ_anx_Y5 + DMQ_dep_Y5 +
                           AgeAtTesting_Y5+DEM_01_Y5, family=quasibinomial(link="logit"), 
                         na.action = na.omit, data=coventure_all)

motives_alc_Y5.glm_sum = summary(motives_alc_Y5.glm)
motives_alc_Y5.glm_sum

#Multicolinearity
vif(motives_alc_Y5.glm)

#Get the CI
confint.default(motives_alc_Y5.glm)

#Get the exponentiated coefficients
exp(coef(motives_alc_Y5.glm))
exp(cbind(OR = coef(motives_alc_Y5.glm), confint.default(motives_alc_Y5.glm)))

# Logistic regression predicting at-risk from low risk youths using all 5 subscales from the Cannabis Motives 
# Questionnaire
motives_can_Y5.glm = glm(PSYCHOTIC_Abn_Y5 ~ CMQ_soc_Y5 + CMQ_enh_Y5 + CMQ_conf_Y5 + CMQ_anx_Y5 + CMQ_dep_Y5 +
                           DEM_01_Y5+AgeAtTesting_Y5, family = quasibinomial(link="logit"), na.action = na.omit, 
                         data=coventure_all)

motives_can_Y5.glm_sum = summary(motives_can_Y5.glm)
motives_can_Y5.glm_sum

vif(motives_can_Y5.glm)

#Get the CI
confint.default(motives_can_Y5.glm)

#Get the exponentiated coefficients
exp(coef(motives_can_Y5.glm))
exp(cbind(OR = coef(motives_can_Y5.glm), confint.default(motives_can_Y5.glm)))



#########################################################################
## SENSITIVITY ANALYSES: GLM MODELS WITH SCHOOLS AS CLUSTER-VARIABLE ####
#########################################################################
# Logistic regression predicting at-risk from low risk youths using all 5 subscales from the Modified Drinking 
# Motives Questionnaire - Revised
DMQ_school <-miceadds::glm.cluster(data=coventure_all, formula=PSYCHOTIC_Abn_Y5 ~ DMQ_soc_Y5 + 
                                     DMQ_enh_Y5+DMQ_conf_Y5+DMQ_anx_Y5+DMQ_dep_Y5+AgeAtTesting_Y5+DEM_01_Y5,
                                   cluster="School_Y5", family="binomial")

summary(DMQ_school)
coef(DMQ_school)

#Get the CI
confint.default(DMQ_school)

#Get the exponentiated coefficients
DMQ_cluster_exp<-exp(coef(DMQ_school))
DMQ_cluster_exp
exp(cbind(OR = coef(DMQ_school), confint.default(DMQ_school)))

# Logistic regression predicting at-risk from low risk youths using all 5 subscales from the Cannabis Motives 
# Questionnaire
CMQ_school <-miceadds::glm.cluster(data=coventure_all, formula=PSYCHOTIC_Abn_Y5 ~ CMQ_soc_Y5 + 
                                     CMQ_enh_Y5+CMQ_conf_Y5+CMQ_anx_Y5+CMQ_dep_Y5+AgeAtTesting_Y5+DEM_01_Y5,
                                   cluster="School_Y5", family="binomial")

summary(CMQ_school)
coef(CMQ_school)

#Get the CI
confint.default(DMQ_school)

#Get the exponentiated coefficients
CMQ_cluster_exp<-exp(coef(CMQ_school))
CMQ_cluster_exp
exp(cbind(OR = coef(CMQ_school), confint.default(CMQ_school)))



