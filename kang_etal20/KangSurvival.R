######################################
# Survival analysis on Galleria data #
###################################### 

# check R version
version

data<-read.table("KangGalleriaData.tsv",header=T)

# note: here death = time_last_seen_alive (in hours), experiments = replicates (which came from different Galleria batches)

# Including controls, 525 Galleria, 75 for each of 5 treatments and 2 controls
nrow(data)
table(data$treatment)

# Mean survival time (hrs) - 50C had lowest survival, then other fungal treatments, then controls
sort(tapply(data$death,data$treatment,mean))

# Drop the controls ie reduce to a dataset of 375 Galleria (75 for each of 5 treatments)
data2<-data[data$treatment!="PBS"&data$treatment!="37",]
nrow(data2)
attach(data2)

treatment<-droplevels(treatment)
table(treatment)
sort(tapply(death,treatment,mean))
#       50       Fe       MM       Zn     NaCl 
# 22.41333 29.61333 31.05333 31.85333 32.48000 

########################################################################
# Set up the GLMM to estimate differences between treatments (Table 3) #
########################################################################

library(lme4)
glmmodel<-glmer(death~treatment+(1|experiment),family=Gamma)
summary(glmmodel)
#Fixed effects:
#               Estimate Std. Error t value Pr(>|z|)    
#(Intercept)    0.046545   0.005073   9.176  < 2e-16 ***
#treatmentFe   -0.010486   0.003600  -2.913 0.003584 ** 
#treatmentMM   -0.011983   0.003539  -3.386 0.000709 ***
#treatmentNaCl -0.013331   0.003486  -3.824 0.000131 ***
#treatmentZn   -0.012754   0.003508  -3.635 0.000278 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Note that releveling so that summary shows results for exp2 and Fe first does not change the results of the model
ftreat<-relevel(treatment,"Fe")
expB<-relevel(experiment,"exp2")
glmmodel2<-glmer(death~ftreat+(1|expB),family=Gamma)
summary(glmmodel2)

# After accounting for random host batch effects, there are differences in survival times among the 5 spore treatments 
# (Likelihood ratio test, Chisq=19.1, df=4, P=0.0008). 
drop1(glmmodel,test="Chi")

# Get predicted means for a table of different survival times
tapply(predict(glmmodel,type="response"),treatment,mean)
#       50       Fe       MM     NaCl       Zn 
# 22.39427 29.50834 30.92310 32.32158 31.70734 

# If we round to the nearest hour, then the model predictions and the observed averages are the same:
#       50       Fe       MM     NaCl       Zn 
# 	22 	 30 	  31	 32 	    32


# These stay the same regardless of input order for treatment or expts 
tapply(predict(glmmodel2,type="response"),ftreat,mean)


#Despite substantial differences among experimental replicates in overall Galeria survival times.
# The 50C spore treatment always resulted in lowest survival times:
tapply(predict(glmmodel,type="response"),list(experiment,treatment),mean)

# complicated, but you could figure out the colors for normalised shading that would visually show 50 is always the lowest value
#           50       Fe       MM     NaCl       Zn
#exp1 24.42905 32.84141 34.54011 36.22650 35.48494
#exp2 19.56349 24.61229 25.55414 26.46564 26.06766
#exp3 17.49153 21.42016 22.13002 22.81035 22.51410
#exp4 24.89687 33.69252 35.48280 37.26488 36.48066
#exp5 25.59044 34.97534 36.90845 38.84050 37.98932

# Model checking: errors look well-spaced and random
plot(glmmodel)


####################################################
## CONTROL ANALYSIS WITH A COX MIXED EFFECTS MODEL #
####################################################

# Using coxme version 2.2-16 
# used in the survival analysis of Chrostek et al 2013 (Luis Teixeira).. 
# cox proportional hazards corrects for the fact that not all Af-injected Galleria were dead at the end of the experiment
# though actually this is probably not much of a problem because 93% were dead at the end
library(coxme)
table(status)
348/(348+27)

# without the random effect of experiment
coxT<-coxph(Surv(death,status)~treatment)
# with the random effect of experiment
coxTE<-coxme(Surv(death,status)~treatment + (1|experiment))

# There is a significant effect of experiment (probably Galleria batch) Likelihood ratio test: df =1, 2.127e-10
# this approach follows the example by the authors of coxme
anova(coxT,coxTE)

# Likelihood ratio test shows a significant effect of treatment. df=4, P=0.002619
coxTE2<-coxme(Surv(death,status)~1 + (1|experiment))
anova(coxTE,coxTE2)

# Note that the conversion to a response feature of the predict function doesn't work with coxme (only lp or risk - not really what I want):
tapply(predict(coxTE,type="response"),treatment,mean)

summary(coxTE)

# It is not easy to extract model predictions from coxme (for showing differences among treatments), so coxme will be the control analysis
# In the coxme.pdf documentation it says:
# "For any model which can be fit by both lmekin and lmer, the latter routine would normally be prefered due to a much wider selection of post-fit tools for residuals, prediction and plotting."

