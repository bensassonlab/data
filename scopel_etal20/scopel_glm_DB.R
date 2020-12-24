# GLM analyses and others by Douda Bensasson 
# Analysis of Aug31 data with SSD1 info, Mixed Origins and Mosaic Beer included
# excluded haploids, and 2 strains with clade assignments that we can't recapitulate
# using the latest phylogeny to go through the clade simplification 

dataAug31<-read.delim("peter_newAug31.txt",header=T,strip.white=T)
data453<-read.csv("peter_scopel453strains.csv",header=T)

attach(dataAug31)

###########################
# HETEROZYGOSITY ANALYSIS #
###########################

table(Aneuploidy_type,Ploidy)
tapply(heterozygosity[Ploidy==2],Aneuploidy_type[Ploidy==2],median,na.rm = TRUE)

#Even when controlling for ploidy effects by examining only diploids, cells with any type of aneuploidy were
# more heterozygous (median 0.05%, N=145) than euploids (median 0.04%, N=642; Wilcoxon test, P=0.0002).
# euploids 0.04% heterozygous N=642
median(heterozygosity[Ploidy==2&Aneuploidy_type=="Euploid"],na.rm=T)
length(heterozygosity[Ploidy==2&Aneuploidy_type=="Euploid"])
# aneuploids 0.05% heterozygous N=145
median(heterozygosity[Ploidy==2&Aneuploidy_type!="Euploid"],na.rm=T)
length(heterozygosity[Ploidy==2&Aneuploidy_type!="Euploid"])
wilcox.test(heterozygosity[Ploidy==2&Aneuploidy_type=="Euploid"],heterozygosity[Ploidy==2&Aneuploidy_type!="Euploid"])

table(Aneuploidy_type,chr1only)
#although the higher heterozygosity was not significant for the few strains that exclusively lost Chr 2-16 (median 0.08%, N=10; Wilcoxon test, P = 0.3). 
# loss only chr2-16 0.08% heterozygous N=10
median(heterozygosity[Ploidy==2&Aneuploidy_type=="Loss"],na.rm=T)
length(heterozygosity[Ploidy==2&Aneuploidy_type=="Loss"])
wilcox.test(heterozygosity[Ploidy==2&Aneuploidy_type=="Euploid"],heterozygosity[Ploidy==2&Aneuploidy_type=="Loss"])

#More specifically, heterozygosity was elevated for strains that exclusively gained Chr 2-16 (median 0.05%, N=108; Wilcoxon test, P = 0.0001),
# gain only chr2-16 0.08% heterozygous N=108
median(heterozygosity[Ploidy==2&Aneuploidy_type=="Gain"],na.rm=T)
length(heterozygosity[Ploidy==2&Aneuploidy_type=="Gain"])
wilcox.test(heterozygosity[Ploidy==2&Aneuploidy_type=="Euploid"],heterozygosity[Ploidy==2&Aneuploidy_type=="Gain"])

###################################################################
# Are there difference among the clades in heterozygosity levels? #
###################################################################

# gather up all the het data and clade data for diploids that are not M. Mosaic (N=641)
het<-heterozygosity[complete.cases(dataAug31$Clades)&complete.cases(dataAug31$heterozygosity)&Ploidy==2&Clades!="M. Mosaic"]
clade<-Clades[complete.cases(dataAug31$Clades)&complete.cases(dataAug31$heterozygosity)&Ploidy==2&Clades!="M. Mosaic"]
length(het)
length(clade)
table(clade)
# There are major differences that generally correlate with what we know about aneuploidy 
# with notable high heterozygosity / low aneuploidy exceptions e.g. Brazilian bioethanol, Mexican Agave, West African cocoa
tapply(het,clade,median)
par(mar=c(8,4,4,2),mfrow=c(1,1))
plot(clade,het,las=2)
par(mar=c(5,4,4,2))
hetclademodel<-lm(het~clade)
nullhetmodel<-lm(het~1)
par(mfrow=c(5,4),mar=c(1,1,1,1))
tapply(het,clade,hist)
par(mfrow=c(2,2),mar=c(5,4,4,2))
# The model checking - the qqplot looks bad - weird effects at low and high heterozygosity
# It could be because of pure homozygotes, and interspecies hybridization
plot(hetclademodel)
# There is a strong correlation between heterozygosity and clade (R^2=0.8) 
# but it could be an overestimate because of the broken model assumptions at high and low heterzygosity
anova(hetclademodel,nullhetmodel,test="F")
summary(hetclademodel)
# Given the broken normality assumption, the non-parametric Kruskal Wallis test is better
# Yes there are significant differences among clades (chi-square=356, df=25, P<1x10-16)
kruskal.test(het,clade)


##########################
## MIXED ORIGINS EXAMPLE #
##########################

# there are 72 Mixed origins strains
# most of these (38) are aneuploid
table(droplevels(Clades[Clades=="8. Mixed origin"]),Aneuploidy_type[Clades=="8. Mixed origin"])
sum(table(droplevels(Clades[Clades=="8. Mixed origin"]),Aneuploidy_type[Clades=="8. Mixed origin"]))
72-34

# 963 strains with Clade info
sum(table(Clades,Aneuploidy_type))
table(Clades,Aneuploidy_type)

# 891 strain excluding Mixed Origins
table(Clades,Aneuploidy_type,exclude="8. Mixed origin",useNA="no")
sum(table(Clades,Aneuploidy_type,exclude="8. Mixed origin",useNA="no"))

# 171 aneuploids
table(Clades[Aneuploidy_type!="Euploid"],Aneuploidy_type[Aneuploidy_type!="Euploid"],useNA="no",exclude="8. Mixed origin")
sum(table(Clades[Aneuploidy_type!="Euploid"],Aneuploidy_type[Aneuploidy_type!="Euploid"],useNA="no",exclude="8. Mixed origin"))

# Mixed origin strains are more likely aneuploid Fisher's exact test (P=1e-9)
fisher.test(matrix(c(171,38,891-171,72-38),nrow=2))

## POLYPLOIDY AND 8. MIXED ORIGIN
# all 72 of the Mixed Origins strains have polyploidy information
# 29 are polyploid

polyploidy<-Ploidy
levels(polyploidy)
polyploidy<-Ploidy
levels(polyploidy)
levels(polyploidy)[1]<-"haploid"
levels(polyploidy)[2]<-"diploid"
levels(polyploidy)[3:6]<-"polyploid"
table(polyploidy)


table(droplevels(Clades[Clades=="8. Mixed origin"]),polyploidy[Clades=="8. Mixed origin"])
sum(table(droplevels(Clades[Clades=="8. Mixed origin"]),polyploidy[Clades=="8. Mixed origin"]))

# 963 strains with Clade and polyploidy info
sum(table(Clades,polyploidy))
table(Clades,polyploidy)

# 891 strain excluding Mixed Origins
table(Clades,polyploidy,exclude="8. Mixed origin",useNA="no")
sum(table(Clades,polyploidy,exclude="8. Mixed origin",useNA="no"))

# 60 polyploids
table(Clades[polyploidy=="polyploid"],polyploidy[polyploidy=="polyploid"],useNA="no",exclude="8. Mixed origin")
sum(table(Clades[polyploidy=="polyploid"],polyploidy[polyploidy=="polyploid"],useNA="no",exclude="8. Mixed origin"))

# Mixed origin strains are more likely polyploid Fisher's exact test (P=6e-14)
fisher.test(matrix(c(60,29,891-60,72-29),nrow=2))

## HETEROZYGOSITY AND 8. MIXED ORIGIN
# heterozygosity is higher for mixed origin (P<2e-16)
summary(heterozygosity[Clades=="8. Mixed origin"])
summary(heterozygosity[Clades!="8. Mixed origin"])
wilcox.test(heterozygosity[Clades=="8. Mixed origin"],heterozygosity[Clades!="8. Mixed origin"])



###############
## GLM MODELS #
###############

# Remove HO deletions, 1 strain without aneuploidy info and Strope data from Peter et al
# Peter et al assigned aneuploidy calls and clade classifications to 962 strains; 
nrow(dataAug31)
data2<-dataAug31[complete.cases(dataAug31$Aneuploidies)&complete.cases(dataAug31$Clades),]
nrow(data2)

# 812 of these were assigned to clades and 150 were mosaics
nrow(data2[data2$Clades=="M. Mosaic",])
data3<-data2[data2$Clades!="M. Mosaic",]
nrow(data3)

# drop 86 HO deletion strains and 45 monosporic derivatives from Strope et al
nrow(data3[data3$HO.deletion!="no",])
nrow(data3[data3$Reference=="2",])
data4<-data3[data3$HO.deletion=="no"&data3$Reference!="2",]
nrow(data4)

# Only compare gain of chr2-16 and euploids
data5<-data4[data4$Aneuploidy_type=="Euploid"|data4$Aneuploidy_type=="Gain",]
nrow(data5)

# Drop haploids, CBS382 and CBS1593 because our phylogenetic analysis shows they should not be in the 7. Mosaic Beer clade and 24. Asian Islands
data6<-data5[data5$Isolate.name!="CBS382"&data5$Isolate.name!="CBS1593"&data5$Ploidy!=1,]
nrow(data6)

# this leaves analysis of 621 strains from 26 clades
nrow(data6)
table(droplevels(data6$Clades),droplevels(data6$Ploidy))
levels(droplevels(data6$Clades))

attach(data6)
table(Clades,Ploidy)
sum(table(Clades,Ploidy))
# 520 euploids, 101 with gain
table(droplevels(Aneuploidy_type))


# Set up a binary response variable (gain: TRUE or FALSE)
gain<-Aneuploidy_type=="Gain"
gain

# Set up a simple polyploidy explanatory variable (polyploidy: diploid or polyploid)
polyploidy<-Ploidy
levels(polyploidy)
polyploidy<-Ploidy
levels(polyploidy)
levels(polyploidy)[1]<-"haploid"
levels(polyploidy)[2]<-"diploid"
levels(polyploidy)[3:6]<-"polyploid"
table(polyploidy)


#####################################################
## ANALYSIS OF Clades * polyploidy * heterozygosity #
#####################################################

# Start with the full model

# Does converge 
model1<-glm(gain~Clades*polyploidy*heterozygosity,binomial)
# deviance looks good
summary(model1)
anova(model1,test="Chi")

# no signfiicant effect of Clades:polyploidy:heterozygosity but it was close (df=1, 0.05017)
model2<-update(model1, ~. -  Clades:polyploidy:heterozygosity)
anova(model1,model2,test="Chi")
anova(model2,test="Chi")
# no signfiicant effect of Clades:heterozygosity (df=24,P=0.563) 
model3<-update(model2, ~. -  Clades:heterozygosity)
anova(model2,model3,test="Chi")
anova(model3,test="Chi")
# no signfiicant effect of Clades:heterozygosity (df=7,P=0.7) 
model4<-update(model3, ~. -  Clades:polyploidy)
anova(model3,model4,test="Chi")
anova(model4,test="Chi")
# no signfiicant effect of polyploidy:heterozygosity (df=1,P=0.8187) 
model5<-update(model4, ~. - polyploidy:heterozygosity)
anova(model4,model5,test="Chi")
anova(model5,test="Chi")
# no signfiicant effect of heterozygosity (df=1,P=0.807) 
model6<-update(model5, ~. - heterozygosity)
anova(model5,model6,test="Chi")
anova(model6,test="Chi")
# no signfiicant effect of polyploidy (df=1,P=0.2501) 
model7<-update(model6, ~. - polyploidy)
anova(model6,model7,test="Chi")
# effect of clades (df=25,P=9.488e-11)
model8<-update(model7, ~. - Clades)
anova(model7,model8,test="Chi")

# explains 17.9% of deviance AIC=504.54. This simple clade model (model7) is the same as cemodel8 below.
anova(model7,test="Chi")
summary(model7)
(551.48-452.54)/551.48



##########################################################################
## ANALYSIS OF Ecological.origins * Clades * heterozygosity * polyploidy #
##########################################################################

# Start with the full model
# Does NOT converge
cemodel1<-glm(gain~Clades*Ecological.origins,binomial)
# Clades:Ecological.origins won't work
cemodel1<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Clades:Ecological.origins,binomial)
cemodel1<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Clades:heterozygosity + Ecological.origins:heterozygosity + Clades:polyploidy +  Ecological.origins:polyploidy + heterozygosity:polyploidy,binomial)
cemodel1<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Clades:heterozygosity + Ecological.origins:heterozygosity + Clades:polyploidy +  Ecological.origins:polyploidy,binomial)


# This does converge!
cemodel1<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + heterozygosity:polyploidy,binomial)
# now step simplifies to Ecological.origins and heterozygosity
step(cemodel1)

# The fullest model has 4 of the 6 possible 2-way interactions (missing is Clades:polyploidy + Clades:Ecological.origins
cemodel1<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Clades:heterozygosity + Ecological.origins:heterozygosity + Ecological.origins:polyploidy + heterozygosity:polyploidy,binomial)
anova(cemodel1,test="Chi")

# no signfiicant effect of  Ecological.origins:heterozygosity  (df=21,P= 0.7125)
cemodel2<-update(cemodel1, ~. -  Ecological.origins:heterozygosity)
anova(cemodel1,cemodel2,test="Chi")
anova(cemodel2,test="Chi")
# no signfiicant effect of  Clades:heterozygosity (df=24,P= 0.5191)
cemodel3<-update(cemodel2, ~. -  Clades:heterozygosity)
anova(cemodel2,cemodel3,test="Chi")
anova(cemodel3,test="Chi")
# no signfiicant effect of  Ecological.origins:polyploidy (df=11,P= 0.1689)
cemodel4<-update(cemodel3, ~. -  Ecological.origins:polyploidy)
anova(cemodel3,cemodel4,test="Chi")
anova(cemodel4,test="Chi")
# no signfiicant effect of  heterozygosity:polyploidy  (df=1,P=0.7352)
cemodel5<-update(cemodel4, ~. -  heterozygosity:polyploidy)
anova(cemodel4,cemodel5,test="Chi")
anova(cemodel5,test="Chi")
# no signfiicant effect of  heterozygosity  (df=1,P= 0.5465)
cemodel6<-update(cemodel5, ~. -  heterozygosity)
anova(cemodel5,cemodel6,test="Chi")
anova(cemodel6,test="Chi")
# no signfiicant effect of  polyploidy  (df=1,P= 0.1584)
# THE MODEL WITH Clades + Ecological.origins
# and explains 0.1793936 vs  0.04696816 of deviance with this input order
cemodel7<-update(cemodel6, ~. -  polyploidy)
anova(cemodel6,cemodel7,test="Chi")
summary(cemodel7)
anova(cemodel7,test="Chi")
98.932/551.48
25.902/551.48
# with the other input order, ecology explains 0.163908 vs 0.06245376 for the other input order
cemodel7ii<-glm(gain~Ecological.origins+Clades,binomial)
anova(cemodel7ii,test="Chi")
90.392/551.48
34.442/551.48
# no signfiicant effect of  Ecological.origins (df=22,P=0.2559)
cemodel8<-update(cemodel7, ~. -   Ecological.origins)
anova(cemodel7,cemodel8,test="Chi")
anova(cemodel8,test="Chi")
# effect of clades (df=25,P=9.488e-11)
cemodel8ii<-update(cemodel8, ~. -  Clades)
anova(cemodel8,cemodel8ii,test="Chi")

# dropping Clades from the model sooner:
# no signfiicant effect of  polyploidy  (df=1,P= 0.1362)
cemodel9<-update(cemodel5, ~. -  polyploidy)
anova(cemodel5,cemodel9,test="Chi")
anova(cemodel9,test="Chi")
# no signfiicant effect of  Clades  (df=25,P= 0.4769) 
cemodel10<-update(cemodel9, ~. -  Clades)
anova(cemodel9,cemodel10,test="Chi")
anova(cemodel10,test="Chi")
# effect of Ecological origins (deviance=14.75357%, df=22,P=9.653e-09)
cemodel10ii<-update(cemodel10, ~. -  Ecological.origins)
anova(cemodel10,cemodel10ii,test="Chi")
81.363/551.48
# effect of heterozygosity (deviance=1.782875%,df=1,P=0.002)
cemodel10iii<-update(cemodel10, ~. -  heterozygosity)
anova(cemodel10,cemodel10iii,test="Chi")
9.8322/551.48

# dropping clade: deviance=6.245376%, df 25, P =0.0988 and leaving in Ecological origins .. this is no longer significantly worse!
model1noclade<-glm(gain~Ecological.origins,binomial)
anova(cemodel7,model1noclade,test="Chi")
34.442/551.48

# AIC=522.64 for clade + ecology
summary(cemodel7)
# AIC=507.09 for ecology, explains 16.39044% deviance
summary(model1noclade)
(551.48-461.09)/551.48
# AIC=504.54 for Clades only (= same model as above), explains 17.9% deviancee, is not worse than ecology+clade (deviance=4.696816%,df=22,P=0.2559)
summary(cemodel8)
(551.48-452.54)/551.48
anova(cemodel7,cemodel8,test="Chi")
25.902/551.48
# AIC=499.25 for Ecological.origins + heterozygosity, explains 18.2% deviance, is not worse than ecology+clade (df=24,P=0.4272)
summary(cemodel10)
(551.48-451.25)/551.48
anova(cemodel7,cemodel10,test="Chi")

# Note that heterozygosity is a good predictor of high chr gain with notable excepions (bioethanol, agave) 
# could it also be an indicator of recent hybridization?
sort(tapply(heterozygosity,Clades,median))

"          19. Malaysian        1. Wine/European           17. Taiwanese                16. CHNI 
           0.0002166650            0.0002199695            0.0002814110            0.0002831800 
   4. Mediterranean oak       24. Asian islands               15. CHNII             2. Alpechin 
           0.0002886450            0.0003200165            0.0003221560            0.0003237630 
   22. Far East Russian               20. CHN V              14. CHNIII  23. North American oak 
           0.0003412545            0.0003437055            0.0003533330            0.0003552160 
         21. Ecuadorean       18. Far East Asia                25. Sake   13. African palm wine 
           0.0003731420            0.0003811940            0.0004018040            0.0007282975 
 26. Asian fermentation 10. French Guiana human          7. Mosaic beer  12. West African cocoa 
           0.0007467600            0.0009877790            0.0015532155            0.0017086130 
        5. French dairy         6. African beer 3. Brazilian bioethanol        9. Mexican agave 
           0.0022034600            0.0023124185            0.0025062580            0.0025498590 
           11. Ale beer         8. Mixed origin 
           0.0030242990            0.0035998945 
"

######################
# Simplifying Clades #
######################

table(droplevels(Clades))
length(table(droplevels(Clades)))
sum(table(droplevels(Clades)))
table(ssd1geno)
sum(table(ssd1geno))


# can't simplify to ssd1geno is much worse than clades (df=23,P=5.87e-10)
ssd1model<-glm(gain~ssd1geno,binomial)
anova(cemodel8,ssd1model,test="Chi")
# and it can't replace ecology (df=20,P=1.817e-07)
ssd1model2<-glm(gain~ssd1geno+heterozygosity,binomial)
anova(cemodel10,ssd1model2,test="Chi")

##########
# Using the maximum likelihood tree to decide _a priori_ contrasts

# 2. Alpechin is embedded within the 1. Wine European clade
# I can combine 1. Wine and 2. Alpechin (df=1,P=0.2576)
# I see no other clades within clades
clades2<-Clades
levels(clades2)[1]<-"Clades1_2"
levels(clades2)[12]<-"Clades1_2"
levels(clades2)
scmodel2<-glm(gain~clades2,binomial)
anova(cemodel8,scmodel2,test="Chi")

# 3. Brazilian bioethanol is the sister group and combines with Clades1_2 (df=1,P=0.6895) 
clades3<-clades2
levels(clades3)[1]<-"Clades1_2_3"
levels(clades3)[19]<-"Clades1_2_3"
levels(clades3)
scmodel3<-glm(gain~clades3,binomial)
anova(scmodel2,scmodel3,test="Chi")

# DOESN'T Simplify Here
# 7. Mosaic Beer is the sister group and is sig. diff than Clades1_2_3 (df=1,P=0.03672) 
clades4<-clades3
levels(clades4)[1]<-"Clades1to3_7"
levels(clades4)[22]<-"Clades1to3_7"
levels(clades4)
scmodel4<-glm(gain~clades4,binomial)
anova(scmodel3,scmodel4,test="Chi")

# 8. Mixed Origin and 11. Ale Beer are sister taxa and are NOT SIG DIFF (df=1,P=0.5401)
# It does not make sense to try to combine this group with Clades 1to3_7 
# because there was gain heterogeneity within that group
clades5<-clades3
levels(clades5)[23]<-"Clades8_11"
levels(clades5)[3]<-"Clades8_11"
levels(clades5)
scmodel5<-glm(gain~clades5,binomial)
anova(scmodel3,scmodel5,test="Chi")

# 5. French dairy and 6 African Beer are sister taxa and are NOT SIG DIFF (df=1,P=0.8466)
# It does not make sense to try to combine all the groups tested so far 
# because there was gain heterogeneity 
clades6<-clades5
levels(clades6)[20:21]<-"Clades5_6"
levels(clades6)
scmodel6<-glm(gain~clades6,binomial)
anova(scmodel5,scmodel6,test="Chi")


# 9. Mexican agave and 10. French Guiana human are sister taxa and are NOT SIG DIFF (df=1,P=0.1313)
# It does not make sense to try to combine all the groups tested so far 
# because there was gain heterogeneity 
# for the same reason, 12. West African cacao can't combine with anything else
clades7<-clades6
levels(clades7)[22]<-"Clades9_10"
levels(clades7)[2]<-"Clades9_10"
levels(clades7)
scmodel7<-glm(gain~clades7,binomial)
anova(scmodel6,scmodel7,test="Chi")


# 25. Sake and 26. Asian fermentation are sister taxa and are SIG DIFF (df=1,P=5.51e-07)
# It does not make sense to try to combine 24. Asian islands with 25 and 26 
# because there was gain heterogeneity 
clades8<-clades7
levels(clades8)[17:18]<-"Clades25_26"
levels(clades8)
scmodel8<-glm(gain~clades8,binomial)
anova(scmodel7,scmodel8,test="Chi")

# 23. North American oak and 22. Far East Russian are sister taxa and are NOT SIG DIFF (df=1,no P: deviance is tiny?)
# It does not make sense to try to combine with 24. Asian islands, 25 and 26 
# because there was gain heterogeneity 
# Or to bring in the next one out 13. African Palm wine
clades9<-clades7
levels(clades9)[14:15]<-"Clades22_23"
levels(clades9)
scmodel9<-glm(gain~clades9,binomial)
anova(scmodel7,scmodel9,test="Chi")

# 20. Ecuadorean and 21. CHNV are sister taxa and are NOT SIG DIFF (df=1,P=0.5154)
clades10<-clades9
levels(clades10)[12:13]<-"Clades20_21"
levels(clades10)
scmodel10<-glm(gain~clades10,binomial)
anova(scmodel9,scmodel10,test="Chi")

# 19. Malaysian and 20_21. are sister taxa and are NOT SIG DIFF (df=1,P=0.3773)
# It does not make sense to combine more groups because there was gain heterogeneity
clades11<-clades10
levels(clades11)[11:12]<-"Clades19to21"
levels(clades11)
scmodel11<-glm(gain~clades11,binomial)
anova(scmodel10,scmodel11,test="Chi")
# AIC: 493.88 ie better than the ecology + heterozygosity model now
summary(scmodel11)

# Combine all clades with the ancestral proportion of chr2-16 gains (df=13,P= 0.4908,AIC=480.33), just 5 proportions
clades12<-clades11
levels(clades12)[1:2]<-"ancestral"
levels(clades12)[3:12]<-"ancestral"
levels(clades12)[4:5]<-"ancestral"
levels(clades12)
scmodel12<-glm(gain~clades12,binomial)
anova(scmodel11,scmodel12,test="Chi")
anova(scmodel12,test="Chi")
summary(scmodel12)
(551.48-470.33)/551.48
# not signfiicantly worse than all 26 clades (deviance = - 3.22514%, df=21, P=0.6625)
anova(cemodel8,scmodel12,test="Chi")
# model predictions
#    ancestral     Clades8_11       25. Sake      Clades5_6 7. Mosaic beer      M. Mosaic 
#    0.08686441     0.38983051     0.59459459     0.26666667     0.37500000             NA 
table(droplevels(clades12),gain)
tapply(predict(scmodel12,type="response"),clades12,mean)


##########
# Using the TreeMix tree (6 edges) to decide _a priori_ contrasts

# Wine, Alpechin and Brazilian bioethanol all combine into Clades1_2_3 as above 
anova(scmodel3,test="Chi")

# 8. Mixed Origin and 11. Ale Beer combine as above
anova(scmodel3,scmodel5,test="Chi")

# 7. Mosaic Beer is the sister group and is not sig. diff than Clades8_11 (df=1,P=0.9355) 
tmclades6<-clades5
levels(tmclades6)[3]<-"Clades7_8_11"
levels(tmclades6)[22]<-"Clades7_8_11"
levels(tmclades6)
tmodel6<-glm(gain~tmclades6,binomial)
anova(scmodel5,tmodel6,test="Chi")

# Wine-Alpechin-Bioethanol (Clades1_2_3) are sig diff than the high aneuploidy group MosaicBeer-MixedOrigin-AleBeer (Clades7_8_11) (df=1,P=2.721e-08)
tmclades7<-tmclades6
levels(tmclades7)[1]<-"Clades1to3_7_8_11"
levels(tmclades7)[3]<-"Clades1to3_7_8_11"
levels(tmclades7)
tmodel7<-glm(gain~tmclades7,binomial)
anova(tmodel6,tmodel7,test="Chi")

# therefore MediteranneanOak, African_beer and French_dairy don't get combined with the others (because there would be aneuploidy rate heterogeneity

# 9. Mexican agave and 10. French Guiana human are sister taxa and are NOT SIG DIFF (df=1,P=0.1313) as above
# It does not make sense to try to combine all the groups tested so far 
# because there was gain heterogeneity 
tmclades8<-tmclades6
levels(tmclades8)[22]<-"Clades9_10"
levels(tmclades8)[2]<-"Clades9_10"
levels(tmclades8)
tmodel8<-glm(gain~tmclades8,binomial)
anova(tmodel6,tmodel8,test="Chi")


# 25. Sake and 26. Asian fermentation are sister taxa and are SIG DIFF (df=1,P=5.51e-07) as above
# It does not make sense to try to combine 24. Asian islands, North American oak or Far East Russian with 25 and 26 
# because there was gain heterogeneity 
tmclades9<-tmclades8
levels(tmclades9)[17:18]<-"Clades25_26"
levels(tmclades9)
tmodel9<-glm(gain~tmclades9,binomial)
anova(tmodel8,tmodel9,test="Chi")

# 12. West African cocoa and 13. African palm wine human are sister taxa and are NOT SIG DIFF (df=1,P=0.08587) 
# It does not make sense to try to combine all the groups tested so far 
# because there was gain heterogeneity 
tmclades10<-tmclades8
levels(tmclades10)[4:5]<-"Clades12_13"
levels(tmclades10)
tmodel10<-glm(gain~tmclades10,binomial)
anova(tmodel8,tmodel10,test="Chi")

# 20. CHN V and 21. Ecuadorean are sister taxa and are NOT SIG DIFF  (df=1,P=0.5154) as above 
# It does not make sense to try to combine all the groups tested so far 
# because there was gain heterogeneity 
tmclades11<-tmclades10
levels(tmclades11)[11:12]<-"Clades20_21"
levels(tmclades11)
tmodel11<-glm(gain~tmclades11,binomial)
anova(tmodel10,tmodel11,test="Chi")


# Combine all clades with the ancestral proportion of chr2-16 gains (df=16,P=0.1691,AIC=487.27), just 3 proportions 
tmclades12<-tmclades11
levels(tmclades12)[1:2]<-"ancestral"
levels(tmclades12)[3:13]<-"ancestral"
levels(tmclades12)[4:7]<-"ancestral" # merge in 5. French dairy and 6. French dairy now because they do not form a monophyletic group
levels(tmclades12)
tmodel12<-glm(gain~tmclades12,binomial)
anova(tmodel11,tmodel12,test="Chi")
anova(tmodel12,test="Chi")
summary(tmodel12)
# not as good as the above model for AIC or deviance explained (12.7312%)
(551.48-481.27)/551.48
# not signfiicantly worse than all 26 clades (deviance = - 3.22514%, df=23, P=0.1897)
anova(cemodel8,tmodel12,test="Chi")

# .. but it is worse than the model that recognized French Dairy and African Beer as high frequency groups (df=2,P=0.004218)
anova(scmodel12,tmodel12,test="Chi")





######################
# ECOLOGICAL ORIGINS #
######################

# note there are 23 ecological categories before and after filtering
table(Ecological.origins)
length(table(Ecological.origins))
sum(table(Ecological.origins))
table(dataAug31$Ecological.origins)
length(table(dataAug31$Ecological.origins))
sum(table(dataAug31$Ecological.origins))
table(Ecological.origins)

######################
# Simplifying Ecological.origins using the Peter et al simplification in Figure 1 and Figure S1B
# Note that they leave "Unknown" as it is, and the 2 Lab strains are in gray, so I think they stay as they are

# Strains from "Human" and "Human, clinical" are not significantly different (df=1,P=0.6154)
ecology2<-Ecological.origins
levels(ecology2)[10:11]<-"Human"
levels(ecology2)
semodel2<-glm(gain~ecology2 + heterozygosity,binomial)
anova(cemodel10,semodel2,test="Chi")

# Strains from Wild ("Tree","Nature","Fruit","Soil","Insect","Flower","Water") are not significantly different (df=6,P=0.0592)
ecology3<-ecology2
levels(ecology3)[8:9]<-"Wild"
levels(ecology3)[11]<-"Wild"
levels(ecology3)[12]<-"Wild"
levels(ecology3)[15:16]<-"Wild"
levels(ecology3)[16]<-"Wild"
levels(ecology3)
semodel3<-glm(gain~ecology3 + heterozygosity,binomial)
anova(semodel2,semodel3,test="Chi")

# But Domesticated strains ("Wine","Beer","Sake","Bakery","Fermentation","Palm wine","Industrial","Distillery","Dairy","Ethanol","Cider","Probiotic") are very different and cannot all be combined (df=11,P=4.799e-09)
ecology4<-ecology3
levels(ecology4)[1:7]<-"Domesticated"
levels(ecology4)[4]<-"Domesticated"
levels(ecology4)[5:7]<-"Domesticated"
levels(ecology4)[6]<-"Domesticated"
levels(ecology4)
semodel4<-glm(gain~ecology4 + heterozygosity,binomial)
anova(semodel3,semodel4,test="Chi")

# NB "Human" and "Wild" ecologies do not have differing prevalence of chr gain (df=1,P=0.7146)
ecology5<-ecology3
levels(ecology5)
levels(ecology5)[8:9]<-"Human_Wild"
levels(ecology5)
semodel5<-glm(gain~ecology5 + heterozygosity,binomial)
anova(semodel3,semodel5,test="Chi")

# The simplest model with Peter ecology + heterozygosity has less good AIC (497.63) than the simplest clade model (480.33)
summary(semodel3)
# The simple ecology + heterozygosity model explains 15.9% of deviance
(551.48-463.63)/551.48
anova(semodel3,test="Chi")

#The simple clade model explains 14.7% (4 df for clades) but with 12 fewer parameter estimates (16 df for ecology+het model)
summary(scmodel12)
(551.48-470.33)/551.48
anova(scmodel12,test="Chi")
anova(cemodel8,scmodel12,test="Chi")

#####################################################
## ANALYSIS OF ecology * polyploidy * heterozygosity #
#####################################################

# Start with the full model

# Does converge 
emodel1<-glm(gain~Ecological.origins*polyploidy*heterozygosity,binomial)
# deviance looks good
summary(emodel1)
anova(emodel1,test="Chi")

# no signfiicant effect of Ecological.origins:polyploidy:heterozygosity (df=7, P=0.1793)
emodel2<-update(emodel1, ~. -  Ecological.origins:polyploidy:heterozygosity)
anova(emodel1,emodel2,test="Chi")
anova(emodel2,test="Chi")

# no signfiicant effect of Ecological.origins:heterozygosity (df=21,P= 0.4893)
emodel3<-update(emodel2, ~. -  Ecological.origins:heterozygosity)
anova(emodel2,emodel3,test="Chi")
anova(emodel3,test="Chi")

# no signfiicant effect of Ecological.origins:polyploidy (df=11, P=0.2456)
emodel4<-update(emodel3, ~. -  Ecological.origins:polyploidy)
anova(emodel3,emodel4,test="Chi")
anova(emodel4,test="Chi")

# no signfiicant effect of polyploidy:heterozygosity (df=1, P=0.9127)
emodel5<-update(emodel4, ~. -   polyploidy:heterozygosity)
anova(emodel4,emodel5,test="Chi")
anova(emodel5,test="Chi")

# signfiicant effect of heterozygosity! (df=1, P=0.04203)
emodel6<-update(emodel5, ~. -   heterozygosity)
anova(emodel5,emodel6,test="Chi")

# no signfiicant effect of polyploidy (df=1, P=0.5504)
emodel7<-update(emodel5, ~. -   polyploidy)
anova(emodel5,emodel7,test="Chi")
anova(emodel7,test="Chi")

# significant effect of heterozygosity! (df=1, P=0.001715)
# Now this model with Ecological.origins and heterozygosity looks the best from AIC (499.25)
emodel8<-update(emodel7, ~. -   heterozygosity)
anova(emodel7,emodel8,test="Chi")
summary(emodel7)


#########################################################################
# Do conclusions change when looking only at 453 strains?               #
# after dropping close relatives (pairwise genetic distance < 0.000007) #
#########################################################################

table(data453$Clades,data453$Ploidy)

attach(data453)
table(Clades,Ploidy)

# Set up a binary response variable and polyploidy
gain<-Aneuploidy_type=="Gain"
gain
polyploidy<-as.factor(Ploidy)
levels(polyploidy)
levels(polyploidy)[1]<-"diploid"
levels(polyploidy)[2:4]<-"polyploid"
table(polyploidy)


# Start with the full model

# does not converge
model453i<-glm(gain~Clades*Ecological.origins,binomial)
model453i<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Clades:Ecological.origins,binomial)
model453i<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Clades:heterozygosity + Ecological.origins:heterozygosity + Clades:polyploidy +  Ecological.origins:polyploidy + heterozygosity:polyploidy + Ecological.origins:heterozygosity:polyploidy,binomial)

# does converge with the 3-way interaction Clades:heterozygosity:polyploidy and all 2-way interactions except Clades:Ecological.origins
model453i<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Clades:heterozygosity + Ecological.origins:heterozygosity + Clades:polyploidy +  Ecological.origins:polyploidy + heterozygosity:polyploidy + Clades:heterozygosity:polyploidy,binomial)
# The 3-way interaction is not significant (df=4,P=0.1527)
model453ii<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Clades:heterozygosity + Ecological.origins:heterozygosity + Clades:polyploidy +  Ecological.origins:polyploidy + heterozygosity:polyploidy,binomial)
anova(model453i,model453ii,test="Chi")
anova(model453ii,test="Chi")
# ns Clades:heterozygosity (df=23, P=0.4536)
model453iii<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Ecological.origins:heterozygosity + Clades:polyploidy +  Ecological.origins:polyploidy + heterozygosity:polyploidy,binomial)
anova(model453ii,model453iii,test="Chi")
anova(model453iii,test="Chi")
# ns Ecological.origins:heterozygosity  (df=20, P=0.1512)
model453iv<-glm(gain ~ Clades + Ecological.origins + heterozygosity + polyploidy + Clades:polyploidy +  Ecological.origins:polyploidy + heterozygosity:polyploidy,binomial)
anova(model453iii,model453iv,test="Chi")
anova(model453iv,test="Chi")
# significant Clades:polyploidy (df=4, P=0.03193)
model453v<-update(model453iv, ~. -  Clades:polyploidy)
anova(model453iv,model453v,test="Chi")
# significant Ecological.origins:polyploidy effects (df=7,P=0.007473)
model453vi<-update(model453iv, ~. -  Ecological.origins:polyploidy)
anova(model453iv,model453vi,test="Chi")
# significant heterozygosity:polyploidy effects (df=1,P=0.01169)
model453vii<-update(model453iv, ~. -  heterozygosity:polyploidy )
anova(model453iv,model453vii,test="Chi")
# The same result using step
model453istep<-step(model453i)
anova(model453istep,test="Chi")

# How about if we start with Ecology not clade? Can we drop clade then? 
# No: can't drop ecology or clade because of Ecological.origins:polyploidy and Clades:polyploidy
model453e<-glm(gain ~ Ecological.origins + heterozygosity + polyploidy + Clades + Clades:heterozygosity + Ecological.origins:heterozygosity + Clades:polyploidy +  Ecological.origins:polyploidy + heterozygosity:polyploidy + Clades:heterozygosity:polyploidy,binomial)
model453estep<-step(model453e)
anova(model453estep,test="Chi")

############################
## ANALYSING ONLY DIPLOIDS #
############################

# there are only 31 polyploids in this filtered analysis and all the significant interactions involve polyploidy 
# can we focus better on ecology vs genetic clade by dropping polyploids for this analysis of independent lineages?

data453D<-data453[data453$Ploidy==2,]
attach(data453D)
# there are only 2 Ale beer strains and 3 African beer strains left
table(Clades,Ploidy)

# Set up a binary response variable and polyploidy
gain<-Aneuploidy_type=="Gain"


# Start with the full model

# does not converge
model453Di<-glm(gain~Clades*Ecological.origins,binomial)
model453Di<-glm(gain ~ Clades + Ecological.origins + heterozygosity + Clades:Ecological.origins,binomial)

# this is the maximal model because it does converge
model453Di<-glm(gain ~ Clades + Ecological.origins + heterozygosity + Clades:heterozygosity + Ecological.origins:heterozygosity,binomial)
anova(model453Di,test="Chi")
# no effect of Clades:heterozygosity (df=23,P=0.6109)
model453Dii<-update(model453Di, ~. -  Clades:heterozygosity)
anova(model453Di,model453Dii,test="Chi")
# no effect of Ecological.origins:heterozygosity (df=20,P=0.1421)
model453Diii<-update(model453Dii, ~. -  Ecological.origins:heterozygosity)
anova(model453Dii,model453Diii,test="Chi")
anova(model453Diii,test="Chi")
# no effect of Ecological.origins (df=21,P=0.4784) - note that 1 ecology is missing
model453Div<-update(model453Diii, ~. -  Ecological.origins)
anova(model453Diii,model453Div,test="Chi")
anova(model453Div,test="Chi")
# no effect of heterozygosity (df=1,P=0.9545)
model453Dv<-update(model453Div, ~. -  heterozygosity)
anova(model453Div,model453Dv,test="Chi")
# sig diffs among clades (df=25,P=0.001381) - note that all clades are retained
anova(model453Dv,test="Chi")
# Clade only - AIC: 345.6, 14.92568% deviance
summary(model453Dv)
(345.11-293.60)/345.11


# Can we drop Clade instead of Ecology? Yes
model453D2i<-glm(gain ~ Ecological.origins + heterozygosity + Clades + Ecological.origins:heterozygosity + Clades:heterozygosity,binomial)
anova(model453D2i,test="Chi")
# no effect of Clades:heterozygosity (df=23,P=0.6109) - unchanged
model453D2ii<-update(model453D2i, ~. -  Clades:heterozygosity)
anova(model453D2i,model453D2ii,test="Chi")
# no effect of Ecological.origins:heterozygosity (df=20,P=0.1421) - unchanged
model453D2iii<-update(model453D2ii, ~. -  Ecological.origins:heterozygosity)
anova(model453D2ii,model453D2iii,test="Chi")
anova(model453D2iii,test="Chi")
# no effect of Clades (df=25,P=0.6097)
model453D2iv<-update(model453D2iii, ~. -  Clades)
anova(model453D2iii,model453D2iv,test="Chi")
anova(model453D2iv,test="Chi")
# no effect of heterozygosity (df=1,P=0.4835)
model453D2v<-update(model453D2iv, ~. -  heterozygosity)
anova(model453D2iv,model453D2v,test="Chi")
# sig diffs among ecologies (df=21,P=0.0004623)
anova(model453D2v,test="Chi")
# Ecology only - AIC: 339.85, 14.27371% deviance
summary(model453D2v)
(345.11-295.85)/345.11


# Simplifying Clade using the original phylogeny

# 2. Alpechin is embedded within the 1. Wine European clade
# I can combine 1. Wine and 2. Alpechin (df=1,P=0.2476)
# I see no other clades within clades
clades2<-Clades
levels(clades2)[1]<-"Clades1_2"
levels(clades2)[12]<-"Clades1_2"
levels(clades2)
scmodelD2<-glm(gain~clades2,binomial)
anova(model453Dv,scmodelD2,test="Chi")

# 3. Brazilian bioethanol is the sister group and combines with Clades1_2 (df=1,P=0.3715) 
clades3<-clades2
levels(clades3)[1]<-"Clades1_2_3"
levels(clades3)[19]<-"Clades1_2_3"
levels(clades3)
scmodelD3<-glm(gain~clades3,binomial)
anova(scmodelD2,scmodelD3,test="Chi")

# DOES Simplify Here
# 7. Mosaic Beer (5 euploids:2 aneuploids) is the sister group and is no longer sig. diff than Clades1_2_3, 200:23 (df=1,P=0.1892) 
table(clades3,gain)
clades4<-clades3
levels(clades4)[1]<-"Clades1to3_7"
levels(clades4)[22]<-"Clades1to3_7"
levels(clades4)
scmodelD4<-glm(gain~clades4,binomial)
anova(scmodelD3,scmodelD4,test="Chi")

# 8. Mixed origin is the sister group to 11. Ale beer and they are NOT SIG DIFF (df=1,P=0.1854)
clades5<-clades4
levels(clades5)[3]<-"Clades8_11"
levels(clades5)[22]<-"Clades8_11"
levels(clades5)
scmodelD5<-glm(gain~clades5,binomial)
anova(scmodelD4,scmodelD5,test="Chi")

# Clades8_11 are the sister group to Clades1to3_7 and they are ARE SIG DIFF (df=1,P=0.001263)
clades6<-clades5
levels(clades6)[1]<-"Clades1to3_7to8_11"
levels(clades6)[3]<-"Clades1to3_7to8_11"
levels(clades6)
scmodelD6<-glm(gain~clades6,binomial)
anova(scmodelD5,scmodelD6,test="Chi")

# 5. French dairy and 6. African beer are sister taxa and they are NOT SIG DIFF (df=1,P=0.6773)
clades7<-clades5
levels(clades7)[20:21]<-"Clades5_6"
levels(clades7)
scmodelD7<-glm(gain~clades7,binomial)
anova(scmodelD5,scmodelD7,test="Chi")

# 9. Mexican Agave and 10. French Guiana are sister taxa  and are NOT SIG DIFF (df=1,P=0.09533)
clades8<-clades7
levels(clades8)[2]<-"Clades9_10"
levels(clades8)[21]<-"Clades9_10"
levels(clades8)
scmodelD8<-glm(gain~clades8,binomial)
anova(scmodelD7,scmodelD8,test="Chi")

# 25. Sake and 26. Asian Fermentation are sister taxa  and are SIG DIFF (df=1,P=7.076e-05)
clades9<-clades8
levels(clades9)[17:18]<-"Clades25_26"
levels(clades9)
scmodelD9<-glm(gain~clades9,binomial)
anova(scmodelD8,scmodelD9,test="Chi")

# 22. North American Oak and 23. Far East Russian are sister taxa  and are NOT SIG DIFF (df=1,P=1) : Note: just one strain
clades10<-clades8
levels(clades10)[14:15]<-"Clades22_23"
levels(clades10)
scmodelD10<-glm(gain~clades10,binomial)
anova(scmodelD8,scmodelD10,test="Chi")

# 20. Ecuadorian and 21. CHNV are sister taxa  and are NOT SIG DIFF (df=1,P=0.4305)
clades11<-clades10
levels(clades11)[12:13]<-"Clades20_21"
levels(clades11)
scmodelD11<-glm(gain~clades11,binomial)
anova(scmodelD10,scmodelD11,test="Chi")

# clades 20 and 21 are sister taxa to 19 Malaysian and are NOT SIG DIFF (df=1,P=0.4915)
clades12<-clades11
levels(clades12)[11:12]<-"Clades19to21"
levels(clades12)
scmodelD12<-glm(gain~clades12,binomial)
anova(scmodelD11,scmodelD12,test="Chi")

# NOTE: that 18. Far East Asia has split in two, but this makes no difference to the simplification

# AIC: 337.26
anova(scmodelD12,test="Chi")
summary(scmodelD12)

# Combine all clades with the ancestral proportion of chr2-16 gains
# start with  5. French dairy and 6. African beer NOT SIG DIFF FROM WINE ETC (df=1,P=0.07227)
clades13<-clades12
levels(clades13)[1]<-"ancestral"
levels(clades13)[17]<-"ancestral"
scmodelD13<-glm(gain~clades13,binomial)
anova(scmodelD12,scmodelD13,test="Chi")
clades14<-clades13
levels(clades14)[1:2]<-"ancestral"
levels(clades14)[3:12]<-"ancestral"
levels(clades14)[4:5]<-"ancestral"
levels(clades14)
scmodelD14<-glm(gain~clades14,binomial)
# combining all the remaining ancestral clades (df=13,P=0.4983)
anova(scmodelD13,scmodelD14,test="Chi")
anova(scmodelD14,test="Chi")
summary(scmodelD14)
# AIC: 324.86, but only 7.60627% deviance
(345.11-318.86)/345.11

# not signfiicantly worse than all 26 clades (deviance = - 3.22514%, df=23, P=0.3373)
anova(model453Dv,scmodelD14,test="Chi")
# model predictions
# ancestral Clades8_11   25. Sake 
# 0.1072386  0.3437500  0.5294118 
table(droplevels(clades14),gain)
tapply(predict(scmodelD14,type="response"),clades14,mean)


