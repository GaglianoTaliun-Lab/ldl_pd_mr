library(MendelianRandomization)
library(metafor)
par(mfrow=c(2,2)) 
par(mar=c(4,2,0,2)+.1)

#dat.csv will require minimally beta.exposure, beta.outcome, se.exposure, se.outcome
#Read .csv file
dat<-read.csv("pathway/dat.csv")

#If file has effect alleles (effect_allele.outcome and effect_alle.exposure), make sure effect alleles match. Otherwise change sign of beta (from Prof. Sarah's MR code):
dat$effect_allele.outcome<-as.factor(dat$effect_allele.outcome)
dat$effect_allele.exposure<-as.factor(dat$effect_allele.exposure)
lev2 <- unique( c( levels(dat$effect_allele.outcome), levels(dat$effect_allele.exposure) ) )
dat$effect_allele.outcome <- factor(dat$effect_allele.outcome, levels=lev2)
dat$effect_allele.exposure <- factor(dat$effect_allele.exposure, levels=lev2)
dat$effect_allele.exposure<-gsub(" ", "",dat$effect_allele.exposure, fixed = TRUE)
dat$beta.exposure[dat$effect_allele.exposure!=dat$effect_allele.outcome]<-dat$beta.exposure[dat$effect_allele.exposure!=dat$effect_allele.outcome] * -1

#Individual analysis. Various variables including models can be modified. Consult guidebook/package for more details

#First set MRInputObject:
MRInputObject <- mr_input(dat$beta.exposure, dat$se.exposure, dat$beta.outcome, dat$se.outcome)

#IVW (default model: fixed if 3 snps or less, random if otherwise)
mr_ivw(MRInputObject)

###results for LDLPD-female.csv:
#Inverse-variance weighted method
#(variants uncorrelated, random-effect model)

#Number of Variants : 6 

#------------------------------------------------------------------
# Method Estimate Std Error  95% CI       p-value
#    IVW    0.012     0.262 -0.503, 0.526   0.965
#------------------------------------------------------------------
#Residual standard error =  1.253 
#Heterogeneity test statistic (Cochran's Q) = 7.8513 on 5 degrees of freedom, (p-value = 0.1646). I^2 = 36.3%. 

#Changing model; e.g. fixed-effects model:
mr_ivw(MRInputObject, model="fixed")

###results for LDLPD-female.csv:
#Inverse-variance weighted method
#(variants uncorrelated, fixed-effect model)

#Number of Variants : 6 

#------------------------------------------------------------------
# Method Estimate Std Error  95% CI       p-value
#    IVW    0.012     0.209 -0.399, 0.422   0.956
#------------------------------------------------------------------
#Residual standard error =  1.253 
#Residual standard error is set to 1 in calculation of confidence interval by fixed-effect assumption.
#Heterogeneity test statistic (Cochran's Q) = 7.8513 on 5 degrees of freedom, (p-value = 0.1646). I^2 = 36.3%. 


#Egger (model always uses random effects)
mr_egger(MRInputObject)

###results for LDLPD-female.csv:
#MR-Egger method
#(variants uncorrelated, random-effect model)

#Number of Variants =  6 

#------------------------------------------------------------------
#      Method Estimate Std Error  95% CI        p-value
#    MR-Egger   -1.500     0.619 -2.714, -0.287   0.015
# (intercept)    0.090     0.035  0.022,  0.158   0.009
#------------------------------------------------------------------
#Residual Standard Error :  0.529 
#Residual standard error is set to 1 in calculation of confidence interval when its estimate is less than 1.
#Heterogeneity test statistic = 1.1189 on 4 degrees of freedom, (p-value = 0.8913)
#I^2_GX statistic: 94.9%

###results for LDLPD-male.csv:
#MR-Egger method
#(variants uncorrelated, random-effect model)
#
#Number of Variants =  6 
#
#------------------------------------------------------------------
#      Method Estimate Std Error  95% CI       p-value
#    MR-Egger   -0.297     0.665 -1.601, 1.006   0.655
# (intercept)    0.013     0.034 -0.054, 0.080   0.708
#------------------------------------------------------------------
#Residual Standard Error :  1.422 
#Heterogeneity test statistic = 8.0937 on 4 degrees of freedom, (p-value = 0.0882)
#I^2_GX statistic: 95.8%



#Median (default weighted; can also choose "penalized" or "simple")
mr_median(MRInputObject)
#To change whether weighted or not, ex:
mr_median(MRInputObject, weighting = "simple")
#Example file:
mr_median(MRInputObject)
Weighted median method 

#Number of Variants : 6 
#------------------------------------------------------------------
#                 Method Estimate Std Error  95% CI       p-value
# Weighted median method    0.102     0.294 -0.474, 0.679   0.727


#Leave one out analysis (uses ivw)
mr_loo(MRInputObject)


#Plots

#Forest (default method is ivw)
mr_forest(MRInputObject)
#To change method, ex:
mr_forest(MRInputObject, methods="egger")

#Funnel (vertical dashed line at ivw estiamte)
mr_funnel(MRInputObject)

#Scatter plot. Options for lines are ivw or mr_egger. 
#By default: ivw shown and directed to web interactive version:
mr_plot(MRInputObject)

#To change to egger, use line = "egger"
mr_plot(MRInputObject, line= "egger")
#To not go to web interactive version:
mr_plot(MRInputObject, interactive=FALSE)
#To combine:
mr_plot(MRInputObject, line="egger", interactive=FALSE)

#To generate pdf
pdf("pathway/name.pdf")
mr_funnel(MRInputObject) #insert plot command here, this is an example
dev.off()

#Sources:
# https://rdrr.io/cran/MendelianRandomization/api/
# Guidebook: https://cran.r-project.org/web/packages/MendelianRandomization/MendelianRandomization.pdf
