#R version 4.0.5
#TwoSampleMR version 0.5.6
#Mendelian Randomization version 0.5.1

library(TwoSampleMR)

#dat.csv will require minimally SNP (rsID), beta.exposure, beta.outcome, se.exposure, se.outcome, effect_allele.exposure, effect_allele.outcome, exposure, id.exposure (source name), outcome, id.outcome, and mr_keep (=TRUE)

#Read .csv file
dat<-read.csv("Moriginalbeta.csv")
#To verify file/print:
dat

#Make sure effect alleles match. Otherwise change sign of beta:
dat$effect_allele.outcome<-as.factor(dat$effect_allele.outcome)
dat$effect_allele.exposure<-as.factor(dat$effect_allele.exposure)
lev2 <- unique( c( levels(dat$effect_allele.outcome), levels(dat$effect_allele.exposure) ) )
dat$effect_allele.outcome <- factor(dat$effect_allele.outcome, levels=lev2)
dat$effect_allele.exposure <- factor(dat$effect_allele.exposure, levels=lev2)
dat$effect_allele.exposure<-gsub(" ", "",dat$effect_allele.exposure, fixed = TRUE)
dat$beta.exposure[dat$effect_allele.exposure!=dat$effect_allele.outcome]<-dat$beta.exposure[dat$effect_allele.exposure!=dat$effect_allele.outcome] * -1

#To see MR methods that can be used:
mr_method_list()
#Default methods: IVW, MR-Egger, Weighted Median, Simple Mode, Weighted Mode

#To view default parameters (parameters used unless otherwise specified):
default_parameters()


#Summary MR results
res<-mr(dat)
res
#To precise method, take xmethod, ymethod, etc. from method list:
res<-mr(dat, method_list=c("xmethod", "ymethod"))
#Ex: 
#res<-mr(dat, method_list=c("mr_uwr", "mr_ivw_radial")
#res

#Note: TwoSampleMR uses bootstrapping to calculate SEs for the weighted median, simple mode and weighted mode methods. It is therefore normal for the standard errors and p-values of those methods
#to differ between each run of the command.


#By individual analysis
#General format: take xmethod from mr_method_list
#res<-xmethod(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)

#Examples
#Inverse variance weighted
mr_ivw(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)

#MR-Egger
mr_egger_regression(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)


#Weighted median
mr_weighted_median(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)


#Weighted mode
mr_weighted_mode(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)

#Analysis by single SNPs 
#Default methods: single_method = "mr_wald_ratio", all_method = c("mr_ivw", "mr_egger_regression")
res_single<-mr_singlesnp(dat)
res_single
#Only all_method can be precised:
#res_single<-mr_singlesnp(dat, all_method = c("xmethod"))


#res_single<-mr_singlesnp(dat, all_method = c("mr_weighted_mode"))
#res_single


#Leave-one-out Analysis 
#Default method: method = mr_ivw
res_loo<-mr_leaveoneout(dat)
res_loo
#To specify method:
#res_loo<-mr_leaveoneout(dat, method=xmethod)



#Plots
#Scatter plot (lines with methods used in res)
pdf("m-standard-scatter.pdf")
mr_scatter_plot(res, dat)
dev.off()

#Forest
pdf("m-standard-forest.pdf")
mr_forest_plot(res_single) #or other plot function
dev.off()

#Funnel plot
pdf("m-standard-funnel.pdf")
mr_funnel_plot(res_single)
dev.off()

#Leave-one-out
pdf("m-standard-loo.pdf")
mr_leaveoneout_plot(res_loo)
dev.off()


#Other tests:
#For heterogeneity (Q) by Inverse variance weighted method and MR-Egger method


#Testing for horizontal pleiotropy using MR-Egger
mr_pleiotropy_test(dat)
#id.exposure   id.outcome outcome exposure egger_intercept         se
#1         UKB Blauwendraat      PD    LDL-C     0.006480583 0.01018367
#       pval
#1 0.5325446


#Testing for horizontal pleiotropy using MR-PResso
library(MRPRESSO)
mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure", SdOutcome="se.outcome", SdExposure="se.exposure", OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data = dat)
#Note: the Outlier Test was not done because the Global Test p-value was >0.05 (not significant).
#$`MR-PRESSO results`$`Global Test`$Pvalue
#[1] 0.369

#Sources:
#TwoSampleMR: https://rdrr.io/github/MRCIEU/TwoSampleMR/api/
#User-friendly guide to TwoSampleMR: https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html 
#MR-Base GWAS: https://gwas.mrcieu.ac.uk/
#MR-Presso: https://rdrr.io/github/rondolab/MR-PRESSO/man/
        
#I2GX
library('MendelianRandomization')
MRInputObject <- mr_input(dat$beta.exposure, dat$se.exposure, dat$beta.outcome, dat$se.outcome)
mr_egger(MRInputObject)
#MR-Egger method
#(variants uncorrelated, random-effect model)

#Number of Variants =  20 

#------------------------------------------------------------------
#      Method Estimate Std Error  95% CI       p-value
#    MR-Egger   -0.159     0.164 -0.480, 0.163   0.334
# (intercept)    0.006     0.010 -0.013, 0.026   0.525
#------------------------------------------------------------------
#Residual Standard Error :  1.101 
#Heterogeneity test statistic = 21.8264 on 18 degrees of freedom, (p-value = 0.2398)
#I^2_GX statistic: 98.6%
