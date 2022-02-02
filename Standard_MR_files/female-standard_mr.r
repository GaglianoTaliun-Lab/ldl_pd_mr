#TwoSampleMR package installation on r/4.0.5
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)

#dat.csv will require minimally SNP (rsID), beta.exposure, beta.outcome, se.exposure, se.outcome, effect_allele.exposure, effect_allele.outcome, exposure, id.exposure (source name), outcome, id.outcome, and mr_keep (=TRUE)

#Read .csv file
dat<-read.csv("Foriginalbeta.csv")
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
res<-mr(dat, method_list=c("mr_uwr", "mr_ivw_radial")
res

#Note: TwoSampleMR uses bootstrapping to calculate SEs for the weighted median, simple mode and weighted mode methods. It is therefore normal for the standard errors and p-values of those methods
to differ between each run of the command.


#By individual analysis
#General format: take xmethod from mr_method_list
res<-xmethod(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)

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
pdf("f-standard-scatter.pdf")
scatter<-mr_scatter_plot(res, dat)
dev.off()

#Forest
pdf("f-standard-forest.pdf")
mr_forest_plot(res_single) #or other plot function
dev.off()

#Funnel plot
pdf("f-standard-funnel.pdf")
funnel<-mr_funnel_plot(res_single)
dev.off()

#Leave-one-out
pdf("f-standard-loo.pdf")
leaveoneout<-mr_leaveoneout_plot(res_loo)
dev.off()


#Other tests:
#For heterogeneity (Q) by Inverse variance weighted method and MR-Egger method


#Testing for horizontal pleiotropy using MR-Egger
mr_pleiotropy_test(dat)
#id.exposure   id.outcome outcome exposure egger_intercept          se
#1         UKB Blauwendraat      PD    LDL-C    0.0007451164 0.009181155
#       pval
#1 0.9359385


#Testing for horizontal pleiotropy using MR-PResso
library(MRPRESSO)
mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure", SdOutcome="se.outcome", SdExposure="se.exposure", OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data = dat)
#Note: the Outlier Test was not done because the Global Test p-value was >0.05 (not significant).
#$`MR-PRESSO results`$`Global Test`$Pvalue
#[1] 0.225
        
#directionality test
outcome_dat<-read_outcome_data("Foriginalbeta.csv", beta_col="beta.outcome", se_col="se.outcome", effect_allele_col="effect_allele.outcome", other_allele_col="other_allele.outcome", eaf_col="eaf.outcome", pval_col="pval.outcome", sep=",")
exposure_dat<-read_exposure_data("Foriginalbeta.csv", beta_col="beta.exposure", se_col="se.exposure", effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure", eaf_col="eaf.outcome", pval_col="pval.exposure", sep=",")
harmon<-harmonise_data(exposure_dat, outcome_dat, action = 2)
harmon$samplesize.exposure<-184689
harmon$ncase.outcome<-7384
harmon$ncontrol.outcome<-20330
harmon$prevalence.outcome<-0.01 #1% prevalence

harmon$r.outcome<-get_r_from_lor(
  harmon$beta.outcome,
  harmon$eaf.outcome,
  harmon$ncase.outcome,
  harmon$ncontrol.outcome,
  harmon$prevalence.outcome,
  model = "logit",
  correction = FALSE
)
directionality_test(harmon)


#Using GWAS datasets from MR-Base (https://gwas.mrcieu.ac.uk/)
library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- extract_instruments(c('ieu-id1'))
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ieu-id2'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)

#Sources:
#TwoSampleMR: https://rdrr.io/github/MRCIEU/TwoSampleMR/api/
#User-friendly guide to TwoSampleMR: https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html 
#MR-Base GWAS: https://gwas.mrcieu.ac.uk/
#MR-Presso: https://rdrr.io/github/rondolab/MR-PRESSO/man/
