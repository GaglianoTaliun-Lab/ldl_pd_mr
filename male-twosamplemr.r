library(TwoSampleMR)

dat<-read.csv("LDLPD-male.csv")

LDLPD<-dat

#Make sure effect alleles match. Otherwise change sign of beta for exposure
dat$effect_allele.outcome<-as.factor(dat$effect_allele.outcome)
dat$effect_allele.exposure<-as.factor(dat$effect_allele.exposure)
lev2 <- unique( c( levels(dat$effect_allele.outcome), levels(dat$effect_allele.exposure) ) )
dat$effect_allele.outcome <- factor(dat$effect_allele.outcome, levels=lev2)
dat$effect_allele.exposure <- factor(dat$effect_allele.exposure, levels=lev2)
dat$effect_allele.exposure<-gsub(" ", "",dat$effect_allele.exposure, fixed = TRUE)
dat$beta.exposure[dat$effect_allele.exposure!=dat$effect_allele.outcome]<-dat$beta.exposure[dat$effect_allele.exposure!=dat$effect_allele.outcome] * -1

res<-mr(dat)
res

#id.exposure   id.outcome outcome exposure                    method nsnp
#1         UKB Blauwendraat      PD    LDL-C                  MR Egger    6
#2         UKB Blauwendraat      PD    LDL-C           Weighted median    6
#3         UKB Blauwendraat      PD    LDL-C Inverse variance weighted    6
#4         UKB Blauwendraat      PD    LDL-C               Simple mode    6
#5         UKB Blauwendraat      PD    LDL-C             Weighted mode    6
#            b        se      pval
#1 -0.29739773 0.6651009 0.6779132
#2 -0.24434399 0.2691345 0.3639373
#3 -0.07285402 0.2628024 0.7816113
#4 -0.28293531 0.5513231 0.6296790
#5 -0.29167642 0.3255234 0.4112993

#Inverse variance weighted
mr_ivw(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
#$b
#[1] -0.07285402
#
#$se
#[1] 0.2628024
#
#$pval
#[1] 0.7816113
#
#$nsnp
#[1] 6
#
#$Q
#[1] 8.377937
#
#$Q_df
#[1] 5
#
#$Q_pval
#[1] 0.1366002

#Weighted median
mr_weighted_median(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
#$b
#[1] -0.244344
#
#$se
#[1] 0.2646133
#
#$pval
#[1] 0.3557987
#
#$Q
#[1] NA
#
#$Q_df
#[1] NA
#
#$Q_pval
#[1] NA
#
#$nsnp
#[1] 6

 
#Weighted mode
mr_weighted_mode(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
#$b
#[1] -0.2916764
#
#$se
#[1] 0.3152381
#
#$pval
#[1] 0.3972926
#
#$nsnp
#[1] 6

#Analysis by single SNPs 
#Default methods: single_method = "mr_wald_ratio", all_method = c("mr_ivw", "mr_egger_regression")
res_single<-mr_singlesnp(dat)
res_single

#exposure outcome id.exposure   id.outcome samplesize
#1    LDL-C      PD         UKB Blauwendraat         NA
#2    LDL-C      PD         UKB Blauwendraat         NA
#3    LDL-C      PD         UKB Blauwendraat         NA
#4    LDL-C      PD         UKB Blauwendraat         NA
#5    LDL-C      PD         UKB Blauwendraat         NA
#6    LDL-C      PD         UKB Blauwendraat         NA
#7    LDL-C      PD         UKB Blauwendraat         NA
#8    LDL-C      PD         UKB Blauwendraat         NA
#                              SNP           b        se          p
#1                      rs10062361  0.86464646 0.5030303 0.08563704
#2                      rs11206510 -1.00220264 0.5660793 0.07665587
#3                      rs12714264 -0.29514563 0.2951456 0.31731051
#4                       rs2073547 -0.36394558 0.9217687 0.69296544
#5                       rs2495477  0.63084112 0.5677570 0.26652053
#6                       rs4341893 -0.13818182 0.8363636 0.86877288
#7 All - Inverse variance weighted -0.07285402 0.2628024 0.78161135
#8                  All - MR Egger -0.29739773 0.6651009 0.67791317

#Leave-one-out Analysis 
#Default method: method = mr_ivw
res_loo<-mr_leaveoneout(dat)
res_loo

#  exposure outcome id.exposure   id.outcome samplesize        SNP           b
#1    LDL-C      PD         UKB Blauwendraat         NA rs10062361 -0.25528306
#2    LDL-C      PD         UKB Blauwendraat         NA rs11206510  0.06433277
#3    LDL-C      PD         UKB Blauwendraat         NA rs12714264  0.12679762
#4    LDL-C      PD         UKB Blauwendraat         NA  rs2073547 -0.05801266
#5    LDL-C      PD         UKB Blauwendraat         NA  rs2495477 -0.17602786
#6    LDL-C      PD         UKB Blauwendraat         NA  rs4341893 -0.06876354
#7    LDL-C      PD         UKB Blauwendraat         NA        All -0.07285402
#         se         p
#1 0.2281529 0.2631776
#2 0.2499933 0.7969178
#3 0.3779016 0.7372247
#4 0.2993291 0.8463255
#5 0.2796017 0.5289782
#6 0.3027639 0.8203310
#7 0.2628024 0.7816113

#Plots
#Scatter plot (lines with methods used in res)
pdf("scatter-male.pdf")
mr_scatter_plot(res, dat)
dev.off()

#Forest
pdf("forest-male.pdf")
mr_forest_plot(res_single)
dev.off()

#Funnel plot
pdf("funnel-male.pdf")
mr_funnel_plot(res_single)
dev.off()

#Leave-one-out
pdf("leaveoneout-male.pdf")
mr_leaveoneout_plot(res_loo)
dev.off()

#Testing for horizontal pleiotropy using MR-Egger
mr_pleiotropy_test(dat)
#  id.exposure   id.outcome outcome exposure egger_intercept         se
#1         UKB Blauwendraat      PD    LDL-C      0.01283231 0.03423931
#       pval
#1 0.7268459

#Testing for horizontal pleiotropy using MR-PResso
library(MRPRESSO)
mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure", SdOutcome="se.outcome", SdExposure="se.exposure", OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data = dat)

#$`Main MR results`
#       Exposure       MR Analysis Causal Estimate        Sd     T-stat
#1 beta.exposure               Raw     -0.07285402 0.2628024 -0.2772198
#2 beta.exposure Outlier-corrected              NA        NA         NA
#    P-value
#1 0.7927074
#2        NA
#
#$`MR-PRESSO results`
#$`MR-PRESSO results`$`Global Test`
#$`MR-PRESSO results`$`Global Test`$RSSobs
#[1] 12.68694
#
#$`MR-PRESSO results`$`Global Test`$Pvalue
#[1] 0.221

#directionality test
outcome_dat<-read_outcome_data("LDLPD-male.csv", beta_col="beta.outcome", se_col="se.outcome", effect_allele_col="effect_allele.outcome", other_allele_col="other_allele.outcome", eaf_col="eaf.outcome", pval_col="pval.outcome", sep=",")
exposure_dat<-read_exposure_data("LDLPD-male.csv", beta_col="beta.exposure", se_col="se.exposure", effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure", eaf_col="eaf.outcome", pval_col="pval.exposure", sep=",")
harmon<-harmonise_data(exposure_dat, outcome_dat, action = 2)
