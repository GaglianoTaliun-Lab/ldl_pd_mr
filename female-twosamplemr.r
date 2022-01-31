library(TwoSampleMR)

#dat.csv will require minimally SNP (rsID), beta.exposure, beta.outcome, se.exposure, se.outcome, effect_allele.exposure, effect_allele.outcome, exposure, id.exposure (source name), outcome, id.outcome, and mr_keep (=TRUE)

#Read .csv file
dat<-read.csv("pathway/dat.csv")
#To verify file/print:
dat
#Example file: 
LDLPD<-read.csv("pathway/LDLPD-female.csv")
LDLPD
  X        SNP chrpos.exposure effect_allele.exposure other_allele.exposure beta.exposure se.exposure pval.exposure
1 1 rs10062361      5:74565153                      T                     C      0.066717   0.0039256      9.97e-65
2 2 rs11206510      1:55496039                      C                     T     -0.050126   0.0041043      2.73e-34
3 3 rs12714264      2:21265518                      T                     A     -0.096566   0.0048150      2.27e-89
4 4  rs2073547      7:44582331                      G                     A      0.040959   0.0041382      4.31e-23
5 5  rs2495477      1:55518467                      G                     A     -0.037018   0.0033441      1.80e-28
6 6  rs4341893      2:21135577                      G                     A     -0.040996   0.0033505      2.07e-34
  chrpos.outcome effect_allele.outcome other_allele.outcome eaf.outcome beta.outcome pval.outcome se.outcome        Wald
1     5:74565153                     T                    C      0.2222       0.0061      0.82800     0.0282  0.09143097
2     1:55496039                     T                    C      0.8088       0.0090      0.75630     0.0292  0.17954754
3     2:21265518                     A                    T      0.8601      -0.0658      0.05318     0.0340 -0.68139925
4     7:44582331                     A                    G      0.7227      -0.0467      0.13150     0.0310  1.14016455
5     1:55518467                     A                    G      0.6069       0.0196      0.46950     0.0270  0.52947215
6     2:21135577                     A                    G      0.3247       0.0277      0.27850     0.0256  0.67567568
    Waldvar        lab id.exposure exposure   id.outcome outcome mr_keep
1 0.1786591 rs10062361         UKB    LDL-C Blauwendraat      PD    TRUE
2 0.3393436 rs11206510         UKB    LDL-C Blauwendraat      PD    TRUE
3 0.1239679 rs12714264         UKB    LDL-C Blauwendraat      PD    TRUE
4 0.5728286  rs2073547         UKB    LDL-C Blauwendraat      PD    TRUE
5 0.5319877  rs2495477         UKB    LDL-C Blauwendraat      PD    TRUE
6 0.3899393  rs4341893         UKB    LDL-C Blauwendraat      PD    TRUE


#Make sure effect alleles match. Otherwise change sign of beta:
dat$effect_allele.outcome<-as.factor(dat$effect_allele.outcome)
dat$effect_allele.exposure<-as.factor(dat$effect_allele.exposure)
lev2 <- unique( c( levels(dat$effect_allele.outcome), levels(dat$effect_allele.exposure) ) )
dat$effect_allele.outcome <- factor(dat$effect_allele.outcome, levels=lev2)
dat$effect_allele.exposure <- factor(dat$effect_allele.exposure, levels=lev2)
dat$effect_allele.exposure<-gsub(" ", "",dat$effect_allele.exposure, fixed = TRUE)
dat$beta.exposure[dat$effect_allele.exposure!=dat$effect_allele.outcome]<-dat$beta.exposure[dat$effect_allele.exposure!=dat$effect_allele.outcome] * -1
#Result for example file:
  X        SNP chrpos.exposure effect_allele.exposure other_allele.exposure beta.exposure se.exposure pval.exposure
1 1 rs10062361      5:74565153                      T                     C      0.066717   0.0039256      9.97e-65
2 2 rs11206510      1:55496039                      C                     T      0.050126   0.0041043      2.73e-34
3 3 rs12714264      2:21265518                      T                     A      0.096566   0.0048150      2.27e-89
4 4  rs2073547      7:44582331                      G                     A     -0.040959   0.0041382      4.31e-23
5 5  rs2495477      1:55518467                      G                     A      0.037018   0.0033441      1.80e-28
6 6  rs4341893      2:21135577                      G                     A      0.040996   0.0033505      2.07e-34
  chrpos.outcome effect_allele.outcome other_allele.outcome eaf.outcome beta.outcome pval.outcome se.outcome        Wald
1     5:74565153                     T                    C      0.2222       0.0061      0.82800     0.0282  0.09143097
2     1:55496039                     T                    C      0.8088       0.0090      0.75630     0.0292  0.17954754
3     2:21265518                     A                    T      0.8601      -0.0658      0.05318     0.0340 -0.68139925
4     7:44582331                     A                    G      0.7227      -0.0467      0.13150     0.0310  1.14016455
5     1:55518467                     A                    G      0.6069       0.0196      0.46950     0.0270  0.52947215
6     2:21135577                     A                    G      0.3247       0.0277      0.27850     0.0256  0.67567568
    Waldvar        lab id.exposure exposure   id.outcome outcome mr_keep
1 0.1786591 rs10062361         UKB    LDL-C Blauwendraat      PD    TRUE
2 0.3393436 rs11206510         UKB    LDL-C Blauwendraat      PD    TRUE
3 0.1239679 rs12714264         UKB    LDL-C Blauwendraat      PD    TRUE
4 0.5728286  rs2073547         UKB    LDL-C Blauwendraat      PD    TRUE
5 0.5319877  rs2495477         UKB    LDL-C Blauwendraat      PD    TRUE
6 0.3899393  rs4341893         UKB    LDL-C Blauwendraat      PD    TRUE


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
#Example file (default mode):
res<-mr(LDLPD)
res
#  id.exposure   id.outcome outcome exposure                    method nsnp           b        se       pval
#1         UKB Blauwendraat      PD    LDL-C                  MR Egger    6 -1.50012998 0.6190943 0.07251925
#2         UKB Blauwendraat      PD    LDL-C           Weighted median    6  0.10252434 0.2936191 0.72695818
#3         UKB Blauwendraat      PD    LDL-C Inverse variance weighted    6  0.01152343 0.2624287 0.96497558
#4         UKB Blauwendraat      PD    LDL-C               Simple mode    6  0.36802608 0.5233763 0.51332918
#5         UKB Blauwendraat      PD    LDL-C             Weighted mode    6  0.17865945 0.4519793 0.70893539
#Note: TwoSampleMR uses bootstrapping to calculate SEs for the weighted median, simple mode and weighted mode methods. It is therefore normal for the standard errors and p-values of those methods
to differ between each run of the command.


#By individual analysis
#General format: take xmethod from mr_method_list
res<-xmethod(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)

#Examples
#Inverse variance weighted
mr_ivw(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
#Example file:
mr_ivw(LDLPD$beta.exposure, LDLPD$beta.outcome, LDLPD$se.exposure, LDLPD$se.outcome)
#$b
#[1] 0.01152343
#
#$se
#[1] 0.2624287
#
#$pval
#[1] 0.9649756
#
#$nsnp
#[1] 6
#
#$Q
#[1] 7.851282
#
#$Q_df
#[1] 5
#
#$Q_pval
#[1] 0.1646245

#MR-Egger
mr_egger_regression(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
#Example file:
mr_egger_regression(LDLPD$beta.exposure, LDLPD$beta.outcome, LDLPD$se.exposure, LDLPD$se.outcome)
#$b
#[1] -1.50013
#
#$se
#[1] 0.6190943
#
#$pval
#[1] 0.07251925
#
#$nsnp
#[1] 6
#
#$b_i
#[1] 0.09018586
#
#$se_i
#[1] 0.03475798
#
#$pval_i
#[1] 0.06038538
#
#$Q
#[1] 1.118915
#
#$Q_df
#[1] 4
#
#$Q_pval
#[1] 0.8912596
#
#$mod
#
#Call:
#lm(formula = b_out ~ b_exp, weights = 1/se_out^2)
#
#Weighted Residuals:
#       1        2        3        4        5        6 
# 0.56732 -0.20515 -0.32719  0.57929 -0.55756 -0.03854 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  0.09019    0.01838   4.906  0.00801 **
#b_exp       -1.50013    0.32744  -4.581  0.01017 * 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.5289 on 4 degrees of freedom
#Multiple R-squared:  0.8399,	Adjusted R-squared:  0.7999 
#F-statistic: 20.99 on 1 and 4 DF,  p-value: 0.01017


#$dat
#    b_out    b_exp    se_exp se_out flipped
#1  0.0061 0.066717 0.0039256 0.0282   FALSE
#2  0.0090 0.050126 0.0041043 0.0292   FALSE
#3 -0.0658 0.096566 0.0048150 0.0340   FALSE
#4  0.0467 0.040959 0.0041382 0.0310    TRUE
#5  0.0196 0.037018 0.0033441 0.0270   FALSE
#6  0.0277 0.040996 0.0033505 0.0256   FALSE

#Weighted median
mr_weighted_median(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)
#$b
#[1] 0.1025243
#
#$se
#[1] 0.2783118
#
#$pval
#[1] 0.7125904
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
#[1] 0.1786594
#
#$se
#[1] 0.4252419
#
#$pval
#[1] 0.69183
#
#$nsnp
#[1] 6

#Analysis by single SNPs 
#Default methods: single_method = "mr_wald_ratio", all_method = c("mr_ivw", "mr_egger_regression")
res_single<-mr_singlesnp(dat)
res_single
#Only all_method can be precised:
res_single<-mr_singlesnp(dat, all_method = c("xmethod"))
#Example file (with default methods):
res_single<-mr_singlesnp(LDLPD)
res_single
#exposure outcome id.exposure   id.outcome samplesize                             SNP           b        se          p
#1    LDL-C      PD         UKB Blauwendraat         NA                      rs10062361  0.09143097 0.4226809 0.82874451
#2    LDL-C      PD         UKB Blauwendraat         NA                      rs11206510  0.17954754 0.5825320 0.75791556
#3    LDL-C      PD         UKB Blauwendraat         NA                      rs12714264 -0.68139925 0.3520908 0.05295421
#4    LDL-C      PD         UKB Blauwendraat         NA                       rs2073547  1.14016455 0.7568544 0.13195128
#5    LDL-C      PD         UKB Blauwendraat         NA                       rs2495477  0.52947215 0.7293749 0.46788418
#6    LDL-C      PD         UKB Blauwendraat         NA                       rs4341893  0.67567568 0.6244512 0.27923864
#7    LDL-C      PD         UKB Blauwendraat         NA All - Inverse variance weighted  0.01152343 0.2624287 0.96497558
#8    LDL-C      PD         UKB Blauwendraat         NA                  All - MR Egger -1.50012998 0.6190943 0.07251925
#Example file with weighted mode as all_method:
res_single<-mr_singlesnp(LDLPD, all_method = c("mr_weighted_mode"))
res_single
#  exposure outcome id.exposure   id.outcome samplesize                 SNP           b        se          p
#1    LDL-C      PD         UKB Blauwendraat         NA          rs10062361  0.09143097 0.4226809 0.82874451
#2    LDL-C      PD         UKB Blauwendraat         NA          rs11206510  0.17954754 0.5825320 0.75791556
#3    LDL-C      PD         UKB Blauwendraat         NA          rs12714264 -0.68139925 0.3520908 0.05295421
#4    LDL-C      PD         UKB Blauwendraat         NA           rs2073547  1.14016455 0.7568544 0.13195128
#5    LDL-C      PD         UKB Blauwendraat         NA           rs2495477  0.52947215 0.7293749 0.46788418
#6    LDL-C      PD         UKB Blauwendraat         NA           rs4341893  0.67567568 0.6244512 0.27923864
#7    LDL-C      PD         UKB Blauwendraat         NA All - Weighted mode  0.17865945 0.4487005 0.70693714


#Leave-one-out Analysis 
#Default method: method = mr_ivw
res_loo<-mr_leaveoneout(dat)
res_loo
#To specify method:
res_loo<-mr_leaveoneout(dat, method=xmethod)
#Example file (default method):
res_loo<-mr_leaveoneout(LDLPD)
res_loo
#  exposure outcome id.exposure   id.outcome samplesize        SNP           b        se         p
#1    LDL-C      PD         UKB Blauwendraat         NA rs10062361 -0.01447501 0.3367586 0.9657148
#2    LDL-C      PD         UKB Blauwendraat         NA rs11206510 -0.01341608 0.3125068 0.9657569
#3    LDL-C      PD         UKB Blauwendraat         NA rs12714264  0.39088407 0.2605184 0.1335086
#4    LDL-C      PD         UKB Blauwendraat         NA  rs2073547 -0.08205520 0.2542248 0.7468722
#5    LDL-C      PD         UKB Blauwendraat         NA  rs2495477 -0.03501416 0.2953868 0.9056425
#6    LDL-C      PD         UKB Blauwendraat         NA  rs4341893 -0.07264361 0.2850432 0.7988380
#7    LDL-C      PD         UKB Blauwendraat         NA        All  0.01152343 0.2624287 0.9649756
#Example file (using MR-Egger):
res_loo<-mr_leaveoneout(LDLPD, method=mr_egger_regression)
#res_loo
#  exposure outcome id.exposure   id.outcome samplesize        SNP          b        se          p
#1    LDL-C      PD         UKB Blauwendraat         NA rs10062361 -1.6461225 0.6571428 0.08732345
#2    LDL-C      PD         UKB Blauwendraat         NA rs11206510 -1.5088658 0.6203153 0.09312859
#3    LDL-C      PD         UKB Blauwendraat         NA rs12714264 -0.8452433 1.1704858 0.52241318
#4    LDL-C      PD         UKB Blauwendraat         NA  rs2073547 -1.3938083 0.6404871 0.11777165
#5    LDL-C      PD         UKB Blauwendraat         NA  rs2495477 -1.6845228 0.6764068 0.08845289
#6    LDL-C      PD         UKB Blauwendraat         NA  rs4341893 -1.5097887 0.6539863 0.10416500
#7    LDL-C      PD         UKB Blauwendraat         NA        All -1.5001300 0.6190943 0.07251925


#Plots
#Scatter plot (lines with methods used in res)
scatter<-mr_scatter_plot(res, dat)

#Forest
forest<-mr_forest_plot(res_single)

#Funnel plot
funnel<-mr_funnel_plot(res_single)

#Leave-one-out
leaveoneout<-mr_leaveoneout_plot(res_loo)

#To save as pdf:
pdf("pathway/name.pdf")
mr_forest_plot(res_single) #or other plot function
dev.off()


#Other tests:
#For heterogeneity (Q) by Inverse variance weighted method and MR-Egger method
mr_heterogeneity(dat)
#Example file:
mr_heterogeneity(LDLPD)
  id.exposure   id.outcome outcome exposure                    method        Q Q_df    Q_pval
1         UKB Blauwendraat      PD    LDL-C                  MR Egger 1.118915    4 0.8912596
2         UKB Blauwendraat      PD    LDL-C Inverse variance weighted 7.851282    5 0.1646245


#Testing for horizontal pleiotropy using MR-Egger
mr_pleiotropy_test(dat)
#Example file:
mr_pleiotropy_test(LDLPD)
 id.exposure   id.outcome outcome exposure egger_intercept         se       pval
1         UKB Blauwendraat      PD    LDL-C      0.09018586 0.03475798 0.06038538


#Testing for horizontal pleiotropy using MR-PResso
library(MRPRESSO)
mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure", SdOutcome="se.outcome", SdExposure="se.exposure", OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data = dat)
#Example file:
mr_presso(BetaOutcome="beta.outcome", BetaExposure="beta.exposure", SdOutcome="se.outcome", SdExposure="se.exposure", OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data = LDLPD)
$`Main MR results`
       Exposure       MR Analysis Causal Estimate        Sd     T-stat   P-value
1 beta.exposure               Raw      0.01152343 0.2624287 0.04391071 0.9666753
2 beta.exposure Outlier-corrected              NA        NA         NA        NA

$`MR-PRESSO results`
$`MR-PRESSO results`$`Global Test`
$`MR-PRESSO results`$`Global Test`$RSSobs
[1] 14.09026

$`MR-PRESSO results`$`Global Test`$Pvalue
[1] 0.175
#Note: the Outlier Test was not done because the Global Test p-value was >0.05 (not significant).

#directionality test
outcome_dat<-read_outcome_data("LDLPD-female.csv", beta_col="beta.outcome", se_col="se.outcome", effect_allele_col="effect_allele.outcome", other_allele_col="other_allele.outcome", eaf_col="eaf.outcome", pval_col="pval.outcome", sep=",")
exposure_dat<-read_exposure_data("LDLPD-female.csv", beta_col="beta.exposure", se_col="se.exposure", effect_allele_col="effect_allele.exposure", other_allele_col="other_allele.exposure", eaf_col="eaf.outcome", pval_col="pval.exposure", sep=",")
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
#id.exposure id.outcome exposure outcome snp_r2.exposure snp_r2.outcome
#1      ZL8A9s     8VyBqt exposure outcome      0.00654467    0.000373162
#  correct_causal_direction steiger_pval
#1                     TRUE           NA


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

