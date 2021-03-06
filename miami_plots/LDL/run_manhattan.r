source('manhattan.plot-blue.R') #see https://github.com/GaglianoTaliun-Lab/useful_code_genetic_analysis/blob/main/Plotting/manhattan.plot-blue.R

#female
prep<-read.table("30780_irnt.gwas.imputed_v3.female.varorder.txt", h=T, as.is=T) #the summary statistics, three required columns: chromsosome, position and p-value
dd<-subset(prep, prep$pval>0)
png("LDL-FemaleManhattan.png", width=950, height=500)
print(manhattan.plot(dd$CHR, dd$BP, dd$pval, sig.level=5e-8)) #where chr, pos and pvalue are the appropriate column headers in dd
dev.off()

#similar for MALE
prepmale<-read.table("30780_irnt.gwas.imputed_v3.male.varorder.txt", h=T, as.is=T) #the summary statistics, three required columns: chromsosome, position and p-value
dat<-subset(prepmale, prepmale$pval>0)
png("LDL-MaleManhattan.png", width=950, height=500)
print(manhattan.plot(dat$CHR, dat$BP, dat$pval, sig.level=5e-8)) #where chr, pos and pvalue are the appropriate column headers in dd
dev.off()

