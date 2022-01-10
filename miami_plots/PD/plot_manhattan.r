source('manhattan.plot-blue.R') #see https://github.com/GaglianoTaliun-Lab/useful_code_genetic_analysis/blob/main/Plotting/manhattan.plot-blue.R
dd<-read.table("FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt", h=T, as.is=T)
png("FemaleManhattan.png", width=950, height=500)
print(manhattan.plot(dd$CHR, dd$BP, dd$P, sig.level=5e-8))
dev.off()

#similar for MALE
