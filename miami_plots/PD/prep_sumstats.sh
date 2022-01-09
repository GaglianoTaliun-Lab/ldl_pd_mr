##PD sumstats downloaded from IPDGC downloads page under 'Male and Female specific summary stats Blauwendraat et al 2021'

#need to replace first column MarkerName (format chr:pos) into two separate columns called CHR and BP, respectively
#first unzip
gzip -d FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz
#split chr:bp in MarkerName column into two tab-separated columns
perl -p -i -e 's/:/\t/g' FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt
#now just need to change the header from MarkerName to two new headers: CHR and BP
sed -i 's/MarkerName/CHR:BP/' FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt
perl -p -i -e 's/:/\t/g' FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt

#EasyStrata doesn't like dashes in column names, so let's rename the P-value header
sed -i 's/P-value/P/' FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt

#repeat for male-only sumstats file
