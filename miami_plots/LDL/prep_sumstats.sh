##LDL sumstats downloaded from Ben Neale UK Biobank Round 2 https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=178908679

#let's copy the first column 'variant' to create 4 additional columns CHR BP REF ALT
#first unzip
zless -S 30780_irnt.gwas.imputed_v3.female.varorder.tsv.bgz > 30780_irnt.gwas.imputed_v3.female.varorder.tsv
#copy first column into a new file and split into 4 new columns 
cut -f 1 30780_irnt.gwas.imputed_v3.female.varorder.tsv > female-col1
perl -p -i -e 's/:/\t/g' female-col1

#now just need to change the header from variant to four new headers
sed -i 's/variant/CHR:BP:REF:ALT/' female-col1
perl -p -i -e 's/:/\t/g' female-col1

#merge 4 new columns to original sumstats file
#remove NaN columns
#remove chrX
paste female-col1 30780_irnt.gwas.imputed_v3.female.varorder.tsv | grep -v NaN | awk '$1!=X'> 30780_irnt.gwas.imputed_v3.female.varorder.txt

#clean up
rm female-col1 30780_irnt.gwas.imputed_v3.female.varorder.tsv  

#repeat for male-only sumstats file
