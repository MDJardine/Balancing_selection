# simple command to extract data from any number of IHA sheets
# data sheets must have ending _raw.xls

for filename in *_raw.xls; do echo $filename; Rscript excel2other.R $filename output/out_$filename; done


## these files can then be imported and tinkered with in R
## long ter would like to merge these data sets automatically in R but beyond skills atm
## doesn't matter as onyl haev 1 experiment but if usign more in the future it is somehting to think about
