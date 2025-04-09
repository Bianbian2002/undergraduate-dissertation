plink2 --bfile ${wd}/bfile/hg19/${pop}/${pop}_1kg_maf_0.01 \
--keep-allele-order \
--extract /gpfs/gibbs/pi/zhao/zb222/Data/snplist/Imputed/ukbAFR_imputed.txt \
--make-bed \
--out ${wd}/bfile/hg19/${pop}/${pop}_${snplist}_maf_0.01