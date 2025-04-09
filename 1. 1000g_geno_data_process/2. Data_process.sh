wd="/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all"
module load PLINK/2

# Step 1. Transform pfile to bfile contatinting ~20M 1KG SNP
plink2 --pfile ${wd}/all_hg38 'vzs' \
--allow-extra-chr \
--chr 1-22 \
--max-alleles 2 \
--extract ${wd}/bfile/snplist_1kg_20M.txt \
--make-bed \
--out ${wd}/bfile/hg38/all_variants/all_hg38_1kg_20M

# Step 2. Keep only "ACTG" SNPs.
plink2 --bfile ${wd}/bfile/hg38/all_variants/all_hg38_1kg_20M \
--snps-only \
--make-bed \
--out ${wd}/bfile/hg38/all_variants/all_hg38_1kg_20M_snponly

# Step 3. Remove "A-T","C-G","T-A","G-C" pairs
module load PLINK/1
qcdir="/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg38"
refdir="/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg38/all_variants"
refname="all_hg38_1kg_20M_snponly"

awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" \
|| $5$6 == "AT" || $5$6 == "TA") {print $2}' \
$refdir/$refname.bim > \
$qcdir/$refname.ac_gt_snps

plink --bfile $refdir/$refname \
--exclude $qcdir/$refname.ac_gt_snps \
--make-bed \
--out $refdir/$refname.no_ac_gt_snps

# Step 4. Use liftOver to convert hg38 to hg19 position
# Generate input for liftover
awk '{print "chr" $1, $4 -1, $4, $2 }' $refdir/$refname.no_ac_gt_snps.bim > \
$refdir/1kg_phase3_hg38.tolift

# Mapping
module load liftOver
liftOver $refdir/1kg_phase3_hg38.tolift ${wd}/bfile/hg38/hg38ToHg19.over.chain \
${wd}/bfile/hg19/all_variants/liftover_output_hg19 ${wd}/bfile/hg19/all_variants/unmapped

# Step 5. Update the original bfile with new variants and position
cd ${wd}/bfile/hg19/all_variants

# Ectract mapped variants 
awk '{print $4}' liftover_output_hg19 > 1kg_CGRCh37.snps

# Ectract updated positions 
awk '{print $4, $3}' liftover_output_hg19 > 1kg_CGRCh37.pos
sort 1kg_CGRCh37.pos | uniq > 1kg_CGRCh37.pos.unique
# Remove duplicated variants with same rsid but different BP
awk '!seen[$1]++' 1kg_CGRCh37.pos.unique > 1kg_CGRCh37.pos.unique.bp

# Update bfile
module load PLINK/2
plink2 --bfile ${refdir}/${refname}.no_ac_gt_snps \
--rm-dup exclude-all 'list' \
--exclude unmapped \
--make-bed \
--out ${wd}/bfile/hg19/all_variants/rm_dup_all/all_hg19_1kg_20M_unupdated

## Alternative choice: use the --rm-dup force-first option, then update all position.
## Step 5a. Use the --rm-dup force-first option
# plink2 --bfile ${refdir}/${refname}.no_ac_gt_snps \
# --rm-dup force-first \
# --exclude unmapped \
# --make-bed \
# --out ${wd}/bfile/hg19/all_variants/keep_first/all_hg19_1kg_20M_unupdated

plink2 --bfile ${wd}/bfile/hg19/all_variants/keep_first/all_hg19_1kg_20M_unupdated \
--update-map 1kg_CGRCh37.pos.unique.bp \
--make-bed \
--out ${wd}/bfile/hg19/all_variants/keep_first/all_hg19_1kg_20M

plink2 --bfile ${wd}/bfile/hg19/all_variants/rm_dup_all/all_hg19_1kg_20M_unupdated \
--update-map 1kg_CGRCh37.pos.unique.bp \
--make-bed \
--out ${wd}/bfile/hg19/all_variants/rm_dup_all/all_hg19_1kg_20M

# Step 6. Sorted variants by CHR and BP
plink2 --bfile ${wd}/bfile/hg19/all_variants/rm_dup_all/all_hg19_1kg_20M \
--sort-vars \
--make-pgen \
--out ${wd}/bfile/hg19/all_variants/rm_dup_all/all_hg19_1kg_20M_sorted \

# Step 7. Extract dataset for each population
cd ${wd}/bfile/hg19
for pop in EUR EAS SAS AFR AMR; do
    sample_list="/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/${pop}.txt"
plink2 --pfile ${wd}/bfile/hg19/all_variants/rm_dup_all/all_hg19_1kg_20M_sorted \
--keep-allele-order \
--keep $sample_list \
--make-bed \
--out ${pop}/${pop}_1kg
done

# Step 8. 0.01 maf filter in each population 
for pop in EUR EAS AFR SAS AMR; do
plink2 --bfile ${wd}/bfile/hg19/${pop}/${pop}_1kg \
--keep-allele-order \
--maf 0.01 \
--make-bed \
--out ${wd}/bfile/hg19/${pop}/${pop}_1kg_maf_0.01
done