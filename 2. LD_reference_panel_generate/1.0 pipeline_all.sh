wd="/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all"

# For hapmap3 panel
Reference="mine_1KG_hapmap3"
snplist="hapmap3"

# For hapmap3+mega panel
Reference="mine_1KG_hapmap3_mega_union"
snplist="hapmap3_mega_union"

# For 20M 1kg panel
Reference="mine_1KG_20M"
snplist="1kg_20M"

# For Annotation+hm3 panel
Reference="BaselineLD_v2.3_hm3"
snplist="hm3+LDv2.3_union"

# Step 0. Generate dir
mkdir /gpfs/gibbs/pi/zhao/zb222/reference/${Reference}
cd /gpfs/gibbs/pi/zhao/zb222/reference/${Reference}
for pop in eur eas afr sas amr; do
mkdir ldblk_1kg_${pop}
done

mkdir ldblk
cd /gpfs/gibbs/pi/zhao/zb222/reference/${Reference}/ldblk
for pop in EUR EAS AFR SAS AMR; do
mkdir ${pop}
done

cd ..
mkdir snplist_ldblk
cd /gpfs/gibbs/pi/zhao/zb222/reference/${Reference}/snplist_ldblk
for pop in EUR EAS AFR SAS AMR; do
mkdir ${pop}
done

# Step 1. Extract corresponding SNP list
module load PLINK/2
for pop in EUR EAS AFR SAS AMR; do
plink2 --bfile ${wd}/bfile/hg19/${pop}/${pop}_1kg_maf_0.01 \
--keep-allele-order \
--extract /gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/snplist_${snplist}.txt \
--make-bed \
--out ${wd}/bfile/hg19/${pop}/${pop}_${snplist}_maf_0.01
done

# Step 2. Generate SNP info file and multi_info file
/gpfs/gibbs/pi/zhao/zb222/reference/MEGA/generate_snpinfo.sh "${Reference}" "${snplist}"

module load R/4.2.0-foss-2020b
Rscript /gpfs/gibbs/pi/zhao/zb222/reference/MEGA/generate_multi_snpinfo.R "${Reference}" "${snplist}"

# Step 3. Divide LD blocks

#!/bin/bash
#SBATCH --partition=day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --out=/gpfs/gibbs/pi/zhao/zb222/Log/div_blk_${snplist}%j.log

module load R
Rscript /gpfs/gibbs/pi/zhao/zb222/reference/MEGA/div_blk.R "${Reference}" "${snplist}"

# Step 4. Calculate LD for each block
/gpfs/gibbs/pi/zhao/zb222/reference/MEGA/calc_ld.sh "${Reference}" "${snplist}"

# Step 5. Write LD blocks to h5py file
ancestries=("AFR" "EAS" "SAS" "EUR" "AMR")
for pop in "${ancestries[@]}"; do
    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=${pop}_writeblk
#SBATCH --out=/gpfs/gibbs/pi/zhao/zb222/Log/write_blk_${pop}%j.log
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL

module load miniconda
conda activate deep_learning

python /gpfs/gibbs/pi/zhao/zb222/reference/MEGA/write_ldblk.py ${pop} ${Reference}

EOT
done