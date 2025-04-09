# Function 1.(generate_snpinfo.sh) Generate snp_info for each population 

#!/bin/bash
#SBATCH --job-name=generate_snpinfo
#SBATCH --output=/gpfs/gibbs/pi/zhao/zb222/Log/snpinfo%A_%j.log
#SBATCH --partition=day
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --array=0-4

# Load module if necessary
module load PLINK/1
# Choose current reference panel (e.g. hapmap3 SNP list)
REFERENCE=$1
SNPLIST=$2
# Define an array of populations
POPULATIONS=("AFR" "EAS" "SAS" "EUR" "AMR")
POP=${POPULATIONS[$SLURM_ARRAY_TASK_ID]}
POP_lower=${POP,,}
# Directory containing the PLINK files
INPUT_DIR="/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19"
OUTPUT_DIR="/gpfs/gibbs/pi/zhao/zb222/reference/${REFERENCE}/ldblk_1kg_${POP_lower}"
# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# File paths
BFILE="${INPUT_DIR}/${POP}/${POP}_${SNPLIST}_maf_0.01"
OUTPUT_FILE="${OUTPUT_DIR}/snpinfo_1kg_hm3"
FRQ_FILE="${OUTPUT_DIR}/maf_${POP}.frq"

# Calculate MAF using PLINK
plink --bfile ${BFILE} --allow-extra-chr --freq --out ${OUTPUT_DIR}/maf_${POP}

# Create output file and write header
echo -e "CHR\tSNP\tBP\tA1\tA2\tMAF" > ${OUTPUT_FILE}

# Merge MAF data with .bim data
awk 'NR==FNR{a1[$2]=$3; a2[$2]=$4; maf[$2]=$5; next; next} {print $1 "\t" $2 "\t" $4 "\t" a1[$2] "\t" a2[$2] "\t" maf[$2]}' ${FRQ_FILE} ${BFILE}.bim >> ${OUTPUT_FILE}
cp ${OUTPUT_FILE} /gpfs/gibbs/pi/zhao/zb222/reference/${REFERENCE}/snpinfo_1kg_hm3_${POP}

rm ${FRQ_FILE}


# Function 2.(generate_multi_snpinfo.R) Generate multi_snpinfo 

args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
wd=paste0("/gpfs/gibbs/pi/zhao/zb222/reference/",as.character(args[1]))
snplist=as.character(args[2])
setwd(wd) 
options(scipen = 999)
# Create a union of all SNPs
library(dplyr)

# Example populations data frames (to be replaced with actual read operations)
snpinfo_1kg_AFR <- read.table(paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19/AFR/AFR_",snplist,"_maf_0.01.bim"), header = F)[,-3]
snpinfo_1kg_AMR <- read.table(paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19/AMR/AMR_",snplist,"_maf_0.01.bim"), header = F)[,-3]
snpinfo_1kg_EAS <- read.table(paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19/EAS/EAS_",snplist,"_maf_0.01.bim"), header = F)[,-3]
snpinfo_1kg_EUR <- read.table(paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19/EUR/EUR_",snplist,"_maf_0.01.bim"), header = F)[,-3]
snpinfo_1kg_SAS <- read.table(paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19/EAS/EAS_",snplist,"_maf_0.01.bim"), header = F)[,-3]
df <- full_join(snpinfo_1kg_AFR,snpinfo_1kg_AMR,by=c("V1","V2","V4","V5","V6"))
df <- df %>% full_join(snpinfo_1kg_EAS,by=c("V1","V2","V4","V5","V6"))
df <- df %>% full_join(snpinfo_1kg_EUR,by=c("V1","V2","V4","V5","V6"))
df <- df %>% full_join(snpinfo_1kg_SAS,by=c("V1","V2","V4","V5","V6"))
# Function to merge and calculate frequencies and flip status
colnames(df) <- c("CHR","SNP","BP","A1","A2")

maf_AFR <- read.table("snpinfo_1kg_hm3_AFR",header = T)
maf_AMR <- read.table("snpinfo_1kg_hm3_AMR",header = T)
maf_EAS <- read.table("snpinfo_1kg_hm3_EAS",header = T)
maf_EUR <- read.table("snpinfo_1kg_hm3_EUR",header = T)
maf_SAS <- read.table("snpinfo_1kg_hm3_SAS",header = T)

for (pop in c("AFR","AMR","EAS","EUR","SAS")){
  maf_pop <- get(paste0("maf_", pop))
  df <- full_join(df,maf_pop,by=c("CHR","SNP","BP"),suffix = c("","_pop"))
  df <- df %>% 
    mutate(
      FRQ = case_when(
        A1 == A1_pop & A2 == A2_pop ~ MAF ,
        A1 == A2_pop & A2 == A1_pop ~ 1 - MAF ,
        TRUE ~ NA_real_
      )
    )
  names(df)[names(df) == "FRQ"] <- paste0("FRQ_", pop)
  df <- df %>%
    select(-A1_pop, -A2_pop, -MAF)
}
df[is.na(df)]=0
p <- nrow(df)

for (pop in c("AFR","AMR","EAS","EUR","SAS")){
  df <- cbind(df, setNames(data.frame(rep(1, p)), paste0("FLP_", pop)))
}
df <- df[order(df$CHR,df$BP),]
write.table(df, "snpinfo_mult_1kg_hm3", sep = "\t", row.names = FALSE, quote = F)

# Function 3.(div_blk.R) Divide LD blocks
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
Reference=as.character(args[1])
snplist=as.character(args[2])
setwd("/gpfs/gibbs/pi/zhao/zb222/reference/MEGA")

for (pop in c("EUR", "AFR", "EAS", "SAS", "AMR")){
	if (pop == "EUR"){
	a = read.table("EUR_block.txt", header=T)
	b = read.table(paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19/",pop,"/",pop, "_",snplist,"_maf_0.01.bim"))
	}else if (pop %in% c("SAS", "EAS")){
		a = read.table("EAS_block.txt", header=T)
		b = read.table(paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19/",pop,"/",pop, "_",snplist,"_maf_0.01.bim"))
	}else if(pop %in% c("AFR", "AMR")){
		a = read.table("AFR_block.txt", header=T)
		b = read.table(paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19/",pop,"/",pop, "_",snplist,"_maf_0.01.bim"))
	}


size = c()
chrom = c()
n = 1
for (i in 1:nrow(a)) {
    chr = as.numeric(sub("chr","",a[i,1]))
    idx = which(b$V4 >= as.numeric(a[i,2]) & b$V4 < as.numeric(a[i,3]) & b$V1 == chr)
    if (length(idx) > 0) {
	write.table(b[idx,2], file=paste0("/gpfs/gibbs/pi/zhao/zb222/reference/",Reference,"/snplist_ldblk/",pop,"/",n), 
		    row.names=F, col.names=F, quote=F, append=F)
	size[n] = length(idx)
	chrom[n] = chr
	n = n + 1
    }
}

write.table(size, file=paste0("/gpfs/gibbs/pi/zhao/zb222/reference/",Reference,"/ldblk/",pop,"/blk_size"), 
	row.names=F, col.names=F, quote=F, append=F)

write.table(chrom, file=paste0("/gpfs/gibbs/pi/zhao/zb222/reference/",Reference,"/ldblk/",pop,"/blk_chr"), 
	row.names=F, col.names=F, quote=F, append=F)

}

# Function 4.(calc_ld.sh) Calculate LD for each block in each population
#!/bin/bash
Reference=$1
snplist=$2

ancestries=("AFR" "EAS" "SAS" "EUR" "AMR") 
for pop in "${ancestries[@]}"; do
  # Define input and output directories
  input_dir="/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg19"
  output_dir="/gpfs/gibbs/pi/zhao/zb222/reference/${Reference}/ldblk/${pop}"
  snp_list_dir="/gpfs/gibbs/pi/zhao/zb222/reference/${Reference}/snplist_ldblk/${pop}"

  # Submit the processing script as a job
  sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=${pop}_panel_${Reference}
#SBATCH --output=/gpfs/gibbs/pi/zhao/zb222/Log/${pop}_panel_${Reference}%j.log
#SBATCH --partition=day
#SBATCH --requeue
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL

/gpfs/gibbs/pi/zhao/zb222/reference/MEGA/process_pop.sh "${pop}" "${input_dir}" "${snp_list_dir}" "${output_dir}" "${snplist}"

EOT
done
# Function 5.(process_pop.sh) Run the LD calculation within each population.
#!/bin/bash
# Usage: process_pop.sh <pop> <input_dir> <snp_list_dir> <output_dir> <snplist>

pop=$1
input_dir=$2
snp_list_dir=$3
output_dir=$4
snplist=$5

module load PLINK/1.9b_6.21-x86_64

# Iterate through each file in the snp_list_dir
for file in "$snp_list_dir"/*; do
  # Extract the block number from the filename
  blk=$(basename "$file")

  # Define the extract file and output file
  extract_file="$file"  # Use the full path directly
  output_file="${output_dir}/ldblk${blk}_1kg"

  # Run PLINK command, no keep-allel-order
  plink --bfile "${input_dir}/${pop}/${pop}_${snplist}_maf_0.01" \
        --keep-allele-order \
        --extract "${extract_file}" \
        --r square \
        --out "${output_file}"
done

# Function 6.(write_ldblk.py) Write LD blocks to h5py file
#!/usr/bin/python
"""
Read LD blocks and write as hdf5
"""

import numpy as np
import os
import h5py
import sys

pop = sys.argv[1]
Reference = sys.argv[2]

os.chdir('/gpfs/gibbs/pi/zhao/zb222/reference/' + Reference)

OUT_DIR = './ldblk_1kg_' + pop.lower()
os.makedirs(OUT_DIR, exist_ok=True)

with open('./ldblk/'+pop+'/blk_chr') as ff:
        blk_chr = [int(line.strip()) for line in ff]
with open('./ldblk/'+pop+'/blk_size') as ff:
        blk_size = [int(line.strip()) for line in ff]

n_blk = len(blk_chr)

n_chr = max(blk_chr)

for chrom in range(1,n_chr+1):
    print('... parse chomosome %d ...' % chrom)
    chr_name = OUT_DIR + '/ldblk_1kg_chr' + str(chrom) + '.hdf5'
    hdf_chr = h5py.File(chr_name, 'w')
    blk_cnt = 0
    for blk in range(n_blk):
        if blk_chr[blk] == chrom:
            if blk_size[blk] > 0:
                blk_file = './ldblk/' + pop + '/ldblk' + str(blk+1) + '_1kg.ld'
                with open(blk_file) as ff:
                    ld = [[float(val) for val in (line.strip()).split()] for line in ff]
                print('blk %d size %s' % (blk+1, np.shape(ld)))
                snp_file = './snplist_ldblk/' + pop + '/' + str(blk+1)
                with open(snp_file) as ff:
                    snplist = [line.strip() for line in ff]
            else:
                ld = []; snplist = []
            blk_cnt += 1
            hdf_blk = hdf_chr.create_group('blk_%d' % blk_cnt)
            hdf_blk.create_dataset('ldblk', data=np.array(ld), compression="gzip", compression_opts=9)
            hdf_blk.create_dataset('snplist', data=snplist, compression="gzip", compression_opts=9)

