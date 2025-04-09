# Down load 1KG geno data (2022-08-04 Byrska-Bishop et al. build 38, 3202 samples, contigs unphased) 
# from PLINK2 website "https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg"

wd="/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all"

# Step 1. Obtain sample IID for each population.
pedigree_info <- read.table("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/all_hg38.psam", header=F, sep="\t")
colnames(pedigree_info) <- c("IID", "PAT", "MAT", "SEX", "SuperPop","Population")
snp_info <- list()
for (pop in c("EUR","EAS","SAS","AFR","AMR")){
  snp_info[[pop]] <- pedigree_info[which(pedigree_info$SuperPop==pop),]
  writeLines(snp_info[[pop]]$IID,paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/",pop,".txt"))
}

# Step 2. Get SNPs that have maf >= 0.01 in at least one population.
module load PLINK/2

# Split dataset into different populations
for pop in EUR EAS AFR SAS AMR; do
sample_list="${wd}/${pop}.txt"
plink2 --pfile ${wd}/all_hg38 'vzs' \
--allow-extra-chr \
--max-alleles 2 \
--keep-allele-order \
--keep $sample_list \
--make-bed \
--out ${wd}/bfile/hg38/${pop}/${pop}_1kg
done

# Extract SNPs that have maf greater or equal than 0.01
for pop in EUR EAS AFR SAS AMR; do
plink2 --bfile ${wd}/bfile/hg38/${pop}/${pop}_1kg \
--allow-extra-chr \
--keep-allele-order \
--maf 0.01 \
--make-bed \
--out ${wd}/bfile/hg38/${pop}/${pop}_1kg_maf_0.01
done

# Step 3. Obtain the union of SNPs across five populations 
pop_list = c("EUR", "AFR", "AMR", "EAS", "SAS")
snp_info_maf_0.01 = list()
for (pop in pop_list){
    snp_info_maf_0.01[[pop]] <- read.table(paste0("/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/hg38/",pop,"/",pop,"_1kg_maf_0.01.bim"), header = F)[,2]
}

snp_info_maf_0.01 <- lapply(snp_info_maf_0.01, as.character)
snp_union = Reduce(union, snp_info_maf_0.01)

writeLines(snp_union, "/gpfs/gibbs/pi/zhao/zb222/Data/geno_data_plink2_all/bfile/snplist_1kg_20M.txt")
