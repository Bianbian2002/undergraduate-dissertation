#!/bin/bash
#SBATCH --job-name=calculate_EUR_PRS_%j
#SBATCH --output=/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/method_calculate/get_PRS/Log/3.17/calculate_EUR_PRS_%A_%a.log
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mem=30G
#SBATCH -p scavenge,day,week
#SBATCH -t 12:00:00
#SBATCH -c 1

module load PLINK/2

SCRIPT_DIR="/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/method_calculate/get_PRS/EUR"
OUT_DIR="/gpfs/gibbs/pi/zhao/zb222/GWAS_subsample/Project3_331/result/PRS_score/LD_v2.3_panel"

# Read parameters from job_params.txt
IFS=',' read -r trait pop k param_phi < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SCRIPT_DIR}/job_params.txt)

echo "Processing: Trait: $trait, Population: $pop, k: $k, param_phi: $param_phi"

# Define input and output files
INPUT_FILE="/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/beta_collect/LD_v2.3_panel/trait_${trait}/tar_${pop}/baselineLD_v2.3_Annot_${k}/${trait}_${pop}_JointPRS_baselineLD_v2.3_Annot_${k}_combined.txt"
OUTPUT_PREFIX="${OUT_DIR}/trait_${trait}/tar_${pop}/baselineLD_v2.3_Annot_${k}/${trait}_${pop}_JointPRS_baselineLD_v2.3_Annot_${k}"

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_PREFIX")"

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop}_10K \
       --double-id \
       --threads 1 \
       --score ${INPUT_FILE} 2 4 6 \
       --out ${OUTPUT_PREFIX}


echo "This is the end of the job."