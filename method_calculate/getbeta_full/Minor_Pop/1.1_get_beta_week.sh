#!/bin/bash
#SBATCH --job-name=beta_minor_pop_cal_week
#SBATCH --partition=week
#SBATCH --requeue
#SBATCH --output=/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/method_calculate/getbeta_full/Log/3.25/beta_minor_pop_cal_week_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mem=30G
#SBATCH --time=30:00:00

module load miniconda
# source /vast/palmer/apps/avx2/software/miniconda/24.3.0/etc/profile.d/conda.sh
conda activate JointPRS

cd /gpfs/gibbs/pi/zhao/lx94/JointPRS_extend/data/

declare -A sample_sizes=(
    ["HDL,EUR"]=875625 ["HDL,EAS"]=117296 ["HDL,AFR"]=89220 ["HDL,SAS"]=33953 ["HDL,AMR"]=39602
    ["LDL,EUR"]=832800 ["LDL,EAS"]=81417 ["LDL,AFR"]=86256 ["LDL,SAS"]=33658 ["LDL,AMR"]=35885
    ["TC,EUR"]=918113 ["TC,EAS"]=144819 ["TC,AFR"]=90873 ["TC,SAS"]=34135 ["TC,AMR"]=41222
    ["logTG,EUR"]=853550 ["logTG,EAS"]=82722 ["logTG,AFR"]=87939 ["logTG,SAS"]=34023 ["logTG,AMR"]=37273
    ["Height","EUR"]=251955 ["Height","EAS"]=159095 ["Height","AFR"]=49781
    ["BMI","EUR"]=233634 ["BMI","EAS"]=158284 ["BMI","AFR"]=49335
    ["T2D","EUR"]=193440 ["T2D","EAS"]=118493 ["T2D","AFR"]=38919
    ["BrC","EUR"]=245620 ["BrC","EAS"]=27138 ["BrC","AFR"]=7434
)

OUT_DIR="/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/beta/LD_v2.3_panel"
mkdir -p "$OUT_DIR"
REF_DIR="/gpfs/gibbs/pi/zhao/zb222/reference/BaselineLD_v2.3_hm3"  # Union of all BaselineLD_v2.3 annotation snps and hm3 snps
SCRIPT_DIR="/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/method_calculate/getbeta_full/Minor_Pop"
RUN_SCRIPT="/gpfs/gibbs/pi/zhao/zb222/JointPRS_original/JointPRS.py"
aux_pop="EUR"  # define the auxiliary population

IFS=',' read -r trait pop k chr param_phi < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SCRIPT_DIR}/job_params_week.txt)
echo "Processing: Trait: $trait, Target Population: $pop, Chromosome: $chr, phi: $param_phi"

mkdir -p "$OUT_DIR/trait_${trait}/tar_${pop}/baselineLD_v2.3_Annot_${k}"

sample_size=${sample_sizes["$trait,$pop"]}

start=$(date +%s)

python $RUN_SCRIPT \
--ref_dir=$REF_DIR \
--bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal \
--sst_file=/gpfs/gibbs/pi/zhao/zb222/Data/summary_data/PRS-CSx/BaselineLD_v2.3_hm3/${trait}_${pop}/auxiliary_pop_${aux_pop}/${chr}_JointPRS_baselineLD_v2.3_Annot_${k}.txt,/gpfs/gibbs/pi/zhao/zb222/Data/summary_data/PRS-CSx/BaselineLD_v2.3_hm3/${trait}_${pop}/${chr}_JointPRS_baselineLD_v2.3_Annot_${k}.txt \
--n_gwas=${sample_sizes["$trait","$aux_pop"]},${sample_size} \
--rho_cons=1,1,1,1,1 \
--chrom=${chr} \
--pop=${aux_pop},${pop} \
--out_dir=$OUT_DIR/trait_${trait}/tar_${pop}/baselineLD_v2.3_Annot_${k} \
--out_name=${trait}_${aux_pop}_${pop}_JointPRS_baselineLD_v2.3_Annot_${k}

end=$(date +%s)
duration=$((end - start))
duration_hours=$(echo "$duration / 3600" | bc -l)
# Print the duration in hours, formatted to 2 decimal places
printf "Job took %.2f hours to complete.\n" $duration_hours

echo "End of the Job."