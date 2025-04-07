#!/bin/bash
#SBATCH --job-name=submit_EUR_cal
#SBATCH --output=/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/method_calculate/getbeta_full/Log/3.16/submit_EUR_cal_%j.log
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mem=1G
#SBATCH -p scavenge,day,week
#SBATCH -t 0:10:00
#SBATCH -c 1 

SCRIPT_DIR="/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/method_calculate/getbeta_full/EUR"
OUT_DIR="/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/beta/LD_v2.3_panel"
PARAMS_FILE="${SCRIPT_DIR}/job_params.txt"
selected_annot_dir="/gpfs/gibbs/pi/zhao/zb222/ldsc/pipeline/estimate_enrichment/3.choose_annot/result/selected_annot/binary/baselineLD_v2.3_hm3_mega"

# Define trait groups
declare -A trait_groups
trait_groups[group1]="HDL LDL TC logTG"
trait_groups[group2]="Height BMI T2D BrC"
trait_groups[group3]="HDL LDL TC logTG Height BMI T2D BrC"  # Add more groups as needed

# Define population groups
declare -A pop_groups
pop_groups[group1]="EUR"
pop_groups[group2]="EUR"  # AFR
pop_groups[group3]="EUR EAS AFR"    # Add more groups as needed

# Specify which groups to process
groups_to_run=(group1 group2)  # Modify this array to choose which groups to process e.g. (group1 group2)

# Initialize parameter file
> "$PARAMS_FILE"

# Function to check if file exists
check_file() {
    local file="$1"
    if [[ -f "$file" && -s "$file" ]]; then
        return 0
    else
        return 1
    fi
}


# Generate combinations
for group in "${groups_to_run[@]}"; do
    echo "Processing group: $group"
    
    # Get traits and populations for this group
    traits=(${trait_groups[$group]})
    populations=(${pop_groups[$group]})
    
    echo "Traits: ${traits[*]}"
    echo "Populations: ${populations[*]}"
    
    for trait in "${traits[@]}"; do
        for pop in "${populations[@]}"; do
            selected_annot_file="${selected_annot_dir}/selected_annot_${trait}_${pop}.txt"
            K=$(wc -l < "$selected_annot_file")
            echo "There are $K annotations selected for trait:${trait}; pop:${pop}"
            for k in $(seq 1 $K); do
                mkdir -p "$OUT_DIR/trait_${trait}/tar_${pop}/baselineLD_v2.3_Annot_${k}"
                for chr in {1..22}; do
                    for param_phi in auto; do
                        output_file="$OUT_DIR/trait_${trait}/tar_${pop}/baselineLD_v2.3_Annot_${k}/${trait}_${pop}_JointPRS_baselineLD_v2.3_Annot_${k}_${pop}_pst_eff_a1_b0.5_phi${param_phi}_chr${chr}.txt"
                        # Check if output file exists, if not, append to job_params.txt
                        if ! check_file "$output_file"; then
                            echo "${trait},${pop},${k},${chr},${param_phi}" >> "$PARAMS_FILE"
                        fi
                    done
                done
            done
        done
    done
done


# Count the number of jobs
job_count=$(wc -l < "$PARAMS_FILE")
echo "Total number of jobs: $job_count"

# Get the number of jobs to run (default to all if not specified)
n_jobs=${1:-$job_count}

# Submit the array job
sbatch --array=1-$n_jobs ${SCRIPT_DIR}/1.1_get_beta.sh

echo "Submitted beta calculation array job with $n_jobs tasks."