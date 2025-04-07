#!/bin/bash
#SBATCH --job-name=submit_minor_pop_PRS_cal
#SBATCH --output=/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/method_calculate/get_PRS/Log/3.27/submit_minor_pop_PRS_cal_%j.log
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mem=1G
#SBATCH -p scavenge,day,week
#SBATCH -t 0:10:00
#SBATCH -c 1 

# Define directories
SCRIPT_DIR="/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/method_calculate/get_PRS/Minor_Pop"
OUTPUT_FILE="${SCRIPT_DIR}/job_params.txt"
selected_annot_dir="/gpfs/gibbs/pi/zhao/zb222/ldsc/pipeline/estimate_enrichment/3.choose_annot/result/selected_annot/binary/baselineLD_v2.3_hm3_mega"

# Define trait and populations
declare -A trait_groups
trait_groups[group1]="logTG"  # "HDL LDL TC logTG"
trait_groups[group2]="Height BMI T2D BrC"

declare -A pop_groups
pop_groups[group1]="SAS AFR AMR"
pop_groups[group2]="EAS AFR"  # AFR

groups_to_run=(group1) 

# Generate combinations
> "$OUTPUT_FILE"
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
                    for param_phi in auto; do
                        echo "${trait},${pop},${k},${param_phi}" >> "$OUTPUT_FILE"
                    done
            done
        done
    done
done


# Count the number of jobs
job_count=$(wc -l < "$OUTPUT_FILE")
echo "Total number of jobs: $job_count"

# Get the number of jobs to run (default to all if not specified)
n_jobs=${1:-$job_count}

# Submit the array job
sbatch --array=1-$n_jobs ${SCRIPT_DIR}/1.1_do_calculate_PRS.sh

echo "Submitted PRS calculation array job with $n_jobs tasks."