#!/bin/bash
#SBATCH --job-name=concatenate_minor_pop
#SBATCH --output=/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/method_calculate/getbeta_full/Log/3.27/concatenate_minor_pop_%j.log
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mem=10G
#SBATCH -p scavenge,day,week
#SBATCH -t 6:00:00
#SBATCH -c 1 

# Define directories
IN_DIR="/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/beta/LD_v2.3_panel"
OUT_DIR="/gpfs/gibbs/pi/zhao/zb222/Anotation_pipeline/result/beta_collect/LD_v2.3_panel"
selected_annot_dir="/gpfs/gibbs/pi/zhao/zb222/ldsc/pipeline/estimate_enrichment/3.choose_annot/result/selected_annot/binary/baselineLD_v2.3_hm3_mega"

# Define trait and populations
declare -A trait_groups
trait_groups[group1]="HDL LDL TC logTG"
trait_groups[group2]="Height BMI T2D BrC"
trait_groups[group3]="HDL LDL TC logTG Height BMI T2D BrC"  # Add more groups as needed

# Define population groups
declare -A pop_groups
pop_groups[group1]="EAS SAS AFR AMR"
pop_groups[group2]="EAS AFR"
pop_groups[group3]="EUR EAS AFR"    # Add more groups as needed

# Specify which groups to process
groups_to_run=(group1 group2)  # Modify this array to choose which groups to process e.g. (group1 group2)

aux_pop="EUR"
# Function to combine files
combine_files() {
    local pop=$1
    local k=$2
    local out_file=$3

    # Create header
    # echo -e "rsID\tA1\tauto" > "$out_file"
    echo -e "CHR\tSNP\tBP\tA1\tA2\tBETA" > "$out_file"

    # Combine files
    for chr in {1..22}; do
        input_file="${IN_DIR}/trait_${trait}/tar_${pop}/baselineLD_v2.3_Annot_${k}/${trait}_${aux_pop}_${pop}_JointPRS_baselineLD_v2.3_Annot_${k}_${pop}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"
        tail -n +2 "$input_file" >> "$out_file"
    done
}

# Main loop
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
                # Create output directory
                out_dir="${OUT_DIR}/trait_${trait}/tar_${pop}/baselineLD_v2.3_Annot_${k}"
                mkdir -p "$out_dir"

                # Define output file
                out_file="${out_dir}/${trait}_${pop}_JointPRS_baselineLD_v2.3_Annot_${k}_combined.txt"

                # Combine files
                combine_files "$pop" "$k" "$out_file"

                echo "Combined file created: $out_file"
            done
        done
    done
done
echo "All files have been combined."