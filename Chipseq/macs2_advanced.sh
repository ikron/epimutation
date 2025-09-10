#!/bin/bash -l
#SBATCH --job-name=macs3
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END
#SBATCH --mail-user=marvilla@jyu.fi

module load macs/2.2.7.1

output=/scratch/project_2000350/genomics/chipseq/macs
control_bed=$output/Control_filterdup.bed

# -------------------------
# Define your samples
# Format: sample_name:extsize
# -------------------------
samples=(
    "anc1_matA_rep1:240"
    "anc1_matA_rep2:288"
    "anc2_matA_rep1:292"
    "anc1_mata_rep1:266"
    "anc2_mata_rep1:273"
    "anc2_matA_rep2:263"
    "G40_L31:271"
    "G5_L11:276"
    "G8_L11:291"
    "G10_L11:288"
    "G20_L11:260"
    "G40_L11:285"
    "G5_L5:283"
    "G8_L5:295"
    "G10_L5:271"
    "G20_L5:283"
    "G40_L5:298"
    "G5_L25:287"
    "G8_L25:259"
    "G10_L25:268"
    "G20_L25:251"
    "G40_L25:283"
    "G5_L31:290"
    "G8_L31:273"
    "G10_L31:274"
    "G20_L31:266"
)

# -------------------------
# Generate ChIP coverage tracks
# -------------------------
for sample in "${samples[@]}"; do
    name="${sample%%:*}"
    extsize="${sample##*:}"
    bed_file="$output/${name}.bed"
    output_bdg="$output/${name}.pileup.bdg"

    echo "Running pileup for $name"
    macs2 pileup -i "$bed_file" -o "$output_bdg" -f BED --extsize "$extsize"
done

# -------------------------
# Build local bias from track control
# -------------------------
declare -A d_bg_extsizes=(
    ["anc1_matA_rep1"]=120 ["anc1_matA_rep2"]=144 ["anc2_matA_rep1"]=146
    ["anc1_mata_rep1"]=133 ["anc2_mata_rep1"]=136 ["anc2_matA_rep2"]=131
    ["G40_L31"]=135 ["G5_L11"]=138 ["G8_L11"]=145 ["G10_L11"]=144
    ["G20_L11"]=130 ["G40_L11"]=142 ["G5_L5"]=141 ["G8_L5"]=147
    ["G10_L5"]=135 ["G20_L5"]=141 ["G40_L5"]=149 ["G5_L25"]=143
    ["G8_L25"]=129 ["G10_L25"]=134 ["G20_L25"]=125 ["G40_L25"]=141
    ["G5_L31"]=145 ["G8_L31"]=136 ["G10_L31"]=137 ["G20_L31"]=133
)

for name in "${!d_bg_extsizes[@]}"; do
    extsize="${d_bg_extsizes[$name]}"
    output_bdg="$output/${name}_d_bg.bdg"
    
    echo "Running d_bg pileup for $name"
    macs2 pileup -i "$control_bed" -B --extsize "$extsize" -f BED -o "$output_bdg"
done

# -------------------------
# SLOCAL background (1k)
# -------------------------
macs2 pileup -i "$control_bed" -B --extsize 500 -f BED -o "$output/1k_bg.bdg"

declare -A norm_factors_1k=(
    ["G40_L31"]=0.271 ["anc1_matA_rep1"]=0.240 ["anc1_matA_rep2"]=0.288
    ["anc2_matA_rep1"]=0.292 ["anc1_mata_rep1"]=0.266 ["anc2_mata_rep1"]=0.273
    ["anc2_matA_rep2"]=0.263 ["G5_L11"]=0.276 ["G8_L11"]=0.291
    ["G10_L11"]=0.288 ["G20_L11"]=0.260 ["G40_L11"]=0.285 ["G5_L5"]=0.283
    ["G8_L5"]=0.295 ["G10_L5"]=0.271 ["G20_L5"]=0.283 ["G40_L5"]=0.298
    ["G5_L25"]=0.287 ["G8_L25"]=0.259 ["G10_L25"]=0.268 ["G20_L25"]=0.251
    ["G40_L25"]=0.283 ["G5_L31"]=0.290 ["G8_L31"]=0.273 ["G10_L31"]=0.274
    ["G20_L31"]=0.266
)

for name in "${!norm_factors_1k[@]}"; do
    factor="${norm_factors_1k[$name]}"
    macs2 bdgopt -i "$output/1k_bg.bdg" -m multiply -p "$factor" -o "$output/${name}_1k_bg_norm.bdg"
done

# -------------------------
# LLCOAL background (10k)
# -------------------------
macs2 pileup -i "$control_bed" -B --extsize 5000 -f BED -o "$output/10k_bg.bdg"

declare -A norm_factors_10k=(
    ["G40_L31"]=0.0271 ["anc1_matA_rep1"]=0.0240 ["anc1_matA_rep2"]=0.0288
    ["anc2_matA_rep1"]=0.0292 ["anc1_mata_rep1"]=0.0266 ["anc2_mata_rep1"]=0.0273
    ["anc2_matA_rep2"]=0.0263 ["G5_L11"]=0.0276 ["G8_L11"]=0.0291
    ["G10_L11"]=0.0288 ["G20_L11"]=0.0260 ["G40_L11"]=0.0285 ["G5_L5"]=0.0283
    ["G8_L5"]=0.0295 ["G10_L5"]=0.0271 ["G20_L5"]=0.0283 ["G40_L5"]=0.0298
    ["G5_L25"]=0.0287 ["G8_L25"]=0.0259 ["G10_L25"]=0.0268 ["G20_L25"]=0.0251
    ["G40_L25"]=0.0283 ["G5_L31"]=0.0290 ["G8_L31"]=0.0273 ["G10_L31"]=0.0274
    ["G20_L31"]=0.0266
)

for name in "${!norm_factors_10k[@]}"; do
    factor="${norm_factors_10k[$name]}"
    macs2 bdgopt -i "$output/10k_bg.bdg" -m multiply -p "$factor" -o "$output/${name}_10k_bg_norm.bdg"
done

# -------------------------
# Compare 1k vs 10k
# -------------------------
for name in "${!norm_factors_1k[@]}"; do
    macs2 bdgcmp -m max -t "$output/${name}_1k_bg_norm.bdg" -c "$output/${name}_10k_bg_norm.bdg" -o "$output/${name}_1k_10k_bg_norm.bdg"
done

