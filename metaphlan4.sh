#!/bin/bash -l
#SBATCH --job-name=metaphlan_analysis
#SBATCH -t 120:00:00
#SBATCH --ntasks-per-node 5
#SBATCH -n 8    ## This is the number of cores that you can specify as --p-threads, and they should be the same
#SBATCH --mail-type=ALL
#SBATCH --mail-user=olateju.idowu@ucdconnect.ie


# Set up
input_dir=/scratch/12355656/Benchmarking/Generator/Metagenomes/Bacillus
output_dir=/scratch/12355656/Benchmarking/Simulated_Metagenomes/metaphlan_results/Bacillus
database=/scratch/12355656/Benchmarking/Simulated_Metagenomes/metaphlan_databases
in_ext=_R1.fastq.gz
out_ext=.txt
bowtie_ext=.bowtie2.bz2

# Loop over input files
for files in "${input_dir}"/*f2*_R1.fastq.gz; do

        out_basename="$(basename "${files}" "${in_ext}")"
        out_path="${output_dir}"/"${out_basename}${out_ext}"
        # Do the analysis
        metaphlan "${files}" --bowtie2out "${output_dir}"/"${out_basename}${bowtie_ext}" --input_type fastq \
        --unclassified_estimation --index mpa_vOct22_CHOCOPhlAnSGB_202212 --bowtie2db "${database}" \
        -o "${out_path}"
done
