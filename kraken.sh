#!/bin/bash -l
#SBATCH --job-name=kraken2_wastewater
#SBATCH -t 100:00:00
#SBATCH --ntasks-per-node 5
#SBATCH -N 4    ## Only set this if you need more than one node for your job, otherwise this line is not required at all
#SBATCH --mail-type=ALL
#SBATCH --mail-user=olateju.idowu@ucdconnect.ie

# give input and output paths
input_dir=/scratch/12355656/Benchmarking/Generator/genome
output_dir=/scratch/12355656/Benchmarking/Simulated_Metagenomes/kraken2_results/

# set file extensions
r1_extension=_R1.fastq
r2_extension=_R2.fastq

# give database path
db_path=/scratch/12355656/Benchmarking/Simulated_Metagenomes/k2_standard_08gb_20230314

# set input and output file paths and names
for r1_path in "${input_dir}"/*"${r1_extension}"; do
  r2_path="${r1_path/${r1_extension}/${r2_extension}}"
  output_basename="$(basename "${r1_path}" "${r1_extension}")"
  output_path="${output_dir}/out/${output_basename}_out.tsv"
  report_path="${output_dir}/report/${output_basename}_report.tsv"
  unclassified_path="${output_dir}/Unclassified/${output_basename}.unclassified#.fastq"
  
  echo "Running kraken2 on ${r1_path} and ${r2_path}..."
  echo "Output file path: ${output_path}"
  echo "Report file path: ${report_path}"
  echo "Unclassified reads file path: ${unclassified_path}"
  
  # run kraken2
  kraken2 --db "${db_path}" \
    --memory-mapping \
    --output "${output_path}" \
    --report "${report_path}" \
    --threads 2 \
    --paired "${r1_path}" "${r2_path}" \
    --unclassified-out "${unclassified_path}"
    
  echo "Done."
done
