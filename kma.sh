#!/bin/bash -l
#SBATCH --job-name=kma_Yersinia
#SBATCH -t 120:00:00
#SBATCH -n 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=me@email.ie

# set input and output paths
input_dir=Benchmarking/Generator/Metagenomes/Yersinia
output_dir=Benchmarking/Simulated_Metagenomes/kma_results/Yersinia
kma=Benchmarking/Metagenomes/GalwayWW/kma/kma

# set file extensions
r1_extension=_R1.fastq
r2_extension=_R2.fastq

# set database path
db_path=Benchmarking/compress_ncbi_nt/ncbi_nt

# loop over input files
for r1_path in "${input_dir}"/*"${r1_extension}"; do
  r2_path="${r1_path/${r1_extension}/${r2_extension}}"
  output_basename="$(basename "${r1_path}" "${r1_extension}")"
  
  ## run kma
  "${kma}" -ipe "${r1_path}" "${r2_path}" -o "${output_basename}" -t_db "${db_path}" -1t1 -mem_mode -and -apm f -ef -matrix
done
