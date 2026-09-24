#!/bin/bash -l
#SBATCH --job-name=blast_patho
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node 5
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=me@email.ie

# set input and output paths
input_dir=Benchmarking/Simulated_Metagenomes/FASTA
output_dir=Benchmarking/Simulated_Metagenomes/BLAST_OUT/Patho
db_dir=Benchmarking/Simulated_Metagenomes/BLAST_DB/all_genomes

# Set file extensions
in_ext=.fasta
out_ext=.tab

# Loop over input files
for files in "${input_dir}"/*; do
        out_basename="$(basename "${files}" "${in_ext}")"
        out_path="${output_dir}"/"${out_basename}${out_ext}"
        
        # Do the analysis
        blastn -query "${files}" \
        -db "${db_dir}" \
        -num_threads 8 \
        -perc_identity 98 \
        -qcov_hsp_perc 100 \
        -out "${out_path}" \
        -outfmt 7

        echo "Done."
done
