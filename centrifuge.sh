#!/bin/bash -l
#SBATCH --job-name=centrifuge
#SBATCH -t 200:00:00
#SBATCH --ntasks-per-node 5
#SBATCH -N 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=me@email.ie

## Centrifuge commands
./centrifuge/centrifuge -x centrifuge_DB/p_compressed --sample-sheet centrifuge_sample-sheet.tsv
