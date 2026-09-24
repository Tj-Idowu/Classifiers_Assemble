#!/bin/bash -l
#SBATCH --job-name=CLARK_zero
#SBATCH -t 150:00:00
#SBATCH --nodelist=sonicmem3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=me@email.ie

## Clean up
./clean.sh

## Update Taxonomy
./updateTaxonomy.sh


## Create and index the database
./set_targets.sh CLARK_DB bacteria --species

## Classify the samples
./classify_metagenome.sh -P ../zero_F.txt ../zero_R.txt -R ../zero_results.txt
