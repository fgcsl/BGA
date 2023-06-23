#!/usr/bin/bash

echo -e "\nCreating checkm env and installing CheckM and Assembly-Stats..."
mamba create -n checkm -c bioconda checkm-genome -y
# https://data.ace.uq.edu.au/public/CheckM_databases
# While installing checkm with bioconda, databse is auto installed
mamba install -c bioconda assembly-stats -n checkm -y
