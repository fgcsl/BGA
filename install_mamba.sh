#!/bin/bash

echo "Downloading Mambaforge..."
cd Downloads
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
echo "Installing Mambaforge..."
chmod +x Mambaforge-Linux-x86_64.sh
./Mambaforge-Linux-x86_64.sh
echo -e "\nAdding Bioconda Channel..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible
echo -e "\nMamba installed\n"
