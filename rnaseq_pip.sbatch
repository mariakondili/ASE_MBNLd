#!/usr/bin/env bash
#SBATCH -A mbnl_dct
#SBATCH -p fast
#SBATCH --cpus-per-task=20
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


module load nextflow/20.04.1

cd /shared/projects/mbnl_dct/mouse_seq/Nextflow_STAR/

### ---> Change of Command to run Nextflow since 2021 !!!

#1 / config file must be in same dir as .nf file ,named "nextflow.config"
# 2/ .nf script should be named "main.nf"
###  mkdir rnaseq_pip
###  mv config_pip.conf rnaseq_pip/nextflow.config
###  mv rnaseq_nxf_pip.nf  rnaseq_pip/main.nf
# 3/ above files must be in a dir named after the project or Pipeline ,for this one : rnaseq_pip/

# config="/shared/home/mkondili/MBNL_DCT/mouse_analysis/NXF_files/config_pip.conf"
nextflow  run -w PIP_out  rnaseq_pip
