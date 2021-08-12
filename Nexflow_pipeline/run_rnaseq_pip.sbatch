#!/usr/bin/env bash
#SBATCH -A mbnl_dct
#SBATCH -p fast
#SBATCH --cpus-per-task=20
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


module load nextflow/20.04.1

cd /mbnl_dct/mouse_seq/Nextflow_pipeline/

##> All used files must be in a dir named after the project or Pipeline ,here :
###  mkdir rnaseq_pip
#1/ config file must be in same dir as .nf file ,named "nextflow.config"
#2/ .nf script should be named "main.nf"

###  mv rnaseq_pip.conf    rnaseq_pip/nextflow.config
###  mv rnaseq_nxf_pip.nf  rnaseq_pip/main.nf

mkdir -p PIP_out/ # workspace with temp.files of Nextflow processes 
nextflow  run -w PIP_out  rnaseq_pip
