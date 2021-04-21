#!/usr/bin/env bash
#SBATCH -A mbnl_dct
#SBATCH -p fast
#SBATCH --mem-per-cpu 40GB
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


module load nextflow/20.04.1
cd /home/MBNL_DCT/
mkdir -p rnaseq_pip
# run : -c <config_file> run -w <workspace>  <nextflow_pip_script.nf>
# Changed in 2021: nextflow run -w nxf_workspace rnaseq_pip/
#..where rnaseq_pip/ contains : nextflow.config , main.nf, and fastq_samplesheet.tsv
# nxf_workspace : directory will be created by Nextflow, and keep tmp files of run.
# attention: main.nf and nextflow.config should always be named this way.

nextflow run -w nxf_workspace rnaseq_pip
