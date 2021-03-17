#!/usr/bin/env bash
#SBATCH -A mbnl_dct
#SBATCH -p fast
#SBATCH --mem-per-cpu 40GB
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


module load nextflow/20.04.1

# run : -c <config_file> run -w <workspace>  <nextflow_pip_script.nf>
nextflow -c /home//MBNL_DCT/config_pip.conf  \
            run -w /mouse_seq/NXF_STAR/PIP_out \
            /mouse_analysis/rnaseq_nxf_pip.nf
