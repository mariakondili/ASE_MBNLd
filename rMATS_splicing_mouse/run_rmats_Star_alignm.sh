#!/usr/bin/env bash
#SBATCH -A mbnl_dct
#SBATCH --cpus-per-task=6
#SBATCH -p fast
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load conda
module load python/3.7
module load rmats/4.1.0

rmats_dir="/software/miniconda/envs/rmats-4.1.0/bin"
cores=6
gtf="/genomes/Annotation/mm10/gencode.vM25.annotation.gtf"
workDir="/projects/mbnl_dct/mouse_seq/rMATS_splicingCounts/STAR_alignment"

outDir1="${workDir}/MBNLd_vs_Saline_Splicing/"
outDir2="${workDir}/AAVGFP_vs_Saline_Splicing/"

tmpDir1="${outDir1}/tmp/"
tmpDir2="${outDir2}/tmp/"

mkdir -p ${outDir1}; mkdir -p ${tmpDir1};
mkdir -p ${outDir2}; mkdir -p ${tmpDir2};


group1="${workDir}/Saline_STAR_bam.txt"
group2="${workDir}/MBNLdecoy_STAR_bam.txt"
group3="${workDir}/AAVGFP_STAR_bam.txt"


python ${rmats_dir}/rmats.py --b1 ${group1} --b2  ${group2} \
--gtf ${gtf} --readLength 120 \
--od ${outDir1} --tmp ${tmpDir1} \
-t "paired"  --libType fr-firststrand \
--nthread ${cores} --tstat 6 \
--task both


python ${rmats_dir}/rmats.py --b1 ${group1} --b2  ${group3} \
--gtf ${gtf} --readLength 120 \
--od ${outDir2} --tmp ${tmpDir2} \
-t "paired"  --libType fr-firststrand \
--nthread ${cores} --tstat 6 \
--task both

##>run with:
#$ sbatch run_rmats.sh
#> follow with:
#$  sacct -j <jobID>
