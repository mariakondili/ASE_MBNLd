#!/usr/bin/env bash
#SBATCH -A mbnl_dct
#SBATCH --cpus-per-task=16
#SBATCH -p fast
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load conda
module load python/3.7
module load rmats/4.1.0

rmats_dir="/shared/ifbstor1/software/miniconda/envs/rmats-4.1.0/bin"
cores=16
gtf="/shared/home/mkondili/genomes/Annotation/mm10/gencode.vM25.annotation.gtf"
workDir="/shared/projects/mbnl_dct/mouse_seq/rMATS_splicingCounts/"
outDir="${workDir}/Saline_vs_AAvGfp_Splicing_firststrand/"
tmpDir="${outDir}/tmp/"
mkdir -p ${outDir}
mkdir -p ${tmpDir}


group1="${workDir}/Saline_bam.txt"
group2="${workDir}/AAVGFP_bam.txt"

#for readLength (of mapped reads), find out from bam file :
# samtools view -F 4 file.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n | uniq -c
# 75-150 bp length

python ${rmats_dir}/rmats.py --b1 ${group1} --b2  ${group2} \
--gtf ${gtf} --readLength 120 \
--od ${outDir} --tmp ${tmpDir} \
-t "paired"  --libType fr-firststrand \
--nthread ${cores} --tstat ${cores} \
--task both

##>run with:
#$ sbatch run_rmats.sh
#> follow with:
#$  sacct -j <jobID>
