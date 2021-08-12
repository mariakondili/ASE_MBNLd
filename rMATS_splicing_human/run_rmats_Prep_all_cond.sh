#!/usr/bin/env bash
#SBATCH -A mbnl_dct
#SBATCH --cpus-per-task=8
#SBATCH -p fast
#SBATCH -o slurm.%j.out

module load conda

module load rmats/4.1.1
rmats_dir="/shared/ifbstor1/software/miniconda/envs/rmats-4.1.1/bin"  # module location
# rmats_dir="/shared/home/mkondili/.conda/pkgs/rmats-4.1.1-py37h2140d24_0/bin"  #local installation

cores=8
#gtf="/shared/home/mkondili/genomes/Annotation/hg19/Homo_sapiens.GRCh37.87.CHR.gtf"
gtf="/shared/bank/homo_sapiens/hg19/gff/gencode.v19.annotation.gtf"
workDir="/shared/projects/mbnl_dct/human_cellline/rMATS_SplicingCounts"
outDir_prep="${workDir}/Prep_Ctrl_vs_DM1_vs_MBNLdecoy"


mkdir -p ${outDir_prep}

group1="${workDir}/Control_bam.txt"
group2="${workDir}/DM1_bam.txt"
group3="${workDir}/MBNLdecoy_bam.txt"

#for readLength (of mapped reads), find out from bam file :
# samtools view -F 4 file.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n | uniq -c
# 75-150 bp length

###--- PREP CTRL-vs-DM1 ----###
mkdir -p ${outDir_prep}/tmp_dm1/
python ${rmats_dir}/rmats.py --b1 ${group1} --b2  ${group2} \
--gtf ${gtf} --readLength 101 \
--od ${outDir_prep} --tmp ${outDir_prep}/tmp_dm1/ \
-t "paired"  --libType fr-firststrand \
--nthread ${cores} --tstat ${cores} \
--task prep


###--- PREP CTRL-vs-MBNld ----###
mkdir -p ${outDir_prep}/tmp_dct/
python ${rmats_dir}/rmats.py --b1 ${group1} --b2  ${group3} \
--gtf ${gtf} --readLength 101 \
--od ${outDir_prep} --tmp ${outDir_prep}/tmp_dct/  \
-t "paired"  --libType fr-firststrand \
--nthread ${cores} --tstat ${cores} \
--task prep
