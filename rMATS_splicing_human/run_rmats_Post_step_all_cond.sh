#!/usr/bin/env bash
#SBATCH -A mbnl_dct
#SBATCH --cpus-per-task=8
#SBATCH -p fast
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load conda
module load rmats/4.1.1

rmats_dir="/software/miniconda/envs/rmats-4.1.1/bin"  # module location
# rmats_dir="/shared/home/mkondili/.conda/pkgs/rmats-4.1.1-py37h2140d24_0/bin"  #local installation
## /shared/ifbstor1/software/miniconda/envs/rmats-4.1.1/rMATS/cp_with_prefix.py
cores=8
gtf="/home/mkondili/genomes/Annotation/hg19/Homo_sapiens.GRCh37.87.CHR.gtf"
# gtf="/shared/bank/homo_sapiens/hg19/gff/gencode.v19.annotation.gtf"
workDir="/projects/mbnl_dct/human_cellline/rMATS_SplicingCounts"

outDir_prep="${workDir}/Prep_Ctrl_vs_DM1_vs_MBNLd"
outDir_post="${workDir}/Post_Ctrl_vs_DM1_vs_MBNLd"


mkdir -p ${outDir_prep}/combo_tmp/
mkdir -p ${outDir_post}

group123="${workDir}/CTRL_and_DM1_and_MBNLd_bam.txt" # one line "," sep,with all replicates of all conditions in the order of the filename

## Add a prefix to files and add them to a common tmp folder for POST task.
#> call cp_with_prefix.py
#python ${rmats_dir}/../rMATS/cp_with_prefix.py "Ctrl_vs_DM1_"   ${outDir_prep}/combo_tmp/  ${outDir_prep}/tmp_dm1/*.rmats
#python ${rmats_dir}/../rMATS/cp_with_prefix.py "Ctrl_vs_MBNLd_" ${outDir_prep}/combo_tmp/  ${outDir_prep}/tmp_dct/*.rmats
#> Deleted  in /combo_tmp/ : _0.rmats, _1.rmats, _2.rmats because refer to Ctrl samples that are already analysed in Ctrl_vs_DM1.


###----POST-treatment --------###
python ${rmats_dir}/rmats.py --b1 ${group123}  \
        --gtf ${gtf} --readLength 101     \
        --od ${outDir_post} --tmp ${outDir_prep}/combo_tmp/  \
        -t "paired"  --libType fr-firststrand   \
        --nthread ${cores} --tstat ${cores} \
        --task post


###------ STATS ------###
mkdir -p ${outDir_post}/Ctrl_vs_DM1_output
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py --new-output-dir ${outDir_post}/Ctrl_vs_DM1_output \
                    --old-output-dir ${outDir_post} --group-1-indices 0,1,2 --group-2-indices 3,4,5

python ${rmats_dir}/rmats.py --od ${outDir_post}/Ctrl_vs_DM1_output --tmp ${outDir_post}/Ctrl_vs_DM1_output/tmp --task stat

mkdir -p ${outDir_post}/Ctrl_vs_MBNLd_output
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py --new-output-dir ${outDir_post}/Ctrl_vs_MBNLd_output \
        --old-output-dir  ${outDir_post}  --group-1-indices 0,1,2 --group-2-indices 6,7,8 #--> should be bam of MBNLd last after Ctrl,DM1.

python ${rmats_dir}/rmats.py --od ${outDir_post}/Ctrl_vs_MBNLd_output --tmp ${outDir_post}/Ctrl_vs_MBNLd_output/tmp  --task stat

mkdir -p ${outDir_post}/DM1_vs_MBNLd_output
python ${rmats_dir}/../rMATS/rMATS_P/prepare_stat_inputs.py --new-output-dir ${outDir_post}/DM1_vs_MBNLd_output \
        --old-output-dir ${outDir_post} --group-1-indices 3,4,5 --group-2-indices 6,7,8

python ${rmats_dir}/rmats.py --od ${outDir_post}/DM1_vs_MBNLd_output --tmp ${outDir_post}/DM1_vs_MBNLd_output/tmp --task stat
