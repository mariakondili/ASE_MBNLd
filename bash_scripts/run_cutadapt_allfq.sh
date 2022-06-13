#!/bin/bash
#SBATCH -A mbnl_dct
#SBATCH -p fast
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH -o slurm.trim.%j.out
#SBATCH -e slurm.trim.%j.err

##> Can verify with FASTQC what adapter is found in the file before applying trimming tool.
module load fastqc/0.11.9
# fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
#            [-c contaminant file] seqfile1 .. seqfileN

input_dir=/mbnl_dct/public_Data_SRA/Fastq

module load cutadapt/2.10

####-----------------------CUTADAPT----------------------------------------###


mkdir -p ${input_dir}/Trimmed/

sampleslist=$(ls ${input_dir}/)

for d in ${sampleslist[@]}; do
    r1=${input_dir}/${d}/*_1.fastq
    r2=${input_dir}/${d}/*_2.fastq

    # Universal adapter Illumina to be removed = AGATCGGAAGAG
    cutadapt --cores 16 -B AGATCGGAAGAG --minimum-length 60 --report  minimal  \
    --output ${input_dir}/Trimmed/${d}_1_trimmed.fastq \
    --paired-output  ${input_dir}/Trimmed/${d}_2_trimmed.fastq \
    -Z ${r1} ${r2}

done

##-------------- FAST_QC --------------------##

## run for a pair of fastq after trimming ,to verify adapters have been removed
fastqc_dir=${input_dir}/../Quality_Check/
mkdir -p ${fastqc_dir}

java_dir=/software/miniconda/envs/fastqc-0.11.9/bin/java
#java_dir to be found with command "$ which java", when fastqc is loaded
fastqc -o ${fastqc_dir} -j ${java_dir} \
${input_dir}/Trimmed/SRR11548475_1_trimmed.fastq  \
${input_dir}/Trimmed/SRR11548475_2_trimmed.fastq
