#!/bin/bash
#SBATCH -A mbnl_dct
#SBATCH -p fast
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH -o slurm.align.%j.out
#SBATCH -e slurm.align.%j.err


####------------------------ ALIGNMENT --------------------------------------------###

# Usage:
#   hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]
#>Source doc : https://research.csc.fi/rnaseq-alignment
#>manual     :  https://daehwankimlab.github.io/hisat2/manual/

## INFO :
# The alignment score for a paired-end alignment equals the sum of the alignment scores of the individual mates
#-S <hit>: File to write SAM alignments to. By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console
# -x <hisat2-idx>: The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.ht2 / etc

module load  hisat2/2.1.0
module load samtools

index_dir=/mbnl_dct/annotation_files/hisat2_indexes
data_dir=/mbnl_dct/public_Data_SRA/Fastq
alignm_dir=${data_dir}/Trimmed/HISAT2_Alignment
mkdir -p ${alignm_dir}

## for dir in $(ls ${data_dir}); --> if each pair of .fq is in a folder with sample_name

for s in ${data_dir}/SRR* ; do
    ## s="/../../../SRR11548476_HSA_2"
    sid=$(basename ${s})

    r1=${data_dir}/Trimmed/${sid}_1_trimmed.fastq
    r2=${data_dir}/Trimmed/${sid}_2_trimmed.fastq

    ## unzipped=${r1%".gz"} #--> the filename without ".gz" suffix


    hisat2-align-s --wrapper basic-0 \
    -p 16 -x  ${index_dir}/index_GRCh38 \
    --rna-strandness RF \
    --pen-cansplice 0 \
    --pen-noncansplice 12 \
    --pen-canintronlen G,-8.0,1.0 \
    --pen-noncanintronlen G,-8.0,1.0 \
    --known-splicesite-infile ${index_dir}/GRCh38_splice_sites.txt \
    --min-intronlen 20 \
    --max-intronlen 500000 \
    --summary-file ${alignm_dir}/summary_${id}.txt \
    -1 ${rz1} -2 ${rz2} | samtools view -bS - > ${alignm_dir}/${sid}.bam

    ##> Pipe the output of hisat2 to samtools directly, so no need to create a ".sam" file.
    ## the "-" is for stdout of previous command that gets input in samtools.
    ##> otherwise :
    ## samtools view -bS ${alignm_dir}/${sample_id}.sam >  ${alignm_dir}/${sample_id}.bam

    if [ -s  ${alignm_dir}/${sid}.bam ]; then  ## verify File is created
        ## Compress FASTQ to gain splace
        gzip ${r1}
        gzip ${r2}
    fi;
done
