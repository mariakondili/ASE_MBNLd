#!/bin/bash

##> script to run in Linux
samplelist="/shared/home/mkondili/MBNL_DCT/mouse_analysis/fastq_samplesheet_mouse.tsv"
outDir="/shared/projects/mbnl_dct/mouse_seq/Aligned"
data_dir="/shared/projects/mbnl_dct/mouse_seq/fastq/"
refGenome="/shared/home/mkondili/genomes/Indexes/STAR_2.7.5/mm10"
ref_Gtf="/shared/home/mkondili/genomes/Annotation/mm10/gencode.vM25.annotation.gtf" ## latest, otherwise : vM19
cpus=8

## Change "-"(hyphen) to "_"(underscore),which allows better manip of bash cmd
#  for f in `ls  ${data_dir}*.fastq.gz` ;do
#    new=`echo $f |tr '-' '_'`
#    mv $f $new
# done
## Read the samplelist one-by-one line :
IFS=$'\n'
for line in $(cat $samplelist);do
  # extract Read.1,Read.2 ,split by space :
  r1=$(echo $line | cut -d ' ' -f1)
  r2=$(echo $line | cut -d ' ' -f2)
  suffix="_R1.fastq.gz"
  nm=$(basename $r1 $suffix)
  mkdir -p ${outDir}/${nm}
  echo "R1=${r1}";
  echo "R2=${r2}";
  STAR --runMode  alignReads --runThreadN ${cpus} --genomeDir ${refGenome} --sjdbGTFfile ${ref_Gtf} --readFilesIn ${r1}  ${r2} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outDir}/${nm}/${nm}_   --readFilesCommand zcat
  echo  "=> Alignment of ${nm} completed\n"
done;
echo "=> Alignment totally completed !\n"
