#!/bin/bash
#SBATCH -A mbnl_dct
#SBATCH -p fast
#SBATCH -n 1
#SBATCH --mem=160G #!!! Attention ! Really needs 160G, and not "GB" to run !
#SBATCH -o slurm.idx_hisat.%j.out
#SBATCH -e slurm.idx_hisat.%j.err


####------------------------ ALIGNMENT --------------------------------------------###
module load  hisat2/2.1.0

# Usage:
#  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]
#>Source doc : https://research.csc.fi/rnaseq-alignment
#>manual     :  https://daehwankimlab.github.io/hisat2/manual/

## INFO :
# The alignment score for a paired-end alignment equals the sum of the alignment scores of the individual mates
#-S <hit>: File to write SAM alignments to. By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console
# -x <hisat2-idx>: The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.ht2 / etc

### 1/ Create Index
hsap_genome=/bank/homo_sapiens/GRCh38/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa
hsap_annot=/bank/homo_sapiens/GRCh38/gtf/Homo_sapiens.GRCh38.101.gtf

index_dir=/projects/mbnl_dct/annotation_files/hisat2_indexes
mkdir -p ${index_dir}

# hisat2-build ${hsap_genome} hisat-indexes/hsap_GRCh38
#or :
hisat2_extract_splice_sites.py ${hsap_annot} >  ${index_dir}/GRCh38_splice_sites.txt
hisat2_extract_exons.py ${hsap_annot}        >  ${index_dir}/GRCh38_exons.txt
hisat2-build  --ss ${index_dir}/GRCh38_splice_sites.txt --exon ${index_dir}/GRCh38_exons.txt ${hsap_genome} ${index_dir}/index_GRCh38

## Files that will be created :

# $ ll -h /projects/mbnl_dct/annotation_files/hisat2_indexes
# total 42G
# -rw-rw----+ 1 user user 7.1M Apr  3 11:35 GRCh38_exons.txt
# -rw-rw----+ 1 user user 8.4M Apr  3 11:35 GRCh38_splice_sites.txt
# -rw-rw----+ 1 user user 3.5G Apr  2 20:00 index_GRCh38.10.rf
# -rw-rw----+ 1 user user 361M Apr  2 19:54 index_GRCh38.11.rf
# -rw-rw----+ 1 user user 1.7G Apr  3 15:13 index_GRCh38.1.ht2
# -rw-rw----+ 1 user user 3.4G Apr  2 20:00 index_GRCh38.1.rf
# -rw-rw----+ 1 user user 704M Apr  3 15:13 index_GRCh38.2.ht2
# -rw-rw----+ 1 user user 3.9G Apr  2 20:01 index_GRCh38.2.rf
# -rw-rw----+ 1 user user  12K Apr  3 11:36 index_GRCh38.3.ht2
# -rw-rw----+ 1 user user 3.7G Apr  2 20:01 index_GRCh38.3.rf
# -rw-rw----+ 1 user user 703M Apr  3 11:36 index_GRCh38.4.ht2
# -rw-rw----+ 1 user user 3.0G Apr  2 20:00 index_GRCh38.4.rf
# -rw-rw----+ 1 user user 1.7G Apr  3 15:51 index_GRCh38.5.ht2
# -rw-rw----+ 1 user user 3.2G Apr  2 20:00 index_GRCh38.5.rf
# -rw-rw----+ 1 user user 716M Apr  3 15:51 index_GRCh38.6.ht2
# -rw-rw----+ 1 user user 2.7G Apr  2 19:59 index_GRCh38.6.rf
# -rw-rw----+ 1 user user  14M Apr  3 11:36 index_GRCh38.7.ht2
# -rw-rw----+ 1 user user 4.4G Apr  2 20:02 index_GRCh38.7.rf
# -rw-rw----+ 1 user user 2.7M Apr  3 11:36 index_GRCh38.8.ht2
# -rw-rw----+ 1 user user 4.0G Apr  2 20:02 index_GRCh38.8.rf
# -rw-rw----+ 1 user user 4.4G Apr  2 20:01 index_GRCh38.9.rf
