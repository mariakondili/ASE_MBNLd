#!/usr/bin/env nextflow

//
// ~~~ Subject :
// Align sequences of mouse for MBNL_DCT, 3 conditions: Saline, AAV, Mbnl-Delta
// Count genes, and Exons for DESeq2, and DEXSeq analysis
// ~~~

// Author: Maria Kondili
//Date : 18 Nov 20


workflow.onComplete {
    //any worlflow property can be used here
    if ( workflow.success ) {
        println "Pipeline Complete"
    }
    println "Command line: $workflow.commandLine"
}

workflow.onError {
    println "Oops .. something went wrong"
}

// when cutadapt done: use this samplesheet for alignment
// params.samplelist = "/shared/home/mkondili/MBNL_DCT/mouse_analysis/NXF_files/cutadapt_samplesheet_mouse.tsv"
params.samplelist  = "/shared/home/mkondili/MBNL_DCT/mouse_analysis/NXF_files/fastq_samplesheet_mouse.tsv"
params.data_dir  = "/shared/projects/mbnl_dct/mouse_seq/fastq/"
params.cutadapt_dir = "/shared/projects/mbnl_dct/mouse_seq/Nextflow_STAR/cut_adapted"
params.outDir= "/shared/projects/mbnl_dct/mouse_seq/Nextflow_STAR/Aligned_Encode"
params.countsDir= "/shared/projects/mbnl_dct/mouse_seq/Nextflow_STAR/HTSeq_counts_Encode"
//params.dexseqDir   = "/shared/projects/mbnl_dct/mouse_seq/NXF_STAR/DEXSeq_counts"
params.refGenome = "/shared/bank/mus_musculus/mm10/star-2.7.5a"
params.ref_Gtf = "/shared/home/mkondili/genomes/Annotation/mm10/gencode.vM25.annotation.gtf"
//params.dexseq_gtf="/shared/home/mkondili/genomes/Annotation/mm10/Galaxy_produced_DEXSeq_annot_mm10.gtf"
params.cpus = 20

// create pairs as variables
samplist = file("${params.samplelist}")
reader = samplist.newReader()
// read the samplesheet and create array with two reads to Input in next Process
pair=[]
samplist.withReader {
    String line
    while ( line = reader.readLine() ) {
        String r1 = line.split(" ")[0]  // Problem when using "\t" separator !
        println "read_1 = ${r1}"
        String r2 = line.split(" ")[1]
        println "read_2 = ${r2}"
        String bn =  r1.split("${params.data_dir}")[1].split("_R1")[0]
        // String bn =  r1.split("${params.data_dir}")[1].split("/")[1].split("_R1")[0]
        println "Base-name = ${bn}"
        pair.add([bn,r1,r2])
    }
}

inputChannel  = Channel.fromList(pair)



process Cut_Adapters {

         cpus "${params.cpus}"
         memory "40G"
         module "cutadapt/2.10"


         input:
         tuple val(bn),path(r1),path(r2) from inputChannel

         output:
         tuple val(bn),path("${bn}_R1_cutadapt.fastq.gz"),path("${bn}_R2_cutadapt.fastq.gz") into cutadChannel

         shell:
         """
         module li ;
         mkdir -p !{params.cutadapt_dir}/!{bn}/

         cutadapt --cores !{params.cpus} \
         -B AGATCGGAAGAG \
         --minimum-length 75 --error-rate 0.1 \
         --report  minimal \
         --output !{bn}_R1_cutadapt.fastq.gz \
         -p       !{bn}_R2_cutadapt.fastq.gz \
         -Z !{r1} !{r2}

         cp !{bn}_R1_cutadapt.fastq.gz !{params.cutadapt_dir}/!{bn}/
         cp !{bn}_R2_cutadapt.fastq.gz !{params.cutadapt_dir}/!{bn}/

         """

}



// mkdir -p /shared/projects/mbnl_dct/mouse_seq/cut_adapted/  \
// cutadapt --cores 8 -B AGATCGGAAGAG --minimum-length 75 --report  minimal  \
// --o /shared/projects/mbnl_dct/mouse_seq/cut_adapted/1051_R1_cutadapt.fastq.gz \
// -p /shared/projects/mbnl_dct/mouse_seq/cut_adapted/1051_R2_cutadapt.fastq.gz \
// -Z \
// /shared/projects/mbnl_dct/mouse_seq/fastq/1051_R1.fastq \
// /shared/projects/mbnl_dct/mouse_seq/fastq/1051_R2.fastq


process STAR_Alignment {

      cpus "${params.cpus}"
      memory "40G"
      module "star/2.7.5a:perl/5.26.2"

      input:
      tuple val(bn),path(f1),path(f2) from cutadChannel
      // tuple val(bn),path(f1),path(f2) from inputChannel

      output:
      tuple path("${bn}/${bn}_Aligned.sortedByCoord.out.bam"),val(bn) into align1Channel
      // tuple path("${bn}/${bn}_Aligned.sortedByCoord.out.bam"),val(bn) into align2Channel


      shell:
      """
      mkdir -p !{params.outDir}/!{bn};
      echo "SampleID=", !{bn};
      STAR --runMode  alignReads \
           --runThreadN !{params.cpus} \
           --genomeDir !{params.refGenome}  \
           --outFilterMultimapNmax 20 \
           --alignSJoverhangMin 8 \
           --alignSJDBoverhangMin 1 \
           --outFilterMismatchNmax 999 \
           --outFilterMismatchNoverReadLmax 0.04 \
           --alignIntronMin 20 \
           --alignIntronMax 1000000 \
           --alignMatesGapMax 1000000 \
           --readFilesIn !{f1}  !{f2}  \
           --outFileNamePrefix !{bn}/!{bn}_  \
           --outSAMtype BAM SortedByCoordinate \
           --readFilesCommand zcat;

       cp  -r !{bn}/*  !{params.outDir}/!{bn}/ ;

       """
}

// --outFilterType BySJout : reduces the number of ”spurious” junctions
// --outFilterMultimapNmax 20 : max  number  of  multiple  alignments  allowed  for  a  read:  if  exceeded,  the  read  is  considered unmapped
// --alignSJoverhangMin 8 : minimum overhang for unannotated junctions
// --alignSJDBoverhangMin 1 : minimum overhang for annotated junctions
// --outFilterMismatchNmax 999 : maximum number of mismatches per pair, large number switches off this filter
// --outFilterMismatchNoverReadLmax 0.04 : max number of mismatches per pair relative to read length:  for 2x100b,
//                         max number of mis-matches is 0.04*200=8 for the paired read
// --alignIntronMin 20  : minimum intron length
// --alignIntronMax 1000000 : maximum intron length
// --alignMatesGapMax 1000000 : maximum genomic distance between mates

process HTSeq_counts_on_Gene {

    	memory "40G"
    	module "samtools/1.10:htseq/0.12.4"


    	input:
    	tuple path(bam),val(bn) from align1Channel

    	output:
    	path("${bn}_counts_gene.txt") into htseqChannel

    	shell:
        """
        mkdir -p !{params.countsDir};

        samtools index !{bam} ;

        if [ -e !{bam}.bai ] ; then echo "Index created"; fi

        htseq-count  -r pos  -s reverse  --type gene --idattr gene_id \
        --mode union --minaqual 10  \
        --counts_output  !{bn}_counts_gene.txt \
        !{bam}  !{params.ref_Gtf} ;

        cp  !{bn}_counts_gene.txt  !{params.countsDir}/

        """
}

/*
process DEXSeq_counts_Exons {

         memory "40G"
         module "samtools/1.10:python/3.7:htseq/0.11.2"

         input:
         tuple path(bam),val(bn) from align2Channel

         output:
         path("${bn}_dexseq_counts.txt") into dexseqChannel

         shell:
         """
         mkdir -p !{params.dexseqDir}/
         samtools index !{bam} ;
         if [ -e !{bam}.bai ] ; then echo "Index created"; fi

         /shared/home/mkondili/Tools_Packages/dexseq_count.py \
         --format bam --paired yes --stranded reverse \
         --minaqual 10 --order pos  \
         !{params.dexseq_gtf} !{bam}  !{bn}_dexseq_counts.txt

         cp !{bn}_dexseq_counts.txt  !{params.dexseqDir}/
         ## Usage: python dexseq_count.py [options] <flattened_gff_file> <alignment_file> <output_file>

         """
         //  /shared/home/mkondili/Tools_Packages/dexseq_count.py \
         //  --format bam --paired yes --stranded reverse --minaqual 10 --order pos
         //  /shared/home/mkondili/genomes/Annotation/mm10/Galaxy_produced_DEXSeq_annot_mm10.gtf \
         // /shared/projects/mbnl_dct/mouse_seq/STAR_Aligned/1051/1051_Aligned.sortedByCoord.out.bam \
         // /shared/projects/mbnl_dct/mouse_seq/STAR_DEXSeq_counts/1051_dexseq_counts.txt
}

*/

// in python:
// gff_file = "/shared/home/mkondili/genomes/Annotation/mm10/Galaxy_produced_DEXSeq_annot_mm10.gtf"
// sam_file = "/shared/projects/mbnl_dct/mouse_seq/STAR_Aligned/1051/1051_Aligned.sortedByCoord.out.bam"
// out_file = "/shared/projects/mbnl_dct/mouse_seq/STAR_DEXSeq_counts/1051_dexseq_counts.txt"
// stranded = True
// reverse = True
// is_PE = True
// minaqual = 10
// order = "pos"


// --outReadsUnmapped = Fastx:
// will output unmapped and partially mapped
//  (i.e.  mapped only onemate of a paired end read) reads into separate file(s)
// Unmapped.out.mate1(2), formatted the sameway as input read files
// (i.e.  FASTQ or FASTA).
// Appended to the read name line are tag to indicatemapping status of the read mates:
// 00:  mates were not mapped;
// 10:  1st mate mapped, 2nd unmapped
// 01:  1st unmapped, 2nd mapped


// -- ENCODE options ----
// An example of ENCODE standard options for long RNA-seq pipeline is given below:
// --outFilterType BySJout : reduces the number of ”spurious” junctions
// --outFilterMultimapNmax 20 : max  number  of  multiple  alignments  allowed  for  a  read:  if  exceeded,  the  read  is  considered unmapped
// --alignSJoverhangMin 8 : minimum overhang for unannotated junctions
// --alignSJDBoverhangMin 1 : minimum overhang for annotated junctions
// --outFilterMismatchNmax 999 : maximum number of mismatches per pair, large number switches off this filter
// --outFilterMismatchNoverReadLmax 0.04 : max number of mismatches per pair relative to read length:  for 2x100b,
//                         max number of mis-matches is 0.04*200=8 for the paired read
// --alignIntronMin 20  : minimum intron length
// --alignIntronMax 1000000 : maximum intron length
// --alignMatesGapMax 1000000 : maximum genomic distance between mates
