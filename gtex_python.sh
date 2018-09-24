#!/bin/bash

#PBS -N gtex_p
#PBS -l nodes=1:ppn=12
#PBS -l gres=localhd:20
#PBS -l mem=72gb
#PBS -l vmem=72gb
#PBS -l walltime=30:00:00
#PBS -e /home/acelik/lp/acelik_files/MH/rnaseq/gtex_p.err
#PBS -o /home/acelik/lp/acelik_files/MH/rnaseq/gtex_p.out


###############
## this script takes 3 arguments using the -v option of qsub
## they are namesd variables
## outdir= where the outputs should go
## readdir= where the reads are
## sample_id= what the samples are names
## folders will be created using the outdir and the sample_id variables under outdir
## if you are submitting things in a loop make sure the names do no clash
##############


module load rsem/1.2.22
module load star/2.5.4b
module load picard-tools/2.18.0
module load java/1.8.0_161
module load trimmomatic/0.32
module load python/3.5.6
module load samtools
module load bowtie 


adapters="/hpf/tools/centos6/mugqic-pipelines/source/resource/software/trimmomatic/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa"
star_index="/home/acelik/lp/acelik_files/gtex/reference/star_reference"
rsem_reference="/home/acelik/lp/acelik_files/gtex/reference/rsem_reference/rsem_reference"
genome_fasta="/home/acelik/lp/acelik_files/gtex/reference/Homo_sapiens_assembly19.fasta"
genes_gtf="/home/acelik/lp/acelik_files/gtex/reference/gencode.v19.annotation.patched_contigs.gtf"

picard="/hpf/tools/centos6/picard-tools/2.18.0/picard.jar"
rnaseqqc="/hpf/tools/centos6/rna-seqc/1.1.8/RNA-SeQC_v1.1.8.jar"
trimmomatic="/hpf/tools/centos6/trimmomatic/0.32/trimmomatic-0.32.jar"

pyscripts=/home/acelik/lp/acelik_files/gtex/gtex-pipeline/rnaseq/src


##############################################
################# TODO #######################
##############################################

## in a later version need another python script to 
## parse steps
## also, similar to genpiples needs to take sample name, libtype, reads1 and reads2 
## if readtype is single end. 
## need another script that reads a files generates qsub files and submits them if specified
## all files in the readset files needs to be checked before scritp generation. 

fastq1=$readsdir/${sample_id}_R1.fastq.gz
fastq2=$readsdir/${sample_id}_R2.fastq.gz

echo "read files are $fastq1 $fastq2"
echo "starting adapter trimming"

mkdir -p $outdir/trim/$sample_id

java -jar $trimmomatic PE \
  -threads 12 \
  -phred33 \
  $fastq1 \
  $fastq2 \
  $outdir/trim/${sample_id}/${sample_id}_R1_pair_trimmed.fastq.gz \
  $outdir/trim/${sample_id}/${sample_id}_R1_single_trimmed.fastq.gz \
  $outdir/trim/${sample_id}/${sample_id}_R2_pair_trimmed.fastq.gz \
  $outdir/trim/${sample_id}/${sample_id}_R2_single_trimmed.fastq.gz \
  ILLUMINACLIP:${adapters}:2:30:15 \
  TRAILING:30 \
  MINLEN:32 2> $outdir/trim/${sample_id}_trim.log

CODE=$?

if [ $CODE -eq 0 ]
then
echo "trimmonatic done"
else
echo "there is some error in trimmomatic see $outdir/trim/${sample_id}_trim.log"
fi

echo "starting star alignment"
mkdir -p $outdir/star/$sample_id

python3 $pyscripts/run_STAR.py --index $star_index \
   --fastq $outdir/trim/${sample_id}/${sample_id}_R1_pair_trimmed.fastq.gz,$outdir/trim/${sample_id}/${sample_id}_R2_pair_trimmed.fastq.gz \
   --prefix $sample_id \
   -o $outdir/star/$sample_id \
   --annotation_gtf $genes_gtf \
   -t 12 2>$outdir/star/${sample_id}.log

CODE=$?

if [ $CODE -eq 0 ]
then
echo "STAR alignment done"
else
echo "there is some error in alignment see $outdir/star/${sample_id}"
fi


echo "starting Picard MarkDuplicates"
mkdir -p $outdir/picard/$sample_id

java -jar $picard MarkDuplicates \
   INPUT=$outdir/star/$sample_id/${sample_id}.Aligned.sortedByCoord.out.bam \
   OUTPUT=$outdir/picard/$sample_id/${sample_id}.Aligned.sortedByCoord.out.md.bam \
   METRICS_FILE=$outdir/picard/${sample_id}.metrics.txt 2> $outdir/picard/${sample_id}.log

samtools index $outdir/picard/$sample_id/${sample_id}.Aligned.sortedByCoord.out.md.bam

#python3 $pyscripts/run_MarkDuplicates.py \
#        $outdir/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam \
#        ${sample_id}.Aligned.sortedByCoord.out.patched.md \
#        --output_dir $outdir/picard/$sample_id \
#        --jar $picard 2> $outdir/picard/${sample_id}.log

if [ $CODE -eq 0 ]
then
echo "Mark Duplicates done"
else
echo "there is some error in markduplicates see $outdir/picard/${sample_id}"
fi

echo "starting RSEM"
mkdir -p $outdir/rsem/$sample_id

rsem-calculate-expression -p 12 --bowtie2 --forward-prob 0.0 \
   --no-bam-output --calc-ci --estimate-rspd \
   --paired-end $outdir/trim/${sample_id}/${sample_id}_R1_pair_trimmed.fastq.gz $outdir/trim/${sample_id}/${sample_id}_R2_pair_trimmed.fastq.gz \
   $rsem_reference $outdir/rsem/$sample_id 2>$outdir/rsem/${sample_id}.rsem.log

#convert-sam-for-rsem $outdir/star/${sample_id}/${sample_id}.Aligned.toTranscriptome.out.bam \
#   $outdir/star/${sample_id}/${sample_id}.converted

#python3 $pyscripts/run_RSEM.py \
#  -o $outdir/rsem/$sample_id \
#  --rsem_ref_dir $rsem_reference \
#  --input_file $outdir/star/$sample_id/${sample_id}.converted \
#  --prefix $sample_id

if [ $CODE -eq 0 ]
then
echo "RSEM alignment done"
else
echo "there is some error in alignment see $outdir/rsem/${sample_id}"
fi

echo "starting QC"
mkdir -p $outdir/qc/$sample_id

echo "switching java version"
module unload java/1.8.0_161
module load rna-seqc/1.1.8

java -jar $rnaseqqc \
-n 1000 -s $sample_id,$outdir/picard/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.md.bam,${sample_id} \
-t $genes_gtf \
-r $genome_fasta \
-o $outdir/qc/${sample_id} \
-noDoc -strictMode 2>$outdir/qc/${sample_id}.qc.log
 
#python3 $pyscripts/run_rnaseqc.py \
#    $outdir/picard/$sample_id/${sample_id}.Aligned.sortedByCoord.out.md.bam \
#    ${genes_gtf} \
#    ${genome_fasta} \
#    ${sample_id} \
#    --jar $rnaseqqc
#    --output_dir $outdir/qc/$sample_id 2> $outdir/qc/${sample_id}.log

if [ $CODE -eq 0 ]
then
echo "QC done"
else
echo "there is some error in QC see $outdir/qc/${sample_id}"
fi


