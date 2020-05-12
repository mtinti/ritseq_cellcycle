#!/bin/bash 
#$ -l local_free=50G
#$ -l h_vmem=4G
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 8

set -e
set -u 
set -o pipefail

experiment='CM03_H_GG2M_S7'
genome='tb927_3'
base_fastq='S7_'


function timer()
{
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')

        if [[ -z "$stime" ]]; then stime=$etime; fi

        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        printf '%d:%02d:%02d' $dh $dm $ds
    fi
}


echo 'copy data in $TMPDIR'
##variables

tmr=$(timer)

mkdir -p $TMPDIR'/genomes/'$genome'/' && cp -avr  'genomes/'$genome'/' $TMPDIR'/genomes/' 
mkdir -p $TMPDIR'/'$experiment'/data/' && cp -avr  $experiment'/data/' $TMPDIR'/'$experiment'/'
mkdir -p $TMPDIR'/FastQ_Screen_Genomes/' && cp -avr  'FastQ_Screen_Genomes/' $TMPDIR
path_genome_index=$TMPDIR'/genomes/'$genome'/'$genome
path_trascriptome_index=$TMPDIR'/trascriptome/'$genome'/'$genome
printf 'Elapsed time: %s\n' $(timer $tmr) 

echo 'folder $TMPDIR'
ls -l $TMPDIR

path_out=$TMPDIR'/'$experiment'/'
mkdir $path_out'fastqc'
#mkdir $path_out'fastp'
mkdir $path_out'fastqs'
mkdir $path_out'qual_out'
mkdir $path_out'qual_out/'$experiment
ls $path_out -l
##qc of fastq 
echo 'run 0.1' 


fastq_1=$path_out'data/'$base_fastq'1.fastq.gz'
fastq_2=$path_out'data/'$base_fastq'2.fastq.gz'


echo 'find size of fastq file'
inFileLength=$(echo $(zcat $fastq_1 | wc -l)/4|bc)
echo 'inFileLength all: '$inFileLength


echo 'run fastqc' 
#fastqc  $fastq_1 -o $path_out'fastqc'
#fastqc  $fastq_2 -o $path_out'fastqc'

#remove problematic reads
#echo 'run fastp' 
#out_1=$path_out'data/'$base_fastq'f1.fastq.gz'
#out_2=$path_out'data/'$base_fastq'f2.fastq.gz'
#fastp -i $fastq_1 -I $fastq_2 -o $out_1 -O $out_2 -h $path_out'fastp/fastp_'$base_fastq'.html' -j $path_out'fastp/'$base_fastq'_fastp.json'

echo 'run bowtie2'
(bowtie2 --very-sensitive-local -p 8 -x $path_genome_index -1 $fastq_1 -2 $fastq_2) \
2>$path_out'qual_out/'$experiment'/'$experiment'.log' | samtools view -bSu | \
samtools sort -@ 8 -o $path_out$base_fastq'sorted.bam'

samtools index $path_out$base_fastq'sorted.bam'

#rm $out_1
#rm $out_2


unset DISPLAY
echo 'run qualimap bamqc' 
qualimap bamqc --java-mem-size=4G -bam $path_out$base_fastq'sorted.bam' \
-outdir $path_out'qual_out/'$experiment \
-outfile $experiment'.bamqc.html' -outformat 'HTML'

echo 'run qualimap rnaseq' 
qualimap rnaseq --java-mem-size=4G -bam $path_out$base_fastq'sorted.bam'  \
-gtf $path_genome_index'.gtf' \
-outdir $path_out'qual_out/'$experiment \
-pe -outfile $experiment'.rnaseq.html' -outformat 'HTML'


## sorting because picard MarkDuplicates
## needs name sorted to filter out secondary aligment
##maybe not
#echo 'run samsort by name' 
#samtools sort -o $path_out$base_fastq'sorted.bam' -n -@ 8 $path_out$base_fastq'sorted.bam'
#rm $path_out$base_fastq'sorted.bam'
#rm $path_out$base_fastq'sorted.bam.bai'

bedtools genomecov -ibam $path_out$base_fastq'sorted.bam' \
-bg -pc > $path_out$base_fastq'sorted_pc_nodepup.bed'
gzip $path_out$base_fastq'sorted_pc_nodepup.bed'


bedtools genomecov -ibam $path_out$base_fastq'sorted.bam' \
-bg > $path_out$base_fastq'sorted_nodepup.bed'
gzip $path_out$base_fastq'sorted_nodepup.bed'




echo 'run picard MarkDuplicates' 
picard MarkDuplicates \
I=$path_out$base_fastq'sorted.bam' \
O=$path_out$base_fastq'sorted_dedup.bam' \
M=$path_out'qual_out/'$experiment'/'$experiment'.metrics.txt' \
REMOVE_DUPLICATES=true \



## sorting becouse the extract_barcode scrips
## needs reference sorted outputs and index of bam file present
echo 'run samsort by refrence' 
samtools sort -@ 8 -o $path_out$base_fastq'sorted_dedup.bam' $path_out$base_fastq'sorted_dedup.bam'
samtools index $path_out$base_fastq'sorted_dedup.bam' 

#rm $path_out$base_fastq'sorted_name.bam'
#rm $path_out$base_fastq'sorted_name_dedup.bam'




echo 'extract unmapped' 
mkdir $path_out'unmap/'
samtools view -bu -f 12 -F 256 $path_out$base_fastq'sorted_dedup.bam' | samtools sort -n -@ 8 -o $path_out$base_fastq'unmap.bam'
bamToFastq -i $path_out$base_fastq'unmap.bam' -fq $path_out'unmap/'$base_fastq'unmap_1.fastq' -fq2 $path_out'unmap/'$base_fastq'unmap_2.fastq'

echo 'run fastq_screen' 
sed -i "s@tmp_folder@$TMPDIR@g" $TMPDIR'/FastQ_Screen_Genomes/fastq_screen2.conf'

fastq_screen --nohits $path_out'unmap/'$base_fastq'unmap_1.fastq' $path_out'unmap/'$base_fastq'unmap_2.fastq' \
--conf $TMPDIR'/FastQ_Screen_Genomes/fastq_screen2.conf' --outdir $path_out'qual_out/'$experiment --aligner bowtie2

rm $path_out$base_fastq'unmap.bam'
rm $path_out'unmap/'$base_fastq'unmap_2.fastq'
gzip $path_out'unmap/'$base_fastq'unmap_1.fastq'


echo 'extract properly paired'
samtools view -b -f 2 -@ 8 -o $path_out$base_fastq'sorted_dedup_.bam' $path_out$base_fastq'sorted_dedup.bam'
rm $path_out$base_fastq'sorted_dedup.bam'
mv $path_out$base_fastq'sorted_dedup_.bam' $path_out$base_fastq'sorted_dedup.bam'

#samtools sort -@ 8 -o $path_out$base_fastq'sorted_dedup.bam' $path_out$base_fastq'sorted_dedup.bam'
#samtools index $path_out$base_fastq'sorted_dedup.bam' 


echo 'run extract_barcode' 
python mylib/extract_barcodes_def2.py \
$path_out$base_fastq'sorted_dedup.bam' GTGAGGCCTCGCGA TCGCGAGGCCTCAC 1 True


bedtools genomecov -ibam $path_out$base_fastq'sorted_dedup_F_plus_R.bam' \
-bg > $path_out$base_fastq'sorted_dedup_F_plus_R.bed'
bedtools genomecov -ibam $path_out$base_fastq'sorted_dedup_WO.bam' \
-bg > $path_out$base_fastq'sorted_dedup_WO.bed'
bedtools genomecov -ibam $path_out$base_fastq'sorted_dedup.bam' \
-bg > $path_out$base_fastq'sorted_dedup.bed'

gzip $path_out$base_fastq'sorted_dedup_F_plus_R.bed'
gzip $path_out$base_fastq'sorted_dedup_WO.bed'
gzip $path_out$base_fastq'sorted_dedup.bed'


#bedtools genomecov -ibam $path_out$base_fastq'sorted_dedup_F_plus_R.bam' \
#-bg -pc > $path_out$base_fastq'sorted_dedup_F_plus_Rpc.bed'
#bedtools genomecov -ibam $path_out$base_fastq'sorted_dedup_WO.bam' \
#-bg -pc > $path_out$base_fastq'sorted_dedup_WOpc.bed'
#bedtools genomecov -ibam $path_out$base_fastq'sorted_dedup.bam' \
#-bg -pc > $path_out$base_fastq'sorted_dedup_pc.bed'

#gzip $path_out$base_fastq'sorted_dedup_F_plus_Rpc.bed'
#gzip $path_out$base_fastq'sorted_dedup_WOpc.bed'
#gzip $path_out$base_fastq'sorted_dedup_pc.bed'




echo 'run featureCounts' 
#echo 'run 11'
#B BothEndsMapped
#C exclude chimeric (reads map on different chr)
#M count multi mapping
featureCounts -p -B -C -M -T 8 -t transcript -g gene_id -a $path_genome_index'.gtf' \
-o $path_out'counts.txt' \
$path_out$base_fastq'sorted_dedup.bam' \
$path_out$base_fastq'sorted_dedup_F_plus_R.bam' \
$path_out$base_fastq'sorted_dedup_WO.bam' \
$path_out$base_fastq'sorted.bam'

cp $path_out'counts.txt.summary' $path_out'qual_out/'$experiment'/'

echo 'copy results' 
rm -r $TMPDIR'/'$experiment'/data'
mkdir -p $experiment'/res' && cp -avr $TMPDIR'/'$experiment $experiment'/res'



