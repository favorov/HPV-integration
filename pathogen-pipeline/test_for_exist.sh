#!/bin/bash
name=$1

if [ ! -n "$1" ]
then
 echo Parameter 1 must be name prefix to be expanede to fasta names by add _1.fq.qz and _2.fq.gz
 exit
fi 

fasta1=${name}_1.fq.gz
fasta2=${name}_2.fq.gz

if [ ! -f ${fasta1} ]
then
 echo file ${fasta1} does not ezist  
 exit
fi 

if [ ! -f ${fasta2} ]
then
 echo file ${fasta2} does not ezist  
 exit
fi

export pathogen=~/pathogen_pipeline_2014
hg19=${pathogen}/genomes/hg19_bwa_indexed/Homo_sapiens_assembly19.fa

${pathogen}/align_pair_fastas_bwa.pl --fasta $hg19 --fastq-01 ${fasta1} --fastq-02 ${fasta2} --out ${name}-aligned || exit
#aligned
echo filtering out the mapped ...
samtools view -bh -f 4 ${name}-aligned.bam > ${name}-unmapped.bam || exit
echo done...
#mapped filtered out
${pathogen}/pipeline.bwa.v8.viral.v012913_bam.pl ${name}-unmapped.bam || exit
#unmapped aligned on viral genome, look for results in *ex.so 

mv ${name}-unmapped-f.c.sam.sort.bam.graph.gt4c.ex.so ${name}-pathogen-report.txt
rm ${name}-unmapped-f.c.sam.sort.bam.bai
mkdir -p presence-test-junk
mv ${name}-unmapped-s* presence-test-junk
mv ${name}-unmapped-m* presence-test-junk
mv ${name}-unmapped-f* presence-test-junk
mv ${name}-unmapped.xa* presence-test-junk
