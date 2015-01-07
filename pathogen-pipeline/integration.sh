#!/bin/bash
#here, we suppose that test_for_exist script is already competed,
#and we have the bam file with the name 
#here, there are two parameters
# $1 is the name of the fasta without _1 and _2. It is kind of id of the experinemt and the previouos script generated ${name}_aligned.bam
# in the instructions from Misha, the ${name}_aligned.bam is called original.bam
# $2 is full name (with path and extension) of the reference viral FASTA file, we suppose it is indexed by bwa in its folder (bwa index -a is myvirus.fasta)


if [ ! -n "$1" ]
then
 echo Parameter 1 must be name prefix to be expanede to original bam file by adding '_aligned.bam'  
 exit
fi 

name=$1

original_bam=${name}-aligned.bam

if [ ! -f ${original_bam} ]
then
 echo file ${original_bam} does not ezist  
 exit
fi 

if [ ! -n "$2" ]
then
 echo Parameter 2 must be name the full name \(with path and extension\) reference viral genome  
 exit
fi 

virus=$2

if [ ! -f ${virus} ]
then
 echo file ${virus} does not ezist  
 exit
fi 

export pathogen=~/pathogen_pipeline_2014

unmapped_mate_of_mapped=${name}_umom.bam

mapped_mate_of_unmapped=${name}_moum.bam

if (! ( [ -f $unmapped_mate_of_mapped ] && [ -s $unmapped_mate_of_mapped ] )) 
then
	echo filter out unmapped mate of mapped
	samtools view -bh -f 4 -F 264 $original_bam > $unmapped_mate_of_mapped 
	#-f include all
	#-F excluse (any?)
	# 264=octal 108
	#0x4 segment unmapped
	#0x8 next segment in the template unmapped
	#0x100 secondary alignment
	echo done...
fi


if (! ( [ -f $mapped_mate_of_unmapped ] && [ -s $mapped_mate_of_ummapped ] )) 
then
	echo filter out mapped mate of ummapped
	samtools view -bh -f 8 -F 260 $original_bam > $mapped_mate_of_unmapped 
	#-f include all
	#-F excluse (any?)
	# 264=octal 108
	#0x4 segment unmapped
	#0x8 next segment in the template unmapped
	#0x100 secondary alignment
	echo done...
fi

#Map unmapped ends to HPV16 genome
#the result is ${name}.s.bam
if (! ( [ -f ${name}.s.bam ] && [ -s ${name}.s.bam ] )) 
then
	echo map ummapped to virus
	${pathogen}/e.map_se_bam2genome.cl.pl --fasta ${virus} --bam ${unmapped_mate_of_mapped} --out ${name}
	echo done...
fi



