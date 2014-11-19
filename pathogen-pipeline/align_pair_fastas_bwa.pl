#!/usr/bin/perl -w
#was e.map_pe_fastq_gz2genome.cl.pl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;


my ($input_fq_01_file, $input_fq_02_file, $input_fasta_file, $sai_1_file, $sai_2_file, $output_bam_file, $out_prefix);
GetOptions ("fasta=s" => \$input_fasta_file,
			"fastq-01=s" => \$input_fq_01_file,
			"fastq-02=s" => \$input_fq_02_file,
            "out=s" => \$out_prefix,
            );
			
$output_bam_file = $out_prefix.".bam";
my $output_bam_unsorted_file_prefix = $out_prefix.".unsorted";
my $output_bam_unsorted_file = $output_bam_unsorted_file_prefix.".bam";

if (-f $output_bam_file and -s $output_bam_file and -f "$output_bam_file.bai" and -s "$output_bam_file.bai")
{
	print "$output_bam_file $output_bam_file.bai exist, we do not run fasta alignment.\nIf it is an error, remove them manually before start.\n";
	exit 0;
}

$sai_1_file = $out_prefix.".01.sai";
$sai_2_file = $out_prefix.".02.sai";

if (not -f $sai_1_file or not -s $sai_1_file)
{
	my $fasta1_run="bwa aln -l 32 -k 2 -t 4 $input_fasta_file -1 $input_fq_01_file > $sai_1_file";
	system ($fasta1_run) == 0 
		or die "system $fasta1_run died: $?";
}

if (not -f $sai_2_file or not -s $sai_2_file)
{
	my $fasta2_run="bwa aln -l 32 -k 2 -t 4 $input_fasta_file -2 $input_fq_02_file > $sai_2_file";
	system ($fasta2_run) == 0 
		or die "system $fasta2_run died: $?";
}

if (not -f $output_bam_unsorted_file or not -s $output_bam_unsorted_file)
{
	my $fasta_combine_run="bwa sampe $input_fasta_file $sai_1_file $sai_2_file $input_fq_01_file $input_fq_02_file | samtools view -Shb - > $output_bam_unsorted_file";
	system ($fasta_combine_run) == 0 
		or die "system $fasta_combine_run died: $?";
}

if (not -f $output_bam_file or not -s $output_bam_file)
{
	print "samtools sort\n";
	my $samtools_sort_run="samtools sort $output_bam_unsorted_file $out_prefix";
	system ($samtools_sort_run) == 0
		or die "system $samtools_sort_run died: $?";
}

if (not -f "$output_bam_file.bai" or not -s "$output_bam_file.bai")
{
	print "samtools index\n";
	my $samtools_index_run="samtools index $output_bam_file";
	system ($samtools_index_run) == 0 
		or die "system $samtools_index_run died : $?";
}

system "rm -f $sai_1_file $sai_2_file $output_bam_unsorted_file";
