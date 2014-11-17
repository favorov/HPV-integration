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

$sai_1_file = $out_prefix.".01.sai";
$sai_2_file = $out_prefix.".02.sai";

system "bwa aln -l 32 -k 2 -t 4 $input_fasta_file -1 $input_fq_01_file > $sai_1_file";
system "bwa aln -l 32 -k 2 -t 4 $input_fasta_file -2 $input_fq_02_file > $sai_2_file";
system "bwa sampe $input_fasta_file $sai_1_file $sai_2_file $input_fq_01_file $input_fq_02_file| samtools view -Shb - > $output_bam_unsorted_file" or die;
system "samtools sort $output_bam_unsorted_file $out_prefix" or die;
system "samtools index $output_bam_file" or die;
system "rm -f $sai_1_file $sai_2_file $output_bam_unsorted_file" or die;
