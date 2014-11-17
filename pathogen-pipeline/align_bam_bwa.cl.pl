#!/usr/bin/perl -w
#was e.map_se_bam2genome.cl.pl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;

my ($input_bam_file, $input_fasta_file, $prefix, $sai_1_file, $output_bam_file);
GetOptions ("fasta=s" => \$input_fasta_file,
			"bam=s" => \$input_bam_file,
			"out=s" => \$prefix,
            );
			
$output_bam_file = $prefix.".bam";

my $output_bam_unsorted_file_prefix = $prefix.".unsorted";
my $output_bam_unsorted_file = $output_bam_unsorted_file_prefix.".bam";

my $output_bam_for_sort_file = $prefix.".s";
my $output_bam_sort_file = $prefix.".s.bam";

$sai_1_file = $prefix.".01.sai";

system "bwa aln -l 32 -k 2 -t 4 $input_fasta_file -b $input_bam_file > $sai_1_file";
system "bwa samse $input_fasta_file $sai_1_file $input_bam_file | samtools view  -Shb -F4 - > $output_bam_unsorted_file" or die;
system "samtools sort $output_bam_unsorted_file $prefix" or die;
system "samtools index $output_bam_file" or die;
system "rm -f $sai_1_file $output_bam_unsorted_file" oe die;
