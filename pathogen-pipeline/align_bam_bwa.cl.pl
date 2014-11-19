#!/usr/bin/perl -w
#was e.map_se_bam2genome.cl.pl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;

my ($input_bam_file, $input_fasta_file, , $out_prefix, $sai_1_file, $output_bam_file);
GetOptions ("fasta=s" => \$input_fasta_file,
			"bam=s" => \$input_bam_file,
			"out=s" => \$out_prefix,
            );
			
$output_bam_file = $out_prefix.".bam";

my $output_bam_unsorted_file_prefix = $out_prefix.".unsorted";
my $output_bam_unsorted_file = $output_bam_unsorted_file_prefix.".bam";

my $output_bam_file = $out_prefix.".bam";

if (-f $output_bam_file and -s $output_bam_file and -f "$output_bam_file.bai" and -s "$output_bam_file.bai")
{
	print "$output_bam_file $output_bam_file.bai exist, we do not run fasta alignment.\nIf it is an error, remove them manually before start.\n";
	exit 0;
}

$sai_1_file = $out_prefix.".01.sai";

#bwa aln -l 32 -k 2 -t 4 $input_fasta_file -b $input_bam_file > $sai_1_file
#bwa samse $input_fasta_file $sai_1_file $input_bam_file | samtools view  -Shb -F4 - > $output_bam_unsorted_file
#samtools sort $output_bam_unsorted_file $out_prefix
#samtools index $output_bam_file

if (not -f $sai_1_file or not -s $sai_1_file)
{
	my $fasta1_run="bwa aln -l 32 -k 2 -t 4 $input_fasta_file -b $input_bam_file > $sai_1_file";
	system ($fasta1_run) == 0 
		or die "system $fasta1_run died: $?";
}


if (not -f $output_bam_unsorted_file or not -s $output_bam_unsorted_file)
{
	my $fasta_combine_run="bwa samse $input_fasta_file $sai_1_file $input_bam_file | samtools view  -Shb -F4 - > $output_bam_unsorted_file";
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

system "rm -f $sai_1_file $output_bam_unsorted_file";
