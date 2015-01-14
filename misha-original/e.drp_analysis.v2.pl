#!/usr/bin/perl -w
use strict;
use warnings;
use feature 'say';
use Getopt::Long;

### 
### 

my ($input_hg19_sam_file, $output_prefix, @raw_data, @drp_array, $n, $threshold);
GetOptions ("sam=s" => \$input_hg19_sam_file,
			"out=s" => \$output_prefix,
            "thresh=i" => \$threshold,);

my $output_drp_file = $output_prefix.".drp";


open my ($input_sam), '<', $input_hg19_sam_file or die "Cannot open fam file: $!";
open my ($output_drp), '>>', $output_drp_file or die "Cannot open file for writing: $!";

@raw_data = <$input_sam>;
my $size_raw_data = @raw_data;
close $input_sam;
my $i;
$n = 0;
my $real_clusters = 0;
my $count = 0;
for ($i=0; $i<$size_raw_data; $i++) {
	
	if ($i>0) {
		my $last_full_line = $raw_data[$i-1];
		chomp($last_full_line);
		my @last_sam_array = split /\t/, $last_full_line;
		my $current_full_line = $raw_data[$i];
		chomp($current_full_line);
		my @current_sam_array = split /\t/, $current_full_line;
		my $difference = abs($current_sam_array[3] - $last_sam_array[3]);
	
		if (($current_sam_array[2] eq $last_sam_array[2]) and ($difference < 100)) {
			if ($n==0) {
				push(@drp_array, $last_full_line);
			}
			push(@drp_array, $current_full_line);
			$n=1;
		}
		else {
			if (@drp_array) {
				$count++;
				my $size_drp_array = @drp_array;
				if ($size_drp_array >= $threshold) {
                    $real_clusters++;
					my @first_line_drp_array = split /\t/, $drp_array[0];
					my @last_line_drp_array = split /\t/, $drp_array[-1];
					say $output_drp "### CLUSTER number $count; # of reads: $size_drp_array ###\n### CLUSTER: $first_line_drp_array[2]:$first_line_drp_array[3]-$last_line_drp_array[3] ###END";
					foreach(@drp_array) {
						say $output_drp "$_";
					}
					say $output_drp "\n";
				}
				@drp_array = ();
				$n=0;
			}
		}
	}
}
#print $output_drp "$count";
say "Number of drp clusters: $real_clusters";
say "Discordant reads threshold: $threshold";
close $output_drp;	
	
