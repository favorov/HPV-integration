#!/usr/bin/perl -w

use strict;
use warnings;
use feature 'say';

my $result = $ARGV[0];
my $summary_result = $result.".s";
my $output_sort_file = $result.".so";
open my ($unput_file), '<', $result or die "Cannot open input file: $!";
open my ($output_file), '>', $summary_result or die "Cannot open file for writing: $!";
while (<$unput_file>) {
        print $output_file "$_" if (/^>gi/);
}
close $unput_file;
close $output_file;
system "sort -rn -t \$'\\t' -k3,3 $summary_result > $output_sort_file";
#system "sort -r -t'\t' -k3,3 $output_file > $output_sort_file";
