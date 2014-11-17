#!/usr/bin/perl -w
use strict;
use warnings;
use feature 'say';
use Readonly;

=head1
    01/08/14
        1) Because /files doesn't work anymore with LSF I copied the script, fasta, and index files (and bwa) to /groups 

    02/11/13
        1) I introduce check procedure for the muliple hits.

    01/28/13
        1) Use updated viral FASTA file with added 24 Human coxsackieviruses
        /files/Genetics/seidman/michael/genomes/viruses/current_work_viral_version/viral.genomes.mp205.v012813.fasta

    01/16/13
        1) Use 'fastq.gz' input file instead of 'fq'.
        2) Use the latest BWA version bwa-0.6.2, /files/Genetics/seidman/michael/save/bwa.
        3) 10 same nucleotides in a row filtering is excluded.
# 03/29/11
# BWA reports multiple hits in XA tags. Take them into analysis but only if
# they relate to another genomes not the same genome but another position
# HOWEVER! IGV will show only one randomly chosen by BWA hit 

# 04/07/11
# Now .so file will include genomes not only with 4 or more covered regions but also results   with 200bp or more coverage (not to miss 1-3 large contiguous regions)
# I chose 200bp though I'm not sure that it is the best threshold

# It would be great to know the length of covered regions shared with other genomes (in the case of multiple hits. The only problem is to distinguish mh from different genomes and from the different positions in the same genome).

# Use unmasked genomes! There is no repeats noise on the level of more than 2 covered regions.

# 04/12/11
# Reads that could be equally good mapped on to genomes of different species (not strains) would be excluded
# from the coverage analysis though retained in the brief and full summaries and in multgenomes file.
# I don't want to look on the genomes with reads mapped only on to conservative ribosomal RNA regions
# that are common for numerous species

# 04/14/11
# Use new version of bacterial genome reference: v041311
# complete genomes + complete sequences (excluding all E.coli except one strain and excluding repeated genomes from 4 species)

#- coverage threshold 0 (instead of 4)

=cut

## viral reference
#Readonly my $fastaFile => "/files/Genetics/seidman/michael/tasha_TCGA_unmapped/bwa_index/bacterial_genomes/RM_907.WedMar91925582011/bacterial.genomes.complete.masked.fa";
#Readonly my $fastaFile => "/groups/seidman/michael/pathogen_pipeline_2014/viral.genomes.mp205.v012813.fasta";
Readonly my $fastaFile => "/home/sasha/pathogen_pipeline_2014/genomes/viral.genomes/viral.genomes.mp205.v012813.fasta";
#Readonly my $faiFile => "/groups/seidman/michael/pathogen_pipeline_2014/viral.genomes.mp205.v012813.fasta.fai";
Readonly my $faiFile => "/home/sasha/pathogen_pipeline_2014/genomes/viral.genomes/viral.genomes.mp205.v012813.fasta.fai";
# .genome file for bedtools: genome short name - genome length
#Readonly my $genomeFile => "/files/Genetics/seidman/michael/tasha_TCGA_unmapped/bwa_index/bacterial_genomes/RM_907.WedMar91925582011/bacterial.genomes.complete.masked.genome";
#Readonly my $genomeFile => "/groups/seidman/michael/pathogen_pipeline_2014/viral.genomes.mp205.v012813.genome";
Readonly my $genomeFile => "/home/sasha/pathogen_pipeline_2014/genomes/viral.genomes/viral.genomes.mp205.v012813.genome";
#value for separate normalization
Readonly my $normSepValue => 1_000;
#value for general normalization (both genome length and read's number)
Readonly my $normGenValue => 1_000_000;
#value for coverage filter
my $coverageThreshold = 0;
#value for general coverage normalization 
Readonly my $normCovValue => 10_000_000;

(my $sec,my $min,my $hour,my $mday,my $mon,my $year,my $wday, my $yday,my $isdst)=localtime(time);
printf ("Start: %02d:%02d:%02d\n", $hour,$min,$sec);

my $fastqfname = $ARGV[0];
my $whole_sam_name = $fastqfname;
#$whole_sam_name =~ s/fq$/sam/;
$whole_sam_name =~ s/gz$/sam/;
my $sai_name = $whole_sam_name;
$sai_name =~ s/sam$/sai/;
## name for sam file with only mapped reads for the analysis
my $sam_name = $whole_sam_name;
$sam_name =~ s/sam$/vm.sam/;
my @name = split(/\./, $sam_name);
my $namePrefix = $name[0];
my $userSummaryFile = $namePrefix."-s";
my $additionBrief = ".b.sum";
my $summaryBriefFile = $userSummaryFile.$additionBrief;
my $additionFull = ".f.sum";
my $summaryFullFile = $userSummaryFile.$additionFull;
my $filtered_sam = $namePrefix."-f.sam";
my $final_clean_sam_file = $namePrefix."-f.c.sam";
my $reportFile = $namePrefix."-r";
my $log_file = $namePrefix.".log";
my $mult_hits_genomes_file = $namePrefix."-multgenomes";

#system "/files/Genetics/seidman/michael/save/bwa/bwa aln -l 32 -k 2 -t 4 $fastaFile $fastqfname > $sai_name";
system "bwa aln -l 32 -k 2 -t 4 $fastaFile $fastqfname > $sai_name";
#system "/files/Genetics/seidman/michael/save/bwa/bwa samse -n 1000 $fastaFile $sai_name $fastqfname > $whole_sam_name ";
system "bwa samse -n 1000 $fastaFile $sai_name $fastqfname > $whole_sam_name ";
system "samtools view -S -F 4 $whole_sam_name > $sam_name";

#open(FQFILE, "$fastqfname") || die "Can't open file: $!";
#while (<FQFILE>) {}
#my $fqsize = $.;
#my $readsnumber = $fqsize/4;
#close(FQFILE);
my $readsnumber = 100000; # just random number for trial purposes. I don't need it.

open my ($sam_file), '<', $sam_name or die "Cannot open sam file: $!";
my @uniqueHits = <$sam_file>;
close $sam_file;
my $confirmedUniqueSize = @uniqueHits;

my $percentConfirmedAligned = sprintf("%.2f%%", $confirmedUniqueSize*100/$readsnumber);

open my ($report), '>', $reportFile or die "Can't open report file for writing: $!";
say $report ("Original database file: $fastaFile\n");
say $report ("Original reads file: $fastqfname\n");
say $report ("Initial total number of reads: $readsnumber ");
say $report ("Total number of aligned reads: $confirmedUniqueSize ($percentConfirmedAligned) \n");
close $report;

my @arrayGenomeName; #contains full names of bacterial genomes from .fasta
my %genomePartFull; #hash Short name [AC_000014.1] -> Full name
my %fullNameAndHits; #hash Full name -> corresponding hits;
my %clean_full_name_and_hits; # hash of arrays: Full name -> corresponding hits;
# only for full summary report; hits without parsing multiple hits

open(INFOGENFILE, "$fastaFile") || die "Can't open file: $!";
while (<INFOGENFILE>) {
	if (/>gi/)  {
		push @arrayGenomeName, $_;
		my @separateGI = split /\|/;
		$genomePartFull{$separateGI[3]} = $_; #fill array Short name -> Full name
	}
}
close(INFOGENFILE);

my $sizeArrayGenomeName = @arrayGenomeName;	
say "Number of viral references: $sizeArrayGenomeName \n";

### I exclude filtering for RNAseq data.
my @filteredHits; #hits without 10 in a row
my $count_repeats = 0;
foreach (@uniqueHits) {
	my @hitString = split /\t/;
	#push @filteredHits, $_ if $hitString[4] !~ /A{10}|T{10}|C{10}|G{10}/;
	if ($hitString[9] !~ /A{10}|T{10}|C{10}|G{10}/) {
		push @filteredHits, $_;
	}
	else {
		push @filteredHits, $_;
        $count_repeats++;
	}
}
open(NOVEL, ">$filtered_sam") || die "Can't open file for writing: $!";
print NOVEL (@filteredHits);
close NOVEL;
say "Reads with 10 same bases in a row were filtered out: $count_repeats";
open $report, '>>', $reportFile or die "Can't open report file for writing: $!";
say $report "Number of filtered out reads with 10 same bases in a row: $count_repeats";
close $report;

($sec, $min, $hour, $mday, $mon, $year, $wday,  $yday, $isdst)=localtime(time);
printf ("Filtering 10 in a row: completed: %02d:%02d:%02d\n", $hour,$min,$sec);
my $sizeFilteredHits = @filteredHits;

# I need to filter out such reads that could be mapped on to different species:
# e.g., Pseudomonas stutzeri and Pseudomonas syringae.
# but I want to keep reads that could be mapped on to different strains:
# e.g., Staphylococcus epidermidis RP62A and Staphylococcus epidermidis ATCC 12228
# it will help me to get rid of conservative ribosome RNA sequences
# In order to do that I create a new array of useful reads, write it in a new SAM file that would be analyzed by genomeCoverageBed
my @reads_for_final_clean_sam_file;
my $reads_number_mulhit = 0;
my $xa_x0_test = $namePrefix.".xa_x0_test";
open my ($xa_result), '>', $xa_x0_test or die "Cannot open file for writing: $!";
open my ($mult_output), '>', $mult_hits_genomes_file or die "Cannot open file for writing: $!";
foreach (@filteredHits) {
		#an array for details of i hit from SAM file
		my $current_sam_line = $_;
		my @separateUHStringArray = split /\|/; # to extract part name (NC...)
		# fill the hash Full name -> corresponding main hit
		push @{$fullNameAndHits{$genomePartFull{$separateUHStringArray[3]}}}, $_;
		push @{$clean_full_name_and_hits{$genomePartFull{$separateUHStringArray[3]}}}, $_;  
		my $ocurred_genomes = $separateUHStringArray[3];
		# analyze multiple hits in XA tags
		my @sam_full_string = split /\t/;
		my $best_hits_number = $1 if ($current_sam_line =~ /X0:i:(\S+)/);
		
		#my @best_hits_number_line = split /:/, $sam_full_string[13]; #look on X0 tag
		#my $best_hits_number = $best_hits_number_line[2]; #take number of equally best hits
# in BWA parameters I use to report XA only if there are =< 1000 hits. If >1000 XA would be undef!
		if ($best_hits_number < 2 or $best_hits_number > 50) {
			push @reads_for_final_clean_sam_file, $current_sam_line; # instead of $_ - just more informative
		}
		# if multiple best hits
		else {
			print $mult_output "\n$sam_full_string[9]\n$genomePartFull{$separateUHStringArray[3]}";
			$reads_number_mulhit++;
			my %first_two_genome_names;
			my $current_full_name = $genomePartFull{$separateUHStringArray[3]};
			my @full_name = split / /, $current_full_name;
			my $first_two_words = $full_name[1]." ".$full_name[2];
			$first_two_genome_names{ $first_two_words } = 1; # just any value, I'm interested only in key
			my $cycle = 1;
			my $hit_in_a_row = 3;
			my $xa_string = $1 if ($current_sam_line =~ /XA:Z:(\S+)/);
            say $xa_result "BHN: $best_hits_number\n$current_sam_line" if (!$xa_string);
			#if ($sam_full_string[19] and $sam_full_string[19] =~ /XA/) {
			#	$xa_string = $sam_full_string[19];
			#}
			#elsif ($sam_full_string[18] and $sam_full_string[18] =~ /XA/) {
			#	$xa_string = $sam_full_string[18];
			#}
			while ($cycle < $best_hits_number) {
				my @multiple_hits = split /\|/, $xa_string; # look on XA tag
                if ($multiple_hits[$hit_in_a_row]) {
				$current_full_name = $genomePartFull{$multiple_hits[$hit_in_a_row]};
				@full_name = split / /, $current_full_name;
				$first_two_words = $full_name[1]." ".$full_name[2];
				$first_two_genome_names{ $first_two_words } = 1;
				push @{$fullNameAndHits{$genomePartFull{$multiple_hits[$hit_in_a_row]}}}, $current_sam_line;
				print $mult_output "$genomePartFull{$multiple_hits[$hit_in_a_row]}";
                }
				$cycle++;
				$hit_in_a_row += 4;
			}
			my $key_number = keys %first_two_genome_names;
			if ($key_number < 2) {
				push @reads_for_final_clean_sam_file, $current_sam_line; # the same genome or just different strains
			}
		}
}

close $mult_output;
close $xa_result;
open my ($final_sam_input), '>', $final_clean_sam_file or die "Cannot open file for writing: $!";
print $final_sam_input (@reads_for_final_clean_sam_file);
close $final_sam_input;
open $report, '>>', $reportFile or die "Can't open report file for writing: $!";
my $percent_mulhit = sprintf("%.2f%%", $reads_number_mulhit*100/$confirmedUniqueSize);
say $report "Number of reads with multiple hits: $reads_number_mulhit ($percent_mulhit)";
close $report;

#process .genome file
my %genomeHash; # Short name [NC_013590.2] -> genome length
open(GENOMEFILE, "$genomeFile") || die "Can't open file: $!";
while (<GENOMEFILE>) {
	chomp;
	(my $tmp0, my $tmp1, my $tmp2, my $tmpShortName, my $tmpTab, my $tmpGenomeLength) = split /[\||\t]/;
	$genomeHash{$tmpShortName} = $tmpGenomeLength;
}
close(GENOMEFILE);

# I changed format of full summary file using array with clean hits (one genome or diff. strains of the same species)
# otherwise it will be too big and useless in the case of parsing all multiple hits
open(NOVEL, ">$summaryFullFile") || die "Can't open file for writing: $!";
foreach my $string (keys %clean_full_name_and_hits) {
    say NOVEL ("$string@{$clean_full_name_and_hits{$string}}"); 
}

open(NOVEL, ">$summaryBriefFile") || die "Can't open file for writing: $!";
say NOVEL ("Hits\tHits/1000bp of genome\tHits/1000 initial reads\tNormalized hits\tGenome size\tVirus\n");
foreach my $string (sort { @{$fullNameAndHits{$b}} <=> @{$fullNameAndHits{$a}} } keys %fullNameAndHits) {
  my $sizeHitsArray = @{$fullNameAndHits{$string}};
  my @tmpName = split /\|/, $string;
  # number of mapped reads per 1000 bases of viral genome
  my $hitsNormPerGenome = "";
  if (defined($genomeHash{$tmpName[3]})) {
  $hitsNormPerGenome = sprintf("%.2f", ($sizeHitsArray*$normSepValue)/$genomeHash{$tmpName[3]});
 }
 # number of mapped reads per 1000 reads in the initial .fastq file
  my $hitsNormPerReads = "";
  $hitsNormPerReads = sprintf("%.4f", ($sizeHitsArray*$normSepValue)/$readsnumber);
  my $normGenHits = sprintf("%.4f", ($hitsNormPerGenome*$normGenValue)/$readsnumber);
  print NOVEL ("$sizeHitsArray\t$hitsNormPerGenome\t$hitsNormPerReads\t$normGenHits\t$genomeHash{$tmpName[3]}\t$string");
}

close NOVEL;
($sec, $min, $hour, $mday, $mon, $year, $wday,  $yday, $isdst)=localtime(time);
printf ("Redistribution of hits per each genome: completed: %02d:%02d:%02d\n", $hour,$min,$sec);	

#my $samFile = $filtered_sam.".sam";
#my $final_clean_sam_file = $namePrefix."-f.c.sam";
my $bamFile = $final_clean_sam_file.".bam";
my $sortedBamFile = $final_clean_sam_file.".sort";
my $sBFforIndex = $sortedBamFile.".bam";
my $graphCoverageFile = $sBFforIndex.".graph";
my $gt4CoverageFile = $graphCoverageFile.".gt4";

#---------------------------------------------------------------------------------------
if (-e $final_clean_sam_file) {
#without header -> so we need to use indexed ref for sam2bam convertation
#samtools faidx ref.fasta
system "samtools view -bt $faiFile  $final_clean_sam_file > $bamFile";
system "samtools sort $bamFile $sortedBamFile";
system "samtools index $sBFforIndex";
system "genomeCoverageBed -bg -ibam $sBFforIndex -g $genomeFile > $graphCoverageFile";
}
($sec, $min, $hour, $mday, $mon, $year, $wday,  $yday, $isdst)=localtime(time);
printf ("Sam->Bam->Graph: completed: %02d:%02d:%02d\n", $hour,$min,$sec);

my @gt4Array;
open (NOVEL, "$graphCoverageFile") || die "Can't open file: $!";
while (<NOVEL>) {
	my @string = split /\t/;
	push @gt4Array, $_ if $string[3] > $coverageThreshold;
}
close(NOVEL);
open (NEW, ">$gt4CoverageFile");
say NEW (@gt4Array);
close(NEW);

open GRAPHFILE, "$gt4CoverageFile" or die "Can't open file: $!";
my @graphArray = <GRAPHFILE>;
my $sizeGraphArray = @graphArray;
my $size = $sizeGraphArray-1;

my $outGraphFile = $gt4CoverageFile."c";

my @unitedGraphArray;
my $unitedString;
my $finalUnitedString;
my $startNum;
my $d = 0;
my $n = 0;
my $m = 0;
my %regionHashArray;
my $genome = "";
my $region = 0;
my $previousGenome = "";
my $previousRegion;

for ($d=0; $d<$sizeGraphArray-1; $d++) {
		my $singleString = $graphArray[$d];
		my @singleStringArray = split(/\t/, $singleString);
		my $nextSingleString = $graphArray[$d+1];
		my @nextSingleStringArray = split(/\t/, $nextSingleString);
		
		
		## only the first and one hit for one genome
		if (($singleStringArray[0] ne $nextSingleStringArray[0]) and ($singleStringArray[0] ne $previousGenome)) {
			$unitedString = $singleString;
			chomp($unitedString);
			my @genomeFullName = split /\|/, $singleStringArray[0];
			my $single_region = $singleStringArray[2]-$singleStringArray[1];
			my $procentOfGenome = sprintf("%.2f%%", ($single_region*100)/$genomeHash{$genomeFullName[3]});
			my $part = $genomePartFull{$genomeFullName[3]};
			chomp($part);
			my $finalUnitedString = $unitedString."\n".$part."\t".$single_region."\t".$procentOfGenome."\n\n";
			push(@unitedGraphArray, $finalUnitedString);
			#$previousGenome = $singleStringArray[0];
			$previousGenome = "";
		}
		## the last hit for current genome
		if (($singleStringArray[0] ne $nextSingleStringArray[0]) and ($singleStringArray[0] eq $previousGenome)) {
			$unitedString = $singleString;
			chomp($unitedString);
			my @genomeFullName = split /\|/, $singleStringArray[0];
			$region += $singleStringArray[2]-$singleStringArray[1];
			my $procentOfGenome = sprintf("%.2f%%", ($region*100)/$genomeHash{$genomeFullName[3]});
			my $part = $genomePartFull{$genomeFullName[3]};
			chomp($part);
			my $finalUnitedString = $unitedString."\n".$part."\t".$region."\t".$procentOfGenome."\n\n";
			push(@unitedGraphArray, $finalUnitedString);
			$previousGenome = $singleStringArray[0];
			$region = 0;
		}
		## non-consecutive hit for the same genome
		if (($singleStringArray[0] eq $nextSingleStringArray[0]) and ($singleStringArray[2] != $nextSingleStringArray[1])) {
			$unitedString = $singleString;
			#chomp($unitedString);
			$region += $singleStringArray[2]-$singleStringArray[1];
			push(@unitedGraphArray, $unitedString);
			$previousGenome = $singleStringArray[0];
		}
		## consecutive hits for same genome
		if (($singleStringArray[0] eq $nextSingleStringArray[0]) and ($singleStringArray[2] == $nextSingleStringArray[1])) {
			## choose max cov value: current or next)
			if ($singleStringArray[3] <= $nextSingleStringArray[3]) {
				$graphArray[$d+1] = "$singleStringArray[0]\t$singleStringArray[1]\t$nextSingleStringArray[2]\t$nextSingleStringArray[3]";
				$unitedString = "$singleStringArray[0]\t$singleStringArray[1]\t$nextSingleStringArray[2]\t$nextSingleStringArray[3]\n";
				$n++;
				$m=0;
			}
			else {
				$graphArray[$d+1] = "$singleStringArray[0]\t$singleStringArray[1]\t$nextSingleStringArray[2]\t$singleStringArray[3]";
				$unitedString = "$singleStringArray[0]\t$singleStringArray[1]\t$nextSingleStringArray[2]\t$singleStringArray[3]\n";
				$n++;
				$m=0;
			}
			$previousGenome = $singleStringArray[0];
		}
}

open (OUTGRAPH, ">$outGraphFile");
print OUTGRAPH (@unitedGraphArray);
close OUTGRAPH;
($sec, $min, $hour, $mday, $mon, $year, $wday,  $yday, $isdst)=localtime(time);
printf ("Coverage analysis: completed: %02d:%02d:%02d\n", $hour,$min,$sec);

my $summaryCoverage = $final_clean_sam_file.".cov";
my $tmp = $final_clean_sam_file.".cov.tmp";
system "grep '>gi' $outGraphFile > $tmp";

open (SUM, ">$summaryCoverage") or die "Can't open file: $!";
open (NOVEL, "$tmp") or die "Can't open file: $!";
while (<NOVEL>) {
	chomp;
	my @fullString = split /\t/;
	chop $fullString[2];
	my $normPercentCoverage = sprintf("%.3f", ($fullString[2]*$normCovValue)/$readsnumber);
	say SUM ("$_\t$normPercentCoverage");
}
close NOVEL;
close SUM;
my $summaryCoverageSort = $summaryCoverage.".sort";
system "sort -n -r -t \$'\t' -k 4,4 $summaryCoverage > $summaryCoverageSort";
system "rm -f $summaryCoverage $tmp";
system "mv $summaryCoverageSort $summaryCoverage";

## extract info from the gt4c files for the genomes with at least 4 covered regions.

my $sum_gt3_file = $outGraphFile.".ex";
my @cov_genome;
my %hash_gt4;
open my ($summary_cov), '<', $outGraphFile or die "Cannot open gt4c file: $!";
while (<$summary_cov>) {
	next if /^\s*$/; #next if blank line
	if (/^gi/) {
		push @cov_genome, $_;
	}
	if (/^>gi/) {
		my $size_cov_genome = @cov_genome;
	    my @sum_gi_line  = split /\t/;
        chomp(my $covered_space = $sum_gi_line[1]);
		$hash_gt4{$_} = [ @cov_genome ] if ($size_cov_genome > 3 or $covered_space >=200);
		@cov_genome = ();
	}
}
close $summary_cov;
open my ($summary_gt3_cov), '>', $sum_gt3_file or die "Cannot open file for writing: $!";
$" = ''; # change default list separator. Actually, a bad idea, but for the moment...
for my $genome ( keys %hash_gt4 ) {
    say $summary_gt3_cov "$genome@{ $hash_gt4{$genome} }";
}
close $summary_gt3_cov;
system "$ENV{'pathogen'}"."/ex_summary.pl $sum_gt3_file";
system "rm -f $whole_sam_name $sai_name $sam_name $final_clean_sam_file $reportFile $bamFile $sBFforIndex $graphCoverageFile $gt4CoverageFile $outGraphFile $summaryCoverage $sum_gt3_file ";

		
