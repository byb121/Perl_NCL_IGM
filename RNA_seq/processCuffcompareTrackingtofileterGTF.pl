#!/usr/bin/perl
use strict;
use warnings;

my ($trcaking_file, $gtf_file) = @ARGV;
my $output_gtf = $gtf_file."Transcripts60perSamples.gtf";


my %tcons_counts;
my %tcons_counts_NOF;
my %tcons_counts_OA;

open TRACKING, $trcaking_file or die "cannot open file $trcaking_file.\n";
while (my $line =<TRACKING>) {
	chomp $line;
	my @columns = split(/\t/, $line);
	my $tcons_id = $columns[0];
	if(exists $tcons_counts{$tcons_id}) {
		print "Error: duplicate TCONS is in tracking file $trcaking_file\n";
		exit;
	}
	my $cuff_gene_id = $columns[1];
	my $ref_gene_id = $columns[2];
	my $tcons_code = $columns[3];
	#my @expression;
	#for (my $i=4;$i< scalar @columns;$i++){
	#	push @expression, $columns[$i];
	#}
	$tcons_counts{$tcons_id} = 0;
	$tcons_counts_NOF{$tcons_id} = 0;
	$tcons_counts_OA{$tcons_id} = 0;
	for (my $i=4;$i< scalar @columns;$i++){
		if($columns[$i] ne "-" ) {
			if (exists $tcons_counts{$tcons_id}) {
				$tcons_counts{$tcons_id} += 1;
			} 
			
			if($i >= 4 && $i <= 9 ) {
				if (exists $tcons_counts_NOF{$tcons_id}) {
					$tcons_counts_NOF{$tcons_id} += 1;
				} 
			}
			
			if($i >= 10 && $i <= 19 ) {
				if (exists $tcons_counts_OA{$tcons_id}) {
					$tcons_counts_OA{$tcons_id} += 1;
				}
			}
		}
	}	
}
close(TRACKING);

open GTF, $gtf_file or die "cannot open file $gtf_file.\n";
open OUTPUT, ">$output_gtf" or die "cannot open the file $output_gtf to write.\n";
while (my $line =<GTF>) {
	chomp $line;
	my @columns = split(/\t/, $line);
	my @tags =  split(/;/, $columns[8]);
	foreach my $element (@tags) {
		if ($element =~ m/^\stranscript\_id\s\"(.+)\"/) {
			if($tcons_counts_NOF{$1} >= 4 || $tcons_counts_OA{$1} >= 6 ) {
				print OUTPUT $line."\n";
			}
		}		 
	}
}
close(OUTPUT);
close(GTF);

exit;
