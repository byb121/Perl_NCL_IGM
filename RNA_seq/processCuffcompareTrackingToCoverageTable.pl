#!/usr/bin/perl
use strict;
use warnings;

my ($trcaking_file) = @ARGV;
my $output_table = $trcaking_file.".coverageTable.txt";


open TRACKING, $trcaking_file or die "cannot open file $trcaking_file.\n";
open TABLE, ">$output_table" or die "cannot open file $output_table.\n";

while (my $line =<TRACKING>) {
	chomp $line;
	my @columns = split(/\t/, $line);
	my $tcons_id = $columns[0];
	
	my $cuff_gene_id = $columns[1];
	my $ref_gene_id = $columns[2];
	my $tcons_code = $columns[3];
	
	if ($ref_gene_id ne "-") {
		my @temp = split(/\|/, $ref_gene_id);
		print TABLE  $tcons_id."\t".$cuff_gene_id."\t".$temp[0]."\t".$temp[1]."\t".$tcons_code."\t";
	} else {
		print TABLE  $tcons_id."\t".$cuff_gene_id."\t"."-"."\t"."-"."\t".$tcons_code."\t";
	}
	
	for (my $i=4;$i< scalar @columns;$i++){
		if($columns[$i] ne "-" ) {
			my @words = split ( /\|/, $columns[$i]);
			if ($i == (scalar @columns -1)) {
				print TABLE $words[6]."\n";
			} else {
				print TABLE $words[6]."\t";
			}
		} else {
			if ($i == (scalar @columns -1)) {
				print TABLE "-"."\n";
			} else {
				print TABLE "-"."\t";
			}
		}
	}
}

close TABLE;
close TRACKING;

exit;
