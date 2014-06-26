#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

###############################################
# Input: exon list with postions and gene names on each row
# Output: bed format with ID on each row end similar to 
###############################################

my ($Exon_list_file)  = @ARGV;
my %Exon_List;

open ExonList, $Exon_list_file or die "ExonList2Bed.pl: Cannot open ".$Exon_list_file."\n";
while (my $line =<ExonList>) {
	chomp $line;
	#print $line."\n";
	$line =~ s/^\s//;
	if ($line =~ m/^chr/i) {
		my @words = split(/\s+/,$line);
		my $chr=$words[0];
		$chr =~ tr/C/c/;
		$chr =~ s/\W//g;
		my $start = $words[1];
		$start =~ s/\W//g;
		my $end = $words[2];
		$end =~ s/\W//g;
		my $strand = "+";
		
		if ($end <= $start) {
			$strand = "-";
			$start = $words[2];
			$start =~ s/\W//g;
			$end = $words[1];
			$end =~ s/\W//g;
		}
		$start -= 1;
		my $gene_name = $words[3];
		$gene_name =~ s/\W//g;
		my $length = $end-$start;
		$Exon_List{$chr."\t".$start."\t".$end} = $chr.":".$start."-".$end.":".$gene_name."\t".$length."\t".$strand; 
	}
}
close ExonList;

my $output = $Exon_list_file.".bed";
open OUTPUT,  ">$output" or die "ExonList2Bed.pl: Cannot open ".$output." to output\n";
foreach my $exon (sort keys %Exon_List) {
	print OUTPUT $exon."\t".$Exon_List{$exon}."\n";
}
close OUTPUT;
exit;






