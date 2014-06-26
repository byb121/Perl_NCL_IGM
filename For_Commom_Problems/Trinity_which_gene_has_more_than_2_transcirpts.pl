#!/usr/bin/perl
use strict;
use warnings;

my ($fasta_file)= @ARGV;

my %gene_transcript_count;

open FASTA, $fasta_file or die "cannot open file $fasta_file.\n";
while (my $line =<FASTA>) {
	chomp $line;
	if($line =~ m/^>/){
		my @header = split(' ',$line);
		my @g_g_seq = split('_',$header[0]);
		my $key = $g_g_seq[0]."_".$g_g_seq[1];
		if (exists $gene_transcript_count{$key}) {
			$gene_transcript_count{$key} =  $gene_transcript_count{$key} + 1;
		} else {
			$gene_transcript_count{$key} =1;
		}
	}
}
close FASTA;

foreach my $key (keys %gene_transcript_count) {
	if ($gene_transcript_count{$key} >= 2) {
		print $key."\n";
	}
}
 
exit;