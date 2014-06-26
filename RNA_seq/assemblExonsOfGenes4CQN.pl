#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my $fasta_file = "mart_export.txt";

my %gene;
my %exon_start;
my %exon_end;
my %exon_seq;

my $gene_name = "";
my $exon_name = "";
my $seq = "";

open FASTA, $fasta_file or die "Cannot open the file $fasta_file";
while (my $line= <FASTA>) {
	chomp $line;
	if ($line =~ m/^\>/){
		
		if ($seq ne "" && $exon_name ne ""){
			$exon_seq{$exon_name} = $seq;
		}
		
		my @elements = split(/\|/, $line);
		$gene_name = $elements[0].'|'.$elements[1].'|'.$elements[2].'|'.$elements[3];
		$exon_name = $elements[4];
		if(exists $gene{$gene_name}) {
			$gene{$gene_name} = $gene{$gene_name}.".".$exon_name;
		} else {
			$gene{$gene_name} = $exon_name;
		}
		$exon_start{$exon_name} = $elements[5];
		$exon_end{$exon_name} = $elements[6];
		$seq = "";
		
	} else {
		$seq = $seq.$line;
	}
}
$exon_seq{$exon_name} = $seq;

close(FASTA);

my %gene_assembled_seq;
foreach my $key (keys %gene) {
	print "assembl exons of gene $key\n";
	my @exons = split (/\./, $gene{$key});
	my $start_compare = 0;
	my $end_compare = 0;
	my @assemble_seq;
	my %sub_exon_start;
	foreach my $exon (@exons){
		$sub_exon_start{$exon} = $exon_start{$exon};
	}
	
	@exons = sort {$sub_exon_start{$a} <=> $sub_exon_start{$b}} keys %sub_exon_start;
	for (my $i = 0; $i < scalar(@exons); $i++){
		if($i==0){
			push @assemble_seq, $exon_seq{$exons[$i]};
		} else {
			if ($exon_start{$exons[$i]} <= $exon_end{$exons[$i-1]}){
				print "overlap exons found $exons[$i] and $exons[$i-1]\n";
				if($exon_end{$exons[$i]} <= $exon_end{$exons[$i-1]}) {
					push @assemble_seq, "";
					$exon_end{$exons[$i]} = $exon_end{$exons[$i-1]};
				} else {
					my $extra_length = $exon_end{$exons[$i]} - $exon_end{$exons[$i-1]};
					my $extra_start = length($exon_seq{$exons[$i]}) - $extra_length; 
					my $extra_bases=substr($exon_seq{$exons[$i]}, $extra_start, $extra_length);
					push @assemble_seq,$extra_bases;
				}
			} else {
				push @assemble_seq, $exon_seq{$exons[$i]};
			}
		} 
	}
	
	#transfer the array into just on seq
	my $final_seq;
	foreach my $ele (@assemble_seq){
		$final_seq = $final_seq.$ele;
	}
	
	$gene_assembled_seq{$key} = $final_seq;
}

my $output_fasta = "mart_export.txt.assembled_exons.fasta";
open OUTPUT, ">$output_fasta" or die "Can not open the file to write $output_fasta";
foreach my $key (keys %gene_assembled_seq) {
	print OUTPUT $key."\n";
	print OUTPUT substr($gene_assembled_seq{$key}, 0, 60, '')."\n" while (length($gene_assembled_seq{$key}));
}
close OUTPUT;

exit;
