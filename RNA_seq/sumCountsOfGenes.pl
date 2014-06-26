#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my ($file) = @ARGV;
my $output_file=$file."_summed_by_gene.txt";
my %lines;
my %new_counts;

open INPUT, $file or die "Cannot open the file $file";
while (my $line= <INPUT>) {
	chomp $line;
	my @elements = split("\t", $line);
	my $temp_value = $elements[1];
	for (my $i=2; $i < scalar(@elements); $i++){
		$temp_value = $temp_value."\t";
		$temp_value = $temp_value.$elements[$i];
	}
	$lines{$elements[0]} = $temp_value;
	#print $elements[0]."\n";
	#print $temp_value."\n";
}
close INPUT;

foreach my $key (keys %lines) {
	my @elements = split(/\|/, $key); 
	my $gene_name = $elements[0];
	print $gene_name."\n";
	if ( ! exists $new_counts{$gene_name}) {
		$new_counts{$gene_name} = $lines{$key};
	} else {
		my @counts_1 = split("\t", $new_counts{$gene_name});
		my @counts_2 = split( "\t", $lines{$key});
		my $temp="";
		for (my $i = 0; $i < scalar(@counts_1); $i++){
			my $sum = $counts_1[$i] + $counts_2[$i];
			$temp = $temp.$sum."\t";
		}
		$temp =~ s/\t$//;
		$new_counts{$gene_name} = $temp;
	}
}

open OUTPUT, ">$output_file" or die "Cannot open the file $output_file";
foreach my $key (keys %new_counts) {
	print OUTPUT $key."\t".$new_counts{$key}."\n";
}
close OUTPUT;

exit;