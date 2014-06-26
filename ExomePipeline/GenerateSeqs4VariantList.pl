#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


my ($v_list) = @ARGV;
my $output_file = $v_list."_fasta.txt";

my $hg19_file="/home/yaobo/Documents/GenomeData/HumanGenome/hg19/hg19.fa";
my %hg19;

print "Start to read in reference file.\n";
open HG19, "$hg19_file" or die "Can not open the file $hg19_file\n";
my $string = "";
my $current_chr = "";
while (my $line = <HG19> ) {
	chomp $line;
	if ($line =~ m/^\#/) {
		next;
	} elsif( $line =~ m/^\>/) {
		$line =~ s/^\>//;
		if ($current_chr eq "") {
			$current_chr = $line;
			print "Reading in $current_chr...\n";
		} else {
			$hg19{$current_chr} = $string;
			$string = "";
			$current_chr = $line;
			print "Reading in $current_chr...\n";
		}
	} else {
		$string = $string.$line;
	}
}
close HG19;
print "Start to read in Variant list from $v_list.\n";

my @output;
open VARS, "$v_list" or die "Cannot open the file $v_list\n";
while (my $line = <VARS>) {
	chomp $line;
	if ($line =~ m/\#/) {
		next;
	} elsif ($line =~ m/^chr/) {
		my @elements = split ("\t", $line);
		my $chr = $elements[0];
		my $pos = $elements[1];
		my $upper_length;
		my $down_length;
		my $down_start;
		
		print "Variant found: $chr\t$pos\n";
		if ( ($pos-1-200) < 0 ) {
			$down_length = $pos-1;
			$down_start = 0;
		} else {
			$down_length = 200;
			$down_start = $pos-1-200;
		}
		
		if ($pos + 200 -1 >= (length $hg19{$chr}) -1 ) {
			$upper_length =  (length $hg19{$chr}) - 1 - $pos;
		} else {
			$upper_length = 200;
		}
		my $down_200 = substr $hg19{$chr}, $down_start , $down_length; # string in perl is 0-based
		my $upper_200 = substr $hg19{$chr}, $pos, $upper_length; #
		my $on_the_pos = substr $hg19{$chr}, $pos -1, 1;
		push @output, ">".$chr."_".$pos."\n".$down_200."\n".$on_the_pos."\n".$upper_200."\n";
	}
}
close (VARS);

open OUTPUT, ">$output_file" or die "Cannot open the output file $output_file\n";
print OUTPUT @output;
close OUTPUT;

exit;



