#!/usr/bin/perl 
use strict;
use warnings;

my ($file) = @ARGV;
my $output = $file."_position_extracted";

open INPUT, $file or die "Can not open file $file";
open OUTPUT, ">$output" or die "Can not open file $output";

while (my $line = <INPUT>) {
	if ($line =~ m/^rs/) {
		chomp $line;
		my @words = split(/\s\|\s/, $line);
		print OUTPUT $words[0]."\t";
	} elsif ($line =~ m/^CTG.+GRCh37/){
		chomp $line;
		my @words = split(/\s\|\s/, $line);
		my @tmp1 = split ("=", $words[2]);
		print OUTPUT $tmp1[0].$tmp1[1]."\t";
		my @tmp2 = split ("=", $words[3]);
		if ($tmp2[1] =~ /\d+/) {
			print OUTPUT $tmp2[1]."\t";
		}
	} elsif($line =~ m/^\n/) {
		print OUTPUT "\n";
	} else {
		next;
	}
}
close(INPUT);
close(OUTPUT);
exit;
