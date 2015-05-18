#!/usr/bin/perl
use strict;
use warnings;

my ($file)=@ARGV;
my $output=$file."MQ.extracted.txt";
open IN, "$file" or die "Cannot open the file\n";
open OUT, ">$output";
while (my $line = <IN> ) {
	if($line !~ m/\#/) {
		if ($line =~ /.*MQ\=(\d+\.*\d+)\W/) {
			print OUT $1."\n";
		}
	}
}
close OUT;
close IN;
exit;
