#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($in_file) = @ARGV;

my $output_file = $in_file.".intervals";
open BED, "$in_file" or die "Cannot open the file $in_file\n";
open OUTPUT, ">$output_file" or die "Cannot open the file to output $output_file\n";
while (my $line = <BED> ) {
	chomp $line;
	if ($line =~ m/^\#/) {
		print OUTPUT $line."\n";
		next;
	} else {
		my @words = split("\t", $line);
		my $chr = $words[2];
		my $pos = $words[3];
		print OUTPUT $chr.":".$pos."\n";
	}
}

close BED;
close OUTPUT;
exit;