#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
use Math::CDF qw(pbinom);

my $input_file = "";

open INPUT, $input_file or die "Cannot open the file $input_file.";
while <INPUT> {
	chomp;
	my @fields = split /\t/;
	my $chr = $fields[0];
	my $coor = $fields[1]
	my $a = $fields[2];
	my $a = $fields[3];
	if ($fields[9] ne "N/A") {
		my $ref = $fields[9];
	}elsif ($fields[9] ne "N/A"){
		
	}elsif ($fields[9] ne "N/A"){
		
	} else {
		print "No sufficient coverage on this posion "
	}
	
	
}