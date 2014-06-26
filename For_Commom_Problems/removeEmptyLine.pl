#!/usr/bin/perl
use strict;
use warnings;

my ($file) = @ARGV;
my @output;
my $count=0;

open INPUT, $file or die "Cannot open the file $file";
while (my $line= <INPUT>) {
	if ($line !~ m/^\n/) {
		push @output, $line;
	} else {
		$count += 1;
		print "remove the line: ".$line;
	}
}

close INPUT;

my $output_file = $file."_emptyLineRemoved.txt";
open OUTPUT, ">$output_file" or die "Cannot open file to outout in $output_file";
print OUTPUT @output;
close OUTPUT;

print "Total removed line number: $count \n";
exit;