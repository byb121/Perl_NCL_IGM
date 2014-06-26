#!/usr/bin/perl
use strict;
use warnings;


my ($file1) = @ARGV;
print "##############################################\n";
print "Warning: To make the script worl correctly, P-value needs to be the second column in the file\n";
print "and the first column has NO blank\n";
print "##############################################\n";
print " reading from $file1\n";

my %hash1;
my $header;
my $outputFile = $file1."_DuplicateIlluminaIDRemoved";
my @output;
open FILE1, $file1 or die "cannot open file $file1";
while (my $line =<FILE1>) {
	if ($line =~ /value/i) {
		$header = $line;
	}
	my @words = split(/\t/,$line);
	if ($line =~ /.+(ILMN_\d+).+/){
		if (!exists $hash1{$1}) {
			$hash1{$1} = $line;
		}else {
			print "Duplicate entry found in $file1: ".$line."\n";
			my @words2 = split(/\t/,$hash1{$1});
			if ($words[1]<$words2[1]){
				$hash1{$1} = $line;
			}
		}
	}
}
close FILE1;

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT $header;
foreach my $key (keys %hash1) {
	print OUTPUT $hash1{$key};
}
close OUTPUT;

print "Done!\n";

exit;
