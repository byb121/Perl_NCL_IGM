#!/usr/bin/perl
use strict;
use warnings;

my ($file1) = @ARGV;
print " reading from $file1\n";

my %hash1;

my $outputFile = $file1."_DuplicateRemoved";
my @output;
open FILE1, $file1 or die "cannot open file $file1";
while (my $line =<FILE1>) {
        chomp $line;
        my @words = split(/\t/,$line);
        if (!exists $hash1{$words[0]}) {
        	$hash1{$words[0]} = 0;
        	push @output, $line."\n";
        }else {
                print "Duplicate entry found in $file1: ".$line."\n";
        }
}
close(FILE1);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);
exit;
