#!/usr/bin/perl
use strict;
use warnings;

my ($file1,$file2) = @ARGV;
print " reading from $file1\n";


my %hash1;
my @output;
my $outputFile = $file1."_capitalizedGeneNames";

open FILE1, $file1 or die "cannot open file $file1";
while (my $line =<FILE1>) {
        chomp $line;
        my @words = split(/\t/,$line);
        my $tmp = $words[0];
        $tmp =~ tr/[a-z]/[A-Z]/;
        push @output, $tmp."\t".$line."\n";
}
close(FILE1);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);
exit;
