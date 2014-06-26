#!/usr/bin/perl
use strict;
use warnings;

my ($file1,$file2) = @ARGV;
print " reading from $file1\n";
print " reading from $file2\n";

my %hash1;
my %hash2;

my $outputFile = $file1."_sharedGenes";

open FILE1, $file1 or die "cannot open file $file1";
while (my $line =<FILE1>) {
        chomp $line;
        $line =~ tr/[a-z]/[A-Z]/;
        $line =~ s/\W//g;
        #my @words = split(/\s+/,$line);
        if (!exists $hash1{$line}) {
                        $hash1{$line} = 0;
                        
        }else {
                print "Duplicate entry found in $file1: ".$line."\n";
        }
}
close(FILE1);
my @output;
open FILE2, $file2 or die "cannot open file $file2";
while (my $line =<FILE2>) {
        chomp $line;
        $line =~ tr/[a-z]/[A-Z]/;
        $line =~ s/\W//g;
        #my @words = split(/\s+/,$line);
        if (!exists $hash1{$line}) {
                        $hash2{$line} = 0;
			print $line."\n";
        }else {
                #print "Duplicate entry found in $file2 : ".$line."\n";
                push @output, $line."\n";
        }
}
close(FILE2);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);
exit;
