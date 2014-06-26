#!/usr/bin/perl
use strict;
use warnings;

my ($file1,$file2)= @ARGV;

print "Reading from $file1.\n";
my %gene_list;

open FILE1, $file1 or die "cannot open file $file1.\n";

while (my $line =<FILE1>) {
    chomp $line;
    #my @words = split(/\t/,$line);
    #print $words[0]."\n";
    print $line."\n";
    $line =~ s/\W//g;
    $gene_list{$line} = 1;
}

close(FILE1);

my $outputFile = $file2."_extratedEntry.txt";
my @output;


print "Reading from $file2.\n";
open FILE2, $file2 or die "cannot open file $file2.\n";

while (my $line =<FILE2>) {
    chomp $line;
    my @words = split(/\t/,$line);
    
    $words[0] =~ s/\W//g;
    
    print $words[0]."\n";
    if(exists $gene_list{$words[0]}) {
    	push @output, $line."\n";
    }
}
close(FILE2);


open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);

print "Done!\n";
exit;
