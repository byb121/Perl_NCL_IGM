#!/usr/bin/perl
use strict;
use warnings;

my ($file)= @ARGV;

print "Reading from $file.\n";

my $outputFile = $file."_GeneNameTailsRemoved.txt";
my @output;

open FILE, $file or die "cannot open file $file.\n";

while (my $line =<FILE>) {
    chomp $line;
    if( $line =~ m/^(\w+)\.\d+.*/){
    	my $gene_name = $1;
    	push @output, $gene_name."\t".$line."\n";
    } else {
    	print $line."\n";
    }
    if( $line =~ m/^(\w+\-\w+)\.\d+\t.*/){
    	my $gene_name = $1;
    	push @output, $gene_name."\t".$line."\n";
    }
}

close(FILE);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);

print "Done!\n";
exit;
