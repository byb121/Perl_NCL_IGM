#!/usr/bin/perl
use strict;
use warnings;

my ($file)= @ARGV;

print "Reading from $file.\n";

my $outputFile = $file."_headers_shoterned.txt";
my @output;

my %gene_hash;

open FILE, $file or die "cannot open file $file.\n";

while (my $line =<FILE>) {
    chomp $line;
    if($line =~ m/^>(\w+\|\d+\|.*\|)\s.*\((.*)\),\s(.*)/){
		my $refseq = $1;
	
		if (exists $gene_hash{$refseq}) {
			print "Error: duplicated refseq entry for $refseq.\n";
		} else {
			$gene_hash{$refseq} = 1;
			push @output,">".$refseq."\n";
		}
		
    } else {
    	push @output,$line."\n";
    }
}

close(FILE);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);
print "Done!\n";

exit;
