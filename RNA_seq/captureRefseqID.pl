#!/usr/bin/perl
use strict;
use warnings;

my ($file)= @ARGV;

print "Reading from $file.\n";

my $outputFile = $file."_headers_processed.txt";
my @output;

my %gene_hash;
my %count_hash;

open FILE, $file or die "cannot open file $file.\n";

while (my $line =<FILE>) {
    chomp $line;
    if($line =~ m/^>(\w+\|\d+\|.*\|)\s.*\((.*)\),\s(.*)/){
		my $refseq = $1;
		my $symbol = $2;
		my $type = $3;
		
		if (exists $count_hash{$symbol}) {
			$count_hash{$symbol} += 1;
		} else {
			$count_hash{$symbol} = 1;
		}
		
		if (exists $gene_hash{$refseq}) {
			print "Error: duplicated refseq entry for $refseq.\n";
		} else {
			$gene_hash{$refseq} = 1;
			my $temp = $symbol.".".$count_hash{$symbol};
			push @output,$refseq."\t".$temp."\t".$type."\n"; 
		}
		
    } else {
    	print "$line is not a match.\n";
    }
}

close(FILE);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);
print "Done!\n";

exit;
