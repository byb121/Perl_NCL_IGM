#!/usr/bin/perl
use strict;
use warnings;

#my $temp = '/home/data/AllChromosomes/hg18.fa';
my $temp = "shengqide xiaoren";
print "temp is ", $temp, "\n";
my $string =  readGenomeFile($temp);
print "string is ",$string, "\n";
exit;

sub readGenomeFile { #this will only work with a single file contain all chromosomes in fasta format
        my ($file) = @_;
	print @_,"\n";
        print $file,"1212\n";
        my $chr = "chr";
        my $seq = "seq";
        return "done!!!\n";
}

