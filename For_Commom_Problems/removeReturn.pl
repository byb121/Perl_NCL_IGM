#!/usr/bin/perl
use strict;
use warnings;

my ($file)= @ARGV;

print "Reading from $file.\n";

my $outputFile = $file."_return_Removed.txt";
my @output;

open FILE, $file or die "cannot open file $file.\n";

while (my $line =<FILE>) {
    
    chomp $line;
    
    if($line =~ m/^>/){
        push @output, "\n".$line." ";
    } else {
    	push @output, $line;
    }
    
}

close(FILE);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);

print "Done!\n";
exit;