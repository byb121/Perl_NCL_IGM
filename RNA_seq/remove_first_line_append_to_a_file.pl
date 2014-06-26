#!/usr/bin/perl
use strict;
use warnings;

my ($file, $outputFile)= @ARGV;

print "Reading from $file.\n";

my @output;

open FILE, $file or die "cannot open file $file.\n";

my $line_number = 1;

while (my $line =<FILE>) {
    
    if($line_number != 1){
    	push @output, $line;
    }
    
    $line_number += 1;
    
}

close(FILE);

open OUTPUT, ">>$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);

print "Done!\n";
exit;