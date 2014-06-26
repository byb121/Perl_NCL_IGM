#!/usr/bin/perl
use strict;
use warnings;

my ($file, $the_word)= @ARGV;

print "Reading from $file and the word $the_word within double quotes and a tab delimit will be inserted to the first line.\n";

my $outputFile = $file."_a_word_inserted.txt";
my @output;

open FILE, $file or die "cannot open file $file.\n";

my $line_number = 1;

while (my $line =<FILE>) {
    
    if($line_number == 1){
        push @output, '"'.$the_word.'"'."\t".$line;
    } else {
    	push @output, $line;
    }
    
    $line_number += 1;
    
}

close(FILE);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);

print "Done!\n";
exit;