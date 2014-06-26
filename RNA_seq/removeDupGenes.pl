#!/usr/bin/perl
use strict;
use warnings;

my ($file)= @ARGV;

print "Reading from $file.\n";

my $outputFile = $file."_DuplicateGenesRemoved";
my @output;

open FILE, $file or die "cannot open file $file.\n";

while (my $line =<FILE>) {
    chomp $line;
    my @words = split(/\t/,$line);
    if($words[0] !~ m/_dup\d+/){
        push @output, $line."\n";
    }
}

close(FILE);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);

print "Done!\n";
exit;