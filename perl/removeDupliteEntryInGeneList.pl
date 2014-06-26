#!/usr/bin/perl
use strict;
use warnings;

my $file = $ARGV[0];
my $colum_number = $ARGV[1];

print "Reading from $file, duplicated value on column $colum_number will be removed from the file! \n";
print "If this is not the column that you meant to, it's too late to change!\n";

my %hash;

my $outputFile = $file."_DuplicateRemoved.txt";
my @output;

open FILE, $file or die "cannot open file $file.\n";

while (my $line =<FILE>) {
    chomp $line;
    my @words = split(/\t/,$line);
	
	#if (@words < 12) {
	#	next;
	#}

	if (!exists $hash{$words[$colum_number-1]}) {
		print $words[$colum_number-1]."\n";
		#print $words[0]."\n";
        $hash{$words[$colum_number-1]} = 0;
        push @output, $line."\n";
    }else {
        print "Duplicate entry found in $file: ".$line."\n";
    }
}

close(FILE);

open OUTPUT, ">$outputFile" or die "cannont open $outputFile";
print OUTPUT @output;
close(OUTPUT);

print "Done!\n";
exit;
