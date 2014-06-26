#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my ($input, $output) = @ARGV;

my $dir = cwd();

$output = $dir."/".$output;
my @output;

open INPUT, $input or die "cannot open file $input.\n";

while (my $line =<INPUT>) {
	chomp $line;
	if( $line =~ m/\texon\t/ || $line =~ m/\ttranscript\t/ || $line =~ m/\tCDS\t/ || $line =~ m/\tstart\_codon\t/ || $line =~ m/\tstop\_codon\t/ || $line =~ m/\tgene\t/ ){
		push @output, $line."\n";
	}
}

close INPUT;

open OUTPUT, ">$output" or die "cannot open file $output.\n";
print OUTPUT @output;
close OUTPUT;

exit;