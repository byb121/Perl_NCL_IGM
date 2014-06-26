#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my ($pattern, $input, $output) = @ARGV;

my $dir = cwd();
$output = $dir."/".$output;

my @output;

print "\n";
print "####################################################################\n";
print "#   This script will lines matching a pattern from the input file  #\n"; 
print "#            Usage: .pl 'perl_pattern' input_file_path             #\n";
print "####################################################################\n";
print "\n";
print "\n";

print "Reading file: $input\n";

print "Removed lines:\n";
open INPUT, $input or die "cannot open file $input.\n";
while (my $line =<INPUT>) {
	chomp $line;
	if( $line !~ m/($pattern)/ ){
		push @output, $line."\n";
	} else {
		print $line."\n";
	}
}
close INPUT;

open OUTPUT, ">$output" or die "cannot open file $output.\n";
print OUTPUT @output;
close OUTPUT;

exit;

