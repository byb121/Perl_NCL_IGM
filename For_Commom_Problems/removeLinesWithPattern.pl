#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my ($pattern, $input, $output) = @ARGV;

my $dir = cwd();

$output = $dir."/".$output;

open INPUT, $input or die "cannot open file $input.\n";
open OUTPUT, ">$output" or die "cannot open file $output.\n";

while (my $line =<INPUT>) {
	chomp $line;
	if( $line =~ m/($pattern)/ ){
		next;
	} else {
		print OUTPUT $line."\n";
	}
}

close INPUT;
close OUTPUT;

exit;