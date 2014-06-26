#!/usr/bin/perl
use strict;
use warnings;

my $exon_seq = "AAAAAACAAAAA";
my $motif = "AAAA";
my $motif_count = 0;

while ($exon_seq =~ m/(?=($motif))/g) {
	print "start: @-\n";
	my $start = $-[1];
	print "start: $start\n";
	print "end: @+\n";
	print "end: $+[1]\n";
	my $string = substr $exon_seq, $start, (length $motif)+3;
	print "Substring: ".$string."\n";
	print $motif_count."\n";
	$motif_count += 1;
}

exit;