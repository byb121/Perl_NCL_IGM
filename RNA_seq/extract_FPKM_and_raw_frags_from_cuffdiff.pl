#!/usr/bin/perl
use strict;
use warnings;

my ($group_tracking_file, $diff_file) = @ARGV;

my $frags_output = $group_tracking_file."_merged_with_diff_frags.txt";
my $FPKM_output = $group_tracking_file."_merged_with_diff_FPKM.txt";
my @frags_output;
my @FPKM_output;

my %group_frags_hash;
my %group_FPKM_hash;
my %diff_hash;

open GROUP_INFO, $group_tracking_file or die "cannot open file $group_tracking_file.\n";
while (my $line =<GROUP_INFO>) {
	my @words = split("\t", $line);
	if (exists $group_frags_hash{$words[0]}) {
		$group_frags_hash{$words[0]} = $group_frags_hash{$words[0]}."\t".$words[3];
	} else {
		$group_frags_hash{$words[0]} = $words[3];
	}
}
close(GROUP_INFO);

open GROUP_INFO, $group_tracking_file or die "cannot open file $group_tracking_file.\n";
while (my $line =<GROUP_INFO>) {
	my @words = split("\t", $line);
	if (exists $group_FPKM_hash{$words[0]}) {
		$group_FPKM_hash{$words[0]} = $group_FPKM_hash{$words[0]}."\t".$words[6];
	} else {
		$group_FPKM_hash{$words[0]} = $words[6];
	}
}
close(GROUP_INFO);

open DIFF_INFO, $diff_file or die "cannot open file $diff_file.\n";
while (my $line =<DIFF_INFO>) {
	chomp $line;
	my @words = split("\t", $line);
	if (exists $diff_hash{$words[0]}) {
		next;
	} elsif ($words[6] eq "OK") {
		$diff_hash{$words[0]} = $line;
	} else {
		next;
	}
}
close(DIFF_INFO);


open FRAGS_OUTPUT, ">$frags_output" or die "cannot open the file $frags_output";
for my $key ( keys %diff_hash ) {
	print FRAGS_OUTPUT $diff_hash{$key}."\t".$group_frags_hash{$key}."\n";
}
close (FRAGS_OUTPUT);

open FPKM_OUTPUT, ">$FPKM_output" or die "cannot open the file $FPKM_output";
for my $key ( keys %diff_hash ) {
	print FPKM_OUTPUT $diff_hash{$key}."\t".$group_FPKM_hash{$key}."\n";
}
close (FPKM_OUTPUT);


