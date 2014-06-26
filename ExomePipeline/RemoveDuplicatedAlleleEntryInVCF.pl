#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $in_vcf;

my $Results=GetOptions("in=s"=>\$in_vcf);
my $out="$in_vcf.processed.vcf";
my $rejected_v_file="$in_vcf.rejected.vcf";

my @output;
my @rejected;
open INPUT, $in_vcf or die "Cannot open $in_vcf\n";
line: foreach my $Line (<INPUT>) {
	chomp $Line;
	if($Line =~ m/^#/){
		push @output, $Line."\n";
	} else {
		my @ele = split("\t", $Line);
		my @Vs = split(",", $ele[4]);
		if (scalar @Vs > 1) {
			my %V;
			foreach my $v (@Vs) {
				if (exists $V{$v}) {
					push @rejected, $Line."\n";
					next line;
				} else {
					$V{$v} = 1;
				}
			}
			push @output, $Line."\n";
		} else {
			push @output, $Line."\n";
		}
	}	
}

close INPUT;

open(OUT, ">$out") or die "Cannot open file \"$out\" to write to!\n";
print OUT @output;
close OUT;

open(REJ, ">$rejected_v_file") or die "Cannot open file \"$rejected_v_file\" to write to!\n";
print REJ @rejected;
close REJ;

exit;