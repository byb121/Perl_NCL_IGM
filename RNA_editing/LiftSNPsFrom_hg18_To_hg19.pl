#!/usr/bin/perl
use strict;
use warnings;
#### Yaobo ####
my $VCF_hg19="/users/a5907529/lustre/Yaobo/GenomeData/GATK_bundle/dbsnp_137.hg19.4GATK.vcf";
my ($input, $output) = @ARGV;

my %V_hg19;
open HG19, $VCF_hg19 or die "Can not open the file $VCF_hg19";
while (my $line = <HG19>) {
	if($line =~ m/^\#/) {
		next;
	} else {
		chomp $line;
		my @words = split("\t", $line);
		my $chr = $words[0];
		my $pos = $words[1];
		if( ! exists $V_hg19{$words[2]}) {
			$V_hg19{$words[2]} = $chr."_".$pos;
		} else {
			print "Duplicate entry of $words[2] found on line $line.\n";
		}
	}
}
close HG19;

open OUTPUT, ">$output" or die "Can not open the output file $output";
open HG18, "$input" or die "Cannot open the file $input";
while (my $line = <HG18>) {
	if ( $line !~ /^rs/) {
		next;
	} else {
		chomp $line;
		my @words = split ("\t", $line);
		if (exists $V_hg19{$words[0]}) {
			my @temp = split ("_", $V_hg19{$words[0]});
			print OUTPUT $words[0]."\t".$words[1]."\t".$temp[0]."\t".$temp[1]."\t".$words[4]."\n";
		} else {
			print "one entry $words[0] is not found in hg19.\n"
		}
	}
}
close HG18;
close OUTPUT;
exit;