#!/usr/bin/perl
use strict;
use warnings;

my $omni = "HumanOmni5-4-v1-0-D-auxilliary-file.txt";
my $chipSNP = "SNP_Table.txt";
my $dbsnp141 = "../dbSNP141_GRCh38.simplified.txt";
my @output;

my (%om, %ch, %db);

open DB, "$dbsnp141" or die "Cannot open the file $dbsnp141\n";
print "reading in $dbsnp141.\n";
while (my $line = <DB> ) {
	chomp $line;
	my @eles = split("\t", $line);
	if ($eles[0] eq "Index") {
		next;
	} else {
		if ( ! exists $db{$eles[2]}) {
			$db{$eles[2]} = $line;
		} 
	}
}
close DB;

print "reading in $omni.\n";
open OM, "$omni" or die "Cannot open the file $omni";
while (my $line = <OM> ) {
	chomp $line;
	my @eles = split("\t", $line);
	if ($eles[0] eq "Name") {
		next;
	} else {
		if ( ! exists $om{$eles[0]}) {
			$om{$eles[0]} = $eles[1];
		} else {
			print "Error, duplicated name $eles[0] in $omni\n";
			exit;
		}
	}
}
close OM;

push @output, "Index_chip\tName_chip\tChr_chip_h19\tPos_chip_h19\tSNP_chip\tA\tB\tRsID_chip\tChr_GRCh38\tPos_GRCh38\tID_dbsnp141\tRef_dbsnp141\tAlt_dbsnp141\n";


print "reading in $chipSNP.\n";
open CHIP, "$chipSNP" or die "Cannot open the file $chipSNP";
while (my $line = <CHIP> ) {
	chomp $line;
	my @eles = split("\t", $line);
	if ($eles[0] eq "Index") {
		next;
	} else {
		$eles[4] =~ s/\[|\]//g;
		my @ab = split("/", $eles[4]);
		my $string = $line."\t".$ab[0]."\t".$ab[1]."\t".$om{$eles[1]};
		if (exists $db{$om{$eles[1]}}) {
			$string = $string."\t".$db{$om{$eles[1]}}."\n";
		} else {
			$string = $string."\t"."NA"."\t"."NA"."\t"."NA"."\t"."NA"."\t"."NA"."\n";
		}
		push @output, $string;
	}
}
close OM;

my $output = "$chipSNP.processed.txt";
open OUT, ">$output" or die "Cannot open the file $output\n";
print OUT @output;
close OUT;
exit;
