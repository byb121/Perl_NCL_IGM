#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $in_vcf;
my $outDir;

my $Results=GetOptions("in=s"=>\$in_vcf, "outDir=s"=>\$outDir);
my @samples;
my %calls;

my %SampleNames;
my @vcfHeader;

# read in samples names and vcf headers
open INPUT, $in_vcf or die "Cannot open $in_vcf\n";
line: foreach my $Line (<INPUT>) {
	chomp $Line;
	if($Line =~ m/^#CHR/){
		my @linesplit = split (/\t/,$Line);
		my $temp_string="";
		for(my $i=0;$i<scalar(@linesplit);$i++) {
			if($i >= 9) {
				if (! exists $SampleNames{$linesplit[$i]} ) {
					$SampleNames{$linesplit[$i]} = 1;
					push @samples, $linesplit[$i];
				} else {
					print "VCF fatal error: duplicated sample name\n";
					exit;
				}
			} else {
				$temp_string = $temp_string."\t".$linesplit[$i];
			}
		}
		$temp_string =~ s/^\t//;
		push @vcfHeader, $temp_string;
		last;
	} else {
		if($Line =~ m/^#/){
			push @vcfHeader, $Line."\n";
		#} else {
			#my @linesplit = split (/\t/,$Line);
			#for(my $i=9;$i<scalar(@linesplit);$i++) {
				#$calls{$linesplit[0]}{$linesplit[1]}{$linesplit[2]}{$linesplit[3]}{$linesplit[4]}{$linesplit[5]}{$linesplit[6]}{$linesplit[7]}{$linesplit[8]}{$samples[$i-9]} = $linesplit[$i]; 
			#}
		}
	}
}
close INPUT;

@samples=("117-3_S21","117-4_S22",
"117-5_S23","117-6_S24",
"117-7_S25","142-3_S26",
"142-4_S27","14-4_S15",
"14-5_S16","14-6_S17",
"15-3_S18","15-4_S19",
"15-5_S20","252-3_S28",
"252-4_S29","D62458_S1",
"D62459_S6","D68317_S7",
"D68323_S8","D68324_S14",
"D68325_S9","D70368_S10",
"D70419_S11","D70420_S12"
,"D70421_S13","D75944_S4",
"D80618_S5","D80620_S2",
"D80624_S3");


$outDir =~ s/\/$//;
for (my $i=0; $i< scalar @samples; $i++) {
	my $output_file = $outDir."/".$samples[$i].".vcf";
	open OUTPUT, ">$output_file" or die "Cannot open the file to write $output_file\n";
	open INPUT, $in_vcf or die "Cannot open $in_vcf\n";
	my @output;
	print OUTPUT @vcfHeader;
	print OUTPUT "\t$samples[$i]\n";
	line: foreach my $Line (<INPUT>) {
		chomp $Line;
		if($Line !~ m/^#/){
			my @linesplit = split (/\t/,$Line);
			if ($linesplit[$i+9] !~ m/^\./ && $linesplit[$i+9] !~ m/^0\/0/) {
				print OUTPUT $linesplit[0]."\t".$linesplit[1]."\t".$linesplit[2]."\t".$linesplit[3]."\t".
										$linesplit[4]."\t".$linesplit[5]."\t".$linesplit[6]."\t".$linesplit[7]."\t".$linesplit[8]."\t".$linesplit[$i+9]."\n";
			}
		}
	}
	close INPUT;
	close OUTPUT;
}

exit;



