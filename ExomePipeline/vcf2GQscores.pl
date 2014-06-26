#!/usr/bin/perl
use strict;
use warnings;

print "\n\n";
print "#################################################################\n";
print "take a input vcf file, and will only read the first sample column\n"; 
print "and output a new vcf file with 3 columns expended:\n";
print "genotype call quality scores (GQ),\n";
print "reference coverage\n";
print "variants covareage\n";
print "#################################################################\n";

my ($vcf) = @ARGV;

my @output;
open VCF, "$vcf" or die "Cannot open the file $vcf, check if it exists\n";
while ( my $line = <VCF> ) {
	chomp $line;
	if ($line =~ m/^\#/) {
		push @output, $line."\n";
	} else {
		my ($ref_cov, $var_cov, $GQ, $PL);
		my @temp = split("\t", $line);
		
		my $format = $temp[8];
		my @format_fields = split(":", $format);
		my $GQ_index;
		my $AD_index;
		my $GT_index;
		my $PL_index;
		for (my $i=0;$i<scalar(@format_fields);$i++){
			if ($format_fields[$i] =~ m/GQ/) {
				$GQ_index = $i;
				#print $line."$GQ_index\n";
			}
			if ($format_fields[$i] =~ m/GT/) {
				$GT_index = $i;
				#print $line."$GT_index\n";
			}
			if ($format_fields[$i] =~ m/AD/) {
				$AD_index = $i;
				#print $line."$AD_index\n";
			}
			if ($format_fields[$i] =~ m/PL/) {
				$PL_index = $i;
				#print $line."$AD_index\n";
			}
		}
		
		if (defined $GQ_index && defined $AD_index && defined $GT_index) {
			my $sample = $temp[9];
			if ($sample !~ m/\.\/\./) {
				my @fields = split(":", $sample);
				my $GT = $fields[$GT_index];
				$GQ = $fields[$GQ_index];
				my $AD = $fields[$AD_index];
				$PL="";				
				#print $line."\t$GT\t$GQ\t$AD\n";
				if ($GT =~ m/0\/1|1\/1/) {
					my @cov = split(",", $AD);
					$ref_cov = $cov[0];
					$var_cov = $cov[1];
					my @PL_temp = split(",", $fields[$PL_index]);
					foreach my $PL_tt (@PL_temp) {
						$PL=$PL."\t".$PL_tt;
					}
					$PL =~ s/^\t//;
				} else {
					$ref_cov = "NA";
					$var_cov = "NA";
					$GQ = "NA";
					$PL = "NA\tNA\tNA";
				}
			} else {
				$ref_cov = "NA";
				$var_cov = "NA";
				$GQ = "NA";
				$PL = "NA\tNA\tNA";
			}
		} else {
			$ref_cov = "NA";
			$var_cov = "NA";
			$GQ = "NA";
			$PL = "NA\tNA\tNA";
		}
		
		push @output, $line."\t".$ref_cov."\t".$var_cov."\t".$GQ."\t".$PL."\n";
	}
}
my $output_file = $vcf.".GQextracted.vcf";
open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
print OUTPUT @output;
close OUTPUT;

exit;
