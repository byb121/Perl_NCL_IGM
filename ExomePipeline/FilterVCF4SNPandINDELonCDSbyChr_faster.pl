#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

print "#################\n";
print "Warning: on the loci where multi alternative variants are observed, indel will overwrite snp\n";
print "#################\n";

my $reference="/users/a5907529/lustre/Yaobo/GenomeData/gencode.v14.annotation.gtf";
my $vcf;
my $feature_in="exon";
my $variant_type="BOTH";
my $output_file = "NO"; #output file name
my $INchr;

my $Results=GetOptions("vcf=s"=>\$vcf, "feature=s"=>\$feature_in, "variant_type=s"=>\$variant_type,  "output=s"=>\$output_file,  "chr=s"=>\$INchr);

unless ($variant_type eq "BOTH" || $variant_type eq "SNP" || $variant_type eq "INDEL") {
	print "Unrecognise variant_type. Accept values are SNP, INDEL and BOTH. Aborted! Retry!!\n";
	exit;
}

#readin features
my @Features;
my @Feature_starts;
my @Feature_ends;

my %features_avail;
open REF, "$reference" or die "Cannot open the ref file to read $reference\n";
while (my $line = <REF>) {
	if($line =~ m/^\#/) {
		next;
	} else {
		chomp $line;
		my @words = split ("\t", $line);
		my $chr = $words[0];
		if ($chr ne $INchr){
			next;
		} else {
			my $feature = $words[2];
			push @Features, $feature;
			if (!exists $features_avail{$feature}) {
				$features_avail{$feature} = 0
			}
			my $start = $words[3];
			push @Feature_starts, $start;
			my $end = $words[4];
			push @Feature_ends, $end;
		}
	}
}
close REF;

print scalar @Features."\n";
print scalar @Feature_starts."\n";
print scalar @Feature_ends."\n";

print "Available features are: \n";
my $green_light = "red"; # decide if keep going or not
foreach my $feature (keys %features_avail) {
	print $feature."\n";
	if ($feature_in eq $feature){
		$green_light = "green";
	}
}
if ($green_light eq "green") {
	print "Enquiry feature $feature_in is one of the above, counting.....\n";
} else {
	print "Enquiry feature $feature_in is NOT one of the above. Aborted!\n";
	exit;
}

#sort feature starts;
my @idx = sort { $Feature_starts[$a] <=> $Feature_starts[$b] } 0 .. $#Feature_starts;
@Features = @Features[@idx];
@Feature_starts = @Feature_starts[@idx];
@Feature_ends = @Feature_ends[@idx];

print scalar @Features."\n";
print scalar @Feature_starts."\n";
print scalar @Feature_ends."\n";

my $search_start_idx=0;

open OUTPUT, ">$output_file" or die "Cannot open the file $output_file to output.\n";
open VCF, "$vcf" or die "Cannot open the vcf file to read $vcf\n";
vcf_loop: while (my $line = <VCF>) {
	if($line =~ m/^\#/) {
		print OUTPUT $line."\n";
		next vcf_loop;
	} else {
		chomp $line;
		my @elements  = split ("\t", $line);
		
		my $chr;
		if ( $elements[0] eq "MT" ){
			$chr = "chrM";
		} else {
			$chr = "chr".$elements[0];
		}
		
		if ($chr ne $INchr) {
			next vcf_loop;
		}
				
		my $pos=$elements[1];
		my $R=$elements[3];
		my $A=$elements[4];
		my $V_type;
		
		if ($R =~ m/^[atgcATGC]$/ && $A =~ m/^[atgcATGC]$/) { # confirm it's a SNP not INDEL
			$V_type = "SNP";
			#print $V_type."\n";
			#$total_SNPs_count += 1;
		} else {
			$V_type = "INDEL";
			#print $V_type."\n";
			#$total_INDELs_count += 1;
		}
		
		my @features_on;
		#print "This is the start of the feature hash loop.\n";
		my $next_search_idx;
		start_loop: for (my $i=$search_start_idx;$i<scalar @Feature_starts;$i++ )  {
			if ($pos >= $Feature_starts[$i] && $pos <= $Feature_ends[$i]) {
				unless ($next_search_idx) {
					$next_search_idx = $i;
					$search_start_idx = $next_search_idx;
				}
				push @features_on, $Features[$i];
			} elsif ($pos < $Feature_starts[$i]) {
				last start_loop;
			}
		}
		#print "This is the end of the feature hash loop.\n";
		
		if (scalar @features_on != 0 ) {
			foreach my $feature_found (@features_on) {
				if ($variant_type eq "BOTH" && $feature_found eq $feature_in ) {
					print OUTPUT $line."\n";
					next vcf_loop;
				} elsif  ($variant_type eq "INDEL" || $variant_type eq "SNP") {
					if ($V_type eq $variant_type && $feature_found eq $feature_in) {
						print OUTPUT $line."\n";
						next vcf_loop;
					}
				}
			}
		}
	}
}
close VCF;
close OUTPUT;

exit;

