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
my $output = "NO"; #decide if output the lines of vcf that overlap with the feature searched

my $Results=GetOptions("vcf=s"=>\$vcf, "feature=s"=>\$feature_in, "variant_type=s"=>\$variant_type, "output=s"=>\$output );

unless ($variant_type eq "BOTH" || $variant_type eq "SNP" || $variant_type eq "INDEL") {
	print "Unrecognise variant_type. Accept values are SNP, INDEL and BOTH. Aborted! Retry!!\n";
	exit;
}

#readin features
my %Features;
my %features_avail;
open REF, "$reference" or die "Cannot open the ref file to read $reference\n";
while (my $line = <REF>) {
	if($line =~ m/^\#/) {
		next;
	} else {
		chomp $line;
		my @words = split ("\t", $line);
		my $chr = $words[0];
		my $feature = $words[2];
		if (!exists $features_avail{$feature}) {
			$features_avail{$feature} = 0
		}
		my $start = $words[3];
		my $end = $words[4];
		$Features{$chr}{$start}{$end}{$feature} = 0;
	}
}
close REF;

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


#read VCF positions and count number of SNPs and INDELs

my $total_SNPs_count = 0;
my $total_INDELs_count = 0;
my $SNP_on_feature_count = 0;
my $INDEL_on_feature_count = 0;
my %count_snp_overlap_with_feature; #this is used to count how many indel/snps overlapped with the feature.
my %count_indel_overlap_with_feature;

open VCF, "$vcf" or die "Cannot open the vcf file to read $vcf\n";
while (my $line = <VCF>) {
	if($line =~ m/^\#/) {
		next;
	} else {
		chomp $line;
		my @elements  = split ("\t", $line);
		my $chr=$elements[0];
		my $pos=$elements[1];
		my $R=$elements[3];
		my $A=$elements[4];
		my $V_type;
		
		if ($R =~ m/^[atgcATGC]$/ && $A =~ m/^[atgcATGC]$/) { # confirm it's a SNP not INDEL
			$V_type = "SNP";
			$total_SNPs_count += 1;
		} else {
			$V_type = "INDEL";
			$total_INDELs_count += 1;
		}
		
		my @features_on;
		print "this is the feature hash loop start \n";
		start_loop: foreach my $start ( sort {$a<=>$b} keys %{$Features{$chr}} ) {
			print $start."\n";
			foreach my $end (keys %{$Features{$chr}{$start}}) {
				if ( $pos >= $start && $pos <= $end) {
					foreach my $feature_found ( keys %{$Features{$chr}{$start}{$end}} ) {
						# print "for $chr $pos found: $chr $start $end $feature_found\n";
						push @features_on, $feature_found;
					}
				} elsif ($pos < $start) {
					last start_loop;
				}
			}
		}
		print "this is the feature hash loop end \n";
		
		if (scalar @features_on != 0 ) {
			foreach my $feature_found (@features_on) {
				if ($V_type eq "SNP" && $feature_found eq $feature_in) {
					$SNP_on_feature_count += 1;
					print "Found a SNP on the feature: $feature_in: on $chr\t $pos\n";
					if (!exists $count_snp_overlap_with_feature{$chr."_".$pos}) {
						
						 $count_snp_overlap_with_feature{$chr."_".$pos} = 0;
					}
				} elsif ($V_type eq "INDEL" && $feature_found eq $feature_in){
					$INDEL_on_feature_count += 1;
					 print "Found a indel on the feature: $feature_in: on $chr\t $pos\n";
					if (!exists $count_indel_overlap_with_feature{$chr."_".$pos}) {
						 $count_indel_overlap_with_feature{$chr."_".$pos} = 0;
					}
				}
			}
		}
		
		
	}
}
close VCF;


print "Count number: ";
if ($variant_type eq "BOTH") {
	my $number_1 = $SNP_on_feature_count+$INDEL_on_feature_count; 
	my $number_2 = scalar( keys %count_snp_overlap_with_feature) + scalar(keys %count_indel_overlap_with_feature);
	print $number_1." features overlap with the variants\n";
	print $number_2." variants overlap with the feature";
} elsif ($variant_type eq "SNP") {
	my $number_1 = $SNP_on_feature_count; 
	my $number_2 = scalar( keys %count_snp_overlap_with_feature);
	print $number_1." features overlap with the SNPs\n";
	print $number_2." SNPs overlap with the feature";
} elsif ($variant_type eq "INDEL") {
	my $number_1 = $INDEL_on_feature_count; 
	my $number_2 = scalar( keys %count_indel_overlap_with_feature);
	print $number_1." features overlap with the INDELs\n";
	print $number_2." INDELs overlap with the feature";
}	
print "\n";

#output
if ($output =~ m/Y|y/) {
	my $output_file = $vcf."_overlap_with_".$feature_in.".vcf";
	open OUTPUT, ">$output_file" or die "Cannot open the file $output_file to output.\n";
	open VCF, "$vcf" or die "Cannot open the vcf file to read $vcf\n";
	while (my $line = <VCF>) {
		if($line =~ m/^\#/) {
			print OUTPUT $line;
			next;
		} else {
			chomp $line;
			my @elements  = split ("\t", $line);
			my $chr=$elements[0];
			my $pos=$elements[1];
			if ($variant_type eq "BOTH") {
				if (exists($count_snp_overlap_with_feature{$chr."_".$pos}) || exists($count_indel_overlap_with_feature{$chr."_".$pos}) ) {
					print OUTPUT $line."\n";
				}
			} elsif ($variant_type eq "SNP") {
				if (exists $count_snp_overlap_with_feature{$chr."_".$pos} ) {
					print OUTPUT $line."\n";
				}
			} elsif ($variant_type eq "INDEL") {
				if (exists $count_indel_overlap_with_feature{$chr."_".$pos} ) {
					print OUTPUT $line."\n";
				}
			}	
		}
	}
	close VCF;
	close OUTPUT;
}


exit;
