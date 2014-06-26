#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $vcf_files;
my $reference_file="/users/a5907529/lustre/Yaobo/GenomeData/GATK_bundle/ucsc.hg19.4GATK.fasta";
my $output_file_AMR="AMR_RAF_table.txt";
my $output_file_ASN="ASN_RAF_table.txt";
my $output_file_AFR="AFR_RAF_table.txt";
my $output_file_EUR="EUR_RAF_table.txt";

my $Results=GetOptions('vcfs=s' => \$vcf_files);

print "##############################################################\n";
print "Warning: To save memory, there's no hash used for variants, meaning:\n";
print "               There's no way to tell if a duplicate exists.\n";
print "Warning: ONLY SNPs are used, INDELs are IGNORED.\n";
print "Warning: Variants with different ref base between vcf and hg19 will be filtered.\n";
print "##############################################################\n";

#get the reference genome
my %hg19;

print "Start to read in reference file.\n";
open HG19, "$reference_file" or die "Can not open the file $reference_file\n";
my $string = "";
my $current_chr = "";
while (my $line = <HG19> ) {
	chomp $line;
	if ($line =~ m/^\#/) {
		next;
	} elsif( $line =~ m/^\>/) {
		$line =~ s/^\>//;
		if ($current_chr eq "") {
			$current_chr = $line;
			print "Reading in $current_chr...\n";
		} else {
			$hg19{$current_chr} = $string;
			$string = "";
			$current_chr = $line;
			print "Reading in $current_chr...\n";
		}
	} else {
		$string = $string.$line;
	}
}
close HG19;
$hg19{$current_chr} = $string;
print "Reference is imported\n";


#Read VCFs and Get Frequency for postions recorded
#get the VCF files
#my %AMR_Freq;
#my %AFR_Freq;
#my %ASN_Freq;
#my %EUR_Freq;


open OUTPUT_AMR, ">$output_file_AMR" or die "Cannot open file $output_file_AMR to output. \n";
open OUTPUT_ASN, ">$output_file_ASN" or die "Cannot open file $output_file_ASN to output. \n";
open OUTPUT_AFR, ">$output_file_AFR" or die "Cannot open file $output_file_AFR to output. \n";
open OUTPUT_EUR, ">$output_file_EUR" or die "Cannot open file $output_file_EUR to output. \n";

my @vcfs=split(',', $vcf_files);
VCF_LOOP: foreach my $vcf (@vcfs) {
	print "Reading vcf file: $vcf\n";
	my @sample_select_index;
	open VCF, $vcf or die "Cannot open the vcf file $vcf\n";
	LINE_LOOP: while (my $line = <VCF> ) {
		if ($line =~ /^\#/) {
			 next;
		} else {
			chomp $line;
			my @words = split (/\t/, $line);
			if ($words[6] ne "PASS") { # only variants that passed the filter
				next LINE_LOOP;
			}
			if ($words[3] =~ m/^(A|T|G|C)$/i && $words[4] =~ m/^(A|T|G|C)$/i) {
				# compare the ref base to the reference
				my ($chr, $pos, $rs, $ref, $var);
				if ($words[0] =~ m/MT/) { ###################could be different on different sets
					$chr="chrM";
				} else {
					$chr="chr".$words[0];
				}
				$pos = $words[1];
				$rs = $words[2];
				$ref = $words[3];
				$var = $words[4];
				#my $hg19_ref = ;
				if (uc($ref) eq uc(substr($hg19{$chr}, $pos-1, 1))) {
					##############################start to count frequency of samples
					my $info = $words[7];
					my @temp_2 = split (";", $info);
					my ($amr, $asn, $afr, $eur);
					for (my $i=0;$i<scalar @temp_2; $i++) {
						if ($temp_2[$i] =~ m/AMR\_AF\=/) {
							$amr = $temp_2[$i];
							$amr =~ s/AMR\_AF\=//;
						} elsif ($temp_2[$i] =~ m/ASN\_AF\=/) {
							$asn = $temp_2[$i];
							$asn =~ s/ASN\_AF\=//;
						} elsif ($temp_2[$i] =~ m/AFR\_AF\=/) {
							$afr = $temp_2[$i];
							$afr =~ s/AFR\_AF\=//;
						} elsif ($temp_2[$i] =~ m/EUR\_AF\=/) {
							$eur = $temp_2[$i];
							$eur =~ s/EUR\_AF\=//;
						} 
					}
					
					#if (exists $AMR_Freq{$chr}{$pos}{$rs}{$ref}{$var} && $amr) {
						#print "Error: Duplicate position info on $chr: $pos. Last reported kept.\n"
					#} else {
						if ($amr) {
							print OUTPUT_AMR "$rs\t$ref".'/'."$var\t$chr\t$pos\t".sprintf("%.3f", 1-$amr)."\n";
						} #else {
						#	$AMR_Freq{$chr}{$pos}{$rs}{$ref}{$var} = 1;
						#}
					#}
					
					#if (exists $ASN_Freq{$chr}{$pos}{$rs}{$ref}{$var} && $asn) {
						#print "Error: Duplicate position info on $chr: $pos. Last reported kept.\n"
					#} else {
						if($asn) {
							print OUTPUT_ASN "$rs\t$ref".'/'."$var\t$chr\t$pos\t".sprintf("%.3f", 1-$asn)."\n";
							#$ASN_Freq{$chr}{$pos}{$rs}{$ref}{$var} = sprintf("%.3f", 1-$asn);
						} #else {
						#	$ASN_Freq{$chr}{$pos}{$rs}{$ref}{$var} = 1;
						#}
					#}
					
					#if (exists $AFR_Freq{$chr}{$pos}{$rs}{$ref}{$var} && $afr) {
					#	print "Error: Duplicate position info on $chr: $pos. Last reported kept.\n"
					#} else {
						if ($afr) {
							print OUTPUT_AFR "$rs\t$ref".'/'."$var\t$chr\t$pos\t".sprintf("%.3f", 1-$afr)."\n";
							#$AFR_Freq{$chr}{$pos}{$rs}{$ref}{$var} = sprintf("%.3f", 1-$afr);
						} #else {
						#	$AFR_Freq{$chr}{$pos}{$rs}{$ref}{$var} = 1;
						#}
					#}
					
					#if (exists $EUR_Freq{$chr}{$pos}{$rs}{$ref}{$var} && $eur) {
					#	print "Error: Duplicate position info on $chr: $pos. Last reported kept.\n"
					#} else {
						if($eur) {
							print OUTPUT_EUR "$rs\t$ref".'/'."$var\t$chr\t$pos\t".sprintf("%.3f", 1-$eur)."\n";
							#$EUR_Freq{$chr}{$pos}{$rs}{$ref}{$var} = sprintf("%.3f", 1-$eur);
						} #else {
						#	$EUR_Freq{$chr}{$pos}{$rs}{$ref}{$var} = 1;
						#}
					#}
					
				} else {
					print "Different ref from hg19. REF in vcf is: $ref, REF in hg19 is:".substr($hg19{$chr}, $pos-1, 1)."\n";
					next LINE_LOOP;
				}
			}		
		}
	}
	close VCF;
}

close OUTPUT_AMR;
close OUTPUT_ASN;
close OUTPUT_AFR;
close OUTPUT_EUR;


#print output;
#open OUTPUT, ">$output_file_AMR" or die "Cannot open file $output_file_AMR to output. \n";
#foreach my $chr (keys %AMR_Freq) {
#	foreach my $pos (sort  {$a <=> $b} keys %{$AMR_Freq{$chr}}) {
#		foreach my $rs (keys %{$AMR_Freq{$chr}{$pos}}) {
#			foreach my $ref (keys %{$AMR_Freq{$chr}{$pos}{$rs}}) {
#				foreach my $var (keys %{$AMR_Freq{$chr}{$pos}{$rs}{$ref}}) {
#					#print OUTPUT "$chr\t$pos\t$rs\t$ref\t$var\t".$AMR_Freq{$chr}{$pos}{$rs}{$ref}{$var}."\n";
#					print OUTPUT "$rs\t$ref".'/'."$var\t$chr\t$pos\t".$AMR_Freq{$chr}{$pos}{$rs}{$ref}{$var}."\n";
#				}
#			}
#		}
#	}
#} 
#close OUTPUT;

#open OUTPUT, ">$output_file_ASN" or die "Cannot open file $output_file_ASN to output. \n";
#foreach my $chr (keys %ASN_Freq) {
#	foreach my $pos (sort  {$a <=> $b} keys %{$ASN_Freq{$chr}}) {
#		foreach my $rs (keys %{$ASN_Freq{$chr}{$pos}}) {
#			foreach my $ref (keys %{$ASN_Freq{$chr}{$pos}{$rs}}) {
#				foreach my $var (keys %{$ASN_Freq{$chr}{$pos}{$rs}{$ref}}) {
#					#print OUTPUT "$chr\t$pos\t$rs\t$ref\t$var\t".$ASN_Freq{$chr}{$pos}{$rs}{$ref}{$var}."\n";
#					print OUTPUT "$rs\t$ref".'/'."$var\t$chr\t$pos\t".$ASN_Freq{$chr}{$pos}{$rs}{$ref}{$var}."\n";
#				}
#			}
#		}
#	}
#} 
#close OUTPUT;

#open OUTPUT, ">$output_file_AFR" or die "Cannot open file $output_file_AFR to output. \n";
#foreach my $chr (keys %AFR_Freq) {
#	foreach my $pos (sort  {$a <=> $b} keys %{$AFR_Freq{$chr}}) {
#		foreach my $rs (keys %{$AFR_Freq{$chr}{$pos}}) {
#			foreach my $ref (keys %{$AFR_Freq{$chr}{$pos}{$rs}}) {
#				foreach my $var (keys %{$AFR_Freq{$chr}{$pos}{$rs}{$ref}}) {
#					#print OUTPUT "$chr\t$pos\t$rs\t$ref\t$var\t".$AFR_Freq{$chr}{$pos}{$rs}{$ref}{$var}."\n";
#					print OUTPUT "$rs\t$ref".'/'."$var\t$chr\t$pos\t".$AFR_Freq{$chr}{$pos}{$rs}{$ref}{$var}."\n";
#				}
#			}
#		}
#	}
#} 
#close OUTPUT;

#open OUTPUT, ">$output_file_EUR" or die "Cannot open file $output_file_EUR to output. \n";
#foreach my $chr (keys %EUR_Freq) {
#	foreach my $pos (sort  {$a <=> $b} keys %{$EUR_Freq{$chr}}) {
#		foreach my $rs (keys %{$EUR_Freq{$chr}{$pos}}) {
#			foreach my $ref (keys %{$EUR_Freq{$chr}{$pos}{$rs}}) {
#				foreach my $var (keys %{$EUR_Freq{$chr}{$pos}{$rs}{$ref}}) {
#					#print OUTPUT "$chr\t$pos\t$rs\t$ref\t$var\t".$EUR_Freq{$chr}{$pos}{$rs}{$ref}{$var}."\n";
#					print OUTPUT "$rs\t$ref".'/'."$var\t$chr\t$pos\t".$EUR_Freq{$chr}{$pos}{$rs}{$ref}{$var}."\n";
#				}
#			}
#		}
#	}
#} 
#close OUTPUT;

exit;

