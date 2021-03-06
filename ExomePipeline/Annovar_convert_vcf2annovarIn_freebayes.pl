#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

#Program to convert vcf from GATK pipeline to correct format for Annovar
# YX edits
# change file writer append to create new '>>' to '>'
# so when the program stops unexpectedly, new files will be generated
# corrected genotype call for multiple alternative alleles

my $file ="";

my $Results=GetOptions("infile=s"=>\$file);

my $outfile = $file."_avinput.tsv";
open(OUT, ">$outfile") || die "Cannot open file \"$outfile\" to write to!\n";
my $Line;
open INPUT, $file or die "Cannot open $file";
line: foreach $Line (<INPUT>){
	if($Line=~/^#/){next line;}
	chomp $Line;
	my @linesplit = split (/\t/,$Line);
	my $chr = $linesplit[0];
	$chr =~ s/chr//;
	my $pos_1 = $linesplit[1];
	my $pos_2 = $pos_1;
	my $ref = $linesplit[3];
	my $var = $linesplit[4];
	my $format = $linesplit[8];
	
	my $Sample_Call="";
	for(my $i=9;$i<scalar(@linesplit);$i++) {
		$Sample_Call=$Sample_Call."\t".$linesplit[$i];
	}
	$Sample_Call =~ s/^\t//;
	
	my @vars=();
	if($var!~/,/){push(@vars,$var);}
	if($var=~/,/){
		@vars=split(/,/,$var);
	}
	
	if($ref=~/,/){print "Warnning: \n$Line\n";} #check for multiple ref alleles!
	foreach my $v (@vars){
		my $ref_now = $ref;
		my $v_now = $v;
		my $pos_1_now = $pos_1;
		my $pos_2_now = $pos_2;
		my $Sample_Call_Processed = AddGenoTypeToSampleCalls_CompondHet($v, $var, $format, $Sample_Call);
		my $out_var=$v;
		for (my $i=0;$i<scalar @vars;$i++) {
			if($v ne $vars[$i]) {
				$out_var=$out_var.",".$vars[$i];
			}
		}
		
		if ($ref_now =~ m/^[atgcATGC]$/) { # singlge base ref
			#insertion
			if(length($v_now)>length($ref_now)){
				#my @ref_elements=split("", $ref_now);
				my @var_chars=split(//, $v_now);
				if ($ref_now eq $var_chars[0]) {
					$ref_now="-";$v_now=~s/^.//;$pos_1_now++;$pos_2_now=$pos_1_now;
				}
			} #remove first 'ref' base if it is the same as first base of the variant
		} else { # multiple bases ref
			my @var_chars=split(//, $v_now);
			my @ref_chars=split(//, $ref_now);
			if(length($v_now) >= length($ref_now)){ #insertion
				Compare1:	for(my $i=0;$i<scalar @ref_chars;$i++) {
					if($ref_chars[$i] eq $var_chars[$i]) {
						if ($i != scalar @ref_chars-1) {
							my @temp_chars = split(//, $ref_now);
							shift @temp_chars;
							$ref_now=join '', @temp_chars;
						} else {
							$ref_now="-";
						}
						my @temp_chars = split(//, $v_now);
						shift @temp_chars;
						$v_now=join '', @temp_chars;
						$pos_1_now++;
					} else {
						last Compare1;
					}
				}
				$pos_2_now = $pos_1_now + length($ref_now) - 1;
			} else { #deletion
				Compare2:	for(my $i=0;$i<scalar @var_chars;$i++) {
					if($var_chars[$i] eq $ref_chars[$i]) {
						if ($i != scalar @var_chars -1) {
							my @temp_chars = split(//, $v_now);
							shift @temp_chars;
							$v_now=join '', @temp_chars;
						} else {
							$v_now="-";
						}
						my @temp_chars = split(//, $ref_now);
						shift @temp_chars;
						$ref_now=join '', @temp_chars;
						$pos_1_now++;
					} else {
						last Compare2;
					}
				}
				$pos_2_now = $pos_1_now + length($ref_now) - 1;
			}
		}
		
		print OUT "$chr\t$pos_1_now\t$pos_2_now\t$ref_now\t$v_now\t$linesplit[0]\t$linesplit[1]\t$linesplit[2]\t".
           "$linesplit[3]\t$out_var\t$linesplit[5]\t$linesplit[6]\t$linesplit[7]\t$linesplit[8]\t$Sample_Call_Processed\n";
	}
}
close OUT;
close INPUT;

exit;


sub AddGenoTypeToSampleCalls_CompondHet {
	my ($first_V, $variants, $format, $sample_call) = @_;
	my @format_fields = split(":", $format);
	my $GT_index;
	for (my $i=0;$i<scalar @format_fields ;$i++){
		if ($format_fields[$i] =~ m/GT/) {
			$GT_index = $i;
		}		
	}
	my $sample_call_processed=""; # vcf_sample1 \t geno_sample1 \t vcf_sample2 \t geno_sample2 ....
	my @sample_call_split = split("\t", $sample_call);
	for(my $i=0;$i<scalar(@sample_call_split);$i++){
		my $sample = $sample_call_split[$i];
		if ($sample !~ m/\.\/\./) {
			#print $sample."\n";
			my @fields = split(":", $sample);
			my $GT = $fields[$GT_index];
			if( $first_V eq $variants) {
				$sample_call_processed = $sample_call_processed."\t".$sample;
			} else {
				#determine which alternative allele is moved to the front
				my @V = split("," , $variants);
				my $V_index;
				for(my $h=0;$h<scalar @V;$h++) {
					if($first_V eq $V[$h]) {
						$V_index = $h+1;
					}
				}
				
				my %temp_hash;
				$temp_hash{$V_index} = 1;
				
				my $j = 2;
				for(my $h=1;$h<=scalar @V;$h++) {
					if ($h != $V_index) {
						$temp_hash{$h} = $j;
						$j+=1;
					}
				}
				
				my $left_number;
				my $right_number;
				if ($fields[$GT_index] =~ m/(\d)\/(\d)/ ) {
					my $A_left = $1;
					my $A_right = $2;
					#print "left $A_left\n";
					#print "right $A_right\n";
					if ($A_left == $A_right && $A_left ==0 ) {
						$left_number = "0";
						$right_number = "0";
					} elsif ( $A_left ==0 ) {
						$left_number = 0;
						$right_number = $temp_hash{$A_right};
					} else {
						$left_number = $temp_hash{$A_left};
						$right_number = $temp_hash{$A_right};
					}
				}
				###replace sample GT numbers
				$sample = $left_number."/".$right_number;
				for(my $h=1;$h < scalar @fields;$h++) {
					$sample=$sample.":".$fields[$h];
					my $teno = scalar @fields;
				}
				$sample_call_processed = $sample_call_processed."\t".$sample;
			}
		} else {
			$sample_call_processed = $sample_call_processed."\t".$sample;
		}
	}
	$sample_call_processed =~ s/^\t//;
	return $sample_call_processed;
}