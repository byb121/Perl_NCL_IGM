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
	###correct chr name with no chr prefix
	if ($linesplit[0] =~ /^x$/) {$linesplit[0] = "X";}
	if ($linesplit[0] =~ /^y$/) {$linesplit[0] = "Y";}
	if ($linesplit[0] =~ /^MT$/) {$linesplit[0] = "M";}
	if ($linesplit[0] =~ /^m$/) {$linesplit[0] = "M";}
	$linesplit[0] =~ s/\.\d+$//;
	$linesplit[1] =~ s/\.\d+$//;
	if ($linesplit[0] !~ /^chr/) {$linesplit[0] = "chr".$linesplit[0];}
	
	my $chr = $linesplit[0];
	$chr =~ s/chr//;
	my $pos_1 = $linesplit[1];
	my $pos_2 = $pos_1;
	my $ref = $linesplit[3];
	my $var = $linesplit[4];
	my $format = $linesplit[8];
	if ($chr =~ /^M$/) {$chr = "MT";}
	
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
		if ($Sample_Call_Processed =~ m/1\// || $Sample_Call_Processed =~ m/\/1/) { # if the first variant is still the interest of the genotype call in any sample
			print OUT "$chr\t$pos_1_now\t$pos_2_now\t$ref_now\t$v_now\t$linesplit[0]\t$linesplit[1]\t$linesplit[2]\t".
           	"$linesplit[3]\t$out_var\t$linesplit[5]\t$linesplit[6]\t$linesplit[7]\t$linesplit[8]\t$Sample_Call_Processed\n";
		}
	}
}
close OUT;
close INPUT;

exit;

sub AddGenoTypeToSampleCalls_CompondHet {
	my ($first_V, $variants, $format, $sample_call) = @_;
	#print "$first_V, $variants, $format, $sample_call\n";
	my $gene_type_call_qual = 13; ##### genotype call quality cut off
	my @format_fields = split(":", $format);
	my $HT_index;
	#my $GQ_index;
	my $AD_index;
	for (my $i=0;$i<scalar @format_fields ;$i++){
		if ($format_fields[$i] =~ m/HT/) {
			$HT_index = $i;
		}
		if ($format_fields[$i] =~ m/AD/) {
			$AD_index = $i;
		}		
	}
	my $sample_call_processed=""; # vcf_sample1 \t geno_sample1 \t vcf_sample2 \t geno_sample2 ....
	my @sample_call_split = split("\t", $sample_call);
	for(my $i=0;$i<scalar(@sample_call_split);$i++){
		my $sample = $sample_call_split[$i];
		if ($sample !~ m/\.\/\./) {
			my @fields = split(":", $sample);
			my $HT = $fields[$HT_index];
			my $AD;
			if(defined $AD_index) {
				$AD = $fields[$AD_index];
			}
			
			if( $first_V eq $variants) {
				$sample_call_processed = $sample_call_processed."\t".$sample;
			} else {
				#determine which alternative allele is moved to the front
				#print "$sample\n";
				my @V = split("," , $variants);
				my @HT_array;
				
				@HT_array = split('/', $HT);
				
				my $V_index;
				my $the_first_V_AD_index; 
				for(my $h=0;$h<scalar @V;$h++) {
					if($first_V eq $V[$h]) {
						$V_index = $h+1;
						$the_first_V_AD_index = $V_index;
					}
				}
				#print "new first V index: $the_first_V_AD_index\n"; 
				my @new_HT_array;
				foreach my $ht_number (@HT_array) {
					#print "here is $ht_number: ";
					if ($ht_number == 0) {
						#print "0\n";
						push @new_HT_array, 0;
					} elsif ($ht_number < $the_first_V_AD_index) {
						#print "+1\n";
						push @new_HT_array, $ht_number+1;
					} elsif ($ht_number > $the_first_V_AD_index) {
						#print "no change\n";
						push @new_HT_array, $ht_number;
					} else {
						#print "to 1\n";
						push @new_HT_array, 1;
					}
				}	
				
				$sample="";
				for (my $i=0; $i < scalar @fields; $i++) {
					if ($i != $HT_index) {
						$sample=$sample.":".$fields[$i];
					} else {
						my $temp_1="";
						foreach my $temp_2 (@new_HT_array) {
							$temp_1=$temp_1."/".$temp_2;
						}
						$temp_1 =~ s/^\///;
						$sample=$sample.":".$temp_1;
					}
				}
				$sample =~ s/^\://;
				
				$sample_call_processed = $sample_call_processed."\t".$sample;
			}
		} else {
			$sample_call_processed = $sample_call_processed."\t".$sample;
		}
	}
	$sample_call_processed =~ s/^\t//;
	return $sample_call_processed;
}


