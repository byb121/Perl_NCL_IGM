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
				for(my $i=0;$i<scalar @var_chars;$i++) {
					
				}
				$ref_now="-";
				for (my $i=1;$i<=length($ref);$i++){
					$v_now=~s/^.//;
					$pos_1_now++;
				}
				$pos_2_now = $pos_1_now;
			} elsif(length($v_now)==length($ref_now)) { #in fact a single base SNP on the first base of ref
				$ref_now = substr $ref,0,1;
				$v_now = substr $v,0,1;
			} else { #deletion
				for(my $i=1;$i <= length($v_now);$i++) {
					$ref_now=~s/^.//;
					$pos_1_now++;
				}
				$v_now="-";
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