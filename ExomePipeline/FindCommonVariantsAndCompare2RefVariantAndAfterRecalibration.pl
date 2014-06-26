#!/usr/bin/perl
use strict;
use warnings;

my $folder = "/users/a5907529/lustre/Kidney_20130211/A1569_Fastq/SensSpec/Compare_B4_AF_Recali";
my $ref_var_file="/users/a5907529/lustre/scripts/SensSpecScript/VariantList_hg19.txt";
my $af_recali_vcf="/users/a5907529/lustre/Kidney_20130211/A1569_Fastq/UnifiedtypeCaller_All/unifiedgenotyper_recali_SNP_INDEL.vcf";
my $output = "/users/a5907529/lustre/Kidney_20130211/A1569_Fastq/SensSpec/Compare_B4_AF_Recali/Summary.txt";

my %All_Variants;
my %Hetereo_count;
my %Homo_count;
my %Ref_count;

opendir (DIR, $folder) or die "Cann't find the directory";
while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./);
	if ($file =~ m/\.onRef.vcf$/){ ## This line controls the name pattern of files that are read in
		my $input_file = $folder."/".$file;
		open INPUT, $input_file or die "cannot open file $input_file.\n";
		print "reading ".$input_file."\n";
		while (my $line =<INPUT>) {
			if($line =~ m/^\#/) {
				next;
			} else {
				chomp $line;
				my @elements  = split ("\t", $line);
				my $chr=$elements[0];
				my $pos=$elements[1];
				my $R=$elements[3];
				my $A=$elements[4];
				if (exists $All_Variants{$chr}{$pos}{$R}{$A} ) {
					$All_Variants{$chr}{$pos}{$R}{$A} = $All_Variants{$chr}{$pos}{$R}{$A} + 1;
				} else {
					$All_Variants{$chr}{$pos}{$R}{$A} = 0
				}
				my $FORMAT=$elements[8];
				my $Sample_Call=$elements[9];
				my $GT = AddGenoTypeToSampleCalls($FORMAT, $Sample_Call);
				if ($GT eq "R/R") {
					if (exists $Ref_count{$chr}{$pos}){
						$Ref_count{$chr}{$pos} = $Ref_count{$chr}{$pos}  + 1;
					} else {
						$Ref_count{$chr}{$pos} = 1;
					}
				} elsif ($GT eq "R/V") {
					if (exists $Hetereo_count{$chr}{$pos}){
						$Hetereo_count{$chr}{$pos} = $Hetereo_count{$chr}{$pos}  + 1;
					} else {
						$Hetereo_count{$chr}{$pos} = 1;
					}
				} elsif ($GT eq "V/V") {
					if (exists $Homo_count{$chr}{$pos}){
						$Homo_count{$chr}{$pos} = $Homo_count{$chr}{$pos}  + 1;
					} else {
						$Homo_count{$chr}{$pos} = 1;
					}
				} else {
					print "It's wrong wrong wrong.....\n";
				}
			}
		}
		close INPUT;
	}
}
close DIR;

my %Filter_info;
open RECALI_VCF, "$af_recali_vcf" or die "Cannot open the file.\n";
while (my $line = <RECALI_VCF>) {
	if($line =~ m/^\#/) {
		next;
	} else {
		chomp $line;
		my @elements  = split ("\t", $line);
		my $chr=$elements[0];
		my $pos=$elements[1];
		my $filter=$elements[6];
		if (exists $All_Variants{$chr} && exists $All_Variants{$chr}{$pos}) {
				$Filter_info{$chr}{$pos} = $filter;
		}
	}
}
close  RECALI_VCF;

my %Ref_Var_Info;
open REF_VAR, "$ref_var_file" or die "Cannot open the file $ref_var_file\n";
while (my $line = <REF_VAR>) {
	if($line =~ m/^\#/) {
		next;
	} else {
		chomp $line;
		my @elements  = split ("\t", $line);
		my $chr=$elements[2];
		my $pos=$elements[3];
		my $change=$elements[1];
		my $freq=$elements[4];
		my $rs=$elements[0];
		if (exists $All_Variants{$chr} && exists $All_Variants{$chr}{$pos}) {
			$Ref_Var_Info{$chr}{$pos} = "$change\t$freq\t$rs";
		}
	}
}
close REF_VAR;

open OUTPUT, ">$output" or die "Cannot open the file $output\n";
print OUTPUT "Chr\tPos\tRef\tVar\tRef_Change\tRef_freq\tDB_snp\tFilter_Info\tRR_count\tRV_Count\tVV_Count\n";
foreach my $chr (keys %All_Variants) {
	foreach my $pos (sort {$a <=> $b} keys %{$All_Variants{$chr}} ) {
		foreach my $ref (keys %{$All_Variants{$chr}{$pos}} ) {
			foreach my $var (keys %{$All_Variants{$chr}{$pos}{$ref}} ) {
				my $RR_count = "0";
				my $RV_count = "0";
				my $VV_count = "0";
				if (exists $Ref_count{$chr}{$pos}) { 
					$RR_count = $Ref_count{$chr}{$pos}; 
				}
				if (exists $Hetereo_count{$chr}{$pos}) {
					$RV_count = $Hetereo_count{$chr}{$pos};
				}
				if (exists $Homo_count{$chr}{$pos}) {
					$VV_count = $Homo_count{$chr}{$pos};
				}
				print OUTPUT "$chr\t$pos\t$ref\t$var\t$Ref_Var_Info{$chr}{$pos}\t$Filter_info{$chr}{$pos}\t$RR_count\t$RV_count\t$VV_count\n";
			}
		}
	}
}
close OUTPUT;

exit;

sub AddGenoTypeToSampleCalls {
	my ($format, $sample_call) = @_;
	my $gene_type_call_qual = 13;
	my $GT;
	
	my @format_fields = split(":", $format);
	my $GT_index;
	my $GQ_index;
	for (my $i=0;$i<scalar(@format_fields);$i++){
		if ($format_fields[$i] =~ m/GT/) {
			$GT_index = $i;
		}
		if ($format_fields[$i] =~ m/GQ/) {
			$GQ_index = $i;
		}
	}
	
	if ($sample_call !~ m/\.\/\./) {
		my @fields = split(":", $sample_call);
		$GT = $fields[$GT_index];
		my $GQ = $fields[$GQ_index];
		if ($GQ < $gene_type_call_qual) {
			$GT = "R/R";
		} elsif ($GT =~ m/0\/0/) {
			$GT = "R/R";
		} elsif ($GT =~ m/0\/[123456789]/) {
			$GT = "R/V";
		} elsif ($GT =~ m/[123456789]\/[123456789]/) {
			$GT = "V/V";
		}
	}
	return $GT;
}
			

