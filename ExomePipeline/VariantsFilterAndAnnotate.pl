#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#input are vcfs
# will looking for same genotype in 'in' samples but different from ANY of the 'difffrom' samples

my $output="Filtered_v";
my $in;
my $difffrom;
my $help;

usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, "in=s"=>\$in, "difffrom=s"=>\$difffrom, 
, 'output=s' => \$output) || defined $help );

unless (defined $in && defined $difffrom && defined $output) {
	die "You have not supplied all parameters\n\n";
}

my @in_vcfs=split(",", $in);
my @in_sample_names;
foreach my $vcf (@in_vcfs) {
	my $bash_output=`grep '#CHR' $vcf`;
	my @bash_output_split = split("\t", $bash_output);
	print "Looking for sample names:\n";
	for(my $i=9;$i<scalar @bash_output_split;$i++){
		print "Found: ".$bash_output_split[$i]."\n";
		push @in_sample_names, $bash_output_split[$i];
	}
}

my @difffrom_vcfs=split(",", $difffrom);
my @diff_sample_names;
foreach my $vcf (@difffrom_vcfs) {
	my $bash_output=`grep '#CHR' $vcf`;
	my @bash_output_split = split("\t", $bash_output);
	print "Looking for sample names:\n";
	for(my $i=9;$i<scalar @bash_output_split;$i++){
		print "Found: ".$bash_output_split[$i]."\n";
		push @diff_sample_names, $bash_output_split[$i];
	}
}

#combine vcfs
#`module load apps/gatk/2.7.2/noarch`;

my $combined_vcf = "temp_combined.vcf";
my $combine_command = 'java -jar $GATKDIR/GenomeAnalysisTK.jar -T CombineVariants -R /users/a5907529/data/GATK/bundle2.5/hg19/ucsc.hg19.fasta'; 
$combine_command = $combine_command.' -genotypeMergeOptions REQUIRE_UNIQUE';
$combine_command = $combine_command." -o $combined_vcf";
foreach my $vcf (@in_vcfs) { 
	$combine_command = $combine_command." --variant $vcf";
}
foreach my $vcf (@difffrom_vcfs) { 
	$combine_command = $combine_command." --variant $vcf";
}
`$combine_command`;

#select variants
my @output;
open VCF, "$combined_vcf" or die "$combined_vcf does not exist or GATK failed to combine vcfs\n";
LOOP_VCF: while (my $line = <VCF> ) {
	chomp $line;
	if($line =~ m/^\#/) {
		push @output, $line."\n";
		next;
	} else {
		my @elements  = split ("\t", $line);
		my $format=$elements[8];
		my $GT_index;
		my $GQ_index;
		my @format_fields = split(":", $format);
		my $Sample_Call="";
		for(my $i=9;$i<scalar(@elements);$i++) {
			$Sample_Call=$Sample_Call."\t".$elements[$i];
		}
		$Sample_Call =~ s/^\t//;
		my $genotype_string = AddGenoTypeToSampleCalls($format, $Sample_Call);
		my @genotypes = split("\t", $genotype_string);
		#start comparison
		if (scalar @in_sample_names > 1) {
			for (my $i=1;$i<scalar @in_sample_names;$i++) {
				if ($genotypes[$i] ne $genotypes[0]) {
					next LOOP_VCF;
				}
			}
		}
		
		if (scalar @diff_sample_names > 1) {
			for (my $i=scalar @in_sample_names;$i<scalar @genotypes;$i++) {
				if ($genotypes[$i] eq $genotypes[0]) {
					next LOOP_VCF;					
				} 
			}
		} else {
			if ($genotypes[scalar(@genotypes)-1] eq $genotypes[0]) {
					next LOOP_VCF;					
			}
		}
		
		push @output, $line."\n";
	}
}

close (VCF);

my $output_vcf = $output.".vcf";
my $output_xlsx = $output.".xlsx";

open OUTPUT, ">$output_vcf" or die "Cannot open file $output_vcf to output. \n";
print OUTPUT @output;
close OUTPUT;

`rm $combined_vcf*`;

# To annotated filtered results
`perl /users/data/Files_HG/vcf_annotation_november2013/VCF_2_annotated_xlsx_20131204.pl --vcf $output_vcf --out $output_xlsx`;

exit;

sub AddGenoTypeToSampleCalls {
	my ($format, $sample_call) = @_;
	my $gene_type_call_qual = 13; ##### genotype call quality cut off
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
	my $sample_call_processed=""; # vcf_sample1 \t geno_sample1 \t vcf_sample2 \t geno_sample2 ....
	my @sample_call_split = split("\t", $sample_call);
	for(my $i=0;$i<scalar(@sample_call_split);$i++){
		my $sample = $sample_call_split[$i];
		if ($sample !~ m/\.\/\./) {
			my @fields = split(":", $sample);
			my $GT = $fields[$GT_index];
			my $GQ = $fields[$GQ_index];
			if ($GQ < $gene_type_call_qual) {
				#$GT = "R/R";
				$GT = "NA"; #if low quality make a null call "NA"
			} elsif ($GT =~ m/0\/0/) {
				$GT = "R/R";
			} elsif ($GT =~ m/0\/[123456789]/) {
				$GT = "R/V";
			} elsif ($GT =~ m/([123456789])\/([123456789])/) {
				my $left = $1;
				my $right = $2;
				if($left!=$right) {
					$GT = "V1/V2"; 
				} else {
					$GT = "V/V";
				}
			}
			$sample_call_processed = $sample_call_processed."\t".$GT;
		} else {
			$sample_call_processed = $sample_call_processed."\t"."NA";
		}
	}
	$sample_call_processed =~ s/^\t//;
	return $sample_call_processed;
}

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: VariantsFilterAndAnnotate.pl \n";
    print "--in vcfs saparated with ',' looking for common variants in these vcfs.\n";
    print "--difffrom vcfs saparated with ',' looking for variants is different from ANY of the variants in these files.\n"; 
    print "--output file name of output. vcf file has the filtered variants, XLSX file has the annotated results.\n";
    return(1);
}