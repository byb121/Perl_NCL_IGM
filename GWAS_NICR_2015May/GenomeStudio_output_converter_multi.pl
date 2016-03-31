#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $output_prefix;
my $genotype_files; #can be a list of GenomeStudio output genotype files, need to be separated with ",".
my $SNP_files;  #can be a list of GenomeStudio output SNP info files, need to be separated with ",".
my $sample_files;  #can be a list of GenomeStudio output Sample info files, need to be separated with ",".
my $gender_file; #file to provide gender information if it is not available in sample file, otherwise will take the gender column in the Sample file.
my $sample_list_file; #optional file to choose an subset of samples, only the first column (first line ignored) will be taken as the sample IDs.
my $phenotype = 1;
my $help;

my $datestring;

if ( @ARGV < 4 || ! GetOptions('help|?' => \$help, "geno=s"=>\$genotype_files, "snp=s"=>\$SNP_files, 'sample=s' => \$sample_files, 
				'sex=s' => \$gender_file, 's_list=s' => \$sample_list_file, 'pheno=s' => \$phenotype, 'out=s' => \$output_prefix) 
				|| defined $help ) {
	usage(); exit;
}

#check file existence
$datestring = localtime();
print STDOUT "[$datestring]: Starting...\n";
my @genotype_file_array = split(",", $genotype_files);
foreach my $genotype_file (@genotype_file_array) {
	if (-e $genotype_file) {
		print STDOUT "-- Found $genotype_file\n";
	} else {
		print STDOUT "$genotype_file does not exist\n exit\n";
		exit;
	}
}

my @snp_file_array = split(",", $SNP_files);
foreach my $SNP_file (@snp_file_array) {
	if (-e $SNP_file) {
		print STDOUT "-- Found $SNP_file\n";
	} else {
		print STDOUT "$SNP_file does not exist\n exit\n";
		exit;
	}
}

my @sample_file_array = split(",", $sample_files);
foreach my $sample_file (@sample_file_array) {
	if (-e $sample_file) {
		print STDOUT "-- Found $sample_file\n";
	} else {
		print STDOUT "$sample_file does not exist\n exit\n";
		exit;
	}
}

my %sex;
if (defined $gender_file) {
	if (-e $gender_file) {
		print STDOUT "-- Found $gender_file\n";
		print STDOUT "-- Reading in gender information....\n";
		open SEX, "$gender_file" or die "Cannot open the file $gender_file.\n";
		my $line_count = 1;
		while (my $line=<SEX>) {
			if ($line_count == 1 ) {
				$line_count += 1;
				next;
			} else {
				chomp $line;
				my @temp_line = split("\t", $line);
				$temp_line[1] =~ s/\W//g;
				if (! exists $sex{$temp_line[0]}) {
					$sex{$temp_line[0]} = $temp_line[1];
				} else {
					print STDOUT "Error (030): Duplicated sample ID: $temp_line[0] is found in the gender file $gender_file\nexit!\n";
					exit;
				}
				$line_count += 1;
			}
		}
	} else {
		print STDOUT "$gender_file does not exist\n exit\n";
		exit;
	}
}

#search ".GType" column in genotype file to get sample ID and count sample column number and row number
$datestring = localtime();
print STDOUT "[$datestring]: Start to count rows and columns...\n";
my %sample_ID_index;
my $genotype_row=0;
my $genotype_geno_column=0;
open GENO, "$genotype_file" or die "Can not open the file $genotype_file.\n";
while (my $line=<GENO>) {
	$genotype_row += 1;
	if ($genotype_row == 1) {
		chomp $line; #print $line."\n";
		my @temp_split = split("\t", $line);
		foreach (my $i=0;$i<scalar @temp_split;$i++) {
			#print $temp_split[$i]."\n";
			if ($temp_split[$i] =~ m/\.GType.*/) {
				$genotype_geno_column += 1;
				my $temp_word = $temp_split[$i];  $temp_word =~ s/\.GType.*//; #print ":$genotype_geno_column:\t:$temp_word;;;;"."\n";
				if (!exists $sample_ID_index{$temp_word}) {
					$sample_ID_index{$temp_word} = $i;
				} else {
					print STDOUT "Error (001): Duplicated sample ID is found: $temp_word\n exit!\n";
					exit;
				}
				if (defined $gender_file) {
					if (!exists $sex{$temp_word}) {
						print STDOUT "Error (040): Sample ID: $temp_word is not found in the gender file provided\n exit!\n";
						exit;
					}
				}
			}
		}
	}
}
close GENO;

#count snp table row number
my $snp_row=0;
open SNP, "$SNP_file" or die "Can not open the file $SNP_file.\n";
while (my $line=<SNP>) {
	$snp_row += 1;
}
close SNP;

#count sample table row number
my $sample_row=0;
open SAM, "$sample_file" or die "Can not open the file $sample_file.\n";
while (my $line=<SAM>) {
	$sample_row += 1;
}
close SAM;

#genotype table column number should be the same as sample file row numbers
if($genotype_row != $snp_row) {
	print STDOUT "Error (002): Number of genotype file rows is not equal to number of SNP file rows.\nexit!\n";
	exit;
}
#SNP table row need to be the same as genotypes file row
#print "$sample_row\t$genotype_geno_column\t$snp_row\n";
if($genotype_geno_column != $sample_row-1) {
	print STDOUT "Error (003): Number of genotype file sample columns is not equal to number of sample file rows.\nexit!\n";
	exit;
}


my %geno;
# a sample list could be supplied
my @sample_list;
# a snp list could be supplied
my @snp_list;
#read in sample file;
$datestring = localtime();
print STDOUT "[$datestring]: Start to read sample info ...\n";
if (! @sample_list) {
	open SAM, "$sample_file" or die "Can not open the file $sample_file.\n";
	my $line_count = 1;
	my $id_index;
	my $sex_index;
	while (my $line=<SAM>) {
		if ($line_count == 1) { #locate sample ID column
			my @temp_line = split("\t", $line);
			for (my $i=0;$i<scalar @temp_line;$i++) {
				if($temp_line[$i] eq "Sample ID") {
					$id_index = $i;
				} elsif ($temp_line[$i] eq "Gender") {
					$sex_index = $i;
				}
			}
			unless (defined $id_index) {
				print STDOUT "Error (004): \"Sample ID\" column was not found in the sample file $sample_file.\nexit!\n";
				exit;
			}
			if (! defined $gender_file) {
				unless (defined $sex_index) {
					print STDOUT "Error (005): \"Gender\" column was not found in the sample file $sample_file.\nexit!\n";
					exit;
				}
			}
		} else {
			#print $line;
			my @temp_line = split("\t", $line);
			#print "hulalall:$temp_line[$sex_index]:$temp_line[$id_index]:\n";
			if (! defined $gender_file) {
				if (! exists $sex{$temp_line[$id_index]}) {
					$sex{$temp_line[$id_index]} = $temp_line[$sex_index];
				} else {
					print STDOUT "Warning (001): ".$temp_line[$id_index]." found more than 1 time in the sample files.\nexit.\n";
				}
			}
			if (exists $sample_ID_index{$temp_line[$id_index]}) {
				push @sample_list, $temp_line[$id_index];
			} else {
				print STDOUT "Error (006): ".$temp_line[$id_index]." does not exists in genotype files.\nexit.\n";
				exit;
			}
		}
		$line_count += 1;
	}
	close SAM;
	print STDOUT "-- Found ".scalar @sample_list." samples from $sample_file\n";
}

foreach my $sample (@sample_list) {
	if (! exists $sex{$sample} ) {
		print STDOUT "Error (041): $sample is not found in the gender info.\n exit!\n";
		exit;
	}
}

$datestring = localtime();
print STDOUT "[$datestring]: Start to read snp info ...\n";
if (! @snp_list) {
	open SNP, "$SNP_file" or die "Can not open the file $SNP_file.\n";
	my $line_count = 1;
	my ($chr_index, $name_index, $pos_index, $snp_index);
	while (my $line=<SNP>) {
		if ($line_count == 1) { #locate sample ID column
			my @temp_line = split("\t", $line);
			for (my $i=0;$i<scalar @temp_line;$i++) {
				#print ";$temp_line[$i];;\n";
				if($temp_line[$i] eq "Chr") {
					$chr_index = $i;
					#print ";chr;\t$chr_index;\n";
				}
				if($temp_line[$i] eq "Name") {
					$name_index = $i;
				}
				if($temp_line[$i] eq "Position") {
					$pos_index = $i;
				}
				if($temp_line[$i] eq "SNP") {
					$snp_index = $i;
				}
			}
			unless (defined $chr_index) {
				print STDOUT "Error (007): \"Chr\" column was not found in the snp file $SNP_file.\nexit!\n";
				exit;
			}
			unless (defined $name_index) {
				print STDOUT "Error (008): \"Name\" column was not found in the snp file $SNP_file.\nexit!\n";
				exit;
			}
			unless (defined $pos_index) {
				print STDOUT "Error (009): \"Position\" column was not found in the snp file $SNP_file.\nexit!\n";
				exit;
			}
			unless (defined $snp_index) {
				print STDOUT "Error (010): \"SNP\" column was not found in the snp file $SNP_file.\nexit!\n";
				exit;
			}
		} else {
			my @temp_line = split("\t", $line);
			my @temp_line_SNP = split('/', $temp_line[$snp_index]);
			my $A = $temp_line_SNP[0]; $A =~ s/^\[//;
			my $B = $temp_line_SNP[1]; $B =~ s/\]$//;
			push @snp_list, $temp_line[$chr_index].'@'.$temp_line[$name_index].'@'.$temp_line[$pos_index].'@'.$A.'@'.$B;
		}
		$line_count += 1;
	}
	close SNP;
	print "-- Found ".scalar @snp_list." snps from $SNP_file\n";
}

$datestring = localtime();
print STDOUT "[$datestring]: Start to read in genotypes from $genotype_file ...\n";
open GENO, "$genotype_file" or die "Can not open the file $genotype_file.\n";
my $line_count = -1;
while (my $line=<GENO>) {
	if ($line_count == -1) { #locate sample ID column
		$line_count += 1;
		next;
	} else {
		my @temp_split = split("\t", $line);
		my @temp_snp_split = split('@', $snp_list[$line_count]);
		foreach my $sample (@sample_list) {
			$temp_split[$sample_ID_index{$sample}] =~ s/\W//g;
			if ($temp_split[$sample_ID_index{$sample}] eq "AA" ) {
				$geno{$temp_snp_split[0]}{$temp_snp_split[2]}{$temp_snp_split[1]}{$sample} = $temp_snp_split[3]." ".$temp_snp_split[3];
			} elsif ($temp_split[$sample_ID_index{$sample}] eq "AB" ) {
				$geno{$temp_snp_split[0]}{$temp_snp_split[2]}{$temp_snp_split[1]}{$sample} = $temp_snp_split[3]." ".$temp_snp_split[4];
			} elsif ($temp_split[$sample_ID_index{$sample}] eq "BB" ) {
				$geno{$temp_snp_split[0]}{$temp_snp_split[2]}{$temp_snp_split[1]}{$sample} = $temp_snp_split[4]." ".$temp_snp_split[4];
			} elsif ($temp_split[$sample_ID_index{$sample}] eq "NC" ) {
				$geno{$temp_snp_split[0]}{$temp_snp_split[2]}{$temp_snp_split[1]}{$sample} = "0 0";
			} else {
				print STDOUT "Error (020): Can not recognise the genotype ".$temp_split[$sample_ID_index{$sample}].".\nexit!\n";
			}
		}
		$line_count += 1;
	}
	if ($line_count%100000 == 0) {
		print STDOUT "---- processed $line_count lines.\n";
	}
}
close GENO;

#output
$datestring = localtime();
print STDOUT "[$datestring]: Start to produce output...\n";
my $tped_output = $output_prefix.".tped"; open TPED, ">$tped_output" or die "Can not open the file to output.\n";
my $tfam_output = $output_prefix.".tfam"; open TFAM, ">$tfam_output" or die "Can not open the file to output.\n";

foreach my $chr ( sort {$a cmp $b} keys %geno) {
	foreach my $pos ( sort {$a<=>$b} keys %{$geno{$chr}}) {
		foreach my $snp_name ( keys %{$geno{$chr}{$pos}}) {
			my $temp_string = "$chr\t$snp_name\t0\t$pos\t";
			foreach my $sample (@sample_list) {
				$temp_string = $temp_string.$geno{$chr}{$pos}{$snp_name}{$sample}."\t";
			}
			$temp_string =~ s/\t$/\n/;
			print TPED $temp_string;
		}
	}
}
close TPED;

foreach my $sample (@sample_list) {
	print TFAM "$sample\t$sample\t0\t0\t".$sex{$sample}."\t$phenotype\n";
}
close TFAM;
print STDOUT "All done!\n";

exit;

#subroutines

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "You have not supplied all parameters I need, please check the usage below.\n";
    print "\nusage: VCF_2_annotated_excel_20131120.pl \n";
    print "--vcf input vcf file (of a single sample or a family;\n";
    print "--InterestedGenes file of interested gene names list (optional); Format: Ensembl gene IDs on the 1st column.\n"; 
    print "--out output excel file name (optional, only when required);\n";
    print "--CNV CNV result file, output of HG's Annotate_CNVs_combine_multiple_files.pl, sample names must be consistent with the vcf.\n\n";
	return(1);
	exit;
}
