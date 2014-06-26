#!/usr/bin/perl
use strict;
use warnings;
#use DBI;
#use DBD::mysql;
use LWP::Simple;
use Getopt::Long;

print "\n";
print "#################################################################\n";
print "# Author is Yaobo. The script should be faster than SQL version #\n";
print "#################################################################\n";
print "\n";
print "\n";
print "\n";
print "\n";
my($input_file) = @ARGV;

my $output_file_no_matched = "$input_file.dbSNP135.noSQL.filtered.txt";
my $output_file_matched = "$input_file.dbSNP135.noSQL.Matches.txt";
my $final_output = "$input_file.dbSNP135.genotype.filtered.txt";

my $SNP_file = "/users/a5907529/GenomeData/hg19/snp135.txt";

print "Variants data will be taken from: $input_file.\n";
print "Filtered result will be in: $output_file_no_matched.\n";
print "Matched result will be in: $output_file_matched.\n";
print "The script will take the first 3 colums as chromosome_name, position and Allele1/Allele2.\n";
print "Comparing to SNPs recorded in $SNP_file. Make sure that the file exists."."\n";

#Comparing_2_DB_SNPfile($input_file, $SNP_file, $output_file_no_matched, $output_file_matched);	

my $filteredresult = DetermineIfRNAediting($output_file_no_matched);

print "#################################################\n";
print "########### output final results ################\n";
print "#################################################\n";

open FINAL, ">$final_output";
print FINAL @$filteredresult;
close FINAL;

exit;
	
sub Comparing_2_DB_SNPfile {
	my ($input_file, $SNP_file, $output_file_no_matched, $output_file_matched)=@_;
	
	print "\n";
	print "Input data is from: $input_file.\n";
	print "Variants in the file will be compared to SNP recorded in $SNP_file.\n";
	print "\n";
	
	my %query_variants;
	
	open VARIANTS, $input_file or die "Cannot open $input_file"; 
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			#if($het_hom=~/heterozygous/ and $SplitLine[8]>85){next A_loop;}
			#if($het_hom=~/homozygous/ and $SplitLine[8]<=85){next A_loop;}
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			my $Mut=$SplitLine[2];
			#my $ComplimentMutant=$SplitLine[5];
			#$ComplimentMutant=~tr/ACGT/TGCA/;
			if (exists $query_variants{$Chr."_".$Pos}) {
				next;
				print "Warning: Ignored a duplicated entry found for ".$Chr."_".$Pos." in query file: $input_file.\n"
			} else {
				$query_variants{$Chr."_".$Pos} = $Mut;	
			}
		}
	}
	close VARIANTS;
	
	my %dbsnp_matches;
	open SNP_DB, $SNP_file or die "Cannot open $SNP_file"; 
	while (my $line = <SNP_DB>){  
		chomp $line;
		my @SplitLine=split(/\t/, $line);
		my $Chr=$SplitLine[1];
		my $Pos=$SplitLine[2]+1;
		my $Mut=$SplitLine[9];
		my $rs=$SplitLine[4];
		my $ref_base=$SplitLine[8];
		my $strand=$SplitLine[6];
		 
		#my $ComplimentMutant=$SplitLine[5];
		#$ComplimentMutant=~tr/ACGT/TGCA/;
		
		if (exists $query_variants{$Chr."_".$Pos}) {
			$dbsnp_matches{$Chr."_".$Pos} = $rs."\t".$Mut."\t".$ref_base."\t".$strand;
		} else {
			next;
		}
	}
	close SNP_DB;
	
	my @matched;
	my @filtered;
	
	push @matched, "chrom"."\t"."position"."\t"."genotype"."\t"."mapped gene"."\t"."ensembl ID"."\t"."heterozygous in OA"."\t"."heterozygous in NOF"."\t"."total heterozygous"."\t";
	push @matched, "rs"."\t"."recorded mutation"."\t"."ucsc ref base"."\t"."strand"."\n";
	push @filtered, "chrom"."\t"."position"."\t"."genotype"."\t"."mapped gene"."\t"."ensembl ID"."\t"."heterozygous in OA"."\t"."heterozygous in NOF"."\t"."total heterozygous"."\n";
	
	open VARIANTS, $input_file or die "Cannot open $input_file"; 
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			#my $Mut=$SplitLine[2];
			#my $ComplimentMutant=$SplitLine[5];
			#$ComplimentMutant=~tr/ACGT/TGCA/;
			if (exists $dbsnp_matches{$Chr."_".$Pos}) {
				push @matched, $line."\t".$dbsnp_matches{$Chr."_".$Pos}."\n";
			} else {
				push @filtered, $line."\n";
			}
		}
	}
	close VARIANTS;
	
	print "###########################################\n";
	print "########### output results ################\n";
	print "###########################################\n";
	print "$output_file_matched\n";
	print "$output_file_no_matched\n";
	
	open MATCHED, ">$output_file_matched";
	print MATCHED @matched;
	close MATCHED;
	
	open FILTERED, ">$output_file_no_matched";
	print FILTERED @filtered;
	close FILTERED;
	
	#return($output_file_matched, $output_file_no_matched);
	
}

sub DetermineIfRNAediting {
	
	my ($input_file)=@_;
		
	print "\n";
	print "Input data is from: $input_file.\n";
	print "Variants in the file will be filtered for RNA-eidting.\n";
	print "Possible RNA-editing events is define in the script as following:.\n";
	print "A-to-I: A/T, A/C, T/A, T/G.\n";
	print "C-to-U: C/T, G/A.\n";
	print "\n";
	
	my @filtered;
	
	open VARIANTS, $input_file or die "Cannot open $input_file"; 
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			my $Mut=$SplitLine[2];
			my $ComplimentMutant=$Mut;
			$ComplimentMutant=~tr/ACGT/TGCA/;
			
			#Determine if the Mut type is an possible RNA-editing event.
			if ($Mut =~ m/^[Aa]\/[Tt]|^[Aa]\/[Cc]|^[Tt]\/[Aa]|^[Tt]\/[Gg]|^[Cc]\/[Tt]|^[Gg]\/[Aa]/ | $ComplimentMutant =~ m/^[Aa]\/[Tt]|^[Aa]\/[Cc]|^[Tt]\/[Aa]|^[Tt]\/[Gg]|^[Cc]\/[Tt]|^[Gg]\/[Aa]/) {
				push @filtered, $line."\n";
			}
		}
	}
	close VARIANTS;
	
	return \@filtered;
}

	