#!/usr/bin/perl
use strict;
use warnings;

my($input_file) = @ARGV;
my $SNP_file = "/users/a5907529/GenomeData/hg19/snp135.txt";

my $snp_output = $input_file."_SNP135_tested.txt";
my $editing_output = $input_file."_Editing_tested.txt";

Comparing_2_DB_SNPfile($input_file, $snp_output, $SNP_file);
DetermineIfRNAediting($snp_output, $editing_output);

exit;

sub Comparing_2_DB_SNPfile {
	my ($input_file, $output, $SNP_file)=@_;
	
	print "\n";
	print "Input data is from: $input_file.\n";
	print "Variants in the file will be compared to SNP recorded in $SNP_file.\n";
	print "An extra column will be added to the file to indicate if it's a recored SNP.\n";
	print "Column values includes: 'Recorded' and 'Noval'.\n";
	print "\n";
	
	my %dbsnp;
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
		
		if (!exists $dbsnp{$Chr."_".$Pos}) {
			$dbsnp{$Chr."_".$Pos} = 1;
		} else {
			$dbsnp{$Chr."_".$Pos} += 1;
		}
	}
	close SNP_DB;
	

	my @output;
	



		
		
	open VARIANTS, $input_file or die "Cannot open $input_file";
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^Chromosome/) {
			push @output, $line."\t"."SNP135";
		} elsif ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			#if($het_hom=~/heterozygous/ and $SplitLine[8]>85){next A_loop;}
			#if($het_hom=~/homozygous/ and $SplitLine[8]<=85){next A_loop;}
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			my $Ref=$SplitLine[3];
			my $Var=$SplitLine[4];
			my $Mut=$Ref."/".$Var;
			#my $ComplimentMutant=$SplitLine[5];
			#$ComplimentMutant=~tr/ACGT/TGCA/;
			if (exists $dbsnp{$Chr."_".$Pos}) {
				push @output, $line."\t"."Recorded";
			} else {
				push @output, $line."\t"."Novel";
			}
		}
	}
	close VARIANTS;
	
	open OUTPUT, ">$output" or die "Cannot open the file $output";
	foreach my $line (@output) {
		print OUTPUT $line."\n";
	}
	close OUTPUT;

	print "###########################################\n";
	print "########### output results ################\n";
	print "###########################################\n";
	print "$output\n";

	#return($output_file_matched, $output_file_no_matched);
	
}

sub DetermineIfRNAediting {
	
	my ($input_file, $output)=@_;
		
	print "\n";
	print "Input data is from: $input_file.\n";
	print "Variants in the file will be tested for RNA-eidting.\n";
	print "Possible RNA-editing events is define in the script as following:.\n";
	print "A-to-I: A/T, A/C, T/A, T/G.\n";
	print "C-to-U: C/T, G/A.\n";
	print "An extra column will be added to the file to indicate if it's a Editing event.\n";
	print "Column Values are: Y - yes; N - no.\n";
	print "\n";
	
	my @output;
	
	open VARIANTS, $input_file or die "Cannot open $input_file"; 
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^Chromosome/) {
			push @output, $line."\t"."RNA_Editing";
		} elsif ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			my $Ref=$SplitLine[3];
			my $Var=$SplitLine[4];
			my $Mut=$Ref."/".$Var;
			my $ComplimentMutant=$Mut;
			$ComplimentMutant=~tr/ACGT/TGCA/;
			
			#Determine if the Mut type is an possible RNA-editing event.
			if ($Mut =~ m/^[Aa]\/[Tt]|^[Aa]\/[Cc]|^[Tt]\/[Aa]|^[Tt]\/[Gg]|^[Cc]\/[Tt]|^[Gg]\/[Aa]/ | $ComplimentMutant =~ m/^[Aa]\/[Tt]|^[Aa]\/[Cc]|^[Tt]\/[Aa]|^[Tt]\/[Gg]|^[Cc]\/[Tt]|^[Gg]\/[Aa]/) {
				push @output, $line."\t"."Y";
			} else {
				push @output, $line."\t"."N";
			}
		}
	}
	close VARIANTS;
	
	open OUTPUT, ">$output" or die "Cannot open the file $output";
	foreach my $line (@output) {
		print OUTPUT $line."\n";
	}
	close OUTPUT;
}