#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use LWP::Simple;
use Getopt::Long;


print "\n";
print "############################################################\n";
print "# Author is not Yaobo. The script was adapted from Helen's #\n";
print "############################################################\n";
print "\n";
print "\n";

my($input_file) = @ARGV;
my $output_file_no_matched = "$input_file.dbSNP135.filtered.txt";
my $output_file_matched = "$input_file.dbSNP135.Matches.txt";

print "Variants data will be taken from: $input_file.\n";
print "Filtered result will be in: $output_file_no_matched.\n";
print "Matched result will be in: $output_file_matched.\n";
print "The script will take the first 3 colums as chromosome_name, position and Allele1/Allele2.\n";

#CALL: /users/nhrg/AnalysisScripts/VarScan2PolyPhenWithCCDS_hg19.pl --geno "homozygous"  --pat "patientname" --lane "lane_no." --file "snps_on_target.txt"
#my $patient="Patient";
#my $lane="Lane";
#my $vs_file="vs_file"; #probableSNPs file

#my $het_hom="all"; #SNP type
#my $Results=GetOptions("geno=s"=>\$het_hom, "pat=s"=>\$patient, "lane=s"=>\$lane, "file=s"=>\$vs_file);

#my $dbSNPmatchfile = MySQLforVariants($CurrentDir, $vs_file, $patient, $lane, $het_hom);

my ($dbSNPmatchfile, $filteredResult)  = Comparing2dbSNPs($input_file, $output_file_no_matched, $output_file_matched);
                        
#my $dbSNPmatchfile = "AlldbSNPmatches_Lane4_CuB_new_filter_all.txt";                                    

#FindRareVariants($CurrentDir, $dbSNPmatchfile, $vs_file, $patient, $lane, $het_hom);

#my $plhh=$patient."_".$lane."_".$het_hom;
#`rm ReformattedVariantFile_$plhh.txt`;

exit;

#################### SUBROUTINES ##########################

sub Comparing2dbSNPs{
	#my ($CurrentDir, $probableSNPs, $patient, $lane, $het_hom)=@_;
	#my $IDplh=$patient."_".$lane."_".$het_hom;
	my ($input_file, $output_file_no_matched, $output_file_matched)=@_;
	
	print "\n";
	print "Warning: No ".'","'. " is allowed in the input file: $input_file, it will cause error in MySQL query.\n";
	print "\n";
	
	my @OutputFile;
	
	open ALLSNPS, $input_file or die "Cannot open $input_file"; 
	#my $header=<ALLSNPS>;
	
	A_loop: while (my $line = <ALLSNPS>){  
		chomp $line;
		if ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			#if($het_hom=~/heterozygous/ and $SplitLine[8]>85){next A_loop;}
			#if($het_hom=~/homozygous/ and $SplitLine[8]<=85){next A_loop;}
			my $Chr=$SplitLine[0]; 
			my $Pos=$SplitLine[1]-1; #dbSNP is 0-based. VCF file is 1-based
			my $Mut=$SplitLine[2];
			#my $ComplimentMutant=$SplitLine[5];
			#$ComplimentMutant=~tr/ACGT/TGCA/;
			push @OutputFile, $Chr.",".$Pos.",".$Mut.",".$line."\n";#only forward strand mut listed in file, no (rev strand) complement listed
		}
	}
	close ALLSNPS;
	
	my $ReformattedOutput = "$input_file"."_Comparing2dbSNP135.temp.txt";
	open REFORMAT, ">$ReformattedOutput";
	print REFORMAT @OutputFile;
	close REFORMAT;
	
	print "intermidiate file is generated with reformated SNPs info from $input_file.\n";
	
	
	##MySQL comparison.
	
	#Connect the MySQL DB
	print"Connecting the database....\n";
	my $ds="DBI:mysql:Data:headnode";
	my $user="nhrg"; 	# Connect as the user
	my $dbh=DBI->connect($ds,$user)||die "Cant connect to MySQL";
	print"Connection is made.\n";
	
	
	print"Prepaing tables......";
	
	# create a table to contain  variants/SNP included in the input file - reformated in the previous block.
	# Format: col1: chromosome name             col2: position              col3: ref_base/var_base            col4: the whole line of the input data    
	my $LocatedSNPsTable = $dbh -> prepare("CREATE TABLE LocatedSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(25), col4 VARCHAR(1000))");
	$LocatedSNPsTable -> execute();
	print "Table for input data is created!\n";
	
	print "Loading data to the table ........";
	# To load the reformated input file data into the table just created.
	my $LoadLocatedSNPsTable = $dbh -> prepare("LOAD DATA LOCAL INFILE '$ReformattedOutput' INTO TABLE LocatedSNPs FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n'");
	$LoadLocatedSNPsTable -> execute();
	print "Successfully loaded.\n";
	
	print "Indexing the first 4 columns...........";
	#index table to prepare the query
	my $Index1 = $dbh -> prepare("CREATE INDEX col_index1 ON LocatedSNPs (col1)");
	my $Index2 = $dbh -> prepare("CREATE INDEX col_index2 ON LocatedSNPs (col2)");
	my $Index3 = $dbh -> prepare("CREATE INDEX col_index3 ON LocatedSNPs (col3)");
	my $Index4 = $dbh -> prepare("CREATE INDEX col_index3 ON LocatedSNPs (col4)");
	
	$Index1 -> execute();
	$Index2 -> execute();
	$Index3 -> execute();
	$Index4 -> execute();
	print "Done!\n";
	
	print"Comparing input data to dbSNP 135.............................";
	my $AllMatchingSNPsTable = $dbh->prepare("CREATE TABLE AllMatchingSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(25), col4 CHAR(15), col5 CHAR(10), col6 CHAR(54), col7 CHAR(10))");
	$AllMatchingSNPsTable -> execute();
	
	my $LoadMatchingSNPsTable1 = $dbh -> prepare(
	"INSERT INTO AllMatchingSNPs
	SELECT LocatedSNPs.col1,LocatedSNPs.col2,LocatedSNPs.col3,snp135.name,snp135.observed,snp135.alleleFreqs,snp135.alleles
	FROM snp135, LocatedSNPs WHERE snp135.chrom=LocatedSNPs.col1 AND snp135.chromStart=LocatedSNPs.col2"); #output all chr pos matches 
	$LoadMatchingSNPsTable1 -> execute();
	print"Done!\n";
	
	print "Filtering input data to remove dbSNP 135 records........";
	# select variants/SNPs which are contained in the table with reformated data and dbSNP135 based on the location only.
	#output matched/overlaped hits.
	my $RemoveMatchedSNPs = $dbh -> prepare("DELETE FROM LocatedSNPs
	WHERE
	EXISTS(SELECT chrom FROM snp135
		WHERE snp135.chrom=LocatedSNPs.col1 AND snp135.chromStart=LocatedSNPs.col2)"); #remove rows which have matches in dbSNP135 
	$RemoveMatchedSNPs -> execute();
	print"Done!\n";
	
	print"Generating result.......";
	
	my $AllMatchingSNPsStatement = "SELECT col1,col2,col3,col4,col5,col6,col7 FROM AllMatchingSNPs";
	my $AllData = $dbh -> prepare($AllMatchingSNPsStatement);
	$AllData -> execute();
	my $AllDataContent = $AllData -> fetchall_arrayref();
	
	my $Filtered_Input = "SELECT col4 FROM LocatedSNPs";
	my $FilteredData = $dbh -> prepare($Filtered_Input);
	$FilteredData -> execute();
	my $FilteredDataContent = $FilteredData -> fetchall_arrayref();	

	print"Done!\n";
	
	
	
	print"Cleaning tables and intermidiate files......";
	#####Removing all the tables form MySQL.
	my $RemoveAllTable = $dbh -> prepare("DROP TABLE AllMatchingSNPs");
	$RemoveAllTable -> execute();
	my $RemoveLocatedTable = $dbh -> prepare("DROP TABLE LocatedSNPs");
	$RemoveLocatedTable -> execute();
	#`rm $ReformattedOutput`;
	print"Done!\n";
	
	
	my @Temp;
	my %RemoveDuplicates;
	my @RemoveDuplicatesArray;
	my $DuplicateCounter = 0;
	foreach my $Line(@$AllDataContent){
		push @Temp, "@$Line\n";
	}
	foreach my $TempLine(@Temp){
		if(!exists $RemoveDuplicates{$TempLine}){ ##There will be some duplicates as in matching to dbSNP, some of the reverse compliments are the same-such as C/G (Leads to G/C,C/G,GC).
		$RemoveDuplicates{$TempLine} = $DuplicateCounter;
		push @RemoveDuplicatesArray, $TempLine;
		$DuplicateCounter++;
		}
	}
	
	print "###########################################\n";
	print "########      output results     ##########\n";
	print "###########################################\n";
	
	print "$output_file_matched\n";
	my $AllMatchesOutput = $output_file_matched;
	open ALL, ">$AllMatchesOutput";
	print ALL @RemoveDuplicatesArray;
	#print ALL @Temp;
	close ALL;
	
	print "$output_file_no_matched\n";
	my $FilteredOutput = $output_file_no_matched;
	open FILTERED, ">$FilteredOutput";
	my @temp;
	foreach my $line(@$FilteredDataContent){
		push @temp, "@$line\n";
	}
	print FILTERED @temp; 
	close FILTERED;
	
	return($AllMatchesOutput, $FilteredOutput);
	
}

sub Comparing_2_DB_SNPfile {
	my ($input_file, $SNP_file, $output_file_no_matched, $output_file_matched)=@_;
	
	print "\n";
	print "Input data if from: $input_file.\n";
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
	open SNP_DB, $input_file or die "Cannot open $SNP_file"; 
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
	
	push @matched, "chrom"."\t"."position"."genotype"."\t"."mapped gene"."\t"."ensembl ID"."\t"."heterozygous in OA"."\t"."heterozygous in NOF"."\t"."total heterozygous"."\t";
	push @matched, "rs"."\t"."recorded mutation"."\t"."ucsc ref base"."\t"."strand"."\n";
	push @filtered, "chrom"."\t"."position"."genotype"."\t"."mapped gene"."\t"."ensembl ID"."\t"."heterozygous in OA"."\t"."heterozygous in NOF"."\t"."total heterozygous"."\n";
	
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
	
	return($output_file_matched, $output_file_no_matched);
	
}


sub FindRareVariants{
	###PUT INTO CODE###
	###Are allele freqs defined? No=next line.
	####For each dbSNP match, does Var-allele or Complement of Var-allele match any of the observed alleles with freqs? No=next line. N.B sometimes observed but no freqs present
	###If more than one allele freq. due to complement, choose rarest allele freq, if rare (MAF<0.01) =next line.
	###Common SNP match found = add to %dbSNP.
	###################
	
	my ($CurrentDir, $dbSNPmatchfile, $probableSNPs, $patient, $lane, $het_hom)=@_;
	
	my %dbSNP;
	my @NotIndbSNP;
	my $common_dbSNPcount=0;
	open DBSNP, $dbSNPmatchfile or die "cannot open $dbSNPmatchfile";
	snploop: while (<DBSNP>) {
		my $ln=$_;
		chomp $ln;
		my @line = split(/\s+/, $ln);
		my $v = $line[2];
		my $cv = $line[2];
		$cv=~tr/ACGT/TGCA/;
		if(!defined $line[5]){next snploop;}#don't inlcude SNPs with undefined MAFs in the dbSNPs to filter list
		my @MAFs = split(/,/, $line[5]);
		my @Alleles = split(/,/,$line[6]);
		my $numAlle=@Alleles;
		my %Alls=(); 
		for(my $a=0; $a<$numAlle; $a++){
			$Alls{$Alleles[$a]}=$MAFs[$a];
		}
		if(!exists $Alls{$v} and !exists $Alls{$cv}){next snploop;}#variant allele not observed
		my $lowest=1;
		if(exists $Alls{$v}){$lowest=$Alls{$v};}
		if(exists $Alls{$cv} and $Alls{$cv}<$lowest){$lowest=$Alls{$cv};}
		if($lowest<=0.01){next snploop;}#don't include rare dbSNP alleles 
		if(!exists $dbSNP{$line[0]."_".$line[1]."_".$v}) {
			$dbSNP{$line[0]."_".$line[1]."_".$v}="indbSNP";
			$common_dbSNPcount++;
		}
	}
	print "$patient $lane $het_hom:\nCommon dbSNPs $common_dbSNPcount\n";
	close(DBSNP);
	
	open ALLSNPS, $probableSNPs or die "cannot open $probableSNPs";
	my $head=<ALLSNPS>;
	var_loop: while (<ALLSNPS>) {
		chomp $_;
		my @Line = split(/\t/, $_);
		if($het_hom=~/heterozygous/ and $Line[8]>85){next var_loop;}
		if($het_hom=~/homozygous/ and $Line[8]<85){next var_loop;}
		my $newstart = $Line[3]; #0-based in dbSNP matches and SnpsOnTarget!!!!
		if (!exists $dbSNP{$Line[1]."_".$newstart."_".$Line[5]}) {
			push @NotIndbSNP, $Line[1]."\t".$newstart."\t".$Line[4]."\/".$Line[5]."\t".$Line[8]."\n";
		}
	}
	my $NonMatchOutput="VariantsNotMatchingdbSNP_MAF-0.01_".$patient."_".$lane."_".$het_hom.".txt";
	open NONMATCH, ">$NonMatchOutput";
	print NONMATCH @NotIndbSNP;
	close NONMATCH;

}


sub MySQLforVariants{
	#my ($CurrentDir, $probableSNPs, $patient, $lane, $het_hom)=@_;
	#my $IDplh=$patient."_".$lane."_".$het_hom;
	my ($input_file, $output_file)=@_;
	my @OutputFile;
	
	open ALLSNPS, $input_file or die "Cannot open $input_file"; 
	#my $header=<ALLSNPS>;
	
	A_loop: while(<ALLSNPS>){
		chomp $_;
		my @SplitLine=split(/\t/, $_);
		#if($het_hom=~/heterozygous/ and $SplitLine[8]>85){next A_loop;}
		#if($het_hom=~/homozygous/ and $SplitLine[8]<=85){next A_loop;}
		my $Chr=$SplitLine[0];
		my $Pos=$SplitLine[1];
		my $Mut=$SplitLine[2];
		#my $ComplimentMutant=$SplitLine[5];
		#$ComplimentMutant=~tr/ACGT/TGCA/;
		push @OutputFile, $Chr.",".$Pos.",".$Mut."\n";#only forward strand mut listed in file, no (rev strand) complement listed
	}
	close ALLSNPS;
	
	my $ReformattedOutput = $input_file."_Comparing2dbSNP135.temp.txt";
	open REFORMAT, ">$ReformattedOutput";
	print REFORMAT @OutputFile;
	close REFORMAT;
	
	print "intermidiate file is generated with reformated SNPs info from $input_file.\n";
	
	
	##MySQL comparison.
	print"Connecting the database....\n";
	my $ds="DBI:mysql:Data:headnode";
	my $user="nhrg";
	my $dbh=DBI->connect($ds,$user)||die "Cant connect to MySQL";
	print"Connection is made.\n";
	
	#####Generate a table containing all the exons which should be covered by our baits.
	print"Querying the databse......\n";
	my $AllMatchingSNPsTable = $dbh->prepare("CREATE TABLE AllMatchingSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(25), col4 CHAR(15), col5 CHAR(10), col6 CHAR(54), col7 CHAR(10))");
	$AllMatchingSNPsTable -> execute();
	my $LocatedSNPsTable = $dbh -> prepare("CREATE TABLE LocatedSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(25))");
	$LocatedSNPsTable -> execute();
	my $LoadLocatedSNPsTable = $dbh -> prepare("LOAD DATA LOCAL INFILE '$ReformattedOutput' INTO TABLE LocatedSNPs FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n'");
	$LoadLocatedSNPsTable -> execute();
	my $Index1 = $dbh -> prepare("CREATE INDEX col_index1 ON LocatedSNPs (col1)");
	my $Index2 = $dbh -> prepare("CREATE INDEX col_index2 ON LocatedSNPs (col2)");
	my $Index3 = $dbh -> prepare("CREATE INDEX col_index3 ON LocatedSNPs (col3)");
	$Index1 -> execute();
	$Index2 -> execute();
	$Index3 -> execute();
	
	print"Comparing........";
	my $LoadMatchingSNPsTable1 = $dbh -> prepare(
	"INSERT INTO AllMatchingSNPs
	SELECT LocatedSNPs.col1,LocatedSNPs.col2,LocatedSNPs.col3,snp135.name,snp135.observed,snp135.alleleFreqs,snp135.alleles
	FROM snp135, LocatedSNPs WHERE snp135.chrom=LocatedSNPs.col1 AND snp135.chromStart=LocatedSNPs.col2"); #output all chr pos matches 
	$LoadMatchingSNPsTable1 -> execute();
	print"Done!\n";
	
	print"Generating result.......";
	my $AllMatchingSNPsStatement = "SELECT col1,col2,col3,col4,col5,col6,col7 FROM AllMatchingSNPs";
	my $AllData = $dbh -> prepare($AllMatchingSNPsStatement);
	$AllData -> execute();
	my $AllDataContent = $AllData -> fetchall_arrayref();
	print"Done!\n";
	
	print"Cleaning tables and intermidiate files......";
	#####Removing all the tables form MySQL.
	my $RemoveAllTable = $dbh -> prepare("DROP TABLE AllMatchingSNPs");
	$RemoveAllTable -> execute();
	my $RemoveLocatedTable = $dbh -> prepare("DROP TABLE LocatedSNPs");
	$RemoveLocatedTable -> execute();
	`rm $ReformattedOutput`;
	print"Done!\n";
	
	
	my @Temp;
	my %RemoveDuplicates;
	my @RemoveDuplicatesArray;
	my $DuplicateCounter = 0;
	foreach my $Line(@$AllDataContent){
		push @Temp, "@$Line\n";
	}
	foreach my $TempLine(@Temp){
		if(!exists $RemoveDuplicates{$TempLine}){ ##There will be some duplicates as in matching to dbSNP, some of the reverse compliments are the same-such as C/G (Leads to G/C,C/G,GC).
		$RemoveDuplicates{$TempLine} = $DuplicateCounter;
		push @RemoveDuplicatesArray, $TempLine;
		$DuplicateCounter++;
		}
	}
	
	print "###########################################\n";
	print "########### output results ################\n";
	print "###########################################\n";
	print "$output_file\n";
	my $AllMatchesOutput = $output_file;
	open ALL, ">$AllMatchesOutput";
	print ALL @RemoveDuplicatesArray;
	#print ALL @Temp;
	close ALL;
	return($AllMatchesOutput);
}

        
