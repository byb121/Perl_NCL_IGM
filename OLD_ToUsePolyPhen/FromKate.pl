
#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use DBD::mysql;
use Math::CDF qw(pbinom); #Generate probabilities and quantiles from several statistical probability functions 
## why needs to calculate this probability
use LWP::Simple; #simple procedural interface to LWP -Library for WWW in Perl

my $CurrentDir=`pwd`;
chomp $CurrentDir;

my $helens = "helens_NOMATCH.txt";
my $dbSNPmatchfile = ParseHelenNoMatchFile($helens);
#my $dbSNPmatchfile = MySQLforVariants_yaobo($CurrentDir, $probableSNPs);
my $NovelVariantsWithCCDS = ExtractNovelVariantRefseqs_Yaobo($CurrentDir, $dbSNPmatchfile);

my $CCDS = "/home/kjs/CCDSref/CCDS.20090327.txt";
#print $probableSNPs."\n";
my ($VariantcDNAfile) = ConvertVariantsTocDNA($CCDS, $NovelVariantsWithCCDS);
my $PolyPhenInput = IdentifyNonSynonymousChanges($VariantcDNAfile);

#RunPolyPhen($CurrentDir);
#ParsePolyPhenOutput($CurrentDir);


####################       SUBROUTINES           ##########################
sub ParseHelenNoMatchFile {# This is to take helen's no mathch file and parse it to the format that can be used with Katies sub ExtractNovelVariantRefseqs
	 my ($helens) = @_;
	 my @Output;
	 open HELEN, $helens or die "Cannot open $helens";
	 while(<HELEN>){
	 	my @SplitLine = split(/\s+/, $_);
	 	if (@SplitLine != 8) {
	 		print $_."\n";
	 	} else {
	 		my $patient = $SplitLine[0]; chomp $patient;
	 		my $chr = $SplitLine[1]; chomp $chr;
	 		my $gene = $SplitLine[2]; chomp $gene;
	 		my $position_0 = $SplitLine[3]; chomp $position_0;
	 		my $ref = $SplitLine[4]; chomp $ref;
	 		my $var = $SplitLine[5]; chomp $var;
	 		my $total_coverage = $SplitLine[6]; chomp $total_coverage;
	 		my $var_coverage = $SplitLine[7]; chomp $var_coverage;
	 		if (length $var ==2){
	 			$var =~ s/$ref//;
	 			if (length $var == 2) {
	 				my @chars = split(//, $var);
	 				push @Output, $patient."\t".$chr."\t".$gene."\t".$position_0."\t".$ref."/".$chars[0]."\t".$total_coverage."\t".$var_coverage."\n";
	 				push @Output, $patient."\t".$chr."\t".$gene."\t".$position_0."\t".$ref."/".$chars[1]."\t".$total_coverage."\t".$var_coverage."\n";
	 			} else {
	 				push @Output, $patient."\t".$chr."\t".$gene."\t".$position_0."\t".$ref."/".$var."\t".$total_coverage."\t".$var_coverage."\n";
	 			}
	 		} else {
	 			push @Output, $patient."\t".$chr."\t".$gene."\t".$position_0."\t".$ref."/".$var."\t".$total_coverage."\t".$var_coverage."\n";
	 		}
	 	} 
	 }
	 my $OutputFile = $helens."_parsedForKatie";
	 open OUTPUT, ">$OutputFile" or die "Cannot create $OutputFile";
	 print OUTPUT @Output;
	 close OUTPUT;
	 return ($OutputFile);
}


sub FilterVarScanResults{ ### currently formatted for varscan version 1 (on Abarth.ncl.ac.uk). For version 2 will need to change $TotalCoverage, $NumberRefStrands and $NumberVariantStrands column names
        my ($CurrentDir, $input) = @_;
        my @FilteredVarScanResults;

        open INPUT, $input or die "Cannot open $input";
        while(<INPUT>){
        	my @SplitLine=split(' ', $_);
            my $Chrom=$SplitLine[1]; chomp $Chrom;
            my $Position=$SplitLine[3]; chomp $Position;
            my $Ref=$SplitLine[4]; chomp $Ref;
            my $Variant=$SplitLine[5]; chomp $Variant;
            my $TotalCoverage=$SplitLine[6]; chomp $TotalCoverage;
            #my $RefCoverage=$SplitLine[5]; chomp $RefCoverage;
            my $VarCoverage=$SplitLine[7]; chomp $VarCoverage;
            my $RefCoverage = $TotalCoverage-$VarCoverage;
                        #my $NumberRefStrands=$SplitLine[10]; chomp $NumberRefStrands;
                        #my $NumberVariantStrands=$SplitLine[11]; chomp $NumberVariantStrands; 
                        #my $x=$VarCoverage; 
                        #my $n=$TotalCoverage; 
                        #my $p=0.95;             
                        #if(($NumberVariantStrands==2)&&($TotalCoverage>=20)&&($VarCoverage>=5)){
                                #CDF Binomial function below generates cumulative probabilities of having n number of
                                #hits or fewer. Over a certain threshold the CDF Binomial Function just labels 
                                #the binomial probability as zero.
                      #          my $sub_binomCDF=pbinom($x,$n,$p);
                       #         if($sub_binomCDF>=0.1){
                        #                #print $sub_binomCDF."\n";
                         #               push @FilteredVarScanResults, $Chrom." ".$Position." ".$Ref."/".$Variant." ".$RefCoverage."/".$VarCoverage."\n";
                          #      }
                           #     else{
             print $Chrom." ".$Position." ".$Ref."/".$Variant." ".$RefCoverage."/".$VarCoverage."\n";

                              
			}
                
        
        close INPUT;
        
        my $OutputFile="AllSNPs_ParsedOnProbability.txt";
        open OUTPUT, ">$OutputFile" or die "Cannot create $OutputFile";
        print OUTPUT @FilteredVarScanResults;
        close OUTPUT;   
        return ($OutputFile);
}

sub MySQLforVariants{
        my ($CurrentDir, $probableSNPs)=@_;
        
        my %UniqueVariants;
        my $UniqueVariantCounter=0;
        my @OutputFile;
        open ALLSNPS, $probableSNPs or die "Cannot open $probableSNPs"; 
        while(<ALLSNPS>){
                my @SplitLine=split(' ', $_);
                my $Chr=$SplitLine[0];
                my $Pos=$SplitLine[1]-1;
                chomp $Pos;
                if(!exists $UniqueVariants{$Pos}){
                        $UniqueVariants{$Pos}=$UniqueVariantCounter;
                        $UniqueVariantCounter++;
                }
                my @SplitCoverages=split('/', $SplitLine[3]);
                my $RefCoverage=$SplitCoverages[0];
                my $MutCoverage=$SplitCoverages[1];
                my @SplitRefAndMut=split('/', $SplitLine[2]);           
                my $Ref=$SplitRefAndMut[0];
                my $Mut=$SplitRefAndMut[1];
                my $ComplimentMutant=$SplitRefAndMut[1];
                $ComplimentMutant=~tr/ACGT/TGCA/;
                my $ComplimentRef=$SplitRefAndMut[0];
                $ComplimentRef=~tr/ACGT/TGCA/; 
                push @OutputFile, $Chr.",".$Pos.",".$Ref."/".$Mut."\n";
                push @OutputFile, $Chr.",".$Pos.",".$Mut."/".$Ref."\n";
                push @OutputFile, $Chr.",".$Pos.",".$ComplimentRef."/".$ComplimentMutant."\n";
                push @OutputFile, $Chr.",".$Pos.",".$ComplimentMutant."/".$ComplimentRef."\n";
        }       
        close ALLSNPS;
#       print $UniqueVariantCounter;
        my $ReformattedOutput="ReformattedVariantFile.txt";
        open REFORMAT, ">$ReformattedOutput";
        print REFORMAT @OutputFile;
        close REFORMAT;
        
        ##MySQL comparison.     
        my $ds="DBI:mysql:nimblegen_analysis:localhost";
        my $user="kjs";
        my $passwd="kjs010203";
        my $dbh=DBI->connect($ds,$user,$passwd)||die "Cant connect to MySQL";
        
        #####Generate a table containing all the exons which should be covered by our baits.
        my $AllMatchingSNPsTable=$dbh->prepare("CREATE TABLE AllMatchingSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(5), col4 CHAR(15))");
    $AllMatchingSNPsTable->execute();
    my $LocatedSNPsTable=$dbh->prepare("CREATE TABLE LocatedSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(5))");
        $LocatedSNPsTable->execute();
    my $LoadLocatedSNPsTable=$dbh->prepare("LOAD DATA INFILE '$CurrentDir/ReformattedVariantFile.txt' INTO TABLE LocatedSNPs FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n'");
        $LoadLocatedSNPsTable->execute();
        my $Index1=$dbh->prepare("CREATE INDEX col_index1 ON LocatedSNPs (col1)");
        my $Index2=$dbh->prepare("CREATE INDEX col_index2 ON LocatedSNPs (col2)");
        my $Index3=$dbh->prepare("CREATE INDEX col_index3 ON LocatedSNPs (col3)");
        $Index1->execute();
        $Index2->execute();
        $Index3->execute();     
        my $LoadMatchingSNPsTable1=$dbh->prepare("INSERT INTO AllMatchingSNPs 
        SELECT LocatedSNPs.col1,LocatedSNPs.col2,LocatedSNPs.col3,dbsnp_first.rsid   
        FROM dbsnp_first, LocatedSNPs WHERE dbsnp_first.chr=LocatedSNPs.col1 AND dbsnp_first.chr_start=LocatedSNPs.col2 AND dbsnp_first.obs_change=LocatedSNPs.col3");  
        $LoadMatchingSNPsTable1->execute();
        my $LoadMatchingSNPsTable2=$dbh->prepare("INSERT INTO AllMatchingSNPs 
        SELECT LocatedSNPs.col1,LocatedSNPs.col2,LocatedSNPs.col3,dbsnp_second.rsid   
        FROM dbsnp_second, LocatedSNPs WHERE dbsnp_second.chr=LocatedSNPs.col1 AND dbsnp_second.chr_start=LocatedSNPs.col2 AND dbsnp_second.obs_change=LocatedSNPs.col3");  
        $LoadMatchingSNPsTable2->execute();
        my $AllMatchingSNPsStatement="SELECT col1,col2,col3,col4 FROM AllMatchingSNPs";
        my $AllData=$dbh->prepare($AllMatchingSNPsStatement);
        $AllData->execute();
        my $AllDataContent=$AllData->fetchall_arrayref();
                
        #####Removing all the tables form MySQL.
        my $RemoveAllTable=$dbh->prepare("DROP TABLE AllMatchingSNPs");
        $RemoveAllTable->execute();
        my $RemoveLocatedTable=$dbh->prepare("DROP TABLE LocatedSNPs");
        $RemoveLocatedTable->execute();

        my @Temp;
        my %RemoveDuplicates;
        my @RemoveDuplicatesArray;
        my $DuplicateCounter=0;
        foreach my $Line(@$AllDataContent){
                        push @Temp, "@$Line\n";
        }
        my %UniqueRsNumbers;
        my $UniquersNumberCounter=0;
        foreach my $TempLine(@Temp){
                my @SplitLine=split(' ', $TempLine);
                my $rsNumber=$SplitLine[3];
                if(!exists $UniqueRsNumbers{$rsNumber}){ ##This hash is to obtain how many unique rsNumbers there are in the data (i.e. Unique variants matching dbSNP).
                        $UniqueRsNumbers{$rsNumber}=$UniquersNumberCounter;
                        $UniquersNumberCounter++;
                }
                if(!exists $RemoveDuplicates{$TempLine}){ ##There will be some duplicates as in matching to dbSNP, some of the reverse compliments are the same-such as C/G (Leads to G/C,C/G,GC).
                        $RemoveDuplicates{$TempLine}=$DuplicateCounter;
                        push @RemoveDuplicatesArray, $TempLine;
                        $DuplicateCounter++;
                }
        }
        my $AllMatchesOutput="AlldbSNPmatches.txt";
        open ALL, ">$AllMatchesOutput";
        print ALL @RemoveDuplicatesArray;
        close ALL;
        return($AllMatchesOutput);
}

sub MySQLforVariants_Yaobo{
        my ($CurrentDir, $probableSNPs)=@_;
        
        my $db_input = $CurrentDir."/".$probableSNPs;
        
        ##MySQL comparison.     
        my $ds="DBI:mysql:nimblegen_analysis:localhost";
        my $user="kjs";
        my $passwd="kjs010203";
        my $dbh=DBI->connect($ds,$user,$passwd)||die "Cant connect to MySQL";
        
        #####Generate a table containing all the exons which should be covered by our baits.
        my $AllMatchingSNPsTable=$dbh->prepare("CREATE TABLE AllMatchingSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(5), col4 CHAR(15))");
    $AllMatchingSNPsTable->execute();
    my $LocatedSNPsTable=$dbh->prepare("CREATE TABLE LocatedSNPs (col1 CHAR(5), col2 INT(15), col3 CHAR(5))");
        $LocatedSNPsTable->execute();
    my $LoadLocatedSNPsTable=$dbh->prepare("LOAD DATA INFILE '$db_input' INTO TABLE LocatedSNPs FIELDS TERMINATED BY ',' LINES TERMINATED BY '\n'");
        $LoadLocatedSNPsTable->execute();
        my $Index1=$dbh->prepare("CREATE INDEX col_index1 ON LocatedSNPs (col1)");
        my $Index2=$dbh->prepare("CREATE INDEX col_index2 ON LocatedSNPs (col2)");
        my $Index3=$dbh->prepare("CREATE INDEX col_index3 ON LocatedSNPs (col3)");
        $Index1->execute();
        $Index2->execute();
        $Index3->execute();     
        my $LoadMatchingSNPsTable1=$dbh->prepare("INSERT INTO AllMatchingSNPs 
        SELECT LocatedSNPs.col1,LocatedSNPs.col2,LocatedSNPs.col3,dbsnp_first.rsid   
        FROM dbsnp_first, LocatedSNPs WHERE dbsnp_first.chr=LocatedSNPs.col1 AND dbsnp_first.chr_start=LocatedSNPs.col2 AND dbsnp_first.obs_change=LocatedSNPs.col3");  
        $LoadMatchingSNPsTable1->execute();
        my $LoadMatchingSNPsTable2=$dbh->prepare("INSERT INTO AllMatchingSNPs 
        SELECT LocatedSNPs.col1,LocatedSNPs.col2,LocatedSNPs.col3,dbsnp_second.rsid   
        FROM dbsnp_second, LocatedSNPs WHERE dbsnp_second.chr=LocatedSNPs.col1 AND dbsnp_second.chr_start=LocatedSNPs.col2 AND dbsnp_second.obs_change=LocatedSNPs.col3");  
        $LoadMatchingSNPsTable2->execute();
        my $AllMatchingSNPsStatement="SELECT col1,col2,col3,col4 FROM AllMatchingSNPs";
        my $AllData=$dbh->prepare($AllMatchingSNPsStatement);
        $AllData->execute();
        my $AllDataContent=$AllData->fetchall_arrayref();
                
        #####Removing all the tables form MySQL.
        my $RemoveAllTable=$dbh->prepare("DROP TABLE AllMatchingSNPs");
        $RemoveAllTable->execute();
        my $RemoveLocatedTable=$dbh->prepare("DROP TABLE LocatedSNPs");
        $RemoveLocatedTable->execute();

        my @Temp;
        my %RemoveDuplicates;
        my @RemoveDuplicatesArray;
        my $DuplicateCounter=0;
        foreach my $Line(@$AllDataContent){
                        push @Temp, "@$Line\n";
        }
        my %UniqueRsNumbers;
        my $UniquersNumberCounter=0;
        foreach my $TempLine(@Temp){
                my @SplitLine=split(' ', $TempLine);
                my $rsNumber=$SplitLine[3];
                if(!exists $UniqueRsNumbers{$rsNumber}){ ##This hash is to obtain how many unique rsNumbers there are in the data (i.e. Unique variants matching dbSNP).
                        $UniqueRsNumbers{$rsNumber}=$UniquersNumberCounter;
                        $UniquersNumberCounter++;
                }
                if(!exists $RemoveDuplicates{$TempLine}){ ##There will be some duplicates as in matching to dbSNP, some of the reverse compliments are the same-such as C/G (Leads to G/C,C/G,GC).
                        $RemoveDuplicates{$TempLine}=$DuplicateCounter;
                        push @RemoveDuplicatesArray, $TempLine;
                        $DuplicateCounter++;
                }
        }
        my $AllMatchesOutput="AlldbSNPmatches.txt";
        open ALL, ">$AllMatchesOutput";
        print ALL @RemoveDuplicatesArray;
        close ALL;
        return($AllMatchesOutput);
}


sub ExtractNovelVariantRefseqs{
        my ($CurrentDir, $dbSNPmatchfile, $probableSNPs)=@_;
        my %dbSNP;
        my @NotIndbSNP;
        open DBSNP, $dbSNPmatchfile or die "cannot open $dbSNPmatchfile";
        while (<DBSNP>) {
                chomp $_;
                my @line = split(/\s+/, $_);
                if (!exists $dbSNP{$line[0]."_".$line[1]}) {
                        $dbSNP{$line[0]."_".$line[1]}="indbSNP";
                }
        }

        
        open ALLSNPS, $probableSNPs or die "cannot open $probableSNPs";
        while (<ALLSNPS>) {
                chomp $_;
                my @Line = split(/\s+/, $_);
                my $newstart = $Line[1]-1;
                if (!exists $dbSNP{$Line[1]."_".$newstart}) {
                        push @NotIndbSNP, $Line[0]."\t".$newstart."\t".$Line[2]."\n"; #chr."\t".pos."\t".Ref/Variant
                }
        }
        my $NonMatchOutput="VariantsNotMatchingdbSNP.txt";
        open NONMATCH, ">$NonMatchOutput";
        print NONMATCH @NotIndbSNP;
        close NONMATCH;

        ##MySQL is used to get the CCDS id's as well as the CCDS strands for 
        ##all the variants which did not match to dbSNP variants.
        my $ds="DBI:mysql:nimblegen_analysis:localhost";
        my $user="kjs";
        my $passwd="kjs010203";
        my $dbh=DBI->connect($ds,$user,$passwd)||die "Cant connect to MySQL";
        
        my $NonMatchingSNPsPlusCCDSTable=$dbh->prepare("CREATE TABLE NonMatchingSNPsPlusCCDS (col1 CHAR(5), col2 INT(15), col3 CHAR(7), col4 CHAR(25), col5 CHAR(2))");
    $NonMatchingSNPsPlusCCDSTable->execute();
    my $NonMatchingSNPsNoCCDSTable=$dbh->prepare("CREATE TABLE NonMatchingSNPsNoCCDS (col1 CHAR(5), col2 INT(15), col3 CHAR(7))");
        $NonMatchingSNPsNoCCDSTable->execute();
    my $LoadNonMatchingSNPsNoCCDSTable=$dbh->prepare("LOAD DATA INFILE '$CurrentDir/$NonMatchOutput' INTO TABLE NonMatchingSNPsNoCCDS FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n'");
        $LoadNonMatchingSNPsNoCCDSTable->execute();
        my $Index1=$dbh->prepare("CREATE INDEX col_index1 ON NonMatchingSNPsNoCCDS (col1)");
        my $Index2=$dbh->prepare("CREATE INDEX col_index2 ON NonMatchingSNPsNoCCDS (col2)");
        $Index1->execute();
        $Index2->execute();
        my $LoadNonMatchingSNPsPlusCCDSTable=$dbh->prepare("INSERT INTO NonMatchingSNPsPlusCCDS
        SELECT NonMatchingSNPsNoCCDS.col1,NonMatchingSNPsNoCCDS.col2,NonMatchingSNPsNoCCDS.col3,CCDS_20090327.CCDSID,CCDS_20090327.cds_strand   
        FROM CCDS_20090327, NonMatchingSNPsNoCCDS WHERE CCDS_20090327.chrom=NonMatchingSNPsNoCCDS.col1 AND CCDS_20090327.cds_from<=NonMatchingSNPsNoCCDS.col2 AND CCDS_20090327.cds_to>=NonMatchingSNPsNoCCDS.col2");  
        $LoadNonMatchingSNPsPlusCCDSTable->execute();
        my $NonMatchingSNPsPlusCCDSStatement="SELECT col1,col2,col3,col4,col5 FROM NonMatchingSNPsPlusCCDS";
        my $AllData=$dbh->prepare($NonMatchingSNPsPlusCCDSStatement);
        $AllData->execute();
        my $AllDataContent=$AllData->fetchall_arrayref();
        
        #####Removing all the tables form MySQL.
        my $RemoveAllTable=$dbh->prepare("DROP TABLE NonMatchingSNPsNoCCDS");
        $RemoveAllTable->execute();
        my $RemoveLocatedTable=$dbh->prepare("DROP TABLE NonMatchingSNPsPlusCCDS");
        $RemoveLocatedTable->execute();

        my @Temp;
        my %RemoveDuplicates;
        my @RemoveDuplicatesArray;
        my $DuplicateCounter=0;
        foreach my $Line(@$AllDataContent){
                push @Temp, "@$Line\n";
        }       
        foreach my $TempLine(@Temp){
                if(!exists $RemoveDuplicates{$TempLine}){
                        $RemoveDuplicates{$TempLine}=$DuplicateCounter;
                        push @RemoveDuplicatesArray, $TempLine;
                        $DuplicateCounter++;
                }
        }
        my $AllMatchesOutput="AllNonMatchingVariantsPlusCCDS.txt";
        open ALL, ">$AllMatchesOutput";
        print ALL @RemoveDuplicatesArray;
        close ALL;
        
        return($AllMatchesOutput);
}

sub ExtractNovelVariantRefseqs_Yaobo{
        my ($CurrentDir, $dbSNPmatchfile)=@_;
        #my %dbSNP;
        my @NotIndbSNP;
        
        open SNPS, $dbSNPmatchfile or die "cannot open $dbSNPmatchfile";
        while (<SNPS>) {
                chomp $_;
                my @Line = split(/\s+/, $_);
                my $newstart = $Line[3];
                push @NotIndbSNP, $Line[1]."\t".$newstart."\t".$Line[4]."\t".$Line[2]."\n"; #chr."\t".pos."\t".Ref/Variant
        }
        my $NonMatchOutput="VariantsNotMatchingdbSNP.txt";
        open NONMATCH, ">$NonMatchOutput";
        print NONMATCH @NotIndbSNP;
        close NONMATCH;

        ##MySQL is used to get the CCDS id's as well as the CCDS strands for 
        ##all the variants which did not match to dbSNP variants.
        my $ds="DBI:mysql:nimblegen_analysis:localhost";
        my $user="kjs";
        my $passwd="kjs010203";
        my $dbh=DBI->connect($ds,$user,$passwd)||die "Cant connect to MySQL";
        
        my $NonMatchingSNPsPlusCCDSTable=$dbh->prepare("CREATE TABLE NonMatchingSNPsPlusCCDS (col1 CHAR(5), col2 INT(15), col3 CHAR(7), col4 CHAR(25), col5 CHAR(2), col6 CHAR(25))");
    $NonMatchingSNPsPlusCCDSTable->execute();
    my $NonMatchingSNPsNoCCDSTable=$dbh->prepare("CREATE TABLE NonMatchingSNPsNoCCDS (col1 CHAR(5), col2 INT(15), col3 CHAR(7), col4 CHAR(25))");
        $NonMatchingSNPsNoCCDSTable->execute();
    my $LoadNonMatchingSNPsNoCCDSTable=$dbh->prepare("LOAD DATA INFILE '$CurrentDir/$NonMatchOutput' INTO TABLE NonMatchingSNPsNoCCDS FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n'");
        $LoadNonMatchingSNPsNoCCDSTable->execute();
        my $Index1=$dbh->prepare("CREATE INDEX col_index1 ON NonMatchingSNPsNoCCDS (col1)");
        my $Index2=$dbh->prepare("CREATE INDEX col_index2 ON NonMatchingSNPsNoCCDS (col2)");
        $Index1->execute();
        $Index2->execute();
        my $LoadNonMatchingSNPsPlusCCDSTable=$dbh->prepare("INSERT INTO NonMatchingSNPsPlusCCDS
        SELECT NonMatchingSNPsNoCCDS.col1,NonMatchingSNPsNoCCDS.col2,NonMatchingSNPsNoCCDS.col3,CCDS_20090327.CCDSID,CCDS_20090327.cds_strand,NonMatchingSNPsNoCCDS.col4   
        FROM CCDS_20090327, NonMatchingSNPsNoCCDS WHERE CCDS_20090327.chrom=NonMatchingSNPsNoCCDS.col1 AND CCDS_20090327.cds_from<=NonMatchingSNPsNoCCDS.col2 AND CCDS_20090327.cds_to>=NonMatchingSNPsNoCCDS.col2");  
        $LoadNonMatchingSNPsPlusCCDSTable->execute();
        my $NonMatchingSNPsPlusCCDSStatement="SELECT col1,col2,col3,col4,col5, col6 FROM NonMatchingSNPsPlusCCDS";
        my $AllData=$dbh->prepare($NonMatchingSNPsPlusCCDSStatement);
        $AllData->execute();
        my $AllDataContent=$AllData->fetchall_arrayref();
        
        #####Removing all the tables form MySQL.
        my $RemoveAllTable=$dbh->prepare("DROP TABLE NonMatchingSNPsNoCCDS");
        $RemoveAllTable->execute();
        my $RemoveLocatedTable=$dbh->prepare("DROP TABLE NonMatchingSNPsPlusCCDS");
        $RemoveLocatedTable->execute();

        my @Temp;
        my %RemoveDuplicates;
        my @RemoveDuplicatesArray;
        my $DuplicateCounter=0;
        foreach my $Line(@$AllDataContent){
                push @Temp, "@$Line\n";
        }       
        foreach my $TempLine(@Temp){
                if(!exists $RemoveDuplicates{$TempLine}){
                        $RemoveDuplicates{$TempLine}=$DuplicateCounter;
                        push @RemoveDuplicatesArray, $TempLine;
                        $DuplicateCounter++;
                }
        }
        my $AllMatchesOutput="AllNonMatchingVariantsPlusCCDS.txt";
        open ALL, ">$AllMatchesOutput";
        print ALL @RemoveDuplicatesArray;
        close ALL;
        
        return($AllMatchesOutput);
}

sub ConvertVariantsTocDNA{
        ##Create Hashes of Start and End positions for CCDS exons##
        my ($CCDS, $variants) =@_;
        my %StartHash;
        my %StopHash;
        my %newStartHash;
        my %StartHashneg;
        my %StopHashneg;
        my $count;
        open CCDS, $CCDS or die "cannot open $CCDS";
        while (<CCDS>) {
                my @line = split (/\t/, $_);
                if ($line[5] eq "Public" or $line[5] eq "Reviewed, update pending") {
                        $count++;
                        my @starts;
                        my @stops;
                        my $chr = "chr".$line[0];
                        my $CCDSID = $line[4];
                        my $exons = substr($line[9], 1, length($line[9])-2);
                        my @exons = split (/,/, $exons);
                        foreach my $exon (@exons) {
                                my @ends = split(/\-/, $exon);
                                my $start = $ends[0];
                                $start =~ s/^\s+//;
                                my $end = $ends[1];
                                push @starts, $start;
                                push @stops, $end; ## need to check if these numbers are 1 or 0 based and add or subtract 1 accordingly. Variant numbers are 0-based.
                        }
                        if ($line[6] eq '+') {
                                if (!exists $StartHash{$CCDSID."_".$chr}) {
                                        $StartHash{$CCDSID."_".$chr}=[@starts];
                                        $StopHash{$CCDSID."_".$chr}=[@stops];
                                } else {
                                        print $CCDSID."\n";
                                }
                        } elsif ($line[6] eq '-') {
                                my @startsrev = reverse @starts;
                                my @stopsrev = reverse @stops;
                                $StartHashneg {$CCDSID."_".$chr} = [@stopsrev]; 
                                $StopHashneg {$CCDSID."_".$chr} = [@startsrev];
                        }
                                
                                                
                }
        }
        foreach my $CCDS (keys %StartHashneg) {
                for (my $i=0; $i < (scalar @{$StartHashneg{$CCDS}}); $i++) {
                        if ($i == 0) {
                                if (exists $StartHash{$CCDS}) {
                        }
                                $StartHash{$CCDS}=[($StartHashneg{$CCDS}[0] - $StartHashneg{$CCDS}[$i])+1];
                                $StopHash{$CCDS}=[($StartHashneg{$CCDS}[0] - $StopHashneg{$CCDS}[$i])+1];
                        } else {
                        push @{$StartHash{$CCDS}}, ($StartHashneg{$CCDS}[0] - $StartHashneg{$CCDS}[$i])+1;
                        push @{$StopHash{$CCDS}}, ($StartHashneg{$CCDS}[0] - $StopHashneg{$CCDS}[$i])+1;
                        }
                }
        }
### still losing 27 CCDS IDs somewhere along the line ###       

#print $count."\n";
#print scalar (keys(%StartHash))."\n";
#print scalar (keys(%StopHash))."\n";
#print scalar (keys(%StartHashneg))."\n";
#print scalar (keys(%StopHashneg));
        
        foreach my $currentCCDS (keys %StartHash) {
                $newStartHash{$currentCCDS}=[1];
                for (my $i=1; $i < (scalar @{$StopHash{$currentCCDS}}); $i++) { 
                        push @{$newStartHash{$currentCCDS}}, ($StopHash{$currentCCDS}[$i-1]-$StartHash{$currentCCDS}[$i-1]+$newStartHash{$currentCCDS}[$i-1])+1;
                }
        }
        
        my $CurrentDir=`pwd`;
        chomp $CurrentDir;
        my @Output;
        chomp $variants;
        open VARIANTS, $variants or die "cannot open $variants";
        while (my $Var=<VARIANTS>) {
                my @VariantsSplit = split(/\s+/, $Var);
                my $variantID = $VariantsSplit[3]."_".$VariantsSplit[0];
                my $strand = $VariantsSplit[4];
                my $variantpos = "";
                if ($strand eq '+') {
                        $variantpos = $VariantsSplit[1]; 
                } elsif ($strand eq '-') {
                        if (exists $StartHashneg{$variantID}) {
                                $variantpos = $StartHashneg{$variantID}[0] - ($VariantsSplit[1])+1; 
                        }
                }
                if (exists $newStartHash{$variantID}) {
                        my $Found = 0;
                        for (my $i=0; $i < (scalar @{$newStartHash{$variantID}}) && !$Found; $i++) {
                                if ($variantpos >= $StartHash{$variantID}[$i] && $variantpos <= $StopHash{$variantID}[$i]) {
                                        $Found=1;
                                        push @Output, $variantID, "\t", ($variantpos-$StartHash{$variantID}[$i])+$newStartHash{$variantID}[$i], "\t", $VariantsSplit[0], "\t", $VariantsSplit[1]+1, "\t", $VariantsSplit[2]."\t".$strand."\t".$VariantsSplit[5]."\n"; ## currently full output ##
                                }
                        }
                }
        }

        my $Output = "$CurrentDir/VariantcDNApos_CCDS.txt";
        open OUTPUT, ">$Output" or die "Cannot open $Output";
        print OUTPUT @Output;
        close OUTPUT;
        @Output=();
        return($Output);
}       

sub IdentifyNonSynonymousChanges{
        ##Script is from Katie-19 March 2010.
        my ($VariantscDNA) = @_;
        chomp $VariantscDNA;
        my $CurrentDir=`pwd`;
        chomp $CurrentDir;
        my %SequenceHash;
        my @SNPout;
        my @ProteinFastas;

        open VARIANTS, $VariantscDNA or die "cannot open $VariantscDNA";
        while (<VARIANTS>) {
                chomp $_;
                my @variants = split (/\t/ , $_);
                my @idsplit = split (/_/ , $variants[0]);
                my @change = split (/\// , $variants[4]);
                my $strand = $variants[5];
                my $sequence = ExtractFastaSequence ($idsplit[0], $idsplit[1]);

                my $refprotein;
                my $variantprotein;
                
                if ($strand eq '+') {
                        $refprotein = translate_seq($sequence, 1, length($sequence));
                        my $variantsequence = create_variant_sequence ($sequence, $variants[1]-1, $change[0], $change[1]);
                        unless ($variantsequence eq "variant locus not found"){
                                $variantprotein = translate_seq ($variantsequence, 1, length($sequence)); 
                        }
                }

                elsif ($strand eq '-') {
                        $change[0] =~tr/ATCG/TAGC/;
                        $change[1] =~tr/ATCG/TAGC/;
                        $refprotein = translate_seq($sequence, 1, length($sequence));
                        my $variantsequence = create_variant_sequence ($sequence, $variants[1]-1, $change[0], $change[1]);
                        unless ($variantsequence eq "variant locus not found"){
                                $variantprotein = translate_seq ($variantsequence, 1, length($sequence)); 
                        }
                }
                                

                if ($refprotein eq $variantprotein) {
                                print "$idsplit[0]\t$idsplit[1]\t$variants[1]\tsynonymous change\n";
                } else {
                        my $aachange = find_aminoacid_change ($refprotein, $variantprotein);
                        push @SNPout, "$idsplit[0]\t$idsplit[1]\t$variants[1]\t$aachange\t$variants[2]\t$variants[3]\t$variants[4]\t$variants[5]\t$variants[6]\n";    ##\t$patient[0]_$patient[1]_$variants[3]\n";
                }

                if (!exists $SequenceHash{">".$idsplit[0]}) {
                        $SequenceHash{">".$idsplit[0]} = $refprotein;
                }
                
                $refprotein = ();
                $variantprotein = ();
        }
                        
        
                foreach my $keys (keys %SequenceHash) {
                        push @ProteinFastas, "$keys\n$SequenceHash{$keys}\n";
                }
        my $PolyPhenFile = CreatePolyphenInput (\@SNPout);
#       my $SNPsOut = $CurrentDir."/NonsynSNPs.txt";
#       open OUTPUT, ">$SNPsOut" or die "Cannot open $SNPsOut";
#       print OUTPUT @SNPout;
#       close OUTPUT;
        my $ProteinFastaOut = $CurrentDir."/Protein.fas";
        open OUTPUT2, ">$ProteinFastaOut" or die "Cannot open $ProteinFastaOut";
        print OUTPUT2 @ProteinFastas;
        close OUTPUT2;  
        return(0);
}

sub ExtractFastaSequence {
        my ($CCDSID, $chr) = @_;
        my @Fasta;
        my $FastaFile = "/home/kjs/CCDSref/CCDS_20090327.fa";
        open FASTA, $FastaFile or die "cannot open $FastaFile";
        my $extractline=0;
        while (<FASTA>) {
                if ($_ =~ /^>$CCDSID\|Hs36.3\|$chr/) {
                        push @Fasta, $_;
                        $extractline=1;
                } elsif ($_=~/^>/ && $extractline==1) {
                        last;
                } elsif ($extractline == 1) {
                        push @Fasta, $_;
                }
        }
        my $sequence = extract_sequence_from_fasta_data (@Fasta);
        return ($sequence);
}

sub extract_sequence_from_fasta_data { 
        ##Script is from Katie-19 March 2010.
    my(@fasta_file) = @_;
    my $sequence = "";
    foreach my $line (@fasta_file) {

        # discard blank line
        if ($line =~ /^\s*$/) {
            next;

        # discard comment line
        } elsif($line =~ /^\s*#/) {
            next;

        # discard fasta header line
        } elsif($line =~ /^>/) {
            next;

        # keep line, add to sequence string
        } else {
            $sequence .= $line;
        }
    }

    # remove non-sequence data (in this case, whitespace) from $sequence string
    $sequence =~ s/\s//g; 
    return ($sequence);
}

sub translate_seq{ 
        ##Script is from Katie-19 March 2010.
        my ($seq, $start, $end) = @_;
        my $protein;
        unless($end) {
       $end = length($seq);
    }
        return (dnatranslate ( substr ( $seq, $start-1, $end -$start+1)));
}

sub dnatranslate{
        ##Script is from Katie-19 March 2010.
    my($dna) = @_;

    # Initialize variables
    my $protein = '';
        
        # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }

    return ($protein);
}

sub codon2aa{
        ##Script is from Katie-19 March 2010.
        my($codon) = @_;
    $codon = uc $codon;
    my(%genetic_code)=(
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );
    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}

sub create_variant_sequence{
        ##Script is from Katie-19 March 2010.
        my ($sequence, $varpos, $changefrom, $changeto) = @_;
        my @sequence = split (// , $sequence);
        if ($sequence[$varpos] eq $changefrom) {
                $sequence[$varpos] =~ s/$sequence[$varpos]/$changeto/;
                my $variantseq = join ('',@sequence);
                $variantseq =~ s/\s//g;
                return $variantseq;
        }       
        else { 
                print $sequence[$varpos]."\t".$changefrom."\t"."variant locus not found\n";
                return "variant locus not found";
        }       
}

sub find_aminoacid_change{
        ##Script is from Katie-19 March 2010.
        my ($ref, $var) = @_;
        my $match = ($ref ^ $var) =~ /^(\x00*)/;

        substr $ref, 0, $+[1], '' if $match;
        substr $var, 0, $+[1], '' if $match;

        $match = (reverse ($ref) ^ reverse ($var))=~ /^(\x00*)/;
        substr $ref, -$+[1], length ($ref), '' if $match;
        substr $var, -$+[1], length ($var), '' if $match;
        return "$ref\t$var";
}

sub CreatePolyphenInput{
        my ($NonSynSNPs_ref) = @_;
        my $PolyPhenInput = "PolyPhenInput.txt";
        my @PolyPhenArray;
	my @OutputWithDnaPosition;
	my $ExtraOutput = "PolyPhenInput.ExtraInfor.txt";

	
        chomp $CurrentDir;
        my @NonSynSNPs = @$NonSynSNPs_ref;
        foreach my $SNP (@NonSynSNPs) {
                chomp $SNP;
                my @snp = split (/\t/, $SNP);
                my $codonpos = $snp[2]/3;
                my $positioncorrected = $codonpos;
                if ($codonpos=~/(\d)\.(\d)/) {
                        $positioncorrected = round ($codonpos);
                }
		#print "snp0 is ".$snp[0]."\n";
 		#print "snp3 is ".$snp[3]."\n";
		#print "snp4 is ".$snp[4]."\n";
                push @PolyPhenArray,"?\t?\t".$snp[0]."\t".$positioncorrected."\t".$snp[3]."\t".$snp[4]."\n";
		push @OutputWithDnaPosition, $snp[0]."\t".$positioncorrected."\t".$snp[3]."\t".$snp[4]."\t".$snp[5]."\t".$snp[6]."\t".$snp[7]."\t".$snp[8]."\t".$snp[9]."\n";
        }
        open OUTPUT, ">$PolyPhenInput" or die "cannot open $PolyPhenInput";
        print OUTPUT @PolyPhenArray;
        close OUTPUT;
	open OUTPUT2, ">$ExtraOutput" or die "cannot open $ExtraOutput";
        print OUTPUT2 @OutputWithDnaPosition;
        close OUTPUT2;
        return(0);
}

sub round{
        my ($number) = shift;
        return (int ($number + 1));
}

sub RunPolyPhen{
        my $CurrentDir=$_[0];
        chomp $CurrentDir;
        my $PolyPhenDir="/home/data/PolyPhen.1.18";
        chomp $PolyPhenDir;
        `cp Protein.fas $PolyPhenDir`;
        `cp PolyPhenInput.txt $PolyPhenDir`;
        `perl /home/data/PolyPhen.1.18/pph.pl -s Protein.fas PolyPhenInput.txt 1>PolyPhenOutput2.txt 2>PolyPhenLog`;
        #`cp $PolyPhenDir/PolyPhenOutput2.txt $CurrentDir`;
        #`cp $PolyPhenDir/PolyPhenLog $CurrentDir`;
        #`rm $PolyPhenDir/PolyPhenOutput2.txt`;
        #`rm $PolyPhenDir/PolyPhenLog`;
        #`rm $PolyPhenDir/Protein.fas`;
        #`rm $PolyPhenDir/PolyPhenInput.txt`;
}



exit;