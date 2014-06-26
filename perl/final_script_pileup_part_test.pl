#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

#my ($PileupFile,$VarscanFile, $GenoFile, $HomoFile) = @ARGV;
my @GenoFiles;
my @PileupFiles;
#my @VarFiles;

@GenoFiles = ('F08_272_s234_genodata.txt_on_Targets',
'F05_339_s56_genodata.txt_on_Targets',
'F05_232_s78_genodata.txt_on_Targets');

@PileupFiles = ('../s234_bowtie_best_hg18.sam.bam.sorted.bam.pileup_on_GenoChipPositions',
'../s56/s56_bowtie_best_hg18.sam.bam.sorted.bam.pileup_on_GenoChipPositions',
'../s78/s78_bowtie_best_hg18.sam.bam.sorted.bam.pileup_on_GenoChipPositions');

my $quality_cutoff = 15; 
my $cvrg_cutoff = 20;

my $integratedGenoFile = InteGenoData(\@GenoFiles);

foreach my $file (@PileupFiles) {
	#my $filtered_pileup_file = FilterPileupWithPositions($file,$integratedGenoFile);
	my $processed_filtered_pileup_file = PileupParser($file, $quality_cutoff, $cvrg_cutoff);
	Merge2FilesWithGenomePosition($integratedGenoFile, $processed_filtered_pileup_file);
}
exit;

sub Merge2FilesWithGenomePosition {
	#the script will merge two tables in seperate file into the first file
	#first two columns as indentifier to map each row
	#N/A will be added to each column if the genome position is not found in the second file
	my ($file1, $file2) = @_;
	my %hash1;
	my %hash2;
	my $addedColumnNumber = 5;
	
	open FILE2, $file2 or die "Cannot open the file $file2";
	while (my $line = <FILE2>) {
		if ($line =~ m/^chr/) {
			chomp $line;
			my @splitLine = split(/\t/, $line);
			my $chrom = $splitLine[0]; chomp $chrom;
			my $pos =  $splitLine[1]; chomp $pos;
			if(!exists($hash2{$chrom."_".$pos})){
				for (my $i=2; $i<scalar @splitLine; $i++){
					if($i==2){
						$hash2{$chrom."_".$pos}= $splitLine[$i];
					} else {
						$hash2{$chrom."_".$pos} = $hash2{$chrom."_".$pos}."\t".$splitLine[$i];
					}
				}			
			} else {
				print "Warning: Found duplicate entry in $file2 on position $chrom $pos\tAction: Removed!!\n";
			}
			
		}
	}
	close FILE2;

	my @output;
	open FILE1, $file1 or die "Cannot open the file $file1";
	my $line_number = 0;	
	while (my $line = <FILE1>) {
		$line_number += 1;
		if ($line =~ m/^chr/) {
			chomp $line;
			my @splitLine = split(/\t/, $line);
			my $chrom = $splitLine[0]; chomp $chrom;
			my $pos =  $splitLine[1]; chomp $pos;
			if(!exists($hash1{$chrom."_".$pos})){
				$hash1{$chrom."_".$pos} = 0;
				if(exists($hash2{$chrom."_".$pos})) {
					my @temp_array = split(/\t/,$hash2{$chrom."_".$pos});
					if(scalar @temp_array == 5){
						print "scalar @temp_array == 5       at $line_number.\n";
						push @output, $line."\t".$hash2{$chrom."_".$pos}."\n";
					} else {
						my $x = scalar @temp_array;
						print $x."\n";
						foreach my $value (@temp_array) {
							print $value."\n";
						}
						print "scalar @temp_array <= 5       at $line_number.\n";
						my $temp_string;
						for (my $i=scalar @temp_array + 1; $i<=$addedColumnNumber; $i++){
							if ($i == scalar @temp_array + 1) {
								$temp_string = "N/A";
							} else {
								$temp_string = $temp_string."\t"."N/A"; 
							}
						}
						push @output, $line."\t".$hash2{$chrom."_".$pos}."\t".$temp_string."\n";
					}
				} else {
					#print "scalar temp_array == 0       at $line_number.\n";
					my $temp_string ="";
					for (my $i=1; $i<=$addedColumnNumber; $i++){
						$temp_string = $temp_string."\t"."N/A"; 
					}
					push @output, $line.$temp_string."\n";
				}
			} else {
				print "Warning: Found duplicate entry in $file1 on position $chrom $pos\tAction: Removed!!\n";
			}
		}
	}
	close FILE1;
	
	open OUTPUT, ">$file1" or die "Cannot open the file $file1 to write";
	print OUTPUT @output;
	close OUTPUT;
}

sub InteGenoData{
	my ($files) = @_;
	my @GenoFiles = @$files;
	my %GenoHash;
	foreach my $file (@GenoFiles) {
		print "reading data from file $file.....\n";
		open GENO, $file or die "Cannot open the file $file";
		while (my $line = <GENO>) {
			if ($line =~ m/^chr/) {
				chomp $line;
				my @splitLine = split(/\t/, $line);
				my $chrom = $splitLine[0]; chomp $chrom;
				my $pos =  $splitLine[1]; chomp $pos;
				my @alleles = split(/\//, $splitLine[2]);
				if(!exists($GenoHash{$chrom."_".$pos})){
					$GenoHash{$chrom."_".$pos} = $alleles[0]."\t".$alleles[1]."\t".$splitLine[3];			
				}
			}
		}
		close(GENO);
		print "done!\n";
	}

	#To add pheno data of each patient
	my $i = 0;
	foreach my $file (@GenoFiles) {
		open GENO, $file or die "Cannot open the file $file";
		$i += 1;
		print "reading phenotype data from file $file for patient $i.....\n";
		my %temp_hash;
		while (my $line = <GENO>) {
			if ($line =~ m/^chr/) {
				chomp $line;
				my @splitLine = split(/\t/, $line);
				my $chrom = $splitLine[0]; chomp $chrom;
				my $pos =  $splitLine[1]; chomp $pos;
				$temp_hash{$chrom."_".$pos} = $splitLine[4];
			}
		}
		close GENO;
	
		foreach my $key (keys %GenoHash) {
			if(exists $temp_hash{$key}){
				$GenoHash{$key} = $GenoHash{$key}."\t".$temp_hash{$key};
			}else{
				$GenoHash{$key} = $GenoHash{$key}."\t"."no call";
			}
		}
	}
	my $CurrentDir=`pwd`;
	chomp $CurrentDir;
	my $output_file = $CurrentDir."/"."Integrated_Geno_File_of".$i."_patients.txt";

	open OUTPUT, ">$output_file" or die "Cannot open the file $output_file";	
	foreach my $key (keys %GenoHash) {
		my @seperate_key = split("_", $key);
		print OUTPUT $seperate_key[0]."\t".$seperate_key[1]."\t";
		print OUTPUT $GenoHash{$key}."\n";
	}
	close OUTPUT;	
	return $output_file;
}

sub FilterPileupWithPositions{
	my($PileupFile,$GenoFile) = @_;
	my %positionsHash;
	open InputPositions, $GenoFile or die "Cannot open $GenoFile";
	while (my $line = <InputPositions>) {
		chomp $line;
		#print $line."\n";
		if ($line =~ m/^chr/) {
			my @splitLine = split(/\s+/, $line);
			my $chrom = $splitLine[0]; chomp $chrom;
			my $pos =  $splitLine[1]; chomp $pos;
			$positionsHash{$chrom."_".$pos} = 0;
		}
	}
	close(InputPositions);
	my @output;
	
	open PILEUP, $PileupFile or die "Cannot open $PileupFile";
	while (my $line=<PILEUP>) {
		chomp $line;
		if ($line =~ m/chr/) {
			my @splitLine = split(/\s+/, $line);
			my $chrom = $splitLine[0];
			my $pos = $splitLine[1];
			if(exists $positionsHash{$chrom."_".$pos}){
				push @output, $line."\n";
			} 
		}
	}
	close(PILEUP);
	#ouput to file
    my $OutputFile = $PileupFile."_on_GenoChipPositions";
    open OUTPUT, ">$OutputFile" or die "Cannot create $OutputFile";
    print OUTPUT @output;
    close OUTPUT;
    return ($OutputFile);
}

sub  Homo_or_Hetero_varscan{

#filter varscan result with higher coverage and var coverge before using this scipt for better accruracy 
#this will only compare the position of Variant and ignoring the ref/var change.

	my ($var_file) = @_;
	print " reading from $var_file\n";

	my %var_hetero;
	my %var_homo;

	my $output_homo = $var_file."_homo";
	my $output_hetero = $var_file."_hetero";

	open VAR, $var_file or die "cannot open file $var_file";
	while (my $line =<VAR>) {
        if ($line =~ m/^chr/){
        	chomp $line;
            my @words = split(/\s+/,$line);
            $words[6] =~ s/\%//g;
            if ($words[6] >= 80 ) {
            	if (exists $var_hetero{$words[0]."_".$words[1]}) {
            		next;
            	} else {
            		if (exists $var_homo{$words[0]."_".$words[1]})  {
            			$var_hetero{$words[0]."_".$words[1]} = $line;
            			delete $var_homo{$words[0]."_".$words[1]};
            		} else {
            			$var_homo{$words[0]."_".$words[1]} = $line;
            		}
            	}
            } elsif ($words[6] >= 30 && $words[6]< 80){
            	$var_hetero{$words[0]."_".$words[1]} = $line;
            }
        }
	}
	close(VAR);

	open HOMO, ">$output_homo" or die "can not open the file $output_homo";
	for my $key ( keys %var_homo ) {
        print HOMO $var_homo{$key}."\n";
	}
	close HOMO;

	open HETERO, ">$output_hetero" or die "can not open the file $output_hetero";
	for my $key ( keys %var_hetero ) {
        print HETERO $var_hetero{$key}."\n";
	}
	close HETERO;
	return $output_homo;
}

sub PileupParser {
	#modified from the script on the webpage 
	#https://bitbucket.org/galaxy/galaxy-central/src/tip/tools/samtools/pileup_parser.pl
	#die "Usage: pileup_parser.pl <in_file> <qv cutoff> <coverage cutoff> <out_file>\n" unless @ARGV == 4;

	#my $in_file = $ARGV[0];
	#my $quality_cutoff = $ARGV[1]; # phred scale integer
	#my $cvrg_cutoff = $ARGV[2]; # unsigned integer
	#my $out_file = $ARGV[3];
	
	my ($in_file,$quality_cutoff, $cvrg_cutoff) = @_;
	my $out_file = $in_file."_parsed";
	my $ref_base_column = 2; # 0 based
	my $read_bases_column = 4; # 0 based
	my $base_quality_column = 5; # 0 based
	my $cvrg_column = 3; # 0 based
	my $SNPs_only = "NO"; # set to "Yes" to print only positions with SNPs; set to "No" to pring everything
	my $bed = "NO"; #set to "Yes" to convert coordinates to bed format (0-based start, 1-based end); set to "No" to leave as is
	my $coord_column = 1; #0 based 
	my $total_diff = "Yes"; # set to "Yes" to print total number of deviant based
	my $print_qual_bases = "Yes"; #set to "Yes" to print quality and read base columns
	my $invalid_line_counter = 0;
	my $first_skipped_line = "";
	my %SNPs = ('A',0,'T',0,'C',0,'G',0);
	my $above_qv_bases = 0;
	my $SNPs_exist = 0;
	my $out_string = "";
	my $diff_count = 0;
	
	open (IN, "<$in_file") or die "Cannot open $in_file $!\n";
	open (OUT, ">$out_file") or die "Cannot open $out_file $!\n";
	
	while (<IN>) {
		chop;
		next if m/^\#/;
		my @fields = split /\t/;
		next if $fields[ $ref_base_column ] eq "*"; # skip indel lines
		my $read_bases   = $fields[ $read_bases_column ];
		die "Coverage column" . ($cvrg_column+1) . " contains non-numeric values. Check your input parameters as well as format of input dataset." if ( not isdigit $fields[ $cvrg_column ] );
		next if $fields[ $cvrg_column ] < $cvrg_cutoff;
		my $base_quality = $fields[ $base_quality_column ];
		if ($read_bases =~ m/[\$\^\+-]/) {
			$read_bases =~ s/\^.//g; #removing the start of the read segement mark
			$read_bases =~ s/\$//g; #removing end of the read segment mark
			while ($read_bases =~ m/[\+-]{1}(\d+)/g) {
				my $indel_len = $1;
				$read_bases =~ s/[\+-]{1}$indel_len.{$indel_len}//; # remove indel info from read base field
			}
		}
		if ( length($read_bases) != length($base_quality) ) {
			$first_skipped_line = $. if $first_skipped_line eq "";
	        ++$invalid_line_counter;
	        next;
		}
		# after removing read block and indel data the length of read_base 
		# field should identical to the length of base_quality field
		
		my @bases = split //, $read_bases;
		my @qv    = split //, $base_quality;
		for my $base ( 0 .. @bases - 1 ) {
			if ( ord( $qv[ $base ] ) - 33 >= $quality_cutoff and $bases[ $base ] ne '*'){
				++$above_qv_bases;
				if ( $bases[ $base ] =~ m/[ATGC]/i ){
					$SNPs_exist = 1; 
					$SNPs{ uc( $bases[ $base ] ) } += 1;
					$diff_count += 1;
				} elsif ( $bases[ $base ] =~ m/[\.,]/ ) {
					$SNPs{ uc( $fields[ $ref_base_column ] ) } += 1;
				}
			}
		}
		
		my $ref_count = 0;
		
		my $out_string = $fields[0]."\t".$fields[$coord_column]."\t".$above_qv_bases."/".$fields[$cvrg_column]."\t".$fields[$ref_base_column]."/".$SNPs{uc($fields[ $ref_base_column])};
		foreach my $key (keys %SNPs){
			if ($key ne uc($fields[$ref_base_column]) and $SNPs{$key} > 0) {
				$out_string = $out_string."\t".$key."/".$SNPs{$key};
			}
		}
		print OUT $out_string."\n";
		
		$above_qv_bases = 0;
		%SNPs = ('A',0,'T',0,'C',0,'G',0);
	}
	

	print "Skipped $invalid_line_counter invalid line(s) beginning with line $first_skipped_line\n" if $invalid_line_counter > 0;
	close IN;
	close OUT;
	return ($out_file);
}		



