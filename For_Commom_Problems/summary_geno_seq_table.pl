#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Math::CDF qw(pbinom);


sub SummayTableAndCalculateSensiAndSpeci {
	##########################################################
	###### The script will also calculate the sensitivity and specificity of 
	###### variant detection of sequencing data based on the geno chip data
	###### eg: geno chip data is considered as reality
	###### 			Scheme			
	######        Ref	Var	TP	TN	FP	FN
	###### Chip	+	+	Yes			
	###### Seq	-	+				
	###### Chip	+	+				Yes
	###### Seq	+	-				
	###### Chip	+	+	Yes			
	###### Seq	+	+				
	###### Chip	-	+	Yes			
	###### Seq	-	+				
	###### Chip	-	+				Yes
	###### Seq	+	-				
	###### Chip	-	+	Yes			
	###### Seq	+	+				
	###### Chip	+	-			Yes	
	###### Seq	-	+				
	###### Chip	+	-		Yes		
	###### Seq	+	-				
	###### Chip	+	-			Yes	
	###### Seq	+	+										
	###### Chip	+	+				Yes
	###### Seq	-	-				
	###### Chip	+	-		Yes		
	###### Seq	-	-				
	###### Chip	-	+				Yes
	###### Seq	-	-				
	###### "no call" is considered invalid value
	###### When both allele A and allele B are related to ref base
	###### the one equals ref represents the ref base
	###### the other on represents the var base
	##########################################################
	
	
	#my $input_file = "Integrated_Geno_File_of3_patients.txt";
	#my $input_file = "testfile.txt";
	
	my ($input_file) = @_;
	my @output;
	my $output_file = $input_file."_summarized";
	
	my @true_positives = (0, 0 ,0);
	my @true_negatives = (0, 0, 0);
	my @false_positives = (0, 0, 0);
	my @false_negatives = (0, 0, 0);
	
	open OUTPUT, ">$output_file" or die "Cannot open the file $output_file.";
	
	open INPUT, $input_file or die "Cannot open the file $input_file.";
	while (my $line = <INPUT>) {
		chomp $line;
		my @fields = split(/\t/,$line);
		my $chr = $fields[0];
		my $coor = $fields[1];
		my $g_a = $fields[2];
		my $g_b = $fields[3];
		my @geno_homo_hetero = ($fields[5],$fields[6],$fields[7]);
		my @seq_homo_hetero;
		my @most_covered_bases; #array of hashes
		my @total_coverage;
		my $ref;
		my $var;
		my $a_ref;
		my $b_ref;
		my $a_var;
		my $b_var;
		my $output_string;
		
		my @ref_positions = (9,14,19);
		foreach my $ref_position (@ref_positions){
			if ($fields[$ref_position] ne "N/A" && !$ref){
				my @ref_temp = split(/\//, $fields[$ref_position]);
				$ref = $ref_temp[0];
			}
		}
		
		if (!$ref){
			$output_string = $chr."\t".$coor."\t".$g_a."\t".$g_b."\t"."No sufficient coverage"."\t"."No sufficient coverage"."\t"."No sufficient coverage"."\t"."No sufficient coverage";
			next;
		}
		## lines(SNPs)with no sufficient sequencing data coverage on the position will be removed from the SNP file.
		
		$ref = uc($ref);
	
		my $c_ref =  ComplementaryBase($ref); 
		
		$output_string =  $chr."\t".$coor."\t".$g_a."\t".$g_b."\t".$ref;
		
		#to assign the total coverage values
		if ($fields[8] ne "N/A") {
			my @coverage_numbers = split (/\//, $fields[8]);
			push @total_coverage, $coverage_numbers[0];
		} else {
			push @total_coverage,0;
		}
		
		if ($fields[13] ne "N/A") {
			my @coverage_numbers = split (/\//, $fields[13]);
			push @total_coverage, $coverage_numbers[0];
		} else {
			push @total_coverage, 0;
		}
		
		if ($fields[18] ne "N/A") {
			my @coverage_numbers = split (/\//, $fields[18]);
			push @total_coverage, $coverage_numbers[0];
		} else {
			push @total_coverage,0;
		}
		
		#to read in two most covered based on that posions for patients
		my $temp_hash = {};
		for (my $i = 9; $i<= 12; $i++) {
			if ($fields[$i] ne "N/A") {
				my @tmp = split (/\//, $fields[$i]);
				$temp_hash->{uc($tmp[0])} = $tmp[1];
			}
		}
			
		if (scalar(keys %$temp_hash)  == 0) {
			print "$chr $coor: patient 1 No sufficient coverage to define any base. \n";
			push @most_covered_bases, $temp_hash;
		} elsif (scalar(keys %$temp_hash) >0 && scalar(keys %$temp_hash) <=2) {
			push @most_covered_bases, $temp_hash;
		} else {
			my @sorted = sort {$temp_hash->{$b} <=> $temp_hash->{$a}} keys %$temp_hash; 
			my $temp_hash_2 = {};
			$temp_hash_2->{$sorted[0]} = $temp_hash->{$sorted[0]}; 
			$temp_hash_2->{$sorted[1]} = $temp_hash->{$sorted[1]};
			push @most_covered_bases, $temp_hash_2;
		}
		
		$temp_hash = {};
		for (my $i = 14; $i<= 17; $i++) {
			if ($fields[$i] ne "N/A") {
				my @tmp = split (/\//, $fields[$i]);
				$temp_hash->{uc($tmp[0])} = $tmp[1];
			}
		}	
		if (scalar(keys %$temp_hash)  ==0) {
			print "$chr $coor: patient 2 No sufficient coverage to define any base. \n";
			push @most_covered_bases, $temp_hash;
		} elsif (scalar(keys %$temp_hash) >0 && scalar(keys %$temp_hash) <=2) {
			push @most_covered_bases, $temp_hash;
		} else {
			my @sorted = sort {$temp_hash->{$b} <=> $temp_hash->{$a}} keys %$temp_hash; 
			my $temp_hash_2 = {};
			$temp_hash_2->{$sorted[0]} = $temp_hash->{$sorted[0]}; 
			$temp_hash_2->{$sorted[1]} = $temp_hash->{$sorted[1]};
			push @most_covered_bases, $temp_hash_2;
		}
		
		$temp_hash = {};
		for (my $i = 19; $i<= 22; $i++) {
			if ($fields[$i] ne "N/A") {
				my @tmp = split (/\//, $fields[$i]);
				$temp_hash->{uc($tmp[0])} = $tmp[1];
			}
		}	
		if (scalar(keys %$temp_hash)  ==0) {
			print "$chr $coor: patient 3 No sufficient coverage to define any base. \n";
			push @most_covered_bases, $temp_hash;
		} elsif (scalar(keys %$temp_hash) >0 && scalar(keys %$temp_hash) <=2) {
			push @most_covered_bases, $temp_hash;
		} else {
			my @sorted = sort {$temp_hash->{$b} <=> $temp_hash->{$a}} keys %$temp_hash; 
			my $temp_hash_2 = {};
			$temp_hash_2->{$sorted[0]} = $temp_hash->{$sorted[0]}; 
			$temp_hash_2->{$sorted[1]} = $temp_hash->{$sorted[1]};
			push @most_covered_bases, $temp_hash_2;
		} 
		
		
		#determine homo/heter zygous type from the sequencing data
		for(my $i =0 ; $i <3; $i++) {
			if ($total_coverage[$i] == 0 ) {
				push @seq_homo_hetero, "N/A";
			} elsif ( keys %{$most_covered_bases[$i]} == 0) {
				push @seq_homo_hetero, "N/A";
			} elsif ($total_coverage[$i] > 0 && keys %{$most_covered_bases[$i]} > 0) {
				my %new_temp_hash;
				foreach my $key (keys %{$most_covered_bases[$i]}){
					#print "pbinom x is $most_covered_bases[$i]{$key}\n";
					#print "pbinom total is $total_coverage[$i]\n";
					if (pbinom($most_covered_bases[$i]{$key}, $total_coverage[$i], 0.05) > 0.95) {
						$new_temp_hash{$key} = $most_covered_bases[$i]{$key};
					}
				}
				
				$most_covered_bases[$i] = \%new_temp_hash;
	
				my %new_temp_hash_2 = %new_temp_hash;
				if(keys %new_temp_hash_2 == 0) {
					push @seq_homo_hetero, "N/A";
				} elsif (keys %new_temp_hash_2 == 1) {
					if(exists $new_temp_hash_2{$ref} || exists $new_temp_hash_2{$c_ref}) {
						push @seq_homo_hetero, "RR";
					}else {
						push @seq_homo_hetero, "VV";
					}
				} elsif (keys %new_temp_hash_2 == 2) {
					if(exists $new_temp_hash_2{$ref} || exists $new_temp_hash_2{$c_ref}) {
						if (exists $new_temp_hash_2{$ref}) {
							delete $new_temp_hash_2{$ref};
						}else {
							delete $new_temp_hash_2{$c_ref};
						}
						
						if(exists $new_temp_hash_2{$ref} || exists $new_temp_hash_2{$c_ref}) {
							push @seq_homo_hetero, "RR";
						} else {
							push @seq_homo_hetero, "RV";
						}
					} else {
						my @keys = keys %new_temp_hash_2;
						if ($keys[0] ne ComplementaryBase($keys[1])) {
							push @seq_homo_hetero, "RV";
							print "$chr $coor: possible error for reference base at postion.\n";
						}else {
							push @seq_homo_hetero, "VV";
						}
					}
				} else {
					print "$chr $coor :";
					print "new_temp_hash has invalid elements number:";
					print keys %new_temp_hash_2;
					print "\n";
				}
			}
		}
	
		#compare A B bases to ref
		
		if (CompareAwithB($g_a,$ref) eq "No Match"){
			$a_ref = "N/A";
		} else {
			if (CompareAwithB($g_a,$ref) eq "A=B") {
				$a_ref = "A=Ref";
			}else {
				$a_ref = "A=C_Ref";
			}	
		}	
	 	
	 	if (CompareAwithB($g_b,$ref) eq "No Match"){
			$b_ref = "N/A";
		} else {
			if (CompareAwithB($g_b,$ref) eq "A=B") {
				$b_ref = "B=Ref";
			}else {
				$b_ref = "B=C_Ref";
			}	
		}
		
		if (CompareAwithB($g_b,$ref) eq "No Match" && CompareAwithB($g_a,$ref) eq "No Match"){
			print "$chr $coor: What a Ref!!!!!\n";
		}
		
		######################################
		##### determine if the geno chip data is consistant with sequencing data
		###########################################################################
		##### and determine number of true positives, true negatives
		##### ,false positives and false negatives (take geno chip data as reality).
		##### This is used calculate the sensitivity and specificity of Sequencing data to detect variants
		######################################
		
	
		
		for(my $i =0 ; $i <3; $i++) {
			
			if ($geno_homo_hetero[$i] eq "no call"  || keys %{$most_covered_bases[$i]} == 0 ){ 
			# if nocall in genotyping data for the patient or there's no base is confirmed on the position from the seqencing data
				if ($geno_homo_hetero[$i] eq "no call" ) {
					$output_string = $output_string."\t"."no call";
				} else {
					$output_string = $output_string."\t"."No sufficient coverage";
					############################## 17/02/2011 ########################
					##### To handle no coverage but have calls on geno chip
					##### To determine true/false negatives ################
					if ($a_ref ne "N/A" && $b_ref eq "N/A" ) {
						if ($geno_homo_hetero[$i] eq "AA"){
							$true_negatives[$i] += 1;
						} else {#AB, BB
							$false_negatives[$i] += 1;
						}
					} elsif ( $a_ref eq "N/A" && $b_ref ne "N/A" ){
						if  ($geno_homo_hetero[$i] eq "BB") {
							$true_negatives[$i] += 1;
						} else { #AA, AB
							$false_negatives[$i] += 1;
						}
					} else { 
						if ($a_ref  eq "A=Ref") {
							if ($geno_homo_hetero[$i] eq "AA") {
								$true_negatives[$i] += 1;
							} else {
								$false_negatives[$i] += 1;
							}
						} elsif ($b_ref  eq "B=Ref") {
							if ($geno_homo_hetero[$i] eq "BB") {
								$true_negatives[$i] += 1;
							} else {
								$false_negatives[$i] += 1;
							}
						}else{
							# if both Allele A and Allele B are NON related to ref or complemtary ref
							$false_negatives[$i] += 1;
						}
					}
					################################# 17/02/2011 #############################################
				}
				next;
			} elsif (keys %{$most_covered_bases[$i]} == 1) {
			# if one base is confirmed on the position from the seqencing data
				if ($geno_homo_hetero[$i] eq "AA" || $geno_homo_hetero[$i] eq "BB"  ) {
				#geno typing data situations
					if ($a_ref ne "N/A" && $b_ref eq "N/A" ) { 
					# if Allele A is equal to ref or complemtary ref
						if ($geno_homo_hetero[$i] eq "AA") { # AA 
							if ($a_ref  eq "A=Ref"){
								if (exists $most_covered_bases[$i]{$g_a}) {
									$output_string = $output_string."\t"."Consistent";
									$true_negatives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									 $false_positives[$i] += 1;
								}
							} else {
								if (exists $most_covered_bases[$i]{ComplementaryBase($g_a)}) {
									$output_string = $output_string."\t"."Consistent";
									$true_negatives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									$false_positives[$i] += 1;
								}
							}
						} else { #BB
							if ($a_ref  eq "A=Ref") {
								if (exists $most_covered_bases[$i]{$g_b}) {
									$output_string = $output_string."\t"."Consistent";
									$true_positives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									if (exists $most_covered_bases[$i]{$g_a}) {
										$false_negatives[$i] += 1;
									} else {
										$false_positives[$i] += 1;
									}
								}
							} else {
								if (exists $most_covered_bases[$i]{ ComplementaryBase($g_b)}) {
									$output_string = $output_string."\t"."Consistent";
									$true_positives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									if (exists $most_covered_bases[$i]{ComplementaryBase($g_a)}) {
										$false_negatives[$i] += 1;
									} else {
										$false_positives[$i] += 1;
									}
								}
							}
						}
					} elsif ($a_ref eq "N/A" && $b_ref ne "N/A" ) {
					# if Allele B is equal to ref or complemtary ref
						if  ($geno_homo_hetero[$i] eq "AA") { #AA
							if ($b_ref  eq "B=Ref"){
								if (exists $most_covered_bases[$i]{$g_a}) {
									$output_string = $output_string."\t"."Consistent";
									$true_positives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									if (exists $most_covered_bases[$i]{$g_b}) {
										$false_negatives[$i] += 1;
									} else {
										$false_positives[$i] += 1;
									}
								}
							} else {
								if (exists $most_covered_bases[$i]{ ComplementaryBase($g_a)}) {
									$output_string = $output_string."\t"."Consistent";
									$true_positives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									if (exists $most_covered_bases[$i]{ComplementaryBase($g_b)}) {
										$false_negatives[$i] += 1;
									} else {
										$false_positives[$i] += 1;
									}
								}
							}
						} else { #BB
							if ($b_ref  eq "B=Ref"){
								if (exists $most_covered_bases[$i]{$g_b}) {
									$output_string = $output_string."\t"."Consistent";
									$true_negatives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									$false_positives[$i] += 1;
								}
							} else {
								if (exists $most_covered_bases[$i]{ ComplementaryBase($g_b)}) {
									$output_string = $output_string."\t"."Consistent";
									$true_negatives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									$false_positives[$i] += 1;
								}
							}
						}					
					} else { #Either Allele A or Allele B is either equal or complementary to ref
					#treat Allele A and Allele B are non related, one stands for Ref and the other one stands for Var
						if($a_ref eq "A=Ref"){
							if  ($geno_homo_hetero[$i] eq "AA"){
								if (exists $most_covered_bases[$i]{$g_a}) {
									$output_string = $output_string."\t"."Consistent";
									$true_negatives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									$false_positives[$i] += 1;
								}
							} else {#BB
								if (exists $most_covered_bases[$i]{$g_b}) {
									$output_string = $output_string."\t"."Consistent";
									$true_positives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									if(exists $most_covered_bases[$i]{$g_a}){
										$false_negatives[$i] += 1;
									} else {
										$false_positives[$i] += 1;
									}
								}
							}
						}elsif($b_ref eq "B=Ref"){
							if  ($geno_homo_hetero[$i] eq "AA"){
								if (exists $most_covered_bases[$i]{$g_a}) {
									$output_string = $output_string."\t"."Consistent";
									$true_positives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									if(exists $most_covered_bases[$i]{$g_b}){
										$false_negatives[$i] += 1;
									} else {
										$false_positives[$i] += 1;
									}
								}
							} else {#BB
								if (exists $most_covered_bases[$i]{$g_b}) {
									$output_string = $output_string."\t"."Consistent";
									$true_negatives[$i] += 1;
								} else {
									$output_string = $output_string."\t"."Inconsistent";
									$false_positives[$i] += 1;
								}
							}
						} else {
							##########################################################################
							# if both Allele A and Allele B  are NON related to ref or complemtary ref
							#rare but possible
							##but this does not necessrayly mean that var from sequencing data is not confirmed by chip data
							###require further code here. It's not right to treat the situation as N/A
							###########################################################################
							###########################################################################
							#changed code - 16/02/2011
							if(exists $most_covered_bases[$i]{$ref} ) {
								$output_string = $output_string."\t"."Inconsistent";
								$false_negatives[$i] += 1;
							} else {
								if ($most_covered_bases[$i]{$g_a} && $geno_homo_hetero[$i] eq "AA"){
									$output_string = $output_string."\t"."Consistent";
									$true_positives[$i] += 1;
								}elsif($most_covered_bases[$i]{$g_b} && $geno_homo_hetero[$i] eq "BB"){
									$output_string = $output_string."\t"."Consistent";
									$true_positives[$i] += 1;
								}else{
									$output_string = $output_string."\t"."Inconsistent";
									$false_positives[$i] += 1;
								}
							}
						}
					}
				} else { # "AB" appears on geno chip
					$output_string = $output_string."\t"."Inconsistent";
					if ($a_ref ne "N/A" && $b_ref eq "N/A" ){
						if ($a_ref  eq "A=Ref") {
							if (exists $most_covered_bases[$i]{$g_b}) {
								$true_positives[$i] += 1;
							} else{
								if (exists $most_covered_bases[$i]{$ref}) {
									$false_negatives[$i] += 1;
								} else {
									$false_positives[$i] += 1;
								}
							}
						} else {
							if (exists $most_covered_bases[$i]{ComplementaryBase($g_b)}) {
								$true_positives[$i] += 1;
							} else{
								if (exists $most_covered_bases[$i]{$ref}) {
									$false_negatives[$i] += 1;
								} else {
									$false_positives[$i] += 1;
								}
							}
						}
					} elsif ($b_ref ne "N/A" && $a_ref eq "N/A") {
						if ($b_ref  eq "B=Ref") {
							if (exists $most_covered_bases[$i]{$g_a}) {
								$true_positives[$i] += 1;
							} else{
								if (exists $most_covered_bases[$i]{$ref}) {
									$false_negatives[$i] += 1;
								} else {
									$false_positives[$i] += 1;
								}
							}
						} else {
							if (exists $most_covered_bases[$i]{ComplementaryBase($g_a)}) {
								$true_positives[$i] += 1;
							} else{
								if (exists $most_covered_bases[$i]{$ref}) {
									$false_negatives[$i] += 1;
								} else {
									$false_positives[$i] += 1;
								}
							}
						}
					} else {
						if($a_ref ne "N/A"){
							if (exists $most_covered_bases[$i]{$ref} || exists $most_covered_bases[$i]{ComplementaryBase($ref)} ){
								$true_negatives[$i] += 1;
							} else {
								$false_positives[$i] += 1;
							}
						}else {
							if (exists $most_covered_bases[$i]{$ref} || exists $most_covered_bases[$i]{ComplementaryBase($ref)} ){
								$false_negatives[$i] += 1;
							}elsif (exists $most_covered_bases[$i]{$g_a} || exists $most_covered_bases[$i]{$g_b} ){
								$true_positives[$i] += 1;
							} else {
								$false_positives[$i] += 1;
							}
						}
					}
				}
			} elsif (keys %{$most_covered_bases[$i]} == 2) {
				if ($geno_homo_hetero[$i] eq "AB") {
					if($a_ref ne "N/A" && $b_ref eq "N/A" ){
						if($a_ref  eq "A=Ref"){
							if(exists $most_covered_bases[$i]{$g_a} && exists $most_covered_bases[$i]{$g_b}) {
								$output_string = $output_string."\t"."Consistent";
								$true_positives[$i] += 1;
							} else {
								$output_string = $output_string."\t"."Inconsistent";
								if ( exists $most_covered_bases[$i]{$g_a} ||  exists $most_covered_bases[$i]{$g_b}) {
									if(exists $most_covered_bases[$i]{$g_a}){
										$false_positives[$i] += 1;
									}else{
										$true_positives[$i] += 1;
									}
								}else{
									$false_positives[$i] += 1;
								}
							}
						} else {
							if(exists $most_covered_bases[$i]{ComplementaryBase($g_a)} && exists $most_covered_bases[$i]{ComplementaryBase($g_b)}) {
								$output_string = $output_string."\t"."Consistent";
								$true_positives[$i] += 1;
							} else {
								$output_string = $output_string."\t"."Inconsistent";
								if ( exists $most_covered_bases[$i]{ComplementaryBase($g_a)} ||  exists $most_covered_bases[$i]{ComplementaryBase($g_b)}) {
									if(exists $most_covered_bases[$i]{ComplementaryBase($g_a)}){
										$false_positives[$i] += 1;
									}else{
										$true_positives[$i] += 1;
									}
								} else {
									$false_positives[$i] += 1;
								}
							}
						}
					} elsif ($a_ref eq "N/A" && $b_ref ne "N/A" ) {
						if($b_ref  eq "B=Ref"){
							if(exists $most_covered_bases[$i]{$g_a} && exists $most_covered_bases[$i]{$g_b}) {
								$output_string = $output_string."\t"."Consistent";
								$true_positives[$i] += 1;
							} else {
								$output_string = $output_string."\t"."Inconsistent";
								if ( exists $most_covered_bases[$i]{$g_a} ||  exists $most_covered_bases[$i]{$g_b}) {
									if(exists $most_covered_bases[$i]{$g_b}){
										$false_positives[$i] += 1;
									} else{
										$true_positives[$i] += 1;
									}
								} else {
									$false_positives[$i] += 1;
								}
							}
						} else {
							if(exists $most_covered_bases[$i]{ComplementaryBase($g_a)} && exists $most_covered_bases[$i]{ComplementaryBase($g_b)}) {
								$output_string = $output_string."\t"."Consistent";
								$true_positives[$i] += 1;
							} else {
								$output_string = $output_string."\t"."Inconsistent";
								if ( exists $most_covered_bases[$i]{ComplementaryBase($g_a)} ||  exists $most_covered_bases[$i]{ComplementaryBase($g_b)}) {
									if(exists $most_covered_bases[$i]{ComplementaryBase($g_b)}){
										$false_positives[$i] += 1;
									} else{
										$true_positives[$i] += 1;
									}
								} else {
									$false_positives[$i] += 1;
								}
							}
						}
					}else {
						if($a_ref eq "A=Ref"){
							if(exists $most_covered_bases[$i]{$g_a} && exists $most_covered_bases[$i]{$g_b}) {
								$output_string = $output_string."\t"."Consistent";
								$true_positives[$i] += 1;
							} else {
								$output_string = $output_string."\t"."Inconsistent";
								if ( exists $most_covered_bases[$i]{$g_a} ||  exists $most_covered_bases[$i]{$g_b}) {
									if(exists $most_covered_bases[$i]{$g_a}){
										$false_positives[$i] += 1;
									} else{
										$true_positives[$i] += 1;
									}
								} else {
									$false_positives[$i] += 1;
								}
							}
						}elsif($b_ref eq "B=Ref"){
							if(exists $most_covered_bases[$i]{$g_a} && exists $most_covered_bases[$i]{$g_b}) {
								$output_string = $output_string."\t"."Consistent";
								$true_positives[$i] += 1;
							} else {
								$output_string = $output_string."\t"."Inconsistent";
								if ( exists $most_covered_bases[$i]{$g_a} ||  exists $most_covered_bases[$i]{$g_b}) {
									if(exists $most_covered_bases[$i]{$g_b}){
										$false_positives[$i] += 1;
									} else{
										$true_positives[$i] += 1;
									}
								} else {
									$false_positives[$i] += 1;
								}
							}
						} else {
							# if both Allele A and Allele B  are NON related to ref or complemtary ref
							#rare but possible
							$output_string = $output_string."\t"."N/A2222";
							############################### 17/02/2011 ################################# Done here ###############
	
							if ( exists $most_covered_bases[$i]{$g_a} &&  exists $most_covered_bases[$i]{$g_b}) {
								$true_positives[$i] += 1;
							} else {
								$false_positives[$i] += 1;
							}
							############################### 17/02/2011 ################################# Done here ###############
						}
					}				
				} else {
					$output_string = $output_string."\t"."Inconsistent";
					
					######################################################################
					######################################################################
					####   hard bit											 #############
					######################################################################
					######################################################################
					
					if ($a_ref ne "N/A" && $b_ref eq "N/A" ) {
						if ($geno_homo_hetero[$i] eq "AA"){
							$false_positives[$i] += 1;
						} else {#BB
							if ($a_ref  eq "A=Ref") {
								if (exists $most_covered_bases[$i]{$g_b}) {
									$true_positives[$i] += 1;
								}else{
									$false_positives[$i] += 1;
								}
							} else {
								if (exists $most_covered_bases[$i]{ ComplementaryBase($g_b)}) {
									$true_positives[$i] += 1;
								} else {
									$false_positives[$i] += 1;
								}
							}
						}
					} elsif ( $a_ref eq "N/A" && $b_ref ne "N/A" ){
						if  ($geno_homo_hetero[$i] eq "AA") { #AA
							if ($b_ref  eq "B=Ref"){
								if (exists $most_covered_bases[$i]{$g_a}) {
									$true_positives[$i] += 1;
								}else{
									$false_positives[$i] += 1;
								}
							} else {
								if (exists $most_covered_bases[$i]{ ComplementaryBase($g_a)}) {
									$true_positives[$i] += 1;
								}else{
									$false_positives[$i] += 1;
								}
							}
						} else { #BB
							$false_positives[$i] += 1;
						}
					} else { # if both Allele A and Allele B are NON related to ref or complemtary ref
						if ($a_ref  eq "A=Ref") {
							if ($geno_homo_hetero[$i] eq "AA") {
								$false_positives[$i] += 1;
							} else {
								if (exists $most_covered_bases[$i]{$g_b}) {
									$true_positives[$i] += 1;
								}else{
									$false_positives[$i] += 1;
								}
							}
						} elsif ($b_ref  eq "B=Ref") {
							if ($geno_homo_hetero[$i] eq "BB") {
								$false_positives[$i] += 1;
							} else {
								if (exists $most_covered_bases[$i]{$g_a}) {
									$true_positives[$i] += 1;
								}else{
									$false_positives[$i] += 1;
								}
							}
						}elsif( exists $most_covered_bases[$i]{$g_a} || exists $most_covered_bases[$i]{$g_b}) {
							$true_positives[$i] += 1;
						} else {
							$false_positives[$i] += 1;
						}
					}
				}
			}
		}
		
		print OUTPUT $output_string."\n";
	}
	
	close INPUT;
	close OUTPUT;
	
	for(my $j =0 ; $j <3; $j++) {
		print "For patient $j :\n";
		print "\tTrue Positives = ";
		print $true_positives[$j];
		print "\n";
		print "\tFalse Positives = ";
		print $false_positives[$j];
		print "\n";
		print "\tTrue Negatives = ";
		print $true_negatives[$j];
		print "\n";
		print "\tFalse Negatives = ";
		print $false_negatives[$j];
		print "\n";
		print "\n";
		print "\tSensitivity of Exome sequencing is ";
		print $true_positives[$j]/($true_positives[$j] + $false_negatives[$j]);
		print "\n";
		print "\tSpecificity is ";
		print $true_negatives[$j]/($true_negatives[$j] + $false_positives[$j]);
		print "\n";
		print "\n";
		print "\n";
	}
	exit;
}



sub ComplementaryBase {
	my($x) = @_;
	my $c_x = $x;
	$c_x =~ tr/[A,a,T,t,G,g,C,c]/[T,t,A,a,C,c,G,g]/;
	return $c_x;
}

sub CompareAwithB {
	my ($x,$y) = @_;
	if ($x eq  $y) {
		return "A=B";
	}elsif ($x eq ComplementaryBase($y)) {
		return "A=C_B";
	}else {
		return "No Match";
	}
}


