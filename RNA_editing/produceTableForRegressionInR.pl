#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
#use DBI;
#use DBD::mysql;
#use LWP::Simple;
#use Getopt::Long;


my ($folder) = @ARGV;

my $common_heterozygous = Looking_4_common_hereozygous_in_VCF_file_folder_4_regression($folder);

exit;


sub Looking_4_common_hereozygous_in_VCF_file_folder_4_regression {
	
	my ($folder) = @_;

	my $dir = cwd();
	$folder =~ s/\/$//;
	my $files = getFilesInFlowCellLaneOrder($folder);
	my @files = @$files;
	
	my @sample_names;
	my %variant_hash;
	
	my $depth_filter =  0; #locus with seqencing depth lower than the value will be dropped.
	my $gq_quality_filter = 13; # Phred-scaled p-value of the genotype calling is wrong
	
	my @file_list_4_mpileup;
	my %ref_hash;
	
	foreach my $file (@files) {
		next if ($file =~ m/^\./);
		if ($file =~ m/\.flt.vcf$/){ ## This line controls the name pattern of files that are read in
			my $input_file = $folder."/".$file;
			
			my $temp_sample_name = $file;
			$temp_sample_name =~ s/\.var\.flt\.vcf//;
			push @sample_names,$temp_sample_name;
			
			open INPUT, $input_file or die "cannot open file $input_file.\n";
			push @file_list_4_mpileup, $file;
			print "reading ".$input_file."\n";
			while (my $line =<INPUT>) {
				chomp $line;
				if( $line =~ m/^chr/ ){
					my @words = split("\t", $line);
					my $key = $words[0]."_".$words[1];
					
					##filters: include Sequencing depth filter and GenoType call quality filter
					my $depth = $words[5];
					my @vcf_sample = split (":", $words[9]);
					my $gq_quality = $vcf_sample[3];
					
					if($depth >= $depth_filter && $gq_quality >= $gq_quality_filter) {
						#if($vcf_sample[0] =~ m/0\/1/){
							my @chars = split(//, $words[3]);
							if(scalar @chars == 1){ # not to include indel but if a single change variant is detected in other sample
												# the position is still included. coverage of the position on indel samples will be 
												# extracted from pileup file							
								if(! exists $ref_hash{$key}) {
									$ref_hash{$key} = $words[3];
								}	
							
								if (exists $variant_hash{$key}) {
									$variant_hash{$key} += 1;
								} else {
									$variant_hash{$key} = 1;
								}	
							}
						#}

					}
				} 
			}
			close INPUT;
		}
	}
	
	print "\n";
	print "size of variant hash:  " . keys( %variant_hash ) . ".\n";
	print "\n";
	
	my $stupid_count = 0;
	for my $key ( keys %variant_hash ) {
		if ($variant_hash{$key}  == 16) { 		#select those variants that exists in any sample
			$stupid_count += 1;
		}
	}
	print "\n";
	print "size of common variants  existing in all samples:  " . $stupid_count . ".\n";
	print "\n";

	print "\n";
	print "size of ref hash:  " . keys( %ref_hash ) . ".\n";
	print "\n";
	
	
	print "################################################################\n";
	print "#                     Asign SNP to genes                       #\n";
	print "#REF/Var could be biased, check the script it self for detail./#\n";
	print "################################################################\n";
	print "\n";
	#print "######\n";
	#print "Warning: Files will only have the ref and var bases extracted from the last vcf file processed. God blessed all samples have the same ref and var!\n";
	#print "######\n";
	print "\n";
	
	my %annoted_hash;
	my $temp_output_file = $dir."/"."temptemptemp_file_for_bedtools.bed";
	my $temp_output_file_2 = $dir."/"."temptemptemp_file_with_genes_Names.txt";
	
	open OUTPUT, ">$temp_output_file" or die "cannot open file $temp_output_file.\n";
	for my $key ( keys %variant_hash ) {
		if ($variant_hash{$key} >= 1) { 		#select those variants that exists in any sample
			my @words = split("_", $key);
			my $bed_start = $words[1]-1; #bed format is 0-numbered
			print OUTPUT $words[0]."\t".$bed_start."\t".$words[1]."\n";
		}
	}
	close OUTPUT;
	
	print "Executing command: /users/a5907529/biosofts/BEDTools-Version-2.12.0/bin/intersectBed -wo -a $temp_output_file -b ~/lustre/Yaobo/GenomeData/gencode.v10.annotation.gtf > $temp_output_file_2\n";
	`intersectBed -wo -a $temp_output_file -b ~/lustre/Yaobo/GenomeData/gencode.v10.annotation.gtf > $temp_output_file_2`;
	
	open INPUT, $temp_output_file_2 or die "cannot open file $temp_output_file_2.\n";
	while (my $line =<INPUT>) {
		chomp $line;
		if( $line =~ m/^chr/ ){
			my @words = split("\t", $line);
			my $gene_name;
			my $gene_id;
			my $key = $words[0]."_".$words[2];
			if (exists $words[11]) {
				my @attributes = split ("; ", $words[11]);
				#print "word 11 are::";
				#print @attributes;
				#print "\n";
				
				foreach my $attribute (@attributes) {
					if ($attribute =~ m/^gene\_id/){
						$attribute =~ s/gene\_id//;
						$attribute =~ s/\"//g;
						$attribute =~ s/\s//g;
						#print "subtitued gene_id is:";
						#print $attribute;
						#print "\n";
						$gene_id = $attribute;
					}
					
					if ($attribute =~ m/^gene\_name/){
						$attribute =~ s/gene\_name//;
						$attribute =~ s/\"//g;
						$attribute =~ s/\s//g;
						#print "subtitued gene_name is:";
						#print $attribute;
						#print "\n";
						$gene_name = $attribute;
					}
				}
				$annoted_hash{$key} = $gene_name."_".$gene_id;	
			}
		}
	}
	close INPUT;
	
	#print "annotaed_hash length after bedtools is ";
	#$size = keys %annoted_hash;
	#print $size;
	#print "\n";
	
	
	#`rm $temp_output_file`;
	#`rm $temp_output_file_2`;
	
	for my $key ( keys %variant_hash ) {  #select those variants that is heterzygous in at least 1 sample
		if ($variant_hash{$key} >= 1 &&  ! exists $annoted_hash{$key} ) {
				$annoted_hash{$key} = "NA"."_"."NA";
		}
	}
	
	print "\n";
	print "size of Annoted hash:  " . keys( %annoted_hash ) . ".\n";
	print "\n";
	
	#print "annotaed_hash length after NNAA is ";
	#$size = keys %annoted_hash;
	#print $size;
	#print "\n";
	
	##Counting Genotype#
	
	my %genoType_hash;
	#my %ref_var_hash; # this hash will only contain the ref and var from the last vcf file processed. God blessed all samples have the same ref and var!"
	my %homo_count;
	my %hetre_count;
	
	my $file_counts = 0;
	#my %export_2_R;
	
	foreach my $file (@files) {
		next if ($file =~ m/^\./);
		if ($file =~ m/\.flt.vcf$/){ ## This line controls the name pattern of files that are read in
			$file_counts += 1;
			my $input_file = $folder."/".$file;
			open INPUT, $input_file or die "cannot open file $input_file.\n";
			
			print "reading ".$input_file."\n";
			while (my $line =<INPUT>) {
				chomp $line;
				if( $line =~ m/^chr/ ){
					my @words = split("\t", $line);
					my $key = $words[0]."_".$words[1];
					my $key_sample = $words[0]."_".$words[1]."_".$file;
					#$ref_var_hash{$key} = $words[3]."/".$words[4];
					if (exists $annoted_hash{$key}) {
						my $ref_base = $words[3];
						my $var_base = $words[4];
						my @vcf_info = split (";", $words[7]);
						my @vcf_sample = split (":", $words[9]);
						my $ref2var_coverage;
						
						##### filters ###### 
						my $depth = $words[5];
						my $gq_quality = $vcf_sample[3];
						##### filters end #####
						
						if($depth >= $depth_filter && $gq_quality >= $gq_quality_filter) {
						# quality of the variant on this position must pass the filter
						# otherwise consider there no variant and only ref base coverage
						# thus these positions should be passed to mpilup
						
							foreach my $tag(@vcf_info) { # it requires that the vcf file info column has to have a DP4 tag
								if ($tag =~ m/^DP4/) {
									$tag =~ s/^DP4\=//;
									my @numbers = split (",", $tag);
									my $temp1 = $numbers[0] + $numbers[1]; #extract coverage of ref and revers_ref
									my $temp2 = $numbers[2] + $numbers[3]; #extract coverage of var and revers_var
									$ref2var_coverage = $temp1.":".$temp2; #recorde the ratio
									
									my @chars = split(",",$var_base);
									$var_base = $chars[0];
									
									# to resolve multi alternative bases/Variants
									if (scalar @chars == 2 && $vcf_sample[0] =~ m/\d\/2/) {
											$var_base = $chars[1];
									} elsif (scalar @chars == 3 && $vcf_sample[0] =~ m/\d\/2/){
										$var_base = $chars[1];
									} elsif (scalar @chars == 3 && $vcf_sample[0] =~ m/\d\/3/) {
										$var_base = $chars[2];
									}
									
									
									if (!exists $genoType_hash{$key_sample}) {
										$genoType_hash{$key_sample} =$ref_base."_".$var_base."_".$vcf_sample[0]."_".$ref2var_coverage; #genotye (eg 0/1 or 1/1) and base ratio
										#print $key_sample." is the wow ".$genoType_hash{$key_sample}."\n"; # test
									} else {
										print "Error 01: duplicated entry in $file for $key.\n";
										exit;
									}
								} else {
									next;
								}
							}
						} else {
							next;
						}
						
						if (!exists $hetre_count{$key}) {# initial the count to allow value 0
							$hetre_count{$key} = 0;
						}
						if (!exists $homo_count{$key}) {
							$homo_count{$key} = 0;
						}
						
						if ($vcf_sample[0] =~ m/0\/1/) {# start to count
							$hetre_count{$key} += 1;
						} else {
							$homo_count{$key} += 1;						
						}
						
					} else {
						next;
					}
				} 
			}
			
			close INPUT;
			
			#print "size of geneType hash b4 pileup for $file:  " . keys( %genoType_hash ) . ".\n";
			
			#### find those variants in annoted_hash but not in genoType_hash yet, record them in an array
			my @pileup_on_positions;
			$file =~ m/^f\_(\d)\.s\_(\d).*flt\.vcf$/;
			print "the file is $file\n";
			my $bam_file = $folder."/".'f_'.$1.".".'s_'.$2.'.gsnap20120620.sorted.bam';
			
			for my $key (keys %annoted_hash) {
				my $key_sample = $key."_".$file;
				if(!exists $genoType_hash{$key_sample}) {
					push @pileup_on_positions,$key;
				}
			}
			
			my $ref_pileup_hash = getCoverageForPositions(\@pileup_on_positions, $bam_file);
			
			my %ref_pileup_hash = %$ref_pileup_hash;
			
			for my $key (keys %annoted_hash) { #find those variants that not recored in this vcf file but in annoted_hash, assign an value
				my $key_sample = $key."_".$file;
				if(!exists $genoType_hash{$key_sample}) {

					my $ref_count = $ref_pileup_hash{$key};
					my $ref_base = $ref_hash{$key};
					$genoType_hash{$key_sample} = $ref_base."_".'-'."_".'0/0_'."$ref_count".':0';

				} 
			}
		}
	}
	
	my %OA_group = ( 1 => 1, 3 => 1, 4 => 1, 6 => 1, 8 => 1, 9 => 1,  11 => 1, 12 => 1, 15 => 1, 16 => 1);
	my %NOF_group = ( 2 => 1, 5 => 1, 7 => 1, 10 => 1, 13 => 1, 14 => 1);
	my @output;
	my @output_hetereo_freq;
	
	my $total_count = 0;
	my $hetero_in_NOF_only = 0;
	my $hetero_in_OA_only = 0;

	my %tricky_variant_positions;
	my $tricky_variant_positions_count;
	
	foreach my $key ( keys %annoted_hash ) {
		my @words = split("_", $key);
		if ($words[0] =~ m/chrY/){
			next;
		} else {
			my $sample_count = 0;
			my $OA_hetero_count = 0;
			my $NOF_hetero_count = 0;
			my $oa_freq_string = "";
			my $nof_freq_string = "";
			my $oa_ref_freq_sum = 0;
			my $oa_var_freq_sum = 0;
			my $nof_ref_freq_sum = 0;
			my $nof_var_freq_sum = 0;
			my $disease_type = "";
			
			my %var_hash_for_samples_on_one_pos;
			
			############################ test if the var is consistent with the majortity change on the position 
			foreach my $sample (@file_list_4_mpileup){
				my $key_sample = $key."_".$sample;
				my @temp = split("_", $genoType_hash{$key_sample} );
				$var_hash_for_samples_on_one_pos{$key_sample} = $temp[1];
			}
			
			my %var_count_hash;
			
			while ( my ($key_temp, $value) = each %var_hash_for_samples_on_one_pos ) {
				if ($value =~ m/\w/) { # to exclude the - 
					if (! exists $var_count_hash{$value}) {
						$var_count_hash{$value} = 1;
					} else {
						$var_count_hash{$value} += 1;
					}
				}
			}
			
			my $majority_var;
			my $var_freq_sum;
			my $most_frequent_var;

			for my $key_temp (keys %var_count_hash) {
				$var_freq_sum += $var_count_hash{$key_temp};
			}
			
			my $max_temp = 0;
			for my $key_temp (keys %var_count_hash) {
				if ($var_count_hash{$key_temp} > $max_temp) {
					$max_temp = $var_count_hash{$key_temp};
					$most_frequent_var = $key_temp;
				}
			}
			
			for my $key_temp (keys %var_count_hash) {
				if ($var_count_hash{$most_frequent_var} > $var_freq_sum/2) {
					$majority_var = $most_frequent_var;	
				}
			}
			
			if(defined $majority_var){
				for my $key_temp (keys %var_hash_for_samples_on_one_pos) {
					if ($var_hash_for_samples_on_one_pos{$key_temp} =~ m/\w/ && $var_hash_for_samples_on_one_pos{$key_temp} ne $majority_var) {
						my @temp = split("_", $genoType_hash{$key_temp} );
						$genoType_hash{$key_temp} = $temp[0]."_".$temp[1]."_".$temp[2]."_".'-:-';
					}
				}
			} else {
				$tricky_variant_positions_count += 1;
				$tricky_variant_positions{$key} = 1;
			}
			
			foreach my $sample (@file_list_4_mpileup){
				$sample_count += 1;
				my $key_sample = $key."_".$sample;
				
				#print "dfghdfgdf   ".$key_sample."                              ".$genoType_hash{$key_sample}."\n";
				my @temp = split("_", $genoType_hash{$key_sample} );
				my @numbers = split(":", $temp[3]);
				
				if($temp[2] =~ m/0\/1/) {
					if(exists($OA_group{$sample_count})){
						$OA_hetero_count += 1;
					} elsif (exists($NOF_group{$sample_count})) {
						$NOF_hetero_count += 1;
					}else {
						print "Fatal error, Grouping failed!"."\n";
						exit;
					}
				}
				
				if(exists($OA_group{$sample_count})){
					$disease_type = "OA";
					#if ( $oa_freq_string eq "") {
					#$oa_freq_string = $temp[3];
					#} else {
					#	$oa_freq_string = $oa_freq_string.";".$temp[3];
					#}
					#$oa_ref_freq_sum = $oa_ref_freq_sum + $numbers[0];
					#$oa_var_freq_sum = $oa_var_freq_sum + $numbers[1];
				} elsif (exists($NOF_group{$sample_count})) {
					$disease_type = "NOF";
					#$NOF_hetero_count += 1;
					### this is to generate numbers for ref and var frequences of OA and NOF samples seperatedly
					#if($nof_freq_string eq "") {
					#	$nof_freq_string = $temp[3];
					#} else {
					#	$nof_freq_string = $nof_freq_string.";".$temp[3];
					#}
					#$nof_ref_freq_sum = $nof_ref_freq_sum + $numbers[0];
					#$nof_var_freq_sum = $nof_var_freq_sum + $numbers[1];
					###
				}else {
					print "Fatal error, Grouping failed!"."\n";
					exit;
				}
				
				#####
				# here to output result
				#####
				my @temp2 = split("_", $annoted_hash{$key});
	            #"Chromosome", "Postition", "Postition", "Ref", "Var", "Symbol","ENSEMBL_ID", "Sample_name", "Disease", "cRef", "cVar", "Genotype";
				my $output_string = $words[0]."\t".$words[1]."\t".$words[1]."\t".$temp[0]."\t".$temp[1]."\t".$temp2[0]."\t".$temp2[1]."\t".$sample."\t".$disease_type."\t".$numbers[0]."\t".$numbers[1]."\t".$temp[2]."\n";
				push @output,$output_string;
			}
			
			if($OA_hetero_count == 0 && $NOF_hetero_count != 0) {
				$hetero_in_NOF_only += 1;
			} elsif($OA_hetero_count != 0 && $NOF_hetero_count == 0) {
				$hetero_in_OA_only += 1;
			}
		}
	}
	
	print "###########################################\n";
	print "########### output results ################\n";
	print "###########################################\n";
	
		
	my $output_file = $dir."/"."Common_Heterozygous_with_info_filtered_4_regression_newVersion.txt";
	my $column_header = "Chromosome"."\t"."Start"."\t"."End"."\t"."Ref"."\t"."Var"."\t"."Symbol"."\t"."ENSEMBL_ID"."\t"."Sample_name"."\t"."Disease"."\t"."cRef"."\t"."cVar"."\t"."Genotype";
	open OUTPUT, ">$output_file" or die "cannot open file $output_file.\n";
	print OUTPUT "#Total qulified common heterozygous are: $total_count"."\n";
	print OUTPUT "#Heterozygous in NOF only: $hetero_in_NOF_only"."\n";
	print OUTPUT "#Heterozygous in OA only: $hetero_in_OA_only"."\n";
	print OUTPUT $column_header."\n";
	print OUTPUT @output; 
	close OUTPUT;
	
	my $output_file_2 = $dir."/"."tricky_variant_positions.txt";
	open OUTPUT1, ">$output_file_2" or die "Cannot open $output_file_2";
	for my $key (keys %tricky_variant_positions) {
		print OUTPUT1 $key."\n";
	}
	close OUTPUT1;
	
	print "Total qulified common heterozygous are: $total_count"."\n";
	print "Heterozygous in NOF only: $hetero_in_NOF_only"."\n";
	print "Heterozygous in OA only: $hetero_in_OA_only"."\n";
	print "number of tricky variant positions: $tricky_variant_positions_count \n";
	print "Result is in: $output_file\n";
	
	
	return $output_file;
}


sub getCoverageForPositions {
	my ($positions,$bam_file) = @_;
	my $coverage;
	print "Using mpileup to get coverage for $bam_file.\n";


	# To test the script
	my $name = $bam_file;
	$name =~ s/\W/\_/g;
	#
	
	my $temp_mpileup_output = "temp_temp_mpileup_$name.pileup";
	my $temp_mpileup_position_bed = "temp_temp_mpileup_$name.bed";
	
	my @positions = @$positions;
	
	open BED, ">$temp_mpileup_position_bed" or die "Cannot open $temp_mpileup_position_bed.";
	foreach my $key (@positions) {
		my @coordinates = split("_", $key);
		my $chr = $coordinates[0];
		my $pos = $coordinates[1];
		print BED $chr."\t".$pos."\n";
	}
	close BED; 
	print "Executing command: ~/biosofts/samtools-0.1.18/samtools mpileup -f ~/lustre/Yaobo/GenomeData/hg19.fa -l $temp_mpileup_position_bed $bam_file > $temp_mpileup_output\n";
	
	`samtools mpileup -f ~/lustre/Yaobo/GenomeData/hg19.fa -l $temp_mpileup_position_bed $bam_file > $temp_mpileup_output`;
	
	print "mpileup for $bam_file is done!.\n";
	
	my %result_hash;
		
	open PILEUP, $temp_mpileup_output or die "Cannot open $temp_mpileup_output.";
	while (my $line = <PILEUP>){ 
		chomp $line;
		my @words = split("\t", $line);
		my $key = $words[0]."_".$words[1];
		if (! exists $result_hash{$key}) {
			#The following 7 lines are to ignore the coverage caused by reference skips in pileup file
			my @chars = split(//, $words[4]);
			my $count = 0;
			foreach my $char (@chars) {
				if ($char eq '>' | $char eq '<') {
					$count += 1;
				}
			}
			$result_hash{$key} = $words[3]-$count;
		} else {
			print "It's wrong wrong wrong!\n";	
		}
	}
	close PILEUP;
	
	foreach my $key (@positions) {
		if(! exists $result_hash{$key}){
			$result_hash{$key} = 0;
		}
	}
	
	return \%result_hash;
}

sub getFilesInFlowCellLaneOrder {
	my ($folder) = @_;
	$folder =~ s/\/$//;
	#my $dir = cwd();
	my %file_list;
	my @ordered_list;
	
	opendir (DIR, $folder) or die "Cann't find the directory";
	while (my $file = readdir(DIR)) {
		next if ($file =~ m/^\./);
		if ($file =~ m/\.flt.vcf$/){
			$file =~ m/^f\_(\d)\.s\_(\d).*flt\.vcf$/;
			$file_list{$file} = $1."_".$2;
		}
	}
	
	while (keys %file_list){
		my @loop_array = keys %file_list;	
		#print "new loop with ".$loop_array[0]."\n";
		my $smaller_key = $loop_array[0];
		my @temp_1 = split("_", $file_list{$smaller_key});
		
		while (my ($key_1, $value) = each %file_list){
			my @temp_2 = split ("_", $value);
			if($temp_1[0] == $temp_2[0]){
				if($temp_1[1] > $temp_2[1]){
					$smaller_key = $key_1;
					@temp_1 = split("_", $file_list{$smaller_key});
				}				
			} elsif ($temp_1[0] > $temp_2[0]) {
					$smaller_key = $key_1;
					@temp_1 = split("_", $file_list{$smaller_key});
			}
		}
		
		push @ordered_list,$smaller_key;
		#print "delete ".$smaller_key."\n";
		delete $file_list{$smaller_key};
	}
	
	return \@ordered_list;
}