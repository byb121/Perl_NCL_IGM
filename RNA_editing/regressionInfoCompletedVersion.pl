#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
#use DBI;
#use DBD::mysql;
use LWP::Simple;
use Getopt::Long;


my ($folder) = @ARGV;

my $common_heterozygous = Looking_4_common_hereozygous_in_VCF_file_folder_4_regression($folder);

exit;


sub Looking_4_common_hereozygous_in_VCF_file_folder_4_regression {
	
	my ($folder) = @_;
	$folder =~ s/\/$//;
	my $dir = cwd();
	
	my @sample_names;
	my %hetero_hash;
	
	my $depth_filter = 10; #locus with seqencing depth lower than the value will be dropped.
	my $alt_quality_filter = 13; # Phred-scaled p-value of the alternative base calling is wrong
	my $gq_quality_filter = 13; # Phred-scaled p-value of the genotype calling is wrong
	
	my @file_list_4_mpileup;
	my %ref_hash;
	
	opendir (DIR, $folder) or die "Cann't find the directory";
	while (my $file = readdir(DIR)) {
		next if ($file =~ m/^\./);
		if ($file =~ m/\.flt.vcf$/){ ## This line controls the name pattern of files that are read in
			my $input_file = $folder."/".$file;
			
			my $temp_sample_name = $file;
			$temp_sample_name =~ s/\.var\.flt\.vcf//;
			push @sample_names,$temp_sample_name;
			
			open INPUT, $input_file or die "cannot open file $input_file.\n";
			push @file_list_4_mpileup, $input_file."\n";
			print "reading ".$input_file."\n";
			while (my $line =<INPUT>) {
				chomp $line;
				if( $line =~ m/^chr/ ){
					my @words = split("\t", $line);
					my $key = $words[0]."_".$words[1];
					
					##filters: include Sequencing depth filter, ALT quaility filter and GenoType call quality filter
					my $depth = $words[5];
					my $alt_quality;
					my @vcf_info = split (";", $words[7]);
					foreach my $tag(@vcf_info) {
						if ($tag =~ m/^DP\=/) {
							$tag =~ s/^DP\=//;
							$alt_quality = $tag;
						}
					}
					my @vcf_sample = split (":", $words[9]);
					my $gq_quality = $vcf_sample[2];
					
					if($depth >= $depth_filter && $alt_quality >= $alt_quality_filter && $gq_quality >= $gq_quality_filter) {
						if($vcf_sample[0] =~ m/0\/1/){
							if (exists $hetero_hash{$key}) {
								$hetero_hash{$key} += 1;
							} else {
								$hetero_hash{$key} = 1;
							}	
						}
						
						if(! exists $ref_hash{$key}) {
							$ref_hash{$key} = $words[3];
						}
					}
				} 
			}
			close INPUT;
		}
	}
	closedir(DIR);
	
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
	for my $key ( keys %hetero_hash ) {
		if ($hetero_hash{$key} >= 1) { 		#select those variants that exists in any sample
			my @words = split("_", $key);
			my $bed_start = $words[1]-1; #bed format is 0-numbered
			print OUTPUT $words[0]."\t".$bed_start."\t".$words[1]."\n";
		}
	}
	close OUTPUT;
	
	print "Executing command: /users/a5907529/biosofts/BEDTools-Version-2.12.0/bin/intersectBed -wo -a $temp_output_file -b ~/GenomeData/hg19/gencode.v10.annotation.gtf > $temp_output_file_2\n";
	`~/biosofts/BEDTools-Version-2.12.0/bin/intersectBed -wo -a $temp_output_file -b ~/GenomeData/hg19/gencode.v10.annotation.gtf > $temp_output_file_2`;
	
	open INPUT, $temp_output_file_2 or die "cannot open file $temp_output_file_2.\n";
	while (my $line =<INPUT>) {
		chomp $line;
		if( $line =~ m/^chr/ ){
			my @words = split("\t", $line);
			my $gene_name;
			my $gene_id;
			my $key = $words[0]."_".$words[1];
			if (exists $words[11]) {
				my @attributes = split ("; ", $words[11]);
				#print "word 11 are::";
				#print @attributes;
				#print "\n";
				
				foreach my $attribute (@attributes) {
					if ($attribute =~ m/^gene\_id/){
						$attribute =~ s/gene\_id//;
						$attribute =~ s/\"//g;
						#print "subtitued gene_id is:";
						#print $attribute;
						#print "\n";
						$gene_id = $attribute;
					}
					
					if ($attribute =~ m/^gene\_name/){
						$attribute =~ s/gene\_name//;
						$attribute =~ s/\"//g;
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
	
	#`rm $temp_output_file`;
	#`rm $temp_output_file_2`;
	
	for my $key ( keys %hetero_hash ) {
		if ($hetero_hash{$key} >= 1 && ! exists $annoted_hash{$key} ) { 	#select those variants that is heterzygous in at least 1 sample
			$annoted_hash{$key} = "NA"."_"."NA";
		}
	}
	
	
	
	##Counting Genotype#
	
	my %genoType_hash;
	#my %ref_var_hash; # this hash will only contain the ref and var from the last vcf file processed. God blessed all samples have the same ref and var!"
	my %homo_count;
	my %hetre_count;
	
	my $file_counts = 0;
	#my %export_2_R;
	
	opendir (DIR, $folder) or die "Cann't find the directory";
	while (my $file = readdir(DIR)) {
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
						my $alt_quality;
						foreach my $tag(@vcf_info) {
							if ($tag =~ m/^DP\=/) {
								$tag =~ s/^DP\=//;
								$alt_quality = $tag;
							}
						}
						my $gq_quality = $vcf_sample[2];
						##### filters end #####
						
						if($depth >= $depth_filter && $alt_quality >= $alt_quality_filter && $gq_quality >= $gq_quality_filter) {
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
									if (!exists $genoType_hash{$key_sample}) {
										$genoType_hash{$key_sample} =$ref_base."_".$var_base."_".$vcf_sample[0]."_".$ref2var_coverage; #genotye (eg 0/1 or 1/1) and base ratio
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
			#### find those variants in annoted_hash but not in genoType_hash yet, record them in an array
			my @pileup_on_positions;
			$file =~ m/^f\_(\d)\.s\_(\d).*flt\.vcf$/;
			print "the file is $file\n";
			my $bam_file = '/users/a5907529/RNA_seq/gsnap_no_snp_align/f_'.$1.".".'s_'.$2.'.gsnap20120620.sorted.bam';
			
			for my $key (keys %annoted_hash) {
				my $key_sample = $key."_".$file;
				if(!exists $genoType_hash{$key_sample}) {
					push @pileup_on_positions,$key;
				}
			}
			
			#### use the array as positions to generate pileup file
			my $ref_pileup_hash = getCoverageForPositions(\@pileup_on_positions, $bam_file);
			
			my %ref_pileup_hash = %$ref_pileup_hash;
			
			for my $key (keys %annoted_hash) { #find those variants that not recored in this vcf file but in annoted_hash, assign an value
				my $key_sample = $key."_".$file;
				if(!exists $genoType_hash{$key_sample}) {
					################################################################
					print "Check now.\n";
					print $key."\n";
					my $ref_count = $ref_pileup_hash{$key};
					my $ref_base = $ref_hash{$key};
					print $ref_count."\n";
					print $ref_base."\n";
					print "Check oh go.\n";
					$genoType_hash{$key_sample} = $ref_base."_".'-'."_".'0/0_'."$ref_count".':0';####error
					################################################################
				} #else {
					#print "Error: don't know why on $key_sample.\n"; ####error
				#}
			}
			
			
		}
	}
	closedir(DIR);
	
	#########################################################################################################
	# need to adapt this part
	########################## loop the keys in annotated_hash and loop the files in @file will have the same thing
	
	my %OA_group = ( 1 => 1, 3 => 1, 4 => 1, 6 => 1, 8 => 1, 9 => 1,  11 => 1, 12 => 1, 15 => 1, 16 => 1);
	my %NOF_group = ( 2 => 1, 5 => 1, 7 => 1, 10 => 1, 13 => 1, 14 => 1);
	my @output;
	my @output_hetereo_freq;
	
	my $total_count = 0;
	my $hetero_in_NOF_only = 0;
	my $hetero_in_OA_only = 0;
	
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
			
			foreach my $sample (@file_list_4_mpileup){
				$sample_count += 1;
				my $key_sample = $key."_".$sample;
				my @temp = split("_", $genoType_hash{$key_sample});
				my @numbers = split(":", $temp[3]);
				if($temp[2] =~ m/0\/1/) {
					if(exists($OA_group{$sample_count})){
						$disease_type = "OA";
						$OA_hetero_count += 1;
						if ( $oa_freq_string eq "") {
							$oa_freq_string = $temp[3];
						} else {
							$oa_freq_string = $oa_freq_string.";".$temp[3];
						}
						$oa_ref_freq_sum = $oa_ref_freq_sum + $numbers[0];
						$oa_var_freq_sum = $oa_var_freq_sum + $numbers[1];
					} elsif (exists($NOF_group{$sample_count})) {
						$disease_type = "NOF";
						$NOF_hetero_count += 1;
						### this is to generate numbers for ref and var frequences of OA and NOF samples seperatedly
						if($nof_freq_string eq "") {
							$nof_freq_string = $temp[3];
						} else {
							$nof_freq_string = $nof_freq_string.";".$temp[3];
						}
						$nof_ref_freq_sum = $nof_ref_freq_sum + $numbers[0];
						$nof_var_freq_sum = $nof_var_freq_sum + $numbers[1];
						###
					}else {
						print "Fatal error, Grouping failed!"."\n";
						exit;
					}
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
	
	print "Total qulified common heterozygous are: $total_count"."\n";
	print "Heterozygous in NOF only: $hetero_in_NOF_only"."\n";
	print "Heterozygous in OA only: $hetero_in_OA_only"."\n";
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
	print "Executing command: ~/biosofts/samtools-0.1.18/samtools mpileup -f ~/GenomeData/hg19/hg19.fa -l $temp_mpileup_position_bed $bam_file > $temp_mpileup_output\n";
	
	`~/biosofts/samtools-0.1.18/samtools mpileup -f ~/GenomeData/hg19/hg19.fa -l $temp_mpileup_position_bed $bam_file > $temp_mpileup_output`;
	
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

	