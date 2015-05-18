#!/usr/bin/perl
use strict;
use warnings;

my @genotyping_files = (
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_1_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_2_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_3_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_4_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_5_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_6_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_7_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_8_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_9_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_10_Full_Data_Table.txt_expanded.txt", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_11_Full_Data_Table.txt_expanded.txt",
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_New_Samples_Full_Data_Table.txt_expanded.txt",
"/home/a5907529/WORKING_DATA/Allelic_Exp/GenoChipData/Test/A2311_Test_Full_Data_Table.txt_expanded.txt"
);

my @vcf_files = (
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr1_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr2_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr3_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr4_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr5_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr6_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr7_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr8_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr9_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr10_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr11_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr12_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr13_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr14_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr15_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr16_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr17_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr18_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr19_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr20_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr21_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chr22_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chrX_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chrY_Freebayes.vcf", 
"/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/Free_out_refined/Yaobo_RNAVariants_AA_refinedBAM_chrM_Freebayes.vcf" 
);

my $output_file = "/home/a5907529/WORKING_DATA/Allelic_Exp/Freebayes_RANVarints/processed_genotyping_vcf.txt";
my @require_samples = (
"HDBR375", "HDBR380", "HDBR383", "HDBR701", "HDBR702", "HDBR703", "HDBR705", "HDBR706", "HDBR707", "HDBR708", "HDBR709",
"HDBR710", "HDBR712", "HDBR713", "HDBR714", "HDBR715", "HDBR716", "HDBR717", "HDBR718", "HDBR720", "HDBR721", "HDBR722",
"HDBR723", "HDBR724", "HDBR725", "HDBR726", "HDBR727", "HDBR728", "HDBR729", "HDBR731", "HDBR732", "HDBR733", "HDBR734",
"HDBR735", "HDBR736", "HDBR737", "HDBR738", "HDBR739", "HDBR740", "HDBR741", "HDBR742", "HDBR743", "HDBR744", "HDBR745",
"HDBR746", "HDBR747", "HDBR748", "HDBR749", "HDBR750", "HDBR751", "HDBR752", "HDBR753", "HDBR754", "HDBR755", "HDBR756",
"HDBR757", "HDBR758", "HDBR759", "HDBR760", "HDBR762", "HDBR763", "HDBR765", "HDBR766", "HDBR767", "HDBR768", "HDBR769",
"HDBR770", "HDBR771", "HDBR772", "HDBR774", "HDBR775", "HDBR776", "HDBR777", "HDBR778", "HDBR779", "HDBR780", "HDBR782",
"HDBR783", "HDBR784", "HDBR785", "HDBR786", "HDBR787", "HDBR788", "HDBR790", "HDBR791", "HDBR792", "HDBR793", "HDBR794",
"HDBR795", "HDBR796", "HDBR797", "HDBR798", "HDBR799", "HDBR802", "HDBR803", "HDBR805", "HDBR806", "HDBR807", "HDBR809",
"HDBR810", "HDBR811", "HDBR812", "HDBR813", "HDBR814", "HDBR815", "HDBR816", "HDBR818", "HDBR819", "HDBR820", "HDBR821",
"HDBR822", "HDBR823", "HDBR824", "HDBR825", "HDBR826", "HDBR827", "HDBR828", "HDBR829", "HDBR831", "HDBR832", "HDBR837",
"HDBR838", "HDBR839", "HDBR840", "HDBR841", "HDBR842", "HDBR843", "HDBR844", "HDBR845", "HDBR846", "HDBR848", "HDBR849",
"HDBR850", "HDBR851", "HDBR852", "HDBR853", "HDBR854", "HDBR855", "HDBR856", "HDBR857", "HDBR858", "HDBR860", "HDBR861",
"HDBR862", "HDBR863", "HDBR864", "HDBR866", "HDBR867", "HDBR868", "HDBR869", "HDBR870", "HDBR871", "HDBR872", "HDBR873",
"HDBR874", "HDBR876", "HDBR877", "HDBR878", "HDBR881", "HDBR882", "HDBR883", "HDBR884", "HDBR885", "HDBR886", "HDBR887",
"HDBR889", "HDBR890", "HDBR891", "HDBR892", "HDBR894", "HDBR896", "HDBR897", "HDBR898", "HDBR899", "HDBR965", "HDBR966",
"HDBR967", "HDBR968", "HDBR969"
);

#my @require_samples = (
#"HDBR375", "HDBR380", "HDBR383", "HDBR705", "HDBR706", "HDBR707", "HDBR708", "HDBR709", "HDBR710", "HDBR712", "HDBR713", "HDBR714", 
#"HDBR715", "HDBR716", "HDBR717", "HDBR718", "HDBR720", "HDBR721", "HDBR722", "HDBR723", "HDBR724", "HDBR725", "HDBR726", "HDBR727", 
#"HDBR728", "HDBR729", "HDBR731", "HDBR732", "HDBR733", "HDBR734", "HDBR735", "HDBR736", "HDBR737", "HDBR738", "HDBR739", "HDBR740", 
#"HDBR741", "HDBR742", "HDBR743", "HDBR744", "HDBR745", "HDBR746", "HDBR747", "HDBR748", "HDBR749", "HDBR750", "HDBR751", "HDBR752"
#);

my $count_file = "/sharedlustre/IGM/ymtgroup-work/Temp_HDBR_RNAseq_2014Jul/Analysis_2014Jul/DESeq2/counts_filtered_15MreadsSample.txt";
my $rpkm_file= "/sharedlustre/IGM/ymtgroup-work/Temp_HDBR_RNAseq_2014Jul/Analysis_2014Jul/DESeq2/RPKM.cqn.noScaling_15MreadsSamples.txt";
#my %candidate_V;
#open GENO, $genotyping_files[0] or die "Can not reads the file $genotyping_files[0]\n";
#read in genotyped variants list
#close GENO;

my %candidate_V = ( '3_36993455_G_A' => 'rs1800734_NA' , '3_37012077_A_C,G' => 'rs1799977_MLH1' ); # need to be GRCh38 reference
#my %candidate_V = ( '3_37012077_A_C,G' => 'rs1799977_MLH1' );
my %candidate_V_indexed;  # multiple layer hash: chr pos  = ref_alt_rs_geneName

#######
#Initiation
#######
my %require_samples_hashed;
foreach my $sample (@require_samples) {
	$require_samples_hashed{$sample} = 1;
}

my %vcf_search_results; #Only interested in SNV. chr pos sample => chr\tpos\tref\talt\trefCount\taltCount(s)
	                                            #chr pos 'VCFinputlines' => array of closeby vcf lines
my %genochip_search_results;
my %genecount_search_results;
my %geneRPKM_search_results;
my %seachingGenes;

foreach my $key (keys %candidate_V) {
	my @ele=split("_",$key);
	
	if ($ele[0] =~ /^x$/) {$ele[0] = "X";}
	if ($ele[0] =~ /^y$/) {$ele[0] = "Y";}
	if ($ele[0] =~ /^MT$/) {$ele[0] = "M";}
	if ($ele[0] =~ /^m$/) {$ele[0] = "M";}
	$ele[0] =~ s/\.\d+$//;
	my $chr = "chr".$ele[0];
	if (! exists $candidate_V_indexed{$chr}{$ele[1]}) {
		$candidate_V_indexed{$chr}{$ele[1]} = $ele[2]."_".$ele[3]."_".$candidate_V{$key};
		my $array = [];
		$vcf_search_results{$chr}{$ele[1]}{"VCFinputlines"} = $array;
		foreach my $sample (@require_samples) {
			$vcf_search_results{$chr}{$ele[1]}{$sample} = "NA\tNA\tNA\tNA\tNA\tNA";
			$genochip_search_results{$chr}{$ele[1]}{$sample} = "NA\tNA\tNA";
			my @temp_split = split("_", $candidate_V{$key});
			if (! exists $seachingGenes{$temp_split[1]}) {
				$seachingGenes{$temp_split[1]} = 1;
			}
			if (! exists $genecount_search_results{$sample}{$temp_split[1]}) {
				$genecount_search_results{$sample}{$temp_split[1]} = "NA\tNA";
				$geneRPKM_search_results{$sample}{$temp_split[1]} = "NA\tNA";
			}
		}
	} else {
		print "Duplicated query SNP positions: $chr $ele[1].\n";
	}
}
#######
#Initiation done
#######

my %vcf_samples; # sample names => sample index in vcf, 0 based
FILE: foreach my $vcf (@vcf_files) {
	open VCF, $vcf or die "Can not open the file $vcf\n";
	print "Reading VCF file $vcf\n";
	LINE: while (my $line = <VCF>) {
		if ($line =~ m/^\#CHR/ && ! %vcf_samples) {
			chomp $line;
			my @ele = split("\t", $line);
			for (my $i=9; $i < scalar @ele; $i++) {
				if (! exists $vcf_samples{$ele[$i]} ) {
					if (exists $require_samples_hashed{$ele[$i]}) {
						$vcf_samples{$ele[$i]} = $i;
						print "Found sample $ele[$i] in the vcf on colum $i (0-based count).\n";
					}
				} else {
					print "Error: VCF have duplicated sample names :$ele[$i]\n";
					exit;
				}
			}
			foreach my $sample (keys %require_samples_hashed) {
				if ( ! exists $vcf_samples{$sample}) {
					print "Sample $sample is not found in the vcfs;\n";
					exit;
				}
			}
		} elsif($line =~ m/^\#/) {
			next LINE;
		} else {  # for the sake of the speed amd memory, only process 20bp up and down stream of any targeted variants
 			chomp $line;
			my @ele = split("\t", $line);
			my $chr = $ele[0];
			my $pos = $ele[1];
			if (exists $candidate_V_indexed{$chr}) {
				foreach my $q_pos ( sort {$a <=> $b} keys %{$candidate_V_indexed{$chr}} ) {
					if( $pos >= $q_pos-20 && $pos <= $q_pos) {
						print "push in a line\n$line\n";
						push @{$vcf_search_results{$chr}{$q_pos}{"VCFinputlines"}}, $line;
					}
				}
				
				#start to process VCF lines
				if (eof) {
					print "Reached to the end of the file:\n";
					foreach my $q_pos (sort {$a <=> $b} keys %{$vcf_search_results{$chr}}) {
						my @lines = @{$vcf_search_results{$chr}{$q_pos}{"VCFinputlines"}};
						if (scalar @lines >= 0) {
							QLINE: foreach my $q_line (@lines) {
								print "processing line:$q_line\n";
								my @q_ele = split("\t", $q_line);
								############################## vcf filtering creterion should be placed here ####
								if ($q_ele[5] >= 20) { #### filter 1 #### low qual filter
									my $l_chr = $q_ele[0];
									my $l_pos = $q_ele[1];
									my $ref = $q_ele[3];
									my @alts = split(",", $q_ele[4]);
									my @l_format_fields = split( ":", $q_ele[8]);
									my $stretch_length = length $ref;
									foreach my $alt_temp (@alts) {
										if (length $alt_temp > $stretch_length) {
											$stretch_length = length $alt_temp;
										}
									}
									if ($q_pos > ($l_pos + $stretch_length -1) ) {
										next QLINE;
									}
									my ($GT_index, $GQ_index, $RO_index, $AO_index, $DP_index);
									for (my $i=0;$i<scalar(@l_format_fields);$i++){
										if ($l_format_fields[$i] =~ m/GQ/) {
											$GQ_index = $i+1;
											print "GQ_index is $GQ_index \n";
										}
										if ($l_format_fields[$i] =~ m/GT/) {
											$GT_index = $i+1;
											print "GT_index is $GT_index \n";
										}
										if ($l_format_fields[$i] =~ m/DP/) {
											$DP_index = $i+1;
											print "DP_index is $DP_index \n";
										}
										if ($l_format_fields[$i] =~ m/RO/) {
											$RO_index = $i+1;
											print "RO_index is $RO_index \n";
										}
										if ($l_format_fields[$i] =~ m/AO/) {
											$AO_index = $i+1;
											print "AO_index is $AO_index \n";
										}
									}
									
									if ($GT_index && $GQ_index && $RO_index && $AO_index && $DP_index) {
									SAMPLE:	foreach my $sample (keys %vcf_samples) {
											if ($sample =~ m/^\./) {
												next SAMPLE;
											}
											my @sample_cols = split(":", $q_ele[$vcf_samples{$sample}]);
											if (scalar @sample_cols == scalar @l_format_fields) { #### filter 2 #### no call filter
												my @bases = split("/", $sample_cols[$GT_index-1]);
												my ($alt, $alt_cov);
												if ($bases[1] == 0) {
													$alt = $alts[0]; $alt_cov = 0;
												} else {
													$alt = $alts[$bases[1] - 1];
													my @aos = split(",", $sample_cols[$AO_index-1]);
													$alt_cov = $aos[$bases[1] - 1];
												}
												#######################################
												# work out true location and ref alt
												########
												#compare ref and alt length
												my @ref_bases = split("", $ref);
												my @alt_bases = split("", $alt);
												print "ref length: ".scalar @ref_bases."\n"; 
												print "var length: ".scalar @alt_bases."\n"; 
												if (scalar @ref_bases == scalar @alt_bases) {#same length
													my $this_pos = $l_pos;
													my $this_ref = $ref;
													my $this_alt = $alt;
													for (my $i=0;$i <= scalar @ref_bases; $i++){ #remove same bases
														if ($ref_bases[$i] eq $alt_bases[$i]) {
															$this_pos += 1;
															$this_ref =~ s/^.//;
															$this_alt =~ s/^.//
														} else {
															last;
														}
													}
													my $type = "SNV";
													if (scalar @ref_bases > ($this_pos - $l_pos + 1)) {
														for(my $i = ($this_pos - $l_pos + 1);$i < scalar @ref_bases;$i++) {
															if ($ref_bases[$i] ne $alt_bases[$i]) {
																$type="REP";
																last;
															}
														}
													}
													
													if ($type eq "SNV" && $this_pos == $q_pos) {
														print "SNV found\n";
														if($vcf_search_results{$chr}{$q_pos}{$sample} eq "NA\tNA\tNA\tNA\tNA\tNA") {
															$vcf_search_results{$chr}{$q_pos}{$sample} =
														 $l_chr."\t".$this_pos."\t".$this_ref."\t".$this_alt."\t".$sample_cols[$RO_index-1]."\t".$alt_cov;
														}
													} elsif ($type ne "SNV") {
														print "Not SNV\n";
														if ($q_pos >= $this_pos && $q_pos <= ($this_pos + length $this_ref -1) ) {
															if($vcf_search_results{$chr}{$q_pos}{$sample} eq "NA\tNA\tNA\tNA\tNA\tNA") {
																$this_ref = substr $this_ref, $q_pos-$this_pos, 1;
																$this_alt = substr $this_alt, $q_pos-$this_pos, 1;
																$vcf_search_results{$chr}{$q_pos}{$sample} =
																 $l_chr."\t".$q_pos."\t".$this_ref."\t".$this_alt."\t".$sample_cols[$RO_index-1]."\t".$alt_cov;
															}
														} else {
															next SAMPLE;
														}
													}
												} elsif (scalar @ref_bases > scalar @alt_bases) { #different lengths - deletion
													my $this_pos = $l_pos;
													my $this_ref = $ref;
													my $this_alt = $alt;
													for (my $i=0;$i <= scalar @alt_bases; $i++){ #remove same bases
														if ($ref_bases[$i] eq $alt_bases[$i]) {
															$this_pos += 1;
															$this_ref =~ s/^.//;
															$this_alt =~ s/^.//
														} else {
															last;
														}
													}
													if ($q_pos >= $this_pos && $q_pos <= ($this_pos + length $this_ref -1) ) {
														if($vcf_search_results{$chr}{$q_pos}{$sample} eq "NA\tNA\tNA\tNA\tNA\tNA" ) {
															$this_ref = substr $this_ref, $q_pos-$this_pos, 1; 
															if (length $this_alt >=  $q_pos-$this_pos + 1) {
																$this_alt = substr $this_alt, $q_pos-$this_pos, 1;
															} else {
																$this_alt = "-";
															}
															$vcf_search_results{$chr}{$q_pos}{$sample} =
															$l_chr."\t".$q_pos."\t".$this_ref."\t".$this_alt."\t".$sample_cols[$RO_index-1]."\t".$alt_cov;
														}
													} else {
														next SAMPLE;
													}
												} else { #different lengths - insertion
													my $this_pos = $l_pos;
													my $this_ref = $ref;
													my $this_alt = $alt;
													for (my $i=0;$i <= scalar @ref_bases; $i++){ #remove same bases
														if ($ref_bases[$i] eq $alt_bases[$i]) {
															$this_pos += 1;
															$this_ref =~ s/^.//;
															$this_alt =~ s/^.//
														} else {
															last;
														}
													}
													if ($q_pos >= $this_pos && $q_pos <= ($this_pos + length $this_alt -1) ) {
														if($vcf_search_results{$chr}{$q_pos}{$sample} eq "NA\tNA\tNA\tNA\tNA\tNA" ) {
															$this_alt = substr $this_alt, $q_pos-$this_pos, 1; 
															if (length $this_ref >=  $q_pos-$this_pos + 1) {
																$this_ref = substr $this_ref, $q_pos-$this_pos, 1;
															} else {
																$this_ref = "-";
															}
															$vcf_search_results{$chr}{$q_pos}{$sample} =
															$l_chr."\t".$q_pos."\t".$this_ref."\t".$this_alt."\t".$sample_cols[$RO_index-1]."\t".$alt_cov;
														}
													} else {
														next SAMPLE;
													}
												}
											}
										}
									}
								}
							}
						}
						delete $vcf_search_results{$chr}{$q_pos}{"VCFinputlines"};
					}
				}
			} else {
				close VCF;
				next FILE;
			}
		}
	}
	close VCF; 
}

# reading genochip files
print "Start to read genochip files..\n";
GENOFILE: foreach my $geno_file (@genotyping_files) {
	my %sample_index;
	print "Reading genofile $geno_file....\n";
	open GENO, "$geno_file" or die "Cannot open the genotyping file $geno_file\n";
	GENOLINE: while (my $line = <GENO>) {
		chomp $line;
		if ($line =~ m/^Index\_chip/) {
			my @ele = split("\t", $line);
			for (my $i=0; $i < scalar @ele; $i++){
				if ( $ele[$i] =~ m/A2311[\_|\-](\d+.*)\.GType/ ) {
					my $sample_name = $1;
					$sample_name =~ s/^0//;
					$sample_name =~ s/\_\d$//;
					$sample_name = "HDBR".$sample_name;
					if (exists $require_samples_hashed{$sample_name}) {
						$sample_index{$sample_name} = $i;
						print "Found sample $sample_name\n";
					}
				}
			}
			unless (%sample_index) {next GENOFILE;}
		} else {
			my @ele = split("\t", $line);
			if ($ele[8] =~ /^x$/) {$ele[8] = "X";}
			if ($ele[8] =~ /^y$/) {$ele[8] = "Y";}
			if ($ele[8] =~ /^MT$/) {$ele[8] = "M";}
			if ($ele[8] =~ /^m$/) {$ele[8] = "M";}
			$ele[0] =~ s/\.\d+$//;
			my $chr = "chr".$ele[8];
			my $pos = $ele[9];
			if (exists $candidate_V_indexed{$chr}) {
				if (exists $candidate_V_indexed{$chr}{$pos} && $ele[5] !~ m/N/ && $ele[6] !~ m/N/) {
					my $geno_A = $ele[5];
					my $geno_B = $ele[6];
					my $vcf_ref = $ele[11];
					my @vcf_alts = split(",",$ele[12]);
					##Determine A B in respect the strand of the ref and alt
					my ($temp_A, $temp_B);
					if ($geno_A eq $vcf_ref) {
						foreach my $vcf_alt (@vcf_alts) {
							if ($geno_B eq $vcf_alt) {
								print "Genochip AB match to the same strand as Ref and Alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								$temp_A = $vcf_ref;
								$temp_B = $vcf_alt;
								last;
							}
						}
						if ( (! defined $temp_A) || ( ! defined $temp_B) ) {
							if (scalar @vcf_alts == 2) {
								if ( ($geno_A eq $vcf_alts[0] && $geno_B eq $vcf_alts[1])	
								 || ($geno_A eq $vcf_alts[1] && $geno_B eq $vcf_alts[0]) ) {
								 	print "Genochip AB match to the same strand of muliti alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	$temp_A = $geno_A;
									$temp_B = $geno_B;
								 } elsif ( ($geno_A eq Complementary($vcf_alts[0]) && $geno_B eq Complementary($vcf_alts[1]) )	
								 || ($geno_A eq Complementary($vcf_alts[1]) && $geno_B eq Complementary($vcf_alts[0]) ) ) {
								 	print "Genochip AB match to the reverse strand of muliti alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	$temp_A = Complementary($geno_A);
									$temp_B = Complementary($geno_B);
								 } else {
								 	print "No match AB to ref alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	next GENOLINE;
								 }
							} else {
								print "No match AB to ref alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								next GENOLINE;
							}
						}
					} elsif (IfComplementary($geno_A, $vcf_ref)) {
						foreach my $vcf_alt (@vcf_alts) {
							if (IfComplementary($geno_B, $vcf_alt)) {
								$temp_A = $vcf_ref;
								$temp_B = $vcf_alt;
								print "Ginochip A is the reverse strand ref: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								last;
							}
						}
						if ( (! defined $temp_A) || ( ! defined $temp_B) ) {
							if (scalar @vcf_alts == 2) {
								if ( ($geno_A eq $vcf_alts[0] && $geno_B eq $vcf_alts[1])	
								 || ($geno_A eq $vcf_alts[1] && $geno_B eq $vcf_alts[0]) ) {
								 	print "Genochip AB match to the same strand of muliti alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	$temp_A = $geno_A;
									$temp_B = $geno_B;
								 } elsif ( ($geno_A eq Complementary($vcf_alts[0]) && $geno_B eq Complementary($vcf_alts[1]) )	
								 || ($geno_A eq Complementary($vcf_alts[1]) && $geno_B eq Complementary($vcf_alts[0]) ) ) {
								 	print "Genochip AB match to the reverse strand of muliti alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	$temp_A = Complementary($geno_A);
									$temp_B = Complementary($geno_B);
								 } else {
								 	print "No match AB to ref alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	next GENOLINE;
								 }
							} else {
								print "No match AB to ref alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								next GENOLINE;
							}
						}
					} elsif ($geno_B eq $vcf_ref) {
						foreach my $vcf_alt (@vcf_alts) {
							if ($geno_A eq $vcf_alt) {
								$temp_A = $vcf_alt;
								$temp_B = $vcf_ref;
								print "Ginochip B is the same strand ref: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								last;
							}
						}
						if ( (! defined $temp_A) || ( ! defined $temp_B) ) {
							if (scalar @vcf_alts == 2) {
								if ( ($geno_A eq $vcf_alts[0] && $geno_B eq $vcf_alts[1])	
								 || ($geno_A eq $vcf_alts[1] && $geno_B eq $vcf_alts[0]) ) {
								 	print "Genochip AB match to the same strand of muliti alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	$temp_A = $geno_A;
									$temp_B = $geno_B;
								 } elsif ( ($geno_A eq Complementary($vcf_alts[0]) && $geno_B eq Complementary($vcf_alts[1]) )	
								 || ($geno_A eq Complementary($vcf_alts[1]) && $geno_B eq Complementary($vcf_alts[0]) ) ) {
								 	print "Genochip AB match to the reverse strand of muliti alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	$temp_A = Complementary($geno_A);
									$temp_B = Complementary($geno_B);
								 } else {
								 	print "No match AB to ref alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	next GENOLINE;
								 }
							} else {
								print "No match AB to ref alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								next GENOLINE;
							}
						}
					} elsif (IfComplementary($geno_B, $vcf_ref) ) {
						foreach my $vcf_alt (@vcf_alts) {
							if (IfComplementary($geno_A, $vcf_alt)) {
								$temp_A = $vcf_alt;
								$temp_B = $vcf_ref;
								print "Ginochip B is the reverse strand ref: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								last;
							}
						}
						if ( (! defined $temp_A) || ( ! defined $temp_B) ) {
							if (scalar @vcf_alts == 2) {
								if ( ($geno_A eq $vcf_alts[0] && $geno_B eq $vcf_alts[1])	
								 || ($geno_A eq $vcf_alts[1] && $geno_B eq $vcf_alts[0]) ) {
								 	print "Genochip AB match to the same strand of muliti alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	$temp_A = $geno_A;
									$temp_B = $geno_B;
								 } elsif ( ($geno_A eq Complementary($vcf_alts[0]) && $geno_B eq Complementary($vcf_alts[1]) )	
								 || ($geno_A eq Complementary($vcf_alts[1]) && $geno_B eq Complementary($vcf_alts[0]) ) ) {
								 	print "Genochip AB match to the reverse strand of muliti alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	$temp_A = Complementary($geno_A);
									$temp_B = Complementary($geno_B);
								 } else {
								 	print "No match AB to ref alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								 	next GENOLINE;
								 }
							} else {
								print "No match AB to ref alt: A: $geno_A; B: $geno_B; Ref: $vcf_ref;  Alt: @vcf_alts;\n";
								next GENOLINE;
							}
						}
					} else {
						next GENOLINE;
					}
					###################################################################
					GENOSAMPLE: foreach my $sample (keys %sample_index) {
						print "process sample: $sample ";
						my $geno_call = $ele[$sample_index{$sample}];
						my $geno_call_score = $ele[$sample_index{$sample}+1];
						if ($geno_call !~ m/D/ && $geno_call !~ m/I/ && $geno_call !~ m/NC/) {
							$geno_call =~ s/A/$temp_A/g;
							$geno_call =~ s/B/$temp_B/g;
							my @geno_call_split = split("",$geno_call);
							$genochip_search_results{$chr}{$pos}{$sample} = 
								$geno_call_split[0]."\t".$geno_call_split[1]."\t".$geno_call_score;
							print "call: $geno_call\n";
						} else {
							$genochip_search_results{$chr}{$pos}{$sample} = 
								"NC/DI"."\t"."NC/DI"."\t"."NA";
								print "call: NC/DI\n";
						}
					}
				} else {
					next GENOLINE;
				}
			} else {
				next GENOLINE;
			}
		}
	}
	close GENO;
}

print "Start to read the count file ..\n";
open COUNT, "$count_file" or die "Cannot open the genotyping file $count_file\n";
my %count_samples_index;
while (my $line = <COUNT> ) {
	chomp $line;
	if ($line =~ m/^GeneName/) {
		my @ele=split("\t", $line);
		for(my $i=1; $i < scalar @ele; $i++) {
			my $sample;
			if ($ele[$i] =~ m/\.(\w+)$/) {
				$sample = $1;
			}
			if (exists $require_samples_hashed{$sample}){
				$count_samples_index{$sample} = $i;
			}
		}
	} else {
		my @ele=split("\t", $line);
		if (exists $seachingGenes{$ele[0]}) {
			foreach my $sample (keys %require_samples_hashed) {
				$genecount_search_results{$sample}{$ele[0]} = $ele[0]."\t".$ele[$count_samples_index{$sample}];
			}
		}
	}
}
close COUNT;

print "Start to read the RPKM file ..\n";
open RPKM, "$rpkm_file" or die "Cannot open the genotyping file $rpkm_file\n";
my %rpkm_samples_index;
while (my $line = <RPKM> ) {
	chomp $line;
	if ($line =~ m/^GeneName/) {
		my @ele=split("\t", $line);
		for(my $i=1; $i < scalar @ele; $i++) {
			my $sample;
			if ($ele[$i] =~ m/\.(\w+)$/) {
				$sample = $1;
			}
			if (exists $require_samples_hashed{$sample}){
				$rpkm_samples_index{$sample} = $i;
			}
		}
	} else {
		my @ele=split("\t", $line);
		if (exists $seachingGenes{$ele[0]}) {
			foreach my $sample (keys %require_samples_hashed) {
				$geneRPKM_search_results{$sample}{$ele[0]} = $ele[0]."\t".$ele[$count_samples_index{$sample}];
			}
		}
	}
}
close RPKM;



###start to output file
print "Output result to file....\n";
open OUTPUT, ">$output_file" or die "Cannot open the file to output $output_file\n";
#print the header first
print OUTPUT "Sample\tChromosome\tPos(GRCh38)\tRef(dbSNP141)\tAlt(dbSNP141)\tdbSNP_ID\tGeneName\tChromosome".
"\tPos\tRef(RNAseq)\tAlt(RNAseq)\tRefCount(RNAseq)\tAltCount(RNAseq)".
"\tAllele1(Genochip)\tAllele2(Genochip)\tScore(Genochip)\tGeneName\tCount(RNAseq)\tGeneName\tRPKM(RNAseq)\n";

#my %candidate_V_indexed;  # multiple layer hash: chr pos  = ref_alt_rs_geneName
foreach my $chr (keys %candidate_V_indexed) {
	foreach my $pos (sort {$a <=> $b} keys %{$candidate_V_indexed{$chr}}) {
		foreach my $sample (sort {$a cmp $b} keys %require_samples_hashed) {
			my @ele=split("_", $candidate_V_indexed{$chr}{$pos});
			my $temp_string = "$sample\t$chr\t$pos\t$ele[0]\t$ele[1]\t$ele[2]\t$ele[3]\t".$vcf_search_results{$chr}{$pos}{$sample};
			$temp_string = $temp_string."\t".$genochip_search_results{$chr}{$pos}{$sample};
			$temp_string = $temp_string."\t".$genecount_search_results{$sample}{$ele[3]};
			$temp_string = $temp_string."\t".$geneRPKM_search_results{$sample}{$ele[3]};
			print OUTPUT $temp_string."\n";
		}
	}
}
close OUTPUT;

print "Done!\n";
exit;

sub IfComplementary {
	my ($a, $b) = @_;
	$a =~ tr/ATGC/TACG/;
	if ($a eq $b) {
		return 1;
	} else {
		return 0;
	}
}
sub Complementary {
	my ($a) = @_;
	$a =~ tr/ATGC/TACG/;
	return $a;
}


