#!/usr/bin/perl;
use strict;
use warnings;
use Getopt::Long;
use Excel::Writer::XLSX;

my ($sampleNames, $sampleColumns, $input_file, $output_excel, $help);
my $CNV_file="";

usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, 'sampleNames=s' => \$sampleNames, 'sampleColumns=s' => \$sampleColumns, 'in=s' => \$input_file, 'out=s' => \$output_excel, 'CNV=s' => \$CNV_file ) || defined $help );

#Change the file place if needed, its just for link CNVs to Gene Names
#my $genefile="/users/a5907529/lustre/scripts/AfterDindel/Ensembl_Genes_R67.txt";
my $genefile="Ensembl_Genes_R67.txt";
my $en_genes_ref;
my %en_genes;
my $insterested_gene_file; #This file should have gene vs ensembl id just for CNV data exists
if ($CNV_file ne "") {
	$en_genes_ref = GetGeneCoords($genefile);
	%en_genes = %$en_genes_ref;
	$insterested_gene_file = "Genes_PID_05NOV2013_YAOBO.txt"; ######## may need to change if the file is in diff place
}

unless (defined $sampleNames) {
	die "You have not supplied sample names using --sampleNames\n";
}

unless (defined $sampleColumns) {
	die "You have to specify which column(s) have sample genotype using --sampleColumns\n";
}

unless (defined $input_file) {
	die "You have not supplied a input file with --in\n";
}

unless (defined $output_excel) {
	die "You have not supplied the output file name with --out\n";
}

my @sample_names = split(",", $sampleNames);
print "Looking for sample:\n";
foreach my $sample (@sample_names) {
	print $sample."\n";
}

my @sample_columns = split(",", $sampleColumns);
print "On columns (0-based):\n";
foreach my $column (@sample_columns) {
	print $column."\n";
}

if(scalar @sample_columns != scalar @sample_names) {
	die "Number of sample names and columns are inconsistent.\n";
}

my @output_all;
my @output_filtered;
my @output_possibleHits;
my @very_rare_interested_vars;  #to filter very rare < 0.01 variants on interested genes
my @output_Xlinked;
my @output_ARmodelHits;
my @output_CNV;
my @output_exp_compound;
my %CNV_regions;
my %gene_name_count; #it's for selecting Xlinked rare variants



# To find extended compound hetero: multiple variant hits and CNVs on the same gene 
my %all_gene_id; # to record gene_id => gene_name
my %CNV_on_samples; # to record all CNVs in a new manner [Gene][sample]=cnv_id
my %Variants_on_samples; # to record all rare variants in a new manner [Gene][sample]=chr_v-position_genoType_alleleFreq(ESP)_alleleFreq(1000G)_alleleFreq(in house);

if ($CNV_file ne "") {
	print "CNV file provided:$CNV_file\n";
	open CNV, $CNV_file or die "Cannot open the CNV file $CNV_file\n";
	my  %CNV_ids;
	while (my $line=<CNV>) {
		chomp $line;
	 	if ($line eq "") {
			next;
		}
		# remove \t within double quotes
		my @temp=split("\t", $line);
		my @elements;
		my $temp_string="";
		my $extension_switch = "OFF";
		foreach my $ele (@temp) {
			chomp $ele;
			if ($extension_switch eq "OFF") {
				if($ele =~ m/^".*"$/) {
					push @elements, $ele;
				} elsif ($ele =~ m/^"/) {
					$extension_switch = "ON";
					$temp_string=$ele;
				} else {
					push @elements, $ele;
				}
			} else {
				if ($ele =~ m/"$/) {
					$extension_switch = "OFF";
					push @elements, $temp_string." ".$ele;
				 } else {
					$temp_string=$temp_string." ".$ele;
				}
			}
		}
		if (scalar @elements < 21) {
			my $length=scalar @elements;
			for (my $hh = 0;$hh<(21-$length);$hh++) {
				push @elements, "NA";
			}
		}
		
		if($elements[0] =~ m/SampleID/) {
			next;
		}
		#print scalar @elements."\n";
		my $id = $elements[0];
		if (ElementInArray(\@sample_names, $id) eq "YES") { # if the ID is in sample names
			my $temp_temp_string="";
			foreach my $cell (@elements) {
				#print $cell."\n";
				$temp_temp_string = $temp_temp_string."\t".$cell;
			}
			$temp_temp_string =~ s/^\t//;
			push @output_CNV, $temp_temp_string;
			my $start=$elements[5];
			my $end=$elements[6];
			$elements[7] =~ s/"//g;
			my $chr="chr".$elements[7];
			$CNV_regions{$id}{$chr}{$start}{$end} = $elements[8] ;
			unless( exists $CNV_ids{$id} ) {
				$CNV_ids{$id} =0;
				print "The file contains info for $id.\n";
			}
		}
	}
}
close CNV;

open INPUT, $input_file or die "Cannot open the file $input_file.\n";
while (my $line=<INPUT>) {
	if ($line !~ m/^Func\t/) {
		my @elements = split("\t", $line);
		if ($elements[26] ne "unknown") { # insert "unknown" to make every row in the same length
			splice @elements, 26, 0, 'unknown'; 
		}
		
		if($elements[2] eq "synonymous SNV") { # eliminate synonymous SNV
			next;
		}
		
		my $sample_string="";
		my $rr_na_sum=0;
		for (my $i=0;$i<scalar @sample_names;$i++) {
			my $on_CNV="NO";
			start_loop: foreach my $start ( sort {$a<=>$b} keys %{$CNV_regions{$sample_names[$i]}{$elements[27]}} ) {
				foreach my $end (keys %{$CNV_regions{$sample_names[$i]}{$elements[27]}{$start}}) {
					if ( $elements[28] >= $start && $elements[28] <= $end) {
						$on_CNV=$CNV_regions{$sample_names[$i]}{$elements[27]}{$start}{$end};
						last start_loop;
					}
				}
			}
			$sample_string=$sample_string.$elements[$sample_columns[$i]]."\t".$elements[$sample_columns[$i]+1]."\t".$on_CNV."\t";
			if($elements[$sample_columns[$i]+1] eq 'R/R' || $elements[$sample_columns[$i]+1] eq 'NA') {
				$rr_na_sum+=1;
			}
		}
		
		my $length = scalar @elements -1;
		
		#Change 1000genome and ESP6500 var freq empty value to "0"
		if($elements[6] eq "" ) {
			$elements[6] = 0;
		}
		if($elements[7] eq "" ) {
			$elements[7] = 0;
		}
		#Change in-house Freq 'NA' to '0'
		if ($elements[$length-6] eq "NA") {
			$elements[$length-6] = 0;
		}
		my $string_1 = $elements[27]."\t".$elements[28]."\t".$elements[29]."\t".$elements[30]."\t"
		.$elements[31]."\t".$elements[35]."\t".$elements[32]."\t".$elements[0]."\t".$elements[$length-9]."\t".$elements[$length-8]
		."\t".$elements[$length-7]."\t".$elements[2]."\t".$elements[3]."\t".$elements[4]."\t".$elements[5]."\t".$elements[6]."\t"
		.$elements[7]."\t".$elements[$length-6]."\t".$elements[8]."\t".$elements[9]."\t".$elements[10]."\t".$elements[11]."\t".$elements[12]."\t".
		$elements[13]."\t".$elements[14]."\t".$elements[15]."\t".$elements[16]."\t".$elements[17]."\t".$elements[18]."\t".$elements[19]."\t".$elements[20]."\t"
		.$elements[33]."\t".$elements[36]."\t"; #re-arrange columns
		
		
		$elements[$length-2] =~ s/\'/\"/g;
		$elements[$length-1] =~ s/\'/\"/g;
		$elements[$length] =~ s/\'/\"/g;
		
		my $string_2 = $elements[34]."\t".$elements[$length-5]."\t".$elements[$length-4]."\t".$elements[$length-3]."\t".$elements[$length-2]
		."\t".$elements[$length-1]."\t".$elements[$length]; #re-arrange columns
		
		if($rr_na_sum < scalar @sample_names) { # RR and NA are removed
			push @output_all, $string_1.$sample_string.$string_2;
			
			# for filtered spreadsheet

			if( $elements[6] <= 0.05 && $elements[7] <= 0.05 ) {
				push @output_filtered, $string_1.$sample_string.$string_2;
				
				unless(exists $gene_name_count{$elements[$length-9]}) { # could be deleted
					$gene_name_count{$elements[$length-9]} = 0;
				}
				
				# for deleterious hits spreadsheet
				my $number_deleterious_predictions=0;
				if($elements[13] eq "D") {
					$number_deleterious_predictions+=1;
				}
				if($elements[15] eq "D" || $elements[15] eq "P") {
					$number_deleterious_predictions+=1;
				}
				if($elements[17] eq "D") {
					$number_deleterious_predictions+=1;
				}
				if($elements[19] eq "A" || $elements[19] eq "D") {
					$number_deleterious_predictions+=1;
				}
				
				my $number_CompoundHet=0;
				foreach my $sample_col (@sample_columns) {
					if($elements[$sample_col+1] eq "V1/V2") {
						$number_CompoundHet +=1;
					}
				}
				
				for( my $i=0; $i< scalar @sample_columns;$i++) {
					my $sample_col = $sample_columns[$i];
					if($elements[$sample_col+1] eq "V1/V2") {
						$number_CompoundHet +=1;
					}
					
					# extend the hash to record variants for samples and genes for %Variants_on_samples: $variants_on_samples{sample_name}{gene_id}{gene_name}=chr_pos_genoType
					unless (exists $all_gene_id{$elements[$length-8]}) { $all_gene_id{$elements[$length-8]} = $elements[$length-9]."_".$elements[$length-7]; }
					if ( exists $Variants_on_samples{$sample_names[$i]}{$elements[$length-8]} ) {
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-8]} = 
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-8]}.";".
						$elements[27]."_".$elements[28]."_".$elements[$sample_col+1]."_".sprintf( "%.4f",$elements[6])."_".sprintf( "%.4f",$elements[7])."_".sprintf( "%.4f",$elements[$length-6]);
					} else {
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-8]} = 
						$elements[27]."_".$elements[28]."_".$elements[$sample_col+1]."_".sprintf( "%.4f",$elements[6])."_".sprintf( "%.4f",$elements[7])."_".sprintf( "%.4f",$elements[$length-6]);				
					}
				}
				
				#for rare variants on a same gene
				if($number_CompoundHet < scalar @sample_names) {
					$gene_name_count{$elements[$length-9]} =$gene_name_count{$elements[$length-9]}+1;
				} #retired
				
				if ($elements[$length-7] eq "YES" && $number_deleterious_predictions >=2) {  #those that were predicted as “Deleterious” by at least two predictors
					push @output_possibleHits, $string_1.$sample_string.$string_2;
				} elsif ($elements[$length-7] eq "NO" && $number_deleterious_predictions ==4) { #those that were predicted as “Deleterious” by all 4 predictors
					push @output_possibleHits, $string_1.$sample_string.$string_2;
				} elsif($elements[35] eq "SNP" && $number_CompoundHet >=1) { #Compound heterozygous “V1/V2”* of single base polymorphisms
					push @output_possibleHits, $string_1.$sample_string.$string_2;
				} elsif($elements[$length-7] eq "YES" && $elements[6] <= 0.01 && $elements[7] <= 0.01) { # those that have AAF less than 0.01 in both 1000Genomes project and ESP6500 
					push @output_possibleHits, $string_1.$sample_string.$string_2;
				}
				
				if ($elements[$length-7] eq "YES" && $elements[6] <= 0.01 && $elements[7] <= 0.01) {
					push @very_rare_interested_vars, $string_1.$sample_string.$string_2;
				}
				
				#for Xlinked rare variants
				if ($elements[27] =~ m/chrx/i) {
					push @output_Xlinked, $string_1.$sample_string.$string_2;
				}
				
				# AR model Hits
				my $number_Homo=0;
				foreach my $sample_col (@sample_columns) {
					if($elements[$sample_col+1] eq "V/V") {
						$number_Homo +=1;
					}
				}
				if ($number_Homo == scalar @sample_names) {
					if ( $elements[27] !~ m/chr[XYMxym]/ ) {
						push @output_ARmodelHits, $string_1.$sample_string.$string_2;
					}
				}
			}
		}
	}
}
close(INPUT);


###################################### working here #################################################################
# record CNVs into the hash %CNV_on_samples
if ($CNV_file ne "") {
	####
	# need to read in the whole list of insterested genes, preferred in ENSEMBL id
	my %insterested_genes;
	open IN_GENES, $insterested_gene_file or die "Cannot open the file $insterested_gene_file.\n";
	while (my $line=<IN_GENES>) {
		my @temp_haha=split("\t", $line);
		if (! exists $insterested_genes{$temp_haha[3]}) {
			$insterested_genes{$temp_haha[4]} = 1;
		}
	}
	close (IN_GENES);
	####
	foreach my $sample (keys %CNV_regions) {
		foreach my $chr (keys %{$CNV_regions{$sample}}) {
			foreach my $start_CNV (keys %{$CNV_regions{$sample}{$chr}}) {
				foreach my $end_CNV (keys %{$CNV_regions{$sample}{$chr}{$start_CNV}}) {
					start_loop: foreach my $start_gene ( sort {$a<=>$b} keys %{$en_genes{$chr}} ) {
						foreach my $end_gene (keys %{$en_genes{$chr}{$start_gene}}) {
							if ( $start_gene <= $end_CNV && $end_gene >= $start_CNV ) {
								my @temp = split("\t", $en_genes{$chr}{$start_gene}{$end_gene});
								my $gene_name = $temp[0];
								my $gene_id = $temp[1];
								unless (exists $all_gene_id{$gene_id}) {
									if ( ! exists $insterested_genes{$gene_id}) {
										$all_gene_id{$gene_id} = $gene_name."_"."NO"; 
									} else {
										$all_gene_id{$gene_id} = $gene_name."_"."YES"; 
									}
								}
								if(exists $CNV_on_samples{$sample}{$gene_id}) {
									$CNV_on_samples{$sample}{$gene_id}=$CNV_on_samples{$sample}{$gene_id}.";".$CNV_regions{$sample}{$chr}{$start_CNV}{$end_CNV};
								} else {
									$CNV_on_samples{$sample}{$gene_id} = $CNV_regions{$sample}{$chr}{$start_CNV}{$end_CNV};
								}
							} elsif ( $end_gene < $start_CNV) {
								next start_loop;
							} elsif ($start_gene > $end_CNV) {
								last start_loop;
							}
						}
					}
				}
			}
		}
	}
}


# output in excel format

 
print "Creating excel output...\n";
my @header = ('Chr', 'Start', 'Strand', 'Ref', 'Obs', 'SNPorINDEL', 'VariantCall_quality', 'Func', 'GeneName', 'GeneID', 'isInterested', 'ExonicFunc',	
'AAChange', 'Conserved', 'SegDup', 'ESP6500_ALL', '1000g2012feb_ALL', 'InHouseMAF', 'dbSNP135', 'AVSIFT', 'LJB_PhyloP', 'LJB_PhyloP_Pred', 'LJB_SIFT', 
'LJB_SIFT_Pred', 'LJB_PolyPhen2', 'LJB_PolyPhen2_Pred', 'LJB_LRT', 'LJB_LRT_Pred', 'LJB_MutationTaster', 'LJB_MutationTaster', 'Pred LJB_GERP++', 'Filter',	 'FORMAT')	;

foreach my $sample (@sample_names) {
	push @header, $sample;
	push @header, $sample.".anno";
	push @header, $sample.".cnv";
}
my @header_2 = ('INFO', 'GoTerm', 'WikiGene_Description', 'MIM_Gene_Description', 'GeneCard Link', 'OMIM Link', 'Uniprot Link');
push @header, @header_2;

my $workbook = Excel::Writer::XLSX->new($output_excel);

#format settings
my $format_header = $workbook->add_format();
$format_header->set_bold();

#my $format_possibleHits_Y = $workbook->add_format();
#$format_possibleHits_Y->set_color('lime');
#$format_possibleHits_Y->set_bg_color('green');

#my $format_possibleHits_N = $workbook->add_format();
#$format_possibleHits_N->set_color('red');

#my $format_possibleHits_C = $workbook->add_format();
#$format_possibleHits_C->set_bg_color('yellow');


my $worksheet_all = $workbook->add_worksheet('ALL');
$worksheet_all->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_all;$i++) {
	my @row = split("\t", $output_all[$i]);
	$worksheet_all->write_row($i+1,0,\@row);
}

my $worksheet_filtered = $workbook->add_worksheet('Filtered');
$worksheet_filtered->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_filtered;$i++) {
	my @row = split("\t", $output_filtered[$i]);
	$worksheet_filtered->write_row($i+1,0,\@row);
}

my $worksheet_xlinked = $workbook->add_worksheet('XLinked');
$worksheet_xlinked->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_Xlinked;$i++) {
	my @row = split("\t", $output_Xlinked[$i]);
	$worksheet_xlinked->write_row($i+1,0,\@row);
}

my $worksheet_possibleHits = $workbook->add_worksheet('Deleterious_Hits');
$worksheet_possibleHits ->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_possibleHits;$i++) {
	my @row = split("\t", $output_possibleHits[$i]);
	$worksheet_possibleHits->write_row($i+1,0,\@row);
}

my $worksheet_rareInterestedHits = $workbook->add_worksheet('Rare_Interested_Hits');
$worksheet_rareInterestedHits ->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @very_rare_interested_vars;$i++) {
	my @row = split("\t", $very_rare_interested_vars[$i]);
	$worksheet_rareInterestedHits->write_row($i+1,0,\@row);
}

my $worksheet_ARmodelHits = $workbook->add_worksheet('Homozygous_Hits');
$worksheet_ARmodelHits ->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_ARmodelHits;$i++) {
	my @row = split("\t", $output_ARmodelHits[$i]);
	$worksheet_ARmodelHits->write_row($i+1,0,\@row);
}

if ($CNV_file ne "") {
	my @CNV_header=('SampleID', 'start.p', 'end.p', 'type', 'nexons', 'start', 'end', 'chromosome', 'id', 'BF', 'reads.expected', 'reads.observed', 	
	'reads.ratio', 'Conrad.hg19', 'exons.hg19', 'Averagecontrols_DupDel_acrossCNV', 	
	'min_controls', 'max_controls', 'FullGenesNames', 'GO-terms', 'OMIM');
	my $worksheet_CNV = $workbook->add_worksheet('CNV');
	$worksheet_CNV -> write_row(0,0,\@CNV_header, $format_header);
	for (my $i=0;$i< scalar @output_CNV;$i++) {
		my @row = split("\t", $output_CNV[$i]);
		$worksheet_CNV -> write_row($i+1,0,\@row);
	}
}

############################# output compound heterozygous sheet ###################################################
my @compound_header = ('Gene_ID', 'Gene_Name', 'IsInterested');
foreach my $sample (@sample_names) {
	push @compound_header, $sample."_vars";
	push @compound_header, $sample."_CNV";
}
my $worksheet_compound = $workbook->add_worksheet('Compound_Heterozygous');
$worksheet_compound -> write_row(0,0,\@compound_header, $format_header);
my $number_of_wide_cols = (scalar @sample_names)*2+2;
$worksheet_compound -> set_column( 0, 0 , 19 );
$worksheet_compound -> set_column( 1, 1 , 13 );
$worksheet_compound -> set_column( 2, 2 , 13 );
$worksheet_compound -> set_column( 3, $number_of_wide_cols , 40 );


my %gene_sample_CNV_count; #number of how many CNVs on the gene
my %gene_sample_VAR_count; #number of how many snp/indel on the gene

foreach  my $sample (@sample_names) {
	foreach my $gene_id (keys %{$CNV_on_samples{$sample}}) {
		my @temp = split(";",$CNV_on_samples{$sample}{$gene_id});
		if(exists $gene_sample_CNV_count{$sample}{$gene_id}) {
			$gene_sample_CNV_count{$sample}{$gene_id} = $gene_sample_CNV_count{$sample}{$gene_id} + scalar @temp;
		} else {
			$gene_sample_CNV_count{$sample}{$gene_id} = scalar @temp;
		}
	}
}

foreach  my $sample (@sample_names) {
	foreach my $gene_id (keys %{$Variants_on_samples{$sample}}) {
		$gene_sample_VAR_count{$sample}{$gene_id} = 0;
		my @temp = split(";",$Variants_on_samples{$sample}{$gene_id});
		# Remove 'R/R' 'NA' and same position variants
		my %positions;
		foreach my $var_info (@temp) {
			if ($var_info !~ m/\_R\/R\_/ && $var_info !~ m/NA\_/) {
				#print "var is $var_info \n";
				my @temp_temp = split("_", $var_info);
				#print @temp_temp."\n";
				if (!exists $positions{$temp_temp[1]}) {
					$positions{$temp_temp[1]} = 0;
					if(exists $gene_sample_VAR_count{$sample}{$gene_id}) {
						#print "+1\n";
						$gene_sample_VAR_count{$sample}{$gene_id} = $gene_sample_VAR_count{$sample}{$gene_id} + 1;
					} 
				}
			}
		}
		#testing here
		#print $sample."_".$gene_id."\t".$Variants_on_samples{$sample}{$gene_id}." Count result ".$gene_sample_VAR_count{$sample}{$gene_id}."\n";
	}
}

# ready to output to excel sheet
my $excel_line_index=1;
my $compond_format= $workbook->add_format(text_wrap => 1);
foreach my $gene_id (keys %all_gene_id) {
	my @output;
	my $output_light = "red";
	push  @output, $gene_id;
	my @temp_yuanyuan = split("_", $all_gene_id{$gene_id});
	push  @output, $temp_yuanyuan[0];
	push  @output, $temp_yuanyuan[1];
	my $row_height =1;
	foreach my $sample (@sample_names) {
		if( exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_var = split(";", $Variants_on_samples{$sample}{$gene_id});
			my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
			my $temp = "";
			for(my $i=0;$i<scalar @gene_sample_var;$i++) {
				if ($i == scalar @gene_sample_var -1 ) {
					$temp = $temp.$gene_sample_var[$i];
				} else {
					$temp = $temp.$gene_sample_var[$i]."\n";
				}
			}
			push @output, $temp;
			$temp = "";
			for(my $i=0;$i<scalar @gene_sample_cnv;$i++) {
				if ($i == scalar @gene_sample_cnv -1 ) {
					$temp = $temp.$gene_sample_cnv[$i];
				} else {
					$temp = $temp.$gene_sample_cnv[$i]."\n";
				}
			}
			push @output, $temp;
			if ($gene_sample_CNV_count{$sample}{$gene_id} + $gene_sample_VAR_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
			
			# decided the largest row_height for excel cell
			if (scalar @gene_sample_cnv >= scalar @gene_sample_var) {
				if (scalar @gene_sample_cnv > $row_height) {
					$row_height = scalar @gene_sample_cnv;
				}
			} else {
				if (scalar @gene_sample_var > $row_height) {
					$row_height = scalar @gene_sample_var;
				}
			}
			
		} elsif (exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
			push @output, 'NA';
			my $temp = "";
			for(my $i=0;$i<scalar @gene_sample_cnv;$i++) {
				if ($i == scalar @gene_sample_cnv -1 ) {
					$temp = $temp.$gene_sample_cnv[$i];
				} else {
					$temp = $temp.$gene_sample_cnv[$i]."\n";
				}
			}
			push @output, $temp;
			if ($gene_sample_CNV_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
			if (scalar @gene_sample_cnv > $row_height) {
				$row_height = scalar @gene_sample_cnv;
			}
		} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
			my @gene_sample_var = split(";", $Variants_on_samples{$sample}{$gene_id});
			my $temp = "";
			for(my $i=0;$i<scalar @gene_sample_var;$i++) {
				if ($i == scalar @gene_sample_var -1 ) {
					$temp = $temp.$gene_sample_var[$i];
				} else {
					$temp = $temp.$gene_sample_var[$i]."\n";
				}
			}
			push @output, $temp;
			push @output, 'NA';
			if ($gene_sample_VAR_count{$sample}{$gene_id} >= 2) {
				$output_light = "green";
			}
			
			if (scalar @gene_sample_var > $row_height) {
				$row_height = scalar @gene_sample_var;
			}
		} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
			push @output, 'NA';
			push @output, 'NA';
		}
	}
	
	if ($output_light eq "green") {
		######## adjust row height according to the number of lines in the cell
		my $final_row_height = $row_height*17;
		$worksheet_compound -> set_row ($excel_line_index, $final_row_height);
		$worksheet_compound -> write_row($excel_line_index,0,\@output, $compond_format);
		$excel_line_index += 1;
	}
}
############################# output compound heterozygous sheet  end here ###################################################

############################# output compound heterozygous sheet for Autosomal recessive Model ###################################################
if (scalar @sample_names >= 2) {

	my $worksheet_compound_AR = $workbook->add_worksheet('Compound_Heterozygou_AR_Model');
	$worksheet_compound_AR -> write_row(0,0,\@compound_header, $format_header);
	$worksheet_compound_AR -> set_column( 0, 0 , 19 );
	$worksheet_compound_AR -> set_column( 1, 1 , 13 );
	$worksheet_compound_AR -> set_column( 2, 2 , 13 );
	$worksheet_compound_AR -> set_column( 3, $number_of_wide_cols , 40 );
	
	# ready to output to excel sheet
	$excel_line_index=1;
	foreach my $gene_id (keys %all_gene_id) {
		my @output;
		my $output_light = "red";
		my %gene_sample_count; # to record variants for each sample
		push  @output, $gene_id;
		my @temp_yuanyuan = split("_", $all_gene_id{$gene_id});
		push  @output, $temp_yuanyuan[0];
		push  @output, $temp_yuanyuan[1];
		my $row_height =1;
		foreach my $sample (@sample_names) {
			if( exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
				my @gene_sample_var = split(";", $Variants_on_samples{$sample}{$gene_id});
				my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
				my $temp = "";
				for(my $i=0;$i<scalar @gene_sample_var;$i++) {
					if ($i == scalar @gene_sample_var -1 ) {
						$temp = $temp.$gene_sample_var[$i];
					} else {
						$temp = $temp.$gene_sample_var[$i]."\n";
					}
				}
				push @output, $temp;
				$temp = "";
				for(my $i=0;$i<scalar @gene_sample_cnv;$i++) {
					if ($i == scalar @gene_sample_cnv -1 ) {
						$temp = $temp.$gene_sample_cnv[$i];
					} else {
						$temp = $temp.$gene_sample_cnv[$i]."\n";
					}
				}
				push @output, $temp;
				
				$gene_sample_count{$sample} = $gene_sample_CNV_count{$sample}{$gene_id} + $gene_sample_VAR_count{$sample}{$gene_id};
				
				# decided the largest row_height for excel cell
				if (scalar @gene_sample_cnv >= scalar @gene_sample_var) {
					if (scalar @gene_sample_cnv > $row_height) {
						$row_height = scalar @gene_sample_cnv;
					}
				} else {
					if (scalar @gene_sample_var > $row_height) {
						$row_height = scalar @gene_sample_var;
					}
				}
				
			} elsif (exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
				my @gene_sample_cnv = split(";", $CNV_on_samples{$sample}{$gene_id});
				push @output, 'NA';
				my $temp = "";
				for(my $i=0;$i<scalar @gene_sample_cnv;$i++) {
					if ($i == scalar @gene_sample_cnv -1 ) {
						$temp = $temp.$gene_sample_cnv[$i];
					} else {
						$temp = $temp.$gene_sample_cnv[$i]."\n";
					}
				}
				push @output, $temp;
				$gene_sample_count{$sample} = $gene_sample_CNV_count{$sample}{$gene_id};
				
				if (scalar @gene_sample_cnv > $row_height) {
					$row_height = scalar @gene_sample_cnv;
				}
			} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && exists $Variants_on_samples{$sample}{$gene_id}) {
				my @gene_sample_var = split(";", $Variants_on_samples{$sample}{$gene_id});
				my $temp = "";
				for(my $i=0;$i<scalar @gene_sample_var;$i++) {
					if ($i == scalar @gene_sample_var -1 ) {
						$temp = $temp.$gene_sample_var[$i];
					} else {
						$temp = $temp.$gene_sample_var[$i]."\n";
					}
				}
				push @output, $temp;
				push @output, 'NA';
				$gene_sample_count{$sample} = $gene_sample_VAR_count{$sample}{$gene_id};
				
				if (scalar @gene_sample_var > $row_height) {
					$row_height = scalar @gene_sample_var;
				}
			} elsif(!exists $CNV_on_samples{$sample}{$gene_id} && !exists $Variants_on_samples{$sample}{$gene_id}) {
				push @output, 'NA';
				push @output, 'NA';
				$gene_sample_count{$sample} = 0;
			}
		}
		
		# if each sample has more than 2 hits then output it
		$output_light = "green";
		#print $gene_id."\t";
		foreach my $sample_sample (keys %gene_sample_count) {
			if ($gene_sample_count{$sample_sample} < 2) {
				$output_light = "red";
			}
			#print $sample_sample."_".$gene_sample_count{$sample_sample}.$output_light."\t";
		}
		#print "\n";
		
		if ($output_light eq "green") {
			######## adjust row height according to the number of lines in the cell
			my $final_row_height = $row_height*17;
			$worksheet_compound_AR -> set_row ($excel_line_index, $final_row_height);
			$worksheet_compound_AR -> write_row($excel_line_index,0,\@output, $compond_format);
			$excel_line_index += 1;
		}
	}

}
############################# output compound heterozygous sheet for Autosomal recessive Model Ends Here ###################################################
$workbook->close();

exit;

sub ElementInArray {
	my ($array_ref, $item) = @_;
	my $result="NO";
	my @array_test = @{$array_ref};
	foreach my $element (@array_test) {
		if($element eq $item) {
			$result = "YES";
		}
	}
	return $result;
}

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: wAnnovarOutput_to_Excel.pl [--sampleNames Sample_1,Sample_2,...] [--sampleColumns 0-based_columnNumber_Sample_1,0-based_columnNumber_Sample_2,...] [--in wAnnovar output Tab-delimeted file] [--out output excel file name] [--CNV CNV result file][-help|-?]\n\n";
	return(1);
}

sub GetGeneCoords {
	my ($gene_file) = @_;
	my %gene_coords;
	open INPUT2, $genefile or die "Cannot open $genefile\n";
	while (my $Line = <INPUT2>){
		chomp $Line;
		my @linesplit1 = split(/\t/,$Line);
		if($linesplit1[0] eq 'MT'){
			$linesplit1[0]='M'
		}
		my $chr="chr".$linesplit1[0];
		my $st=$linesplit1[1];
		my $end=$linesplit1[2];
		my $gen=$linesplit1[3]."\t".$linesplit1[4]; # gene name \t gene id
		$gene_coords{$chr}{$st}{$end} = $gen;
	}
	close INPUT2;
	return \%gene_coords;
}
