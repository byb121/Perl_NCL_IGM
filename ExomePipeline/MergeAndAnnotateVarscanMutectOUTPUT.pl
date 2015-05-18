#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Excel::Writer::XLSX;
use Spreadsheet::ParseExcel;

print "\n";
print "\n";
print "############################################################################################\n";
print "############################################################################################\n";
print "\n";

my $varscan_in_files;
my $mutect_in_files;
my $output;
my $help;

usage() if ( @ARGV < 1 || ! GetOptions('help|?' =>\$help, "varscan=s"=>\$varscan_in_files, "mutect=s"=>\$mutect_in_files, 'out=s' => \$output) || defined $help );

unless (defined $varscan_in_files && defined $mutect_in_files && defined $output) {
	die "You have not supplied corret input. Read the help\n\n";
}

my %all_positions;

print "Reading in Varscan file(s)\n";
my @varscan_files=split(",", $varscan_in_files);
my %varscanV;
my @varscan_header;
foreach my $varscanFile (@varscan_files) {
	print "Reading in varscan file $varscanFile....\n";
	open VARSCAN, "$varscanFile" or die "Cannot open the varscan file $varscanFile.\n";
	while (my $line = <VARSCAN>) {
		chomp $line;
		if ($line =~ m/^\#/) {
			next;
		}
		if ($line =~ m/^chrom/ ) {
			if ( ! @varscan_header ) {
				my @temp = split("\t", $line);
				foreach my $ele (@temp) {push @varscan_header, $ele;}
			} else {
				next;
			}
		} else {
			my @temp = split("\t", $line);
			my $chr = $temp[0];
			$chr =~ s/^chr//;
			my $pos = $temp[1];
			my $ref = $temp[2];
			my $var = $temp[3];
			
			my $new_format_ref = VarscanToAnnovarFormat($pos, $ref, $var);
			my @new_format = @{ $new_format_ref };
			if (exists $varscanV{$chr}{$new_format[0]}{$new_format[1]}{$new_format[2]}{$new_format[3]}) {
				print "$chr new pos: $new_format[0] new ref: $new_format[2] new alt: $new_format[3]\n";
				print "Duplicated entry found in varscan file(s): $chr $pos $ref $var;\nexit now, bye!\n";
				exit;
			} else{
				my $string="";
				for(my $i=4;$i<scalar @temp;$i++) {
					$string = $string."\t".$temp[$i];
				}
				$string =~ s/^\t//;
				$varscanV{$chr}{$new_format[0]}{$new_format[1]}{$new_format[2]}{$new_format[3]} = $string;
				$all_positions{$chr}{$new_format[0]}{$new_format[1]}{$new_format[2]}{$new_format[3]} = 1; 
			}
			
		}
	}
	close(VARSCAN);
}
print "Varscan files are done, now let\'s read in mutect file(s)... \n";
print "\n";
print "Reading in mutect file(s)\n";

my @mutect_files=split(",", $mutect_in_files);
my %mutectV;
my @mutect_header;
foreach my $mutectFile (@mutect_files) {
	print "Reading in mutect file $mutectFile....\n";
	open MUTECT, "$mutectFile" or die "Cannot open the mutect file $mutectFile.\n";
	while (my $line = <MUTECT>) {
		chomp $line;
		if ($line =~ m/^\#/) {
			next;
		}
		if ($line =~ m/^contig/) {
			if ( ! @mutect_header) {
				my @temp = split("\t", $line);
				foreach my $ele (@temp) {push @mutect_header, $ele;}
			} else {
				next;
			}
		} else {
			my @temp = split("\t", $line);
			my $chr = $temp[0];
			$chr =~ s/^chr//;
			my $pos = $temp[1];
			my $ref = $temp[3];
			my $var = $temp[4];
			
			if ($ref =~ m/^[atgcATGC]$/ && $var =~ m/^[atgcATGC]$/ && $temp[scalar @temp -1] eq "KEEP") {
				if (exists $mutectV{$chr}{$pos}{$pos}{$ref}{$var}) {
					print "Duplicated entry found in mutect file(s): $chr $pos $ref $var;\nexit now, bye!\n";
					exit;
				} else{
					my $string=$temp[2];
					for(my $i=5;$i<scalar @temp;$i++) {
						$string = $string."\t".$temp[$i];
					}
					$string =~ s/^\t//;
					$mutectV{$chr}{$pos}{$pos}{$ref}{$var} = $string;
					$all_positions{$chr}{$pos}{$pos}{$ref}{$var} = 1; 
				}
			} #else {
			#	print "Since when mutect output indels????\nexit now!!\n";
			#	next;
			#}
		}
	}
	close(MUTECT);
}
print "Mutect files are done. Output Annovar input file now....\n";

my $annovar_input_file = $output.".annovar_input.txt";
open ANNOVAR_IN, ">$annovar_input_file" or die "Cannot open the file to output $annovar_input_file\n";

#print ANNOVAR_IN "chr\tstart\tend\tref\tvar";
#for(my $i=4;$i<scalar @varscan_header;$i++) {
#	print ANNOVAR_IN "\t".$varscan_header[$i];
#}
#print ANNOVAR_IN "\t".$mutect_header[2];
#for(my $i=5;$i<scalar @mutect_header;$i++) {
#	print ANNOVAR_IN "\t".$mutect_header[$i];
#}
#print ANNOVAR_IN "\n";

foreach my $chr ( sort {$a cmp $b} keys %all_positions ) {
	foreach my $start ( sort {$a<=>$b} keys %{$all_positions{$chr}} ) {
		foreach my $end ( keys %{$all_positions{$chr}{$start}} ) {
			foreach my $ref (keys %{$all_positions{$chr}{$start}{$end}} ) {
				foreach my $var (keys %{$all_positions{$chr}{$start}{$end}{$ref}} ) {
					if (exists $varscanV{$chr}{$start}{$end}{$ref}{$var} && exists $mutectV{$chr}{$start}{$end}{$ref}{$var} ) {
						print ANNOVAR_IN "$chr\t$start\t$end\t$ref\t$var\t$varscanV{$chr}{$start}{$end}{$ref}{$var}";
						print ANNOVAR_IN "\t$mutectV{$chr}{$start}{$end}{$ref}{$var}\n";
					} elsif (exists $varscanV{$chr}{$start}{$end}{$ref}{$var}) {
						print ANNOVAR_IN "$chr\t$start\t$end\t$ref\t$var\t$varscanV{$chr}{$start}{$end}{$ref}{$var}";
						for(my $i=5;$i<=scalar @mutect_header;$i++) { #print one more NA than the (header length - 5) because context column is in the firt 5 columns
							print ANNOVAR_IN "\t"."NA";
						}
						print ANNOVAR_IN "\n";
					} elsif (exists $mutectV{$chr}{$start}{$end}{$ref}{$var} ) {
						print ANNOVAR_IN "$chr\t$start\t$end\t$ref\t$var";
						for(my $i=4;$i<scalar @varscan_header;$i++) {
							print ANNOVAR_IN "\t"."NA";
						}
						print ANNOVAR_IN "\t$mutectV{$chr}{$start}{$end}{$ref}{$var}\n";
					} else {
						print "Impossible error. Exit\n";
						exit;
					}
				}
			}
		} 
	}
}
close ANNOVAR_IN;	

################################################################################ run annovar
`sh /sharedlustre/IGM/annovar_2014jul14/run_annovar_tumor.sh $annovar_input_file $annovar_input_file`;
#################################################################################################



#add in house annotation
my $annovar_output_file = "$annovar_input_file.hg19_multianno.txt";
my $inhouse_annotated_file =  AddInHouseAnnotationToAnnovarOutput ($annovar_output_file);

print "start filtering and produce filtered table\n";

my @output_everything;
my @output_filtered;
open INHOUSE, "$inhouse_annotated_file" or die "Cannot open the file $inhouse_annotated_file\n";
while (my $line=<INHOUSE>) {
	if ($line !~ m/^Chr\t/) {
		chomp $line;
		my @elements = split("\t", $line);
		my $length = scalar @elements -1;
		$elements[$length-2] =~ s/\'/\"/g;
		$elements[$length-1] =~ s/\'/\"/g;
		$elements[$length] =~ s/\'/\"/g;
		
		my $string_1 = $elements[0]."\t".$elements[1]."\t".$elements[3]."\t".$elements[4]."\t".$elements[119]."\t".$elements[5]."\t".$elements[120]
		."\t".$elements[121]."\t".$elements[8]."\t".$elements[9]."\t".$elements[30].":".$elements[32].":".$elements[35].":".$elements[38].":"
		.$elements[41].":".$elements[44].":".$elements[47].":".$elements[49]."\t".$elements[56]."\t".$elements[122]."\t".$elements[25]
		."\t".$elements[53]."\t".$elements[54]."\t".$elements[55]."\t".$elements[66]."\t".$elements[67]."\t".$elements[68];
		
		my $variant_caller="";
		if ($elements[69] ne "NA" && $elements[88] ne "NA") {
			$variant_caller="MU,VAR";
		} elsif ($elements[69] ne "NA") {
			$variant_caller="VAR";
		} else {
			$variant_caller="MU";
		}
		
		$string_1 = $string_1."\t".$variant_caller."\t".$elements[69]."\t".$elements[70]."\t".$elements[71]."\t".$elements[72]."\t".$elements[73]
		."\t".$elements[74]."\t".$elements[75]."\t".$elements[76]."\t".$elements[78]."\t".$elements[79]."\t".$elements[114]."\t".$elements[115]
		."\t".$elements[104]."\t".$elements[105]."\t".$elements[94]."\t".$elements[100]."\t".$elements[112]."\t".$elements[80]."\t".$elements[81]
		."\t".$elements[82]."\t".$elements[83]."\t".$elements[84]."\t".$elements[85]."\t".$elements[86]."\t".$elements[87]."\t".$elements[77].
		"\t".$elements[88]."\t".$elements[93]
		."\t".$elements[95]."\t".$elements[96]."\t".$elements[97]."\t".$elements[98]."\t".$elements[99]."\t".$elements[101]."\t".$elements[102]
		."\t".$elements[103]."\t".$elements[106]."\t".$elements[107]."\t".$elements[108]."\t".$elements[109]."\t".$elements[110]."\t".$elements[111]
		."\t".$elements[113]."\t".$elements[116]."\t".$elements[117]."\t".$elements[118]."\t".$elements[20]."\t".$elements[21]."\t".sprintf( "%.4f",$elements[22])
		."\t".sprintf( "%.4f",$elements[23])."\t".sprintf( "%.4f",$elements[24])
		."\t".$elements[127]."\t".$elements[128]."\t".$elements[129]."\t".$elements[123]."\t".$elements[124]
		."\t".$elements[125]."\t".$elements[126]."\t".$elements[10]."\t".$elements[13]."\t".$elements[14]."\t".$elements[15]."\t".$elements[18]
		."\t".$elements[19]
		."\t".$elements[27]."\t".$elements[29]."\t".$elements[30]
		."\t".$elements[31]."\t".$elements[32]."\t".$elements[34]."\t".$elements[35]."\t".$elements[37]
		."\t".$elements[38]."\t".$elements[40]."\t".$elements[41]."\t".$elements[43]."\t".$elements[44]
		."\t".$elements[46]."\t".$elements[47]."\t".$elements[48]."\t".$elements[49]."\t".$elements[50]."\t".$elements[51]."\t".$elements[52];
		
		push @output_everything, $string_1;
		
		if($elements[8] eq "synonymous SNV" && $elements[13] eq "synonymous SNV" && $elements[18] eq "synonymous SNV"  ) { # eliminate synonymous SNV
			next;
		}
		
		if($elements[5] =~ m/(intronic|intergenic|upstream|downstream)/i &&  $elements[10] =~ m/(intronic|intergenic|upstream|downstream)/i 
						&& $elements[15] =~ m/(intronic|intergenic|upstream|downstream)/i  ) { # eliminate intronic or intergenic variants
			next;
		}
		
		if( $elements[56] <= 0.05 && $elements[122] <= 0.05) {
			push @output_filtered, $string_1;
		} else {
			next;
		}
	}
}
close INHOUSE;
print "Reading in inhouse annotation output is done.\n";
print "Output results in excel.\n";

my @header = 
('Chr', 'Start', 'Ref', 'Alt', 'SNPorINDEL', 'Func.ensGene', 'Gene_names', 
'Gene_ID', 'ExonicFunc.ensGene', 'AAChange.ensGene', 'predictions', 'PopFreqMax', 
'maf_GATK', 'snp138', 'clinvar_20140303', 'caddgt10', 'cosmic68', 'gwasCatalog', 
'dgvMerged', 'tfbsConsSites', 'Variant_Caller', 'normal_reads1', 'normal_reads2', 
'normal_var_freq', 'normal_gt', 'tumor_reads1', 'tumor_reads2', 'tumor_var_freq', 
'tumor_gt', 'variant_p_value', 'somatic_p_value', 'MU_n_ref_count', 'MU_n_alt_count', 
'MU_t_ref_count', 'MU_t_alt_count', 'MU_power', 'MU_t_lod_fstar', 'MU_normal_best_gt', 
'tumor_reads1_plus', 'tumor_reads1_minus', 'tumor_reads2_plus', 'tumor_reads2_minus', 
'normal_reads1_plus', 'normal_reads1_minus', 'normal_reads2_plus', 'normal_reads2_minus', 
'somatic_status', 'MU_context', 'MU_covered', 'MU_tumor_power', 'MU_normal_power', 
'MU_total_pairs', 'MU_improper_pairs', 'MU_map_Q0_reads', 'MU_tumor_f', 
'MU_contaminant_fraction', 'MU_contaminant_lod', 'MU_t_ref_sum', 'MU_t_alt_sum', 
'MU_t_ref_max_mapq', 'MU_t_alt_max_mapq', 'MU_t_ins_count', 'MU_t_del_count', 
'MU_init_n_lod', 'MU_n_ref_sum', 'MU_n_alt_sum', 'MU_judgement', 'phastConsElements46way', 
'genomicSuperDups', 'esp6500si_all', '1000g2012apr_all', 'cg69', 'GeneCard Link', 
'OMIM Link', 'Uniprot Link', 'GoTerm', 'WikiGene_Description', 'MIM_Gene_Description', 
'OMIM_Gene_Description', 'Func.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 
'Func.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
'LJB23_SIFT_score_converted', 'LJB23_Polyphen2_HDIV_score', 'LJB23_Polyphen2_HDIV_pred', 
'LJB23_Polyphen2_HVAR_score', 'LJB23_Polyphen2_HVAR_pred', 'LJB23_LRT_score_converted',
 'LJB23_LRT_pred', 'LJB23_MutationTaster_score_converted', 'LJB23_MutationTaster_pred', 
 'LJB23_MutationAssessor_score_converted', 'LJB23_MutationAssessor_pred', 'LJB23_FATHMM_score_converted', 
 'LJB23_FATHMM_pred', 'LJB23_RadialSVM_score_converted', 'LJB23_RadialSVM_pred', 'LJB23_LR_score', 
 'LJB23_LR_pred', 'LJB23_GERP++', 'LJB23_PhyloP', 'LJB23_SiPhy');

my $output_excel = $output.".xlsx";
my $workbook = Excel::Writer::XLSX->new($output_excel);
my $superDupsRow = 68;

#format settings
my $format_header = $workbook->add_format();
$format_header->set_bold();
my $format_segmental_dup = $workbook->add_format(color => 10);

my $worksheet_everything = $workbook->add_worksheet('ALL');
$worksheet_everything->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_everything;$i++) {
	my @row = split("\t", $output_everything[$i]);
	if ($row[$superDupsRow] ne "NA") {
		$worksheet_everything->write_row($i+1,0,\@row, $format_segmental_dup);
	} else {
		$worksheet_everything->write_row($i+1,0,\@row);
	}
}

my $worksheet_filtered = $workbook->add_worksheet('RareExonicVariants');
$worksheet_filtered->write_row(0,0,\@header, $format_header);
for (my $i=0;$i< scalar @output_filtered;$i++) {
	my @row = split("\t", $output_filtered[$i]);
	if ($row[$superDupsRow] ne "NA") {
		$worksheet_filtered->write_row($i+1,0,\@row, $format_segmental_dup);
	} else {
		$worksheet_filtered->write_row($i+1,0,\@row);
	}
}

$workbook->close();



print "\n";



exit;


sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: MergeAndAnnotateVarscanMutectOUTPUT.pl \n";
    print "--varscan input varscan files (multiple files can be separated with \",\" );\n";
    print "--mutect input mutect files (multiple files can be separated with \",\" );\n"; 
    print "--out output excel file name.\n";
    return(1);
}

sub VarscanToAnnovarFormat {
	my ($pos, $ref, $var) = @_;
	if ($ref =~ m/^[atgcATGC]$/ && $var =~ m/^[atgcATGC]$/) {
		my @out = ($pos, $pos, $ref, $var);
		return \@out;
	} else {
		if ($var =~ /^\+/) {
			$var =~ s/^\+//;
			my $start = $pos+1;
			my $end = $pos+1;
			my $ref_new="-";
			my $var_new = $var;
			my @out = ($start, $end, $ref_new, $var_new);
			return \@out;
		} elsif($var =~ /^\-/) {
			$var =~ s/^\-//;
			my $len = length($var);
			my $start = $pos+1;
			my $end = $pos+$len;
			my $ref_new = $var;
			my $var_new = "-";
			my @out = ($start, $end, $ref_new, $var_new);
			return \@out;
		} else {
			print "fatal error of sub VarscanToAnnovarFormat: Cannot process -- $pos, $ref, $var.\n";
		}
	}
}

sub AddInHouseAnnotationToAnnovarOutput {
	#will change NA in MAFs to 0
	#will add SNPorINDEL, inhouse GATK mafs and hyperlinks 
	my ($annover_output) = @_;
	my $output_file ="$annover_output.inhouseAnnotated.txt"; 
	
	my $genefile="/sharedlustre/IGM/annovar_2014jul14/inHouse_db/Ensembl_Genes_GRCh37_75.txt";
	my $sharedfile_GATK="/sharedlustre/IGM/annovar_2014jul14/inHouse_db/InHouse_OnTarget_GATK_MAFs.txt_179exomes";
	my $OMIMfile="/sharedlustre/IGM/annovar_2014jul14/inHouse_db/Ensembl_OMIM_AllGeneNames.txt";
	my $ens_gene_OMIM_Uniprot_Acc = "/sharedlustre/IGM/annovar_2014jul14/inHouse_db/ens_gene_symbol_omim_id_uniprot_id.txt";
	
	print "Require db file: Ensembl_Genes_GRCh37_75.txt\n";
	if (-e $genefile) {
		print "Found $genefile\n";
		print "Consider update when necessary\n";
	} else {
		print "$genefile not exist\n exit\n";
		exit;
	}
	
	print "Require db file: InHouse_OnTarget_GATK_MAFs.txt\n";
	if (-e $sharedfile_GATK) {
		print "Found $sharedfile_GATK\n";
		print "Consider update when necessary\n";
	} else {
		print "$sharedfile_GATK not exist\n exit\n";
		exit;
	}
	
	print "Require db file: Ensembl_OMIM_AllGeneNames.txt\n";
	if (-e $OMIMfile) {
		print "Found $OMIMfile\n";
		print "Consider update when necessary\n";
	} else {
		print "$OMIMfile not exist\n exit\n";
		exit;
	}
	
	print "Require db file: ens_gene_symbol_omim_id_uniprot_id.txt\n";
	if (-e $ens_gene_OMIM_Uniprot_Acc) {
		print "Found $ens_gene_OMIM_Uniprot_Acc\n";
		print "Consider update when necessary\n";
	} else {
		print "$ens_gene_OMIM_Uniprot_Acc not exist\n exit\n";
		exit;
	}
	
	my $en_genes_ref = GetGeneCoords($genefile);
	my %en_genes = %$en_genes_ref;
	my $ens_symbol_map_ref = GetEnsSymbolMap($genefile);
	my %ens_symbol_map = %$ens_symbol_map_ref;
	
	my $inHouse_MAF_GATK_ref = GetInHouseMafGATK($sharedfile_GATK);
	my %inHouse_MAF_GATK = %$inHouse_MAF_GATK_ref;
	my $OMIM_ref = GetOMIManno($OMIMfile);
	my %OMIM = %$OMIM_ref;
	
	my $OmimAcc_ref = GetOmimAcc($ens_gene_OMIM_Uniprot_Acc);
	my %OmimAcc = %$OmimAcc_ref;
	my $UniprotAcc_ref = GetUniprotAcc($ens_gene_OMIM_Uniprot_Acc);
	my %UniprotAcc = %$UniprotAcc_ref;
		
	
	my @output;
	open ANNOVAR_OUTPUT , "$annovar_output_file" or die "Can not open the file $annovar_output_file\n";
	while (my $line = <ANNOVAR_OUTPUT> ) {
		chomp $line;
		if($line =~ m/^Chr/) {
			push @output, $line."\n";
			next;
		} else {
			my @elements  = split ("\t", $line);

			my $annovar_chr=$elements[0];
			my $annovar_pos=$elements[1];
			my $annovar_R=$elements[3];
			my $annovar_A=$elements[4];
			
			my $SNPorINDEL;
			if ($annovar_R =~ m/^[atgcATGC]$/ && $annovar_A =~ m/^[atgcATGC]$/) { # a SNP with just one alternative allele
				$SNPorINDEL = "SNP"; 
			} else {
				$SNPorINDEL="INDEL";
			}
			
			my $genecard_link = "";
			my $gene_name;
			my $ens_id;
			if ($elements[5] eq "intergenic") {
				$ens_id = "NA";
			} elsif ($elements[5] eq "splicing") {
				if ($elements[6] =~ m/^(.*?)\(.*\)$/) {
					$ens_id = $1;
				} elsif ($elements[6] !~ m/(\(|\,)/) {
					$ens_id = $elements[6];
				} elsif ($elements[6] =~ m/(EN\w+\d+?)\,EN.*/) {
					$ens_id = $1;
				} else {
					print "fatal match error on ensembl annotation splicing match;\n";
					print "$line.\n\n";
					exit 2;
				}
			} elsif ($elements[5] eq 'exonic;splicing') {
				if ($elements[6] =~ m/^(.*?)\;/) {
					$ens_id = $1;
				} else {
					$ens_id = $elements[6];
				}
			} else {
				if ($elements[6] ne "NA" &&  $elements[6] !~ m/(\(|\,)/) {
					$ens_id = $elements[6];	
				} else {
					$ens_id = "NA";
				}
			}
			
			if ($ens_id eq "NA") {
				$genecard_link = "NA";
			} else {
				$genecard_link = "=HYPERLINK(\'http://www.genecards.org/cgi-bin/carddisp.pl?id=$ens_id&id_type=ensembl\', \'GeneCard Link\')";
			}
			
			if (exists $ens_symbol_map{$ens_id}) {
				$gene_name = $ens_symbol_map{$ens_id};
			} else {
				$gene_name = "NA";
			}
			
			my $omim_link = "";
			if (exists $OmimAcc{$ens_id}) {
				my $omim_accesion =  $OmimAcc{$ens_id};
				$omim_link = "=HYPERLINK(\'http://omim.org/entry/$omim_accesion\', \'OMIM Link\')";
			} else {
				$omim_link = "NA";
			}
			
			my $uniprot_link = "";
			if (exists $UniprotAcc{$ens_id}) {
				my $uniprot_accession = $UniprotAcc{$ens_id};
				$uniprot_link = "=HYPERLINK(\'http://www.uniprot.org/uniprot/$uniprot_accession\', \'Uniprot Link\')";
			} else {
				$uniprot_link = "NA";
			}
			
			my $omim_anno = "NA\tNA\tNA\tNA";
			if (exists $OMIM{$ens_id}) {
				$omim_anno = $OMIM{$ens_id};
			}
			
			my $maf_GATK = "";
			if (exists $inHouse_MAF_GATK{$annovar_chr}{$annovar_pos}{$annovar_R}{$annovar_A}) {
				$maf_GATK = $inHouse_MAF_GATK{$annovar_chr}{$annovar_pos}{$annovar_R}{$annovar_A};
			} else {
				$maf_GATK = "0";
			}
			
			if($elements[22] eq "NA" || $elements[22] eq "." ) {#1000G
				$elements[22] = 0;
			}
			if($elements[23] eq "NA" || $elements[23] eq ".") { #esp6500
				$elements[23] = 0;
			}
			if($elements[24] eq "NA" || $elements[24] eq "." ) { #cg69
				$elements[24] = 0;
			}
			if($elements[56] eq "NA" ) { #PopFreqMax
				$elements[56] = 0;
			}
			if($elements[54] eq "NA" ) { #CADD gt 10
				$elements[54] = 0;
			}
			
			for (my $i=0;$i< scalar @elements;$i++){ push @output, $elements[$i]."\t";}
			push @output, $SNPorINDEL ."\t".$gene_name."\t".$ens_id."\t".$maf_GATK."\t".
					$omim_anno."\t".$genecard_link."\t".$omim_link."\t".$uniprot_link."\n";
		
		}
	}
	close ANNOVAR_OUTPUT;
	
	open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
	print OUTPUT @output;
	close OUTPUT;
	return $output_file;
}

sub GetOmimAcc {
	my ($sharedfile) = @_;
	my %OmimAcc;
	open INPUT, $sharedfile or die "Cannot open $sharedfile\n";
	while (my $line = <INPUT>) {
		if ($line =~ m/^ENSG/) {
			chomp $line;
			my @elements = split(/\t/, $line ,-1);
			if ($elements[2] ne "") {
				$OmimAcc{$elements[0]} = $elements[2];
			}
		} else {
			next;
		}
	}
	close INPUT;
	return \%OmimAcc;
}

sub GetUniprotAcc {
	my ($sharedfile) = @_;
	my %UniprotAcc;
	open INPUT, $sharedfile or die "Cannot open $sharedfile\n";
	while (my $line = <INPUT>) {
		if ($line =~ m/^ENSG/) {
			chomp $line;
			my @elements = split(/\t/, $line, -1);
			if ($elements[3] ne "") {
				$UniprotAcc{$elements[0]} = $elements[3];
			}
		} else {
			next;
		}
	}
	close INPUT;
	return \%UniprotAcc;
}

sub GetInHouseMafGATK {
	my ($sharedfile) = @_;
	my %IHcontrols;
	open INPUT3, $sharedfile or die "Cannot open $sharedfile\n"; 
	shloop: while (my $Line = <INPUT3>){
		chomp $Line;
		my @linesplit2 = split(/\t/,$Line);
		my $c=$linesplit2[0];
		my $p=$linesplit2[1];
		my $r=$linesplit2[2];
		my $v=$linesplit2[3];
		my $maf=$linesplit2[7];
		$IHcontrols{$c}{$p}{$r}{$v}=$maf;
	}
	close INPUT3;
	return \%IHcontrols;
}

sub GetOMIManno {
	my ($OMIMfile) = @_;
	my %OMIM;
	my %add_OMIM;
	my %seen_already;
	open INPUT, $OMIMfile or die "Cannot open $OMIMfile\n"; 
	while (my $Line = <INPUT>){
		chomp $Line;
		my @linesplit = split(/\t/,$Line);
		my $ens_id = $linesplit[0];
		my $go;
		my $wiki;
		my $MIM;
		my $oMIM;
		if ( exists $linesplit[2] ) {
			if ($linesplit[2] eq "") {
				$go = "NA";
			} else {
				$go = $linesplit[2];
			}
		} else {
			$go = "NA";
		}
			
		
		if ( exists $linesplit[3] ){
			if ($linesplit[3] eq "") {
				$wiki = "NA";
			} else {
				$wiki = $linesplit[3];
			}
		} else {
			$wiki = "NA";
		}
		
		if ( exists $linesplit[4] ){
			if ($linesplit[4] eq "") {
				$MIM = "NA";
			} else {
				$MIM = $linesplit[4];
			}
		} else {
			$MIM = "NA";
		}
		
		if ( exists $linesplit[5] ){
			if ($linesplit[5] eq "") {
				$oMIM = "NA";
			} else {
				$oMIM = $linesplit[5];
			}
		} else {
			$oMIM = "NA";
		}
		
		#group all individual terms per gene
		
		if(!exists $add_OMIM{$ens_id}){
			$add_OMIM{$ens_id}{'go'}="";
			$add_OMIM{$ens_id}{'wiki'}="";
			$add_OMIM{$ens_id}{'MIM'}="";
			$add_OMIM{$ens_id}{'oMIM'}="";
			}
		
		if(!exists $seen_already{$ens_id}{$go}){
		$add_OMIM{$ens_id}{'go'} = $add_OMIM{$ens_id}{'go'}.", ".$go;
		$seen_already{$ens_id}{$go}=0;
		}
		if(!exists $seen_already{$ens_id}{$wiki}){
		$add_OMIM{$ens_id}{'wiki'} = $add_OMIM{$ens_id}{'wiki'}.", ".$wiki;
		$seen_already{$ens_id}{$wiki}=0;
		}
		if(!exists $seen_already{$ens_id}{$MIM}){
		$add_OMIM{$ens_id}{'MIM'} = $add_OMIM{$ens_id}{'MIM'}.", ".$MIM;
		$seen_already{$ens_id}{$MIM}=0;
		}
		if(!exists $seen_already{$ens_id}{$oMIM}){
		$add_OMIM{$ens_id}{'oMIM'} = $add_OMIM{$ens_id}{'oMIM'}.", ".$oMIM;
		$seen_already{$ens_id}{$oMIM}=0;
		}
	}
	close INPUT;
		
	#loop to add all individual terms per ENSID together
	foreach my $e (keys %add_OMIM){
		$add_OMIM{$e}{'go'}=~s/^\,\s//;
		$add_OMIM{$e}{'wiki'}=~s/^\,\s//;
		$add_OMIM{$e}{'MIM'}=~s/^\,\s//;
		$add_OMIM{$e}{'oMIM'}=~s/^\,\s//;
		
		$OMIM{$e} = $add_OMIM{$e}{'go'}."\t".$add_OMIM{$e}{'wiki'}."\t".$add_OMIM{$e}{'MIM'}."\t".$add_OMIM{$e}{'oMIM'};
	}

	return \%OMIM;
}


sub GetGeneCoords {
	my ($genefile) = @_;
	my %gene_coords;
	open INPUT, $genefile or die "Cannot open $genefile\n";
	while (my $Line = <INPUT>){
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
	close INPUT;
	return \%gene_coords;
}

sub GetEnsSymbolMap{
	my ($genefile) = @_;
	my %map;
	open INPUT, $genefile or die "Cannot open $genefile\n";
	while (my $Line = <INPUT>){
		chomp $Line;
		my @linesplit1 = split(/\t/,$Line);
		$map{$linesplit1[4]}=$linesplit1[3]; #ens_id = symbol
	}
	close INPUT;
	return \%map;
}

