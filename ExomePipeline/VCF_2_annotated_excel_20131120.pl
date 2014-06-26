#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Excel::Writer::XLSX;


print "\n";
print "VCF_2_annotated_excel_20131120.pl says\n";
print "############################################################################################\n";
print "# YX: in house annotation part was amended from HG's                                       #\n";
print "# YX: the script will invoke Annovar and produce annotated output in Excel XLSX format     #\n";
print "# YX: input shoud be a single vcf file of a family or single sample                        #\n";                                 
print "# HG edits: to take (standalone) annovar output (with vcf file lines appended after)       #\n";
print "# vcf output (i.e. chr) starts at column 37 (fileline array element 36)                    #\n";
print "# column number of 'vcf' chr will depend on how many annotations are included from annovar #\n";
print "############################################################################################\n";
print "\n";

my $vcf_in;
my $InterestedGenefile=""; #Ensembl gene ids on each line
my $CNV_file="";
my $output_excel=""; #if empty then inHouse annotated file will not converted to excel file
my $help;

usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, "vcf=s"=>\$vcf_in, "InterestedGenes=s"=>\$InterestedGenefile, 'CNV=s' => \$CNV_file, 'out=s' => \$output_excel) || defined $help );

unless (defined $vcf_in && -e $vcf_in) {
	die "You have not supplied input vcf file using --vcf or the file does not exist\n\n";
}

print "Start to run Annovar on $vcf_in\n";
`/users/a5907529/data/Files_HG/vcf_annotation_november2013/annovar_YX.sh $vcf_in`;
print "Annovar is done!\n\n";

my $annovar_output_filename = $vcf_in.".hg19_multianno.txt";
my $vcf = $annovar_output_filename;
my $output_file=$annovar_output_filename."_inHouse_annotated.txt"; #in house annotation script output file

my $input_file;
if ($output_excel ne "") {$input_file=$output_file;} #input file for the script part to convert it to excel file

#### input part for converting excel file part
my $bash_output=`grep '#CHR' $vcf_in`;
chomp $bash_output;
my @sample_names;
my @sample_columns;
my @bash_output_split = split("\t", $bash_output);
my $first_sample_column_index = 47; # change if add or reduced Annovar annotation 0_based
# it is the sample column not the genotype column eg: 0/0:..... not the R/R column

print "Looking for sample names:\n";
for(my $i=9;$i<scalar @bash_output_split;$i++){
	print "Found: ".$bash_output_split[$i]."\n";
	push @sample_names, $bash_output_split[$i];
}

for(my $i=1;$i<=scalar @sample_names;$i++){
	push @sample_columns, $first_sample_column_index;
	$first_sample_column_index += 2;
}


#### Requiring input in house db files to add as annotation
my $genefile="/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/Ensembl_Genes_R67.txt";
my $sharedfile="/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/InHouse_OnTarget_Variants_MAFs.txt_394exomes";
my $sharedfile_GATK="/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/InHouse_OnTarget_GATK_MAFs.txt_297exomes";
my $OMIMfile="/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/Ensembl_OMIM_AllGeneNames.txt";
my $ens_gene_OMIM_Uniprot_Acc = "/users/data/Files_HG/vcf_annotation_november2013/inHouse_db/ens_gene_symbol_omim_id_uniprot_id.txt";

print "\n";
print "Require db file: Ensembl_Genes_R67.txt\n";
if (-e $genefile) {
	print "Found $genefile\n";
	print "Consider update when necessary\n";
} else {
	print "$genefile not exist\n exit\n";
	exit;
}

print "Require db file: InHouse_OnTarget_Variants_MAFs.txt\n";
if (-e $sharedfile) {
	print "Found $sharedfile\n";
	print "Consider update when necessary\n";
} else {
	print "$sharedfile not exist\n exit\n";
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

# An comlumn to tell if the variant is a SNP or INDEL
my $SNPorINDEL;

### Read in all annotations first ###
my %en_genes;
my %inHouse_MAF;
my %inHouse_MAF_GATK;
my %OMIM;
my %isInterestedGenes;
my %OmimAcc;
my %UniprotAcc;

my $en_genes_ref = GetGeneCoords($genefile);
%en_genes = %$en_genes_ref;
my $inHouse_MAF_ref = GetInHouseMaf($sharedfile);
%inHouse_MAF = %$inHouse_MAF_ref;
my $inHouse_MAF_GATK_ref = GetInHouseMafGATK($sharedfile_GATK);
%inHouse_MAF_GATK = %$inHouse_MAF_GATK_ref;
my $OMIM_ref = GetOMIManno($OMIMfile);
%OMIM = %$OMIM_ref;
if ($InterestedGenefile ne "") {
	my $isInterestedGenes_ref = GetIsInterestedGenes($InterestedGenefile);
	%isInterestedGenes = %$isInterestedGenes_ref;
}
my $OmimAcc_ref = GetOmimAcc($ens_gene_OMIM_Uniprot_Acc);
%OmimAcc = %$OmimAcc_ref;
my $UniprotAcc_ref = GetUniprotAcc($ens_gene_OMIM_Uniprot_Acc);
%UniprotAcc = %$UniprotAcc_ref;

my @output;

open VCF, "$vcf" or die "Can not open the file $vcf";
while (my $line = <VCF> ) {
	chomp $line;
	if($line =~ m/^Chr/) {
		push @output, $line."\n";
		next;
	} else {
		my @elements  = split ("\t", $line);
		$SNPorINDEL = "SNP"; # this will add a column to indicate the variant is a SNP/indel; easy for filtering 
		my $chr=$elements[38];
		my $pos=$elements[39];
		my $R=$elements[41];
		my $A=$elements[42];
		my @variants; # chr \t pos \t ref \t v \t SNPorINDEL # to record variant ( variants, if multi alternative alleles ) on the position
		my @ALTs; # record variants in mutation taster format. it's to query in house MAF
		
		my $FORMAT=$elements[46];
		my $Sample_Call="";
		for(my $i=47;$i<scalar(@elements);$i++) {
			$Sample_Call=$Sample_Call."\t".$elements[$i];
		}
		$Sample_Call =~ s/^\t//;
		
		#print "test: chr $chr pos $pos Ref $R Var $A\n";
		if ($R =~ m/^[atgcATGC]$/ && $A =~ m/^[atgcATGC]$/) { # a SNP with just one alternative allele
			push @variants, "$chr\t$pos\t$R\t$A";
			push @ALTs, $A;
		} elsif ($R =~ m/^[atgcATGC]$/) { # insertions
			if($A !~ /\,/){ #just one alternative allele
				my $v = $A;
				if ( length $v > 1 ) {
					$v =~ s/^./\+/;
					$SNPorINDEL="INDEL";
				}
				push @variants, "$chr\t$pos\t$R\t$v\t$SNPorINDEL";
				push @ALTs, $A;
			} else { # multi ALT bases insertions
				my @vars=split(/\,/,$A);
				foreach my $v (@vars) {
					push @ALTs, $v;
					if ( length $v > 1 ) {
						$v =~ s/^./\+/;
						$SNPorINDEL="INDEL";
					}
					push @variants, "$chr\t$pos\t$R\t$v\t$SNPorINDEL";
				}
			}
		} else { ##deletion(s) !!and in some cases of multiple (,) vars insertions!! 
			$SNPorINDEL="INDEL";
			#print "A:     $A\n";
			if ( $A !~ /\,/) {
					my $v = $R;
					if ( length $v > 1 ) {
						$v =~ s/^./\-/;
					}
					push @variants, "$chr\t$pos\t$R\t$v\t$SNPorINDEL";
					push @ALTs, $A;
			} else {
				my @vars=split(/\,/,$A);
				foreach my $v (@vars) {
					if 	(length $v < length $R){
						my $diff=(length $v)-1; 
						my $alt = substr($R,$diff);
						my $ref = substr($R,$diff, 1);
						$alt =~s/^./\-/; 
						$pos = $pos+$diff;
						push @ALTs, $v;
						push @variants, "$chr\t$pos\t$ref\t$alt\t$SNPorINDEL";
					} elsif(length $v > length $R){
						my $diff=(length $R)-1; 
						my $alt = substr($v,$diff);
						my $ref = substr($v,$diff,1); 
						$alt =~s/^./\+/;
						$pos = $pos+$diff;
						push @ALTs, $v;
						push @variants, "$chr\t$pos\t$ref\t$alt\t$SNPorINDEL";					
					} else {
						print "Warning: Equal length of REF and  VAR were found on line: \n".$line."\n Considering it's a SNP on the first base of RAF\n";
						my $alt =  substr($v, 0, 1);
						my $ref = substr($R, 0, 1);
						push @ALTs, $v;
						push @variants, "$chr\t$pos\t$ref\t$alt\tSNP";
						print "Use the position recorded on the vcf line, in house MAF might be wrong. Moving to the next line now.\n";
					}
				}
			}
		}
		
		# find which gene the variant is on ...
		#.... This will only add 1x gene ... some chr/pos match to more than 1x gene (forward & reverse DNA strands) ... 
		# ... possibly better to get annovar annotation first and then add additional info based on most likely deleterious gene according to annovar annotation!!
		my $vcf_pos=$elements[1];
		my $gene_name_ens_id = "";
		my $genecard_link = "";
		start_loop: foreach my $start ( sort {$a<=>$b} keys %{$en_genes{$chr}} ) {
			foreach my $end (keys %{$en_genes{$chr}{$start}}) {
				if ( $vcf_pos >= $start && $vcf_pos <= $end) {
					$gene_name_ens_id = $en_genes{$chr}{$start}{$end};
					last start_loop;
				}
			}
		}
		my $gene_name;
		my $ens_id;
		if ( $gene_name_ens_id ne "" ) {
			my @temp = split("\t", $gene_name_ens_id);
			$gene_name = $temp[0];
			$ens_id = $temp[1];
			$genecard_link = "=HYPERLINK(\'http://www.genecards.org/cgi-bin/carddisp.pl?id=$ens_id&id_type=ensembl\', \'GeneCard Link\')";
		} else {
			$gene_name = "NA";
			$ens_id = "NA";
			$genecard_link = "NA";
		}
		# add OMIM links and Uniprot links
		############# all links added need to have single quotes replaced bu double quotes after wAnnovar ############
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
		
		# find if the gene is in the interested gene list
		my $isInterested = "NO";
		if (exists $isInterestedGenes{$ens_id}) {
			$isInterested = "YES";
		}
		
		#find the corresponding omim annotation
		my $omim_anno = "NA\tNA\tNA\tNA";
		if (exists $OMIM{$ens_id}) {
			$omim_anno = $OMIM{$ens_id};
		}		
		
		#find in house MAF for variants
		my $maf = "";
		my $v = $variants[0];
		my @temp = split("\t", $v);
		my $chr1 = $temp[0];
		my $pos1 = $temp[1];
		my $ref1 = $temp[2];
		my $alt1 = $temp[3];
		if (exists $inHouse_MAF{$chr1}{$pos1}{$ref1}{$alt1} ) {
			$maf = $inHouse_MAF{$chr1}{$pos1}{$ref1}{$alt1};
		} else {
			$maf = "NA";
		}

		#find in house 'GATK' MAF for variants
		my $maf_GATK = "";
		#my $v = $variants[0];
		#my @temp = split("\t", $v);
		#my $chr1 = $temp[0];
		#my $pos1 = $temp[1];
		#my $ref1 = $temp[2];
		#my $alt1 = $temp[3];
		if (exists $inHouse_MAF_GATK{$chr1}{$pos1}{$ref1}{$alt1}) {
			$maf_GATK = $inHouse_MAF_GATK{$chr1}{$pos1}{$ref1}{$alt1};
		} else {		
			$maf_GATK = "NA";
		}
		

		## spliting the line with multi ALTs and asign MAF to each variant
		#if (scalar(@variants) > 1) {
		#	for (my $i=0;$i<scalar(@variants);$i++) {
		#		my $v = $variants[$i];
		#		my $single_maf;
		#		my $single_maf_GATK;
		#		my @temp = split("\t", $v);
		#		my $chr = $temp[0];
		#		my $pos = $temp[1];
		#		my $ref = $temp[2];
		#		my $alt = $temp[3];
		#		my $V_type = $temp[4];
		#		if (exists $inHouse_MAF{$chr}{$pos}{$ref}{$alt} ) {
		#			$single_maf = $inHouse_MAF{$chr}{$pos}{$ref}{$alt};
		#		} else {
		#			$single_maf = "NA";
		#		}
		#		
		#		if (exists $inHouse_MAF_GATK{$chr}{$pos}{$ref}{$alt} ) {
		#			$single_maf_GATK = $inHouse_MAF_GATK{$chr}{$pos}{$ref}{$alt};
		#		} else {
		#			$single_maf_GATK = "NA";
		#		}
		#		#print "$i: $alt $A\n";
		#		my $alter_temp = $ALTs[$i];
		#		for (my $j=0;$j<scalar @ALTs;$j++) {
		#			if ($ALTs[$i] ne $ALTs[$j]) {
		#				$alter_temp = $alter_temp.",".$ALTs[$j];
		#			}
		#		}
		#		my $Sample_Call_processed = AddGenoTypeToSampleCalls_CompondHet($ALTs[$i], $A, $FORMAT, $Sample_Call);
		#		push @output, $line."\t".$Sample_Call_processed."\t".$SNPorINDEL ."\t".$gene_name."\t".$ens_id."\t".$isInterested."\t"
		#		.$single_maf."\t".$single_maf_GATK."\t".$omim_anno."\t".$genecard_link."\t".$omim_link."\t".$uniprot_link."\n";
		#	}
		#} else {
			my $Sample_Call_processed = AddGenoTypeToSampleCalls($FORMAT, $Sample_Call);
			for (my $i=0;$i<=46;$i++){ push @output, $elements[$i]."\t";}
			push @output, $Sample_Call_processed."\t".$SNPorINDEL ."\t".$gene_name."\t".$ens_id."\t".$isInterested."\t"
			.$maf."\t".$maf_GATK."\t".$omim_anno."\t".$genecard_link."\t".$omim_link."\t".$uniprot_link."\n";
		#}
	}
}
close (VCF);

open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
print OUTPUT @output;
close OUTPUT;



################################################
#
#  output to excel part
#
#


if ($output_excel eq "") {print "#####  No excel file is required. ######\nDone!\n"; exit;}

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
#my %gene_name_count; #it's for selecting Xlinked rare variants



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
	if ($line !~ m/^Chr\t/) {
		my @elements = split("\t", $line);
		#if ($elements[26] ne "unknown") { # insert "unknown" to make every row in the same length
		#	splice @elements, 26, 0, 'unknown'; 
		#}
		
		if($elements[7] eq "synonymous SNV" && $elements[11] eq "synonymous SNV" && $elements[15] eq "synonymous SNV"  ) { # eliminate synonymous SNV
			next;
		}
		if($elements[5] !~ m/(^exonic|^splicing)/ &&  $elements[9] !~ m/(^exonic|^splicing)/ && $elements[13] !~ m/(^exonic|^splicing)/  ) { # eliminate intronic or intergenic variants
			next;
		}
		
		my $sample_string="";
		my $rr_na_sum=0;
		for (my $i=0;$i<scalar @sample_names;$i++) {
			if ($CNV_file ne ""){ #when CNV file is available
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
			} else {
				$sample_string=$sample_string.$elements[$sample_columns[$i]]."\t".$elements[$sample_columns[$i]+1]."\t";
			}
			
			if($elements[$sample_columns[$i]+1] eq 'R/R' || $elements[$sample_columns[$i]+1] eq 'NA') {
				$rr_na_sum+=1;
			}
		}
		
		my $length = scalar @elements -1;
		
		#Change 1000genome, ESP6500 and cg69 var freq empty value to "0"
		if($elements[20] eq "NA" ) {#1000G
			$elements[20] = 0;
		}
		if($elements[19] eq "NA" ) { #esp6500
			$elements[19] = 0;
		}
		if($elements[21] eq "NA" ) { #cg69
			$elements[21] = 0;
		}
		
		#Change in-house and GATK Freq 'NA' to '0'
		if ($elements[$length-8] eq "NA") { #inHouse maf
			$elements[$length-8] = 0;
		}
		if ($elements[$length-7] eq "NA") { #inhouse gatk maf
			$elements[$length-7] = 0;
		}
		
		my $string_1 = $elements[38]."\t".$elements[39]."\t".$elements[41]."\t".$elements[42]."\t"
		.$elements[$length-12]."\t".$elements[43]."\t".$elements[5]."\t".$elements[7]."\t".$elements[8]."\t".$elements[9]
		."\t".$elements[11]."\t".$elements[12]."\t".$elements[13]."\t".$elements[15]."\t".$elements[16]."\t".$elements[$length-11]."\t"
		.$elements[$length-10]."\t".$elements[$length-9]."\t".$elements[17]."\t".$elements[18]."\t".$elements[19]."\t".$elements[20]."\t".$elements[21]."\t".
		$elements[$length-8]."\t".$elements[$length-7]."\t".$elements[22]."\t".$elements[23]."\t".$elements[24]."\t".$elements[25]."\t".$elements[26]."\t".$elements[27]."\t"
		.$elements[28]."\t".$elements[29]."\t".$elements[30]."\t".$elements[31]."\t".$elements[32]."\t".$elements[33]."\t".$elements[34]."\t".$elements[35]."\t"
		.$elements[36]."\t".$elements[37]."\t".$elements[46]."\t";
		
		#re-arrange columns
		
		
		$elements[$length-2] =~ s/\'/\"/g;
		$elements[$length-1] =~ s/\'/\"/g;
		$elements[$length] =~ s/\'/\"/g;
		
		my $string_2 = $elements[44]."\t".$elements[45]."\t".$elements[$length-6]."\t".$elements[$length-5]."\t".$elements[$length-4]
		."\t".$elements[$length-3]."\t".$elements[$length-2]."\t".$elements[$length-1]."\t".$elements[$length]; #re-arrange columns
		
		if($rr_na_sum < scalar @sample_names) { # RR and NA are removed
			push @output_all, $string_1.$sample_string.$string_2;
			
			# for filtered spreadsheet

			if( $elements[19] <= 0.05 && $elements[20] <= 0.05 ) {
				push @output_filtered, $string_1.$sample_string.$string_2;
				
				#unless(exists $gene_name_count{$elements[$length-9]}) { # could be deleted
				#	$gene_name_count{$elements[$length-9]} = 0;
				#}
				
				# for deleterious hits spreadsheet
				my $number_deleterious_predictions=0;
				if($elements[25] eq "D" || $elements[25] eq "P") {
					$number_deleterious_predictions+=1;
				}
				if($elements[27] eq "D" || $elements[27] eq "P") {
					$number_deleterious_predictions+=1;
				}
				if($elements[29] eq "D") {
					$number_deleterious_predictions+=1;
				}
				if($elements[31] eq "A" || $elements[31] eq "D") {
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
					unless (exists $all_gene_id{$elements[$length-10]}) { $all_gene_id{$elements[$length-10]} = $elements[$length-11]."_".$elements[$length-9]; }
					if ( exists $Variants_on_samples{$sample_names[$i]}{$elements[$length-10]} ) {
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-10]} = 
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-10]}.";".
						$elements[38]."_".$elements[39]."_".$elements[$sample_col+1]."_".sprintf( "%.4f",$elements[19])."_".sprintf( "%.4f",$elements[20])."_".sprintf( "%.4f",$elements[$length-8]);
					} else {
						$Variants_on_samples{$sample_names[$i]}{$elements[$length-10]} = 
						$elements[38]."_".$elements[39]."_".$elements[$sample_col+1]."_".sprintf( "%.4f",$elements[19])."_".sprintf( "%.4f",$elements[20])."_".sprintf( "%.4f",$elements[$length-8]);				
					}
				}
				
				#for rare variants on a same gene
				#if($number_CompoundHet < scalar @sample_names) {
				#	$gene_name_count{$elements[$length-11]} =$gene_name_count{$elements[$length-11]}+1;
				#} #retired
				
				if ($elements[$length-9] eq "YES" && $number_deleterious_predictions >=1) {  #those that were predicted as “Deleterious” by at least two predictors
					push @output_possibleHits, $string_1.$sample_string.$string_2;
				} elsif ($elements[$length-9] eq "NO" && $number_deleterious_predictions >=3) { #those that were predicted as “Deleterious” by all 4 predictors
					push @output_possibleHits, $string_1.$sample_string.$string_2;
				} elsif($elements[$length-12] eq "SNP" && $number_CompoundHet >=1) { #Compound heterozygous “V1/V2”* of single base polymorphisms
					push @output_possibleHits, $string_1.$sample_string.$string_2;
				} elsif($elements[$length-9] eq "YES" && $elements[19] <= 0.01 && $elements[20] <= 0.01) { # those that have AAF less than 0.01 in both 1000Genomes project and ESP6500 
					push @output_possibleHits, $string_1.$sample_string.$string_2;
				}
				
				if ($elements[$length-9] eq "YES" && $elements[19] <= 0.01 && $elements[20] <= 0.01) {
					push @very_rare_interested_vars, $string_1.$sample_string.$string_2;
				}
				
				#for Xlinked rare variants
				if ($elements[38] =~ m/chrx/i) {
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
					if ( $elements[38] !~ m/chr[XYMxym]/ ) {
						push @output_ARmodelHits, $string_1.$sample_string.$string_2;
					}
				}
			}
		}
	}
}
close(INPUT);


# record CNVs into the hash %CNV_on_samples
if ($CNV_file ne "") {
	####
	# need to read in the whole list of insterested genes in ENSEMBL id
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
									if ( ! exists $isInterestedGenes{$gene_id}) {
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
my @header = 
('Chr','Start','Ref','Obs','SNPorINDEL','VariantCall_quality','Func.ensGene','ExonicFunc.ensGene','AAChange.ensGene','Func.knownGene','ExonicFunc.knownGene',
'AAChange.knownGene','Func.refGene','ExonicFunc.refGene','AAChange.refGene','GeneName','GeneID','isInterested','phastConsElements46way','genomicSuperDups',
'ESP6500_ALL','1000g2012feb_ALL','cg69','InHouseMAF','InHouseMAF_GATK','dbSNP137','LJB2_SIFT','LJB2_PolyPhen2_HDIV','LJB2_PP2_HDIV_Pred','LJB2_PolyPhen2_HVAR',
'LJB2_PolyPhen2_HVAR_Pred','LJB2_LRT','LJB2_LRT_Pred','LJB2_MutationTaster','LJB2_MutationTaster_Pred','LJB_MutationAssessor','LJB_MutationAssessor_Pred',
'LJB2_FATHMM','LJB2_GERP++','LJB2_PhyloP','LJB2_SiPhy','Format');

foreach my $sample (@sample_names) {
	push @header, $sample;
	push @header, $sample.".anno";
	if ($CNV_file ne "") {
		push @header, $sample.".cnv";
	}
}
my @header_2 = ('Filter', 'INFO', 'GoTerm', 'WikiGene_Description', 'MIM_Gene_Description', 'OMIM_Gene_Description', 
'GeneCard Link', 'OMIM Link', 'Uniprot Link');
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
				######### testing here  ##################
				#foreach my $wolegequ (@temp_temp){
				#	print $wolegequ."  ";
				#}
				#print "\n";
				##########################################
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
    print "\nusage: VCF_2_annotated_excel_20131120.pl \n";
    print "--vcf input vcf file (of a single sample or a family;\n";
    print "--InterestedGenes file of interested gene names list (optional); Format: Ensembl gene IDs on the 1st column.\n"; 
    print "--out output excel file name (optional, only when required);\n";
    print "--CNV CNV result file, output of HG's Annotate_CNVs_combine_multiple_files.pl, sample names must be consistent with the vcf.\n\n";
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

sub GetInHouseMaf {
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

sub GetIsInterestedGenes { # one ensembl gene_id on each line
	my ($InterestedGenefile) = @_;
	my %Mito=();
	if (-e $InterestedGenefile) {
		open MF, $InterestedGenefile or die "cannot open $InterestedGenefile";
		while (my $line = <MF>) {
			chomp $line;
			my @genam=split(/\t/,$line);
			if(!exists $Mito{$genam[0]}){
				$Mito{$genam[0]}=0;
			}	
		}
		close MF;
	} else {
		print "No insterested gene list is provided.\n";
	}
	return \%Mito;
}

sub AddGenoTypeToSampleCalls {
	my ($format, $sample_call) = @_;
	my $gene_type_call_qual = 13; ##### genotype call quality cut off
	my @format_fields = split(":", $format);
	my $GT_index;
	my $GQ_index;
	for (my $i=0;$i<scalar(@format_fields);$i++){
		if ($format_fields[$i] =~ m/GT/) {
			$GT_index = $i;
		}
		if ($format_fields[$i] =~ m/GQ/) {
			$GQ_index = $i;
		}
	}
	my $sample_call_processed=""; # vcf_sample1 \t geno_sample1 \t vcf_sample2 \t geno_sample2 ....
	my @sample_call_split = split("\t", $sample_call);
	for(my $i=0;$i<scalar(@sample_call_split);$i++){
		my $sample = $sample_call_split[$i];
		if ($sample !~ m/\.\/\./) {
			my @fields = split(":", $sample);
			my $GT = $fields[$GT_index];
			my $GQ = $fields[$GQ_index];
			if ($GQ < $gene_type_call_qual) {
				#$GT = "R/R";
				$GT = "NA"; #if low quality make a null call "NA"
			} elsif ($GT =~ m/0\/0/) {
				$GT = "R/R";
			} elsif ($GT =~ m/0\/[123456789]/) {
				$GT = "R/V";
			} elsif ($GT =~ m/([123456789])\/([123456789])/) {
				my $left = $1;
				my $right = $2;
				if($left!=$right) {
					$GT = "V1/V2"; 
				} else {
					$GT = "V/V";
				}
			}
			$sample_call_processed = $sample_call_processed."\t".$sample."\t".$GT;
		} else {
			$sample_call_processed = $sample_call_processed."\t".$sample."\t"."NA";
		}
	}
	$sample_call_processed =~ s/^\t//;
	return $sample_call_processed;
}

sub AddGenoTypeToSampleCalls_CompondHet {
	my ($first_V, $variants, $format, $sample_call) = @_;
	my $gene_type_call_qual = 13; ##### genotype call quality cut off
	my @format_fields = split(":", $format);
	my $GT_index;
	my $GQ_index;
	my $AD_index;
	for (my $i=0;$i<scalar @format_fields ;$i++){
		if ($format_fields[$i] =~ m/GT/) {
			$GT_index = $i;
		}
		if ($format_fields[$i] =~ m/GQ/) {
			$GQ_index = $i;
		}
		if ($format_fields[$i] =~ m/AD/) {
			$AD_index = $i;
		}		
	}
	my $sample_call_processed=""; # vcf_sample1 \t geno_sample1 \t vcf_sample2 \t geno_sample2 ....
	my @sample_call_split = split("\t", $sample_call);
	for(my $i=0;$i<scalar(@sample_call_split);$i++){
		my $sample = $sample_call_split[$i];
		if ($sample !~ m/\.\/\./) {
			my @fields = split(":", $sample);
			my $GT = $fields[$GT_index];
			my $GQ = $fields[$GQ_index];
			my $AD;
			if(defined $AD_index) {
				$AD = $fields[$AD_index];
			}
			
			if( $first_V eq $variants) {
				if ($GQ < $gene_type_call_qual) {
					#$GT = "R/R";
					$GT = "NA"; #low qual genotype ... make a null call "NA"
				} elsif ($GT =~ m/0\/0/) {
					$GT = "R/R";
				} elsif ($GT =~ m/0\/[123456789]/) {
					$GT = "R/V";
				} elsif ($GT =~ m/([123456789])\/([123456789])/) {
					my $left = $1;
					my $right = $2;
					if($left!=$right) {
					$GT = "V1/V2"; 
					} else {
					$GT = "V/V";
					}
				}
				$sample_call_processed = $sample_call_processed."\t".$sample."\t".$GT;
			} else {
				#determine which alternative allele is moved to the front
				my @V = split("," , $variants);
				my @AD_array;
				if(defined $AD_index) {
					@AD_array = split(",", $AD);
				}
				my $V_index;
				my $the_first_V_AD_index; 
				for(my $h=0;$h<scalar @V;$h++) {
					if($first_V eq $V[$h]) {
						$V_index = $h+1;
						$the_first_V_AD_index = $V_index+1;
					}
				}
				
				my %temp_hash;
				$temp_hash{$V_index} = 1;
				my $new_AD;
				if(defined $AD_index) {
					#print "sample_".$sample."\n";
					#print "AD_ARRAY_".$AD_array[0]."\n";
					#print "the first index _".$the_first_V_AD_index."\n";
					#print "ad array [index]_".$AD_array[$the_first_V_AD_index-1]."\n";
					if(scalar @AD_array == 1 && $AD_array[0] =~ m/\./ ) {
						$new_AD = '.';
					} else {
						$new_AD = $AD_array[0].",".$AD_array[$the_first_V_AD_index-1];
						for(my $h=1;$h < scalar @AD_array;$h++) {
							if($h != $the_first_V_AD_index-1) {
								$new_AD = $new_AD.",".$AD_array[$h];
							}
						}
					}
					
				}
				
				my $j = 2;
				for(my $h=1;$h<=scalar @V;$h++) {
					if ($h != $V_index) {
						$temp_hash{$h} = $j;
						$j+=1;
					}
				}
				
				#adding genotye
				if ($GQ < $gene_type_call_qual) {
					#$GT = "R/R";
					$GT = "NA"; #low qual null call!
				} elsif ($GT =~ m/0\/0/) {
					$GT = "R/R";
				} elsif ($GT =~ m/0\/[123456789]/) {
					$GT = "R/V";
				} elsif ($GT =~ m/([123456789])\/([123456789])/) {
					my $left = $1;
					my $right = $2;
					if($left!=$right) {
						$GT = "V1/V2"; 
					} else {
						$GT = "V/V";
					}
				}
				
				my $left_number;
				my $right_number;
				if ($fields[$GT_index] =~ m/(\d)\/(\d)/ ) {
					my $A_left = $1;
					my $A_right = $2;
					#print "left $A_left\n";
					#print "right $A_right\n";
					if ($A_left == $A_right && $A_left ==0 ) {
						$left_number = "0";
						$right_number = "0";
					} elsif ( $A_left ==0 ) {
						$left_number = 0;
						$right_number = $temp_hash{$A_right};
					} else {
						$left_number = $temp_hash{$A_left};
						$right_number = $temp_hash{$A_right};
					}
				}
				###replace sample GT numbers
				$sample = $left_number."/".$right_number;
				for(my $h=1;$h < scalar @fields;$h++) {
					if(defined $AD_index) {
						if($h == $AD_index) {
							$sample=$sample.":".$new_AD;
						} else {
							$sample=$sample.":".$fields[$h];
						}
					} else {
						$sample=$sample.":".$fields[$h];
					}
					my $teno = scalar @fields;
					#print "h is $h   total is $teno   fields[h] is $fields[$h]\n";
					#print "$sample\n";
				}
				
				$sample_call_processed = $sample_call_processed."\t".$sample."\t".$GT;
			}
		} else {
			$sample_call_processed = $sample_call_processed."\t".$sample."\t"."NA";
		}
	}
	$sample_call_processed =~ s/^\t//;
	return $sample_call_processed;
}




