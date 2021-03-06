#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

###############################################
#
#  amended from Helen's
###############################################


my $vcf;
my $InterestedGenefile="/users/a5907529/lustre/Kidney_20130211/john_gene_list.txt";


my $ref="/users/a5907529/lustre/Yaobo/GenomeData/GATK_bundle/ucsc.hg19.4GATK.fasta";
my $genefile="/users/a5907529/lustre/scripts/AfterDindel/Ensembl_Genes_R67.txt";
my $sharedfile="/users/a5907529/lustre/scripts/AfterDindel/InHouse_OnTarget_Variants_MAFs.txt";
my $OMIMfile="/users/a5907529/lustre/scripts/AfterDindel/Ensembl_OMIM_AllGeneNames.txt";

my $Results=GetOptions("vcf=s"=>\$vcf, "InterestedGenes=s"=>\$InterestedGenefile);
my $output_file = $vcf.".inHouse_annotated.txt";

### Read in all annotations first ###
my %en_genes;
my %inHouse_MAF;
my %OMIM;
my %isInterestedGenes;

my $en_genes_ref = GetGeneCoords($genefile);
%en_genes = %$en_genes_ref;
my $inHouse_MAF_ref = GetInHouseMaf($sharedfile);
%inHouse_MAF = %$inHouse_MAF_ref;
my $OMIM_ref = GetOMIManno($OMIMfile);
%OMIM = %$OMIM_ref;
if ($InterestedGenefile ne "") {
	my $isInterestedGenes_ref = GetIsInterestedGenes($InterestedGenefile);
	%isInterestedGenes = %$isInterestedGenes_ref;
}


my @output;
open VCF, "$vcf" or die "Can not open the file $vcf";
while (my $line = <VCF> ) {
	chomp $line;
	if($line =~ m/^\#/) {
		push @output, $line."\n";
		next;
	} else {
		my @elements  = split ("\t", $line);
		
		my $chr=$elements[0];
		my $pos=$elements[1];
		my $R=$elements[3];
		my $A=$elements[4];
		my @variants; # chr \t pos \t ref \t v
		my @ALTs;
		
		my $FORMAT=$elements[8];
		my $Sample_Call="";
		for(my $i=9;$i<scalar(@elements);$i++) {
			$Sample_Call=$Sample_Call."\t".$elements[$i];
		}
		$Sample_Call =~ s/^\t//;
		my $Sample_Call_processed = AddGenoTypeToSampleCalls($FORMAT, $Sample_Call);
		
		#print "$chr\t$pos\t$R\t$A\n";
		
		if ($R =~ m/^[atgcATGC]$/ && $A =~ m/^[atgcATGC]$/) { # confirm it's a SNP not INDEL
			push @variants, "$chr\t$pos\t$R\t$A";
			push @ALTs, $A;
		} elsif ($R =~ m/^[atgcATGC]$/) { # insertions
			if($A !~ /\,/){
				my $v = $A;
				if ( length $v > 1 ) {
					$v =~ s/^./\+/;
				}
				push @variants, "$chr\t$pos\t$R\t$v";
				push @ALTs, $A;
			} else { # multi ALT bases insertions
				my @vars=split(/\,/,$A);
				foreach my $v (@vars) {
					push @ALTs, $v;
					if ( length $v > 1 ) {
						$v =~ s/^./\+/;
					}
					push @variants, "$chr\t$pos\t$R\t$v";
				}
			}
		} else { ##deletion(s) !!and in some cases of multiple (,) vars insertions!! 
			if ( $A !~ /\,/) {
					my $v = $R;
					if ( length $v > 1 ) {
						$v =~ s/^./\-/;
					}
					push @variants, "$chr\t$pos\t$R\t$v";
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
						push @variants, "$chr\t$pos\t$ref\t$alt";
					} elsif(length $v > length $R){
						my $diff=(length $R)-1; 
						my $alt = substr($v,$diff);
						my $ref = substr($v,$diff,1); 
						$alt =~s/^./\+/;
						$pos = $pos+$diff;
						push @ALTs, $v;
						push @variants, "$chr\t$pos\t$ref\t$alt";					
					} else {
						print "ERROR: Dangerous line: \n".$line."\n Try to remove the Line before run annotation.\n";
						exit;
					}
				}
			}
		}
		
		my $vcf_pos=$elements[1];
		my $gene_name_ens_id = "";
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
		} else {
			$gene_name = "NA";
			$ens_id = "NA";
		}
		
		my $isInterested = "NO";
		if (exists $isInterestedGenes{$gene_name}) {
			$isInterested = "YES";
		}
		
		my $maf = "";
		foreach my $v (@variants) {
			my @temp = split("\t", $v);
			my $chr = $temp[0];
			my $pos = $temp[1];
			my $ref = $temp[2];
			my $alt = $temp[3];
			if (exists $inHouse_MAF{$chr}{$pos}{$ref}{$alt} ) {
				$maf = $maf.",".$inHouse_MAF{$chr}{$pos}{$ref}{$alt};
			} else {
				$maf = $maf.","."NA";
			}
		}
		$maf =~ s/^\,//;
		
		my $omim_anno = "NA\tNA\tNA";
		if (exists $OMIM{$ens_id}) {
			$omim_anno = $OMIM{$ens_id};
		}
		
		push @output, $elements[0]."\t".$elements[1]."\t".$elements[2]."\t".$elements[3]."\t".$elements[4]."\t".$elements[5]."\t".$elements[6]."\t".$elements[7]."\t".$elements[8]."\t".$Sample_Call_processed."\t".$gene_name."\t".$ens_id."\t".$isInterested."\t".$maf."\t".$omim_anno."\n";

		
		## spliting the line with multi ALTs and asign MAF to each variant
		for (my $i=0;$i<scalar(@variants);$i++) {
			my $v = $variants[$i];
			my $single_maf;
			my @temp = split("\t", $v);
			my $chr = $temp[0];
			my $pos = $temp[1];
			my $ref = $temp[2];
			my $alt = $temp[3];
			if (exists $inHouse_MAF{$chr}{$pos}{$ref}{$alt} ) {
				$single_maf = $inHouse_MAF{$chr}{$pos}{$ref}{$alt};
			} else {
				$single_maf = "NA";
			}
			push @output, $elements[0]."\t".$elements[1]."\t".$elements[2]."\t".$elements[3]."\t".$ALTs[$i]."\t".$elements[5]."\t".$elements[6]."\t".$elements[7]."\t".$elements[8]."\t".$Sample_Call_processed."\t".$gene_name."\t".$ens_id."\t".$isInterested."\t".$single_maf."\t".$omim_anno."\n";
		}
		#### right here here here
		#### to be done
	}
}
close (VCF);

open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
print OUTPUT @output;
close OUTPUT;
exit;

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

sub GetOMIManno { #undone
	my ($OMIMfile) = @_;
	my %OMIM;
	open INPUT, $OMIMfile or die "Cannot open $OMIMfile\n"; 
	while (my $Line = <INPUT>){
		chomp $Line;
		my @linesplit = split(/\t/,$Line);
		my $ens_id = $linesplit[0];
		my $go;
		my $wiki;
		my $MIM;
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
		
		$OMIM{$ens_id} = $go."\t".$wiki."\t".$MIM;
	}
	close INPUT;
	return \%OMIM;
}

sub GetIsInterestedGenes { # one gene name on each line
	my ($InterestedGenefile) = @_;
	my %Mito=();
	if (-e $InterestedGenefile) {
		open MF, $InterestedGenefile or die "cannot open $InterestedGenefile";
		while (my $line = <MF>) {
			chomp $line;
			if(!exists $Mito{$line}){
				$Mito{$line}=0;
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
	my $gene_type_call_qual = 13;
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
				$GT = "R/R";
			} elsif ($GT =~ m/0\/0/) {
				$GT = "R/R";
			} elsif ($GT =~ m/0\/[123456789]/) {
				$GT = "R/V";
			} elsif ($GT =~ m/[123456789]\/[123456789]/) {
				$GT = "V/V";
			}
			$sample_call_processed = $sample_call_processed."\t".$sample."\t".$GT;
		} else {
			$sample_call_processed = $sample_call_processed."\t".$sample."\t"."NA";
		}
	}
	$sample_call_processed =~ s/^\t//;
	return $sample_call_processed;
}