#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

############ current status : Aborted #############

###############################################
#
#  amended from Helen's
# 1. use abs file path
# 2. read in vcf files. it's to avoid the possible mistakes in recording indels positions
# 3. use mpileup
# 4. read in VCF result of mpileup files, it to avoid dealing with qulity scores and gap alignment in pileup files.
#
###############################################

	      #
	   #	
	###############could change it to use GATK reclibrated alignment files
	  #
	    #	

my $projectID="Any_ID";
my $samplePath="/users/a5907529/lustre/scripts/AfterDindel/Yaobo_pipeline/";
my $sampleNames="test";
my $ref="/users/a5907529/lustre/Yaobo/GenomeData/GATK_bundle/ucsc.hg19.4GATK.fasta";
my $Mitofile="";
#$Mitofile is a file contain of genes of insterests. set to an empty file Genes_empty.txt if it's not avail
#my $Mitofile="/users/a5907529/lustre/scripts/AfterDindel/Yaobo_pipeline/DLB_genes.txt";
my $genefile="/users/a5907529/lustre/scripts/AfterDindel/Yaobo_pipeline/Ensembl_Genes_R67.txt";
#my $sharedfile="/users/a5907529/lustre/scripts/AfterDindel/Yaobo_pipeline/InHouse_OnTarget_Variants_MAFs.txt_noDLB";
my $sharedfile="/users/a5907529/lustre/scripts/AfterDindel/Yaobo_pipeline/InHouse_OnTarget_Variants_MAFs.txt";
my $cov_threshold=5; #min cov in pileup for variant to be called, otherwise '0' call

my $Results=GetOptions("samplePath=s"=>\$samplePath, "sampleNames=s"=>\$sampleNames, "covThreshold=i"=>\$cov_threshold, "ref=s"=>\$ref, "candGenes=s"=>\$Mitofile, "ensemblGenes=s"=>\$genefile, "controls=s"=>\$sharedfile, "projectID=s"=>\$projectID);

# Sample names need to split with commars

#Program to produce list of all 'on-target' variants seen in patients (and unaffecteds)
#Pseudo-call 'missing' variants
#Attach Ensembl gene name and in-house MAFs to variants
#output list to file in Annovar input format

#define variables
my @SampleNames = split(",", $sampleNames);
my %Variants;
my %genes;
$samplePath =~ s/\/$//;


#SNV file
foreach my $sample (@SampleNames) {
	my $file = $samplePath."/".$sample."/"."varscan"."/".$sample."what ever you call this output";
	if (-e $file) {
		open INPUT, $file or die "Cannot open ".$file."\n";
		my $head=<INPUT>;
		while (my $Line = <INPUT>){
			chomp $Line;
			my @linesplit1 = split (/\t/,$Line);
			my $chr = $linesplit1[1];
			my $pos_1 = $linesplit1[3]+1;
			my $ref = $linesplit1[4];
			my $var = $linesplit1[5];
			my $geno = 0; #not covered or ref
			
			if($linesplit1[8]<85){
				$geno="R/V";
			} #het
			
			if($linesplit1[8]>=85){
				$geno="V/V";
			} #hom
			
			if(!exists $Variants{$chr}{$pos_1}{$ref}{$var}){
				foreach my $sample_sub (@SampleNames) {
					$Variants{$chr}{$pos_1}{$ref}{$var}{$sample_sub}=0; ##add in 'all', 'het-only', 'hom-only counts'
				}
			}
			
			if(exists $Variants{$chr}{$pos_1}{$ref}{$var}){
				$Variants{$chr}{$pos_1}{$ref}{$var}{$sample}=$geno;
			}
			
			if(!exists $genes{$chr}{$pos_1}){
				$genes{$chr}{$pos_1}[0]='NO_ENSEMBL_GENE_MATCH';
			}		
		}
		print "$file";
		close INPUT;
	} else {
		print "No SNP file can find for sample $sample.\n";
	}
}
	
	
#Indel file
foreach my $sample (@SampleNames) {
	my $file = $samplePath."/".$sample."/"."dindel"."/".$sample."_nodups_FilteredIndels_OnTarget_vcf4_all";
	if (-e $file) {
		open INPUT, $file or die "Cannot open ".$file."\n";
		while (my $Line = <INPUT>){
			chomp $Line;
			my @linesplit1 = split (/\t/,$Line);
			my $chr = $linesplit1[0];
			my $pos_1 = $linesplit1[1];
			my $ref = $linesplit1[2];
			my $var = $linesplit1[3];
			my $geno = 0; #not covered or ref
			if($Line!~/1\/1/){$geno="R/V";} #het allows for 1/1, 1/2, 1/3 etc!!
			if($Line=~/1\/1/){$geno="V/V";} #hom
			
			if(!exists $Variants{$chr}{$pos_1}{$ref}{$var}){
				foreach my $sample_sub (@SampleNames) {
					$Variants{$chr}{$pos_1}{$ref}{$var}{$sample_sub}=0; ##add in 'all', 'het-only', 'hom-only counts'
				}
			}
			
			if(exists $Variants{$chr}{$pos_1}{$ref}{$var}){
					$Variants{$chr}{$pos_1}{$ref}{$var}{$sample}=$geno;
			}
			
			if(!exists $genes{$chr}{$pos_1}){
				$genes{$chr}{$pos_1}[0]='NO_ENSEMBL_GENE_MATCH';
			}
		}
		print "$file\n";
		close INPUT;
	} else {
		print "No Indel file can find for sample $sample.\n";
	}
}


##generate chr pos1 file to then generate smaller variant specific pileup files for each patient
my $out_ch_pos = $samplePath."/"."chr_pos_for_pileup.txt";
open OUTP, ">$out_ch_pos"  or die "Cannot open file \"$out_ch_pos\" to write to!\n";
foreach my $chrm (sort keys %Variants){
	foreach my $pos (sort {$a<=>$b} keys %{$Variants{$chrm}}){
		print OUTP "$chrm\t$pos\n";
	}
}
close OUTP;

floop: foreach my $sample (@SampleNames){
	my $bam_file = $samplePath."/".$sample."/"."picard"."/".$sample."_nodups.sorted.bam";
	      #
	   #	
	###############could change it to use GATK reclibrated alignment files
	  #
	    #	
	my $pexit = $samplePath."/".$sample."/"."picard"."/".$sample.".pileup_varpos";
	if (-e $pexit) {
		next floop;
	} else {
		print "start samtools\n";
		`samtools mpileup -f $ref -l $out_ch_pos $bam_file >  $pexit`;
		#`samtools mpileup -f $ref -l $out_ch_pos $File > $Id.pileup_varpos`
	}
}


#Add pileup file info to distinguish between Reference alleles (3) and Non-covered variants (0) - pseudo-call (0) variants as hom-ref/het/hom-var from pileup

foreach my $sample (@SampleNames){
	my $file = $samplePath."/".$sample."/"."picard"."/".$sample.".pileup_varpos";
	open INPUT, $file or die "Cannot open $file.\n";
	pileup_loop: while (my $Line = <INPUT>){
		chomp $Line;
		my @linesplit1 = split (/\t/,$Line);
		if(!exists $Variants{$linesplit1[0]}{$linesplit1[1]}){
			next pileup_loop;
		}
		my $chr = $linesplit1[0];
		my $pos_1 = $linesplit1[1];
		my $ref = uc($linesplit1[2]);
		my $cov = $linesplit1[3];
		if($cov<$cov_threshold){
			foreach my $v (keys %{$Variants{$chr}{$pos_1}{$ref}}){
				$Variants{$chr}{$pos_1}{$ref}{$v}{$sample}="0"; #make any variants called below the coverage threshold = '0'
				next pileup_loop;
			}
		}
			#split ref/var calls and count - pseudocall ref/het/var
		my @basesplit=split(//,$linesplit1[4]);
		my $bc=0;
		my $a=0;
		my $indel='no';
		foreach my $b (@basesplit){
			$bc++;
			if($b=~/[ACGTacgt]/){$a++;}
			if($b=~/\+/){$indel='ins';}
			if($b=~/\-/){$indel='del';}
		}
		my $vf=0;
		$vf=$a/$bc;
		foreach my $v (keys %{$Variants{$chr}{$pos_1}{$ref}}){ 
			if($v !~/[\+\-]/ and $Variants{$chr}{$pos_1}{$ref}{$v}{$sample}=~/0/ and $vf<=0.25){
				$Variants{$chr}{$pos_1}{$ref}{$v}{$sample}="R/RP";
			} #homozygous ref and not 0=not-covered!!
			#if($v !~/[\+\-]/ and $Variants{$chr}{$pos_1}{$ref}{$v}{$ID}=~/0/ and $vf>0.25 and $vf<0.85){$Variants{$chr}{$pos_1}{$ref}{$v}{$ID}="R/VP";} #call heterozygous based on pileup bases - not called by Varscan/Dindel!!
			#if($v !~/[\+\-]/ and $Variants{$chr}{$pos_1}{$ref}{$v}{$ID}=~/0/ and $vf>=0.85){$Variants{$chr}{$pos_1}{$ref}{$v}{$ID}="V/VP";} #call homozygous variant based on pileup bases - not called by Varscan/Dindel!!
			if($v =~/[\+\-]/ and $indel =~/[id]/ and $Variants{$chr}{$pos_1}{$ref}{$v}{$sample}=~/0/ and $vf<=0.25){
				$Variants{$chr}{$pos_1}{$ref}{$v}{$sample}="R/RP";
			} #homozygous ref and not 0=not-covered!!
			#if($v =~/[\+\-]/ and $indel =~/[id]/ and $Variants{$chr}{$pos_1}{$ref}{$v}{$ID}=~/0/ and $vf>0.25 and $vf<0.85){$Variants{$chr}{$pos_1}{$ref}{$v}{$ID}="R/VP";} #call heterozygous based on pileup bases - not called by Varscan/Dindel!!
			#if($v =~/[\+\-]/ and $indel =~/[id]/ and $Variants{$chr}{$pos_1}{$ref}{$v}{$ID}=~/0/ and $vf>=0.85){$Variants{$chr}{$pos_1}{$ref}{$v}{$ID}="V/VP";} #call homozygous variant based on pileup bases - not called by Varscan/Dindel!!
		}
	}
	print "$file";
	close INPUT;
}


#get per chromosome hash of numerically sorted positions in an array
my %genes_pos_perchr;
my %genes_sortedpos;
my %startfrom;
foreach my $c (keys %genes){
	if(!exists $genes_pos_perchr{$c}){
		$genes_pos_perchr{$c}=0;
	}
	if(!exists $startfrom{$c}){
		$startfrom{$c}=0;
	}
	foreach my $p (sort {$a<=>$b} keys %{$genes{$c}}){
		#if(!exists $genes_pos_perchr{$c}){$genes_pos_perchr{$c}=0;}
		$genes_sortedpos{$c}[$genes_pos_perchr{$c}]=$p;
		$genes_pos_perchr{$c}++;
	}
}
	
#Known Mito, PID, etc. gene?
my %Mito=();
if (-e $Mitofile) {
	open MF, $Mitofile or die "cannot open $Mitofile";
	while (my $fLine = <MF>) {
		chomp $fLine;
		my @Mline = split (/\t/, $fLine);
		if(!exists $Mito{$Mline[3]}){
			$Mito{$Mline[3]}=0;
		}	
	}
	close MF;
} else {
	print "No insterested gene list is provided.\n";
}


#assign genename to each variant in %variants$samplePath
open INPUT2, $genefile or die "Cannot open $genefile\n";
geneloop: while (my $Line = <INPUT2>){
	chomp $Line;
	my @linesplit1 = split(/\t/,$Line);
	if($linesplit1[0] eq 'MT'){$linesplit1[0]='M'}
	my $chr="chr".$linesplit1[0];
	my $st=$linesplit1[1];
	my $end=$linesplit1[2];
	my $gen=$linesplit1[3];
	my $MitoGene="No";
	if(exists $Mito{$gen}){$MitoGene="Yes";}
	for(my $ct=0; $ct<$genes_pos_perchr{$chr}; $ct++){
	#for(my $ct=$startfrom{$chr}; $ct<$genes_pos_perchr{$chr}; $ct++){
		if($genes_sortedpos{$chr}[$ct]>$end){
			$startfrom{$chr}=$ct; next geneloop;
		}
		if($genes_sortedpos{$chr}[$ct]>=$st and $genes_sortedpos{$chr}[$ct]<=$end){
			if($genes{$chr}{$genes_sortedpos{$chr}[$ct]}[0]!~/NO_ENSEMBL_GENE/){
				$genes{$chr}{$genes_sortedpos{$chr}[$ct]}[0]=$genes{$chr}{$genes_sortedpos{$chr}[$ct]}[0].",".$gen; $startfrom{$chr}=$ct;
			}
			if($genes{$chr}{$genes_sortedpos{$chr}[$ct]}[0]=~/NO_ENSEMBL_GENE/){
				$genes{$chr}{$genes_sortedpos{$chr}[$ct]}[0]=$gen; $startfrom{$chr}=$ct;
			}
			if(!exists $genes{$chr}{$genes_sortedpos{$chr}[$ct]}[1]){
				$genes{$chr}{$genes_sortedpos{$chr}[$ct]}[1]=$MitoGene;
			}
			if($genes{$chr}{$genes_sortedpos{$chr}[$ct]}[1] eq "No" and $MitoGene eq "Yes"){
				$genes{$chr}{$genes_sortedpos{$chr}[$ct]}[1]=$MitoGene;
			}
		}
	}
}
close INPUT2;

#Get in-house control MAFs
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

#open output file and print list of variants in Annovar format
#my $outfile2 = $path1."Annovar_In".$fileID."_onTarget_variant_list_inhouseMAF.txt";
my $outfile2 = $samplePath."/"."Annovar_In_".$projectID."_CovThresh_".$cov_threshold."_onTarget_variant_list_inhouseMAF.txt";
open(OUT2, ">$outfile2") || die "Cannot open file \"$outfile2\" to write to!\n";		
foreach my $chrm (sort keys %Variants){
	foreach my $pos (sort {$a<=>$b} keys %{$Variants{$chrm}}){
		foreach my $r (keys %{$Variants{$chrm}{$pos}}){
			foreach my $v (keys %{$Variants{$chrm}{$pos}{$r}}){
				my $mf=0;
				if(exists $IHcontrols{$chrm}{$pos}{$r}{$v}){
					$mf=$IHcontrols{$chrm}{$pos}{$r}{$v};
				}
				if(!defined $genes{$chrm}{$pos}[1]){
					$genes{$chrm}{$pos}[1]="No";
				}
				#convert hg19 mtDNA pos to rCRS?
				#if($chrm eq "chrM" and $pos>=3107){$pos-=1;}
				#if($chrm eq "chrM" and $pos<3107){$pos-=2;}
				##convert chr, pos, ref, var to Annovar format and add to front of current output file line
				my $ch = $chrm;
				$ch =~ s/chr//;
				my $pos_1 = $pos;
				my $pos_2 = $pos_1;
				my $ref = $r;
				my $var = $v;
				if($var=~/^\+$/ or $var=~/^\+\>/ or $var=~/^\+[DEL]/){
					next;
				} #don't put indels with errors into output file
				#insertion
				if($var=~/\+/){
					$ref="-";$var=~s/\+//;
				}
				#deletion
				if($var=~/\-/){
					$var=~s/\-//;$ref=$var;$var="-";my $length=length($ref);$pos_2=$pos_1+$length;$pos_1++;
				}	
				############################################################################################
				print OUT2 "$ch\t$pos_1\t$pos_2\t$ref\t$var\t$chrm\t$pos\t$r\t$v\t$genes{$chrm}{$pos}[0]\t$genes{$chrm}{$pos}[1]\t$mf";
				my $total_alleles=0;
				my $minor_alleles=0;
				foreach my $pid (sort keys %{$Variants{$chrm}{$pos}{$r}{$v}}){
				print OUT2 "\t$Variants{$chrm}{$pos}{$r}{$v}{$pid}";
				if($Variants{$chrm}{$pos}{$r}{$v}{$pid}=~/R\/RP/){
					$total_alleles+=2;
				}
				if($Variants{$chrm}{$pos}{$r}{$v}{$pid}=~/R\/V/){
					$total_alleles+=2;$minor_alleles+=1;
				}
				if($Variants{$chrm}{$pos}{$r}{$v}{$pid}=~/V\/V/){$total_alleles+=2;$minor_alleles+=2;}
				}
				my $case_maf=2;
				if($total_alleles>0){
					$case_maf=$minor_alleles/$total_alleles;
				}
				print OUT2 "\t$case_maf\t$total_alleles\n";
			}
		}
	}
}
close OUT2;
exit;