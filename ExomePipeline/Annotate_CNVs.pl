#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

#Program to annotate ExomeDepth CNV prediction output files

my $inputfile = "/users/nhrg/lustre/AROS/CNV_results/InHouse_CNVs_with_Counts.txt";
my $CurrentDir="/users/nhrg/lustre/AROS/April2013_PFC_Samples_Results/CNVs";
my $fullname_file="/users/nhrg/lustre/hg19/Ensembl_OMIM_AllGeneNames_withChrPos.txt";
chomp $CurrentDir;

my $Results=GetOptions("inpath=s"=>\$CurrentDir, "controls=s"=>\$inputfile, "omim=s"=>\$fullname_file);

my @DirContent=`ls $CurrentDir`;

my %Variants=();
open INPUT, $inputfile or die "Cannot open $inputfile\n";
fileline: foreach my $Line (<INPUT>){
	chomp $Line;
	my @linesplit1 = split (/\t/,$Line);
	my $chr = $linesplit1[0];
	my $exon = $linesplit1[1];
	$Variants{$chr}{$exon}{'all'}=$linesplit1[2]; 
	$Variants{$chr}{$exon}{'dup'}=$linesplit1[3]; 
	$Variants{$chr}{$exon}{'del'}=$linesplit1[4];
	$Variants{$chr}{$exon}{'allBF20'}=$linesplit1[5]; 
	$Variants{$chr}{$exon}{'dupBF20'}=$linesplit1[6]; 
	$Variants{$chr}{$exon}{'delBF20'}=$linesplit1[7];
}#end fileline loop
close INPUT;
	
#Get full genenames into hash and then combine all terms for same chr-pos/(gene)
my %full_genenames;
open INPUT_F, $fullname_file or die "Cannot open $fullname_file\n";
	my $head_e =<INPUT_F>;
	loopf: while (<INPUT_F>){
		   my $Line_f=$_;
		   chomp $Line_f;
		   my @linesplit_f = split(/\t/,$Line_f);
		   my $ch=$linesplit_f[0];
		   if($ch!~/^\d/){next loopf;}##only include numeric chrs ie. not HG183_PATCH,MT,X,Y
		   my $stpos=$linesplit_f[1];
		   my $endpos=$linesplit_f[2];
		   my $ensg=$linesplit_f[4];
		   my $goT="";
		   if(defined $linesplit_f[5]){$goT=$linesplit_f[5];}
		   my $wikiD="";
		   if(defined $linesplit_f[6]){$wikiD=$linesplit_f[6];} ###this has been swapped to short genename to see if it all prints!! 
		   my $mimD="";
		   if(defined $linesplit_f[8]){$mimD=$linesplit_f[8];}
		   if($wikiD=~/\S+/ and !exists $full_genenames{$ch}{$stpos}{$endpos}{'wiki'}{$wikiD}){$full_genenames{$ch}{$stpos}{$endpos}{'wiki'}{$wikiD}=0;}
		   if($goT=~/\S+/ and !exists $full_genenames{$ch}{$stpos}{$endpos}{'go'}{$goT}) {$full_genenames{$ch}{$stpos}{$endpos}{'go'}{$goT}=0;}
		   if($mimD=~/\S+/ and !exists $full_genenames{$ch}{$stpos}{$endpos}{'mim'}{$mimD}){$full_genenames{$ch}{$stpos}{$endpos}{'mim'}{$mimD}=0;}
}
close INPUT_F;
my %combined_names;
foreach my $c (keys %full_genenames){
	foreach my $sp (keys %{$full_genenames{$c}}){
		foreach my $ep (keys %{$full_genenames{$c}{$sp}}){
		
		if(exists $full_genenames{$c}{$sp}{$ep}{'wiki'}){	
	foreach my $w (keys %{$full_genenames{$c}{$sp}{$ep}{'wiki'}}){
		if(exists $combined_names{$c}{$sp}{$ep}{'wiki'}){$combined_names{$c}{$sp}{$ep}{'wiki'}=$combined_names{$c}{$sp}{$ep}{'wiki'}.";".$w;}
		if(!exists $combined_names{$c}{$sp}{$ep}{'wiki'}){$combined_names{$c}{$sp}{$ep}{'wiki'}=$w;}
		}}
		if(exists $full_genenames{$c}{$sp}{$ep}{'go'}){
	foreach my $g (keys %{$full_genenames{$c}{$sp}{$ep}{'go'}}){
		if(exists $combined_names{$c}{$sp}{$ep}{'go'}){$combined_names{$c}{$sp}{$ep}{'go'}=$combined_names{$c}{$sp}{$ep}{'go'}.";".$g;}
		if(!exists $combined_names{$c}{$sp}{$ep}{'go'}){$combined_names{$c}{$sp}{$ep}{'go'}=$g;}
		}}
		if(exists $full_genenames{$c}{$sp}{$ep}{'mim'}){
	foreach my $m (keys %{$full_genenames{$c}{$sp}{$ep}{'mim'}}){
		if(exists $combined_names{$c}{$sp}{$ep}{'mim'}){$combined_names{$c}{$sp}{$ep}{'mim'}=$combined_names{$c}{$sp}{$ep}{'mim'}.";".$m;}
		if(!exists $combined_names{$c}{$sp}{$ep}{'mim'}){$combined_names{$c}{$sp}{$ep}{'mim'}=$m;}
		}}
		
}}}

#Annotate ExomeDepth CNV results files with in-house control CNV counts and OMIM gene info
foreach my $f (@DirContent){
	chomp $f;
	if($f=~/^CNV_\S+.csv$/){
	my $outputfile2 = $CurrentDir."/Annotated_".$f."2";
	unless(open(OUT2,">$outputfile2")){print "Cannot open file \"$outputfile2\" to write to!\n"; exit;}
	my $infile=$CurrentDir."/".$f;
	open INPUT2, $infile or die "Cannot open $infile\n";
	my $head2=<INPUT2>;
	chomp $head2;
	my @headsplit=split(/,/,$head2);
	for(my $hs=0;$hs<scalar(@headsplit);$hs++){
		print OUT2 "$headsplit[$hs]\t";
	}
	print OUT2 "AverageControls_DupDel_acrossCNV\tmin_controls\tmax_controls\tFullGenesNames\tGO-terms\tOMIM\n";
	fileline2: foreach my $Line2 (<INPUT2>){
		chomp $Line2;
		my @linesplit2 = split (/,/,$Line2);
		my $chr = $linesplit2[6];
		$chr=~s/"//g;
		my $pos_st = 0;
		my $pos_end = 0;
		$pos_st = $linesplit2[4];
		$pos_end = $linesplit2[5];
		my $dup_del = $linesplit2[2];
		$dup_del=~s/"//g;
		my $BF = $linesplit2[8];
		my $num_exons = $linesplit2[3];
		my $total_elements=scalar(@linesplit2);
		
		#get lists of all full gene names, wiki gene descriptions and omim disorders for chr-stpos-endpos coordinates
		my $all_genenames="";
		my $all_go="";
		my $all_omims="";
		
		foreach my $st (sort {$a<=>$b} keys %{$combined_names{$chr}}){
			foreach my $ed (sort {$a<=>$b} keys %{$combined_names{$chr}{$st}}){
				if(($ed>=$pos_st and $ed<=$pos_end) or ($st>=$pos_st and $st<=$pos_end 
						or ($pos_st>=$st and $pos_st<=$ed) or ($pos_end>=$st and $pos_end<=$ed))){
					if(exists $combined_names{$chr}{$st}{$ed}{'wiki'}){$all_genenames=$all_genenames.":".$combined_names{$chr}{$st}{$ed}{'wiki'};}
					if(exists $combined_names{$chr}{$st}{$ed}{'go'}){$all_go=$all_go.":".$combined_names{$chr}{$st}{$ed}{'go'};}
					if(exists $combined_names{$chr}{$st}{$ed}{'mim'}){$all_omims=$all_omims.":".$combined_names{$chr}{$st}{$ed}{'mim'};}			
					#exit clauses??
					}
			}}
		
		my $avExonsTotal=0;
		my $min=500;
		my $max=0;
		my $count=0;
	exonloop2:	for(my $e=($total_elements-1);$e>=($total_elements-$num_exons);$e--){
			if(!defined $linesplit2[$e]){next exonloop2;}
			my $exon=$linesplit2[$e];
			$exon=~s/"//g;
			
			#don't include 'exons' which are ExomeDepth file errors
			if($exon=~/^\d/ or $exon=~/^CNVR\d/ or $exon=~/^chr\d/ or $exon=~/^duplication/ 
				or $exon=~/^deletion/ or $exon=~/^NA$/){next exonloop2;}
			######												#######
			if(!exists $Variants{$chr}{$exon}){$count++;$min=0;next exonloop2;}
				$avExonsTotal+=$Variants{$chr}{$exon}{'all'};
				if($min>$Variants{$chr}{$exon}{'all'}){$min=$Variants{$chr}{$exon}{'all'};}
				if($max<$Variants{$chr}{$exon}{'all'}){$max=$Variants{$chr}{$exon}{'all'};}
				$count++;
			}#for exons
			my $avExons=0;
			if($count>0){$avExons=$avExonsTotal/$count;}
			
			for(my $pl=0;$pl<scalar(@linesplit2);$pl++){
			print OUT2 "$linesplit2[$pl]\t";
			}
			print OUT2 "$avExons\t$min\t$max\t$all_genenames\t$all_go\t$all_omims\n";
	}#end fileline loop
	close INPUT2;
	close OUT2;
}#if matches
}#foreach file

exit; 