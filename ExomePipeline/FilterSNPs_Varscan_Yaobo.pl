#!/usr/bin/perl;
use strict;
use warnings;
use Getopt::Long;

####
# Adapted from Helen Griffin's script
#####

my $Targets="/users/a5907529/lustre/Yaobo/GenomeData/TruSeq-Exome-Targeted-Regions.txt";

my $File="filename";
my $Hom_het_all="all"; ##all, homozygous or heterozygous##				
my $extra_bp=500; #!!!!! +/-?bp to each target !!!!! 
my $MinCov=5;	
my $var_freq_more_than=25; #$var_freq > 25
my $strands_min=2; #strands >= 2

my $Results=GetOptions("inFile=s"=>\$File, "targets=s"=>\$Targets, "hom_het_all=s"=>\$Hom_het_all, "mincov=f"=>\$MinCov, "minvarfreq=i"=>\$var_freq_more_than);

my $Filtered_FileOT = Filter_onTarget($Targets, $File, $Hom_het_all, $extra_bp, $MinCov, $var_freq_more_than, $strands_min);

####SUBROUTINES#################################################################################
sub Filter_onTarget{
my ($Targets, $File, $hha, $extra_bp, $MinCov, $var_freq_more_than, $strands_min)=@_;

my $coord;
my %Targetregions;
my %Target_chr;
my $arraycount;
my $num_chr_coords;
my %startfrom;

#open target regions file
open (FILE, $Targets) || die "Target file not found\n";
while (<FILE>) {
	chomp($_);
	my @targsplit = split(/\t/, $_);
	  my $chr = $targsplit[0];
	  my $T_start= $targsplit[1];
	  my $T_end = $targsplit[2];
	  my $T_info = $targsplit[3];
	  $coord = join('_', $chr, $T_start, $T_end);
		push(@{$Targetregions{$coord}}, $T_start, $T_end, $T_info);
		push(@{$Target_chr{$chr}}, $coord);
}
close FILE;

#open outputfile 
my $outputfileN = $File."OnTartgets_filtered";
 unless (open(SNPSN, ">$outputfileN")) {print "Cannot open file \"$outputfileN\" to write to!\n\n"; exit;}
#Read in Varscan.snp file
unless (open(FILE2, $File)) {print "\nCannot open file $File!\n";}
my $header=<FILE2>;
print SNPSN $header;
foreach my $line (<FILE2>) {
	chomp($line);
		my @SplitLine=split(/\t/,$line);
			my $chrom_name=$SplitLine[0]; 
			my $pos=$SplitLine[1]; 
			my $ref=$SplitLine[2]; 
			my $var=$SplitLine[3]; 
			my $reads1=$SplitLine[4]; 
			my $reads2=$SplitLine[5]; 
			my $var_freq=$SplitLine[6]; 
			my $strands1=$SplitLine[7]; 
			my $strands2=$SplitLine[8]; 
			my $qual1=$SplitLine[9]; 
			my $qual2=$SplitLine[10]; 
			my $pvalue=$SplitLine[11]; 
	
	#filter on minimum coverage
	if(($reads1+$reads2)<$MinCov){next;}
	if($reads2<1){next;}
	#chop $var_freq to remove % symbol
	chomp($var_freq);
	chop($var_freq);

	if(!exists $Target_chr{$chrom_name}){next;} #Not all chr names may be present in targets hash!!
	$num_chr_coords = @{$Target_chr{$chrom_name}};
	if(!exists $startfrom{$chrom_name}){$startfrom{$chrom_name}=0;}
	
	if($var_freq > $var_freq_more_than && $strands2 >= $strands_min){
		
	#foreach chr specific target region, print snps within region to file
	loop: for($arraycount = $startfrom{$chrom_name}; $arraycount < $num_chr_coords; $arraycount++){
		
		my $chr_coord = $Target_chr{$chrom_name}[$arraycount];
	 	my $t_start = $Targetregions{$chr_coord}[0] - $extra_bp; #Add bases to either end of target coordinates
		my $t_end = $Targetregions{$chr_coord}[1] + $extra_bp;	
		my $gene_name = $Targetregions{$chr_coord}[2];
		my $t_length = ($t_end+1) - $t_start;

	#convert 1-based position ($pos) to 0-based position ($pos_minus1)
		my $pos_minus1 = $pos - 1;

	  if($pos_minus1 < $t_start){$startfrom{$chrom_name}=$arraycount; last loop;}
	  
	  if($t_start <= $pos_minus1 && $pos_minus1 <= $t_end){
			if($hha=~/homozygous/ and $var_freq>=85){print SNPSN $line."\n";}
			if($hha=~/heterozygous/ and $var_freq<85){print SNPSN $line."\n";}
			if($hha=~/all/){print SNPSN $line."\n";}
			$startfrom{$chrom_name}=$arraycount;
			last loop;
			}
	   	} #arraycount for loop
	  } #varfreq > 25 && strands2 >= 2 loop
}
close FILE2;
close SNPSN;
return($outputfileN);
}

exit;
