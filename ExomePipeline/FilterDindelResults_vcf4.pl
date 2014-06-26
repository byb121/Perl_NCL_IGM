#!/usr/bin/perl;
use strict;
use warnings;
use Getopt::Long;

###############
#
# origin: Helen Griffin. Adopted to produce vcf format
#
###############

#my $Targets="/users/nhrg/lustre/hg19/Agilentwholeexome50Mb_Targetshg19_GENE_sortednonOV.txt";
my $Targets="/users/a5907529/lustre/Yaobo/GenomeData/TruSeq-Exome-Targeted-Regions.txt";

my $File="filename";														
my $Hom_het_all="all"; ##all, hom or het##				
my $qual=10;																		
my $RF=2;
my $RR=2;
my $extra_bp=500; #+/-?bp to each target 

my $Results=GetOptions("inFile=s"=>\$File, "targets=s"=>\$Targets, "hom_het_all=s"=>\$Hom_het_all, "qual=f"=>\$qual, "reads_F=i"=>\$RF, "reads_R=i"=>\$RR);

if($File=~/(\S+\/)(\S+)_variantCalls.vcf4/){
	my $Filepath=$1;
	my $ID=$2;
	my $Filtered_File1 = Filter_cov_strands_qual($Filepath, $File, $Hom_het_all, $ID, $qual, $RF, $RR);
	my $Filtered_File2 = Filter_onTarget($Filepath, $Targets, $Filtered_File1, $Hom_het_all, $ID, $extra_bp);
}

if($File=~/(\S+\/)(\S+)_variantCalls.vcf4/){
	my $Filepath=$1;
	my $ID=$2;
	my $Filtered_File1 = Filter_cov_strands_qual_outVCF($Filepath, $File, $Hom_het_all, $ID, $qual, $RF, $RR);
	my $Filtered_File2 = Filter_onTarget_outVCF($Filepath, $Targets, $Filtered_File1, $Hom_het_all, $ID, $extra_bp);
}

####SUBROUTINES#################################################################################

sub Filter_cov_strands_qual{ ##filter out SNPs/indels not found on both strands and with a minimum coverage of less than 2 per strand
	my ($Filepath, $File, $hha, $ID, $qual, $RF, $RR)=@_;
	my $NumberForwardReads="";
	my $NumberReverseReads="";
	my $TotalCov="";
	my %Hash=();
		my $Out=$Filepath.$ID."_FilteredIndels_vcf4_qual_cov_".$hha;
		open FILE, $File;
		while(<FILE>){
			chomp $_;
			if($_=~/chr\S+/){
				my @SplitLine=split(/\t/, $_);
				my $Chr=$SplitLine[0];
				my $Pos=$SplitLine[1];
				my $Ref=uc($SplitLine[3]);
				my $Variant=uc($SplitLine[4]);
				
				my $R=$Ref;
				my $Ref_firstbase = substr($R,0,1,q{}); # use 1st base of $Ref as ref!!!	
				#print "1st base: $Ref_firstbase\tRef: $Ref\tVariant: $Variant\t Pos_orig: $Pos\t";
				#variant??? add 1 or more!! ins/dels to @Vars
				my @Vars;
				my @Pos_s;
				if(length $Ref==1){ ##insertion(s)
					if($Variant!~/\,/){
						$Variant =~ s/^./\+/;
						push(@Vars,$Variant);
						push(@Pos_s,$Pos);
					}elsif($Variant=~/\,/){
						my @varsplit=split(/\,/,$Variant);
						foreach my $v (@varsplit){
							$v =~ s/^./\+/;
							push(@Vars, $v);
							push(@Pos_s,$Pos);
							}
						}
				}elsif(length $Ref>1){ ##deletion(s) !!and in some cases of multiple (,) vars insertions!! 
					if($Variant!~/\,/){
					my $v = $Ref;
					$v =~ s/^./\-/;
					push(@Vars,$v);
					push(@Pos_s,$Pos);
					}elsif($Variant=~/\,/){
						my @varsplit=split(/\,/,$Variant);
						foreach my $v (@varsplit){ 			##N.B. if length of ref equals length of var this is currently ignored!!!!
							#print "\nv: $v\n";
							#if(length $v == 1){my $v1=$Ref; $v1 =~ s/^./\-/;push(@Vars,$v1); push(@Pos_s,$Pos);} ###not needed as < below does same!!
							if(length $v < length $Ref){my $diff=(length $v)-1; my $v2= substr($Ref,$diff); $v2=~s/^./\-/; push(@Vars,$v2); my $Pos_a = $Pos+$diff; push(@Pos_s,$Pos_a);}
							if(length $v > length $Ref){my $diff2=(length $Ref)-1; my $v3= substr($v,$diff2); $v3=~s/^./\+/; push(@Vars,$v3); my $Pos_b = $Pos+$diff2; push(@Pos_s,$Pos_b);}
							}
						}
					}
				
				#print "Vars:  @Vars\tPos_s: @Pos_s\n";
				
				my $genotype_qual=$SplitLine[9]; #e.g. 0/1:5 or 1/1:90
				my @gt_split=split(/:/,$genotype_qual);
				my $gt_qual=$gt_split[1];
				my @Info=split(/;/, $SplitLine[7]); #e.g. DP=52;NF=5;NR=0;NRS=15;NFS=0;HP=4
				foreach my $Info(@Info){
					if($Info=~/NF=(\S+)/){
						$NumberForwardReads=$1;
					}elsif($Info=~/NR=(\S+)/){
						$NumberReverseReads=$1;
					}elsif($Info=~/DP=(\S+)/){
						$TotalCov=$1;}
				}
				my $c=0;
				foreach my $entry (@Vars){
				my $VariantLength=length $entry;
				  if($gt_qual>=$qual){ ##filter on genotype quality	
					if(($NumberForwardReads>=$RF)&&($NumberReverseReads>=$RR)){ ##filter on minimum variant read coverage for F & R strands 
					   if($hha eq "het" and $genotype_qual!~/1\/1/){ #heterozygous variants only
					   #print "het: $hha and $genotype_qual\n";
					   $Hash{$Chr}{$Pos_s[$c]}{$Ref_firstbase}{$entry}=$NumberForwardReads."\t".$NumberReverseReads."\t".$genotype_qual;
				       }
					   if($hha eq "hom" and $genotype_qual=~/1\/1/){ #homozygous variants only
					   #print "hom: $hha and $genotype_qual\n";
					   $Hash{$Chr}{$Pos_s[$c]}{$Ref_firstbase}{$entry}=$NumberForwardReads."\t".$NumberReverseReads."\t".$genotype_qual;
				       }
					   if($hha eq "all"){ #both hom & het variants
					   #print "all: $hha and $genotype_qual\n";
					   $Hash{$Chr}{$Pos_s[$c]}{$Ref_firstbase}{$entry}=$NumberForwardReads."\t".$NumberReverseReads."\t".$genotype_qual;
				       }
					}#min read cov on both forward & reverse
				  }#genotype quality
				#}#if var length>=2, selects only indels
			$c++;
			}#foreach $entry pos/var
			
			}#if line starts with chr
		}
		open OUT, ">$Out";
		foreach my $c (sort keys %Hash){
			foreach my $p (sort {$a<=>$b} keys %{$Hash{$c}}){
				foreach my $r (keys %{$Hash{$c}{$p}}){
					foreach my $v (keys %{$Hash{$c}{$p}{$r}}){
		print OUT "$c\t$p\t$r\t$v\t$Hash{$c}{$p}{$r}{$v}\n";
}}}}
		close OUT;
		%Hash=();
		close FILE;		
return($Out);
}

sub Filter_cov_strands_qual_outVCF { ##filter out SNPs/indels not found on both strands and with a minimum coverage of less than 2 per strand
	my ($Filepath, $File, $hha, $ID, $qual, $RF, $RR)=@_;
	my $NumberForwardReads="";
	my $NumberReverseReads="";
	my $TotalCov="";
	my %Hash=();
	my $Out=$Filepath.$ID."_FilteredIndels_vcf4_qual_cov_".$hha.".vcf";
	open FILE, $File;
	while(my $line=<FILE>){
		chomp $line;
		if($line=~/chr\S+/){
			my @SplitLine=split(/\t/, $line);
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			my $Ref=uc($SplitLine[3]);
			my $Variant=uc($SplitLine[4]);
			
			my $R=$Ref;
			my $Ref_firstbase = substr($R,0,1,q{}); # use 1st base of $Ref as ref!!!	
			#print "1st base: $Ref_firstbase\tRef: $Ref\tVariant: $Variant\t Pos_orig: $Pos\t";
			#variant??? add 1 or more!! ins/dels to @Vars
			my @Vars;
			my @Pos_s;
			if(length $Ref==1){ ##insertion(s)
				if($Variant!~/\,/){
					$Variant =~ s/^./\+/;
					push(@Vars,$Variant);
					push(@Pos_s,$Pos);
				}elsif($Variant=~/\,/){
					my @varsplit=split(/\,/,$Variant);
					foreach my $v (@varsplit){
						$v =~ s/^./\+/;
						push(@Vars, $v);
						push(@Pos_s,$Pos);
					}
				}
			}elsif(length $Ref>1){ ##deletion(s) !!and in some cases of multiple (,) vars insertions!! 
				if($Variant!~/\,/){
					my $v = $Ref;
					$v =~ s/^./\-/;
					push(@Vars,$v);
					push(@Pos_s,$Pos);
				}elsif($Variant=~/\,/){
					my @varsplit=split(/\,/,$Variant);
					foreach my $v (@varsplit){ 			##N.B. if length of ref equals length of var this is currently ignored!!!!
						#print "\nv: $v\n";
						#if(length $v == 1){my $v1=$Ref; $v1 =~ s/^./\-/;push(@Vars,$v1); push(@Pos_s,$Pos);} ###not needed as < below does same!!
						if(length $v < length $Ref){my $diff=(length $v)-1; my $v2= substr($Ref,$diff); $v2=~s/^./\-/; push(@Vars,$v2); my $Pos_a = $Pos+$diff; push(@Pos_s,$Pos_a);}
						if(length $v > length $Ref){my $diff2=(length $Ref)-1; my $v3= substr($v,$diff2); $v3=~s/^./\+/; push(@Vars,$v3); my $Pos_b = $Pos+$diff2; push(@Pos_s,$Pos_b);}
					}
				}
			}
				
			#print "Vars:  @Vars\tPos_s: @Pos_s\n";
				
			my $genotype_qual=$SplitLine[9]; #e.g. 0/1:5 or 1/1:90
			my @gt_split=split(/:/,$genotype_qual);
			my $gt_qual=$gt_split[1];
			my @Info=split(/;/, $SplitLine[7]); #e.g. DP=52;NF=5;NR=0;NRS=15;NFS=0;HP=4
			foreach my $Info(@Info){
				if($Info=~/NF=(\S+)/){
					$NumberForwardReads=$1;
				}elsif($Info=~/NR=(\S+)/){
					$NumberReverseReads=$1;
				}elsif($Info=~/DP=(\S+)/){
					$TotalCov=$1;}
			}
			my $c=0;
			foreach my $entry (@Vars){
				my $VariantLength=length $entry;
				if($gt_qual>=$qual){ ##filter on genotype quality	
					if(($NumberForwardReads>=$RF)&&($NumberReverseReads>=$RR)){ ##filter on minimum variant read coverage for F & R strands 
						if($hha eq "het" and $genotype_qual!~/1\/1/){ #heterozygous variants only
							#print "het: $hha and $genotype_qual\n";
				   			$Hash{$Chr}{$Pos_s[$c]}{$Ref_firstbase}{$entry}=$line;
			       		}
				   		if($hha eq "hom" and $genotype_qual=~/1\/1/){ #homozygous variants only
				   			#print "hom: $hha and $genotype_qual\n";
				   			$Hash{$Chr}{$Pos_s[$c]}{$Ref_firstbase}{$entry}=$line;
			       		}
				   		if($hha eq "all"){ #both hom & het variants
				   			#print "all: $hha and $genotype_qual\n";
				   			$Hash{$Chr}{$Pos_s[$c]}{$Ref_firstbase}{$entry}=$line;
				   		}
					}#min read cov on both forward & reverse
			  	}#genotype quality
				#}#if var length>=2, selects only indels
				$c++;
			}#foreach $entry pos/var
			
		} #if line starts with chr
	}
	open OUT, ">$Out";
	foreach my $c (sort keys %Hash){
		foreach my $p (sort {$a<=>$b} keys %{$Hash{$c}}){
			foreach my $r (keys %{$Hash{$c}{$p}}){
				foreach my $v (keys %{$Hash{$c}{$p}{$r}}){
					print OUT "$Hash{$c}{$p}{$r}{$v}\n";
				}
			}
		}
	}
	
	close OUT;
	%Hash=();
	close FILE;		
	return($Out);
}

sub Filter_onTarget{ ##only include 'on-target' indels
my ($Filepath, $Targets,$File, $hha, $ID, $extra_bp)=@_;
#create target region hash
open (FILE, $Targets) || die "Target File not found\n";
my %Targetregions;
my %Targets_per_chr;
while (<FILE>) {
	chomp($_);
	(my $chr, my $T_start, my $T_end, my $T_info) = split(/\t/, $_);
	if(!exists $Targets_per_chr{$chr}){$Targets_per_chr{$chr}=0;}
	$Targetregions{$chr}[$Targets_per_chr{$chr}][0]=$T_start-$extra_bp;			#####Add extra onto each end of target!
	$Targetregions{$chr}[$Targets_per_chr{$chr}][1]=$T_end+$extra_bp;			####Add extra onto each end of target!
	$Targetregions{$chr}[$Targets_per_chr{$chr}][2]=$T_info;
	$Targets_per_chr{$chr}++;
}
close FILE; 
###Filter indels that are on target.###
	my @Array=();
		my $Out=$Filepath.$ID."_FilteredIndels_OnTarget_vcf4_".$hha;
		my %Startfrom;
		open FILE, $File;
		targetloop: while(<FILE>){
			chomp $_;
			my @SplitLine=split(/\t/, $_);
			my $Chr=$SplitLine[0];
			if(!exists $Targets_per_chr{$Chr}){next;}
			my $Pos=$SplitLine[1];
			my $vlength=(length $SplitLine[3])-2; ##-1 for +/- and -1 for adding to start pos
			my $epos=$Pos+$vlength;##end of indel pos
			#loop through targets hash
			if(!exists $Startfrom{$Chr}){$Startfrom{$Chr}=0;}
			for(my $count=$Startfrom{$Chr}; $count<$Targets_per_chr{$Chr}; $count++){
				if($epos < $Targetregions{$Chr}[$count][0]){$Startfrom{$Chr} = $count; next targetloop;}
				if(($Pos >= $Targetregions{$Chr}[$count][0]) and ($Pos <= $Targetregions{$Chr}[$count][1])
					or ($epos >= $Targetregions{$Chr}[$count][0]) and ($epos <= $Targetregions{$Chr}[$count][1])){
						push @Array, $SplitLine[0]."\t".$SplitLine[1]."\t".$SplitLine[2]."\t".$SplitLine[3]."\t".$SplitLine[4]."\t".$SplitLine[5]."\t".$SplitLine[6]."\t".$Targetregions{$Chr}[$count][2]."\n";		
						$Startfrom{$Chr} = $count;
						next targetloop;
				}#end of if pos match loop
			}#end of targets hash loop
		}#end of indel file loop	
		open OUT, ">$Out";
		print OUT @Array;
		close OUT;
		@Array=();
		close FILE;
return($Out);	
}

sub Filter_onTarget_outVCF { ##only include 'on-target' indels
	my ($Filepath, $Targets,$File, $hha, $ID, $extra_bp)=@_;
	#create target region hash
	open (FILE, $Targets) || die "Target File not found\n";
	my %Targetregions;
	my %Targets_per_chr;
	while (<FILE>) {
		chomp($_);
		(my $chr, my $T_start, my $T_end, my $T_info) = split(/\t/, $_);
		if(!exists $Targets_per_chr{$chr}){$Targets_per_chr{$chr}=0;}
		$Targetregions{$chr}[$Targets_per_chr{$chr}][0]=$T_start-$extra_bp;			#####Add extra onto each end of target!
		$Targetregions{$chr}[$Targets_per_chr{$chr}][1]=$T_end+$extra_bp;			####Add extra onto each end of target!
		$Targetregions{$chr}[$Targets_per_chr{$chr}][2]=$T_info;
		$Targets_per_chr{$chr}++;
	}
	close FILE; 
	###Filter indels that are on target.###
	my @Array=();
	my $Out=$Filepath.$ID."_FilteredIndels_OnTarget_vcf4_".$hha.".vcf";
	my %Startfrom;
	open FILE, $File;
	targetloop: while(<FILE>){
		chomp $_;
		my @SplitLine=split(/\t/, $_);
		my $Chr=$SplitLine[0];
		if(!exists $Targets_per_chr{$Chr}){next;}
		my $Pos=$SplitLine[1];
		my $vlength=(length $SplitLine[3])-2; ##-1 for +/- and -1 for adding to start pos
		my $epos=$Pos+$vlength;##end of indel pos
		#loop through targets hash
		if(!exists $Startfrom{$Chr}){$Startfrom{$Chr}=0;}
		for(my $count=$Startfrom{$Chr}; $count<$Targets_per_chr{$Chr}; $count++){
			if($epos < $Targetregions{$Chr}[$count][0]){$Startfrom{$Chr} = $count; next targetloop;}
			if(($Pos >= $Targetregions{$Chr}[$count][0]) and ($Pos <= $Targetregions{$Chr}[$count][1])
			or ($epos >= $Targetregions{$Chr}[$count][0]) and ($epos <= $Targetregions{$Chr}[$count][1])){
				push @Array, $_."\n";		
				$Startfrom{$Chr} = $count;
				next targetloop;
			}#end of if pos match loop
		}#end of targets hash loop
	}#end of indel file loop	
	open OUT, ">$Out";
	print OUT @Array;
	close OUT;
	@Array=();
	close FILE;
	return($Out);	
}

exit;