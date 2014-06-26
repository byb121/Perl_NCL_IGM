#!/usr/bin/perl
# Created 15/10/2013
# Edited 23/10/2013,31/10/2013
use warnings;
use strict;
use Getopt::Long;
use List::Util qw(sum);
use List::Util qw(max);
{
	# Parameters
	my $OutputDirectory="";
	my $InputDirectory="";
	my $PileupFile="short.pileup";
	my $OutputFile="test.txt";
	my $MinDepth=10;
	my $MinVarCnt=0;
	my $MinQual=30;
	my $Parameters= GetOptions("pileupfile:s"=>\$PileupFile,
							"outputfile:s"=>\$OutputFile,
							"indir:s"=>\$InputDirectory,
							"outdir:s"=>\$OutputDirectory,
							"mindepth:i"=>\$MinDepth,
							"minvarcnt:i"=>\$MinVarCnt,
							"minqual:i"=>\$MinQual);
	$PileupFile=$InputDirectory.$PileupFile;
	$OutputFile=$OutputDirectory.$OutputFile;
	
	my ($DepthRef)=readPileAndExclude($PileupFile,$OutputFile,$MinDepth,$MinVarCnt,$MinQual);
	exit;
}

sub readPileAndExclude{
	my ($PileupFileName,$OutputFileName,$MinDepth,$MinVarCnt,$MinQual)=@_;
	open PileupFile, $PileupFileName or die "Pileupfile $PileupFileName not found \n";
	my %ChrSize=(Chr1=>5579133, Chr2=>4509021, Chr3=>2452883);
	my @RefChrom;
	my @RefBase;
	my @MutBase;
	my @MutPos;
	my @MutCnt;
	my @RefCnt;
	my %Frequencies;
	while (<PileupFile>) {
		my $Line=$_;
		chomp($Line);
		next if $Line=~/^#/ ; # skip comment lines
		my @Bits=split /\t/,$Line;
		next if !exists $ChrSize{$Bits[0]};
		next if $Bits[1]>$ChrSize{$Bits[0]};
		my $CurrentDepth=$Bits[3];
		next if $CurrentDepth<$MinDepth;
		my ($BestVar,$BestVarCnt,$WTCnt)=processCnts($Bits[5],$Bits[4],$MinQual,$MinVarCnt);
		push(@RefChrom,$Bits[0]);
		push(@RefBase,$Bits[2]);
        push(@MutBase,$BestVar);
		push(@MutPos,$Bits[1]);
		push(@MutCnt,$BestVarCnt);
		push(@RefCnt,$WTCnt);
	}
	close PileupFile;
	open OutputFile, ">$OutputFileName"; 
	for (my $i=0;$i<scalar @MutBase;$i++)  {
	next if (($RefCnt[$i]+$MutCnt[$i])<$MinDepth);
	next if ($RefBase[$i]!~/[ACGTacgt]/);
	print OutputFile $RefChrom[$i]."\t".$MutPos[$i]."\t".uc($RefBase[$i]);
	print OutputFile "\t".uc($MutBase[$i])."\t".$RefCnt[$i]."\t".$MutCnt[$i]."\n";
	}
	close OutputFile;
	return;
}

sub processCnts{
   my($QualsString,$BaseString,$MinQual,$MinVarCnt)=@_;
   # Exclude indels and read begin/end fields
   $BaseString=~s/\^\S//g;
   $BaseString=~s/\$//g;
   $BaseString=~s/\-[0-9]+[ACGTNacgtn]+//g;
   $BaseString=~s/\+[0-9]+[ACGTNacgtn]+//g;
   # Exclude bases below quality
   my @BasesArray=split //,$BaseString;
   my @QualArray=split //,$QualsString;
   my @CleanBases;
   for (my $i=0;$i< scalar @BasesArray; $i++){
	next if ((ord($QualArray[$i])-33)<$MinQual);
	push @CleanBases,$BasesArray[$i];
   }
   my %Counts=('.'=>0,','=>0);
   my @BaseCounts;
   my @Bases=qw(a c g t);
   ++$Counts{$_} for @CleanBases;
   my @UniqueBases = keys(%Counts);
   # Count variant frequencies
   my $BestVarCnt=0;
   my $BestVar=" ";
   for (my $i=0;$i<4;$i++){
	$BaseCounts[$i]=0;
	if (exists($Counts{$Bases[$i]})) {$BaseCounts[$i]=$Counts{$Bases[$i]};}
	if (exists($Counts{uc($Bases[$i])})) {$BaseCounts[$i]+=$Counts{uc($Bases[$i])};}
	if (($BaseCounts[$i]>$BestVarCnt) && ($BaseCounts[$i]>=$MinVarCnt))  {
	    $BestVarCnt=$BaseCounts[$i];
		$BestVar=uc($Bases[$i]);
	}
   }
   # Count wildtype
   my $WTCnt=$Counts{"."}+$Counts{","};
   return($BestVar,$BestVarCnt,$WTCnt)
}