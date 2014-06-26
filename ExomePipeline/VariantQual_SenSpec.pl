#!/usr/bin/perl;
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);	
# use IPC::System::Simple qw(system);
	
my ($help, $vcf, $ref);
my $currentdir=`pwd`;
my $CurrentScript=abs_path($0);
$CurrentScript=~/(.*)\/VariantQual/;
my $scriptdir=$1; #="~/scripts/darren";#### Path to R script ###
chomp $currentdir;

my $home_dir = $ENV{"HOME"};
chomp $home_dir;


usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, 'vcfFile=s' => \$vcf, 'refSNPs=s' => \$ref) || defined $help );

unless (defined $vcf) {
	die "You have not supplied a vcf file using -vcfFile\n";
}

unless (defined $ref) {
	die "You have not supplied a reference list of SNPs using -refSNPs\n";
}

$vcf =~ s/^\~/$home_dir/;
$vcf = abs_path($vcf);

my $input=Filter($vcf, $ref);
r("$input");
exit;

####SUBROUTINES#################################################################################

sub Filter{ 
	my $File=$_[0];
	my $ref=$_[1];
	my %VariantHash=();
	my @Array=();
	
	open FILE, "$File";
	while(<FILE>){
		chomp $_;
		if($_=~/\#/){
			next;
		}else{
			my @SplitLine=split(/\t/, $_);
			my $Chr=$SplitLine[0];
			unless($SplitLine[0]=~/chr\S+/){
				$Chr="chr".$SplitLine[0];
			}
			my $Pos=$SplitLine[1];
			my $Ref=uc($SplitLine[3]);
			my $Variant=uc($SplitLine[4]);
			my @SplitInfo=split(';', $SplitLine[7]);
			my @SplitGT=split(':', $SplitLine[9]);
			my $GT=$SplitGT[0];
			my $Match=$Chr."_".$Pos."_".$Ref;
			if($Variant=~/\S+\,/){
				my @SplitVar=split(',', $Variant);
				foreach my $Var(@SplitVar){
					if($GT eq "1/1"){
						$VariantHash{$Match}="3";
					}elsif($GT eq "0/1"){
						$VariantHash{$Match}="2";
					}
				}
			}else{
				if($GT eq "1/1"){
					$VariantHash{$Match}="3";
				}elsif($GT eq "0/1"){
					$VariantHash{$Match}="2";
				}					
			}

		}
	}	

	my $rinput="$File.Input.txt";
	open I, ">$rinput";
	print I "rsNumber\tChange\tChromosome\tPosition\tFrequency\tSample\n";
	
	open REF, "$ref";
	while(<REF>){
		chomp $_;
		if($_=~/rsNumber\t/){
			next;
		}else{
			my @Split=split('\t', $_);
			my @splitrefvar=split('/', $Split[1]);
			my $VariantMatch=$Split[2]."_".$Split[3]."_".$splitrefvar[0];
			if(exists $VariantHash{$VariantMatch}){
				print I $Split[0]."\t".$Split[1]."\t".$Split[2]."\t".$Split[3]."\t".$Split[4]."\t".$VariantHash{$VariantMatch}."\n";
			}else{		
				print I $Split[0]."\t".$Split[1]."\t".$Split[2]."\t".$Split[3]."\t".$Split[4]."\t"."1"."\n";
			}
		}
	}
	
	close I;
	return($rinput);
}

sub r{
	my $in=$_[0];
	my $invokeR = "R --vanilla --silent --no-save --args";
	my $scriptR = "$scriptdir/Darrens.r > $in.R.log";
	my $Garbage=system("$invokeR $in < $scriptR ");
	system("rm $in $in.R.log");
	return(1);
}

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: AnalysisScript_v1.pl [-vcfFile vcf file] [-refSNPs reference list of SNPs] [-help|-?]\n\n";
	return(1);
}

