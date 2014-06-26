#!/usr/bin/perl;
use strict;
use warnings;
use Getopt::Long;
use lib qw(/users/a5907529/lib/5.16.1
            /users/a5907529/lib/site_perl/5.16.1);
use Statistics::Descriptive;

my $vcf;
my $outputTail="MQandQD.eval";
my $help;



usage() if ( ! GetOptions('help|?' => \$help, 'vcf=s' => \$vcf, 'outputTail=s' => \$outputTail) || defined $help );

my $output_file= $vcf.".".$outputTail;

unless (defined $vcf) {
	die "You have not supplied the input file using --vcf\n";
}

my @QDs_snp;
my @QDs_indel;
my @MQs_snp;
my @MQs_indel;

open VCF, "$vcf" or die "Can not open the file $vcf";
while (my $line = <VCF> ) {
	if($line =~ m/^\#/) {
		next;
	} else {
		my @elements  = split ("\t", $line);
		
		my $ref=$elements[3];
		my $alt=$elements[4];
		
		if ($ref =~ m/^[atgcATGC]$/ && $alt =~ m/^[atgcATGC]$/) { # confirm it's a SNP not INDEL
			if ($elements[7] =~ m/QD\=(\d+\.?\d*)\;.*/) {
					push @QDs_snp, $1;
			}		
			if ($elements[7] =~ m/MQ\=(\d+\.?\d*)\;.*/) {
				push @MQs_snp, $1;
			}
		} else {
			if ($elements[7] =~ m/QD\=(\d+\.?\d*)\;.*/) {
					push @QDs_indel, $1;
			}
			if ($elements[7] =~ m/MQ\=(\d+\.?\d*)\;.*/) {
					push @MQs_indel, $1;
			}
		}
	}
}
close (VCF);

my $QDs_snp_quantiles_ref = QuantilesOfArray(\@QDs_snp);
my $MQs_snp_quantiles_ref = QuantilesOfArray(\@MQs_snp);
my $QDs_indel_quantiles_ref = QuantilesOfArray(\@QDs_indel);
my $MQs_indel_quantiles_ref = QuantilesOfArray(\@MQs_indel);

print "nQDs_snp ".scalar(@QDs_snp)."\n";
print "nMQs_snp ".scalar(@MQs_snp)."\n";
print "nQDs_indel ".scalar(@QDs_indel)."\n";
print "nMQs_indel ".scalar(@MQs_indel)."\n";


open OUTPUT, ">$output_file" or die "Cannot open the file $output_file to ouput."; 
print OUTPUT "Name\tMin\t1st quantils\tMedian\tMean\t3rd Quantile\tMax\n";
print OUTPUT "QDs_snp\t";
foreach my $quan (@$QDs_snp_quantiles_ref) {
	print OUTPUT $quan."\t";
}
print OUTPUT "\n";
print OUTPUT "MQs_snp\t";
foreach my $quan (@$MQs_snp_quantiles_ref) {
	print OUTPUT $quan."\t";
}
print OUTPUT "\n";
print OUTPUT "QDs_indel\t";
foreach my $quan (@$QDs_indel_quantiles_ref) {
	print OUTPUT $quan."\t";
}
print OUTPUT "\n";
print OUTPUT "MQs_indel\t";
foreach my $quan (@$MQs_indel_quantiles_ref) {
	print OUTPUT $quan."\t";
}
print OUTPUT "\n";
close (OUTPUT);

exit;

sub QuantilesOfArray {
	my ($array_ref) = @_;
	my @quantiles;
	my @numbers = @$array_ref;
	my $stat = Statistics::Descriptive::Full->new( );
	$stat->add_data(@numbers);
	push @quantiles,  sprintf( "%.1f", $stat->quantile(0));
	push @quantiles,  sprintf( "%.1f", $stat->quantile(1));
	push @quantiles,  sprintf( "%.1f", $stat->quantile(2));
	push @quantiles,  sprintf( "%.1f", $stat->mean());
	push @quantiles,  sprintf( "%.1f", $stat->quantile(3));
	push @quantiles,  sprintf( "%.1f", $stat->quantile(4));
	return \@quantiles;
}	

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: UnifiedGenoTyper_output_filter_Select_CoveredByNSamplesSites.pl [--vcf input.vcf][--outputTail tail string adding to the input as output file name] [-help|-?]\n\n";
	return(1);
}

