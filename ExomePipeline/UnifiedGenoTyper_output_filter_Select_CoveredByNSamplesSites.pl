#!/usr/bin/perl;
use strict;
use warnings;
use Getopt::Long;

my ($vcf, $intervals, $coverage, $percentage);

my $help;
$coverage = 10;
$percentage = 0.9;

usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, 'vcf=s' => \$vcf, 'intervals=s' => \$intervals, 'coverage=s' => \$coverage, 'percentage=s' => \$percentage) || defined $help );

unless (defined $vcf) {
	die "You have not supplied the input file using --vcf\n";
}

unless (defined $intervals) {
	die "You have not supplied a file to ouput intervals using --intervals\n";
}

my $dp_index=0;
my $total_number_samples = 0;

open INTERVALS, ">$intervals" or die "Cannot open the file $intervals to ouput."; 
open VCF, "$vcf" or die "Can not open the file $vcf";
while (my $line = <VCF> ) {
	if($line =~ m/^\#/) {
		next;
	} else {
		my @elements  = split ("\t", $line);
		if ($elements[8] !~ m/\:DP\:/) {
			print "Error, no DP in INFO field, check input.\n";
			exit;
		} else {
			if( $dp_index == 0 ) { # find out the index of DP in INFO field
				my @format = split (":", $elements[8]);
				for (my $i=0;$i<=$#format;$i++) {
					if($format[$i] eq "DP") {
						$dp_index = $i;
					}
				}
				my $temp = $dp_index+1;
				print "DP is the ".$temp." element of INFO field\n";
			}
			
			if ( $total_number_samples == 0 ) {
				$total_number_samples = $#elements-9+1;
				print "There are $total_number_samples samples in the vcf file.\n";
			}
			
			my $chr  = $elements[0];
			my $pos = $elements[1];
			my $number_passed_coverage = 0;
			for (my $i=9;$i <= $#elements; $i++) {
				if($elements[$i] =~ m/\.\/\./) {
					if ($coverage == 0 ) {
						$number_passed_coverage +=1;
					}
				} else {
					my @temp = split(":", $elements[$i]);
					if ( $temp[$dp_index] >= $coverage ) {
						$number_passed_coverage += 1;
					}
				}
			}
			
			if  ($number_passed_coverage/$total_number_samples >= $percentage ) {
				print INTERVALS $chr.":".$pos."\n";
			} else {
				next;
			}
		}
	}
}

close (VCF);
close (INTERVALS);
exit;

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: UnifiedGenoTyper_output_filter_Select_CoveredByNSamplesSites.pl [--vcf input.vcf][--intervals intervals.output] [--coverage coverage cutoff]  [--percentage percentage of samples has the coverage] [-help|-?]\n\n";
	return(1);
}

