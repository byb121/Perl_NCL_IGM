#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $targets_file;
my $output_file;
my $population="EUR";

my $help;
usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help,  'target=s' => \$targets_file,  'population=s' => \$population, 'output=s' => \$output_file) || defined $help );

print "##########################################################\n";
print "Target regions require the following format for each row: chr\tstart\tend\n";
print "chr: the chromosome name. need to start with \"chr\".\n";
print "start: start position of a region. 0-based. check the bed format for details.\n";
print "end: end position of a region, 1-based. check the bed format for details.\n";
print "##########################################################\n";

unless ($population eq "AMR" || $population eq "AFR" || $population eq "ASN" || $population eq "EUR") {
	print "Unrecognise variant_type. Accept values are AMR, AFR, ASR and EUR. Aborted! Retry!!\n";
	print "For more details: http://www.1000genomes.org/category/frequently-asked-questions/population \n";
	exit;
}

my $population_raf_file;
if( $population eq "AMR") {
	$population_raf_file="/users/a5907529/lustre/1000Geome_Phase1_Calls/AMR_RAF_table.txt";
} elsif ($population eq "AFR") {
	$population_raf_file="/users/a5907529/lustre/1000Geome_Phase1_Calls/AFR_RAF_table.txt";
} elsif ($population eq "ASN") {
	$population_raf_file="/users/a5907529/lustre/1000Geome_Phase1_Calls/ASN_RAF_table.txt";
} elsif ($population eq "EUR") {
	$population_raf_file="/users/a5907529/lustre/1000Geome_Phase1_Calls/EUR_RAF_table.txt";
}

#get the target region
my %Targetregions;
open (TARGET, $targets_file) || die "Target file not found: $targets_file\n";
print "Reading in targets....\n";
while (my $line = <TARGET>) {
	chomp $line;
	my @targsplit = split(/\t/, $line);
	my $chr = $targsplit[0];
	my $start= $targsplit[1]+1;
	my $end = $targsplit[2];
	if ( !exists$Targetregions{$chr}{$start}) {
		$Targetregions{$chr}{$start} = $end; 
	} else {
		if ($end > $Targetregions{$chr}{$start}) {
			$Targetregions{$chr}{$start} = $end;
		}
	}
}
close TARGET;
print "Done.\n";

my %Targets_starts;
my %Targets_ends;

print "Sorting target regions....\n";
foreach my $chr (keys %Targetregions) {
	my @starts;
	my @ends;
	foreach my $start (keys %{$Targetregions{$chr}}) {
		push @starts, $start;
		push @ends, $Targetregions{$chr}{$start};
	}
	my @idx =  sort { $starts[$a] <=> $starts[$b] } 0 .. $#starts;
	@starts = @starts[@idx];
	@ends = @ends[@idx];
	$Targets_starts{$chr} = \@starts;
	$Targets_ends{$chr} = \@ends;
	print "for $chr, there are ". scalar @starts." starts and ".scalar @ends." ends.\n";
}
undef %Targetregions;
print "Done.\n";

#read the raf file
open RAF, $population_raf_file or die "Cannot open the file: $population_raf_file. Please make sure it exists and readable.\n";
open OUTPUT, ">$output_file" or die "Cannnot open the file to ouput; $output_file\n";

my $my_running_chr = "";
my @running_targets_starts;
my @running_targets_ends;
my $running_reached_index = 0;

while (my $line = <RAF> ) {
	chomp $line;
	my @temp = split("\t", $line);
	my $chr = $temp[2];
	my $pos = $temp[3];
	
	if ($my_running_chr eq "") {
		print "Start to run on $chr\n";
		$my_running_chr = $chr;
		@running_targets_starts = @{$Targets_starts{$chr}};
		@running_targets_ends = @{$Targets_ends{$chr}};
		$running_reached_index = 0;
	} else {
		if ( $chr ne $my_running_chr ) {
			print "Currently running on $chr\n";
			$my_running_chr = $chr;
			@running_targets_starts = @{$Targets_starts{$chr}};
			@running_targets_ends = @{$Targets_ends{$chr}};
			$running_reached_index = 0;
		}
	}
	
	start_loop: for ( my $i = $running_reached_index; $i < scalar @running_targets_starts; $i++ ) {
		if ( $pos >= $running_targets_starts[$i] && $pos <= $running_targets_ends[$i]) {
			print OUTPUT $line."\n";
			$running_reached_index = $i;
			last start_loop;
		}
	}
}
close RAF;
close OUTPUT;
exit;

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: ProduceRefSet4SensAndSpec.pl [--target /Path/to/target_region.bed] [--population \"AFR|AMR|ASN|EUR\"] [--output /Path/to/output.txt] [--help|-?]\n\n";
	return(1);
}

