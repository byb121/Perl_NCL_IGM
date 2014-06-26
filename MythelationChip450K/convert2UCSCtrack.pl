#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(ceil floor);

my ($input) = @ARGV;
my $output = $input."_processed4upload.bed";
my @output;

############# 1st line for the track #############
#############     edit as needed     #############
my $first_line = "track name=\"d0 vs d14 Mythelation\" description=\"d0_vs_d14 _MSC_3Samples_each\" visibility=3 itemRgb=\"On\"\n";
push @output,$first_line;

print "\n";
print "#######################################################\n";
print "#   This script will convert the processed chip data  #\n"; 
print "#     into track file to display in Genome Browser    #\n";
print "#######################################################\n";
print "\n";
print "\n";

print "Reading file: $input\n";

open INPUT, $input or die "cannot open file $input.\n";

while (my $line =<INPUT>) {
	if ($line =~ m/^TargetID/) {
		next;
	} else {
		my @words = split( "\t", $line);
		my $chr = "chr".$words[27];
		my $start;
		if($words[28] eq ""){ #remove those lines with no mapping postitions
			next;
		} else {
			$start = $words[28];
		}
		my $end = $words[28] + 50;
		my $name;
		if ($words[40] eq "" && $words[41] eq "") {
			$name = $chr."_".$start;
		} else {
			$name = $words[40]."_".$words[41];
		}
		my $score = floor(($words[7] + $words[8] + $words[9])*1000/3 + 0.5); 
		my $strand;
		if ($words[32] eq "F"){
			$strand = "+";
		} else {
			$strand = "-";
		}
		my $thickStart = $start;
		my $thickEnd = $end;
		################# Define Colour Scheme ####################
		my $colour;
		if ($score >= 600 ) {
			$colour = "255,127,0"; ##Orange
		} elsif ($score > 200 && $score < 600) {
			$colour = "128,0,128"; ##purple
		} elsif ($score > 0 && $score <= 200) {
			$colour = "0,0,205" ; ##bright blue
		} else {
			$colour = "0,0,0"; ##Black
		}
		
		################# for comparison ##################
		my $c_score = -floor($words[1]*1000/3 + 0.5); 
		
		################ Colour Scheme for comparison #################
		my $c_colour;
		my $p_value = $words[5];
		if ($p_value <= 0.05 && $c_score > 0) {
			$c_colour = "0,238,0"; ##Green
		} elsif ($p_value <= 0.05 && $c_score < 0) {
			$c_colour = "238,0,0"; ##Red			
		} else {
			$c_colour = "0,0,0"; ##Black
			
		}
		
		#push @output, $chr."\t".$start."\t".$end."\t".$name."\t".$score."\t".$strand."\t".$thickStart."\t".$thickEnd."\t".$colour."\n";
		push @output, $chr."\t".$start."\t".$end."\t".$name."\t".$c_score."\t".$strand."\t".$thickStart."\t".$thickEnd."\t".$c_colour."\n";
	}
}
close(INPUT);

print "Writing file: $output\n";
open OUTPUT, ">$output" or die "cannot open file $output.\n";
print OUTPUT @output; 
close OUTPUT;

exit;
