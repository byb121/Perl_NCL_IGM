#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $coordinates="";
my $vcf_file = "";
my $Results=GetOptions("coordinates=s"=>\$coordinates, "vcf=s"=>\$vcf_file);

## Resolve the coodinates
# allowing strings
# 1. empty
# 2. chromosomes
# 3. chromosome:pos1-pos2 seperated by commar

my @coor_sets;

if($coordinates ne "") {
	if($coordinates =~ m/\,/) {
		@coor_sets = split(/\,/, $coordinates);
	} else {
		push @coor_sets,$coordinates;
	}
} else {
	print "no coordinates is provided. Aborted\n";
	exit;
}

if( $vcf_file eq "") {
	print "no vcf file is provided. Aborted\n";
	exit;
}

open INPUT, $vcf_file or die "Cannot open $vcf_file.\n";
while (my $Line = <INPUT>){
	if($Line !~ m/^\#/) { # ignore commend lines
		chomp $Line;
		my @linesplit1 = split (/\t/,$Line);
		my $chr = $linesplit1[0];
		my $pos_1 = $linesplit1[1];
		my $refs = $linesplit1[3];
		my $vars = $linesplit1[4];
		my $qual = $linesplit1[5];
		my $filter_status = $linesplit1[6];
		my $info = $linesplit1[8];
		my $sample1 = $linesplit1[10];
		foreach my $chr_pos (@coor_sets) {
			if($chr_pos =~ /\:/) {
				my @chr_pos_pos = split(/\:/, $chr_pos);
				if ($chr eq $chr_pos_pos[0] ){
					if($chr_pos_pos[1] =~ m/\-/) { #chromosome:pos-pos
						my @pos_pos = split(/\-/, $chr_pos_pos[1]);
						if ( $pos_1 >= $pos_pos[0] && $pos_1 <= $pos_pos[1]) {
							print $Line."\n";
						}							
					} else { #chromosome:pos
						if ( $pos_1 == $chr_pos_pos[1] ) {
							print $Line."\n";
						}
					}
				}
			} else { #only chromosome supplied
				if($chr eq $chr_pos) {
					print $Line."\n";
				}
			}
		}
		
	} else {
		print $Line;
	}
}

exit;
	