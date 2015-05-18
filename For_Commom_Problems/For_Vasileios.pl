#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


my $vcf;
my $tab;
my $out;

my $help;
usage() if ( @ARGV < 3 || ! GetOptions('help|?' => \$help, "vcf=s"=>\$vcf, "tab=s"=>\$tab, 
	'out=s' => \$out) || defined $help );

my %AF;
open VCF, "$vcf" or die "Can not open the fasta index file: $vcf\n";
while (my $line = <VCF> ) {
	if ($line =~ m/^\#/) {
		next;
	} else {
		chomp $line;
		my @eles = split(/\t/, $line);
		my $chr = $eles[0];
		my $pos = $eles[1];
		my $ref = $eles[3];
		my $var = $eles[4];
		my @eles2 = split(";", $eles[7]);
		my $found="no";
		
		foreach my $field (@eles2){
			if ($field =~ /^AF/) {
				$found = "yes";
				my $temp = $field; $temp =~ s/^AF=//;
				if (! exists $AF{$chr}{$pos}{$ref}{$var}) {
					$AF{$chr}{$pos}{$ref}{$var} = $temp;
				} else {
					print "Duplicated entry found $line\n Abort!\n";
					exit;
				}
			}
		}
		
		if ($found eq "no") {
			print "VCF does not have AF field on $line\n Abort!\n";
			exit;
		}
	}
}
close VCF;

open TAB, "$tab" or die "Can not open the fasta index file: $tab\n";
open OUT, ">$out" or die "Can not open the fasta index file: $out\n";
while (my $line = <TAB> ) {
	chomp $line;
	my @eles = split ("\t", $line);
	my $chr = $eles[2];
	my $pos = $eles[3];
	my $ref = $eles[5]; if ($ref eq "0") {$ref="N";}
	my $var = $eles[6];
	
	my $out_string="";
	for (my $i=0;$i<scalar @eles;$i++) {
		if ($i == 7) {
			if (exists $AF{$chr}{$pos}{$ref}{$var}) {
				$out_string = $out_string.$AF{$chr}{$pos}{$ref}{$var}."\t";
			} else {
				print "No AF entry found in VCF for tab line: $line\n";
				$out_string = $out_string.$eles[$i]."\t";
			}
		} else {
			$out_string = $out_string.$eles[$i]."\t";
		}
	}
	$out_string =~ s/\t$/\n/;
	print OUT $out_string;
}
close OUT;
close TAB;

print "Done!";

exit;

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "--vcf the vcf file has AF field;\n";
    print "--tab tab delimited text file;\n"; 
    print "--HC output.\n";
    return(1);
}