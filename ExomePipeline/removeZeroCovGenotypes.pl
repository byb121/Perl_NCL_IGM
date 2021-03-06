#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

print "chr preceding chrosome names will be removed as well. \nVariants with no AD tag will be removed.\nLines with no vailed genotype calls will be removed!\n\n";

my $vcf_in;
my $out_vcf;
my $help;

usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, "vcf=s"=>\$vcf_in, "out=s"=>\$out_vcf) || defined $help );

open VCF, "$vcf_in" or die "Cannot open the file $vcf_in or the file does not exist!\n";
open OUT, ">$out_vcf" or die "Cannot write to the file, please  check the path $out_vcf!\n";
LINE: while (my $line = <VCF>){
	if ($line =~ m/^\#/) {
		print OUT $line;
	} else {
		chomp $line;
		$line =~ s/^chr//;
		my @ele = split("\t", $line);
		# find the index of DP tag
		my $format = $ele[8];
		my @format_fields = split(":", $format);
		my $AD_index;
		for (my $i=0;$i<scalar(@format_fields);$i++){
			if ($format_fields[$i] =~ m/AD/) {
				$AD_index = $i;
			}
		}
		#change 0-covered variant to ./. if it is a single sample VCF, discard the 0-covered line
		
		if (defined $AD_index) {
			my @corrected;
			my $newline="";
			for(my $i=9; $i<scalar @ele;$i++){
				my $sample = $ele[$i];
				my @fields = split(":", $sample);
				if ($sample !~ m/\.\/\./) {
					my $AD = $fields[$AD_index];
					if ($AD ne ".") {
						my @covs = split(",", $AD);
						my $all_0 = "Yes";
						COV: foreach my $cov (@covs) {
							if ($cov != 0) {
								$all_0 = "No";
								last COV;
							}
						}
						if ($all_0 eq "Yes") {
							#print "Found one error $sample in line: $line\n";
							push @corrected, './.';
						} else {
							push @corrected, $sample;
						}
					} else {
						push @corrected, './.';
					}
				} else {
					push @corrected, $sample;
				}
			}
			my $discard="yes";
			foreach my $corrected_call (@corrected) {
				if ($corrected_call ne "./.") {
					$discard="no";
					last;
				}
			}
			if ($discard eq "no") {
				for(my $i=0; $i<9;$i++) {
					$newline = $newline.$ele[$i]."\t";
				}
				foreach my $word (@corrected) {
					$newline = $newline.$word."\t";
				}
				$newline =~ s/\t$//;
				print OUT $newline."\n";
			} else {
				print "Discard All 0 covered: $line\n";
			}
		} else {
			#print OUT $line."\n";
			print "no AD tag on the line: $line"."\n";
		}
	}
}
close VCF;
close OUT;
exit;



sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: VCF_2_annotated_excel_20131120.pl \n";
    print "--vcf input vcf file (of a single sample or a family;\n";
    print "--out output vcf file with 0-cov var set to \'./.\'.\n\n";
	return(1);
}