#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

print "\n\n";
print "#################################################################\n";
print "take a input vcf file, and will only read the first sample column\n"; 
print "and output a new fitlered vcf file.\n";
print "#################################################################\n";

#AO
#DP
#QR  cannot be used
#QA
#MQM mean mapping qualtiy of observed alleles
#SAF
#SAR
#coverage of the reference = DP - sum(AO)

my $in_vcf_string;
my $output_file;

my $help;
usage() if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, "in=s"=>\$in_vcf_string, 'out=s' => \$output_file) || defined $help );

my @vcfs = split(",", $in_vcf_string);			

my @header;
my $sampleLine="";
my %samples;
my $vcf1 = $vcfs[0];
my $format;

#Taking headers
open VCF, "$vcf1" or die "Cannot open the file $vcf1, check if it exists\n";
while ( my $line = <VCF> ) {
	chomp $line;
	if ($line =~ m/^\#\#/) {
		push @header, $line."\n";
	} elsif ($line =~ m/^\#CHR/) {
		my @temp = split("\t", $line);
		for (my $i=0;$i < 9;$i++) {
			$sampleLine = $sampleLine."\t".$temp[$i];
		}
		$sampleLine =~ s/^\t//;
	} else {
		last;
	}
}
close VCF;

my %all_variant;
my @all_samples;

foreach my $vcf (@vcfs) {
	print "open file $vcf\n";
	open VCF, "$vcf" or die "Cannot open the file $vcf, check if it exists\n";
	my @vcf_samples;
	while ( my $line = <VCF> ) {
		chomp $line;
		if ($line =~ m/^\#\#/) {
			next;
		} elsif ($line =~ m/^\#CHR/) {
			my @temp = split("\t", $line);
			for (my $i=9;$i < scalar @temp;$i++) {
				if (! exists $samples{$temp[$i]}) {
					print "Found sample $temp[$i]\n";
					$sampleLine = $sampleLine."\t".$temp[$i];
					$samples{$temp[$i]} = 1;
					push @all_samples, $temp[$i];
					push @vcf_samples, $temp[$i];
				} else {
					print "Duplicated samples found in the files $temp[$i]\nexit!\n";
					exit;
				}
			}
		} else {
			my @temp = split("\t", $line);
			#print "read in $line;\n";
			for (my $i=0; $i < scalar @vcf_samples; $i++) {
				if (! $format) {$format = $temp[8];}
				my $sample = $vcf_samples[$i];
				my $sample_line = $temp[$i+9];
				#print "sample line is $sample_line\n";
				if ( ! exists $all_variant{$temp[0]}{$temp[1]}{$temp[3]}{$sample}{$temp[4]}) {
					$all_variant{$temp[0]}{$temp[1]}{$temp[3]}{$sample}{$temp[4]} = $sample_line;
					print "read in $temp[0] $temp[1] $temp[3] $sample $temp[4] $sample_line\n";
				} else {
					print "It is wrong, duplicated line for one sample?\nexit!\n";
					exit;
				}
			}
		}
	}
	close VCF;
}			

my @output;

foreach my $chr (keys %all_variant) {
	foreach my $pos ( sort {$a <=> $b} keys %{$all_variant{$chr}}) {
		foreach my $ref (keys %{$all_variant{$chr}{$pos}}) {
			my %alts_hash;
			my @alts_array;
			my $alt_index=1;
			my $output_string="$chr\t$pos\t.\t$ref\t";
			foreach my $sample (sort {$a cmp $b} keys %{$all_variant{$chr}{$pos}{$ref}}) {
				foreach my $v (keys %{$all_variant{$chr}{$pos}{$ref}{$sample}}) {
					my @v_s = split(",", $v);
					for(my $i = 0; $i < scalar @v_s;$i++) {
						if (! exists $alts_hash{$v}) {
							$alts_hash{$v} = $alt_index;
							push @alts_array, $v;
							$output_string=$output_string."$v,";
							$alt_index += 1;
						}
					}
				}
			}
			$output_string =~ s/\,$//;
			$output_string = $output_string."\t99\tPASS\t.\t$format";
			
			foreach my $sample (@all_samples) {
				if(exists $all_variant{$chr}{$pos}{$ref}{$sample}) {
					my $new_ht_string="";
					foreach my $v (keys %{$all_variant{$chr}{$pos}{$ref}{$sample}}) {
						my @v_s = split(",", $v);
						my @sample_info_temp = split(":", $all_variant{$chr}{$pos}{$ref}{$sample}{$v});
						my @HT_s=split("/", $sample_info_temp[0]);
						foreach my $het_number (@HT_s) {
							if ($het_number =~ m/^0$/) {
								$new_ht_string=$new_ht_string."/".$het_number;
							} else {
								$new_ht_string=$new_ht_string."/".$alts_hash{$v_s[$het_number-1]};
							}
						}
						$new_ht_string =~ s/^\///;
						$new_ht_string=$new_ht_string.":".$sample_info_temp[1].":".$sample_info_temp[2];
					}
					$output_string = $output_string."\t".$new_ht_string;
				} else {
					$output_string = $output_string."\t"."./.";
				}
			}
			push @output, $output_string."\n";
		}
	}
}

print "Output merge variants to $output_file\n";
open OUTPUT, ">$output_file" or die "Cannot open file $$output_file.\n";
print OUTPUT @header;
print OUTPUT $sampleLine."\n";
print OUTPUT @output;
close OUTPUT;
print "Done!\n";		
			
exit;		
			
			
sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "\nusage: VCF_2_annotated_excel_20131120.pl \n";
    print "--in input vcf file (of a single sample or a family), multiple files are seperated by \",\";\n";
    print "--out output file name.\n";
	return(1);
}