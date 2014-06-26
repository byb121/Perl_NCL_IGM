#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Statistics::TTest;

my ($folder) = @ARGV;

$folder =~ s/\/$//;
my @files;
my $dir = cwd();

my %hash;

my $depth_filter = 10; #locus with seqencing depth lower than the value will be dropped.
my $alt_quality_filter = 13; # Phred-scaled p-value of the alternative base calling is wrong
my $gq_quality_filter = 13; # Phred-scaled p-value of the genotype calling is wrong

opendir (DIR, $folder) or die "Cann't find the directory";
while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./);
	if ($file =~ m/\.flt.vcf$/){ ## This line controls the name pattern of files that are read in
		my $input_file = $folder."/".$file;
		open INPUT, $input_file or die "cannot open file $input_file.\n";
		
		print "reading ".$input_file."\n";
		while (my $line =<INPUT>) {
			chomp $line;
			if( $line =~ m/^chr/ ){
				my @words = split("\t", $line);
				my $key = $words[0]."_".$words[1];
				
				##filters: include Sequencing depth filter, ALT quaility filter and GenoType call quality filter
				my $depth = $words[5];
				my $alt_quality;
				my @vcf_info = split (";", $words[7]);
				foreach my $tag(@vcf_info) {
					if ($tag =~ m/^DP\=/) {
						$tag =~ s/^DP\=//;
						$alt_quality = $tag;
					}
				}
				my @vcf_sample = split (":", $words[9]);
				my $gq_quality = $vcf_sample[2];
				
				if($depth >= $depth_filter && $alt_quality >= $alt_quality_filter && $gq_quality >= $gq_quality_filter) {
					if (exists $hash{$key}) {
						$hash{$key} += 1;
					} else {
						$hash{$key} = 1;
					}	
				}
			} 
		}
		close INPUT;
	}
}
closedir(DIR);

print "###########################################\n";
print "########### output results ################\n";
print "###########################################\n";

my $output_file = $dir."/"."CommonSNPs_filtered.txt";
open OUTPUT, ">$output_file" or die "cannot open file $output_file.\n";
for my $key ( keys %hash ) {
	if ($hash{$key} == 16) {
		my @words = split("_", $key);
		print OUTPUT $words[0]."\t".$words[1]."\n";
	}
}
close OUTPUT;

exit;




