#!/usr/bin/perl
use strict;
use warnings;

print "\n";
print "#################################################################\n";
print "#    output a new vcf file fitlered by DP and GQ on samples     #\n";
print "#################################################################\n";

my ($vcf) = @ARGV;
my ($GQ_threshold, $DP_threshold, $sample_percent_threshold);
$GQ_threshold = 20;
$DP_threshold = 10;
$sample_percent_threshold = 100; # fraction of samples failed any of the two thresholds lower than this percentage will be dropped

my @output;
open VCF, "$vcf" or die "Cannot open the file $vcf, check if it exists\n";
my $line_count=0;
Line: while ( my $line = <VCF> ) {
	chomp $line;
	$line_count += 1;
	if($line_count%10000 == 0) {
		print "processed $line_count lines.\n";
	}
	if ($line =~ m/^\#/) {
		push @output, $line."\n";
	} else {
		my @temp = split("\t", $line);
		my $format = $temp[8];
		my @format_fields = split(":", $format);
		my $GQ_index;
		my $DP_index;
		#my $GT_index;
		#my $PL_index;
		for (my $i=0;$i<scalar(@format_fields);$i++){
			if ($format_fields[$i] =~ m/GQ/) {
				$GQ_index = $i;
				#print $line."$GQ_index\n";
			}
			if ($format_fields[$i] =~ m/DP/) {
				$DP_index = $i;
				#print $line."$GT_index\n";
			}
			#if ($format_fields[$i] =~ m/AD/) {
			#	$AD_index = $i;
			#	#print $line."$AD_index\n";
			#}
			#if ($format_fields[$i] =~ m/PL/) {
			#	$PL_index = $i;
			#	#print $line."$AD_index\n";
			#}
		}
		
		if (defined $GQ_index && defined $DP_index) {
			my @all_GQ;
			my @all_DP;
			
			for(my $i=9;$i < scalar @temp; $i++){
				my $sample = $temp[$i];
				#print $sample."\n";
				if ($sample !~ m/\.\/\./) {
					my @fields = split(":", $sample);
					my $GQ = $fields[$GQ_index];
					my $DP = $fields[$DP_index];
					if ($DP !~ m/^\d+$/) {
						next Line;
					}
					if ($GQ !~ m/^\d+$/) {
						next Line;
					}
					#print "GQ: $GQ, DP: $DP\n";
					push @all_GQ, $GQ;
					push @all_DP, $DP;
				} else {
					#print "GQ: 0, DP: 0\n";
					push @all_GQ, 0;
					push @all_DP, 0;
				}
			}
			my $Good_GQ_count = 0;
			my $Good_DP_count = 0;
			foreach my $GQ (@all_GQ) {
				if($GQ >= $GQ_threshold) {
					$Good_GQ_count += 1;
				}
			}
			foreach my $DP (@all_DP) {
				if($DP >= $DP_threshold) {
					$Good_DP_count += 1;
				}
			}
			if($Good_DP_count <= $Good_GQ_count) {
				if($Good_DP_count/(scalar @temp - 9) >= $sample_percent_threshold/100) {
					push @output, $line."\n";
				}
			} else {
				if($Good_GQ_count/(scalar @temp - 9) >= $sample_percent_threshold/100) {
					push @output, $line."\n";
				}
			}
		}
	}
}
my $output_file = $vcf.".GQ_DP_filtered.vcf";
print "Writing output...\n";
open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
print OUTPUT @output;
close OUTPUT;
print "Done!\n";

exit;
