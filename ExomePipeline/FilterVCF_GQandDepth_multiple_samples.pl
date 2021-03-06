#!/usr/bin/perl
use strict;
use warnings;

print "\n";
print "#################################################################\n";
print "#    output a new vcf file fitlered by DP and GQ on samples     #\n";
print "#    NOTE: Any line without DP or GQ info will be skipped!!     #\n";
print "#################################################################\n";

my ($vcf) = @ARGV;
my ($GQ_threshold, $DP_threshold);
$GQ_threshold = 30;
$DP_threshold = 20;
my $good_sample_threshold = 1; # number of samples failed any of the two thresholds lower than this percentage will be dropped

my $output_file = $vcf.".GQ_DP_filtered.vcf";

open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
open VCF, "$vcf" or die "Cannot open the file $vcf, check if it exists\n";
my $line_count=0;
Line: while ( my $line = <VCF> ) {
	chomp $line;
	$line_count += 1;
	if($line_count%10000 == 0) {
		print "processed $line_count lines.\n";
	}
	if ($line =~ m/^\#/) {
		print OUTPUT $line."\n";
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
			my @smaple_info;
			my $Good_Sample_count = 0;
			
			SAMPLE: for(my $i=9;$i < scalar @temp; $i++){
				my $sample = $temp[$i];
				if ($sample !~ m/^\./) {
					my @fields = split(":", $sample);
					my $GQ = $fields[$GQ_index];
					my $DP = $fields[$DP_index];
					if ($DP =~ m/^\.$/) { # HalplotypeCaller output stringe zero covered var
						push @smaple_info, "./.";
						next SAMPLE;
					}
					
					if ($DP !~ m/^\d+/) {
						#print "\n\nDP index:$DP_index    sample: $sample    DP:$DP\n";
						print "Warning: DP is not just numbers. the line is skipped:\n$line\n";
						next Line;
					}
					if ($GQ !~ m/^\d+/) {
						print "Warning: GQ is not just numbers. the line is skipped:\n$line\n";
						next Line;
					}
					if ($GQ >= $GQ_threshold && $DP >= $DP_threshold) {
						$Good_Sample_count += 1;
						push @smaple_info, $sample;
					} else {
						push @smaple_info, "./.";
					}
				} else {
					push @smaple_info, "./.";
				}
			}
			
			if($Good_Sample_count >= $good_sample_threshold) {
				my $out_string="";
				for(my $i=0;$i < 9; $i++){
					$out_string = $out_string.$temp[$i]."\t";
				}
				foreach my $sample (@smaple_info) {
					$out_string = $out_string.$sample."\t";
				}
				$out_string =~ s/\t$/\n/;
				print OUTPUT $out_string;
			}
		}
	}
}
close VCF;
close OUTPUT;

exit;
