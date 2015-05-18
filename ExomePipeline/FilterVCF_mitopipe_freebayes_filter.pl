#!/usr/bin/perl
use strict;
use warnings;

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

#if ($AOs[$i] >= $AO_threshold && $DPs[$i] >= $DP_threshold && $QAs[$i] >= $QA_threshold*$AOs[$i]*0.9 
#			&& $MQMs[$i] >= $MQM_threshold && $SARs[$i] >= $SAR_threshold*$AOs[$i] && $SAFs[$i] >= $SAF_threshold*$AOs[$i]) {	

my $AO_threshold = 1000; # alternative allele depth threshold
my $DP_threshold = 1000; # total depth threshold
my $QA_threshold = 30; # at least 90% of bases should have base quality scores over the threshold.(not doing that exactly:))
my $MQM_threshold = 50; # mean mapping quality threshold;
my $SAR_threshold = 0.2; # at least this fraction of the alternative obervations are on reverse strand
my $SAF_threshold = 0.2; # at least this fraction of the alternative obervations are on forward strand
my $ref_depth_threshold = 0.1; # at least this fraction of the total depth is reference base to consider reference base presents 

my ($vcf) = @ARGV;

my @output;
open VCF, "$vcf" or die "Cannot open the file $vcf, check if it exists\n";
while ( my $line = <VCF> ) {
	chomp $line;
	if ($line =~ m/^\#\#/) {
		push @output, $line."\n";
	} elsif ($line =~ m/^\#CHR/) {
		push @output, '##FORMAT=<ID=HT,Number=1,Type=String,Description="Heterzygosity type. 0-ref 1-alt 2-second alt and so on, sep with "/" if multiple alleles>'."\n";
		push @output, '##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Count for all reported alleles">'."\n";
		push @output, $line."\n";
	} else {
		my ($AO, $DP, $QA, $MQM, $SAF, $SAR, $RP);
		my @temp = split("\t", $line);
		if ( $temp[4] =~ m/\./) {
			next;
		}
		my $info=$temp[7];
		if ($info =~ m/(\;|^)AO\=(.+?)\;/ || $info =~ m/(\;|^)AO\=(.+?)$/) {
			$AO=$2;
		} else {
			print "AO tag does not exsit on row $line\nexit!\n";
			exit;
		}
		if ($info =~ m/(\;|^)DP\=(.+?)\;/ || $info =~ m/(\;|^)DP\=(.+?)$/) {
			$DP=$2;
		} else {
			print "DP tag does not exsit on row $line\nexit!\n";
			exit;
		}
		if ($info =~ m/(\;|^)QA\=(.+?)\;/ || $info =~ m/(\;|^)QA\=(.+?)$/) {
			$QA=$2;
		} else {
			print "QA not does not exsit on row $line\nexit!\n";
			exit;
		}
		if ($info =~ m/(\;|^)MQM\=(.+?)\;/ || $info =~ m/(\;|^)MQM\=(.+?)$/) {
			$MQM=$2;
		} else {
			print "MQM tag does not exsit on row $line\nexit!\n";
			exit;
		}
		if ($info =~ m/(\;|^)SAF\=(.+?)\;/ || $info =~ m/(\;|^)SAF\=(.+?)$/) {
			$SAF=$2;
		} else {
			print "SAF tag does not exsit on row $line\nexit!\n";
			exit;
		}
		if ($info =~ m/(\;|^)SAR\=(.+?)\;/ || $info =~ m/(\;|^)SAR\=(.+?)$/) {
			$SAR=$2;
		} else {
			print "SAR tag does not exsit on row $line\nexit!\n";
			exit;
		}
		
		my $ref = $temp[3];
		my $alt = $temp[4];
		
		my @alts = split(",", $alt);
		my @AOs = split(",", $AO);
		my @DPs = split(",", $DP);
		my $RO = $DPs[0];
		foreach my $temp_AO (@AOs) {
			$RO = $RO-$temp_AO;
		}
		my @QAs = split(",", $QA);
		my @MQMs = split(",", $MQM);
		my @SARs = split(",", $SAR);
		my @SAFs = split(",", $SAF);
		my @kept_alt_index;
		#print "The line info is: $info\n";
		for (my $i=0; $i < scalar @alts; $i++) {
			if ($AOs[$i] >= $AO_threshold && $DPs[$i] >= $DP_threshold && $QAs[$i] >= $QA_threshold*$AOs[$i]*0.9 
			&& $MQMs[$i] >= $MQM_threshold && $SARs[$i] >= $SAR_threshold*$AOs[$i] && $SAFs[$i] >= $SAF_threshold*$AOs[$i]) {
				print "Variant $ref/$alts[$i] passed the filters\n";
				push @kept_alt_index, $i;
			} #else {
				#if ($AOs[$i] < $AO_threshold) {
				#	print "AO Low. ";
				#}
				#if ($DPs[$i] < $DP_threshold) {
			#		print "DP Low. ";
			#	}
			#	if ($QAs[$i] < $QA_threshold*$AOs[$i]*0.9) {
			#		print "QA Low. ";
			#	}
			#	if ($MQMs[$i] < $MQM_threshold) {
			#		print "MQM Low. ";
			#	}
			#	if ($SARs[$i] < $SAR_threshold*$AOs[$i]) {
			#		print "SAR Low. ";
			#	}
			#	if ($SAFs[$i] < $SAF_threshold*$AOs[$i]) {
			#		print "SAF Low. ";
			#	}
			#	print "\n";
			#}
		}
		if (scalar @kept_alt_index < 1) {
			next;
		}
		
		my %output;
		foreach my $v_index (@kept_alt_index) {
			my $new_entry_ref = ProcessFreeBayesVCF($temp[1], $ref, $alts[$v_index]);
			my @new_entry = @$new_entry_ref;
			if (! exists $output{$temp[0]}{$new_entry[0]}{$new_entry[1]}{$new_entry[2]} ) {
				my @temp1 = ($AOs[$v_index], $DPs[$v_index], $QAs[$v_index], $MQMs[$v_index], $SARs[$v_index], $SAFs[$v_index]);
				$output{$temp[0]}{$new_entry[0]}{$new_entry[1]}{$new_entry[2]} = \@temp1; 
			} else {
				print "Bug a~ Da bug!\nexit!\n";
				exit;
			}
		}
		foreach my $chr (keys %output) {
			print "chr: $chr ";
			foreach my $pos ( sort {$a<=>$b} keys %{$output{$chr}}) {
				print "pos: $pos ";
				foreach my $ref (keys %{$output{$chr}{$pos}}) {
					print "ref: $ref ";
					my $output_string = "$chr\t$pos\t\.\t$ref\t";
					my $var_string = "";
					my $AO_string="";
					my $DP_string="";
					my $QA_string="";
					my $MQM_string="";
					my $SAR_string="";
					my $SAF_string="";
					my $sample_HT_string = "";
					my $sample_AD_string = ""; # heterozygosity type 0/1/2 would mean ref and first and second alternative alleles.
					if ($RO >= $ref_depth_threshold*$DPs[0]) {
						$sample_HT_string = "0";
						$sample_AD_string = "$RO";
					}
					my $count = 1;
					foreach my $var (keys %{$output{$chr}{$pos}{$ref}}) {
						print "var: $var ";
						my $temp_array_ref = $output{$chr}{$pos}{$ref}{$var};
						my @temp_array = @$temp_array_ref;
						$var_string = $var_string.",".$var;
						$AO_string = $AO_string.",".$temp_array[0];
						$DP_string = $DP_string.",".$temp_array[1];
						$QA_string = $QA_string.",".$temp_array[2];
						$MQM_string = $MQM_string.",".$temp_array[3];
						$SAR_string = $SAR_string.",".$temp_array[4];
						$SAF_string = $SAF_string.",".$temp_array[5];
						$sample_HT_string = $sample_HT_string."/".$count;
						$count += 1;
					}
					print "\n";
					$var_string =~ s/^\,//;
					$AO_string =~ s/^\,//;
					$DP_string =~ s/^\,//;
					$QA_string =~ s/^\,//;
					$MQM_string =~ s/^\,//;
					$SAR_string =~ s/^\,//;
					$SAF_string =~ s/^\,//;
					$sample_HT_string =~ s/^\///;
					$sample_AD_string = $sample_AD_string.",".$AO_string;
					$sample_AD_string =~ s/^\,//;
					my $sample_DP_string = $DPs[0];
					$output_string = $output_string."$var_string\t$temp[5]\tPASS\tAO=$AO_string;DP=$DP_string;QA=$QA_string;".
					"MQM=$MQM_string;SAR=$SAR_string;SAF=$SAF_string\tHT:AD:DP\t$sample_HT_string:$sample_AD_string:$sample_DP_string\n";
					push @output, $output_string;
				}
			}
		}
	}
}
my $output_file = $vcf.".filtered.vcf";
open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
print OUTPUT @output;
close OUTPUT;

exit;


sub ProcessFreeBayesVCF {
	my ($pos, $ref, $alt) = @_;
	print "to process $pos, $ref, $alt\n";
	if ($alt eq $ref) {
		print "How wrong freebayes can go!! Pos: $pos alt: $alt ref: $ref\nexit!\n";
		exit;
	}
	
	my @return_array;
	if ($ref =~ m/^[atgcATGC]$/) { # ref is a single base 	
		@return_array = ($pos, $ref, $alt);
	} else { ##deletion or whatever. Freebayes ....
		my @alt_chars = split("", $alt);
		my @ref_chars = split("", $ref);
		for (my $i=0; $i < scalar @alt_chars; $i++) {
			if ($alt_chars[$i] eq $ref_chars[$i] && length $alt > 1 && length $ref > 1) {
				$alt =~ s/^.//;
				$ref =~ s/^.//;
				$pos += 1;
			} else {
				last;
			}
		}
		if (length $alt == length $ref) {
			for (my $i=scalar @alt_chars - 1; $i >=0 ; $i--) {
				if ($alt_chars[$i] eq $ref_chars[$i]) {
					$alt =~ s/.$//;
					$ref =~ s/.$//;
				} else {
					last;
				}
			}
		}
		@return_array = ($pos, $ref, $alt);
	}
	print "after process: $pos, $ref, $alt\n";
	return \@return_array;
}








