#!/usr/bin/perl
use strict;
use warnings;

my ($input) = @ARGV;
my @input_files=split(",", $input);

my $missing_rate_threshold = 0.01;
my $MAF_theshold = 0.1;
my $Hardy_win_threshold = 0.0000001; # 1e-7

#QC on SNPs
my %failed_SNP;

print "Start to do QC on genotyping files\n";
my %missed_call;
my %A_sum;
my %B_sum;
my %number_sample_NC;
my %AA_sum;
my %AB_sum;
my %BB_sum;

print "QC first:\n";
foreach my $input_file (@input_files) {
	print "Reading file: $input_file... ";
	my @sample_columns;
	my @sample_names;
	open IN, $input_file or die "Cannot open the file to read: $input_file\n";
	while (my $line = <IN>) {
		chomp $line;
		if ($line =~ m/^Index/) {
			my @ele = split("\t", $line);
			for( my $i = 0; $i < scalar @ele; $i++) {
				if ($ele[$i] =~ m/(.+)\.GType/) {
					push @sample_columns, $i;
					push @sample_names, $1;
					if (! exists $number_sample_NC{$1}) {
						$number_sample_NC{$1} = 0;
					} else {
						print "Error: repeated sample $1, exit!\n";
						exit;
					}
				}
			}
		} else {
			my @ele = split("\t", $line);
			for (my $i=0; $i<scalar @sample_columns; $i++) {
				if (! exists $A_sum{$ele[14]}) {
					$A_sum{$ele[14]} = 0;
					$B_sum{$ele[14]} = 0;
					$missed_call{$ele[14]} = 0;
					$AA_sum{$ele[14]} = 0;
					$AB_sum{$ele[14]} = 0;
					$BB_sum{$ele[14]} = 0;
			
					if ($ele[$sample_columns[$i]] =~ m/AA/) {
						$A_sum{$ele[14]} = $A_sum{$ele[14]} + 2;
					} elsif($ele[$sample_columns[$i]] =~ m/AB/) {
						$A_sum{$ele[14]} = $A_sum{$ele[14]} + 1;
						$B_sum{$ele[14]} = $B_sum{$ele[14]} + 1;
					} elsif($ele[$sample_columns[$i]] =~ m/BB/) {
						$B_sum{$ele[14]} = $B_sum{$ele[14]} + 2;
					} elsif ($ele[$sample_columns[$i]] =~ m/NC/) {
						$missed_call{$ele[14]} = $missed_call{$ele[14]} + 1;
						$number_sample_NC{$sample_names[$i]} = $number_sample_NC{$sample_names[$i]} + 1;
					}
					
				} else {
					if ($ele[$sample_columns[$i]] =~ m/AA/) {
						$A_sum{$ele[14]} = $A_sum{$ele[14]} + 2;
						$AA_sum{$ele[14]} = $AA_sum{$ele[14]} + 1;
					} elsif($ele[$sample_columns[$i]] =~ m/AB/) {
						$A_sum{$ele[14]} = $A_sum{$ele[14]} + 1;
						$B_sum{$ele[14]} = $B_sum{$ele[14]} + 1;
						$AB_sum{$ele[14]} = $AB_sum{$ele[14]} + 1;
					} elsif($ele[$sample_columns[$i]] =~ m/BB/) {
						$B_sum{$ele[14]} = $B_sum{$ele[14]} + 2;
						$BB_sum{$ele[14]} = $BB_sum{$ele[14]} + 1;
					} elsif ($ele[$sample_columns[$i]] =~ m/NC/) {
						$missed_call{$ele[14]} = $missed_call{$ele[14]} + 1;
						$number_sample_NC{$sample_names[$i]} = $number_sample_NC{$sample_names[$i]} + 1;
					}
				}
				
				
			}
		}
	}
	close IN;
	print "Done!\n";
}

print "Calculating QC numbers.. \n";
my $number_of_samples = scalar (keys %number_sample_NC);
my $hw_failed_count=0;
foreach my $snp (keys %A_sum) {
	my $A_MAF = $A_sum{$snp}/(2*$number_of_samples);
	#print $snp."\tMAF: ".$A_MAF."\n";
	my $B_MAF = $B_sum{$snp}/(2*$number_of_samples);
	#print $snp."\tMAF: ".$B_MAF."\n";
	my $MAF;
	if ($A_MAF <= $B_MAF) {
		$MAF = $A_MAF;
	} else {
		$MAF = $B_MAF;
	}
	
	if ($MAF < $MAF_theshold) {
		$failed_SNP{$snp} = 1;
		#print $snp."\tMAF: ".$MAF."\n";
	} else {
		my $missing_rate = $missed_call{$snp}/$number_of_samples;
		if ($missing_rate > $missing_rate_threshold) {
			$failed_SNP{$snp} = 1;
			#print $snp."\tNocal rate: ".$missing_rate."\n";
		} else {
			my $hwe = snphwe($AB_sum{$snp}, $AA_sum{$snp}, $BB_sum{$snp});
			print "$AB_sum{$snp}, $AA_sum{$snp}, $BB_sum{$snp}; HWE: $hwe\n";
			if ($hwe <= $Hardy_win_threshold) {
				$failed_SNP{$snp} = 1;
				$hw_failed_count += 1;
			}
		}
	}
}
print "QC results:\n";
print "Number of samples: $number_of_samples\n";
my $number_of_failed_snp = scalar (keys %failed_SNP);
my $number_of_snp = scalar (keys %A_sum);
print "Number of SNPs in total: $number_of_snp\n";
print "Number of SNPs failed QC: $number_of_failed_snp\n";
print "Number of SNPs failed HWE: $hw_failed_count\n";
print "no call rate for each sample:\n";
foreach my $sample (sort { $number_sample_NC{$b} <=> $number_sample_NC{$a} } keys %number_sample_NC) {
	my $sample_no_call_rate = $number_sample_NC{$sample}/$number_of_snp;
	print $sample.": ".$sample_no_call_rate."\n";
}

print "\nPress Y to continue, other key to exit:";
my $continue =<STDIN>;
chomp $continue;
print "Continue word caught:$continue\n";
if ($continue !~ m/^Y/i) {
	exit;
}

foreach my $input_file (@input_files) {
	print "Reading file: $input_file\n";
	my $outputfile = $input_file."_processed_for_MatrixEQTL.txt";
	
	my @out_columns=('9','10','12','13','15'); #1-based column numbers
	my @out_column_names=('Chr_GRCh38','Pos_GRCh38','Pos_GRCh38','Alt_dbsnp141','Name');
	
	open IN, $input_file or die "Cannot open the file to read: $input_file\n";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ m/^Index/) {
			my @ele = split("\t", $line);
			for( my $i = 0; $i < scalar @ele; $i++) {
				if ($ele[$i] =~ m/(.+)\.GType/) {
					push @out_columns, $i+1;
					push @out_column_names, $1;
				}
			}
			last;
		} else {
			print "Error: No header found. Exit.\n";
			exit;
		}
	}
	close IN;
	
	#Start to process
	open IN, $input_file or die "Cannot open the file to read: $input_file\n";
	open OUT, ">$outputfile" or die "Cannot open the file to ouput: $outputfile\n";
	my $out_string="";
	foreach my $name (@out_column_names) {
		$out_string=$out_string.$name."\t";
	}
	$out_string =~ s/\t$//;
	print OUT $out_string."\n";
	
	while (my $line = <IN>){
		if ($line =~ m/^Index/) {
			next;
		} else {
			chomp $line;
			$out_string="";
			my @ele = split("\t", $line);
			
			if (exists $failed_SNP{$ele[14]}) {
				next;
			}
			
			foreach my $index (@out_columns) {
				if ($ele[$index-1] =~ m/AA/) {
					$out_string=$out_string."0"."\t";
				} elsif ($ele[$index-1] =~ m/AB/) {
					$out_string=$out_string."1"."\t";
				} elsif ($ele[$index-1] =~ m/BB/) {
					$out_string=$out_string."2"."\t";
				} elsif ($ele[$index-1] =~ m/NC/) {
					$out_string=$out_string."NA"."\t";
				} else {
					$out_string=$out_string.$ele[$index-1]."\t";
				}
			}
			
			$out_string =~ s/\t$//;
			print OUT $out_string."\n";
		}
	}
	close IN;
	close OUT;
	print "Done!\n";
}
exit;

sub snphwe {
    my $obs_hets = shift;
    my $obs_hom1 = shift;
    my $obs_hom2 = shift;

    if($obs_hom1 < 0 || $obs_hom2 < 0 || $obs_hets <0) {
	return(-1);
    }

    # rare homozygotes
    my $obs_homr;

    # common homozygotes
    my $obs_homc;
    if($obs_hom1 < $obs_hom2) {
	$obs_homr = $obs_hom1;
	$obs_homc = $obs_hom2;
    } else {
	$obs_homr = $obs_hom2;
	$obs_homc = $obs_hom1;
    }

    # number of rare allele copies
    my $rare_copies = 2 * $obs_homr + $obs_hets;

    # total number of genotypes
    my $genotypes = $obs_homr + $obs_homc + $obs_hets;

    if($genotypes <= 0) {
	return(-1);
    }
    
    # Initialize probability array
    my @het_probs;
    for(my $i=0; $i<=$rare_copies; $i++) {
	$het_probs[$i] = 0.0;
    }

    # start at midpoint
    my $mid = int($rare_copies * (2 * $genotypes - $rare_copies) / (2 * $genotypes));

    # check to ensure that midpoint and rare alleles have same parity
    if(($rare_copies & 1) ^ ($mid & 1)) {
	$mid++;
    }
    
    my $curr_hets = $mid;
    my $curr_homr = ($rare_copies - $mid) / 2;
    my $curr_homc = $genotypes - $curr_hets - $curr_homr;

    $het_probs[$mid] = 1.0;
    my $sum = $het_probs[$mid];
    for($curr_hets = $mid; $curr_hets > 1; $curr_hets -= 2) {
	$het_probs[$curr_hets - 2] = $het_probs[$curr_hets] * $curr_hets * ($curr_hets - 1.0) / (4.0 * ($curr_homr + 1.0) * ($curr_homc + 1.0));
	$sum += $het_probs[$curr_hets - 2];

	# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
	$curr_homr++;
	$curr_homc++;
    }

    $curr_hets = $mid;
    $curr_homr = ($rare_copies - $mid) / 2;
    $curr_homc = $genotypes - $curr_hets - $curr_homr;
    for($curr_hets = $mid; $curr_hets <= $rare_copies - 2; $curr_hets += 2) {
	$het_probs[$curr_hets + 2] = $het_probs[$curr_hets] * 4.0 * $curr_homr * $curr_homc / (($curr_hets + 2.0) * ($curr_hets + 1.0));
	$sum += $het_probs[$curr_hets + 2];
	
	# add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
	$curr_homr--;
	$curr_homc--;
    }

    for(my $i=0; $i<=$rare_copies; $i++) {
	$het_probs[$i] /= $sum;
    }

    # alternate p-value calculation for p_hi/p_lo
#    my $p_hi = $het_probs[$obs_hets];
#    for(my $i=$obs_hets+1; $i<=$rare_copies; $i++) {
#	$p_hi += $het_probs[$i];
#    }
#    
#    my $p_lo = $het_probs[$obs_hets];
#    for(my $i=$obs_hets-1; $i>=0; $i--) {
#	$p_lo += $het_probs[$i];
#    }
#
#    my $p_hi_lo;
#    if($p_hi < $p_lo) {
#	$p_hi_lo = 2 * $p_hi;
#    } else {
#	$p_hi_lo = 2 * $p_lo;
#    }

    # Initialise P-value 
    my $p_hwe = 0.0;

    # P-value calculation for p_hwe
    for(my $i = 0; $i <= $rare_copies; $i++) {
	if($het_probs[$i] > $het_probs[$obs_hets]) {
	    next;
	}
	$p_hwe += $het_probs[$i];
    }
    
    if($p_hwe > 1) {
	$p_hwe = 1.0;
    }

    return($p_hwe);
}


