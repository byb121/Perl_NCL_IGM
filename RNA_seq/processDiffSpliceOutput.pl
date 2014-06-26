#!/usr/bin/perl
use strict;
use warnings;

my @files = @ARGV;

my $output_file = $files[0].".compared.txt";
my %all_sig_hits_hash;

my $genefile="/users/a5907529/lustre/scripts/AfterDindel/Ensembl_Genes_R67.txt";

my $en_genes_ref = GetGeneCoords($genefile);
my %en_genes = %$en_genes_ref;

foreach my $input_file (@files) {
	print "Process the file $input_file\n";
	open INPUT, $input_file or die "cannot open file $input_file.\n";
	while (my $line = <INPUT>) {
		if ($line =~ m/^chr.+yes/) {
			chomp $line;
			my @words = split("\t", $line);
			my $chr = $words[0];
			my $start = $words[1];
			my $end = $words[2];
			my $cate = $words[3];
			my $stat_diff = $words[4];
			my $sqrt_JSD = $words[5];
			my $cov_1 = $words[6];
			my $cov_2 = $words[7];
			
			LOOP_I: for( my $i=0;$i<=100;$i++) {
				if (exists $all_sig_hits_hash{$chr}{$start+$i}) {
					for( my $j=0;$j<=100;$j++) {
						if (exists $all_sig_hits_hash{$chr}{$start+$i}{$end+$j} || exists $all_sig_hits_hash{$chr}{$start+$i}{$end-$j} ) {
							if (exists $all_sig_hits_hash{$chr}{$start+$i}{$end+$j}) {
								my @temp = split("\t", $all_sig_hits_hash{$chr}{$start+$i}{$end+$j});
								$stat_diff = $stat_diff + $temp[4];
								$sqrt_JSD = $sqrt_JSD + $temp[5];
								$cov_1 = $cov_1 + $temp[6];
								$cov_2 = $cov_2 + $temp[7];
								$temp[8] += 1;
								 $all_sig_hits_hash{$chr}{$start+$i}{$end+$j} = "$chr\t$start\t$end\t$cate\t$stat_diff\t$sqrt_JSD\t$cov_1\t$cov_2\t".$temp[8];	
								print "Found similar entry: $chr: $start+$i $end+$j  oringinal $start $end\n";
							} else {
								my @temp = split("\t", $all_sig_hits_hash{$chr}{$start+$i}{$end-$j});
								$stat_diff = $stat_diff + $temp[4];
								$sqrt_JSD = $sqrt_JSD + $temp[5];
								$cov_1 = $cov_1 + $temp[6];
								$cov_2 = $cov_2 + $temp[7];
								$temp[8] += 1;
								 $all_sig_hits_hash{$chr}{$start+$i}{$end-$j} = "$chr\t$start\t$end\t$cate\t$stat_diff\t$sqrt_JSD\t$cov_1\t$cov_2\t".$temp[8];	
								print "Found similar entry: $chr: $start+$i $end-$j  oringinal $start $end\n";
							} 
							last LOOP_I;
						}
					}
				} elsif (exists $all_sig_hits_hash{$chr}{$start-$i}) {
					for( my $j=0;$j<=100;$j++) {
						if (exists $all_sig_hits_hash{$chr}{$start-$i}{$end+$j} || exists $all_sig_hits_hash{$chr}{$start-$i}{$end-$j} ) {
							if (exists $all_sig_hits_hash{$chr}{$start-$i}{$end+$j}) {
								my @temp = split("\t", $all_sig_hits_hash{$chr}{$start-$i}{$end+$j});
								$stat_diff = $stat_diff + $temp[4];
								$sqrt_JSD = $sqrt_JSD + $temp[5];
								$cov_1 = $cov_1 + $temp[6];
								$cov_2 = $cov_2 + $temp[7];
								$temp[8] += 1;
								 $all_sig_hits_hash{$chr}{$start-$i}{$end+$j} = "$chr\t$start\t$end\t$cate\t$stat_diff\t$sqrt_JSD\t$cov_1\t$cov_2\t".$temp[8];	
								print "Found similar entry: $chr: $start-$i $end+$j  oringinal $start $end\n";
							} else {
								my @temp = split("\t", $all_sig_hits_hash{$chr}{$start-$i}{$end-$j});
								$stat_diff = $stat_diff + $temp[4];
								$sqrt_JSD = $sqrt_JSD + $temp[5];
								$cov_1 = $cov_1 + $temp[6];
								$cov_2 = $cov_2 + $temp[7];
								$temp[8] += 1;
								 $all_sig_hits_hash{$chr}{$start-$i}{$end-$j} = "$chr\t$start\t$end\t$cate\t$stat_diff\t$sqrt_JSD\t$cov_1\t$cov_2\t".$temp[8];	
								print "Found similar entry: $chr: $start-$i $end-$j  oringinal $start $end\n";
							}
							last LOOP_I;
						}
					}
				} elsif ($i == 100) {
					#print "new record\n";
					$all_sig_hits_hash{$chr}{$start}{$end} = "$chr\t$start\t$end\t$cate\t$stat_diff\t$sqrt_JSD\t$cov_1\t$cov_2\t1";
				} else {
					next;
				}
			}
			
		} else {
			next;
		}
	}
	close INPUT;
}

open OUTPUT, ">$output_file" or die "Cannot open the file $output_file to output\n";
foreach my $chr ( keys %all_sig_hits_hash) {
	foreach my $start ( sort {$a<=>$b} keys %{$all_sig_hits_hash{$chr}}) {
		foreach my $end (keys %{$all_sig_hits_hash{$chr}{$start}}) {
			my @temp = split("\t", $all_sig_hits_hash{$chr}{$start}{$end});
			if ($temp[8] == 4) {
				my $chr = $temp[0];
				my $start = $temp[1];
				my $end = $temp[2];
				my $cate = $temp[3];
				my $stat_diff = $temp[4]/$temp[8];
				my $sqrt_JSD = $temp[5]/$temp[8];
				my $cov_1 = $temp[6]/$temp[8];
				my $cov_2 = $temp[7]/$temp[8];
				
				#annotation:
				my $gene_name_ens_id = "";
				start_loop: foreach my $gene_start ( sort {$a<=>$b} keys %{$en_genes{$chr}} ) {
					foreach my $gene_end (keys %{$en_genes{$chr}{$gene_start}}) {
						if ( $start >= $gene_start && $start <= $gene_end) {
							$gene_name_ens_id = $en_genes{$chr}{$gene_start}{$gene_end};
							last start_loop;
						} elsif ( $end >= $gene_start && $end <= $gene_end) {
							$gene_name_ens_id = $en_genes{$chr}{$gene_start}{$gene_end};
							last start_loop;
						}
					}
				}
				
				my @temp_2 = split("\t", $gene_name_ens_id );
				if( scalar @temp_2 < 2) {
					$temp_2[0] = "N/A";
					$temp_2[1] = "N/A";
				}
					
				print OUTPUT  "$chr\t$start\t$end\t$cate\t$stat_diff\t$sqrt_JSD\t$cov_1\t$cov_2\t".$temp[8]."\t".$temp_2[0]."\t".$temp_2[1]."\n";
			}
		}
	}
}

#open INPUT, $files[0] or die "cannot open file $files[0].\n";
#while (my $line = <INPUT>) {
#	if ($line =~ m/^chr.+yes/) {
#		chomp $line;
#		my @words = split("\t", $line);
#		my $chr = $words[0];
#		my $start = $words[1];
#		my $end = $words[2];
#		
#		if (exists $all_sig_hits_hash{$chr}{$start}{$end} && $all_sig_hits_hash{$chr}{$start}{$end} >= 3) {
#			print OUTPUT $line."\n";
#		}
#	}
#}
#close INPUT;
close OUTPUT;

exit;

sub GetGeneCoords {
	my ($gene_file) = @_;
	my %gene_coords;
	open INPUT2, $genefile or die "Cannot open $genefile\n";
	while (my $Line = <INPUT2>){
		chomp $Line;
		my @linesplit1 = split(/\t/,$Line);
		if($linesplit1[0] eq 'MT'){
			$linesplit1[0]='M'
		}
		my $chr="chr".$linesplit1[0];
		my $st=$linesplit1[1];
		my $end=$linesplit1[2];
		my $gen=$linesplit1[3]."\t".$linesplit1[4]; # gene name \t gene id
		$gene_coords{$chr}{$st}{$end} = $gen;
	}
	close INPUT2;
	return \%gene_coords;
}


#my @chromosomes = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
#									"chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY");

#foreach my $chr (@chromosomes) {
#	my %asm_hash;
#	my %splice_hash;
#	my @output;
#	open ASM, $asm_file or die "cannot open file $asm_file.\n";
#	while (my $line =<ASM>) {
#		if ($line =~ m/^$chr/) {
#			chomp $line;
#			my @words = split(/\t/, $line);
#			my $start = $words[3];
#			my $end = $words[4];
#			$words[8] =~ m/gene\_id\s\"(\w+)\"\;\stranscript\_id\s\"(\w+)\"/;
#			my $gene_id = $1;
#			my $transcript_id = $2;
#			
#		} else {
#			next;
#		}
#	}
#	close ASM;
	
#	open SPLICE, $splice_file or die "cannot open file $splice_file.\n";
#	while (my $line =<SPLICE>) {
#		if ($line =~ m/^$chr/) {
#			chomp $line;
#			my @words = split(/\t/, $line);
#			
#		} else {
#			next;
#		}
#	}
#	close SPLICE;
#	open OUTPUT, ">$output_table" or die "cannot open file $output_table.\n";
#	close OUTPUT;
#}

