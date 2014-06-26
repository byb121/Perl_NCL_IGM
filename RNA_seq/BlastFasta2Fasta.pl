#!/usr/bin/perl
use strict;
use warnings;

my ($fasta_file) = @ARGV;
my $output_file = $fasta_file.".clustered_table.txt";
my $min_identify_percent = 0.95;
my $min_match_length_percent = 0.95;
my $max_match_length_percent = 1.05;
my $max_mismatch_count_percent = 0.01;
my $max_gap_open_count_percent = 0.01;
my $max_e_value =1e-60;


#output each seq to a fasta file
my %total_fasta;
my $id;

`makeblastdb -in $fasta_file -dbtype nucl -title total_temp_db -out temp`;

#make the blast db
open TOTAL, $fasta_file or die "Can not open the file $fasta_file";
while (my $line = <TOTAL> ){
	if ( $line =~ m/^\>/) {
		chomp $line;
		$id = $line;
		$id =~ s/^\>//;
		print "Reading seq for $id...\n";
	} else {
		if ( !exists $total_fasta{$id}) {
			$total_fasta{$id} = $line;
		} else {
			$total_fasta{$id} = $total_fasta{$id}.$line;
		}
	}
}
close TOTAL;

my $temp_file="temp_1.fasta";
my %clusters;
my $cluster_number = 1;

foreach my $id (keys %total_fasta) {
	if ( ! exists $clusters{$id}) {
		my %this_cluster;
		my $OA_count = 0;
		my $NOF_count = 0;
		open TEMP_1, ">$temp_file" or die "Cannot open the file $temp_file";
		print TEMP_1 ">".$id."\n";
		print TEMP_1 $total_fasta{$id};
		close TEMP_1;
		my $blast_result_file="blast.output";
		#blast the fasta file to the total fasta file 
		`blastn -query $temp_file -db temp -evalue 0.1 -outfmt 6 -out $blast_result_file -max_target_seqs 100`;
		open BLAST_RESULT, $blast_result_file or die "Can not open the file $blast_result_file\n";
		#gather blast result
		#determine similar seq
		
		my @temp = split ("_", $id);
		my $query_length = $temp[10];
		my $query_dis_type = $temp[2];
		$this_cluster{$id} = 1;
		HIT_LOOP: while (my $line = <BLAST_RESULT> ){
			chomp $line;
			my @temp_2 =  split ("\t", $line);
			my $match_id = $temp_2[1];
			my @temp_3 = split("_", $temp_2[1]);
			my $match_dis_type = $temp_3[2];
			my $match_length = $temp_3[10];
			my $e_value = $temp_2[10];
			if ($e_value > $max_e_value) {
				last HIT_LOOP; 
			}
			my $aln_length_percent = $temp_2[2];
			my $match_length_percent = $match_length/$query_length;
			my $mismatch_count_percent = $temp_2[4]/$query_length;
			my $gap_open_count_percent = $temp_2[5]/$query_length;
			if ($aln_length_percent >= $min_identify_percent && $match_length_percent >= $min_match_length_percent && $match_length_percent <= $max_match_length_percent
						&& $mismatch_count_percent <= $max_mismatch_count_percent && $gap_open_count_percent <=  $max_gap_open_count_percent) {
							$this_cluster{$match_id} = 1;
						}
		}
		close BLAST_RESULT;
		
		my $this_cluster_number;
		FIND_CLUSTER: foreach my $clustered_id (keys %this_cluster) {
			if (exists $clusters{$clustered_id}) {
				$this_cluster_number = $clusters{$clustered_id};
				last FIND_CLUSTER;
			}
		}
		unless ($this_cluster_number) {
			$this_cluster_number=$cluster_number;
			$cluster_number += 1;
		}
		
		foreach my $clustered_id (keys %this_cluster) {
			$clusters{$clustered_id} = $this_cluster_number;
		}
	} else {
		next;
	}
}

my %OA_match_count;
my %OA_match_sample;
my %NOF_match_count;
my %NOF_match_sample;

foreach my $value (values %clusters) {
	$OA_match_count{$value} = 0;
	$NOF_match_count{$value} = 0;
	foreach my $key (keys %clusters) {
		if ($clusters{$key} == $value) {
			my  @temp = split("_", $key);
			my $dis_type = $temp[2];
			my $sample = $temp[0]."_".$temp[1];
			if ( $dis_type eq "OA") {
				if ( !exists $OA_match_sample{$value}{$sample} ) {
					$OA_match_sample{$value}{$sample} = 1;
					$OA_match_count{$value} = $OA_match_count{$value} +1;
				}
			} else {
				if ( !exists $NOF_match_sample{$value}{$sample} ) {
					$NOF_match_sample{$value}{$sample} = 1;
					$NOF_match_count{$value} = $NOF_match_count{$value} +1;
				}
			}
		}
	}
}

open OUTPUT, ">$output_file" or die "Can not open the file to ouput $output_file\n";
foreach my $key (keys %clusters) {
	my $value = $clusters{$key};
	print OUTPUT $key."\t".$value."\t".$OA_match_count{$value}."\t".$NOF_match_count{$value}."\n";
}

exit;