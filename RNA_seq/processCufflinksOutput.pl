#!/usr/bin/perl
use strict;
use warnings;

my @files = (
"f1_s2_ucsc_anno/transcripts.gtf",
"f1_s5_ucsc_anno/transcripts.gtf",
"f1_s7_ucsc_anno/transcripts.gtf",
"f2_s2_ucsc_anno/transcripts.gtf",
"f2_s5_ucsc_anno/transcripts.gtf",
"f2_s6_ucsc_anno/transcripts.gtf",
"f1_s1_ucsc_anno/transcripts.gtf",
"f1_s3_ucsc_anno/transcripts.gtf",
"f1_s4_ucsc_anno/transcripts.gtf",
"f1_s6_ucsc_anno/transcripts.gtf",
"f1_s8_ucsc_anno/transcripts.gtf",
"f2_s1_ucsc_anno/transcripts.gtf",
"f2_s3_ucsc_anno/transcripts.gtf",
"f2_s4_ucsc_anno/transcripts.gtf",
"f2_s7_ucsc_anno/transcripts.gtf",
"f2_s8_ucsc_anno/transcripts.gtf");

my %trancripts_expression_hash;



for(my $i=0;$i<=15;$i++){
	my $input_file = $files[$i];
	open INPUT, $input_file or die "cannot open file $input_file.\n";
	while (my $line =<INPUT>) {
		chomp $line;
		if ($line =~ m/\ttranscript\t.+transcript\_id\s\"uc/){
			
			my $fpkm;
			my $cov;
			my $status;
			my $id;
			
			my @elements_1 = split(/\t/, $line);
			
			if (exists $elements_1[8]){
				my @elements_2 = split ("; ", $elements_1[8]);
				foreach my $attributes (@elements_2){
					if ($attributes =~ m/^FPKM/){
						my @temp = split(/\s/, $attributes);
						$temp[1] =~ s/\"//g;
						$fpkm = $temp[1];
					} elsif ($attributes =~ m/^cov/){
						my @temp = split(/\s/, $attributes);
						$temp[1] =~ s/\"//g;
						$cov = $temp[1];
					} elsif ($attributes =~ m/^full\_read/){
						my @temp = split(/\s/, $attributes);
						$temp[1] =~ s/\"//g;
						$status = $temp[1];
					} elsif ($attributes =~ m/^transcript\_id/){
						my @temp = split(/\s/, $attributes);
						$temp[1] =~ s/\"//g;
						$id = $temp[1];
					}
				}
				
				if ($fpkm && $status && $cov && $id) {
					if (exists $trancripts_expression_hash{$id}) {
						my @temp_temp = split(/\t/, $trancripts_expression_hash{$id});
						my $length = scalar @temp_temp;
						if ($length < 3*$i-3) {
							for (my $j=1;$j<=3*$i-$length-3;$j++){
								$trancripts_expression_hash{$id} = $trancripts_expression_hash{$id}."\t";
							}
						}
						$trancripts_expression_hash{$id} = $trancripts_expression_hash{$id}."\t".$fpkm."\t".$cov."\t".$status;
					} else {
						$trancripts_expression_hash{$id} = "";
						for(my $j=0;$j<$i;$j++){
							$trancripts_expression_hash{$id} = $trancripts_expression_hash{$id}."\t"."\t"."\t";
						}
						$trancripts_expression_hash{$id} = $trancripts_expression_hash{$id}."\t".$fpkm."\t".$cov."\t".$status;
					}
				} else {
					print "Error: one or more of id/FPKM/cov/status is lost.\n";
					exit;
				}
			} else {
				print "Error: no 8th elements on the line.\n";
				exit;
			}
		}
	}
	print $input_file." is done!\n";
	close (INPUT);
}

############ map the id to refseq #################

#reading ucsc -> refseq id
my %id_map;
my $id_map_file = "/users/a5907529/GenomeData/hg19/ucsc2refseq.txt";
open UCSC_REF, $id_map_file or die "cannot open file $id_map_file.\n";
while (my $line =<UCSC_REF>) {
	chomp $line;
	if ($line =~ m/^uc/) {
		my @words = split(/\t/, $line);
		$id_map{$words[0]} = $words[1];
	}
}
close(UCSC_REF);

#reading refseq -> symbol
my %symbol_map;
my $ref_symbol_map = "/users/a5907529/RNA_seq/cufflinks_gsnap_v1/all.counts.txt";
open REF_SYM, $ref_symbol_map or die "cannot open file $ref_symbol_map.\n";
while (my $line =<REF_SYM>) {
	chomp $line;
	if ($line =~ m/^gi/) {
		my @words = split(/\t/, $line);
		my @temp = split(/\|/, $words[0]);
		my @temp_2 = split(/\./, $temp[3]);
		$symbol_map{$temp_2[0]} = $words[1];
	}
}
close(REF_SYM);

#################output_part###########################
my $output_file  ="temp.txt";

open OUTPUT, ">$output_file" or die "cannot open file $output_file.\n";
foreach my $key (keys %trancripts_expression_hash) {
	print OUTPUT $key."\t".$id_map{$key}."\t".$symbol_map{$id_map{$key}}."\t".$trancripts_expression_hash{$key}."\n";
}
close OUTPUT;

exit;


