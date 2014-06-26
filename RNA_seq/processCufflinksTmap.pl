#!/usr/bin/perl
use strict;
use warnings;

my @files = (
"f1_s2_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f1_s5_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f1_s7_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f2_s2_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f2_s5_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f2_s6_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f1_s1_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f1_s3_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f1_s4_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f1_s6_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f1_s8_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f2_s1_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f2_s3_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f2_s4_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f2_s7_ucsc_anno/ALL_compare.transcripts.gtf.tmap",
"f2_s8_ucsc_anno/ALL_compare.transcripts.gtf.tmap");

my %trancripts_expression_hash;

for(my $i=0;$i<=15;$i++){
	my $input_file = $files[$i];

	open INPUT, $input_file or die "cannot open file $input_file.\n";
	while (my $line =<INPUT>) {
		chomp $line;
		if ($line =~ m/\tj\t/){
			my $len;
			my @elements_1 = split(/\t/, $line);
			my $id = $elements_1[1];
			if (exists $trancripts_expression_hash{$id}) {
				my @temp = split("_", $trancripts_expression_hash{$id});
				if ($temp[1] <= $i) {
					$temp[0] = $temp[0] + 1;
					$temp[1] = $i + 1;
					$trancripts_expression_hash{$id} = $temp[0]."_".$temp[1]; 
				} 
			} else {
				my $temp = $i+1;
				$trancripts_expression_hash{$id} = "1"."_".$temp;
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
my $output_file  ="temp_tmap.txt";

open OUTPUT, ">$output_file" or die "cannot open file $output_file.\n";
foreach my $key (keys %trancripts_expression_hash) {
	my @temp = split("_",$trancripts_expression_hash{$key});
	if ($temp[0] == 16) {
		print OUTPUT $key."\t".$id_map{$key}."\t".$symbol_map{$id_map{$key}}."\t".$trancripts_expression_hash{$key}."\n";
	}
}
close OUTPUT;

exit;


