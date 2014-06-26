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

my $depth_filter = 0; #locus with seqencing depth lower than the value will be dropped.
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

print "################################################################\n";
print "#                     Asign SNP to genes                       #\n";
print "#REF/Var could be biased, check the script it self for detail./#\n";
print "################################################################\n";
print "\n";
print "######\n";
print "Warning: Files will only have the ref and var bases extracted from the last vcf file processed. God blessed all samples have the same ref and var!\n";
print "######\n";
print "\n";

my %annoted_hash;
my $temp_output_file = $dir."/"."temptemptemp_file.txt";
my $temp_output_file_2 = $dir."/"."temptemptemp_file_with_genes_Names.txt";

open OUTPUT, ">$temp_output_file" or die "cannot open file $temp_output_file.\n";
for my $key ( keys %hash ) {
	if ($hash{$key} == 16) {
		my @words = split("_", $key);
		print OUTPUT $words[0]."\t".$words[1]."\t".$words[1]."\n";
	}
}
close OUTPUT;

`/users/a5907529/biosofts/BEDTools-Version-2.12.0/bin/intersectBed -wo -a $temp_output_file -b ~/GenomeData/hg19/gencode.v10.annotation.gtf > $temp_output_file_2`;

open INPUT, $temp_output_file_2 or die "cannot open file $temp_output_file_2.\n";
while (my $line =<INPUT>) {
	chomp $line;
	if( $line =~ m/^chr/ ){
		my @words = split("\t", $line);
		my $gene_name;
		my $gene_id;
		my $key = $words[0]."_".$words[1];
		if (exists $words[11]) {
			my @attributes = split ("; ", $words[11]);
			#print "word 11 are::";
			#print @attributes;
			#print "\n";
			
			foreach my $attribute (@attributes) {
				if ($attribute =~ m/^gene\_id/){
					$attribute =~ s/gene\_id//;
					$attribute =~ s/\"//g;
					#print "subtitued gene_id is:";
					#print $attribute;
					#print "\n";
					$gene_id = $attribute;
				}
				
				if ($attribute =~ m/^gene\_name/){
					$attribute =~ s/gene\_name//;
					$attribute =~ s/\"//g;
					#print "subtitued gene_name is:";
					#print $attribute;
					#print "\n";
					$gene_name = $attribute;
				}
			}
			$annoted_hash{$key} = $gene_name."_".$gene_id;	
		}
	}
}
close INPUT;

`rm $temp_output_file`;
`rm $temp_output_file_2`;

##Counting Genotype#

my %genoType_hash;
my %ref_var_hash; # this hash will only contain the ref and var from the last vcf file processed. God blessed all samples have the same ref and var!"
my %homo_count;
my %hetre_count;
my %export_2_R;

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
				$ref_var_hash{$key} = $words[3]."/".$words[4];
				if (exists $annoted_hash{$key}) {
					my @vcf_info = split (";", $words[7]);
					my @vcf_sample = split (":", $words[9]);
					my $ref2var_coverage;
					foreach my $tag(@vcf_info) { # it requires that the vcf file info column has to have a DP4 tag
						if ($tag =~ m/^DP4/) {
							$tag =~ s/^DP4\=//;
							my @numbers = split (",", $tag);
							my $temp1 = $numbers[0] + $numbers[1]; #extract coverage of ref and revers_ref
							my $temp2 = $numbers[2] + $numbers[3]; #extract coverage of var and revers_var
							$ref2var_coverage = $temp1.":".$temp2; #recorde the ratio
							if (exists $genoType_hash{$key}) {
								$genoType_hash{$key} = $genoType_hash{$key}."\t".$vcf_sample[0]."_".$ref2var_coverage; #genotye (eg 0/1 or 1/1) and base ratio
							} else {
								$genoType_hash{$key} = $vcf_sample[0]."_".$ref2var_coverage;
							}
							
							if (exists $export_2_R{$key}) {
								$export_2_R{$key} = $export_2_R{$key}."\t".$words[3]."\t".$temp1."\t".$words[4]."\t".$temp2; 
								#ref_base	ref_freq	var_base	var_freq
							} else {
								$export_2_R{$key} = $words[3]."\t".$temp1."\t".$words[4]."\t".$temp2;
							}
							
						} else {
							
							next;
						}
					}
					
					if (!exists $hetre_count{$key}) {# initial the count to allow value 0
						$hetre_count{$key} = 0;
					}
					if (!exists $homo_count{$key}) {
						$homo_count{$key} = 0;
					}
					
					if ($vcf_sample[0] =~ m/0\/1/) {# start to count
						$hetre_count{$key} += 1;
					} else {
						$homo_count{$key} += 1;						
					}
					
				} else {
					next;
				}
			} 
		}
		close INPUT;
	}
}
closedir(DIR);


##### To calculate the average imbalance
###Choose SNP that are hetreozygous in at least 3 samples in OA and NOF group
###Remove homozygous at the position
###Return avg(Nr/(Nr+Nv))

my %OA_group = ( 2 => 1, 3 => 1, 4 => 1, 5 => 1, 6 => 1, 7 => 1,  8 => 1, 9 => 1, 11 => 1, 12 => 1);
my %NOF_group = ( 1 => 1, 10 => 1, 13 => 1, 14 => 1, 15 => 1, 16 => 1);
my @output;
my @output_hetereo_freq;

for my $key ( keys %genoType_hash ) {
	my @words = split("_", $key);
	if ($words[0] =~ m/chrY/){
		next;
	} else {
		#print OUTPUT $words[0]."\t".$words[1]."\t".$genoType_hash{$key}."\t".$hetre_count{$key}."\t".$homo_count{$key}."\n";
		my $OA_hetero_count = 0;
		my $NOF_hetero_count = 0;
		my @OA_imbalance_value;
		my @NOF_imbalance_value;
		
		my @sample_genoType = split("\t", $genoType_hash{$key});
		for(my $i=0;$i<=15;$i++){
			my @temp = split("_", $sample_genoType[$i]);
			my @numbers = split(":", $temp[1]); 
			if($temp[0] =~ m/0\/1/) {
				if(exists($OA_group{$i+1})){
					$OA_hetero_count += 1;
					push (@OA_imbalance_value, $numbers[0]/($numbers[0]+$numbers[1]));
				} elsif (exists($NOF_group{$i+1})) {
					$NOF_hetero_count += 1;
					push (@NOF_imbalance_value, $numbers[0]/($numbers[0]+$numbers[1]));
				} else {
					print "Fatal error, Grouping failed!"."\n";
					exit;
				}
			}
		}
		my @words_2 = split("_", $annoted_hash{$key});
		
		push(@output_hetereo_freq, $words[0]."\t".$words[1]."\t".$words_2[0]."\t".$words_2[1]."\t".$OA_hetero_count."\t".$NOF_hetero_count."\n" );
		
		if ($OA_hetero_count >= 3 && $NOF_hetero_count >= 3) {
			my $average_OA;
			my $average_NOF;
			
			my $j = 0;
			foreach my $temp_OA (@OA_imbalance_value) {
				$j += $temp_OA;					
			}
			$average_OA = $j/$OA_hetero_count;
			
			$j = 0;
			foreach my $temp_NOF (@NOF_imbalance_value) {
				$j += $temp_NOF;					
			}
			$average_NOF = $j/$NOF_hetero_count;
			
			my $ttest = new Statistics::TTest;
			$ttest->set_significance(95);
			$ttest->load_data(\@OA_imbalance_value,\@NOF_imbalance_value);
			my $t    = $ttest->t_statistic;
			my $df   = $ttest->df;
			my $prob = $ttest->{t_prob};
			
			my @words_2 = split("_", $annoted_hash{$key});
			my $output_string = $words[0]."\t".$words[1]."\t".$ref_var_hash{$key}."\t".$words_2[0]."\t".$words_2[1]."\t".$OA_hetero_count."\t".$average_OA."\t".$NOF_hetero_count."\t".$average_NOF."\t"."t=$t (df = $df)"."\t".$prob."\n";
			# output format chromosome_name position Ref/Var Gene_name Gene_id Number_of_Hereo_in_OA	OA_Average_ImB	Number_of_Hereo_in_NOF	NOF_Average_ImB	T_test	P_value
			
			push (@output, $output_string);
		}
	}
}

print "###########################################\n";
print "########### output results ################\n";
print "###########################################\n";

my $output_file = $dir."/"."CommonSNPs_filtered.txt";
print "$output_file\n";
open OUTPUT, ">$output_file" or die "cannot open file $output_file.\n";
for my $key ( keys %hash ) {
	if ($hash{$key} == 16) {
		my @words = split("_", $key);
		print OUTPUT $words[0]."\t".$words[1]."\t".$ref_var_hash{$key}."\t".$words[1]."\n";
	}
}
close OUTPUT;

my $output_file_anno = $dir."/"."CommonSNPs_filtered_annotated.txt";
print "$output_file_anno\n";
open OUTPUT, ">$output_file_anno" or die "cannot open file $output_file_anno.\n";
for my $key ( keys %annoted_hash ) {
	my @words_1 = split("_", $key);
	my @words_2 = split("_", $annoted_hash{$key});
	print OUTPUT $words_1[0]."\t".$words_1[1]."\t".$words_2[0]."\t".$words_2[1]."\n";
}
close OUTPUT;


my $genoType_output_file = $dir."/"."CommonSNPs_with_genoType_info_filtered.txt";
print "$genoType_output_file\n";
open OUTPUT, ">$genoType_output_file" or die "cannot open file $genoType_output_file.\n";
for my $key ( keys %genoType_hash ) {
		my @words = split("_", $key);
		my @words_2 = split("_", $annoted_hash{$key});
		print OUTPUT $words[0]."\t".$words[1]."\t".$words_2[0]."\t".$words_2[1]."\t".$genoType_hash{$key}."\t".$hetre_count{$key}."\t".$homo_count{$key}."\n";
}
close OUTPUT;
	
my $imbalance_output_file = $dir."/"."CommonSNPs_with_imbalance_info_filtered.txt";
print "$imbalance_output_file\n";
open OUTPUT, ">$imbalance_output_file" or die "cannot open file $imbalance_output_file.\n";
print OUTPUT @output; 
close OUTPUT;


my $hetereo_freq_output_file = $dir."/"."CommonSNPs_with_hetereo_freq_info_filtered.txt";
print "$hetereo_freq_output_file\n";
open OUTPUT, ">$hetereo_freq_output_file" or die "cannot open file $hetereo_freq_output_file.\n";
print OUTPUT @output_hetereo_freq; 
close OUTPUT;

my $export_2_R_output_file = $dir."/"."CommonSNPs_info_export_2_R_filtered.txt";
print "$export_2_R_output_file\n";
open OUTPUT, ">$export_2_R_output_file" or die "cannot open file $export_2_R_output_file.\n";
for my $key ( keys %export_2_R ) {
		my @words = split("_", $key);
		my @words_2 = split("_", $annoted_hash{$key});
		print OUTPUT $words[0]."\t".$words[1]."\t".$words_2[0]."\t".$words_2[1]."\t".$export_2_R{$key}."\n";
}
close OUTPUT;

exit;




