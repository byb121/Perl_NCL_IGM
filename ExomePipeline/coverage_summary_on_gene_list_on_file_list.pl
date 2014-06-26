#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

###############################################
# Emit vary tables of coverage summary on exons of a gene list
# include average coverage on exons of each sample
# more than a certain times covered percentages of each sample
###############################################

my $Output_prefix="prefix_test";
my $SamplePath="/users/a5907529/lustre/Kidney_20130211/A1569_Fastq/";
my $SampleNames="Sample_1,Sample_2,Sample_3,Sample_4,Sample_5,Sample_6,Sample_7,Sample_8,Sample_D09254,Sample_D09255,Sample_D115953,Sample_D116785,Sample_D118160,Sample_D156602,Sample_D21608,Sample_D21683,Sample_D44659,Sample_D53464,Sample_D53686,Sample_D53687,Sample_D62086,Sample_D62088,Sample_D74357,Sample_D79574";
my $Gene_list_file="/users/a5907529/lustre/Kidney_20130211/A1569_Fastq/john_list.txt";

my $Results=GetOptions("Output_prefix=s"=>\$Output_prefix, "SampleNames=s"=>\$SampleNames, "Gene_list_file=s"=>\$Gene_list_file);

my @Names = split(",", $SampleNames);
$SamplePath =~ s/\/$//;
my $output_file_mean_coverage = $SamplePath."/".$Output_prefix."_mean_coverage_on_exons_of_genes.txt";
my $output_file_1x_coverage = $SamplePath."/".$Output_prefix."_1x_coverage_on_exons_of_genes.txt";
my $output_file_20x_coverage = $SamplePath."/".$Output_prefix."_20x_coverage_on_exons_of_genes.txt";
my $output_file_40x_coverage = $SamplePath."/".$Output_prefix."_40x_coverage_on_exons_of_genes.txt";


my %Gene_List;
open GeneList, $Gene_list_file or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$Gene_list_file."\n";
while (my $line =<GeneList>) {
	chomp $line;
	$line =~ s/^\s//;
	my @words = split(/ +/,$line);
	$Gene_List{$words[0]} = 0;
}
close GeneList;
print "Interests Genes are in file: $Gene_list_file \n";
print "Captured genes are: \n";
foreach my $i_gene (keys %Gene_List) {
	print $i_gene."\n";
}

print AverageCoverageOnExonsOfSamples($SamplePath, \@Names, \%Gene_List, $output_file_mean_coverage);
print FractionOfLeastCoverageOnExonsOfSamples($SamplePath,  \@Names, \%Gene_List, 1, $output_file_1x_coverage);
print FractionOfLeastCoverageOnExonsOfSamples($SamplePath,  \@Names, \%Gene_List, 20, $output_file_20x_coverage);
print FractionOfLeastCoverageOnExonsOfSamples($SamplePath,  \@Names, \%Gene_List, 40, $output_file_40x_coverage);

exit;

sub AverageCoverageOnExonsOfSamples{
	
	my ($SamplePath, $names_ref, $gene_list_ref, $output_file) = @_;
	
	my @Names = @{$names_ref};
	my %Gene_List = %{$gene_list_ref};
	
	my %mean_count_on_exons;
	
	foreach my $sample (@Names) {
		my $file = $SamplePath."/".$sample."/"."bedtools"."/"."coverage_on_targets.txt";
		open INPUT, $file or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$file."\n";
		while (my $line =<INPUT>) {
			if ($line =~ m/^chr/) {
				chomp $line;
				my @words = split(/\t/,$line);
				if ($words[3] =~ m/^chr.+\:.+\:(\w+)$/) {
					if(exists $Gene_List{$1}) {
						if ( exists $mean_count_on_exons{$words[3]}{$sample} ) {
							$mean_count_on_exons{$words[3]}{$sample} = $mean_count_on_exons{$words[3]}{$sample} + $words[6]*$words[9];
						} else {
							$mean_count_on_exons{$words[3]}{$sample} = $words[6]*$words[9];
						}
					}	
				}
			}
		}
		close INPUT;
		print "coverage_summary_on_gene_list_on_file_list.pl: $file\n";
	}
	
	open OUTPUT, ">$output_file" or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$output_file." to output.\n";
	foreach my $exon (sort keys %mean_count_on_exons) {
		print OUTPUT $exon;
		foreach my $sample (@Names) {
				print OUTPUT "\t".sprintf( "%.1f",$mean_count_on_exons{$exon}{$sample});
		}
		print OUTPUT "\n";
	}
	close OUTPUT;
	
	return "Success.\n";
}

sub FractionOfLeastCoverageOnExonsOfSamples {
	
	my ($SamplePath, $names_ref, $gene_list_ref, $least_coverage, $output_file ) = @_;
	
	my @Names = @{$names_ref};
	my %Gene_List = %{$gene_list_ref};
	
	my %fraction_on_exons;
	
	foreach my $sample (@Names) {
		my $file = $SamplePath."/".$sample."/"."bedtools"."/"."coverage_on_targets.txt";
		open INPUT, $file or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$file."\n";
		while (my $line =<INPUT>) {
			if ($line =~ m/^chr/) {
				chomp $line;
				my @words = split(/\t/,$line);
				if ($words[3] =~ m/^chr.+\:.+\:(\w+)$/) {
					if(exists $Gene_List{$1}) {
						if ( $words[6] >= $least_coverage ){
								if ( exists $fraction_on_exons{$words[3]}{$sample} ) {
								$fraction_on_exons{$words[3]}{$sample} = $fraction_on_exons{$words[3]}{$sample} + $words[9];
							} else {
								$fraction_on_exons{$words[3]}{$sample} = $words[9];
							}
						} else {
							if ( !exists $fraction_on_exons{$words[3]}{$sample} ) {
								$fraction_on_exons{$words[3]}{$sample} = 0
							}
						}
					}	
				}
			}
		}
		close INPUT;
		print "coverage_summary_on_gene_list_on_file_list.pl: $file\n";
	}
	
	open OUTPUT, ">$output_file" or die "coverage_summary_on_gene_list_on_file_list.pl: Cannot open ".$output_file." to output.\n";
	foreach my $exon (sort keys %fraction_on_exons) {
		print OUTPUT $exon;
		foreach my $sample (@Names) {
				print OUTPUT "\t".sprintf( "%.1f", $fraction_on_exons{$exon}{$sample}*100 );
		}
		print OUTPUT "\n";
	}
	close OUTPUT;
	return "Success.\n";
}

