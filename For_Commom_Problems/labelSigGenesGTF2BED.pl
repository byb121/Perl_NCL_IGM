#!/usr/bin/perl
use strict;
use warnings;
use Math::Round;

################ Intro ##################
# The script will look for those variants that pass the filters.
# If they are exist in at least half of total sampels, 
# the script will look for their position in Gencode V10.
# 
# Genotype of left variants will be outputed into results file.
#########################################


print "\n";
print "#################################################################\n";
print "#                           Intro                               #\n";
print "#################################################################\n";
print "\n";
print "Author: Yaobo Xu.\n";
print "Input: paths of GTF_file and expression file (exported from genespring):\n";
print "       Attributs of GTF file should contain \"gene_name\" tag.\n";
print "Output: Bed file (Can be upload to Genome Browser to view:\n";
print "        It contains all exons annotated by the supplied GTF file.\n";
print "        Exons of those genes that showed significant regulations in the expression file will be colored.\n";
print "        Balck: no significant regulation.\n";
print "        Red: Significant down-regulation.\n";
print "        Green: Significant up-regulation.\n";
print "########################  intro end  ###########################\n";
print "\n";


my ($gtf_file, $expression_file) = @ARGV;
my $output_bed = $expression_file."processed.bed";

# cutoffs
my $p_value_cutoff = 0.05;
my $fold_change_cutoff = 2;

#Firt line of the bed file to display in USCS Genome Browser
my $first_line = "track name=\"d14 vs d0 Expression\" description=\"d14_vs_d0 _MSC_3Samples_each\" visibility=3 itemRgb=\"On\" useScore=1\n";

my %expression_hash;

print "Processing expression file.......";
open EXPRESSION, $expression_file or die "cannot open file $expression_file.\n";

while (my $line =<EXPRESSION>) {
	chomp $line;
	if ($line =~ m/^\d+/){
		my @words = split("\t", $line);
		my $p_value = $words[1];
		my $regulation = $words[4];
		my $fold_change;
		
		if ($regulation eq "up") {
			$fold_change = $words[3];
		} elsif ($regulation eq "down") {
			$fold_change = 0 - $words[3];
		}else {
			print "Error: uexpected regulation direction found: at line $line\n";
		}
			
		my $symbol = $words[7];
				
		if ($p_value <= $p_value_cutoff && abs($fold_change) >= $fold_change_cutoff){
			if (exists $expression_hash{$symbol}) {
				my @values = split("_",$expression_hash{$symbol});
				if ($p_value < $values[0]){
					$expression_hash{$symbol} = $p_value."_".$fold_change;
				}
			} else {
				$expression_hash{$symbol} = $p_value."_".$fold_change;
			}
		}		
	}
}

close (EXPRESSION);
print "Done!\n";

print "Processing annotation file and outputing bed file ..............";
open GTF, $gtf_file or die "cannot open file $gtf_file.\n";
open BED, ">$output_bed" or die "cannot open file $output_bed";

print BED $first_line;

while (my $line = <GTF>) {
	chomp $line;
	if ($line =~ m/^chr/) {
		my @columns = split("\t", $line);
		if ($columns[2] eq "exon") {
			my @elements = split(";", $columns[8]);
			my $gene_name;
			foreach my $attri (@elements) {
				$attri =~ s/^\W|\W$//g; #remove non-word characters from the both ends of the string
				my @temp = split(" ", $attri);
				if ($temp[0] eq "gene_name") {
					$gene_name = $temp[1];
				}
			}
			
			my $colour = "0,0,0"; #decide the color of the exon block
			my $score = 0;
			if (exists $expression_hash{$gene_name}) {
				my @temp = split ("_", $expression_hash{$gene_name});
				#print $expression_hash{$gene_name}."\n";
				my $p_value = $temp[0];
				my $fold_change = $temp[1];
				$score = round(abs($fold_change));
				#print $score."\n";
				if ($fold_change > 0) {
					$colour = "0,238,0"; ##Green
				}elsif($fold_change < 0) {
					$colour = "238,0,0"; ##Red
				} else {
					print "\nWarning: fold change can not be zero! STOP!\n";
					exit;
				}
			}
			print BED $columns[0]."\t".$columns[3]."\t".$columns[4]."\t".$gene_name."\t".$score."\t".$columns[6]."\t".$columns[3]."\t".$columns[4]."\t".$colour."\n";
			
				
		}
		
	}
}

close (GTF);
close (BED);
print "Done!\n";

print "Output file: $output_bed\n";

exit;



