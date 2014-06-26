#!/usr/bin/perl
use strict;
use warnings;

print "The script will add 3 more columns to the diff file including: \"gene name\", \"gene type\", \"gene status\" and \"length\". \n";
print "Annotation can only take from gencode gtf, other gtf may cause problems to match\n";
print "################################################################################################\n";
print "\n";

my ($diff_file, $gene_anno )= @ARGV;

print "Reading in the annotation file...\n";
my %anno;
open ANNO, "$gene_anno" or die "Cannot open the file $gene_anno\n";
while (my $line = <ANNO>) {
	if ($line =~ m/\#/) {
		next;
	} else {
		my @ele = split("\t", $line);
		#my @temp = split(";", $ele[8]);
		
		my $transcript_id;
		my $gene_type;
		my $gene_name;
		my $gene_status;
		
		if ($ele[8] =~ m/transcript\_id.{1}\"(.+?)\"/) {
			$transcript_id = $1;
			$transcript_id =~ s/\.\d+$//;
		} else {
			next;
		}
		if ($ele[8] =~ m/gene\_type.{1}\"(.+?)\"/) {
			$gene_type = $1;
		} else {
			next;
		}
		if ($ele[8] =~ m/gene\_name.{1}\"(.+?)\"/) {
			$gene_name = $1;
		} else {
			next;
		}
		if ($ele[8] =~ m/gene\_status.{1}\"(.+?)\"/) {
			$gene_status = $1;
		} else {
			next;
		}
		
		#print "$transcript_id: $gene_name\t$gene_type\t$gene_status\n";
		if (!exists $anno{$transcript_id}) {
			$anno{$transcript_id} = "$gene_name\t$gene_type\t$gene_status";
		}
	}
}
close ANNO;
print "Anno reading is done.\n";

print "Reading in the diff files now.....\n";
my @output;
open DIFF, "$diff_file" or die "Cannot open the file $diff_file.\n";
while (my $line = <DIFF>) {
	chomp $line;
	if ($line =~ m/^test/) {
		push @output, $line."\tgene_name\tgene_type\tgene_status\tlength\n";
	} else {
		my @ele = split("\t", $line);
		my $start;
		my $end;
		if ($ele[3] =~ m/chr\w+\:(\d+)\-(\d+)$/) {
			$start = $1;
			$end = $2;
		} else {
			$start = 0;
			$end = 0;
		}
		my $length = $end-$start+1;
		if ( $ele[2] eq "-" ) {
			push @output, $line."\tNA\tNA\tNA\t$length\n";
		} else {
			my @temp = split (",", $ele[2]);
			my $transcript_id = $temp[0];
			$transcript_id =~ s/\.\d+$//;
			if (exists $anno{$transcript_id} ) {
				push @output, $line."\t$anno{$transcript_id}\t$length\n";
			} else {
				push @output, $line."\tNA\tNA\tNA\t$length\n";
			}
		}
	}
}
close DIFF;
print "Reading is done..... \n";

my $output_file = $diff_file.".moreCols";
print "Output new file to $output_file.\n";
open OUTPUT, ">$output_file" or die "Cannot open file $output_file to output. \n";
print OUTPUT @output;
close OUTPUT;

exit;
