#!/usr/bin/perl
use strict;
use warnings;

my ($gtf_file) = @ARGV;
my $csv = $gtf_file.".csv";

my $csv_header = "chr,start,end,gene_id,gene_type,gene_name";

open GTF, $gtf_file or die "cannot open file $gtf_file.\n";
open OUTPUT, ">$csv" or die "cannot open the file $csv to write.\n";
print OUTPUT $csv_header."\n";
print "start to process the gtf file $gtf_file.\n";
print "......\n";
while (my $line =<GTF>) {
	chomp $line;
	if($line =~ m/\tgene\t/) {
		my @columns = split(/\t/, $line);
		my @tags =  split(/;/, $columns[8]);
		my ($gene_id, $gene_name, $gene_type);
		foreach my $element (@tags) {
			if ($element =~ m/^gene\_id\s\"(.+)\"/) {
				$gene_id = $1;
			}elsif ($element =~ m/^\sgene\_biotype\s\"(.+)\"/) {
				$gene_type = $1;
			}elsif ($element =~ m/^\sgene\_name\s\"(.+)\"/) {
				$gene_name = $1;
			}else{
				next;
			}																																	
		}
		if (!defined $gene_id) {
			print "Error, no gene is was found in: $line\n";
		}
		if (!defined $gene_name) {
			$gene_name = "NA";
		}
		if (!defined $gene_type) {
			$gene_type = "NA"		
		}
		print OUTPUT "$columns[0],$columns[3],$columns[4],$gene_id,$gene_type,$gene_name\n";
	}
	
}
close(OUTPUT);
close(GTF);

print "Converting is done.\n";
print "The output is $csv.\n";

exit;
