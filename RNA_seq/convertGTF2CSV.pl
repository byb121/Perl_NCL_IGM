#!/usr/bin/perl
use strict;
use warnings;

my ($gtf_file) = @ARGV;
my $csv = $gtf_file.".csv";

my $csv_header = "chr,source,feature,start,end,score,strand,frame,gene_id,transcript_id,gene_type,gene_status,gene_name,transcript_type,transcript_status,transcript_name";

open GTF, $gtf_file or die "cannot open file $gtf_file.\n";
open OUTPUT, ">$csv" or die "cannot open the file $csv to write.\n";
print OUTPUT $csv_header."\n";
print "start to process the gtf file $gtf_file.\n";
print "......\n";
while (my $line =<GTF>) {
	chomp $line;
	if($line =~ m/^chr/) {
		my @columns = split(/\t/, $line);
		if($columns[8] =~ m/gene\_id.+transcript\_id.+gene\_type.+gene\_status.+gene\_name.+transcript\_type.+transcript\_status.+transcript\_name/) {
			print OUTPUT $columns[0].",".$columns[1].",".$columns[2].",".$columns[3].",".$columns[4].",".$columns[5].",".$columns[6].",".$columns[7]; 
			my @tags =  split(/;/, $columns[8]);
			foreach my $element (@tags) {
				if ($element =~ m/^gene\_id\s\"(.+)\"/) {
					print OUTPUT ",".$1;
				}elsif ($element =~ m/^\stranscript\_id\s\"(.+)\"/) {
					print OUTPUT ",".$1;
				}elsif ($element =~ m/^\sgene\_type\s\"(.+)\"/) {
					print OUTPUT ",".$1;
				}elsif ($element =~ m/^\sgene\_status\s\"(.+)\"/) {
					print OUTPUT ",".$1;
				}elsif ($element =~ m/^\sgene\_name\s\"(.+)\"/) {
					print OUTPUT ",".$1;
				}elsif ($element =~ m/^\stranscript\_type\s\"(.+)\"/) {
					print OUTPUT ",".$1;
				}elsif ($element =~ m/^\stranscript\_status\s\"(.+)\"/) {
					print OUTPUT ",".$1;
				}elsif ($element =~ m/^\stranscript\_name\s\"(.+)\"/) {
					print OUTPUT ",".$1;
				}else{
					next;
				}																																					
			}
		} else {
			print "$line is not cosistant with pattern. ignored.\n";
			next;
		}
		print OUTPUT "\n";
	}
}
close(OUTPUT);
close(GTF);

print "Converting is done.\n";
print "The output is $csv.\n";

exit;
