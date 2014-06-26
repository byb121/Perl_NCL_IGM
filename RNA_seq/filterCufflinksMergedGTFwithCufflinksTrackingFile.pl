#!/usr/bin/perl
use strict;
use warnings;

print "The script will remove those lines in gtf for isoforms not passed certain conditions in tracking file \n";
print "columns in the tracking file used for filtering including (either counts greater than the threshold will be passed):\n";
print "\"q0 count\" >= 1\n";
print "\"q0 status\" == OK\n";
print "\"q1 count\" >= 1\n";
print "\"q1 status\" == OK\n";
print "################################################################################################\n";
print "\n";

my ($gtf_file, $tracking_file )= @ARGV;
my $count_threshold = 1; # can change for different count filtering

my %isoform_ID;
print "Reading in the tracking file $tracking_file...    ";
open TRACKING, "$tracking_file" or die "Cannot open the file $tracking_file to read. \n";
while (my $line = <TRACKING>) {
	chomp $line;
	if ($line =~ m/^tracking/) {
		next;
	} else {
		my @ele = split("\t", $line);
		if ($ele[5] eq "OK" && $ele[10] eq "OK") {
			if ($ele[1] >= $count_threshold || $ele[6] >= $count_threshold) {
				$isoform_ID{$ele[0]} = 1;
			}
		}
	}
}
close TRACKING;
print "DONE!!\n\n";

my @output;
print "Reading in the gtf file $gtf_file...    ";
open GTF, "$gtf_file" or die "Cannot open the file $gtf_file to read.\n";
while (my $line = <GTF>) {
	chomp $line;
	my @ele = split ("\t", $line);
	my $transcript_id;
	if ($ele[8] =~ m/transcript\_id.{1}\"(.+?)\"/) {
		$transcript_id = $1;
	} else {
		next;
	}
	if (exists $isoform_ID{$transcript_id}) {
		push @output, $line."\n";
	}
}
close GTF;
print "DONE!!\n\n";

my $output_file = $gtf_file.".filtered";
print "Printing results to file $output_file...    ";
open OUTPUT, ">$output_file" or die "Cannot open the file to ouput $output_file\n";
print OUTPUT @output;
close OUTPUT;
print "DONE!!\n\n";

exit;