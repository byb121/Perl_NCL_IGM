#!/usr/bin/perl
use strict;
use warnings;

#Version 0.99
#Need to improve the way to output results
#Need to handle the output if the output file is already exist
#because the symbol ">>" used in command line.

my ($folder) = @ARGV;
my @files;

my $output_folder  = "~/arcOGEN/filter_output/";


my @SNP_list = ("rs12901499","rs1052488","rs11556090", "rs11071943", "rs12164949");

opendir (DIR, $folder) or die "Cann't find the directory";
$folder =~ s/\/$//;

while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./);
	my $output_file = $output_folder.$file."_John_2nd_list.txt";
	print $output_file."\n";
	my $input_file = $folder."/".$file;
	print "input file is: $input_file\n";
	print "start to filter.\n";
	foreach my $rs (@SNP_list) {
		`less $input_file | grep -P "$rs\\s+" >> $output_file`;
	}
	print "Done!\n";
}
closedir(DIR);

exit;



