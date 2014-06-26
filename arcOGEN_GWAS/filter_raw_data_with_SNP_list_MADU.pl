#!/usr/bin/perl
use strict;
use warnings;

my ($folder) = @ARGV;
my @files;

my @SNP_list = ("rs9806590","rs7166081");

opendir (DIR, $folder) or die "Cann't find the directory";
while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./);
	print $file."\n";
	my $output_file = $file."_Emma_list.txt";
	print "1\n";
	foreach my $rs (@SNP_list) {
		`less $file | grep -P "$rs\\s+" >> $output_file`;
	}
	print "2\n";
}
closedir(DIR);

exit;



