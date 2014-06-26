#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my ($folder) = @ARGV;
my @files;

my @SNP_list = ("rs2615977","rs1676486","rs9659030");
my $dir = cwd();

opendir (DIR, $folder) or die "Cann't find the directory";
while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./);
	print $file."\n";
	my $input_file = $folder.$file;
	my $output_file = $dir."/".$file."_John_list.txt";
	print "1\n";
	print $output_file."\n";
	foreach my $rs (@SNP_list) {
		`less $input_file | grep -P "$rs\\s+" >> $output_file`;
	}
	print "2\n";
}

closedir(DIR);

exit;



