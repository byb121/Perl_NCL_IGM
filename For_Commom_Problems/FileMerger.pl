#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my ($folder) = @ARGV;

my $dir = cwd();
my $output_file = $dir."/"."merged_file.txt";

$folder =~ s/\/$//;

my %file_list;

opendir (DIR, $folder) or die "Cann't find the directory";
while (my $file = readdir(DIR)) {
	next if ($file =~ m/^\./);
	#if ($file =~ m/Novoalign_transcript_reweighted_chr_\d+_\d+\.txt$/){
	#if ($file =~ m/Novoalign_transcript_chr_\d+_\d+\.txt$/){
	if ($file =~ m/regression_array_test_slice_\d+_\d+\.txt$/){
		#$file =~ m/Novoalign_transcript_reweighted_chr_(\d+)_(\d+)\.txt$/;
		#$file =~ m/Novoalign_transcript_chr_(\d+)_(\d+)\.txt$/;
		$file =~ m/regression_array_test_slice_(\d+)_(\d+)\.txt$/;
		$file_list{$file} = $1;
	}
}

foreach my $key ( sort {$file_list{$a} <=> $file_list{$b}} keys %file_list) {
	print "merging $key into output $output_file....\n";
	my $input_file = $folder."/".$key;
	`cat $input_file >> $output_file`;
}

print "Merging is Done!\n";
print "The output is $output_file\n";
exit;