#!/usr/bin/perl
use strict;
use warnings;

my $input_file_1; #"-f1"
my $input_file_2; #"-f2"
my $rank_list_based_on_column; #"-r int[2345...]]"#if the list is seperated based on fold change, then min and max fold change values will be added in output
my $split_length; #"-l int = 1"
my $split_number; #"-n int = 0"

my ($input_line) = @ARGV;
my @argus = split($input_line, "-");

foreach my $value (@argus) {
	if ($value =~ m/^f1/) {
		if ($input_file_1) {
			print "ERROR: input file one can not be defined twice\n";
			exit;
		} else {
			my @temp = split ($value, /\s/);
			if ($temp[1]) {
				$input_file_1 = $temp[1];
			} else {
				print "ERROR: input file one is absent\n";
				exit;
			}
		}
	}

	if ($value =~ m/^f2/) {
		if ($input_file_2) {
			print "ERROR: input file two can not be defined twice\n";
			exit;
		} else {
			my @temp = split ($value, /\s/);
			if ($temp[1]) {
				$input_file_2 = $temp[1];
			} else {
				print "ERROR: input file two is absent\n";
				exit;
			}
		}
	}
	
	if ($value =~ m/^r/) {
		if ($rank_list_based_on_column) {
			print "ERROR: column of ranking based on can not be defined twice\n";
			exit;
		} else {
			my @temp = split ($value, /\s/);
			if ($temp[1] && isint($temp[1])) {
				$rank_list_based_on_column = $temp[1];
			} else {
				print "ERROR: column of ranking based on is absent or not integer.\n";
				exit;
			}
		}
	}
	
	if ($value =~ m/^l/) {
		if ($split_length) {
			print "ERROR: split length can not be defined twice\n";
			exit;
		} else {
			my @temp = split ($value, /\s/);
			if ($temp[1] && isint($temp[1])) {
				$split_length = $temp[1];
			} else {
				print "ERROR: split length is absent or not integer\n";
				exit;
			}
		}
	}
	
	if ($value =~ m/^n/) {
		if ($split_number) {
			print "ERROR: split number can not be defined twice\n";
			exit;
		} else {
			my @temp = split ($value, /\s/);
			if ($temp[1] && isint($temp[1]) ) {
				$split_number = $temp[1];
			} else {
				print "ERROR: split number is absent or not integer\n";
				exit;
			}
		}
	}
}			 
	
if ($split_number and $split_length) {
	print "ERROR: split length and number can not be defined at same time.\n";
	exit;
}

if (!$input_file_1) {
	print "input file one can not be empty";
	exit;
}

if (!$input_file_2) {
	print "input file two can not be empty";
	exit;
}

if (!$split_length && !$split_number) {
	$split_length = 1;
} else {
	
}



print "-f1"."\t".$input_file_1."-f2"."\t".$input_file_2."-r"."\t".$rank_list_based_on_column."-l"."\t".$split_length."-r"."\t".$split_number."\n";
exit;
