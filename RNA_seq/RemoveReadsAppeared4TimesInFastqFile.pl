#!/usr/bin/perl 
use strict;
use warnings;

my ($fastqFile) = @ARGV;

my $output = $fastqFile."ReadsRepeated2Times";
my %readsNameHash;
my @OUTPUT;

open INPUT, $fastqFile or die "cannot open file $fastqFile";
while (my $line =<INPUT>) {
	if($line !~ m/^@/){
		#print $line;
		my @words = split(/\t/,$line);
		if(!exists $readsNameHash{$words[0]}){
			for (keys %readsNameHash)
			{
				delete $readsNameHash{$_};
			}
			$readsNameHash{$words[0]} = 0;
			push @OUTPUT, $line;
		} else {
			$readsNameHash{$words[0]} = $readsNameHash{$words[0]} + 1;
			if ($readsNameHash{$words[0]} <= 1) {
				push @OUTPUT, $line;
			} else {
				next;
			}
		}
	} else {
		push @OUTPUT, $line;
	}
}

close INPUT;

open OUTPUT, ">$output" or die "Cannot open the file $output";
print OUTPUT @OUTPUT;
close OUTPUT;

exit;