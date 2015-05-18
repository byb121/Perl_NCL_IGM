#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


my $file1;
my $file2;
my $out;

my $help;
usage() if ( @ARGV < 3 || ! GetOptions('help|?' => \$help, "file1=s"=>\$file1, "file2=s"=>\$file2, 
	'out=s' => \$out) || defined $help );

my %f1;
my $file1_col_num=0;
my %all;

print "Read in file1 $file1\n";
open F1, "$file1" or die "Can not open the file: $file1\n";
while (my $line = <F1> ) {
	chomp $line;
	$line =~ s/\r+$//;
	my @eles = split ("\t", $line);
	if (scalar @eles > $file1_col_num) {$file1_col_num = scalar @eles;}
	if(!exists $f1{$eles[2]}{$eles[3]}{$eles[6]}) {
		$f1{$eles[2]}{$eles[3]}{$eles[6]} = $line;
		$all{$eles[2]}{$eles[3]}{$eles[6]} = 1;
	} else {
		print "duplicated entres on $eles[2] $eles[3] in file $file1.\nBad input.\nExit!\n";
		exit;
	}
}
close F1;
print "Done!\n";

print "Read in file2 $file2\n";
my %f2;
open F2, "$file2" or die "Can not open the file: $file2\n";
while (my $line = <F2> ) {
	chomp $line;
	$line =~ s/\r+$//;
	my @eles = split ("\t", $line);	
	if(!exists $f2{$eles[2]}{$eles[3]}{$eles[6]}) {
		$f2{$eles[2]}{$eles[3]}{$eles[6]} = $line;
		$all{$eles[2]}{$eles[3]}{$eles[6]} = 1;
	} else {
		print "duplicated entres on $eles[2] $eles[3] in file $file2.\nBad input.\nExit!\n";
		exit;
	}
	
}
close F2;
print "Done!\n";

print "Generate output: $out\n";
open OUT, ">$out" or die "Can not open the file: $out\n";
foreach my $chr ( sort {$a cmp $b} keys %all ) {
	foreach my $pos ( sort {$a cmp $b} keys %{$all{$chr}} ) {
		foreach my $var (keys %{$all{$chr}{$pos}}) {
			my $out_string="";
			if (exists $f1{$chr}{$pos}{$var} && exists $f2{$chr}{$pos}{$var}) {
				#print "line1: $f1{$chr}{$pos}{$var} line2: $f2{$chr}{$pos}{$var}\n";
				$out_string=$f1{$chr}{$pos}{$var}."\t";
				my @eles = split ("\t", $f1{$chr}{$pos}{$var});
				if (scalar @eles < $file1_col_num) {
					for(my $i=0;$i < $file1_col_num - scalar @eles;$i++) {
						$out_string=$out_string."\t";
					}
				}
				$out_string=$out_string.$f2{$chr}{$pos}{$var}."\n";
			} elsif (exists $f1{$chr}{$pos}{$var}) {
				$out_string=$f1{$chr}{$pos}{$var}."\n";
			} elsif (exists $f2{$chr}{$pos}{$var}) {
				for(my $i=0;$i < $file1_col_num;$i++) {
					$out_string=$out_string."\t";
				}
				$out_string=$out_string.$f2{$chr}{$pos}{$var}."\n";
			} else {
				print "Bug!\nExit!\n";
				exit;
			}
			print OUT $out_string;
		}
	}
}
close OUT;


print "Done!";

exit;

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "--file1 the first file to merge;\n";
    print "--file2 the second file to merge;\n"; 
    print "--out output.\n";
    return(1);
}