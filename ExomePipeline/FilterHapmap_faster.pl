#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($hapfile, $targets_file);
$targets_file = "/home/yaobo/Documents/GenomeData/HumanGenome/hg19/TruSeq-Exome-Targeted-Regions.txt";

my $Results=GetOptions("hapmap=s"=>\$hapfile, "Targets=s"=>\$targets_file);

my $output_file = $hapfile."_filtered_500bp.txt";

my $hg19_file="/home/yaobo/Documents/GenomeData/HumanGenome/hg19/hg19.fa";
my %hg19;

print "Start to read in reference file.\n";
open HG19, "$hg19_file" or die "Can not open the file $hg19_file\n";
my $string = "";
my $current_chr = "";
while (my $line = <HG19> ) {
	chomp $line;
	if ($line =~ m/^\#/) {
		next;
	} elsif( $line =~ m/^\>/) {
		$line =~ s/^\>//;
		if ($current_chr eq "") {
			$current_chr = $line;
			print "Reading in $current_chr...\n";
		} else {
			$hg19{$current_chr} = $string;
			$string = "";
			$current_chr = $line;
			print "Reading in $current_chr...\n";
		}
	} else {
		$string = $string.$line;
	}
}
close HG19;
print "Start to read inhapmap file.\n";


my @filtered_hapmap;
#Filter 1: if the hapmap ref is the same as hg19 ref
print $hapfile."\n";
open HAPMAP, "$hapfile" or die "Can not open the file $hapfile\n";
while (my $line = <HAPMAP>) {
	chomp $line;
	if ($line =~ m/^chromosome/) {
		next;
	} else {
		my @words = split("\t", $line);
		my $chr = $words[0];
		my $pos = $words[1];
		my $change = $words[3];
		my $ref = $words[4];
		my $ref_freq = $words[5];
		if ( uc($ref) eq uc( substr($hg19{$chr}, $pos-1, 1) ) ) {
			push @filtered_hapmap, $line;
		}
	}
}
close HAPMAP;
 
 undef %hg19;  #try to release the memory

#Filter 2: If they are on targeted exomes.
my %targets;
my %targets_start;
my %targets_end;
my @filtered_hapmap_2;
if (-e $targets_file) {
	print "target file is $targets_file\n";
	open TARGETS, "$targets_file" or die "Cannot open the target file $targets_file\n";
	while ( my $line =<TARGETS>) {
		chomp $line;
		if ($line =~ m/^chr/) {
			my @words = split ("\t", $line);
			my $chr = $words[0];
			my $start = $words[1];
			my $end = $words[2];
			$targets{$chr}{$start} = $end;
		}
	}
} else {
	print "No targets file provided!\n";
}
close TARGETS;

#sort targets
foreach my $chr (keys %targets) {
	foreach my $start ( sort {$a<=>$b} keys %{$targets{$chr}} ) {
		if (exists $targets_start{$chr}) {
			push @{$targets_start{$chr}}, $start;
			push @{$targets_end{$chr}}, $targets{$chr}{$start};
		} else {
			my @array_1=($start);
			my @array_2=($targets{$chr}{$start});
			$targets_start{$chr} = [@array_1];
			$targets_end{$chr} = [@array_2];	
		}
	}
}
undef %targets;

# output directly after filtering
open OUT, ">$output_file" or die "Cannot open the file $output_file to output.\n";
foreach my $line (@filtered_hapmap) {
	my @words = split("\t", $line);
	my $chr = $words[0];
	my $pos = $words[1];
	my $change = $words[3];
	my $ref_freq = $words[5];
	my $rs_id = $words[8];
	start_loop: for (my $i=0;$i<scalar @{$targets_start{$chr}};$i++){
		#print "array $chr length ";
		#print scalar @{$targets_start{$chr}};
		#print "\t".$targets_start{$chr}[$i]."\t".$targets_end{$chr}[$i];		
		#print "\n";
		
		if ($pos > $targets_end{$chr}[$i] + 500) {
			next;
		} elsif (  $pos <= ( $targets_end{$chr}[$i] + 500 ) && $pos >= ($targets_start{$chr}[$i] - 500) ) {
			print "found one.\n";
			print OUT  "$rs_id\t$change\t$chr\t$pos\t$ref_freq\n";
			last start_loop;
		} elsif ($pos < ($targets_start{$chr}[$i] - 500) ) {
			last start_loop;
		}
	}
}
close OUT;

exit;

