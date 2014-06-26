#!/usr/bin/perl
use strict;
use warnings;

my ($target_gene_file) = @ARGV;
my $id_map_file = "/home/yaobo/Documents/GenomeData/HumanGenome/hg19/ensembl69_id_map.txt";

my %targets;

print "Reading target file $target_gene_file.\n";
open TARGET_FILE, $target_gene_file or die "cannot open file $target_gene_file.\n";
my $count = 0;
while (my $line =<TARGET_FILE>) {
	if($line =~ m/^ENSG/) {
		chomp $line;
		my @cols = split("\t", $line);
		if (exists $targets{$cols[0]}) {
			print "Duplicated entry found $cols[0], removed!\n";
		}	else {
			my $chr = $cols[1];
			my $start = $cols[2];
			my $end = $cols[3];
			$targets{$cols[0]} = $chr."_".$start."_".$end;
			$count += 1;
		}
	}
}
close(TARGET_FILE);

print "Total target found: $count.\n";

my %genes;
my %symbols;

print "reading gene file $id_map_file.\n";
open MAP, $id_map_file or die "Can not open the file $id_map_file.\n";
while(my $line=<MAP>) {
	chomp $line;
	my @cols = split(",", $line);
	my $chr = $cols[1];
	my $start = $cols[2];
	my $end = $cols[3];
	$genes{$cols[0]} = $chr."_".$start."_".$end;
	if ($cols[4] ne "") {
		$symbols{$cols[0]} = $cols[4];
	}
}
close(MAP);

#Start to search
my %upstream;
my %downstream;
my %overlap;

foreach my $target (keys %targets) {
	
	my $up_shortest_dist = 10000000;
	my $down_shortest_dist = 10000000;
	my @target_coords = split("_",$targets{$target});
	my $target_chr = $target_coords[0];
	my $target_start = $target_coords[1];
	my $target_end = $target_coords[2];
	print "Searching neighbouring genes of $target on chromosome: $target_chr".":".$target_start."-".$target_end."\n";
	foreach my $gene (keys %genes) {
		if($gene eq $target) {
			next;
		} else {
			my @coords = split("_",$genes{$gene});
			my $chr = $coords[0];
			my $start = $coords[1];
			my $end = $coords[2];
			
			if($chr eq $target_chr) {
				if($start <= $target_end && $end >= $target_start) {
					if(exists $overlap{$target}) {
						$overlap{$target} = $overlap{$target}.",".$gene;
					} else {
						$overlap{$target} = $gene;
					}
				} elsif ($end < $target_start ) {
					if(($target_start - $end) < $up_shortest_dist) {
						$up_shortest_dist = $target_start - $end;
						$upstream{$target} = $gene."_".$up_shortest_dist;
					}
				} elsif($start > $target_end) {
					if(($start - $target_end) < $down_shortest_dist) {
						$down_shortest_dist = $start - $target_end;
						$downstream{$target} = $gene."_".$down_shortest_dist;
					}
				}
			} else {
				next;
			}
		}
	}
	
}

#output

my $output = $target_gene_file."neighburing_genes.txt";
open OUTPUT, ">$output" or die "Can not open the file $output.\n";
print OUTPUT "Target_ID"."\t"."upstream_neighbour"."\t"."upstream_neighbour_symbol"."\t"."upstream_distance"."\t"."downstream_neighbour"."\t"."downstream_neighbour_symbol"."\t"."downstream_distance"."\t"."overlapped_genes"."\t"."overlapped_genes_symbol"."\n";
foreach my $target (keys %targets) {
	if (exists $upstream{$target}) {
		my @words = split("_", $upstream{$target});
		my $symbol;
		if(exists $symbols{$words[0]}) {
			$symbol = $symbols{$words[0]};
		} else {
			$symbol = "NA";
		}
		print OUTPUT $target."\t".$words[0]."\t".$symbol."\t".$words[1]."\t";
	} else {
		print OUTPUT $target."\t"."NA"."\t"."NA"."\t"."NA"."\t";
	}
	
	if (exists $downstream{$target}) {
		my @words = split("_", $downstream{$target});
		my $symbol;
		if(exists $symbols{$words[0]}) {
			$symbol = $symbols{$words[0]};
		} else {
			$symbol = "NA";
		}
		print OUTPUT $words[0]."\t".$symbol."\t".$words[1]."\t";
	} else {
		print OUTPUT "NA"."\t"."NA"."\t"."NA"."\t";
	}
	
	if (exists $overlap{$target}) {
		my @words = split(",",$overlap{$target} );
		my $symbol="";
		foreach my $gene (@words) {
			if(exists $symbols{$gene}) {
				if ($symbol eq "") {
					$symbol = $symbols{$gene};
				} else {
					$symbol = $symbol.",".$symbols{$gene};
				}
			} else {
				if ($symbol eq "") {
					$symbol = "NA";
				} else {
					$symbol = $symbol.","."NA";
				}
			}
		}
		print OUTPUT $overlap{$target}."\t".$symbol."\n";
	} else {
		print OUTPUT "NA"."\t"."NA"."\n";
	}
}

close(OUTPUT);

print "\n";
print "#\n";
print "Searching is done, results is in file $output.\n";
print "#\n";
exit;
