#!/usr/bin/perl
use strict;
use warnings;
use lib qw(/users/a5907529/lib/5.16.1
            /users/a5907529/lib/site_perl/5.16.1);
use Statistics::Descriptive;


my ($input, $output) = @ARGV;

print "#####################################################\n";
print "# The script is to summarize the output of coverageBed.   -Yaobo #\n";
print "#####################################################\n";
print "\n";
print " reading from $input\n";

my %mean_count_on_exons;

open INPUT, $input or die "mean_coverage_from_coverageBed.pl: cannot open file $input\n";

while (my $line =<INPUT>) {
	if ($line =~ m/^chr/) {
		chomp $line;
		my @words = split(/\t/,$line);
		if ( exists $mean_count_on_exons{$words[3]} ) {
			$mean_count_on_exons{$words[3]} = $mean_count_on_exons{$words[3]} + $words[6]*$words[9];
		} else {
			$mean_count_on_exons{$words[3]} = $words[6]*$words[9];
		}
	}
}
close INPUT;

open OUTPUT, ">$output" or die "mean_coverage_from_coverageBed.pl: cannont open $output\n";

print OUTPUT "Total number of targeted exons: ";
my $hash_length = keys %mean_count_on_exons;
print OUTPUT $hash_length."\n";

my $count_eq0=0;
my $count_10x=0;
my $count_20x=0;
my $count_30x=0;
my $count_50x=0;
my $count_100x=0;
my $count_200x=0;
my $count_500x=0;

foreach my $exon (keys %mean_count_on_exons) {
	if ( $mean_count_on_exons{$exon} == 0 ) {$count_eq0 += 1;}
	if ( $mean_count_on_exons{$exon} >= 10 ) {$count_10x += 1;}
	if ( $mean_count_on_exons{$exon} >= 20 ) {$count_20x += 1;}
	if ( $mean_count_on_exons{$exon} >= 30 ) {$count_30x += 1;}
	if ( $mean_count_on_exons{$exon} >= 50 ) {$count_50x += 1;}
	if ( $mean_count_on_exons{$exon} >= 100 ) {$count_100x += 1;}
	if ( $mean_count_on_exons{$exon} >= 200 ) {$count_200x += 1;}
	if ( $mean_count_on_exons{$exon} >= 500 ) {$count_500x += 1;}
}


print OUTPUT "Number of exons has zero coverage: ".sprintf( "%.2f", $count_eq0/$hash_length*100)." (".$count_eq0.")"."\n";
print OUTPUT "Number of exons has more than 10x coverage: ".sprintf( "%.2f", $count_10x/$hash_length*100)." (".$count_10x.")"."\n";;
print OUTPUT "Number of exons has more than 20x coverage: ".sprintf( "%.2f", $count_20x/$hash_length*100)." (".$count_20x.")"."\n";;
print OUTPUT "Number of exons has more than 30x coverage: ".sprintf( "%.2f", $count_30x/$hash_length*100)." (".$count_30x.")"."\n";;
print OUTPUT "Number of exons has more than 50x coverage: ".sprintf( "%.2f", $count_50x/$hash_length*100)." (".$count_50x.")"."\n";;
print OUTPUT "Number of exons has more than 100x coverage: ".sprintf( "%.2f", $count_100x/$hash_length*100)." (".$count_100x.")"."\n";;
print OUTPUT "Number of exons has more than 200x coverage: ".sprintf( "%.2f", $count_200x/$hash_length*100)." (".$count_200x.")"."\n";;
print OUTPUT "Number of exons has more than 500x coverage: ".sprintf( "%.2f", $count_500x/$hash_length*100)." (".$count_500x.")"."\n";;

print OUTPUT "Min\t1st quantils\tMedian\tMean\t3rd Quantile\tMax:\n";
my $quantiles_ref = QuantilesOfHash(\%mean_count_on_exons);
my @quantiles = @{ $quantiles_ref };
foreach my $quan (@quantiles) {
	print OUTPUT $quan."\t";
}
print OUTPUT "\n";

close OUTPUT;

print "Done!\n";

exit;

sub QuantilesOfHash {
	my ($hash_ref) = @_;
	my @quantiles;
	my %hash = %{$hash_ref};
	my @numbers = values %hash;
	my $stat = Statistics::Descriptive::Full->new( );
	$stat->add_data(@numbers);
	push @quantiles,  sprintf( "%.1f", $stat->quantile(0));
	push @quantiles,  sprintf( "%.1f", $stat->quantile(1));
	push @quantiles,  sprintf( "%.1f", $stat->quantile(2));
	push @quantiles,  sprintf( "%.1f", $stat->mean());
	push @quantiles,  sprintf( "%.1f", $stat->quantile(3));
	push @quantiles,  sprintf( "%.1f", $stat->quantile(4));
	return \@quantiles;
}	


