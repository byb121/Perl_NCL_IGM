#!/usr/bin/perl
use strict;
use warnings;

sub readGenomeFile { #this will only work with a single file contain all chromosomes in fasta format

	my %human_genome = ();
	my $chr = "";
	my $seq = "";
	
#my $file = @_;
	my $file = $_[0];

	
	open GM_FILE, "<$file" or die $!;
	while(my $line = <GM_FILE>) {
	        chomp $line;
	        if($line =~ /^>chr/) {
	                if ($seq =~ /^$/) {
	                        $chr = $line;
	                        $chr = substr($chr,1);
	                        print "Start to process $chr!\n";
	                } else {
	                        #$length = length($seq);
	                        $human_genome{$chr} = $seq;
	                        print "$chr is done!!\n";
	                        $chr = $line;
	                        $chr = substr($chr,1);
	                        print "Start to process $chr!\n";
	                        $seq = "";
	                }
	        } else {
	                $seq = $seq.$line;
	        }
	}
	$human_genome{$chr} = $seq;
	print "$chr is done!!\n";
	
	close(GM_FILE);
	return %human_genome;
}


sub countRead {
	
	my $human_genome = @_;
	my %genome_count = (); #A better/complex structure can be constructed to store both cound and sequence. it could be done in the future!
	
	while (my ($chr, $seq) = each %$human_genome){
	       my $length = length($seq);
		print "constracting zeros for $chr.....\n";
	        for(my $i = 1; $i <= $length; $i++) {
	                $genome_count{$chr}{$i} = 0;
	                if ($i == 0.1*($length - $length%10)) {
	                	print "10% is done!\n";	
	                } elsif ($i == 0.2*($length - $length%10)) {
	                	print "20% is done!\n";
	                } elsif ($i == 0.3*($length - $length%10)) {
	                	print "30% is done!\n";
	                } elsif ($i == 0.4*($length - $length%10)) {
	                	print "40% is done!\n";
	                } elsif ($i == 0.5*($length - $length%10)) {
	                	print "50% is done!\n";
	                } elsif ($i == 0.6*($length - $length%10)) {
	                	print "60% is done!\n";
	                } elsif ($i == 0.7*($length - $length%10)) {
	                	print "70% is done!\n";
	                } elsif ($i == 0.8*($length - $length%10)) {
	                	print "80% is done!\n";
	                } elsif ($i == 0.9*($length - $length%10)) {
	                	print "90% is done!\n";
	                } elsif ($i == $length) {
	                	print "100% is done!\n";
	                }
	        }
		print "$chr is done!\n";
	}
	
	return (\%genome_count);
	
#	while (my ($key, $value) = each %genome_count) {
#	        my $length = keys %$value;
#		print "calculating the length of $key.\n";
#	        print "$key in length $length!\n";
#	}

}
1;
