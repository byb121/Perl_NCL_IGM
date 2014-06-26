#!/usr/bin/perl
use warnings;


open FILE, "<alignment.txt" or die $!;
my $count = 0;
my @chromosome; #to store which chromosome is the alignment lies
my @position; #to store starting position of the read
my @direction; #to store read direction
my @filter; #to store if a read passes the filter

##to construct the genome
my %human_genome = ();
my $chr = "";
my $seq = "";

open GM_FILE, "</home/data/AllChromosomes/hg18.fa" or die $!;
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

#while (($key, $value) = each(%human_genome)){ 
#        $length = length ($value);
#        print "$key in length $length\n";
#}

##Construt a data structure to store count number

my %genome_count = (); #A better/complex structure can be constructed to store both cound and sequence. it could be done in the future!

while (($chr, $seq) = each %human_genome) {
	$length = length($seq);
	%count_ref = ();
	for($i = 1; $i <= $length; $i++) {
		$count_ref -> {$i} = 0;
	}
	$genome_count{$chr} = $count;
}

while (($key, $value) = each %genome_count){ 
        $length = keys($value);
        print "$key in length $length\n";
}


## Start to read the line
while(my $line = <FILE>) {
	chomp $line;
	@elements = split(/\s+/, $line); #split the line with multi spaces
	if($elements[10] =~ /^chr/) { 
		#store the values
		$chromosome[$count] = $elements[10];
                $position[$count] = $elements[11];
                $direction[$count] = $elements[12];
                $filter[$count] = $elements[15];
		#print "$count: the read is at $chromosome[$count] from $position[$count], direction is $direction[$count]. pass the filter $filter[$count]\n";
		
		if ($direction =~ /F/) {
			for (i=0; i<75; i++) {
			$genome_count{$chromosome}{$position+i} += 1;
			}
		} else {
			for (i=75; i>=0; i--) {
			$genome_count{$chromosome}{$position-i} += 1;
			}
		}
		#############
		#to make names of chromosome in different file same
		$chromosome[$count] = $split(/\./,$chromosome[$count])[0];
		                #read in the scarlar of the chromsome
		$count = $count + 1;
		#read in the scarlar of the chromsome
		#add 1 read to the corespond position on the chromosome

	}
}

close(FILE);

##Start to count each base aligned.
