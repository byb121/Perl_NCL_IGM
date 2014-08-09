#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(sum);


########################################################
###############set up parameteres#######################
########################################################

#read length
my $read_length = 75;
#alignment file. currently only *_export.txt file of CASAVA 1.6 are supported
my $alignment_file = "/home/a5907529/alignment_test_file.txt";
my $targets_file = "/home/a5907529/exons.txt";


################## to hand possible warnings ##########################
my @errors;


##### more information about the alignment #######################

my $reads_passed_filter =0;
my $reads_passed_filter_aligned_on_exons=0; ##### how many aligned read passed the filter and aligned onto exons
my $number_of_exon_have_at_least_one_read_aligned=0;
my $output_file ="/home/a5907529/coverage_newloop.output";


############### to import exome positions #########################

my ($hash_ref, $array_ref) = read_ccds($targets_file); #values are references of a hash and array
for my $key (keys %$hash_ref) {
        my $temp = scalar @{$hash_ref -> {$key}};
        print "$key contains $temp entries:\n";
}

my $array_size = @$array_ref; # the array which contains numbers of counts
print $array_size," is the size of the big array contain counts";
my @counts = @$array_ref;
print "\n";

####  Start to read the alignment file
my $n = 0; # line number
my $tmp_hash = {};       

open FILE, "<$alignment_file" or die $!;
while(my $line = <FILE>) {
	
	chomp $line;
	$n += 1;
	# print "reading line $n......\n";
	my @elements = split(/\s+/, $line); #split the line with multi spaces
	if($elements[10] =~ /^chr/) {  #to define values of each variable
		my @temp = split(/\./, $elements[10]); #to make names of chromosome in different files the same
		my $contig = $temp[0];
		my $position = $elements[11];
		my $direction = $elements[12];
		my $filter = $elements[@elements-1];
		
		# to display progress
		#if($n%5000 == 0){
		#	print "processed $n reads .....\n";
		#find all of exon on the chromosome
			
		#print "$count: the read is at $contig[$count] from $position[$count], direction is $direction[$count]. pass the filter $filter[$count]\n";
		if ($filter =~ /Y/) {
			
			if (exists $hash_ref->{$contig}){ #the read must pass the filter
				
				my $left;
				my $right;
				$reads_passed_filter +=1;
				
				if ($direction =~ /F/) {
					$left = $position;
				} elsif ($direction =~ /R/){
					$left = $position - $read_length + 1;
				}
			
				push @{$tmp_hash->{$contig}}, $left;
			
				if($n%100000 ==0 || eof){
				
					print "processing $n line of the file\n";
CONTIG:				for my $key (keys %$tmp_hash) {
						my @sorted_array = sort{$a <=> $b} @{$tmp_hash->{$key}};
						print "$key gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg\n";
					
EXON:					for(my $i = 0; $i < scalar @{$hash_ref->{$key}};$i++) {
							
							my $switch = 1; #in default the switch will be on
								
							if (scalar @sorted_array == 0) {
								# go the next $key if the array is empty
								last;
							}
							
							###### to deal with reads aligned at same position
						
							#if(scalar @sorted_array > 1) {
							#	for(my $x = 0; $x < (scalar @sorted_array-1); $x++){
							#		if ($sorted_array[$x] == $sorted_array[$x+1]){
							#			$repeats += 1;
							#		} else {
							#			last;
							#		}
							#	}
							#}
						
							#### to decide if the read is aligned on to an exon
							my @words = split(/_/, $hash_ref->{$key}[$i]);
							
							$left = $sorted_array[0];
							
							if ($words[1] < $left){
								next;
							}
							
							while (($left + $read_length - 1) < $words[0] && scalar @sorted_array > 1) {
								shift(@sorted_array);
								$left = $sorted_array[0];
							}
							
							while ($switch == 1){
								
								$left = $sorted_array[0];
								$right = $left + $read_length -1;
								
								if ($left > ($words[0]-$read_length+1) && $left <= $words[0]) {
									if($right <= $words[1] ) {
										print "Read: $n--- Situation 1: Reads: $left - $right on $key. Exon: $words[0] -$words[1] on $key.\n";
										for (my $j = $words[2]; $j <= ($words[2]+$right-$words[0]);$j++) { 
											#do something to the @count
											$counts[$j] += 1;
										}
										$reads_passed_filter_aligned_on_exons +=1;
										$words[3] +=1;
										$hash_ref->{$key}[$i] = $words[0]."_".$words[1]."_".$words[2]."_".$words[3];
										shift(@sorted_array);
									} elsif ($right > $words[1]) {
										print "Read: $n--- Situation 2:  Reads: $left - $right on $key. Exon: $words[0] -$words[1] on $key.\n";
										for (my $j = $words[2]; $j <=($words[2]+$words[1]-$words[0]);$j++) { 
											#do something to the @count
											$counts[$j] += 1;
										}
										$reads_passed_filter_aligned_on_exons +=1;
										$words[3] +=1;
										$hash_ref->{$key}[$i] = $words[0]."_".$words[1]."_".$words[2]."_".$words[3]; 
										shift(@sorted_array);                 
									}
								} elsif ($left > $words[0] && $left <= $words[1]) {
									#read starts inside the exon
									if($right <= $words[1] ) {
										# read is insdie the exon
										print "Read: $n--- Situation 3:  Reads: $left - $right on $key. Exon: $words[0] -$words[1] on $key.\n";
										for (my $j = ($words[2] + $left - $words[0]);$j <=($words[2]+$right-$words[0]);$j++) {
											#do something to the @count
											$counts[$j] += 1;
										}
										$reads_passed_filter_aligned_on_exons +=1;
										$words[3] +=1;
										$hash_ref->{$key}[$i] = $words[0]."_".$words[1]."_".$words[2]."_".$words[3];
										shift(@sorted_array);
									} elsif ($right > $words[1]){
										#read is shorter than the exon
										print "Read: $n--- Situation 4:  Reads: $left - $right on $key. Exon: $words[0] -$words[1] on $key.\n";
										for (my $j = ($words[2] + $left - $words[0]); $j <=($words[2]+$words[1]-$words[0]);$j++) { 
											#do something to the @count
											$counts[$j] += 1;
										}
										$reads_passed_filter_aligned_on_exons +=1;
										$words[3] +=1;
										$hash_ref->{$key}[$i] = $words[0]."_".$words[1]."_".$words[2]."_".$words[3];
										shift(@sorted_array);
									}
								} else {
									$switch = 0;
								}
							} 
						}
					}
					
					$tmp_hash = {};
				}
			} else {
				push (@errors, "WARNING: $contig does NOT exist in reference contigs!\n");
			}
		} 
	}
}
close(FILE);

print "Coverage on whole exome is ",mean(\@counts),"\n";
print "The run is finished with".@errors." error.\n";


###################write a report of sammury ######################

open INFO, ">>$output_file" or dir $!;

print INFO "###########################################\n";
print INFO "#####      RUN RESULT of COVERAGE     #####\n";
print INFO "###########################################\n";
print INFO "\n\n\n\n\n";
print INFO "alignment file is: $alignment_file\n";
print INFO "The size of the big array contain counts: ", scalar @counts."\n";
print INFO "exons position was read from: $targets_file\n";
print INFO "Errors and warnings:\n";

if (!@errors) {
        print INFO "The run is finished with no error.\n";
} else {
        foreach my $line (@errors) {
                print INFO "$line";
        }
}
print INFO "Coverage on whole exome is ",mean(\@counts),"\n";
print INFO "$n reads were processed.\n";
print INFO "$reads_passed_filter read passed the filter.\n";
print INFO "$reads_passed_filter_aligned_on_exons are aligned on exons.\n";

########### summary exon with or without reads aligned onto #################

my @landings;

for my $key (keys %$hash_ref) {
	foreach my $count (@{$hash_ref->{$key}}) {
		my @words = split(/_/, $count);
		push(@landings, $words[3])
	} 
}

my %landing_summary;
my %counts_summary;

foreach my $value (@landings) {
	if (exists $landing_summary{$value}) {
		$landing_summary{$value} += 1;
	}else{
		$landing_summary{$value} = 1;
	}
}

foreach my $value (@counts) {
        if (exists $counts_summary{$value}) {
                $counts_summary{$value}  += 1;
        }else{
                $counts_summary{$value} = 1;
        }
}

print INFO "How many exon	landings\n";
while (my ($key,$value) = each(%landing_summary)) {
	print INFO "$value	$key\n ";
}

print INFO "counts	frequency\n";
while (my ($key,$value) = each(%counts_summary)) {
        print INFO "$value	$key\n";
}


close(INFO);

exit;

sub mean {
        my ($array) = @_;
        return @$array ? sum(@$array) / @$array : 0;
}

sub read_ccds{
    my ($file) = @_; #current surport input file format is defined in sort_CCDS_entry.pl
    my @contig_names;
    my @begains;
    my @ends;

    my %coordinates; # hash of arrays to hold exon positions. key is the chromosome number. value is in "start_end_countIndex" format
    my @counts; #array to hold count numbers

    open FILE, "<$file" or die $!; #current surpor
    print "Start to read CCDS file.\n";
    while(my $line = <FILE>) {
        chomp $line;
        my @elements = split(/\s+/, $line); #split the line with multi spaces

        my $chr = $elements[0];
        #print "Gathering information of $chr....\n";

        ##this large section is just to handle overlaps between CCDSs

        if (exists $contig_names[0]) { # if this is not the first element for all of arrays
        	if ($contig_names[$#contig_names] ne $chr) { # if a new contig is read
        		push (@contig_names,$chr);
        		push (@begains,$elements[1]);
        		push (@ends,$elements[2]);
         #           print "Done!\n";
			} else { # when still on the same contig
				if ($elements[1] == $begains[$#contig_names]) {
                        # if the current exon starts at the same position as the last entry
					if ($elements[2] > $ends[$#contig_names]) {
                        # but this exon longer than the last one, change the ending position of last one
                       	$ends[$#ends] = $elements[2];
					}
				} elsif ($elements[1] < $ends[$#contig_names]) {
					# if the current exon starting position is inside the last exon
					if($elements[2] > $ends[$#contig_names]) {
						# and the exon ends outside of the last exon, it could be combine with last situation
						# change the ending position of the last exon
						$ends[$#ends] = $elements[2];
					}
				} elsif ($elements[1] == $ends[$#contig_names]) {
					# if the current exon starts at the ending position of last exon
					#change the ending position of last exon
					$ends[$#ends] = $elements[2];
				} elsif ($elements[1] > $ends[$#contig_names]) {
					# if the current exon is totally at the downstream of the last exon
					push(@contig_names, $chr);
					push(@begains, $elements[1]);
					push(@ends, $elements[2]);
				}
          #              print "Done!\n";
			}
        } else { 
        	# if this is the first entry for all of the arrays
        	push (@contig_names,$chr);
        	push (@begains,$elements[1]);
        	push (@ends,$elements[2]);
        	#     print "Done!\n";
		}
    }
    close(FILE);

    # to construct a hash of arrays. key is the contig name, value is "start_end_countIndex(_exonName)" format
    # count array is constructed at the same time. to store count numbers.
    print "Start to construct the position hash  and the count array.\n";
    for (my $i=0;$i<scalar @contig_names;$i++) {
        my $tmp_value = $begains[$i]."_".$ends[$i]."_".(scalar @counts)."_0";
        my $tmp_key = $contig_names[$i];
        push (@{$coordinates{$tmp_key}},$tmp_value);

        my $tmp_number = scalar @counts;
        for (my $j = $tmp_number; $j<=($tmp_number + $ends[$i] - $begains[$i]); $j++) {
                $counts[$j] = 0;
        }
        my $length = scalar @contig_names;

        ## this section is to get a percentage of the whole progress.
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
                ##############
    }
    print "Positions of exome on genome is recorded successfully!\n";
    return (\%coordinates,\@counts);
}         



                                             
