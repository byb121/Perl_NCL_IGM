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
my $ccds_file = "/home/a5907529/CCDS.txt";


################## to hand possible warnings ##########################
my @errors;


##### to output more information about the alignment #######################

my $reads_passed_filter =0;
my $reads_passed_filter_aligned_on_exons=0; ##### how many aligned read passed the filter and aligned onto exons
my $number_of_exon_have_at_least_one_read_aligned=0;
my $output_file ="/home/a5907529/coverage_hg18_perl.output";


############### to import exome positions #########################

#return (\@contig_names,\@begains,\@ends,\@index,\@landings,\@counts);
my ($contigs,$begains, $ends, $index,$landings, $counts) = read_ccds($ccds_file); #values are references of a hash and array

print "Coordinates array contains ", scalar @$contigs," entries:\n";
print "The size of the big array contain counts: ", scalar @$counts."\n";

print "\n";

####  Start to read the alignment file
my $n = 0; # line number
        
print "Start to read the alignment file...............\n";

open FILE, "<$alignment_file" or die $!;
while(my $line = <FILE>) {

	chomp $line;
        $n += 1;
	
	# to display progress
	if($n%5000 == 0){
		print "processed $n reads .....\n";
	} 
	# print "reading line $n......\n";
	my @elements = split(/\s+/, $line); #split the line with multi spaces
	if($elements[10] =~ /^chr/) { #to define values of each variable
    	my @temp = split(/\./, $elements[10]); #to make names of chromosome in different files the same
        my $contig = $temp[0];
        my $position = $elements[11];
        my $direction = $elements[12];
        my $filter = $elements[@elements-1];

        #find all of exon on the chromosome

        #print "$count: the read is at $contig[$count] from $position[$count], direction is $direction[$count]. pass the filter $filter[$count]\n";
        if ($filter =~ /Y/) { #the read must pass the filter
        	my $left;
            my $right;
            $reads_passed_filter +=1;
                    
            if ($direction =~ /F/) {
            	$left = $position;
                $right = $position + $read_length -1;
            } elsif ($direction =~ /R/){
            	$left = $position - $read_length + 1;
                $right = $position;
           	}
                    
            ###### start to search coodinates
            for(my $i = 0; $i < scalar @$contigs; $i++) {
               	if($contig eq $$contigs[$i]) {
					if ($left > ($$begains[$i]-$read_length+1) && $left <= $$begains[1]) { # read started u`psteam of the exon
        	       		if($right <= $$ends[$i] ) {
            	   			print "Read: $n--- Situation 1: the read: $left to $right on $contig. Exon: $$begains[$i] to $$ends[$i] on $$contigs[$i].\n";
                            for (my $j = $$index[3]; $j <= ($$index+$right-$$begains[1]);$j++) { #do something to the @count
                            	$$counts[$j] += 1;
                            }
                            $reads_passed_filter_aligned_on_exons +=1;
                            $$landings[$i] +=1;
                         } elsif ($right > $$ends[$i]) {
                         	print "Read: $n--- Situation 1: the read: $left to $right on $contig. Exon: $$begains[$i] to $$ends[$i] on $$contigs[$i].\n";
                            for (my $j = $$index[$i]; $j <=($$index[$i]+$$ends[$i]-$$begains[$i]);$j++) { #do something to the @count
                            	$$counts[$j] += 1;
                            }
                            $reads_passed_filter_aligned_on_exons +=1;
                            $$landings[$i] +=1;
                         }
                   	} elsif ($left > $$begains[$i] && $left <= $$ends[$i]) { #read starts inside the exon
                   		if($right <= $$ends[$i] ) { # read is insdie the exon
                        	print "Read: $n--- Situation 1: the read: $left to $right on $contig. Exon: $$begains[$i] to $$ends[$i] on $$contigs[$i].\n";
                            for (my $j = ($$index[$i] + $left - $$begains[$i]); $j <=($$index[$i]+$right-$$begains[$i]);$j++) { #do something to the @count
                            	$$counts[$j] += 1;
                            }
                            $reads_passed_filter_aligned_on_exons +=1;
                            $$landings[4] +=1;
                        } elsif ($right > $$ends[1]){ #right side of the read is outside the exon
                         	print "Read: $n--- Situation 1: the read: $left to $right on $contig. Exon: $$begains[$i] to $$ends[$i] on $$contigs[$i].\n";
                            for (my $j = ($$index[$i] + $left - $$begains[$i]); $j <= ($$index[$i]+$$ends[$i]-$$begains[$i]);$j++) { #do something to the @count
                            	$$counts[$j] += 1;
                            }
                            $reads_passed_filter_aligned_on_exons +=1;
                           	$$landings[$i] += 1;
                        }
                   	}
                #push (@errors, "WARNING: $contig does NOT exist in reference contigs!\n");
               	}
        	}
        }
	}
}
close(FILE);

print "Coverage on whole exome is ",mean($counts),"\n";
# print "The run is finished with".@errors." error.\n";


###################write a report of sammury ######################

open INFO, ">>$output_file" or dir $!;

print INFO "###########################################\n";
print INFO "#####      RUN RESULT of COVERAGE     #####\n";
print INFO "###########################################\n";
print INFO "alignment file is: $alignment_file\n";
print INFO "exons position was read from: $ccds_file\n";
#print INFO "Errors and warnings:\n";


print INFO "Coverage on whole exome is ",mean($counts),"\n";
print INFO "$n reads were processed.\n";
print INFO "$reads_passed_filter read passed the filter.\n";
print INFO "$reads_passed_filter_aligned_on_exons are aligned on exons.\n";

####################################################
##working on this part
####################################################


my %counts_summary;

foreach my $value (@$landings) {
	if (exists $counts_summary{$value}) {
		my $temp = $counts_summary{$value};
		$counts_summary{$value} = $temp + 1;
	}else{
		$counts_summary{$value} = 1;
	}
}

print INFO "How many exon		have	 	how many reads landed on\n";
while (my ($key,$value) = each(%counts_summary)) {
	print INFO "$value						$key\n ";
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
        } else { # if this is the first entry for all of the arrays
                push (@contig_names,$chr);
                push (@begains,$elements[1]);
                push (@ends,$elements[2]);
           #     print "Done!\n";
        }
        ################
    }
    close(FILE);

    my @counts; #array to hold count numbers
    my @index;
    my @landings;
    

    # to construct a hash of arrays. key is the contig name, value is "start_end_countIndex(_exonName)" format
    # count array is constructed at the same time. to store count numbers.
    print "Start to construct the position arrays and the count array.\n";
    for (my $i=0;$i<scalar @contig_names;$i++) {
        
        my $tmp_value =scalar @counts;
        push (@index,$tmp_value);
        push (@landings, 0);

        for (my $j = $tmp_value; $j<=($tmp_value + $ends[$i] - $begains[$i]); $j++) {
                $counts[$j] = 0;
        }


        ## this section is to get a percentage of the whole progress.
        my $length = scalar @contig_names;
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
    return (\@contig_names,\@begains,\@ends,\@index,\@landings,\@counts);
}  
                                         
