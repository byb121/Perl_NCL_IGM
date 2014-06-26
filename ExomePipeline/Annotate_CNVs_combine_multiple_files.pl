#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

#Program to combine Annotated CNV files

my $outputfile = "CNVs_combined_AROS_April2013_annotated_2nd_18.txt";
my $CurrentDir= "/users/nhrg/lustre/AROS/April2013_PFC_Samples_Results/CNVs";
chomp $CurrentDir;

my $Results=GetOptions("inpath=s"=>\$CurrentDir, "outfile=s"=>\$outputfile);

my @DirContent=`ls $CurrentDir`;
my $CNVfilecount=0;
my $head="head";

my $outputfile2 = $CurrentDir."/".$outputfile;
chomp $outputfile2;
unless(open(OUT2,">$outputfile2")){print "Cannot open file \"$outputfile2\" to write to!\n"; exit;}

my %var_pat_seen=();
my %Variants=();
dircontloop: foreach my $file (@DirContent){
	chomp $file;
	unless($file =~ /NONE_CNV_\S+.csv/){
		if($file =~ /^Annotated_CNV_(\S+).csv2$/){
			my $id=$1;
			#print "$file\n";
			$CNVfilecount++;
			my $file_path = $CurrentDir."/".$file;
			open INPUT, $file_path or die "Cannot open $file_path\n";
			$head=<INPUT>;
			foreach my $Line (<INPUT>){
				chomp $Line;
				my @temp = split (/\t/,$Line);
				my @elements;
				my $temp_string="";
				my $extension_switch = "OFF";
				foreach my $ele (@temp) {
					chomp $ele;
					if ($extension_switch eq "OFF") {
						if($ele =~ m/^".*"$/) {
							push @elements, $ele;
						} elsif ($ele =~ m/^"/) {
							$extension_switch = "ON";
							$temp_string=$ele;
						} else {
							push @elements, $ele;
						}
					} else {
						if ($ele =~ m/"$/) {
							$extension_switch = "OFF";
							push @elements, $temp_string.", ".$ele;
						 } else {
							$temp_string=$temp_string.", ".$ele;
						}
					}
				}
				if (scalar @elements < 20) {
					my $length=scalar @elements;
					for (my $hh = 0;$hh<(20-$length);$hh++) {
						push @elements, "NA";
					}
				}
				
				print OUT2 "$id";
				foreach my $temp_temp (@elements) {
					$temp_temp =~ s/\"//g;
					print OUT2 '@'."$temp_temp"; 
				}		
				print OUT2 "\n";
			}
		close INPUT;
		}
	}

}
close OUT2;


exit; 
