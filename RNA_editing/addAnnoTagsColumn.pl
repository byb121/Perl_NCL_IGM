#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
#use DBI;
#use DBD::mysql;
use LWP::Simple;
use Getopt::Long;


my ($input_file) = @ARGV;

my $outputfile = AddAnnotationTags($input_file);

exit;


sub AddAnnotationTags {
	my ($input) = @_;
	
	print "\n";
	print "############################      intro      #########################\n";
	print "Two columns will be added to the input file.\n";
	print "Added Col 1: GeneticPosition (eg: CDS/exon/intron/intergenic).\n";
	print "Added Col 2: ExtraTags (eg: UTR/start_codon/stop_codon/Selenocysteine).\n";
	print "# Adding Annotated tags as a column to the file................ ";


	my @output;
	my %tags;
	my $output_file = $input."_with_annotated_tags.txt";

	my $dir = cwd();
	my $temp_output_file = $dir."/"."temptemptempAnnotatedtags_file.txt";
	my $temp_output_file_2 = $dir."/"."temptemptempAnnotatedtags_file_with_genes_Names.txt";
	
	open INPUT, $input or die "Cannot open $input";
	open OUTPUT, ">$temp_output_file" or die "Cannot open $temp_output_file";
	while (my $line = <INPUT>) {
		chomp $line;
		if ($line =~ m/^chr/){
			my @words = split("\t", $line);
			print OUTPUT $words[0]."\t".$words[1]."\t".$words[1]."\n";
		}
	} 
	close OUTPUT;
	close INPUT;
	
	`/users/a5907529/biosofts/BEDTools-Version-2.12.0/bin/intersectBed -wo -a $temp_output_file -b ~/GenomeData/hg19/gencode.v10.annotation.gtf > $temp_output_file_2`;

	open ANNO, $temp_output_file_2 or die "Cannot open $temp_output_file_2";
	while (my $line = <ANNO>) {
		chomp $line;
		if( $line =~ m/^chr/ ){
			my @words = split("\t", $line);
			my $key = $words[0]."_".$words[1];
			if (exists $tags{$key}) {
				my %temp_hash;
				my @anno_tags = split(",", $tags{$key});
				for my $value (@anno_tags) {
						$temp_hash{$value} = 1;
				}
				
				if (exists $temp_hash{$words[5]}) {
					next;
				} else {
					$tags{$key} = $tags{$key}.",".$words[5];
				}
			} else {
				$tags{$key} = $words[5];
			}
		}
	}
	close ANNO;
	
	open INPUT, $input or die "Cannot open $input";
	while (my $line = <INPUT>) {
		chomp $line;
		if ($line =~ m/^chr/){
			my @words = split("\t", $line);
			if (exists $tags{$words[0]."_".$words[1]}) {
				my $pos_tag;
				if ($tags{$words[0]."_".$words[1]} =~ m/CDS/ ) {
					$pos_tag = "CDS";
				} elsif ($tags{$words[0]."_".$words[1]} =~ m/exon/) {
					$pos_tag = "exon";
				} elsif ($tags{$words[0]."_".$words[1]} =~ m/transcript/) {
					$pos_tag = "intron";
				} elsif ($tags{$words[0]."_".$words[1]} =~ m/gene/) {
					print "\nWarning: $words[0]:$words[1] : Annotation file does not have a transcript record overlapping the position, "; 
					print "but an annotated gene was found for the posion! Possible false annotation! will consider it's on a intron position.\n";
					$pos_tag = "intron";
				}
				
				my $extra_tag;
				if ($tags{$words[0]."_".$words[1]} =~ m/UTR/ ) {
					if (length $extra_tag == 0) {
						$extra_tag = "UTR";
					} else {
						$extra_tag = $extra_tag.","."UTR";
					}
				} 
				
				if ($tags{$words[0]."_".$words[1]} =~ m/start\_codon/) {
					if (length $extra_tag == 0) {
						$extra_tag = "start_codon";
					} else {
						$extra_tag = $extra_tag.","."start_codon";
					}
				} 
				
				if ($tags{$words[0]."_".$words[1]} =~ m/stop\_codon/) {
					if (length $extra_tag == 0) {
						$extra_tag = "stop_codon";
					} else {
						$extra_tag = $extra_tag.","."stop_codon";
					}
				} 
				
				if ($tags{$words[0]."_".$words[1]} =~ m/Selenocysteine/) {
					if (length $extra_tag == 0) {
						$extra_tag = "Selenocysteine";
					} else {
						$extra_tag = $extra_tag.","."Selenocysteine";
					}
				}
				
				if (length $extra_tag == 0) {
					$extra_tag = "N/A";
				}

				push @output, $line."\t".$pos_tag."\t".$extra_tag."\n";
				
			} else {
				push @output, $line."\t"."Intergenic"."\t"."N/A"."\n";
			}
		}
	}
	close INPUT;
	
	open OUTPUT, ">$output_file" or die "Cannont open $output_file";
	print OUTPUT @output;
	close OUTPUT;
	print "Done!\n";
	print "Output file is: $output_file\n";

	#`rm $temp_output_file`;
	#`rm $temp_output_file_2`;

	return $output_file;
}
