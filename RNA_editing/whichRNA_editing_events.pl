#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
#use DBI;
#use DBD::mysql;
use LWP::Simple;
use Getopt::Long;


my ($input_file) = @ARGV;

my $outputfile = DetermineRNAeditingType($input_file);

exit;

sub DetermineRNAeditingType {
	
	my ($input_file)=@_;
	my @output;
	
	print "\n";
	print "############################      intro      #########################\n";
	print "one column will be added to the input file.\n";
	print "Added Col 1: A-to-I/C-to-U.\n";
	print "\n";
	print "Input data is from: $input_file.\n";
	print "Possible RNA-editing events is define in the script as following:.\n";
	print "A-to-I: A/T, A/C, T/A, T/G.\n";
	print "C-to-U: C/T, G/A.\n";
	print "\n";
	print "# Adding Annotated tags as a column to the file................ ";
	

	
	open VARIANTS, $input_file or die "Cannot open $input_file"; 
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			my $Chr=$SplitLine[0];
			my $Pos=$SplitLine[1];
			my $Mut=$SplitLine[2];
			my $ComplimentMutant=$Mut;
			$ComplimentMutant=~tr/ACGT/TGCA/;
			
			#Determine if the Mut type is an possible RNA-editing event.
			if ($Mut =~ m/^[Aa]\/[Tt]|^[Aa]\/[Cc]|^[Tt]\/[Aa]|^[Tt]\/[Gg]/ | $ComplimentMutant =~ m/^[Aa]\/[Tt]|^[Aa]\/[Cc]|^[Tt]\/[Aa]|^[Tt]\/[Gg]/) {
				push @output, $line."\t"."A-to-I"."\n";
			}elsif ($Mut =~ m/^[Cc]\/[Tt]|^[Gg]\/[Aa]/ | $ComplimentMutant =~ m/^[Cc]\/[Tt]|^[Gg]\/[Aa]/) {
				push @output, $line."\t"."C-to-U"."\n";
			} else {
				print "Warning: the input file contains an non A to I/C to U change on postion on $Chr:$Pos. Check the input!\n";
				exit;
			}
		}
	}
	close VARIANTS;
	
	my $output_file = $input_file."_with_AI_CU_decisions.txt";
	
	open OUTPUT, ">$output_file" or die "Cannont open $output_file";
	print OUTPUT @output;
	close OUTPUT;
	print "Done!\n";
	print "Output file is: $output_file\n";
	
	return $output_file;
	
}

