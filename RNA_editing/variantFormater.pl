#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
#use DBI;
#use DBD::mysql;
use LWP::Simple;
use Getopt::Long;

my ($input_file) = @ARGV;

my $outputfile = VariantFormater($input_file);

exit;

sub VariantFormater {
	
	my ($input_file)=@_;
	my @output;
	
	print "\n";
	print "############################      intro      #########################\n";
	print "The script expect variants in input file in format (tab-delimited):.\n";
	print "chr	pos	ref/var.\n";
	print "\n";
	print "Input data is from: $input_file.\n";
	print "The script will change the format to (tab-delimited): .\n";
	print "chr(a single number without \"chr\")	pos	pos	ref	var\n";
	print "\n";
	print "Formating................ ";
	

	
	open VARIANTS, $input_file or die "Cannot open $input_file"; 
	while (my $line = <VARIANTS>){  
		chomp $line;
		if ($line =~ m/^chr/) {
			my @SplitLine=split(/\t/, $line);
			my $Chr=$SplitLine[0];
			$Chr =~ s/chr//;
			my $Pos=$SplitLine[1];
			my $Mut=$SplitLine[2];
			my @bases = split("/", $Mut);
			
			my $newline = $Chr."\t".$Pos."\t".$Pos."\t".$bases[0]."\t".$bases[1];
			for (my $i=3;$i<=$#SplitLine;$i++) {
				$newline = $newline."\t".$SplitLine[$i];
			}
			push @output, $newline."\n";
		}
	}
	close VARIANTS;
	
	my $output_file = $input_file."_formated.txt";
	
	open OUTPUT, ">$output_file" or die "Cannont open $output_file";
	print OUTPUT @output;
	close OUTPUT;
	print "Done!\n";
	print "Output file is: $output_file\n";
	
	return $output_file;
	
}

