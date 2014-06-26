#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Spreadsheet::ParseExcel;

my $exomiserXLSfile = "/home/yaobo/Downloads/test.xls";

my $ArrayOfExomiserXLS_ref;
my %HashOfExomiserXLS;
my @ArrayOfExomiserXLS_Output;

if ($exomiserXLSfile ne "") {
	$ArrayOfExomiserXLS_ref = ExomiserXLS2Array ($exomiserXLSfile);
	foreach my $exomiserRecord (@{$ArrayOfExomiserXLS_ref}) {
		my @temp = split("\t", $exomiserRecord);
		$HashOfExomiserXLS{$temp[0]}{$temp[1]}{$temp[2]}{$temp[3]} = $exomiserRecord;
		print $HashOfExomiserXLS{$temp[0]}{$temp[1]}{$temp[2]}{$temp[3]}."\n";
	}
}



sub ExomiserXLS2Array {
	my ($xlsFile) = @_;
	my @output;
	#my $xlsFile = "/home/yaobo/Downloads/test.xls";
	my $parser = Spreadsheet::ParseExcel->new();
	my $workbook = $parser->parse($xlsFile);
	
	my @readinArray;
	
	die $parser->error(), ".\n" if ( !defined $workbook );
	
	# Following block is used to Iterate through all worksheets
	# in the workbook and print the worksheet content 
	
	my $i=1; #to read only the first worksheet
	for my $worksheet ( $workbook->worksheets() ) {
		# Find out the worksheet ranges
		my ( $row_min, $row_max ) = $worksheet->row_range();
		my ( $col_min, $col_max ) = $worksheet->col_range();
		my ($gene, $score_no_prior, $score_mouse_exists, $chr, $pos, $ref, $var, $nn_change, $pathpogenicity, $path_score, $mouse_pheno, $omim_anno);
		# $pos is as the wAnnovar format start position
		for my $row ( $row_min+1 .. $row_max ) {
			#print "Row $row:\n";
			for my $col ( $col_min .. $col_max ) {
				# Return the cell object at $row and $col
				my $cell = $worksheet->get_cell( $row, $col ); 
				if ($col == 0 ) {
					$gene = $cell->value();
				} elsif($col == 1) {
					$score_no_prior = $cell->value();
				} elsif($col == 2) {
					$score_mouse_exists = $cell->value();
				}elsif ($col ==3) {
					my @lines = split("\n", $cell->value());
					$nn_change = "";
					for (my $j=0;$j<scalar @lines;$j++) {
						if ($j == 0) {					
							my @temp = split(":", $lines[$j]);
							$chr = $temp[0];
							my @temp2 = split(/\./, $temp[1]);
							if ($temp2[1] =~ m/(\d+)(\w+)\>(\w+)\s\[(.*)\]/) {
								$pos = $1;
								$ref = "$2";
								$var = "$3";
							} elsif ($temp2[1] =~ m/(\d+)(\w+)\>(\-)\s\[(.*)\]/) { #deletion
								$pos = $1;
								$ref = "$2";
								$var = "$3";
							} elsif ($temp2[1] =~ m/(\d+)(\-)\>(\w+)\s\[(.*)\]/) { #insertion
								$pos = $1+1;  ###########################might be a bug of exomiser to record positions
								$ref = "$2";
								$var = "$3";
							} else {
								print "shit, a Ghost!!\n";
							}
						} elsif ($j>2) {
							chomp $lines[$j];
							$nn_change = $nn_change."\n".$lines[$j];
						}
					}
					$nn_change =~ s/^\n//;
				} elsif($col==4) {
					my @lines = split("\n", $cell->value());
					$pathpogenicity = $lines[1];
					if ($lines[2] =~ m/(\d\.\d+)$/) {
						$path_score = $1;
					}
				} elsif($col==5) {
					my @lines = split("\n", $cell->value());
					#print scalar @lines."\n";
					$mouse_pheno = $lines[0];
					$mouse_pheno =~ s/No\sOMIM\sdisease\sentry//;
					#chomp $mouse_pheno;
					if (scalar @lines >= 2) {
						$omim_anno="";
						for (my $j=1;$j<scalar @lines;$j++) {
							#chomp $lines[$j];
							if ($lines[$j] =~ m/^OMIM\:/){
								$lines[$j] =~ s/^OMIM\://;
								$omim_anno = $omim_anno."\n".$lines[$j];
								#chomp $omim_anno;
							}
							$lines[$j] =~ s/^OMIM\://;
							#chomp $omim_anno;
						}
						$omim_anno =~ s/^\n//;
					} else {
						$omim_anno="NA";
					}
				}
			}
			
			#print "$gene, $score_no_prior, $score_mouse_exists, $chr, $pos, $ref, $var, $nn_change, $pathpogenicity, $path_score, $mouse_pheno, $omim_anno\n";
			push @output, "$chr\t$pos\t$ref\t$var\t$gene\t$score_no_prior\t$score_mouse_exists\t$nn_change\t$pathpogenicity\t$path_score\t$mouse_pheno\t$omim_anno"
		}
		$i+=1;
		if ($i >=2 ){last;}
	}
	return \@output;
}



exit;