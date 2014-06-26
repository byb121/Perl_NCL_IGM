#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

###############################################
#
# slightly amended from Helen's
# only change file handler to abs file path
#
###############################################


#Program to append full gene descriptions and OMIM disease information to Annovar output
#Also add full gene names for genes that do not feature in OMIM file

my $path1="/users/a5907529/lustre/scripts/AfterDindel/Yaobo_pipeline/test/";
my $fullname_file="/users/a5907529/lustre/scripts/AfterDindel/Yaobo_pipeline/Ensembl_OMIM_AllGeneNames.txt";
my $file="Annovar_Out_TEST_CovThresh_5_onTarget_variant_list_inhouseMAF.txt";
my $ID="Test_ID";

my $Results=GetOptions("inPath=s"=>\$path1, "inFile=s"=>\$file, "ensemblFile=s"=>\$fullname_file, "projectID=s"=>\$ID);

$path1 =~ s/\/$//;
my $outfile1=$path1."/"."Annovar_Out_".$ID."_onTarget_variant_list_inhouseMAF_OMIM.txt";

open(OUT1, ">$outfile1") || die "Cannot open file \"$outfile1\" to write to!\n";
open INPUT, $path1."/".$file or die "Cannot open ".$path1."/".$file."\n";
my $headF1=<INPUT>;
my $headF2=<INPUT>;
my $head=<INPUT>;
chomp $head;
print OUT1 "$head";
print OUT1 "\tChromosome\tPosition\tReference \(R\)\tVariant \(V\)\tGene\tMito Gene\?\tInHouse MAF\tGenotypes\tGenotype MAF\tAlleles Covered\tWikiGene Description\tAssociated GO Terms\tOMIM Disorder\n";

#Get full genenames into hash
my %full_genenames;
open INPUT_F, $fullname_file or die "Cannot open $fullname_file\n";
	my $head_e =<INPUT_F>;
	loopf: while (my $Line_f=<INPUT_F>){
		   chomp $Line_f;
		   my @linesplit_f = split(/\t/,$Line_f);
		   my $ensg=$linesplit_f[0];
		   my $goT="";
		   if(defined $linesplit_f[2]){$goT=$linesplit_f[2];}
		   my $wikiD="";
		   if(defined $linesplit_f[3]){$wikiD=$linesplit_f[3];}
		   my $mimD="";
		   if(defined $linesplit_f[5]){$mimD=$linesplit_f[5];}
		   if(!exists $full_genenames{$ensg}{'wiki'}{$wikiD}){$full_genenames{$ensg}{'wiki'}{$wikiD}=0;}
		   if(!exists $full_genenames{$ensg}{'go'}{$goT}){$full_genenames{$ensg}{'go'}{$goT}=0;}
		   if(!exists $full_genenames{$ensg}{'mim'}{$mimD}){$full_genenames{$ensg}{'mim'}{$mimD}=0;}
}
close INPUT_F;

my %combined_names;
foreach my $en (keys %full_genenames){
	foreach my $w (keys %{$full_genenames{$en}{'wiki'}}){
		if(exists $combined_names{$en}{'wiki'}){$combined_names{$en}{'wiki'}=$combined_names{$en}{'wiki'}.", ".$w;}
		if(!exists $combined_names{$en}{'wiki'}){$combined_names{$en}{'wiki'}=$w;}
		}
	foreach my $g (keys %{$full_genenames{$en}{'go'}}){
		if(exists $combined_names{$en}{'go'}){$combined_names{$en}{'go'}=$combined_names{$en}{'go'}.", ".$g;}
		if(!exists $combined_names{$en}{'go'}){$combined_names{$en}{'go'}=$g;}
		}
	foreach my $m (keys %{$full_genenames{$en}{'mim'}}){
		if(exists $combined_names{$en}{'mim'}){$combined_names{$en}{'mim'}=$combined_names{$en}{'mim'}.", ".$m;}
		if(!exists $combined_names{$en}{'mim'}){$combined_names{$en}{'mim'}=$m;}
		}
}


loop: while (my $Line=<INPUT>){
	chomp $Line;
	if($Line=~/chr/){
		my @linesplit1 = split(/\t/,$Line);
		my @ENSGs=();
		my @temp_ens=();
		if($linesplit1[1]=~/\;/){	#print "$linesplit1[1]\n";
									@temp_ens = split(/\;/,$linesplit1[1]);
									my %seen=();
									foreach my $t (@temp_ens){
										if(!exists $seen{$t}){
										$seen{$t}=0;
										if($t!~/[\,\;\(\)]/){	#print "$t\n";
																push(@ENSGs,$t);}
										my @t2=();
										if($t=~/\(/){@t2=split(/\(/,$t); 
													 if(!exists $seen{$t2[0]}){	push(@ENSGs,$t2[0]);
																				$seen{$t2[0]}=0;
																				#print "$t2[0]\n";
																				}
													}
										}#if not prev. seen
									}#foreach temp-ensg
								}#split ";"
		if($linesplit1[1]=~/\(/ and $linesplit1[1]!~/\;/){	@temp_ens = split(/\(/,$linesplit1[1]);						
															push(@ENSGs,$temp_ens[0]);
															#print "$temp_ens[0]\n";
															}#split "(" not ";"
								
		##coded for all ENSG split possibilities??!!
		if($linesplit1[1]!~/[\,\;\(\)]/){push(@ENSGs,$linesplit1[1]); 
											#print "$linesplit1[1]\n";
										}
				
		print OUT1 "$Line";
				
		##Print out gene info for each ENSG
		printloope: foreach my $ens (@ENSGs){
			if(!exists $combined_names{$ens}){
												print "$ens\n";
												next printloope;}
				print OUT1 "\t";
				print OUT1 $combined_names{$ens}{'wiki'};
				print OUT1 "\t";
				print OUT1 $combined_names{$ens}{'go'};
				print OUT1 "\t";
				print OUT1 $combined_names{$ens}{'mim'};
			}#foreach annovar comma sep gene
			print OUT1 "\n";
	
	}#if matches chr
}#while file line loop
close INPUT;
close OUT1;

exit;
