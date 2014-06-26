#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $reference_file="/users/a5907529/lustre/Yaobo/GenomeData/GATK_bundle/ucsc.hg19.4GATK.fasta";
my $clip_bed_file; #Bed file format of clip sites
my $motif_file; # list of motifs on binding sites
my $exon_list_file; # list of exons with the format : 'Exon_id	Chromosome	Strand(-1 and 1)	Start_sties	End_sites'
my $output_file= "motif.counted.clip.checked.txt";
my $Results=GetOptions( "exon_list_file=s" => \$exon_list_file, "motif_file=s" => \$motif_file,  "clip_bed_file=s" => \$clip_bed_file, "output=s" => \$output_file);

my $check_clip_sites;
if ($clip_bed_file) { $check_clip_sites = "Y";}
my $check_motifs;
if ($motif_file) { $check_motifs = "Y";}

#read in the hg19 reference file
my %hg19;

print "Start to read in reference file.\n";
open HG19, "$reference_file" or die "Can not open the file $reference_file\n";
my $string = "";
my $current_chr = "";
while (my $line = <HG19> ) {
	chomp $line;
	if ($line =~ m/^\#/) {
		next;
	} elsif( $line =~ m/^\>/) {
		$line =~ s/^\>//;
		if ($current_chr eq "") {
			$current_chr = $line;
			print "Reading in $current_chr...\n";
		} else {
			$hg19{$current_chr} = $string;
			$string = "";
			$current_chr = $line;
			print "Reading in $current_chr...\n";
		}
	} else {
		$string = $string.$line;
	}
}
close HG19;
#put in the last sting into the hash
$hg19{$current_chr} = $string;
print "Reference is imported\n";

#read in clip sites;
my %clip_sites_temp;
my %clip_sites;
if ($check_clip_sites eq "Y") {
	open CLIPS, $clip_bed_file or die "Cannot open the motif file $clip_bed_file\n";
	while (my $line = <CLIPS>) {
		next if ($line =~ m/^\#/);
		if ($line =~ m/^chr/) {
			chomp $line;
			my @temp = split ("\t", $line);
			my $chr = $temp[0];
			my $pos = $temp[2]; # bed file end position is 1-based
			$clip_sites_temp{$chr}{$pos} = 0;
		}
	}
	close CLIPS;
	foreach my $chr (keys %clip_sites_temp) {
		my @temp_1;
		foreach my $pos ( sort {$a <=> $b} keys %{$clip_sites_temp{$chr}} ) {
			push @temp_1, $pos;
		}
		$clip_sites{$chr} = \@temp_1;
	}
	undef %clip_sites_temp;
}


#read in motifs;
my @motifs;
if ($check_motifs eq "Y") {
	open MOTIF, $motif_file or die "Cannot open the motif file $motif_file\n";
	while (my $line = <MOTIF>) {
		next if ($line =~ m/^\#/);
		chomp $line;
		push @motifs, $line;
	}
}
close MOTIF;

#read in the exon file
open OUTPUT, ">$output_file" or die "Cannot open the file to output $output_file\n";
open EXONS, $exon_list_file or die "Cannot open the fie $exon_list_file\n";
my $processed_count = 0;
while (my $line = <EXONS>) {
	next if ($line =~ m/^\#/);
	$processed_count += 1;
	chomp $line;
	my @temp = split ("\t", $line);
	my $chr = $temp[1];
	my $strand = $temp[2];
	my $start = $temp[3];
	my $end = $temp[4];
	if ( $chr !~ m/^chr/ || $start !~ m/^\d+$/ || $end !~ m/^\d+$/  ) {
		print "Warning: Cannot recognize the format of exon on line:\n$line\n";
		next;
	}
	#if ($strand eq "-1") {
	#	my $temp = $start;
	#	$start = $end;
	#	$end = $temp;
	#}
	my $exon_length = $end - $start + 1;
	#get the exon sequence
	my $exon_seq = uc(substr($hg19{$chr}, $start-1, $exon_length));
	print OUTPUT $line;
	
	if (%clip_sites) {
		my $clip_count = 0;
		my @clip_sites_array = @{ $clip_sites{$chr} };
		#print "on $chr, there are total ", scalar @clip_sites_array, " clip sites\n";
		CLIP_LOOP: foreach my $clip_pos (@clip_sites_array) {
			if ( $clip_pos >= $start && $clip_pos <= $end) {
				#print "Clip found: $chr $clip_pos\n";
				$clip_count +=1;
			} elsif ( $clip_pos > $end ) {
				last CLIP_LOOP;
			}
		}
		print OUTPUT "\t".$clip_count;
	}
	
	if (@motifs) { #search and count motifs
		my $motif_count = 0;
		my $motif_total_length = 0;
		
		my @starts;
		my @ends;
		my @to_use_motifs = @motifs;
		foreach my $motif ( @to_use_motifs ) {
			$motif = uc($motif);
			if ($strand eq "-1") {
				$motif =~ tr/[A,T,G,C,]/[T,A,C,G]/;
				$motif = scalar reverse("$motif");
			}
			#print ">exon_seq: \n$exon_seq\n";
			#print ">motif: \n$motif\n";
			while ($exon_seq =~ m/(?=($motif))/g) {
				push @starts, $-[1];
				#print "Start: $-[1]";
				push @ends, $+[1];
				#print " ends: $+[1]\n";
				$motif_count += 1;
			}
		}
		
		if ($motif_count == 1) { 
			#print "1 la haha $motif_total_length;\n";
			$motif_total_length = $ends[0] - $starts[0];
		} elsif ( $motif_count > 1 ) {
			#resolve overlap motif matches
			my @idx = sort { $starts[$a] <=> $starts[$b] } 0 .. $#starts;
			@starts = @starts[@idx];
			@ends = @ends[@idx];
			my @new_starts;
			my @new_ends;
			push @new_starts, $starts[0];
			push @new_ends, $ends[0];
			for (my $i=1;$i<scalar @starts;$i++) {
				if ($starts[$i] >= $new_starts[-1] && $starts[$i] <= $new_ends[-1]) {
					if ($ends[$i] >= $new_ends[-1] ) {
						$new_ends[-1] = $ends[$i];
					}
				} elsif ($starts[$i] > $new_ends[-1] ) {
					push @new_starts, $starts[$i];
					push @new_ends, $ends[$i]; 
				}
			}
			
			#counting length of overlaps
			for (my $i=0;$i<scalar @new_starts;$i++) {
				#print "New start: $new_starts[$i]";
				#print "New end: $new_ends[$i]\n";
				#print "length of motif bases: $motif_total_length\n";
				$motif_total_length = $motif_total_length + $new_ends[$i] - $new_starts[$i];
			}
		}
		
		print OUTPUT "\t".$motif_count."\t".$motif_total_length;
	}
	print OUTPUT "\n";
	if ($processed_count % 1000 == 0 ) {
		print "Processed ".$processed_count." exons.\n";
	}
}
close EXONS;
close OUTPUT;
print "Done!\n";

exit;
