#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $help;
my $sample_table;
my $MitoPi_dir="/users/a5907529/data/MK_MiSeq/mito_pipe_20150703";
my $Annovar_dir="/users/a5907529/data/Files_HG/vcf_annotation_november2013";
my $SEQ_PLATFORM="Illumina";
my $Library_ID_prefix="MK.YX.MitoPi";

if ( @ARGV < 1 || ! GetOptions('help|?' => \$help, "in=s"=>\$sample_table, "MitoPi=s"=>\$MitoPi_dir, 
				'Annova=s' => \$Annovar_dir, 'platform=s' => \$SEQ_PLATFORM, 'LibraryID=s' => \$Library_ID_prefix) 
				|| defined $help ) {
	usage(); exit;
}

$MitoPi_dir =~ s/\/$//;
$Annovar_dir =~ s/\/$//;
my $ref = $MitoPi_dir."/ref/GRCh38.p3.fa";

## 1. sample table in format -> sample_ID fastq_1 fastq_2 output_directory
#open the table
#check input exists
#check if fastqs exit
# check if output dir exist, create if not
# check if lane of a sample duplicated
# read in essential parameters
my %in_paras;
my %sample_output;
if ( -e $sample_table) {
	print "Found input sample table file $sample_table\n";
	open INPUT, "$sample_table" or die "Can not open the file $sample_table to read.\n";
	my %sampleLanes;
	 # used to warn different output place for a same sample
	while (my $line = <INPUT>) {
		chomp $line;
		if ($line =~ /^\#/) { next;} else {
			my @temp = split("\t", $line);
			if (scalar @temp < 5 ) { print "Error: line in $sample_table does not have 5 elements or not tab-delimited.\n"; exit;}
			if ( ! -e $temp[1]) { print "Error: File $temp[1] does not exist.\n"; exit;}
			if ( ! -e $temp[2]) { print "Error: File $temp[1] does not exist.\n"; exit;}
			if ( ! -e $temp[4]) {print "Make the output folder $temp[4]\n"; `mkdir -p $temp[4]`;}
			
			if (! exists $sampleLanes{$temp[0]."\t".$temp[3]}) {
				$sampleLanes{$temp[0]."\t".$temp[3]} = 1;
			} else {
				print "Error: $temp[0] has duplicated lanes.\n"; exit;
			}
			
			if (! exists $sample_output{$temp[0]}) {
				$sample_output{$temp[0]} =$temp[4];
			} elsif ($sample_output{$temp[0]} ne $temp[4]) {
				print "Warning: Sample $temp[0] has more than 1 output places, the first appearance in the input table will be used.\n";
			}
			
			if (! exists $in_paras{$temp[0]} ) {
				$in_paras{$temp[0]} = "$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
			} else {
				$in_paras{$temp[0]} = $in_paras{$temp[0]}.'#'."$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
			}
		}
	}
	close INPUT;
} else {
	print "Error: can not find file $sample_table\n";
	exit;
}

## 2. check if Ref exists
if (-e $ref) {
	print "Found reference file: $ref\n";
} else {
	print "Error: can not locate the reference file $ref\n";
	exit;
}
#check if mito pipe dir exits, check if all required scripts are in place, check if ref seq are in place
my $QC_script = $MitoPi_dir."/modules/QC_redup.sh"; if (! -e $QC_script) {print "Error: Can not locate script: $QC_script\n"; exit;}
my $align_script = $MitoPi_dir."/modules/map_call_filter_annotate.sh"; if (! -e $align_script) {print "Error: Can not locate script: $align_script\n"; exit;}
my $bed = $MitoPi_dir."/modules/mito.bed"; if (! -e $bed) {print "Error: Can not locate bed file: $bed\n"; exit;}
my $detectLanes_script = $MitoPi_dir."/modules/detectSampleLanes.pl"; if (! -e $detectLanes_script) {print "Error: Can not locate script: $detectLanes_script\n"; exit;}

#check if annovar dir exists, check if require scripts are in place
# not in use currently
#submit jobs
# stage 2 jobs need to be hold
print "\nSubmitting jobs....\n\n";
foreach my $sample (keys %in_paras) {
	my $SAMPLE_ID=$sample;
	my $JOB_ID="Mitopi_QC_".$Library_ID_prefix."_".$sample;
	my @temp=split('#', $in_paras{$sample});
	my $OUTPUT_DIR=$sample_output{$sample};
	foreach my $in_paras_set (@temp) {
		my @temp_sets = split("\t",$in_paras_set);
		my $FASTQ1 = $temp_sets[0];
		my $FASTQ2 = $temp_sets[1];
		my $LANE = $temp_sets[2];
		#print "qsub -N $JOB_ID $QC_script $SAMPLE_ID $OUTPUT_DIR $MitoPi_dir $FASTQ1 $FASTQ2 $LANE\n";
		`qsub -N $JOB_ID $QC_script $SAMPLE_ID $OUTPUT_DIR $MitoPi_dir $FASTQ1 $FASTQ2 $LANE`;
	}
	my $JOB_ID_2="Mitopi_Mapping_".$Library_ID_prefix."_".$sample;
	#print "qsub -N $JOB_ID_2 -hold_jid $JOB_ID $align_script $SAMPLE_ID $OUTPUT_DIR $MitoPi_dir $SEQ_PLATFORM $Library_ID_prefix\n";
	`qsub -N $JOB_ID_2 -hold_jid $JOB_ID $align_script $SAMPLE_ID $OUTPUT_DIR $MitoPi_dir $SEQ_PLATFORM $Library_ID_prefix`;
}
print "Done!\nNot sure if it is succuessful though ;)\n";

exit;

sub usage {
    print "Unknown option: @_\n" if ( @_ );
    print "You have not supplied all parameters I need, please check the usage below.\n";
    print "\nusage: ********* not written yet :) ********* \n";
	return(1);
	exit;
}

