#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA` }

use strict;
use warnings;
use Set::Scalar;
use Getopt::Long qw(GetOptionsFromArray);
use Time::HiRes qw(gettimeofday tv_interval);
use Switch;
use File::Spec;
use seq;

use constant { FGS => 0, MGA => 1, ALL => 2 };

my @fasta_files;   ## input DNA reads file(s)
my $param_file;    ## file with preset options
my $gene_method;   ## gene calling methods: FragGeneScan, MetaGeneAnnotator
my $post_method;   ## post gene calling: FragGeneScan, MetaGeneAnnotator
my $min_length;    ## length cutoff of final path
my $ncpu;          ## number of CPUs
my $odir;          ## output directory
my $fgs_complete;  ## FGS complete option (0/1)
my $fgs_train;     ## FGS training option
my $drop_indel;    ## FGS drop indeled prediction
my $mga_genome;    ## MGA genome option
my $kmer;          ## k-mer size
my $pairend;       ## pair-end reads flag
my $profile;       ## profile generation flag
my $alignment;     ## alignment generation flag
my $min_coverage;  ## minimum coverage of trimming graph
my $no_trim;
my $percentile;    
my $distance;
my $min_share;
my $med_coverage;
my $base_depth;
my $verbose;       ## verbose
my $numline;       ## no. of lines in each splitted file
my $lsuffix;       ## length of suffix for splitting file
my $help;

my $ts = [gettimeofday];

## Get program options
getopts();

## Make output directories
my ($orf_dir, $spa_dir, $post_dir);
makedirs();

## Call ORFs from short DNA reads
my ($faa, $ffn);
callOrfs();

## Run SPA with ORFs
runSpa();

## Post-processing
cleanSpa();

my $te = [gettimeofday];
print STDERR "Total elapsed time:", tv_interval( $ts, $te ), " sec\n";

##====================================================================
## Make sub directories 
##====================================================================
sub makedirs
{
	$orf_dir  = "$odir/orf";
	$spa_dir  = "$odir/spa";
	$post_dir = "$odir/post";

	my $cmd = "mkdir -p $orf_dir $spa_dir $post_dir 2> /dev/null";
	system( "$cmd" );
}

##====================================================================
## Call ORFs on DNA reads using FGS or MGA
##====================================================================
sub callOrfs
{

	my $cmd = "call_orfs.pl ";
	$cmd .= "-n $ncpu ";
	$cmd .= "-o $orf_dir ";
	$cmd .= "-L $numline ";
	$cmd .= "-S $lsuffix ";
	$cmd .= "-g $gene_method ";
	$cmd .= "-f $fgs_complete " if $gene_method == FGS;
	$cmd .= "-t $fgs_train "    if $gene_method == FGS;
	$cmd .= "-d $drop_indel "   if $gene_method == FGS;
	$cmd .= "-m $mga_genome "   if $gene_method == MGA;
	$cmd .= "-i " . join(" ", @fasta_files );

	print STDERR $cmd, "\n";
	system( "$cmd" ) == 0 or die "ORF calling failed: $?";

	if ( $gene_method == FGS ) {
		$faa = $drop_indel == 0 ? "$orf_dir/fgs.post.faa" : "$orf_dir/fgs.noindel.faa";
		$ffn = $drop_indel == 0 ? "$orf_dir/fgs.post.ffn" : "$orf_dir/fgs.noindel.ffn";
	} else {
		$faa = "$orf_dir/mga.faa";
		$ffn = "$orf_dir/mga.ffn";
	}
}


sub runSpa
{
	my $cmd = "run_spa.pl ";
	$cmd .= "-o $spa_dir ";
	$cmd .= "-c $min_coverage "  if defined $min_coverage;
	$cmd .= "-p $percentile "    if defined $percentile;
	$cmd .= "-d $distance "      if defined $distance;
	$cmd .= "-r $min_share "     if defined $min_share;
	$cmd .= "-C $med_coverage "  if defined $med_coverage;
	$cmd .= "-b $base_depth "    if defined $base_depth;
	$cmd .= "--profile "         if $profile;
	$cmd .= "--alignment "       if $alignment;
	$cmd .= "--no-trim "         if $no_trim;
	$cmd .= "-V "  if $verbose;
	$cmd .= "-M "  if $pairend;
	$cmd .= "-i $faa";

	print STDERR $cmd, "\n";
	system( "$cmd" ) == 0 or die "Running SPA failed: $?";
}

sub cleanSpa
{
	my $cmd = "clean_spa.pl ";
	$cmd .= "-o $post_dir ";
	$cmd .= "-a $faa ";
	$cmd .= "-d $ffn ";
	$cmd .= "-s $spa_dir/spa.fasta ";
	$cmd .= "-r $spa_dir/spa.place ";
	$cmd .= "-T $post_method ";
	$cmd .= "-m $mga_genome " if ( $post_method == MGA || $post_method == ALL );
	$cmd .= "-w $min_length ";
	
	print STDERR $cmd, "\n";
	system( "$cmd" ) == 0 or die "Post-processing failed: $?";
}



sub getopts
{
	Getopt::Long::Configure("no_ignore_case");
	## print command
	print STDERR join(" ", $0, @ARGV), "\n";

	my $status = GetOptionsFromArray( \@ARGV, 
									  "parameters|p=s"   => \$param_file,
									  "output|o=s"       => \$odir,
									  "input|i=s{,}"     => \@fasta_files, ## at least one file
									  "help|h"           => \$help
									);
	
	if ( $status != 1 ) {
	   usage();
	}
	if ( $help ) {
		usage();
	}

	if ( ! defined $param_file ) {
		print "[ERROR] parameter file must be given\n";
		usage();
	}
	if ( @fasta_files == 0 ) {
		print "[ERROR] Input file(s) must be given\n";
		usage();
	}
	
	loadParameters();
	
	
	$odir = `echo -n \`pwd\`` unless defined $odir;

	## file split options
	$numline = 20000; # 10K FASTA

	## should be sufficient
	$lsuffix = 10; #  # 10 digit suffix

	checkParameters();
}

sub loadParameters
{
	open(IN, "< $param_file") or die $!;
	while ( <IN>) {
		chomp;
		next if /^#/;
		my ( $opt, $val ) = split;
		switch ( $opt ) {
			case "gene"         { $gene_method  = $val; next }
			case "complete"     { $fgs_complete = $val; next }
		    case "train"        { $fgs_train    = $val; next }
			case "dropindel"    { $drop_indel   = $val; next }
			case "mga"          { $mga_genome   = $val; next }
			case "ncpu"         { $ncpu         = $val; next }
			case "kmer"         { $kmer         = $val; next }
			case "paired"       { $pairend      = $val; next }
			case "trim"         { $no_trim      =!$val; next }
			case "min-coverage" { $min_coverage = $val; next }
			case "min-share"    { $min_share    = $val; next }
			case "distance"     { $distance     = $val; next }
			case "base-depth"   { $base_depth   = $val; next }
			case "med-coverage" { $med_coverage = $val; next }
			case "percentile"   { $percentile   = $val; next }
			case "profile"      { $profile      = $val; next }							  
			case "alignment"    { $alignment    = $val; next }							  
			case "post"         { $post_method  = $val; next }							  
			case "length"       { $min_length   = $val; next }							  
			case "odir"         { $odir         = $val; next }							  
			case "verbose"      { $verbose      = $val; next }							  			
		}
	}
	close(IN);		
	$odir = File::Spec->rel2abs($odir); 
}

sub checkParameters
{
	if ( $gene_method != 0  && $gene_method != 1 ) { die "[ERROR] Invalid value for gene finder:\n"; }
	if ( $fgs_complete != 0 && $fgs_complete != 1 ) { die "[ERROR] Invalid value for FGS complete option\n"; }
	if ( $fgs_train ne "complete"   && 
		 $fgs_train ne "sanger_5"   && 
		 $fgs_train ne "sanger_10"  &&
		 $fgs_train ne "454_10"     && 
		 $fgs_train ne "454_30"     &&
		 $fgs_train ne "illumina_5" &&
		 $fgs_train ne "illumina_10" ) { die "[ERROR] Invalid value for FGS train option\n"; }
	if ( $drop_indel != 0 && $drop_indel != 1 ) { die "[ERROR] Invalid value for FGS dropindel option\n"; }
	if ( $mga_genome ne "s" && $mga_genome ne "m" ) { die "[ERROR] Invalid value for MGA genome option\n"; }
	if ( $pairend != 0 && $pairend != 1 ) { die "[ERROR] Invalid for paired end option\n"; }
	if ( $no_trim != 0 && $no_trim != 1 ) { die "[ERROR] Invalid for graph trimming option\n"; }
	if ( $min_coverage < 1 ) { die "[ERROR] Invalid value for min-coverage option\n"; }
	if ( $min_share < 1 ) { die "[ERROR] Invalid value for min-share option\n"; }
	if ( $distance < 1 ) { die "[ERROR] Invalid valude for distance option\n"; }
	if ( $base_depth < 1 ) { die "[ERROR] Invalid value for base-depth option\n"; }
	if ( $med_coverage < 1 ) { die "[ERROR] Invalid value for med-coverage option\n"; }
	if ( $percentile <= 0 || $percentile > 100 ) { die "[ERROR] Invalid value for percentile option\n"; }
	if ( $profile != 0 && $profile != 1 ) { die "[ERROR] Invalid value for profile option\n"; }
	if ( $alignment != 0 && $alignment != 1 ) { die "[ERROR] Invalid value for alignment option\n"; }
	if ( $post_method != 0 && $post_method != 1 && $post_method != 2) { die "[ERROR] Invalid value for post processing method\n"; }
	if ( $min_length < 1) { die "[ERROR] Invalid value for minimum length option\n"; }
	if ( $verbose != 0 && $verbose != 1 ) { die "[ERROR] Invalid value for verbose option\n"; }
}

sub usage
{
	print "Usage:$0 [-o <output directory>] -p <parameter file> -i <FASTA files>\n";
	print "      -p, --parameters : [required] string " . "\tProgram parameter file\n";
	print "      -i, --input      : [required] string " . "\tFASTA file(s) of DNA reads\n";
	print "      -o, --output     : [optional] string"  . "\tOutput directory (default:.)\n";
	print "      -h, --help                           " . "\tPrint this message\n";

	exit;
}
