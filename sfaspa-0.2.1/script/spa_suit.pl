#!/usr/bin/env perl

# Author: Youngik Yang
# Last update: Fri Jan 30 17:35:00 EST 2015
# - Switch module is buggy
# - Replaced with perl switch which requires 5.10.1 or recent

BEGIN { push @INC, `echo -n \$SPA_SCRIPT` }

use v5.14;
use strict;
use warnings;
use Set::Scalar;
use Getopt::Long qw(GetOptionsFromArray);
use Time::HiRes qw(gettimeofday tv_interval);
#use Switch; ## buggy 
use feature "switch";
use File::Spec;
use seq;
use timer;

use constant { FGS => 0, MGA => 1, ALL => 2 };

my @fasta_files;   ## input DNA reads file(s)
my $param_file;    ## file with preset options
my $orf_call;      ## ORF calling flag
my $assembly;      ## assembly flag
my $postprocess;   ## post processing flag
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
my $seed_coverage; ## minimum seed coverge
my $seed_reuse;    ## seed reuse
my $overlap_length;    ## back trace string length
my $read_support;     ## minimum reads support
my $verbose;       ## verbosity
# my $ljoin_pe;      ## pair end reads support for long overlapping path extension
# my $sjoin_pe;      ## pair end reads support for short overlapping path extension
# my $rjoin_pe;      ## pair end reads support for read bridging path extension
my $numline;       ## no. of lines in each splitted file
my $lsuffix;       ## length of suffix for splitting file
my $nparts;        ## no. of suffix array
my $help;


my $ts = time();
timer::printLocalTime();

## Get program options
getopts();

## Make output directories
my ($orf_dir, $spa_dir, $post_dir);
makedirs();

## Set input file names
my ($faa, $ffn);
setFiles();

## Call ORFs from short DNA reads
callOrfs();

## determine no. of suffix array
getParts();

## Run SPA with ORFs
runSpa();

## Post-processing
cleanSpa();

my $te = [gettimeofday];
timer::printElapsed( $ts, "Done", 1);
timer::printLocalTime();

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

sub setFiles
{
	if ( $gene_method == FGS ) {
		$faa = $drop_indel == 0 ? "$orf_dir/fgs.post.faa" : "$orf_dir/fgs.noindel.faa";
		$ffn = $drop_indel == 0 ? "$orf_dir/fgs.post.ffn" : "$orf_dir/fgs.noindel.ffn";
	} else {
		$faa = "$orf_dir/mga.faa";
		$ffn = "$orf_dir/mga.ffn";
	}
}

##====================================================================
## Call ORFs on DNA reads using FGS or MGA
##====================================================================
sub callOrfs
{
	return unless $orf_call;

	print STDERR "\nI. Gene calling\n";

	my $cmd = "call_orfs.pl ";
	$cmd .= "--ncpu $ncpu ";
	$cmd .= "--orfdir $orf_dir ";
	$cmd .= "--line $numline ";
	$cmd .= "--suffix $lsuffix ";
	$cmd .= "--gene $gene_method ";
	$cmd .= "--complete $fgs_complete " if $gene_method == FGS;
	$cmd .= "--train $fgs_train "    if $gene_method == FGS;
	$cmd .= "--dropindel $drop_indel "   if $gene_method == FGS;
	$cmd .= "--mgagenome $mga_genome "   if $gene_method == MGA;
	$cmd .= "--input " . join(" ", @fasta_files );

	system( "$cmd" ) == 0 or die "ORF calling failed: $?";

	timer::printElapsed( $ts, "Peptide sequence generated", 1);
	print STDERR '-' x 80, "\n";
}

sub getParts
{
	print STDERR "\nDetermining number of suffix array partitions ...\n";

	my $cmd = "part -i $faa";
	print STDERR $cmd, "\n";
	system( "$cmd > $spa_dir/npart" ) == 0 or die "Determining # suffix array failed: $?";	
	$nparts = `grep parts $spa_dir/npart | cut -f 2 -d ':'`; chomp $nparts;
	if ( $nparts eq "" || $nparts == 0 ) {
		die "Invalid number of suffix arrays";
	}
	print STDERR "# suffix arrays:", $nparts, "\n";

	timer::printElapsed( $ts, "No. of suffix array partitions calculated", 1);
	print STDERR '-' x 80, "\n";
}

sub runSpa
{
	return unless $assembly ;

	print STDERR "\nII. Assembly\n";

	my $cmd = "run_spa.pl ";
	$cmd .= "--kmer $kmer ";
	$cmd .= "--nparts $nparts ";
	$cmd .= "--ncpus $ncpu ";
	$cmd .= "--spadir $spa_dir ";
	$cmd .= "--seed-coverage $seed_coverage "    if defined $seed_coverage;
	$cmd .= "--overlap-length $overlap_length "  if defined $overlap_length;
	$cmd .= "--read-support $read_support "      if defined $read_support;
	$cmd .= "--seed-reuse "                      if defined $seed_reuse and $seed_reuse;
	# $cmd .= "--ljoin-pe "         if $ljoin_pe;
	# $cmd .= "--sjoin-pe "         if $sjoin_pe;
	# $cmd .= "--rjoin-pe "         if $rjoin_pe;

	$cmd .= "--profile "          if $profile;
	$cmd .= "--alignment "        if $alignment;
	$cmd .= "--verbose $verbose " if $verbose;
	
	$cmd .= "--pairend " if $pairend;
	$cmd .= "--input $faa";

	system( "$cmd" ) == 0 or die "Running SPA failed: $?";
	timer::printElapsed( $ts, "SPA assembly completed");
	print STDERR '-' x 80, "\n";
}

sub cleanSpa
{
	return unless $postprocess ;

	print STDERR "\nIII. Post-processing\n";

	my $cmd = "clean_spa.pl ";
	$cmd .= "--postdir $post_dir ";
	$cmd .= "--orf $faa ";
	$cmd .= "--dna $ffn ";
	$cmd .= "--cpu $ncpu ";
	$cmd .= "--spa $spa_dir/spa.fasta ";
	$cmd .= "--place $spa_dir/spa.place.bin ";
	$cmd .= "--post $post_method ";
	$cmd .= "--mga $mga_genome " if ( $post_method == MGA || $post_method == ALL );
	$cmd .= "--length $min_length ";
	
	system( "$cmd" ) == 0 or die "Post-processing failed: $?";
	timer::printElapsed( $ts, "SPA post-processing completed", 1);
	print STDERR '-' x 80, "\n";
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
		for ( $opt ) {
			when ("call")           { $orf_call     = $val; }
			when ("assemble")       { $assembly     = $val; }
			when ("clean")          { $postprocess  = $val; }
			when ("gene")           { $gene_method  = $val; }
			when ("complete")       { $fgs_complete = $val; }
			when ("train")          { $fgs_train    = $val; }
			when ("dropindel")      { $drop_indel   = $val; }
			when ("mga")            { $mga_genome   = $val; }
			when ("ncpu")           { $ncpu         = $val; }
			when ("kmer")           { $kmer         = $val; }
			when ("paired")         { $pairend      = $val; }
			when ("seed-coverage")  { $seed_coverage = $val; }
			when ("seed-reuse")     { $seed_reuse   = $val; }
			when ("read-support")   { $read_support = $val; }
			when ("overlap-length") { $overlap_length = $val; }
			when ("profile")        { $profile      = $val; }							  
			when ("alignment")      { $alignment    = $val; }							  
			when ("post")           { $post_method  = $val; }							  
			when ("length")         { $min_length   = $val; }							  
			when ("odir")           { $odir         = $val; }							  
			when ("verbose")        { $verbose      = $val; }							
			# when ("ljoin-pe"      { $ljoin_pe     = $val; }
			# when ("sjoin-pe"      { $sjoin_pe     = $val; }
			# when ("rjoin-pe"      { $rjoin_pe     = $val; }
		}
	}
	close(IN);		
	$odir = File::Spec->rel2abs($odir); 
}

sub checkParameters
{
	my $max_cpu = `grep -w "processor" /proc/cpuinfo | wc -l`; chomp $max_cpu;
	if ( $ncpu =~ /max/i ) {$ncpu = $max_cpu; } 
	elsif ( $ncpu =~ /min/i ) { $ncpu = 1; }
	elsif ( $ncpu =~ /([0-9]+)/ ) { $ncpu = $1; $ncpu = $max_cpu if $ncpu > $max_cpu; }
	else { die "[Error] Invalid value for cpu option\n"}
	if ( ! defined $kmer ) { die "[ERROR] Size of k-mer must be given\n"; }
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
	if ( $seed_coverage < 1 ) { die "[ERROR] Invalid value for seed-coverage option\n"; }
	if ( $seed_reuse != 0 && $seed_reuse != 1 ) { die "[ERROR] Invalid value for seed-reuse option\n"; }
	if ( $read_support < 1 ) { die "[ERROR] Invalid value for mininum reads support option\n"; }
	if ( $overlap_length < 1 ) { die "[ERROR] Invalid valude for minimum overlap length option\n"; }
	if ( $profile != 0 && $profile != 1 ) { die "[ERROR] Invalid value for profile option\n"; }
	if ( $alignment != 0 && $alignment != 1 ) { die "[ERROR] Invalid value for alignment option\n"; }
	if ( $post_method != 0 && $post_method != 1 && $post_method != 2) { die "[ERROR] Invalid value for post processing method\n"; }
	if ( $min_length < 1) { die "[ERROR] Invalid value for minimum length option\n"; }
	if ( $verbose < 0 || $verbose > 3 ) { die "[ERROR] Invalid value for verbose option\n"; }
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
