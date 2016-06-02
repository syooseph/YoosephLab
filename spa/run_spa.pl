#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA` }

use strict;
use warnings;
use Set::Scalar;
use Getopt::Long qw(GetOptionsFromArray);
use Time::HiRes qw(gettimeofday tv_interval);

use seq;

use constant { FGS => 0, MGA => 1, ALL => 2 };

my @orf_files;     ## input ORFs reads file(s)
my $spa_dir;       ## output directory
my $kmer;          ## k-mer size
my $pairend;       ## pair-end reads flag
my $extend_read;   ## allow to latch reads to a path
my $profile;       ## profile generation flag
my $alignment;     ## alignment generation flag
my $binary_dump;   ## dump binary output
my $dsuffix;       ## a suffix of dump objects
my $verbose;       ## verbose
my $no_trim;
my $min_coverage;  ## minimum coverage of trimming graph
my $percentile;    
my $distance;
my $min_share;
my $med_coverage;
my $base_depth;
my $help;

## Get program options
getopts();

## Run SPA with ORFs
runSpa();


sub runSpa
{
	my $ts = [gettimeofday];

	if ( ! -d $spa_dir ) {
		system( "mkdir -p $spa_dir 2> /dev/null" );
	} 

	my $graph_file = "$spa_dir/gin";
	my $index_file = "$spa_dir/iin";

	## Preprocessing
	my $t0 = [gettimeofday];
	my $cmd = "prespa -k $kmer -G $graph_file -I $index_file ";
	$cmd .= "-i " . join(" ", @orf_files) . " >& $spa_dir/pre.log";
	print STDERR $cmd, "\n";
	system( "$cmd" ) == 0 or die "Proprocessing failed: $?";
	my $t1 = [gettimeofday];
	print STDERR "SPA preprocessing time:", tv_interval( $t0, $t1 ), " sec\n";	

	## SPA
	$t0 = [gettimeofday];
	$cmd  = "spa -k $kmer -G $graph_file -I $index_file -o $spa_dir ";
	$cmd .= "-c $min_coverage "  if defined $min_coverage;
	$cmd .= "-p $percentile "    if defined $percentile;
	$cmd .= "-d $distance "      if defined $distance;
	$cmd .= "-r $min_share "     if defined $min_share;
	$cmd .= "-C $med_coverage "  if defined $med_coverage;
	$cmd .= "-b $base_depth "    if defined $base_depth;
	$cmd .= "--profile "      if $profile;
	$cmd .= "--alignment "    if $alignment;
	$cmd .= "--verbose "      if $verbose;
	$cmd .= "--no-trim "       if $no_trim;
	$cmd .= "--dump-flag --dump-suffix $dsuffix "  if $binary_dump;
	$cmd .= "--pair-end "     if $pairend;
	$cmd .= "-i " . join(" ", @orf_files) . " >& $spa_dir/spa.log";
	print STDERR $cmd, "\n";
	system( "$cmd" ) == 0 or die "SPA failed: $?";

	$t1 = [gettimeofday];
	print STDERR "SPA assembly time:", tv_interval( $t0, $t1 ), " sec\n";	

	my $te = [gettimeofday];
	print STDERR "SPA running time:", tv_interval( $ts, $te ), " sec\n";

	my $nseq = `grep -c '>' $spa_dir/spa.fasta`; chomp $nseq;
	print STDERR "# Sequences:$nseq", "\n";
	die "No sequence assembled\n" if $nseq == 0;
}

sub getopts
{
	Getopt::Long::Configure("no_ignore_case");

	## print command
	print STDERR join(" ", $0, @ARGV), "\n";

	my $status = GetOptionsFromArray( \@ARGV, 
									  "input|i=s{,}"     => \@orf_files, ## at least one file
									  "pairend|M"        => \$pairend,
									  "kmer|k=i"         => \$kmer,
									  "min-coverage|c=i" => \$min_coverage,
									  "percentile|p=f"   => \$percentile,
									  "distance|d=i"     => \$distance,
									  "min-share|r=i"    => \$min_share,
									  "med-coverage|C=i" => \$med_coverage,
									  "base-depth|b=i"   => \$base_depth,
									  "profile|P"        => \$profile,
									  "alignment|A"      => \$alignment,
									  "dump-flag|D"      => \$binary_dump,
									  "dump-suffix|X=s"  => \$dsuffix,
									  "no-trim|N"        => \$no_trim,
									  "spadir|o=s"       => \$spa_dir,
									  "verbose|V"        => \$verbose,
									  "help|h|?"         => \$help
									);
	
	usage() if $status != 1 || @orf_files == 0;
	
	$kmer = 6   unless defined $kmer;
	$spa_dir = `echo -n \`pwd\`` unless defined $spa_dir;
	$dsuffix = "spa" unless defined $dsuffix;
}

sub usage
{
	print "Usage:$0 [options]\n";
	print "      -i, --input       : [required] string " . "\tFasta-file(s) of ORFs\n";
	print "      -k, --kmer        : [optional] integer" . "\tK-mer size (node size, default:6)\n";
	print "      -M, --pairend     : [optional] flag   " . "\t'Interleaved' pair-end reads\n"; 
	print "      -o, --spadir      : [optional] string " . "\tOutput directory (default:.)\n";
	print "      -N, --no-trim     : [optional] flag"    . "\tSkip graph trimming\n";
	print "      -c, --min-coverage: [optional] integer" . "\tMinimum kmer coverage for trimming graph (default:5)\n";
	print "      -p, --percentile  : [optional] float"   . "\tTop x% percentile of seed kmers (default:50)\n";
	print "      -d, --distance    : [optional] integer" . "\tMaximum distant node in a path to find same reads in a current node (default:20)\n";
	print "      -r, --min-share   : [optional] integer" . "\tMinimum node depth in graph traversal (default:5)\n";
	print "      -C, --med-coverage: [optional] integer" . "\tMedian base depth in consensus (default:10)\n";
	print "      -b, --base-depth  : [optional] integer" . "\tMinimum base coverage in consensus (default:3)\n";
	print "      -P, --profile     : [optional] flag   " . "\tGenerate profile\n";
	print "      -A, --alignment   : [optional] flag   " . "\tGenerate alginment\n";
	print "      -D, --dump        : [optional] flag   " . "\tdump program objects in binary\n";
	print "      -X, --suffix      : [optional] string " . "\tsuffix of dump objects\n";
	print "      -h, --help        : [optional] flag   " . "\tGenerate this message\n"; 

	exit;
}
