#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA_SCRIPT` }

use strict;
use warnings;
use Set::Scalar;
use Getopt::Long qw(GetOptionsFromArray);
#use Time::HiRes qw(gettimeofday tv_interval);
use timer;
use seq;

use constant { FGS => 0, MGA => 1, ALL => 2 };

my @orf_files;     ## input ORFs reads file(s)
my $spa_dir;       ## output directory
my $ncpus;         ## no. of CPUs
my $nparts;        ## no. of suffix arrays
my $kmer;          ## k-mer size
my $pairend;       ## pair-end reads flag
# my $ljoin_pe;
# my $sjoin_pe;
# my $rjoin_pe;
my $profile;       ## profile generation flag
my $alignment;     ## alignment generation flag
my $verbose;       ## verbose
my $seed_coverage;
my $seed_reuse;
my $overlap_length;
my $read_support;
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
	print STDERR "\nGenerating SPA inputs ...\n";
	my $t0 = time();
	my $cmd = "prespa -s $nparts -k $kmer -o $spa_dir ";
	$cmd .= "-i " . join(" ", @orf_files) . " &> $spa_dir/pre.log";
	print STDERR "$cmd\n";
	system( "$cmd" ) == 0 or die "Proprocessing failed: $?";
	timer::printElapsed( $t0, "SPA inputs were generated" );

	## SPA
	print STDERR "\nRunning SPA assembly ...\n";
	$t0 = time();
	$cmd  = "spa --kmer-size $kmer ";
	$cmd .= "--nparts $nparts ";
	$cmd .= "--ncpus $ncpus ";
	$cmd .= "--output-dir $spa_dir ";
	$cmd .= "--seed-coverage $seed_coverage "      if defined $seed_coverage;
	$cmd .= "--overlap-length $overlap_length "         if defined $overlap_length;
	$cmd .= "--read-support $read_support "          if defined $read_support;
	$cmd .= "--seed-reuse $seed_reuse "          if defined $seed_reuse;
	#$cmd .= "--read-bridge-extend 1 " if $read_join;
	# $cmd .= "--ljoin-pe-support 1 "   if $ljoin_pe;
	# $cmd .= "--sjoin-pe-support 1 "   if $sjoin_pe;
	# $cmd .= "--rjoin-pe-support 1 "   if $rjoin_pe;
	$cmd .= "--profile "              if $profile;
	$cmd .= "--alignment "            if $alignment;
	$cmd .= "--verbose $verbose"      if $verbose;
	$cmd .= "--pair-end "             if $pairend;
	$cmd .= "--input " . join(" ", @orf_files) . " &> $spa_dir/spa.log";
	print STDERR "$cmd\n";
	system( "$cmd" ) == 0 or die "SPA failed: $?";
	timer::printElapsed( $t0, "Assembly completed");

	my $nseq = `grep -c '>' $spa_dir/spa.fasta`; chomp $nseq;
	print STDERR "# Assembled sequences:$nseq", "\n";
	die "No sequence was assembled\n" if $nseq == 0;
	
}

sub getopts
{
	Getopt::Long::Configure("no_ignore_case");

	## print command
	print STDERR join(" ", $0, @ARGV), "\n";

	my $status = GetOptionsFromArray( \@ARGV, 
									  "input|i=s{,}"      => \@orf_files, ## at least one file
									  "pairend|M"         => \$pairend,
									  "nparts|s=i"        => \$nparts,
									  "ncpus|n=i"         => \$ncpus,
									  "kmer|k=i"          => \$kmer,
									  "seed-coverage|S=i" => \$seed_coverage,
									  "overlap-length|b=i"   => \$overlap_length,
									  "read-support|r=i"     => \$read_support,
									  "seed-reuse|u"        => \$seed_reuse,
									  #"read-join"         => \$read_join,
									  # "ljoin-pe"          => \$ljoin_pe, 
									  # "sjoin-pe"          => \$sjoin_pe, 
									  # "rjoin-pe"          => \$rjoin_pe, 
									  "profile|P"         => \$profile,
									  "alignment|A"       => \$alignment,
									  "spadir|o=s"        => \$spa_dir,
									  "verbose|V=i"       => \$verbose,
									  "help|h|?"          => \$help
									);
	
	usage() if $status != 1 || @orf_files == 0;
	die "# of suffix array must be given" unless defined $nparts;
	$kmer = 6   unless defined $kmer;
	$spa_dir = `echo -n \`pwd\`` unless defined $spa_dir;
	$ncpus = 1 unless defined $ncpus;
}

sub usage
{
	print "Usage:$0 [options]\n";
	print "      -i, --input         : [required] string " . "\tFasta-file(s) of ORFs\n";
	print "      -s, --naprts        : [required] integer" . "\tNo. of suffix arrays\n";
	print "      -k, --kmer          : [optional] integer" . "\tK-mer size (node size, default:6)\n"; 
	print "      -n, --ncpus         : [optional] integer" . "\tNumber of CPUs (default:1)\n";
	print "      -M, --pairend       : [optional] flag   " . "\t'Interleaved' pair-end reads\n"; 
	print "      -o, --spadir        : [optional] string " . "\tOutput directory (default:.)\n";
	print "      -S, --seed-coverage : [optional] integer" . "\tMinimum seed coverage (default:2)\n";
	print "      -U, --seed-reuse    : [optional] flag" . "\tResue seed kmers found in previous paths? (default:false)\n";
	print "      -b, --overlap-length: [optional] integer" . "\tMininum overlap length during suffix array search (default:15)\n";
	print "      -r, --read-support  : [optional] integer" . "\tMinimum read support (default:5)\n";
	print "      -V, --verbose       : [optional] integer" . "\tVerbosity option (0:quiet, 1:brief, 2: detailed, 3: very wordy)\n";
	print "      -P, --profile       : [optional] flag   " . "\tGenerate profile\n";
	print "      -A, --alignment     : [optional] flag   " . "\tGenerate alginment\n";
	print "      -h, --help          : [optional] flag   " . "\tGenerate this message\n"; 

	exit;
}
