#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA` }

use strict;
use warnings;
use Set::Scalar;
use Getopt::Long qw(:config no_ignore_case);
use Time::HiRes qw(gettimeofday tv_interval);

use seq;

use constant { FGS => 0, MGA => 1, ALL => 2 };

my $post_method;  ## post gene calling: FragGeneScan, MetaGeneAnnotator
my $min_length;   ## length cutoff of final path
my $post_dir;     ## output directory
my $orf_faa;      ## Peptide reads from gene predictor
my $orf_ffn;      ## Matching DNA reads of ORFs from gene predictor
my $spa_faa;      ## ORF file from SPA
my $spa_place;    ## SPA read placement file
my $mga_genome;   ## MGA genome option
my $help;

## Get program options
getopts();

## Post-processing
cleanSpa();

sub cleanSpa
{
	my $ts = [gettimeofday];

	if ( ! -d $post_dir ) {
		system( "mkdir -p $post_dir 2>/dev/null" );
	}

	my $rtran_faa = "$post_dir/spa.rtran.fasta";
	my $rtran_log = "$post_dir/spa.rtran.log";
	my $cmd = "rtran -a $orf_faa -d $orf_ffn -p $spa_place -o $rtran_faa >& $rtran_log";
	print STDERR $cmd, "\n";
	system( "$cmd" ) == 0 or die "Reverse translation failed: $?";
	if ( $post_method == FGS || $post_method == ALL ) {
		$cmd = "run_FragGeneScan.pl -genome=$rtran_faa -out=$post_dir/spa.rtran.fgs -complete=1 -train=complete > $post_dir/spa.rtran.fgs.log 2> /dev/null";
		print STDERR $cmd, "\n";
		system( "$cmd" ) == 0 or  die "FragGeneScan failed: $?";
		$cmd = "find_overcall.pl -p $spa_faa -g $post_dir/spa.rtran.fgs.faa > $post_dir/spa.rtran.fgs.post 2> $post_dir/spa.rtran.fgs.post.log";
		print STDERR $cmd, "\n";
		system( "$cmd" ) == 0 or die "Overcall search failed: $?";
	} 
	if ( $post_method == MGA || $post_method == ALL ) {
		$cmd = "mga_linux_ia64 -$mga_genome $rtran_faa > $post_dir/spa.rtran.mga.out";
		print STDERR $cmd, "\n";
		system( "$cmd" ) == 0 or die "MetaGeneAnnotator failed: $?";
		$cmd = "parse_mga.pl -m $post_dir/spa.rtran.mga.out -s $rtran_faa -o $post_dir/spa.rtran.mga.faa -d $post_dir/spa.rtran.mga.ffn";
		print STDERR $cmd, "\n";
		system( "$cmd" ) == 0 or die "Parsing MGA failed: $?";
		$cmd = "find_overcall.pl -p $spa_faa -g $post_dir/spa.rtran.mga.faa > $post_dir/spa.rtran.mga.post 2> $post_dir/spa.rtran.mga.post.log";
		print STDERR $cmd, "\n";
		system( "$cmd" ) == 0 or die "Overcall search failed: $?";
	}

	my $good;
	if ( $post_method == FGS ) {
		my @Ids = `awk '{if(\$2==1) {print \$1}}' $post_dir/spa.rtran.fgs.post`; chomp @Ids;
		$good = Set::Scalar->new( @Ids );
	} elsif ( $post_method == MGA ) {
		my @Ids = `awk '{if(\$2==1) {print \$1}}' $post_dir/spa.rtran.mga.post`; chomp @Ids;
		$good = Set::Scalar->new( @Ids );
	} else {
		my @Ids = `awk '{if(\$2==1) {print \$1}}' $post_dir/spa.rtran.fgs.post`; chomp @Ids;
		$good = Set::Scalar->new( @Ids );
		@Ids = `awk '{if(\$2==1) {print \$1}}' $post_dir/spa.rtran.mga.post`; chomp @Ids;
		$good = $good * Set::Scalar->new( @Ids ); ## set intersection
	}
	
	my @Ids = `seqlen.pl -i $spa_faa -c | awk '{if(\$2 >= $min_length){ print \$1}}'`; chomp @Ids;
	$good = $good * Set::Scalar->new( @Ids ); ## set intersection

	## Extract good sequences only
	seq::extractFastaBySeqId( $good, "$spa_faa", "$post_dir/spa.post.fasta", 0 );

	my $nseq = `grep -c '>' $post_dir/spa.post.fasta`; chomp $nseq;
	print STDERR "# Sequences:$nseq", "\n";

	my $te = [gettimeofday];
	print STDERR "Postprocessing time:", tv_interval( $ts, $te ), " sec\n";
}

sub getopts
{
	## print command
	print STDERR join(" ", $0, @ARGV), "\n";

	GetOptions( "post|T=i"    => \$post_method,
				"length|w=i"  => \$min_length,
				"orf|a=s"     => \$orf_faa, 
				"dna|d=s"     => \$orf_ffn,
				"spa|s=s"     => \$spa_faa,
				"place|r=s"   => \$spa_place,
				"postdir|o=s" => \$post_dir,
				"mga|m=s"     => \$mga_genome,
				"help|h|?"    => \$help );
	usage() if $help || !$orf_faa || !$orf_ffn || !$spa_faa || !$spa_place;
	
	$post_method = ALL unless defined $post_method;
	$min_length = 60   unless defined $min_length;
	$mga_genome = "m"  unless defined $mga_genome;
	$post_dir = `echo -n \`pwd\`` unless defined $post_dir;
}

sub usage
{
	print "Usage:$0 [options]\n";
	print "      -a, --orf     : [required] string " . "\tORF file from ORF prediction\n";
	print "      -d, --dna     : [required] string " . "\tDNA file from ORF prediction\n";
	print "      -s, --spa     : [required] string " . "\tORF file predicted from SPA\n";
	print "      -r, --place   : [required] string " . "\tRead placement from SPA\n";
	print "      -T, --post    : [optional] integer" . "\tPost-processing gene caller (0:FragGeneScan, 1:MetaGeneAnnotator, 2:Both, default:2)\n";
	print "      -w, --length  : [optional] integer" . "\tMininum length of sequence length (default:60)\n";
	print "      -o, --postdir : [optional] string " . "\tOutput directory (default:.)\n";
	print "      -m, --mga     : [optinoal] string " . "\tMGA genome option (s/m, default:m)\n";
	print "      -h, --help    : [optional] flag   " . "\tGenerate this message\n"; 
	exit;
}
