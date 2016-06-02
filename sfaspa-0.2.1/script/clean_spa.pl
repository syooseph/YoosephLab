#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA_SCRIPT` }

use strict;
use warnings;
use Set::Scalar;
use Getopt::Long qw(:config no_ignore_case);
#use Time::HiRes qw(gettimeofday tv_interval);
use timer;
use seq;

use constant { FGS => 0, MGA => 1, ALL => 2 };

my $ncpus;
my $post_method;  ## post gene calling: FragGeneScan, MetaGeneAnnotator
my $min_length;   ## length cutoff of final path
my $post_dir;     ## output directory
my $orf_faa;      ## Peptide reads from gene predictor
my $orf_ffn;      ## Matching DNA reads of ORFs from gene predictor
my $spa_faa;      ## ORF file from SPA
my $spa_place;    ## SPA read placement file
my $mga_genome;   ## MGA genome option
my $numline;
my $lsuffix;
my $help;


## Get program options
getopts();

## Post-processing
cleanSpa();


sub cleanSpa
{
	my $t0 = time();
	my $ts = time();

	if ( ! -d $post_dir ) {
		system( "mkdir -p $post_dir 2>/dev/null" );
	}

	print STDERR "\nGenerating nucleotide sequences ...\n";
	my $rtran_faa = "$post_dir/spa.nuc.fasta";
	my $rtran_log = "$post_dir/rtran.log";
	my $cmd = "rtran -a $orf_faa -d $orf_ffn -p $spa_place -o $post_dir -n $ncpus &> $rtran_log";
	print STDERR "$cmd\n";
	system( "$cmd" ) == 0 or die "Reverse translation failed: $?";
	timer::printElapsed( $ts, "Nucleotide sequences generated");

	if ( $post_method == FGS || $post_method == ALL ) {
		runFGS($rtran_faa);
	} 
	if ( $post_method == MGA || $post_method == ALL ) {
		runMGA($rtran_faa);
	}
	timer::printElapsed( $ts, "Gene recalled");


	print STDERR "\nSearching sequence overcalls ...\n";
	$ts = time();
	if ( $post_method == FGS || $post_method == ALL ) {
		my $cmd = "find_overcall.pl -p $spa_faa -g $post_dir/rtran.fgs.faa > $post_dir/rtran.fgs.post 2> $post_dir/rtran.fgs.post.log";
		print STDERR "$cmd\n";
		system( "$cmd" ) == 0 or die "Overcall search failed: $?";
	}
	if ( $post_method == MGA || $post_method == ALL ) {
		my $cmd = "find_overcall.pl -p $spa_faa -g $post_dir/rtran.mga.faa > $post_dir/rtran.mga.post 2> $post_dir/rtran.mga.post.log";
		print STDERR "$cmd\n";
		system( "$cmd" ) == 0 or die "Overcall search failed: $?";
	}
	timer::printElapsed( $ts, "Overcalled sequences identified");

	my $good;
	if ( $post_method == FGS ) {
		my @Ids = `awk '{if(\$2==1) {print \$1}}' $post_dir/rtran.fgs.post`; chomp @Ids;
		$good = Set::Scalar->new( @Ids );
	} elsif ( $post_method == MGA ) {
		my @Ids = `awk '{if(\$2==1) {print \$1}}' $post_dir/rtran.mga.post`; chomp @Ids;
		$good = Set::Scalar->new( @Ids );
	} else {
		my @Ids = `awk '{if(\$2==1) {print \$1}}' $post_dir/rtran.fgs.post`; chomp @Ids;
		$good = Set::Scalar->new( @Ids );
		@Ids = `awk '{if(\$2==1) {print \$1}}' $post_dir/rtran.mga.post`; chomp @Ids;
		$good = $good * Set::Scalar->new( @Ids ); ## set intersection
	}
	
	my @Ids = `seqlen.pl -i $spa_faa -c | awk '{if(\$NF >= $min_length){ print \$1}}'`; chomp @Ids;
	$good = $good * Set::Scalar->new( @Ids ); ## set intersection

	## Extract good sequences only
	seq::extractFastaBySeqId( $good, "$spa_faa", "$post_dir/post.fasta", 0 );

	my $nseq = `grep -c '>' $post_dir/post.fasta`; chomp $nseq;
	print STDERR "# Sequences:$nseq", "\n";

	#my $te = [gettimeofday];
	#print STDERR "Postprocessing time:", tv_interval( $ts, $te ), " sec\n";
	timer::printElapsed( $t0, "Post-processing done");
}

sub runFGS
{
	my ($rtran_faa) = @_;

	if ( $ncpus > 1 ) {
		splitFiles($rtran_faa);

		print STDERR "\nRunning FragGeneScan ...\n";
		my $cmd = "ls $post_dir/x* | xargs -P $ncpus -i bash -c 'run_FragGeneScan.pl -genome=\{\} -out=\{\}.fgs -complete=1 -train=complete &> \{\}.log'";
		print STDERR "$cmd\n";

		system( "$cmd") == 0 or  die "FragGeneScan failed: $?";
		system( "cat $post_dir/x*.out > $post_dir/rtran.fgs.out" );
		system( "cat $post_dir/x*.ffn > $post_dir/rtran.fgs.ffn" );
		system( "cat $post_dir/x*.faa > $post_dir/rtran.fgs.faa" );
		#system( "cat $post_dir/x*.log > $post_dir/rtran.fgs.log" );
		system( "rm -f $post_dir/x*" );
	} else {
		print STDERR "\nRunning FragGeneScan ...\n";
		my $cmd = "run_FragGeneScan.pl -genome=$rtran_faa -out=$post_dir/rtran.fgs -complete=1 -train=complete > $post_dir/rtran.fgs.log 2> /dev/null";
		print STDERR "$cmd\n";
		system( "$cmd" ) == 0 or  die "FragGeneScan failed: $?";
	}
}

sub runMGA
{
	my ($rtran_faa) = @_;
	if ( $ncpus > 1 ) {
		splitFiles($rtran_faa);

		print STDERR "Running MetaGeneAnnotator ...\n";
		my $cmd = "ls $post_dir/x* | xargs -P $ncpus -i bash -c 'mga_linux -$mga_genome \{\} > \{\}.out'";
		print STDERR "$cmd\n";
		system( "$cmd" ) == 0 or die "MetaGeneAnnotator failed: $?";

		print STDERR "Parsing MetaGeneAnnotator outputs ...\n";
		$cmd = "ls $post_dir/x* | grep -v out | xargs -P $ncpus -i bash -c 'parse_mga.pl -m \{\}.out -s \{\} -o \{\}.faa -d \{\}.ffn'";
		print STDERR "$cmd\n";
		system( "$cmd" ) == 0 or die "Parsing MGA failed: $?";
		system( "cat $post_dir/x*.out > $post_dir/rtran.mga.out" );
		system( "cat $post_dir/x*.faa > $post_dir/rtran.mga.faa" );
		system( "cat $post_dir/x*.ffn > $post_dir/rtran.mga.ffn" );
		system( "rm -f $post_dir/x*" );
	} else {
		print STDERR "Running MetaGeneAnnotator ...\n";
		my $cmd = "mga_linux -$mga_genome $rtran_faa > $post_dir/rtran.mga.out";
		print STDERR "$cmd\n";
		system( "$cmd" ) == 0 or die "MetaGeneAnnotator failed: $?";

		print STDERR "Parsing MetaGeneAnnotator outputs ...\n";
		$cmd = "parse_mga.pl -m $post_dir/rtran.mga.out -s $rtran_faa -o $post_dir/rtran.mga.faa -d $post_dir/rtran.mga.ffn";
		print STDERR "$cmd\n";
		system( "$cmd" ) == 0 or die "Parsing MGA failed: $?";
	}
}

sub splitFiles
{
	my ($rtran_faa) = @_;
	print STDERR "\nSplitting sequence file ..\n";
	my $cmd = "cat $rtran_faa  | awk 'BEGIN{b=1}{if(/^>/&&b){printf\"%s\\n\",\$1;b=0}else if(/^>/){printf\"\\n%s\\n\", \$1} else {printf\"%s\",\$1}}END{printf\"\\n\"}' | split -l $numline -d -a $lsuffix - $post_dir/x";
	print STDERR "$cmd\n";
	system( "$cmd") == 0 or die "Split file failed: $?";
}

sub getopts
{
	## print command
	print STDERR join(" ", $0, @ARGV), "\n";

	GetOptions( "cpu|n=i"     => \$ncpus,
				"post|T=i"    => \$post_method,
				"length|w=i"  => \$min_length,
				"orf|a=s"     => \$orf_faa, 
				"dna|d=s"     => \$orf_ffn,
				"spa|s=s"     => \$spa_faa,
				"place|r=s"   => \$spa_place,
				"postdir|o=s" => \$post_dir,
				"mga|m=s"     => \$mga_genome,
				"help|h|?"    => \$help );
	usage() if $help || !$orf_faa || !$orf_ffn || !$spa_faa || !$spa_place;
	
	$ncpus       = 1   unless defined $ncpus;
	$post_method = FGS unless defined $post_method;
	$min_length = 60   unless defined $min_length;
	$mga_genome = "m"  unless defined $mga_genome;
	$post_dir = `echo -n \`pwd\`` unless defined $post_dir;

	$numline = 2000   unless defined $numline; # 1K FASTA
	$lsuffix = 10     unless defined $lsuffix; # 5 digit suffix

}

sub usage
{
	print "Usage:$0 [options]\n";
	print "      -a, --orf     : [required] string " . "\tORF file from ORF prediction\n";
	print "      -d, --dna     : [required] string " . "\tDNA file from ORF prediction\n";
	print "      -s, --spa     : [required] string " . "\tORF file predicted from SPA\n";
	print "      -r, --place   : [required] string " . "\tRead placement from SPA\n";
	print "      -n, --cpu     : [optional] integer " . "\t# of CPUs\n";
	print "      -T, --post    : [optional] integer" . "\tPost-processing gene caller (0:FragGeneScan, 1:MetaGeneAnnotator, 2:Both, default:0)\n";
	print "      -w, --length  : [optional] integer" . "\tMininum length of sequence length (default:60)\n";
	print "      -o, --postdir : [optional] string " . "\tOutput directory (default:.)\n";
	print "      -m, --mga     : [optinoal] string " . "\tMGA genome option (s/m, default:m)\n";
	print "      -L, --numline : [optional] integer" . "\tnumber of line of splitted files (default:2000 - 1K fasta)\n";
	print "      -S, --suffix  : [optional] integer" . "\tlength of suffix length of splitted files (default:10)\n";
	print "      -h, --help    : [optional] flag   " . "\tGenerate this message\n"; 
	exit;
}
