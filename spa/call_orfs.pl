#!/usr/bin/env perl

BEGIN { push @INC, `echo -n \$SPA` }

use strict;
use warnings;
use Set::Scalar;
use Getopt::Long qw(GetOptionsFromArray);
use Time::HiRes qw(gettimeofday tv_interval);

use seq;

use constant { FGS => 0, MGA => 1, ALL => 2 };

my @fasta_files;  ## input DNA reads file(s)
my $gene_method;  ## gene calling methods: FragGeneScan, MetaGeneAnnotator
my $ncpu;         ## number of CPUs
my $orf_dir;      ## output directory
my $fgs_complete; ## FGS complete option (0/1)
my $fgs_train;    ## FGS training option
my $drop_indel;   ## FGS drop indeled prediction
my $mga_genome;   ## MGA genome option
my $numline;      ## no. of lines in each splitted file
my $lsuffix;      ## length of suffix for splitting file
my $help;

## Get program options
getopts();

## Call ORFs from short DNA reads
my ($faa, $ffn);
callOrfs();

sub callOrfs
{
	my $t0 = [gettimeofday];

	splitFiles();

	if ( ! -d $orf_dir ) {
		system( "mkdir -p $orf_dir 2> /dev/null" );
	} 
	if ( $gene_method == FGS ) {
		## FGS calling
		my $cmd = "ls $orf_dir/x* | xargs -P $ncpu -i sh -c 'run_FragGeneScan.pl -genome=\{\} -out=\{\}.fgs -complete=$fgs_complete -train=$fgs_train >& \{\}.log'";
		print STDERR $cmd, "\n";
		system( "$cmd") == 0 or  die "FragGeneScan failed: $?";
		system( "cat $orf_dir/x*.out > $orf_dir/fgs.out" );
		system( "cat $orf_dir/x*.ffn > $orf_dir/fgs.ffn" );
		system( "cat $orf_dir/x*.faa > $orf_dir/fgs.faa" );
		system( "cat $orf_dir/x*.log > $orf_dir/fgs.log" );

		## post processing
		$cmd = "ls $orf_dir/x* | awk '/x[0-9]+\$/' |  xargs -P $ncpu -i sh -c 'postfgs -t . -o \{\}.fgs.out -n \{\}.fgs.ffn -a \{\}.fgs.faa \{\} >& \{\}.fgs.post.log'";
		print STDERR $cmd, "\n";
		system( "$cmd" ) == 0 or die "Postprocessing of FGS failed: $?";
		system( "cat $orf_dir/x*.out.post > $orf_dir/fgs.post.out" );
		system( "cat $orf_dir/x*.faa.post > $orf_dir/fgs.post.faa" );
		system( "cat $orf_dir/x*.ffn.post > $orf_dir/fgs.post.ffn" );
		system( "cat $orf_dir/x*.post.log > $orf_dir/fgs.post.log" );

		if ( $drop_indel ) {
			## trim indels
			$cmd = "ls $orf_dir/x*.fgs.out | awk -F/ '{print \$NF}' | cut -f 1 -d '.' | xargs -P $ncpu -i  sh -c 'find_indel.pl -f $orf_dir/\{\}.fgs.out > $orf_dir/\{\}.indel.ids'";
			print STDERR $cmd, "\n";
			system( "$cmd" ) == 0 or die "Finding indel output failed: $?";
			$cmd = "ls $orf_dir/x*.fgs.out | awk -F/ '{print \$NF}' | cut -f 1 -d '.' | xargs -P $ncpu -i  sh -c 'trim_FGS.pl -f $orf_dir/\{\}.fgs.faa.post -i $orf_dir/\{\}.indel.ids -n $orf_dir/\{\}.noindel.faa'";
			print STDERR $cmd, "\n";
			system( "$cmd" ) == 0 or die "Triming indeled faa failed: $?";
			$cmd = "ls $orf_dir/x*.fgs.out | awk -F/ '{print \$NF}' | cut -f 1 -d '.' | xargs -P $ncpu -i  sh -c 'trim_FGS.pl -f $orf_dir/\{\}.fgs.ffn.post -i $orf_dir/\{\}.indel.ids -n $orf_dir/\{\}.noindel.ffn'";
			print STDERR $cmd, "\n";
			system( "$cmd" ) == 0 or die "Triming indeled ffn failed: $?";
			system( "cat $orf_dir/x*ids > $orf_dir/fgs.indel.ids" );
			system( "cat $orf_dir/x*noindel.faa > $orf_dir/fgs.noindel.faa" );
			system( "cat $orf_dir/x*noindel.ffn > $orf_dir/fgs.noindel.ffn" );
		} 
	} else {
		## MGA 
		my $cmd = "ls $orf_dir/x* | xargs -P $ncpu -i sh -c 'mga_linux_ia64 -$mga_genome \{\} > \{\}.out'";
		print STDERR $cmd, "\n";
		system( "$cmd" ) == 0 or die "MetaGeneAnnotator failed: $?";
		
		## Parsing MGA
		$cmd = "ls $orf_dir/x* | grep -v out | xargs -P $ncpu -i sh -c 'parse_mga.pl -m \{\}.out -s \{\} -o \{\}.faa -d \{\}.ffn'";
		print STDERR $cmd, "\n";
		system( "$cmd" ) == 0 or die "Parsing MGA failed: $?";
		system( "cat $orf_dir/x*.out > $orf_dir/mga.out" );
		system( "cat $orf_dir/x*.faa > $orf_dir/mga.faa" );
		system( "cat $orf_dir/x*.ffn > $orf_dir/mga.ffn" );
	}
	system( "rm -f $orf_dir/x*" );


	my $t1 = [gettimeofday];
	print STDERR "ORF called:", tv_interval( $t0, $t1 ), " sec\n";
}

sub splitFiles
{	
	if ( $fasta_files[0] =~ /\.gz$/ ) {
		my $cmd = "zcat @fasta_files | awk 'BEGIN{b=1}{if(/^>/&&b){printf\"%s\\n\",\$1;b=0}else if(/^>/){printf\"\\n%s\\n\", \$1} else {printf\"%s\",\$1}}END{printf\"\\n\"}' | split -l $numline -d -a $lsuffix - $orf_dir/x";
		print STDERR $cmd, "\n";
		system( "$cmd") == 0 or die "Split file failed: $?";
	} else {
		my $cmd = "cat @fasta_files  | awk 'BEGIN{b=1}{if(/^>/&&b){printf\"%s\\n\",\$1;b=0}else if(/^>/){printf\"\\n%s\\n\", \$1} else {printf\"%s\",\$1}}END{printf\"\\n\"}' | split -l $numline -d -a $lsuffix - $orf_dir/x";
		print STDERR $cmd, "\n";
		system( "$cmd") == 0 or die "Split file failed: $?";
	}
}


sub getopts
{
	Getopt::Long::Configure("no_ignore_case");

	## print command
	print STDERR join(" ", $0, @ARGV), "\n";

	my $status = GetOptionsFromArray( \@ARGV, 
									  "input|i=s{,}"  => \@fasta_files, ## at least one file
									  "gene|g=i"      => \$gene_method,
									  "ncpu|n=i"      => \$ncpu,
									  "orfdir|o=s"    => \$orf_dir,
									  "complete|f=i"  => \$fgs_complete,
									  "train|t=s"     => \$fgs_train,
									  "dropindel|d=i" => \$drop_indel,
									  "mgagenome|m=s" => \$mga_genome,
									  "line|L=i"      => \$numline,
									  "suffix|S=i"    => \$lsuffix,
									  "help|h|?"      => \$help
									);
	
	usage() if $status != 1 || @fasta_files == 0;
	
	$gene_method = FGS unless defined $gene_method;
	$drop_indel = 1  unless defined $drop_indel;

	$fgs_complete = 0          unless defined $fgs_complete;
	$fgs_train = "illumina_10" unless defined $fgs_train;
	$mga_genome = "m"          unless defined $mga_genome;

	$ncpu = 1   unless defined $ncpu;
	$orf_dir = `echo -n \`pwd\`` unless defined $orf_dir;

	## file split options
	$numline = 200000 unless defined $numline; # 100K FASTA
	$lsuffix = 5      unless defined $lsuffix; # 5 digit suffix
}

sub usage
{
	print "Usage:$0 [options]\n";
	print "      -i, --input    : [required] string\tfasta-file(s) of DNA reads\n";
	print "      -g, --gene     : [optional] integer\tgene predictor (0:FragGeneScan, 1:MetaGeneAnnotator, default:0)\n";
	print "      -n, --ncpu     : [optional] integer\tnumber of cpus (default:1)\n";
	print "      -o, --orfdir   : [optional] string\toutput directory (default:.)\n";
	print "      -f, --complete : [optional] integer\tFGS complete genome option (0 or 1, default:0)\n";
	print "      -t, --train    : [optional] string\tFGS train option (complete, sanger_5, sanger_10, 454_10, 454_30, illumina_5, or illumina_10, default: illumina_10)\n";
	print "      -d, --dropindel: [optional] integer\tignore FGS ORF prediction with indels (0 or 1, default:1)\n"; 
	print "      -m, --mgagenome: [optinoal] string\tMGA genome option (s/m, default:m)\n";
	print "      -L, --line     : [optional] integer\tnumber of line of splitted files (default:200000 - 100K fasta)\n";
	print "      -S, --suffix   : [optional] integer\tlength of suffix length of splitted files (default:5)\n";

	exit;
}
