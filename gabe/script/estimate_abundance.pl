#!/usr/bin/env perl

##====================================================================
## Author: Youngik Yang
## Date: Fri 2015-07-24 01:55:18 PM
## Description:
## Wrapper script of genome abundance estimation
##====================================================================

BEGIN { push @INC, `echo -n $ENV{GABE}/script` }

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use POSIX;
use Cwd;
use timer;


## GABE Path must be set in environment
my $root        = $ENV{GABE};

## global variables and default values
my $help        = undef;
my $bwa_file    = undef;      ## bwa mapping file
my $pairend     = undef;
my $genome_file = undef;      ## genome information file
my $skip_matrix = undef;      ## skip matrix generatioin?
my $ncpus       = 1;          ## no. CPUs 
my $max_iter    = 10000;      ## maximum EM iterations
my $tolerance   = 1.0e-6;     ## EM convergency tolerance
my $equal_start = 1;          ## Same mixing coefficients at start
my $em_verbose  = undef;      ## verbosity
my $bootstrap   = undef;      ## perform bootstrapping?
my $boot_iter   = 100;        ## how many bootstrapping? (typically 1000)
my $sample_rate = 30;         ## percentage of reads to sample (typically 100%)
my $purge_temp  = 1;          ## do not delete temporary file (bam2bed.txt, bed2num.txt)
my $out_dir     = getcwd;     ## output directory

##---------------------------
## Parse command line options
##---------------------------
parseOpts();

system("mkdir -p $out_dir 2>/dev/null") unless -d $out_dir;

##----------------------
## Generate input matrix
##----------------------
my $ts  = time();
my $matrix = "$out_dir/gabe_matrix";
parseBWA();

##--------------------
## Estinamte abundaces
##--------------------
estimate();

timer::printElapsed( $ts, "Done", 0);


sub parseBWA
{
    return if defined $skip_matrix;

    my $cmd =  "$root/bin/gabe_pre -i $bwa_file -g $genome_file -o $out_dir";
    $cmd .= " -p" if defined $pairend;
    
    print STDERR "\nParsing BWA ...\n";
    print STDERR "$cmd\n";

    my $ret = system("$cmd");
    die "Error while parsing BWA:$ret" unless $ret == 0;
    timer::printElapsed( $ts, "BWA parsing completed", 0);
    
    # if ( $purge_temp ) {
    # 	# print STDERR "Temporay files removed:($bed and $num)\n";
    # }
}

sub estimate
{
    die "$matrix does not exist" unless -e $matrix;
    
    my $genome_file = "$out_dir/gabe_lengths";
    my $rname_file  = "$out_dir/gabe_reads";
    my $log = "$out_dir/gabe_estimate.log";
    my $exe = "$root/bin/gabe_run";
    my $cmd = "$exe -i $matrix -d 0 -l $genome_file -q $rname_file -n $ncpus -e $equal_start -t $tolerance -o $out_dir";
    $cmd .= " -v" if defined $em_verbose;
    $cmd .= " -b -B $boot_iter -P $sample_rate" if defined $bootstrap;
    $cmd .= " >& $log";
    print STDERR "\nEstimating genome abundances ...\n";
    print STDERR "$cmd\n";
    my $ret = system("$cmd");
    die "Error during EM:$ret" unless $ret == 0;
    print STDERR "EM and abundance estimation:$log\n";
    timer::printElapsed( $ts, "Abundance estimated", 0);            
}

sub parseOpts
{
    my $cmd = join(" ", $0, @ARGV);
    
    GetOptions( "help|h|H|?"      => \$help,
		"bwafile|b=s"     => \$bwa_file,
		"pairend|p"       => \$pairend,
		"genomes|g=s"     => \$genome_file,
		"outdir|o=s"      => \$out_dir,
		"cpu|n=i"         => \$ncpus,
		"maxiter|m=i"     => \$max_iter,
		"equalmix|e=i"    => \$equal_start,
		"tolerance|t=f"   => \$tolerance,
		"verbose|v"       => \$em_verbose,
		"bootsrap|B"      => \$bootstrap,
		"bootiter|T=i"    => \$boot_iter,
		"percentage|P=f"  => \$sample_rate,
		"nomatrix|M"      => \$skip_matrix,
		#"remove_temp|r=i" => \$purge_temp
	      );


    die getUsage() if defined $help;

    if ( ! defined $skip_matrix ) {
	if ( ! defined $bwa_file ) {
	    print "BWA file is required\n";
	    die getUsage();
	} else {
	    die "$bwa_file does not exist\n" unless -e $bwa_file;
	}
    }

    if ( ! defined $genome_file ) {
	die "Genome information file is required";
    }

    timer::printLocalTime();
    print STDERR "Machine:$ENV{HOSTNAME}\n";
    print STDERR "Command:$cmd\n";
}

sub getUsage
{
    my $usage = "$0 [options]\n";
    $usage .= "\t#======================================================================\n";
    $usage .= "\t#Key             Type        Feature     Description\n";
    $usage .= "\t#======================================================================\n";    
    $usage .= "\t-h, --help       on/off      Optional    Program usage\n";
    $usage .= "\t-o, --outdir     string      Optional    output direcory (default:.)\n";
    $usage .= "\t#----------------------------------------------------------------------\n";
    $usage .= "\t# Mapping and genome parameters\n";
    $usage .= "\t                                         Each line should be a pair of a genome name and its size.\n";
    $usage .= "\t-b, --bwafile    string      Required    BWA mapping file\n";
    $usage .= "\t-p, --pairend    on/off      Optional    Turn on in case of pair-end reads\n";
    $usage .= "\t                                         It is required for the first run which generates an input matrix.\n";
    $usage .= "\t                                         Once a matrix generated, this option is ignored for re-estimation (See -M option).\n";
    $usage .= "\t-g, --genomes    string      Required    Genome information file (See the format in the following line)\n";
    $usage .= "\t                                         contig-id\tgenome-id\toptional-description\n";
    $usage .= "\t                                         [NOTE] contig-id must be identical to sequence-id in BWA mapping\n";    
    $usage .= "\t-M, --nomatrix   on/off      Optional    Skip matrix generation (default:off)\n";
    $usage .= "\t                                         In order to perform EM with different parameters\n";
    $usage .= "\t#----------------------------------------------------------------------\n";
    $usage .= "\t# Expectation Maximization parameters\n";
    $usage .= "\t-n, --cpu        integer     Optional    Number of CPUs (default:$ncpus)\n";
    $usage .= "\t-m, --maxiter    integer     Optional    Maximum number of EM iterations (default:$max_iter)\n";
    $usage .= "\t-t, --tolerance  number      Optional    EM convergence tolerance (default:$tolerance)\n";
    $usage .= "\t                                         In case of performing bootstrapping later\n";
    $usage .= "\t-e, --equalmix   integer     Optioanl    Same coefficients at start? (default:1)\n";
    $usage .= "\t-v, --verbose    on/off      Optional    Show EM detail (default:off)\n";
    $usage .= "\t#----------------------------------------------------------------------\n";    
    $usage .= "\t# Bootstrapping parameters\n";
    $usage .= "\t# Bootstrapping can be computationally intensive.\n";
    $usage .= "\t-B, --bootstrap  on/off      Optional    Perform bootstrapping (default:false)\n";
    $usage .= "\t-T, --bootiter   integer     Optional    Bootstrapping iteration (default:$boot_iter)\n";
    $usage .= "\t-P, --percentage number      Optional    Bootstrapping sampling rate (default:$sample_rate)\n";
    $usage .= "\t#----------------------------------------------------------------------\n";
    $usage .= "\te.g.: $0 -b bwa_file\n";
    return $usage;
}
