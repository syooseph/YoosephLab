#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';

my $q_dir = shift;	# the directory that contains the query files
my $q_label_prefix = shift;	# the query files' prefix id
my $db = shift;	# the db to search against
my $out_dir = shift;	# the output directory

$q_dir = abs_path($q_dir);
$out_dir = abs_path($out_dir);
$db = abs_path($db);

foreach(<$q_dir/$q_label_prefix*.fa>){
    /.*\/(.*)\.fa/;
    my $file_stem = $1;
    print "running psi-blast for $file_stem\n";
    system "/usr/local/packages/ncbi-blast+-2.2.28/bin/psiblast -query $_ -db $db -evalue 1000 -num_threads 8 -outfmt 7 -num_iterations 3 -max_target_seqs 20000000 >$out_dir/$file_stem.psibst";
}
