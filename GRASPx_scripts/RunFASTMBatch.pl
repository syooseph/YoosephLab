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
    print "running fastm for $file_stem\n";
    system "/home/czhong/Tools/fasta-36.3.6/bin/fastm36 $_ $db -s BL62 -E 1000000 -C 100 -d 10000000000 >$out_dir/$file_stem.fastm";
}
