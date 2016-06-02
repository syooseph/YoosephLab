#!/usr/bin/perl -w
use strict;

chdir "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace" or die "Cannot change directory: $!\n";

foreach(<REF*.fa>){
    my $file_stem = $_;
    $file_stem =~ s/\.fa//g;
    print "running fastm for $file_stem\n";
    system "/home/czhong/Tools/fasta-36.3.6/bin/fastm36 /usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/$_ /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/marine.core.reduc_header.faa -s BL62 -E 1000000 -C 100 -d 10000000000 >/usr/local/depot/projects/SPA/czhong/Grasp/Results/$file_stem.frag.fastm";
}
