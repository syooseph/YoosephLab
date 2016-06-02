#!/usr/bin/perl -w
use strict;

chdir "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace" or die "Cannot change directory: $!\n";

foreach(<NT05HA*.fa>){
    my $file_stem = $_;
    $file_stem =~ s/\.fa//g;
    print "running blastp for $file_stem\n";
    system "/usr/local/packages/ncbi-blast+-2.2.28/bin/blastp -query /usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/$_ -db /usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/saliva.faa -evalue 1000 -outfmt 7 -max_target_seqs 20000000 >/usr/local/depot/projects/SPA/czhong/Grasp/Results/$file_stem.bst";
}
