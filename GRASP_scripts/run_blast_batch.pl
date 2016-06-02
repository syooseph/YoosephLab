#!/usr/bin/perl -w
use strict;

chdir "/usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace" or die "Cannot change directory: $!\n";

foreach(<deh*.fa>){
    my $file_stem = $_;
    $file_stem =~ s/\.fa//g;
    print "running blast for $file_stem\n";
    system "/usr/local/packages/ncbi-blast+-2.2.28/bin/tblastn -query /usr/local/depot/projects/SPA/czhong/Grasp/WorkSpace/$_ -db /usr/local/depot/projects/SPA/czhong/Data/all.core.fna -num_threads 4 -evalue 10 -outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart qend sframe sstart send evalue bitscore\" -max_target_seqs 20000000 >/usr/local/depot/projects/SPA/czhong/Grasp/Results/$file_stem.bst";
}
