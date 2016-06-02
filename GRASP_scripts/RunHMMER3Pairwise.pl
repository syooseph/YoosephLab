#!/usr/bin/perl -w
use strict;

#my $file = shift;  # this is expected to be a kegg gene-pfam file

my $hmmer_exe = "/usr/local/packages/hmmer3/bin/hmmsearch";
my $hmm_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/marine_glycolysis_hmm";
my $target_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Data/Marine_Core_Genomes/Annotation";
my $output_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/core_marine_glycolysis_annot_hmm";

foreach(<$hmm_dir/*.hmm>)  {
  my $query = $_; 
  $query =~ /.*\/(.*)\.hmm/;
  my $query_stem = $1;
  foreach(<$target_dir/*.pep>)  {
    my $target = $_;
    $target =~ /.*\/(.*)\.pep/;
    my $target_stem = $1;
    print "running hmmsearch to search $query_stem against $target_stem ...\n";
    system "$hmmer_exe $query $target >$output_dir/$query_stem-$target_stem.hmm.search";
  }
}
