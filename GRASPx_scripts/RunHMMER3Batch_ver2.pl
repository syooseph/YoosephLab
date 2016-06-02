#!/usr/bin/perl -w
use strict;

my $file = shift;  # this is expected to be a kegg gene-pfam file

my $input = shift;
my $results_dir = shift;
my $hmmer_exe = "/usr/local/packages/hmmer3/bin/hmmsearch";
my $hmm_dir = "/usr/local/depot/projects/SPA/czhong/Data/Pfam_HMM";

open my $IN, "<$file" or die "Cannot open file: $!\n";
my %fam_hash;
while(<$IN>)  {
  chomp;
  /\:(.*)\s+.*\:(.*)/;
  $fam_hash{$2} = 1;
}
close $IN;

foreach(keys %fam_hash)  {
  print "running hmmsearch on $_ ...\n";
  system "$hmmer_exe $hmm_dir/$_.hmm $input >$results_dir/$_.hmm.search";
}
