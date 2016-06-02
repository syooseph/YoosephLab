#!/usr/bin/perl -w
use strict;

my $file = shift;  # this is expected to be a kegg gene-pfam file

my $input_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/HMP_saliva/PSIBLAST/derived";
my $results_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Results/HMP_saliva/PSIBLAST/hmm_results";
my $hmmer_exe = "/usr/local/packages/hmmer3/bin/hmmsearch";
my $hmm_dir = "/usr/local/depot/projects/SPA/czhong/Data/Pfam_HMM";

open my $IN, "<$file" or die "Cannot open file: $!\n";
while(<$IN>)  {
  chomp;
  /\:(.*)\s+.*\:(.*)/;
  my $gene = $1;
  my $pfam = $2;
  if(-e "$input_dir/$gene.asm.faa")  {
    print "running hmmsearch on $_ ...\n";
    system "$hmmer_exe $hmm_dir/$pfam.hmm $input_dir/$gene.asm.faa >$results_dir/$gene-$pfam.hmm.search";
  }
}
close $IN;
