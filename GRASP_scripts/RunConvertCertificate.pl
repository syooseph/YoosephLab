#!/usr/bin/perl -w
use strict;

my $script_dir = "/usr/local/depot/projects/SPA/czhong/Works/GuidedAssemble/Scripts";
my $result_dir = "/usr/local/depot/projects/SPA/czhong/Grasp_stable/Amphora2_marine";

chdir $result_dir or die "Cannot change directory: $!\n";
foreach(<*.fa.aln>)       {
  print "$_\n";
  my $aln_file = $_;
  $aln_file =~ /(.*)\.fa\.aln/;
  my $l_tag = $1;
  # extracting the assembled sequence from the alignment file
  system "perl $script_dir/ConvertCertificateToMFASTA.pl $aln_file >$l_tag.asm.faa";
}
